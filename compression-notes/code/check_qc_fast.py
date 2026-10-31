"""Fast adversarial test of the compression conjecture QC at n = 5, 6.
QC: for every rank-2 C there are a copy s and a unit w with Delta_s(w) = q_n(C) - q_n(C_{s,w}) >= 0,
C_{s,w} = (I - ww*)_s C (I - ww*)_s.  Adversary: minimise F(C) = max_{s,w} Delta_s(w) / ||C||^2 over C = Y Z^*.
For a site s, write the 3x3 blocks X_ab = <a|_s C |b>_s and Q[(a,b),(c,d)] = q_{n-1}(X_ab, X_cd) (9x9). Then for the
projector P = I - ww*:  q_n(C_P) = sum_{jk} q(X'_jk) - (1/2) q(T'),  X' = (P (x) P^T) X,  T' = sum_j X'_jj, and q_n(C)
the same with P = I.  The inner maximisation over w uses these 9x9 formulas (cheap); the outer step uses the
Danskin gradient at the maximiser (autograd with (s, w) fixed)."""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)


def Lpow(C, m):
    """L^{(x)m} on C of shape (3,)*m + (3,)*m; L(X) = X - Tr(X) I / 2."""
    for i in range(m):
        tr = torch.diagonal(C, dim1=i, dim2=m + i).sum(-1)
        eye = torch.eye(3, dtype=C.dtype).reshape([3 if k in (i, m + i) else 1 for k in range(2 * m)])
        C = C - 0.5 * tr.unsqueeze(i).unsqueeze(m + i) * eye
    return C


def site_Q(Y, Z, n, s):
    """9x9 matrix Q[(a,b),(c,d)] = q_{n-1}(X_ab, X_cd) at site s."""
    d1 = 3 ** (n - 1)
    Ys = Y.reshape([3] * n + [2]).movedim(s, 0).reshape(3, d1, 2)
    Zs = Z.reshape([3] * n + [2]).movedim(s, 0).reshape(3, d1, 2)
    X = torch.einsum('aik,bjk->abij', Ys, Zs.conj()).reshape(9, d1, d1)      # X_ab = Y_a Z_b^*
    LX = torch.stack([Lpow(X[r].reshape([3] * (2 * (n - 1))), n - 1).reshape(d1, d1) for r in range(9)])
    return torch.einsum('aij,bij->ab', X.conj(), LX)                         # sesquilinear, conj in first


def qn_from_Q(Q, P):
    """q_n of the compression by the 3x3 matrix P at the site (P = I gives q_n(C))."""
    M = torch.kron(P, P.transpose(0, 1).contiguous())                                                   # acts on (a,b) indices: X' = M X
    S1 = torch.einsum('ra,rb,ab->', M.conj(), M, Q)                          # sum_{jk} q(X'_jk)
    v = P.transpose(0, 1).contiguous().reshape(-1)                                                      # T' = sum_{(c,d)} P_{dc} X_cd
    return (S1 - 0.5 * (v.conj() @ Q @ v)).real


def qn_batch(Q, P):
    """q_n of the compressions by a batch of 3x3 matrices P (B, 3, 3)."""
    B = P.shape[0]
    M = torch.einsum('bja,bkc->bjkac', P, P.transpose(1, 2)).reshape(B, 9, 9)   # (P (x) P^T)
    S1 = torch.einsum('bra,brc,ac->b', M.conj(), M, Q)
    v = P.transpose(1, 2).reshape(B, 9)
    return (S1 - 0.5 * torch.einsum('ba,ac,bc->b', v.conj(), Q, v)).real


def inner_max(Qs, n, starts=64, iters=150, gen=None):
    best = (-np.inf, None, None)
    I3 = torch.eye(3, dtype=torch.complex128)
    for s, Q in enumerate(Qs):
        q0 = qn_from_Q(Q, I3)
        W = torch.randn(starts, 3, 2, generator=gen).requires_grad_()
        opt = torch.optim.Adam([W], lr=0.05)
        for it in range(iters):
            w = torch.complex(W[..., 0], W[..., 1]); w = w / torch.linalg.norm(w, dim=1, keepdim=True)
            P = I3[None] - w[:, :, None] * w.conj()[:, None, :]
            vals = q0 - qn_batch(Q, P)
            opt.zero_grad(); (-vals.sum()).backward(); opt.step()
        with torch.no_grad():
            w = torch.complex(W[..., 0], W[..., 1]); w = w / torch.linalg.norm(w, dim=1, keepdim=True)
            P = I3[None] - w[:, :, None] * w.conj()[:, None, :]
            vals = q0 - qn_batch(Q, P)
            k = int(torch.argmax(vals))
            if vals[k].item() > best[0]: best = (vals[k].item(), s, w[k].clone())
    return best


def run(n, seed, outer=300, kind='random'):
    D = 3 ** n
    gen = torch.Generator().manual_seed(seed)
    if kind == 'random':
        Y = torch.randn(D, 2, dtype=torch.complex128, generator=gen); Z = torch.randn(D, 2, dtype=torch.complex128, generator=gen)
    else:                                                                    # near the zero |0^n><0^n| - (-1)^n |1^n><1^n|
        Y = torch.zeros(D, 2, dtype=torch.complex128); Z = torch.zeros(D, 2, dtype=torch.complex128)
        i0 = 0; i1 = sum(3 ** k for k in range(n))
        Y[i0, 0] = 1; Y[i1, 1] = 1; Z[i0, 0] = 1; Z[i1, 1] = -((-1) ** n)
        Y = Y + 0.15 * torch.randn(D, 2, dtype=torch.complex128, generator=gen) / np.sqrt(D)
        Z = Z + 0.15 * torch.randn(D, 2, dtype=torch.complex128, generator=gen) / np.sqrt(D)
    P = [t.clone().requires_grad_() for t in (Y.real, Y.imag, Z.real, Z.imag)]
    opt = torch.optim.Adam(P, lr=0.01)
    I3 = torch.eye(3, dtype=torch.complex128)
    hist = []
    for it in range(outer):
        Yc, Zc = torch.complex(P[0], P[1]), torch.complex(P[2], P[3])
        C2 = torch.trace((Yc.conj().T @ Yc) @ (Zc.conj().T @ Zc)).real            # ||Y Z^*||^2 = Tr[(Y*Y)(Z*Z)]
        with torch.no_grad():
            Qs = [site_Q(Yc, Zc, n, s).detach() for s in range(n)]
        F, s_star, w_star = inner_max(Qs, n, gen=gen)
        Q = site_Q(Yc, Zc, n, s_star)
        Pw = I3 - torch.outer(w_star, w_star.conj())
        val = (qn_from_Q(Q, I3) - qn_from_Q(Q, Pw)) / C2
        qn = (qn_from_Q(Q, I3) / C2).item()
        opt.zero_grad(); val.backward(); opt.step()
        hist.append((val.item(), qn))
        if it % 50 == 0 or it == outer - 1:
            print(f"  n = {n}, {kind} seed {seed}, it {it}: max Delta/||C||^2 = {val.item():+.4e}, q_n/||C||^2 = {qn:+.4e}", flush=True)
    return min(h[0] for h in hist)


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for seed, kind in [(0, 'random'), (1, 'nearzero'), (2, 'random'), (3, 'nearzero')]:
        m = run(n, seed, outer=int(sys.argv[2]) if len(sys.argv) > 2 else 300, kind=kind)
        print(f"n = {n}, {kind} seed {seed}: smallest max-Delta seen {m:+.4e}", flush=True)
