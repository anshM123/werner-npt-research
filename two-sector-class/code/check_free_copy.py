"""One free copy (copy 0, a full qutrit) plus sector copies 1..m (A_s perp B_s, dims (2,1), (1,2) or (1,1), random bases).
Checks:
 (i)   the identity 2^m q_n(C) = sum_mu q1(X_mu) + sum_nu q1(Y_nu) + 2^m [qt(P) + qt(Q)] + 2 (-1)^m Re B1(X_0, Y_0)
       for C = [[X, P], [Q, Y]] supported on C^3 (x) (A + B), X_mu = Tr_A[(I (x) sigma_mu) X] (full-weight Pauli strings on
       the 2-dim A_s, identity on 1-dim ones), q1(Z) = |Z|^2 - |Tr Z|^2/2, B1 its bilinear form, qt(P) = |P|^2 - |Tr_0 P|^2/2;
 (ii)  Lemma L: |B1(Z,W)|^2 <= psi(Z) psi(W), psi(Z) = |Z|_1^2 - |Tr Z|^2/2 (random and adversarial, 3 x 3);
 (iii) Theorem A: q_n >= 0 when rank X <= 1, rank Y <= 1, rank P, rank Q <= 2 (adversarial minimum of q_n/|C|^2);
 (iv)  Corollary: q_n(C) >= 0 for rank-two C whose range is spanned by one vector of C^3 (x) A and one of C^3 (x) B,
       with arbitrary co-range (adversarial minimum), and q_n(C) - q_n(C P_V) >= 0 there.
usage: python check_free_copy.py n trials seed"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64); cplx = torch.complex128
n, TRIALS, SEED = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); D = 3 ** n; m = n - 1
torch.manual_seed(SEED); rng = np.random.default_rng(SEED)
w9 = torch.eye(3, dtype=cplx).reshape(9) / np.sqrt(3)
def Phi(C):
    X = C.reshape((3,) * (2 * n)).permute([i for s in range(n) for i in (s, n + s)]).reshape((9,) * n)
    for s in range(n):
        tr = torch.tensordot(w9.conj(), X, dims=([0], [s])); X = X - 1.5 * torch.moveaxis(torch.tensordot(w9, tr, dims=0), 0, s)
    return X.reshape([3, 3] * n).permute([2 * s for s in range(n)] + [2 * s + 1 for s in range(n)]).reshape(D, D)
def q(C): return torch.vdot(C.reshape(-1), Phi(C).reshape(-1)).real
def kron_all(ms):
    out = torch.ones(1, 1, dtype=cplx)
    for M in ms: out = torch.kron(out, M)
    return out
PA = [torch.tensor(x, dtype=cplx) for x in ([[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]])]
q1 = lambda Z: torch.linalg.norm(Z) ** 2 - 0.5 * torch.abs(torch.trace(Z)) ** 2
B1 = lambda Z, W: torch.vdot(Z.reshape(-1), W.reshape(-1)) - 0.5 * torch.trace(Z).conj() * torch.trace(W)
psi = lambda Z: torch.linalg.matrix_norm(Z, ord='nuc') ** 2 - 0.5 * torch.abs(torch.trace(Z)) ** 2
def strings(dims):  # full-weight Pauli strings on 2-dim factors, identity on 1-dim factors
    two = [k for k, d in enumerate(dims) if d == 2]
    out = []
    for mus in itertools.product(range(3), repeat=len(two)):
        ms, it = [], iter(mus)
        for d in dims: ms.append(PA[next(it)] if d == 2 else torch.eye(1, dtype=cplx))
        out.append(kron_all(ms))
    return out
def hat(X, Sg, dA):  # Tr_A[(I (x) Sg) X], X on C^3 (x) A
    X4 = X.reshape(3, dA, 3, dA); return torch.einsum('ag,igja->ij', Sg, X4)
def tr0(P, dA, dB): return torch.einsum('iaib->ab', P.reshape(3, dA, 3, dB))
idmax, Lmax, thA, cor, corgap = 0.0, 0.0, 1e9, 1e9, 1e9
for trial in range(TRIALS):
    dims = [[(2, 1), (1, 2), (1, 1)][rng.integers(3)] for _ in range(m)]
    Us = [torch.tensor(np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))[0]) for _ in range(m)]
    EA = kron_all([torch.eye(3, dtype=cplx)] + [Us[s][:, :dims[s][0]] for s in range(m)])
    EB = kron_all([torch.eye(3, dtype=cplx)] + [Us[s][:, dims[s][0]:sum(dims[s])] for s in range(m)])
    dA, dB = EA.shape[1] // 3, EB.shape[1] // 3
    SA = strings([d[0] for d in dims]); SB = strings([d[1] for d in dims])
    IA = torch.eye(dA, dtype=cplx); IB = torch.eye(dB, dtype=cplx)
    # (i) identity on a random member of any rank
    X = torch.randn(3 * dA, 3 * dA, dtype=cplx); Y = torch.randn(3 * dB, 3 * dB, dtype=cplx)
    P = torch.randn(3 * dA, 3 * dB, dtype=cplx); Qm = torch.randn(3 * dB, 3 * dA, dtype=cplx)
    C = EA @ X @ EA.conj().T + EB @ Y @ EB.conj().T + EA @ P @ EB.conj().T + EB @ Qm @ EA.conj().T
    rhs = sum(q1(hat(X, S, dA)) for S in SA) + sum(q1(hat(Y, S, dB)) for S in SB) \
        + 2 ** m * (torch.linalg.norm(P) ** 2 - 0.5 * torch.linalg.norm(tr0(P, dA, dB)) ** 2
                    + torch.linalg.norm(Qm) ** 2 - 0.5 * torch.linalg.norm(tr0(Qm, dB, dA)) ** 2) \
        + 2 * (-1) ** m * B1(hat(X, IA, dA), hat(Y, IB, dB)).real
    idmax = max(idmax, (abs(2 ** m * q(C) - rhs) / torch.linalg.norm(C) ** 2).item())
    # (iii) Theorem A: rank X, rank Y <= 1, rank P, rank Q <= 2
    def memberA(p):
        k = 0
        def take(r, c):
            nonlocal k
            v = torch.complex(p[k:k + r * c], p[k + r * c:k + 2 * r * c]).reshape(r, c); k += 2 * r * c; return v
        x, y, xp, yp = take(3 * dA, 1), take(3 * dA, 1), take(3 * dB, 1), take(3 * dB, 1)
        U1, V1, U2, V2 = take(3 * dA, 2), take(3 * dB, 2), take(3 * dB, 2), take(3 * dA, 2)
        C = EA @ (x @ y.conj().T) @ EA.conj().T + EB @ (xp @ yp.conj().T) @ EB.conj().T \
            + EA @ (U1 @ V1.conj().T) @ EB.conj().T + EB @ (U2 @ V2.conj().T) @ EA.conj().T
        return C / torch.linalg.norm(C)
    npA = 2 * (2 * 3 * dA + 2 * 3 * dB + 2 * 3 * dA * 2 + 2 * 3 * dB * 2)
    # (iv) corollary: range spanned by uA in C^3 (x) A and uB in C^3 (x) B, co-range arbitrary
    def memberC(p):
        uA = EA @ torch.complex(p[:3 * dA], p[3 * dA:6 * dA]); o = 6 * dA
        uB = EB @ torch.complex(p[o:o + 3 * dB], p[o + 3 * dB:o + 6 * dB]); o += 6 * dB
        f = torch.complex(p[o:o + 2 * D], p[o + 2 * D:o + 4 * D]).reshape(D, 2)
        C = torch.stack([uA, uB], 1) @ f.conj().T; return C / torch.linalg.norm(C)
    npC = 6 * dA + 6 * dB + 4 * D
    for kind, member, npar in (('A', memberA, npA), ('C', memberC, npC)):
        for st in range(2):
            p = torch.randn(npar, requires_grad=True)
            opt = torch.optim.LBFGS([p], lr=1, max_iter=1500, tolerance_grad=1e-14, tolerance_change=1e-18, history_size=50, line_search_fn='strong_wolfe')
            def closure():
                opt.zero_grad(); v = q(member(p)); v.backward(); return v
            opt.step(closure)
            with torch.no_grad(): v = q(member(p)).item(); Cc = member(p)
            if kind == 'A': thA = min(thA, v)
            else:
                cor = min(cor, v)
                PV = EA @ EA.conj().T + EB @ EB.conj().T
                corgap = min(corgap, (q(Cc) - q(Cc @ PV)).item())
# (ii) Lemma L
for trial in range(TRIALS * 20):
    Z = torch.randn(3, 3, dtype=cplx); W = torch.randn(3, 3, dtype=cplx)
    if trial % 2: Z = Z @ torch.diag(torch.tensor([1, 1, 0], dtype=cplx))
    Lmax = max(Lmax, (torch.abs(B1(Z, W)) ** 2 / (psi(Z) * psi(W))).item())
p = torch.randn(36, requires_grad=True); opt = torch.optim.Adam([p], lr=0.01)
for it in range(3000):
    opt.zero_grad(); Z = torch.complex(p[:9], p[9:18]).reshape(3, 3); W = torch.complex(p[18:27], p[27:]).reshape(3, 3)
    v = -(torch.abs(B1(Z, W)) ** 2 / (psi(Z) * psi(W))); v.backward(); opt.step()
Lmax = max(Lmax, -v.item())
print(f"n = {n}, {TRIALS} random sector configurations")
print(f"(i)   identity: max relative error {idmax:.1e}")
print(f"(ii)  Lemma L: max |B1|^2/(psi psi) = {Lmax:.6f} (random + adversarial)")
print(f"(iii) Theorem A: adversarial min of q_n/|C|^2 = {thA:+.2e}")
print(f"(iv)  corollary class: adversarial min of q_n/|C|^2 = {cor:+.2e}; min of q_n(C) - q_n(C P_V) = {corgap:+.2e}")
