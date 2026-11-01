"""Section 11, one filtering suffices.  (a) Is every element of K_1(3) ONE filtering of rho plus a separable operator?  Test NPT sums of two or three random
filterings of rho: fit by 1 filtering + 80 products versus k filterings + 80 products.
(b) The saved F3 targets: fit by 1, 2, 3 filterings + 60 products, and the overlap of the fitted filtering's
partial-transpose negative vector (a (x) conj b)|Phi> with the negative eigenvector of X^Gamma.
usage: python check_f3_one_filtering.py   (needs PyTorch; about ten minutes)"""
import itertools
import numpy as np, torch
torch.set_default_dtype(torch.float64)
D = 3; N = 9
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2; Rt = torch.tensor(rho, dtype=torch.complex128)
def filt(A, B, m):
    subsets = [S for r in range(m + 1) for S in itertools.combinations(range(m), r)]
    At, Bt = A.reshape([D] + [3] * m), B.reshape([D] + [3] * m); Y = 0
    for S in subsets:
        Sc = [s for s in range(m) if s not in S]
        pa = [0] + [1 + s for s in S] + [1 + s for s in Sc]
        Ap = At.permute(pa).reshape(D, 3 ** len(S), 3 ** len(Sc)); Bp = Bt.permute(pa).reshape(D, 3 ** len(S), 3 ** len(Sc))
        M1 = torch.einsum('axc,eyc->aexy', Ap, Ap.conj()); M2 = torch.einsum('byc,fxc->bfyx', Bp, Bp.conj())
        Y = Y + (-0.5) ** len(S) * torch.einsum('aexy,bfyx->abef', M1, M2).reshape(N, N)
    return Y
def pt(Y): return Y.reshape(D, D, D, D).transpose(0, 3, 2, 1).reshape(N, N) if isinstance(Y, np.ndarray) else Y.reshape(D, D, D, D).permute(0, 3, 2, 1).reshape(N, N)
def fit(X, J, P, RS, seed=0):
    torch.manual_seed(seed); Xt = torch.tensor(X); best = np.inf
    for r in range(RS):
        p = [(0.5 * torch.randn(*s, dtype=torch.complex128)).requires_grad_() for s in ((J, D, 3), (J, D, 3), (P, D), (P, D))]
        def model():
            V = torch.einsum('ka,kb->kab', p[2], p[3]).reshape(P, N)
            Kt = torch.einsum('jab,jcd->jacbd', p[0], p[1]).reshape(J, N, 9)
            return V.T @ V.conj() + (Kt @ Rt @ Kt.conj().transpose(1, 2)).sum(0)
        o = torch.optim.LBFGS(p, lr=1, max_iter=5000, history_size=50, line_search_fn='strong_wolfe', tolerance_grad=1e-14, tolerance_change=1e-16)
        def cl():
            o.zero_grad(); loss = torch.linalg.norm(model() - Xt) ** 2; loss.backward(); return loss
        o.step(cl)
        with torch.no_grad(): best = min(best, float(torch.linalg.norm(model() - Xt) / torch.linalg.norm(Xt)))
    return best
g = np.random.default_rng(5); done = 0
while done < 6:
    k = 2 + done % 2
    S = 0
    for _ in range(k):
        a = g.standard_normal((3, 3)) + 1j * g.standard_normal((3, 3)); b = g.standard_normal((3, 3)) + 1j * g.standard_normal((3, 3))
        K = np.kron(a, b); S = S + K @ rho @ K.conj().T
    S = (S + S.conj().T) / 2; S /= np.trace(S).real
    m = np.linalg.eigvalsh(pt(S))[0]
    if m > -0.01: continue
    done += 1
    print(f"sum of {k} filterings, PT-min {m:+.4f}: 1 filt + 80 prod {fit(S, 1, 80, 6):.2e}; {k} filt + 80 prod {fit(S, k, 80, 4):.2e}", flush=True)

def fit2(X, J, P, RS, seed=0):
    torch.manual_seed(seed); Xt = torch.tensor(X); best = (np.inf, None)
    for r in range(RS):
        p = [(0.5 * torch.randn(*s, dtype=torch.complex128)).requires_grad_() for s in ((J, D, 3), (J, D, 3), (P, D), (P, D))]
        def model():
            V = torch.einsum('ka,kb->kab', p[2], p[3]).reshape(P, N)
            Kt = torch.einsum('jab,jcd->jacbd', p[0], p[1]).reshape(J, N, 9)
            return V.T @ V.conj() + (Kt @ Rt @ Kt.conj().transpose(1, 2)).sum(0)
        o = torch.optim.LBFGS(p, lr=1, max_iter=5000, history_size=50, line_search_fn='strong_wolfe', tolerance_grad=1e-14, tolerance_change=1e-16)
        def cl():
            o.zero_grad(); loss = torch.linalg.norm(model() - Xt) ** 2; loss.backward(); return loss
        o.step(cl)
        with torch.no_grad():
            rr = float(torch.linalg.norm(model() - Xt) / torch.linalg.norm(Xt))
        if rr < best[0]: best = (rr, [q.detach().numpy().copy() for q in p])
    return best
for tag in ["M3_s31", "M4_s31", "M5_s31", "M6_s34", "M6_s35"]:
    Fl = np.load(f"adv_f3_filters_{tag}.npz"); M = int(tag[1])
    X = filt(torch.tensor(Fl["A"]), torch.tensor(Fl["B"]), M).numpy(); X = (X + X.conj().T) / 2; X /= np.trace(X).real
    ev, U = np.linalg.eigh(pt(X)); z = U[:, 0]
    out = []
    for J in (1, 2, 3):
        rr, pp = fit2(X, J, 60, 6)
        out.append(rr)
        if J == 1:
            a, b = pp[0][0], pp[1][0]
            v = np.kron(a, b.conj()) @ np.eye(3).reshape(9)      # (a (x) conj b)|Phi>
            K = np.kron(a, b); G = K @ rho @ K.conj().T
            share = np.trace(G).real
            align = abs(np.vdot(v / np.linalg.norm(v), z))
    print(f"{tag}: PT-min {ev[0]:+.4f}; fit 1 filt + 60 prod {out[0]:.2e}, 2 filt {out[1]:.2e}, 3 filt {out[2]:.2e}; "
          f"single filtering carries {share:.3f} of the trace; |<v|z>| = {align:.4f}", flush=True)
