"""Section 11: certified F3 instances (search step).  For saved adversary filters A, B : (C^3)^{(x)M} -> C^3, find nonnegative weights and
atoms of K_1(3) (filterings (a (x) b) rho (a (x) b)^* with a, b in M_3, and products) such that X - sum w_i G_i is close
to (eps Tr X / 9) I; check_cone_certificate_verify.py then certifies the remainder separable (Gurvits-Barnum ball, N = 9).
usage: python check_f3_certificate.py M seed eps rounds   (reads adv_f3_filters_M{M}_s{seed}.npz; needs PyTorch, SciPy)
then: python check_cone_certificate_verify.py adv_f3_filters_M{M}_s{seed}.npz f3_cert_atoms_M{M}_s{seed}.npz M 1"""
import sys, itertools, numpy as np, torch
from scipy.optimize import nnls
torch.set_default_dtype(torch.float64)
M, SEED, EPS, ROUNDS = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
D = 3; N = D * D
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2
Rt = torch.tensor(rho, dtype=torch.complex128)
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
def pt(Y): return Y.reshape(D, D, D, D).permute(0, 3, 2, 1).reshape(N, N)
Fl = np.load(f"adv_f3_filters_M{M}_s{SEED}.npz"); A, B = Fl["A"], Fl["B"]
X = filt(torch.tensor(A), torch.tensor(B), M).numpy(); X = (X + X.conj().T) / 2
Xn = X / np.trace(X).real
print(f"F3 target M = {M}, seed {SEED}: PT-min/Tr {np.linalg.eigvalsh(pt(torch.tensor(Xn)).numpy())[0]:+.4f}, "
      f"min eigenvalue/Tr {np.linalg.eigvalsh(Xn)[0]:+.3e}; eps = {EPS:.1e}", flush=True)
T = Xn - EPS / N * np.eye(N)
atoms = []
def add(kind, p, q):
    if kind == 'f':
        K = np.kron(p, q); Mx = K @ rho @ K.conj().T
    else:
        v = np.kron(p, q); Mx = np.outer(v, v.conj())
    Mx = (Mx + Mx.conj().T) / 2
    atoms.append((kind, p.copy(), q.copy(), Mx / np.trace(Mx).real))
def param_fit(J, P, RS, seed):
    torch.manual_seed(seed); Xt = torch.tensor(T, dtype=torch.complex128); best = (np.inf, None)
    for r in range(RS):
        p = [(0.4 * torch.randn(*s, dtype=torch.complex128)).requires_grad_() for s in ((J, D, 3), (J, D, 3), (P, D), (P, D))]
        def model():
            V = torch.einsum('ka,kb->kab', p[2], p[3]).reshape(P, N)
            Kt = torch.einsum('jab,jcd->jacbd', p[0], p[1]).reshape(J, N, 9)
            return V.T @ V.conj() + (Kt @ Rt @ Kt.conj().transpose(1, 2)).sum(0)
        o = torch.optim.LBFGS(p, lr=1, max_iter=6000, history_size=50, line_search_fn='strong_wolfe', tolerance_grad=1e-15, tolerance_change=1e-18)
        def cl():
            o.zero_grad(); loss = torch.linalg.norm(model() - Xt) ** 2; loss.backward(); return loss
        o.step(cl)
        with torch.no_grad():
            rr = float(torch.linalg.norm(model() - Xt) / torch.linalg.norm(Xt))
            if rr < best[0]: best = (rr, [q.detach().numpy().copy() for q in p])
    return best
for (J, P, RS, sd) in [(16, 80, 4, 1), (32, 160, 3, 2)]:
    rr, pp = param_fit(J, P, RS, sd)
    print(f"seed fit {J} + {P}: {rr:.2e}", flush=True)
    for j in range(J): add('f', pp[0][j], pp[1][j])
    for k in range(P): add('p', pp[2][k], pp[3][k])
g = np.random.default_rng(1)
for _ in range(400):
    add('f', g.standard_normal((D, 3)) + 1j * g.standard_normal((D, 3)), g.standard_normal((D, 3)) + 1j * g.standard_normal((D, 3)))
    add('p', g.standard_normal(D) + 1j * g.standard_normal(D), g.standard_normal(D) + 1j * g.standard_normal(D))
def vecr(Mx):
    iu = np.triu_indices(Mx.shape[0], 1)
    return np.concatenate([np.diag(Mx).real, np.sqrt(2) * Mx[iu].real, np.sqrt(2) * Mx[iu].imag])
torch.manual_seed(123)
nT = np.linalg.norm(T)
for rd in range(ROUNDS):
    Amat = np.stack([vecr(a[3]) for a in atoms], 1); coef, _ = nnls(Amat, vecr(T), maxiter=500000)
    Rres = T - sum(c * a[3] for c, a in zip(coef, atoms) if c > 0); nR = np.linalg.norm(Rres) / nT
    if rd % 5 == 0 or nR < 1e-14: print(f"round {rd}: atoms {len(atoms)}, active {int(np.sum(coef > 0))}, residual {nR:.3e}", flush=True)
    if nR < 1e-14: break
    Rt_ = torch.tensor(Rres, dtype=torch.complex128); new = []
    for r in range(24):
        a = torch.randn(D, 3, dtype=torch.complex128, requires_grad=True); b = torch.randn(D, 3, dtype=torch.complex128, requires_grad=True)
        o = torch.optim.Adam([a, b], lr=0.03)
        for it in range(300):
            K = torch.kron(a, b); Yg = K @ Rt @ K.conj().T; v = -torch.trace(Rt_ @ Yg).real / torch.trace(Yg).real
            o.zero_grad(); v.backward(); o.step()
        if -v.item() > 0: new.append(('f', a.detach().numpy(), b.detach().numpy()))
        x = torch.randn(D, dtype=torch.complex128, requires_grad=True); y = torch.randn(D, dtype=torch.complex128, requires_grad=True)
        o = torch.optim.Adam([x, y], lr=0.03)
        for it in range(200):
            w = torch.kron(x, y); v = -(w.conj() @ Rt_ @ w).real / (w.conj() @ w).real
            o.zero_grad(); v.backward(); o.step()
        if -v.item() > 0: new.append(('p', x.detach().numpy(), y.detach().numpy()))
    for kind, p, q in new: add(kind, p, q)
act = np.where(coef > 0)[0]
zf, zp = np.zeros((D, 3), complex), np.zeros(D, complex)
np.savez(f"f3_cert_atoms_M{M}_s{SEED}.npz", kinds=np.array([atoms[i][0] for i in act]),
         P=np.array([atoms[i][1] if atoms[i][0] == 'f' else zf for i in act]),
         Q=np.array([atoms[i][2] if atoms[i][0] == 'f' else zf for i in act]),
         x=np.array([atoms[i][1] if atoms[i][0] == 'p' else zp for i in act]),
         y=np.array([atoms[i][2] if atoms[i][0] == 'p' else zp for i in act]),
         coef=coef[act], eps=EPS, residual=nR)
print(f"saved {len(act)} active atoms; final residual {nR:.3e}", flush=True)
