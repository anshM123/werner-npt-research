"""Section 9: certified membership of the 6 x 6 chain targets in K_2(6) (search step; the check is
check_chain_certificate_verify.py).
Target: X = (A (x) B) rho^{(x)M} (A (x) B)^*, with A, B : (C^3)^{(x)M} -> C^6 the chain filters (the two-copy NPT 6 x 6
filtering from the seed-4 gradient path, padded with e_0 on the extra copies and moved by delta = 0.1 toward random
M-copy filters, generator seed 9).  A, B are float64, hence exact dyadic rationals; they are saved.
Column generation for K_2(6) on the shifted target T = X/Tr X - (eps/36) I: the weights are found for T, so that the
remainder X - sum w_i G_i is close to (eps Tr X/36) I and can be certified separable (check_chain_certificate_verify.py).
Atoms: filterings (a (x) b) rho^{(x)2} (a (x) b)^* (a, b : C^9 -> C^6) and products |x y><x y|, saved by parameters.
usage: python check_chain_certificate.py M eps rounds [seed_atoms.npz | fit_terms.pkl]   (needs PyTorch and SciPy)
The four-copy atoms were found in two steps: eps = 0 (41 rounds, 1296 active atoms; the saved file is then the
seed), and eps = 1e-6 seeded with those atoms (4 rounds).  The three-copy run is seeded with the fit terms
chain_seeds_M3.pkl of check_chain_six.py."""
import sys, itertools, numpy as np, torch
from scipy.optimize import nnls
torch.set_default_dtype(torch.float64)
M, EPS, ROUNDS = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
SEEDFILE = sys.argv[4] if len(sys.argv) > 4 else None
D = 6
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2
R2 = np.kron(rho, rho).reshape([3] * 8).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(81, 81)
R2t = torch.tensor(R2, dtype=torch.complex128)
def filt(A, B, m):
    subsets = [S for r in range(m + 1) for S in itertools.combinations(range(m), r)]
    At, Bt = A.reshape([D] + [3] * m), B.reshape([D] + [3] * m); Y = 0
    for S in subsets:
        Sc = [s for s in range(m) if s not in S]
        pa = [0] + [1 + s for s in S] + [1 + s for s in Sc]
        Ap = At.permute(pa).reshape(D, 3 ** len(S), 3 ** len(Sc)); Bp = Bt.permute(pa).reshape(D, 3 ** len(S), 3 ** len(Sc))
        M1 = torch.einsum('axc,eyc->aexy', Ap, Ap.conj()); M2 = torch.einsum('byc,fxc->bfyx', Bp, Bp.conj())
        Y = Y + (-0.5) ** len(S) * torch.einsum('aexy,bfyx->abef', M1, M2).reshape(D * D, D * D)
    return Y
def pt(Y): return Y.reshape(D, D, D, D).permute(0, 3, 2, 1).reshape(D * D, D * D)
# ---- target filters (same construction as chain_l3six.py and cert_pipeline.py)
torch.manual_seed(4)
A2 = torch.randn(D, 9, dtype=torch.complex128, requires_grad=True); B2 = torch.randn(D, 9, dtype=torch.complex128, requires_grad=True)
opt = torch.optim.Adam([A2, B2], lr=0.01)
for it in range(20000):
    Y = filt(A2, B2, 2); Y = (Y + Y.conj().T) / 2; Y = Y / torch.trace(Y).real
    f = torch.linalg.eigvalsh(pt(Y))[0]
    if f.item() < -0.03: break
    opt.zero_grad(); f.backward(); opt.step()
u = torch.zeros(3 ** (M - 2), dtype=torch.complex128); u[0] = 1
A0 = torch.einsum('ax,y->axy', A2.detach(), u).reshape(D, 3 ** M); B0 = torch.einsum('ax,y->axy', B2.detach(), u).reshape(D, 3 ** M)
g2 = torch.Generator().manual_seed(9)
RA = torch.randn(D, 3 ** M, dtype=torch.complex128, generator=g2); RB = torch.randn(D, 3 ** M, dtype=torch.complex128, generator=g2)
A = (A0 / torch.linalg.norm(A0) + 0.1 * RA / torch.linalg.norm(RA)).numpy(); B = (B0 / torch.linalg.norm(B0) + 0.1 * RB / torch.linalg.norm(RB)).numpy()
np.savez(f"chain_cert_target_M{M}.npz", A=A, B=B)
X = filt(torch.tensor(A), torch.tensor(B), M).numpy(); X = (X + X.conj().T) / 2
Xn = X / np.trace(X).real
print(f"M = {M}: target PT-min/Tr {np.linalg.eigvalsh(pt(torch.tensor(Xn)).numpy())[0]:+.4f}, min eigenvalue/Tr "
      f"{np.linalg.eigvalsh(Xn)[0]:+.3e}; eps = {EPS:.1e}", flush=True)
T = Xn - EPS / (D * D) * np.eye(D * D)
# ---- dictionary
def gen_f(a, b):
    K = np.kron(a, b); Yg = K @ R2 @ K.conj().T; return (Yg + Yg.conj().T) / 2
def gen_p(x, y):
    v = np.kron(x, y); return np.outer(v, v.conj())
atoms = []
def add(kind, p, q):
    Mx = gen_f(p, q) if kind == 'f' else gen_p(p, q)
    atoms.append((kind, p.copy(), q.copy(), Mx / np.trace(Mx).real))
def param_fit(J, P, RS, seed):
    torch.manual_seed(seed); Xt = torch.tensor(T, dtype=torch.complex128); best = (np.inf, None)
    for r in range(RS):
        p = [(0.3 * torch.randn(*s, dtype=torch.complex128)).requires_grad_() for s in ((J, D, 9), (J, D, 9), (P, D), (P, D))]
        def model():
            V = torch.einsum('ka,kb->kab', p[2], p[3]).reshape(P, D * D)
            Kt = torch.einsum('jab,jcd->jacbd', p[0], p[1]).reshape(J, D * D, 81)
            return V.T @ V.conj() + (Kt @ R2t @ Kt.conj().transpose(1, 2)).sum(0)
        o = torch.optim.LBFGS(p, lr=1, max_iter=6000, history_size=50, line_search_fn='strong_wolfe', tolerance_grad=1e-15, tolerance_change=1e-18)
        def cl():
            o.zero_grad(); loss = torch.linalg.norm(model() - Xt) ** 2; loss.backward(); return loss
        o.step(cl)
        with torch.no_grad():
            rr = float(torch.linalg.norm(model() - Xt) / torch.linalg.norm(Xt))
            if rr < best[0]: best = (rr, [q.detach().numpy().copy() for q in p])
    return best
if SEEDFILE and SEEDFILE.endswith(".pkl"):
    import pickle
    for pp in pickle.load(open(SEEDFILE, "rb")):
        for j in range(pp[0].shape[0]): add('f', pp[0][j], pp[1][j])
        for k in range(pp[2].shape[0]): add('p', pp[2][k], pp[3][k])
    print(f"seeded with {len(atoms)} fit terms from {SEEDFILE}", flush=True)
elif SEEDFILE:
    S = np.load(SEEDFILE, allow_pickle=True)
    for k in range(len(S["kinds"])):
        if S["kinds"][k] == 'f': add('f', S["P"][k], S["Q"][k])
        else: add('p', S["x"][k], S["y"][k])
    print(f"seeded with {len(atoms)} atoms from {SEEDFILE}", flush=True)
else:
    for (J, P, RS, sd) in [(16, 80, 3, 1), (32, 160, 2, 2)]:
        rr, pp = param_fit(J, P, RS, sd)
        print(f"seed fit {J} + {P}: {rr:.2e}", flush=True)
        for j in range(J): add('f', pp[0][j], pp[1][j])
        for k in range(P): add('p', pp[2][k], pp[3][k])
g = np.random.default_rng(1)
for _ in range(600):
    add('f', g.standard_normal((D, 9)) + 1j * g.standard_normal((D, 9)), g.standard_normal((D, 9)) + 1j * g.standard_normal((D, 9)))
    add('p', g.standard_normal(D) + 1j * g.standard_normal(D), g.standard_normal(D) + 1j * g.standard_normal(D))
def vecr(Mx):
    iu = np.triu_indices(Mx.shape[0], 1)
    return np.concatenate([np.diag(Mx).real, Mx[iu].real, Mx[iu].imag])
torch.manual_seed(123)
nT = np.linalg.norm(T)
for rd in range(ROUNDS):
    Amat = np.stack([vecr(a[3]) for a in atoms], 1); coef, _ = nnls(Amat, vecr(T), maxiter=500000)
    Rres = T - sum(c * a[3] for c, a in zip(coef, atoms) if c > 0); nR = np.linalg.norm(Rres) / nT
    print(f"round {rd}: atoms {len(atoms)}, active {int(np.sum(coef > 0))}, residual {nR:.3e}", flush=True)
    if nR < 1e-14: break
    Rt_ = torch.tensor(Rres, dtype=torch.complex128); new = []
    for r in range(24):
        a = torch.randn(D, 9, dtype=torch.complex128, requires_grad=True); b = torch.randn(D, 9, dtype=torch.complex128, requires_grad=True)
        o = torch.optim.Adam([a, b], lr=0.03)
        for it in range(400):
            K = torch.kron(a, b); Yg = K @ R2t @ K.conj().T; v = -torch.trace(Rt_ @ Yg).real / torch.trace(Yg).real
            o.zero_grad(); v.backward(); o.step()
        if -v.item() > 0: new.append(('f', a.detach().numpy(), b.detach().numpy()))
        x = torch.randn(D, dtype=torch.complex128, requires_grad=True); y = torch.randn(D, dtype=torch.complex128, requires_grad=True)
        o = torch.optim.Adam([x, y], lr=0.03)
        for it in range(300):
            w = torch.kron(x, y); v = -(w.conj() @ Rt_ @ w).real / (w.conj() @ w).real
            o.zero_grad(); v.backward(); o.step()
        if -v.item() > 0: new.append(('p', x.detach().numpy(), y.detach().numpy()))
    for kind, p, q in new: add(kind, p, q)
act = np.where(coef > 0)[0]
zf, zp = np.zeros((D, 9), complex), np.zeros(D, complex)
np.savez(f"chain_cert_atoms_M{M}.npz", kinds=np.array([atoms[i][0] for i in act]),
         P=np.array([atoms[i][1] if atoms[i][0] == 'f' else zf for i in act]),
         Q=np.array([atoms[i][2] if atoms[i][0] == 'f' else zf for i in act]),
         x=np.array([atoms[i][1] if atoms[i][0] == 'p' else zp for i in act]),
         y=np.array([atoms[i][2] if atoms[i][0] == 'p' else zp for i in act]),
         coef=coef[act], eps=EPS, residual=nR)
print(f"saved {len(act)} active atoms; final residual {nR:.3e}", flush=True)
