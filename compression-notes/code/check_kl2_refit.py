"""KL_2 = L_2(9) check for the adversarial 9 x 9 target (adv_lmD_m2_D9_s0.npy): calibrate K_2(9) fits at several budgets
   on a random element of K_2(9) (3 filterings of rho^2 + 30 products), then fit the adversarial target with the same
   budgets.  usage: python check_kl2_refit.py (after check_lmD_adversary.py 2 9 16 80 60 0)"""
import numpy as np, torch
torch.set_default_dtype(torch.float64)
D = 9
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2
def power(k):
    R = np.eye(1)
    for _ in range(k): R = np.kron(R, rho)
    perm = [2 * i for i in range(k)] + [2 * i + 1 for i in range(k)]
    return R.reshape([3] * (4 * k)).transpose(perm + [2 * k + q for q in perm]).reshape(9 ** k, 9 ** k)
R2 = torch.tensor(power(2), dtype=torch.complex128)
def fit(Y, J, P, IT=5000, RS=2):
    Yt = torch.tensor(Y, dtype=torch.complex128); best = np.inf
    for _ in range(RS):
        p = [(0.3 * torch.randn(*s, dtype=torch.complex128)).requires_grad_() for s in ((J, D, 9), (J, D, 9), (P, D), (P, D))]
        def model():
            V = torch.einsum('ka,kb->kab', p[2], p[3]).reshape(P, D * D)
            Kt = torch.einsum('jab,jcd->jacbd', p[0], p[1]).reshape(J, D * D, 81)
            return V.T @ V.conj() + (Kt @ R2 @ Kt.conj().transpose(1, 2)).sum(0)
        opt = torch.optim.LBFGS(p, lr=1, max_iter=IT, history_size=50, line_search_fn='strong_wolfe', tolerance_grad=1e-15, tolerance_change=1e-18)
        def closure():
            opt.zero_grad(); loss = torch.linalg.norm(model() - Yt) ** 2; loss.backward(); return loss
        opt.step(closure)
        with torch.no_grad(): best = min(best, float(torch.linalg.norm(model() - Yt) / torch.linalg.norm(Yt)))
    return best
g = np.random.default_rng(3)
a = g.standard_normal((3, D, 9)) + 1j * g.standard_normal((3, D, 9)); b = g.standard_normal((3, D, 9)) + 1j * g.standard_normal((3, D, 9))
Yc = sum(np.kron(a[i], b[i]) @ power(2) @ np.kron(a[i], b[i]).conj().T for i in range(3))
xs = g.standard_normal((30, D)) + 1j * g.standard_normal((30, D)); ys = g.standard_normal((30, D)) + 1j * g.standard_normal((30, D))
Yc = Yc + sum(np.outer(np.kron(x, y), np.kron(x, y).conj()) for x, y in zip(xs, ys)); Yc /= np.trace(Yc).real
Ya = np.load("adv_lmD_m2_D9_s0.npy")
for (J, P) in [(16, 80), (32, 200), (48, 400)]:
    print(f"budget {J} + {P}: calibration (element of K_2(9)) {fit(Yc, J, P):.2e} | adversarial target {fit(Ya, J, P):.2e}", flush=True)
