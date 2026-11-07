import numpy as np, sys
from scipy.optimize import minimize
n = int(sys.argv[1])
d = np.load(f'mz{n}.npz'); Mz = {tuple(int(x) for x in k.split('_')): d[k] for k in d.files}
fixed = sum(Mz[(s, 0)] for s in range(1, n + 1)) + 0.5 * sum(Mz[(s, s)] for s in range(1, n // 2 + 1))
pairs = [(s, t) for t in range(1, n) for s in range(t + 1, n + 1) if s + t <= n]
greedy = fixed + sum(Mz[p] for p in pairs)
stack = np.vstack([greedy] + [Mz[p] - Mz[(p[1], p[0])] for p in pairs])
u, sv, vt = np.linalg.svd(stack); r = int(np.sum(sv > 1e-9 * sv[0])); Vp = vt[:r].T
F0 = Vp.T @ fixed @ Vp; Fp = [Vp.T @ Mz[p] @ Vp for p in pairs]; Fm = [Vp.T @ Mz[(p[1], p[0])] @ Vp for p in pairs]
def lmin(ph):
    M = F0 + sum(ph[i]*Fp[i] + (1-ph[i])*Fm[i] for i in range(len(pairs)))
    return np.linalg.eigvalsh((M+M.T)/2)[0]
for c in [0.5, 0.8, 1.0, 1.1, 1.2, 1.3, 1.5, 2.0]:
    print(f"uniform phi={c}: lambda_min = {lmin(np.full(len(pairs), c)):.4e}")
best = None
rng = np.random.default_rng(0)
for trial in range(6):
    x0 = np.ones(len(pairs)) + 0.3*rng.standard_normal(len(pairs))
    res = minimize(lambda p: -lmin(p), x0, method='Nelder-Mead', options={'maxiter': 4000, 'xatol': 1e-6, 'fatol': 1e-8})
    if best is None or res.fun < best.fun: best = res
    print(f"trial {trial}: lambda_min = {-res.fun:.4e}", flush=True)
print("best lambda_min found:", -best.fun, dict(zip(pairs, np.round(best.x, 3))))
