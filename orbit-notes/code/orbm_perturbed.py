import numpy as np, sys, warnings
warnings.filterwarnings("ignore")
sys.argv = ["x", "2", "0", "0"]; exec(open("orbm_descent.py").read().split("for k in range(starts):")[0])
base = "../data/"
rng2 = np.random.default_rng(11)
for name in ("variety_critical_n2_low", "variety_critical_n3_golden", "variety_critical_n3_low"):
    X = np.load(base + name + ".npy"); D = X.shape[0]; n = int(round(np.log(D) / np.log(3)))
    Uu, S, Vh = np.linalg.svd(X); u0 = Uu[:, :2].T.copy()
    fobj = lambda x: m_plane(np.array([filt(x) @ v for v in u0]))
    m0 = fobj(np.zeros(9 * n)); out = []
    for k in range(3):
        x0 = 0.02 * rng2.normal(size=9 * n)
        res = minimize(fobj, x0, method='L-BFGS-B', bounds=[(-8, 8)] * (9 * n), options={'maxiter': 3000, 'gtol': 1e-11, 'ftol': 1e-14})
        uf = np.array([filt(res.x) @ v for v in u0]); out.append(f"{fobj(x0):.6f} -> {res.fun:.6f} (ratio {thin_ratio(uf):.1e})")
    print(f"{name} (n={n}): m at the critical point {m0:.6f}; perturbed descents: " + "; ".join(out), flush=True)
