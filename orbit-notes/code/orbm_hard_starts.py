"""ORB_m descents started from the range planes of hard operators stored in the repository (copywise squeeze minima,
single-copy minima, two-sided minima, Levi-positive and low critical points of the orbit notes; tangent-class copywise
minima and shared-shift points; general shared-shift points)."""
import numpy as np, sys, glob, warnings
warnings.filterwarnings("ignore")
sys.argv = ["x", "2", "0", "0"]; exec(open("orbm_descent.py").read().split("for k in range(starts):")[0])
base = "../../"
files = sorted(glob.glob(base + "orbit-notes/data/*.npy")) + [base + "tangent-class/data/copywise_squeeze_min_n3.npy", base + "tangent-class/data/copywise_squeeze_min_n4.npy",
         base + "tangent-class/data/shared_shift_pd_T_n3.npy", base + "tangent-class/data/common_shift_pd_n3.npy", base + "tangent-class/data/common_shift_t0.40_n3.npy"]
for f in files:
    X = np.load(f); Dd = X.shape[0]; nn = int(round(np.log(Dd) / np.log(3)))
    if nn > 4: continue
    n = nn; D = Dd
    Uu, S, Vh = np.linalg.svd(X); u0 = Uu[:, :2].T.copy()
    fobj = lambda x: m_plane(np.array([filt(x) @ v for v in u0]))
    res = minimize(fobj, np.zeros(9 * n), method='L-BFGS-B', bounds=[(-8, 8)] * (9 * n), options={'maxiter': 3000, 'gtol': 1e-11, 'ftol': 1e-14})
    uf = np.array([filt(res.x) @ v for v in u0])
    print(f"{f.split('/')[-2]}/{f.split('/')[-1]} (n={n}): m {fobj(np.zeros(9 * n)):.5f} -> {res.fun:.6f}; smallest range local-support ratio {thin_ratio(uf):.1e}", flush=True)
