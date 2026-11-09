"""S_N = conv{ pibar_C : rank C <= 2 } (identity-count distributions), support points by autograd, then convex hull.
usage: python polytope.py N ndir restarts seed   (PyTorch; writes supp_torch_N{N}.npy)"""
import sys, itertools
import numpy as np, torch
from scipy.spatial import ConvexHull
torch.set_default_dtype(torch.float64)
N, NDIR, RR, SEED = (int(a) for a in sys.argv[1:5])
torch.manual_seed(SEED); rng = np.random.default_rng(SEED)
D = 3 ** N
subsets = list(itertools.product([0, 1], repeat=N))
sizes = np.array([sum(K) for K in subsets])
sup = np.array([[all(K[i] >= L[i] for i in range(N)) for K in subsets] for L in subsets])
mob = torch.tensor(np.where(sup, (-1.0) ** (sizes[None, :] - sizes[:, None]), 0.0))
agg = torch.tensor(np.array([[1.0 if sizes[j] == k else 0.0 for j in range(len(subsets))] for k in range(N + 1)]))
scale = torch.tensor(3.0 ** -sizes)
def pibar(A, B):
    TA = A.reshape((3,) * N + (2,)); TB = B.reshape((3,) * N + (2,))
    mus = []
    for K in subsets:
        kept = [i for i in range(N) if not K[i]]; tr = [i for i in range(N) if K[i]]
        Ma = TA.permute(*kept, *tr, N).reshape(3 ** len(kept), 3 ** len(tr), 2)
        Mb = TB.permute(*kept, *tr, N).reshape(3 ** len(kept), 3 ** len(tr), 2)
        T = torch.einsum('xai,yai->xy', Ma, Mb.conj())
        mus.append((T.abs() ** 2).sum())
    mu = torch.stack(mus) * scale
    C = A @ B.conj().T
    e = mob @ mu
    return (agg @ e) / (C.abs() ** 2).sum()
def support(g, restarts):
    gt = torch.tensor(g); best = (np.inf, None)
    for r in range(restarts):
        ar = torch.randn(D, 2, 2, requires_grad=True); br = torch.randn(D, 2, 2, requires_grad=True)
        opt = torch.optim.LBFGS([ar, br], lr=1, max_iter=500, tolerance_grad=1e-11, tolerance_change=1e-14,
                                history_size=30, line_search_fn='strong_wolfe')
        def closure():
            opt.zero_grad()
            val = gt @ pibar(torch.view_as_complex(ar), torch.view_as_complex(br))
            val.backward(); return val
        for _ in range(3): opt.step(closure)
        with torch.no_grad():
            p = pibar(torch.view_as_complex(ar), torch.view_as_complex(br)).numpy()
        v = float(g @ p)
        if v < best[0]: best = (v, p)
    return best
P = []
dirs = [(-0.5) ** np.arange(N + 1)]
for _ in range(NDIR):
    v = rng.normal(size=N + 1); v -= v.mean(); dirs.append(v / np.linalg.norm(v))
for i, g in enumerate(dirs):
    val, p = support(g, RR)
    P.append(p)
    if i == 0: print(f"endpoint direction: min = {val:+.3e} at pibar x3^N = {np.round(p * 3 ** N, 3).tolist()}", flush=True)
P = np.array(P); np.save(f"supp_torch_N{N}.npy", P)
hull = ConvexHull(P[:, :-1])
V = P[hull.vertices]
uniq = []
for v in V:
    if not any(np.abs(v - u).max() < 2e-3 for u in uniq): uniq.append(v)
print(f"N={N}: {len(P)} support points, ~{len(uniq)} distinct hull vertices (x 3^N, endpoint value):")
for v in sorted(uniq, key=lambda r: tuple(np.round(r, 4))):
    print("  ", np.round(v * 3 ** N, 3).tolist(), "  endpoint", round(float(v @ (-0.5) ** np.arange(N + 1)), 5))
