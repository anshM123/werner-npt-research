"""Six copies, O(3)-singlet test vectors y. By rotation invariance, compressing copy s onto any plane gives one value
g(G_s y) for each copy s. Maximise Gap1(y) = min_s g(G_s y) - g(y) over the 15-dimensional singlet sector. Where Gap1 > 0
(every one-copy compression raises g), test two-copy compressions."""
import itertools, math, time, sys
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import eigh
n = 6; D = 3 ** n
def Pn(X):
    T = X.reshape((3,) * (2 * n))
    for s in range(n):
        Tm = np.moveaxis(T, [s, n + s], [0, 1])
        tr = np.einsum('ii...->...', Tm)
        out = -0.25 * (Tm + np.swapaxes(Tm, 0, 1))
        for i in range(3):
            out[i, i] += 0.75 * tr
        T = np.moveaxis(out, [0, 1], [s, n + s])
    return T.reshape(D, D)
def g(y):
    return eigh(Pn(np.outer(y, y)), eigvals_only=True, subset_by_index=[0, 0])[0] / (y @ y)
def matchings(points):
    if not points:
        yield []
        return
    p = points[0]
    for i in range(1, len(points)):
        rest = points[1:i] + points[i + 1:]
        for mm in matchings(rest):
            yield [(p, points[i])] + mm
def omega_pairs(pairs):
    T = np.zeros([3] * n)
    for idx in itertools.product(range(3), repeat=n // 2):
        full = [0] * n
        for (s, t), v in zip(pairs, idx):
            full[s] = v; full[t] = v
        T[tuple(full)] = 1
    return T.reshape(-1)
SB = np.array([omega_pairs(m) for m in matchings(list(range(n)))]).T
U, sv, _ = np.linalg.svd(SB, full_matrices=False); Qs = U[:, sv > 1e-8 * sv[0]]
k = Qs.shape[1]
def compress(v, copies, ws):
    T = v.reshape([3] * n)
    for s, w in zip(copies, ws):
        T = np.moveaxis(np.tensordot(np.eye(3) - np.outer(w, w), T, axes=([1], [s])), 0, s)
    return T.reshape(-1)
e3 = np.array([0.0, 0.0, 1.0])
def gap1(c):
    y = Qs @ c; y /= np.linalg.norm(y); g0 = g(y)
    return min(g(compress(y, [s], [e3])) for s in range(n)) - g0, g0
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
t0 = time.time(); found = []
for st in range(int(sys.argv[2]) if len(sys.argv) > 2 else 8):
    c0 = rng.standard_normal(k)
    r = minimize(lambda c: -gap1(c)[0], c0, method='Nelder-Mead', options={'maxfev': 500, 'xatol': 1e-8, 'fatol': 1e-12})
    gp, g0 = gap1(r.x)
    print(f'start {st}: max one-copy gap {gp * 4 ** n:+.4f} x 4^-6 at g = {g0 * 4 ** n:.3f} x 4^-6 [{time.time() - t0:.0f}s]', flush=True)
    if gp > 0:
        found.append((gp, r.x))
found.sort(key=lambda t: -t[0])
pass  # points not kept
print(f'{len(found)} singlet test vectors with every one-copy compression raising g', flush=True)
def unit(th):
    return np.array([np.sin(th[0]) * np.cos(th[1]), np.sin(th[0]) * np.sin(th[1]), np.cos(th[0])])
for gp, c in found[:2]:
    y = Qs @ c; y /= np.linalg.norm(y); g0 = g(y)
    best = np.inf
    for pair in [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]:
        def f(th):
            z = compress(y, list(pair), [e3, unit(th)])
            return g(z) if z @ z > 1e-12 else 1.0
        val = min(minimize(f, rng.uniform(0, np.pi, 2) * [1, 2], method='Nelder-Mead', options={'xatol': 1e-8, 'fatol': 1e-13, 'maxfev': 300}).fun for _ in range(3))
        best = min(best, val)
        print(f'  point with one-copy gap {gp * 4 ** n:+.4f}: pair {pair} best two-copy gap {(val - g0) * 4 ** n:+.4f} x 4^-6 [{time.time() - t0:.0f}s]', flush=True)
