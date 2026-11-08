"""Compressing two copies at once (addendum to Observation 7.25, 8 November 2026).

(a) n = 4, the partner vector y = p2 - c(p1 + p3) of Observation 7.26: every one-copy compression raises g, but
    compressing a pair of copies, each onto a plane, lowers it (best plane normals found by Nelder-Mead).
(b) n = 6, the pair minimising the value over spin singlets (alternating minimisation in the 15-dimensional singlet
    sector): g at its test vector over all inputs, the one-copy compression gap (one direction suffices, by rotation
    invariance) and the best two-copy compression gap for the pair of copies (0, 3).
Floating point."""
import itertools, math
import numpy as np
from scipy.optimize import minimize

def Pn(X, n):
    T = X.reshape((3,) * (2 * n))
    for s in range(n):
        Tm = np.moveaxis(T, [s, n + s], [0, 1])
        tr = np.einsum('ii...->...', Tm)
        out = -0.25 * (Tm + np.swapaxes(Tm, 0, 1))
        for i in range(3):
            out[i, i] += 0.75 * tr
        T = np.moveaxis(out, [0, 1], [s, n + s])
    return T.reshape(3 ** n, 3 ** n)
def g(y, n):
    return np.linalg.eigvalsh(Pn(np.outer(y, y), n))[0] / (y @ y)
def unit(th):
    return np.array([np.sin(th[0]) * np.cos(th[1]), np.sin(th[0]) * np.sin(th[1]), np.cos(th[0])])
def compress(v, n, copies, ws):
    T = v.reshape([3] * n)
    for s, w in zip(copies, ws):
        T = np.moveaxis(np.tensordot(np.eye(3) - np.outer(w, w), T, axes=([1], [s])), 0, s)
    return T.reshape(-1)
def matchings(points):
    if not points:
        yield []
        return
    p = points[0]
    for i in range(1, len(points)):
        rest = points[1:i] + points[i + 1:]
        for mm in matchings(rest):
            yield [(p, points[i])] + mm
def omega_pairs(pairs, n):
    T = np.zeros([3] * n)
    for idx in itertools.product(range(3), repeat=n // 2):
        full = [0] * n
        for (s, t), v in zip(pairs, idx):
            full[s] = v; full[t] = v
        T[tuple(full)] = 1
    return T.reshape(-1)
def best_pair(v, n, pair, rng, starts=6, maxfev=2000):
    def f(th):
        z = compress(v, n, list(pair), [unit(th[:2]), unit(th[2:])])
        return g(z, n) if z @ z > 1e-12 else 1.0
    return min(minimize(f, rng.uniform(0, np.pi, 4) * [1, 2, 1, 2], method='Nelder-Mead',
                        options={'xatol': 1e-8, 'fatol': 1e-13, 'maxfev': maxfev}).fun for _ in range(starts))
rng = np.random.default_rng(0)

print('(a) four copies, partner vector')
n = 4
p1, p2, p3 = (omega_pairs(pr, n) for pr in ([(0, 1), (2, 3)], [(0, 2), (1, 3)], [(0, 3), (1, 2)]))
c = (math.sqrt(6805) - 71) / 36
y = p2 - c * (p1 + p3); y /= np.linalg.norm(y); g0 = g(y, n)
one = g(compress(y, n, [0], [np.array([0.0, 0.0, 1.0])]), n) - g0
print(f'g(y) = {g0 * 256:.6f} x 4^-4; one-copy compression gap {one * 256:+.6f}')
for pair in itertools.combinations(range(n), 2):
    print(f'  pair {pair}: best two-copy compression gap {(best_pair(y, n, pair, rng) - g0) * 256:+.5f} x 4^-4')

print('(b) six copies, minimising singlet pair')
n = 6
SB = np.array([omega_pairs(m, n) for m in matchings(list(range(n)))]).T
U, sv, _ = np.linalg.svd(SB, full_matrices=False); Qs = U[:, sv > 1e-8 * sv[0]]
best = (np.inf, None, None)
for st in range(4):
    yy = Qs @ rng.standard_normal(Qs.shape[1]); yy /= np.linalg.norm(yy); lam_old = 9
    for it in range(300):
        ev, V = np.linalg.eigh(Qs.T @ Pn(np.outer(yy, yy), n) @ Qs); xx = Qs @ V[:, 0]
        ev, V = np.linalg.eigh(Qs.T @ Pn(np.outer(xx, xx), n) @ Qs); yy = Qs @ V[:, 0]; lam = ev[0]
        if abs(lam - lam_old) < 1e-14:
            break
        lam_old = lam
    if lam < best[0]:
        best = (lam, xx, yy)
lam, xx, yy = best
gy = g(yy, n)
print(f'singlet-sector minimum {lam * 4 ** n:.4f} x 4^-6; g(test vector) over all inputs {gy * 4 ** n:.4f} x 4^-6')
one = g(compress(yy, n, [0], [np.array([0.0, 0.0, 1.0])]), n) - gy
print(f'one-copy compression gap {one * 4 ** n:+.4f} x 4^-6')
print(f'pair (0, 3): best two-copy compression gap {(best_pair(yy, n, (0, 3), rng, starts=3, maxfev=800) - gy) * 4 ** n:+.4f} x 4^-6')
