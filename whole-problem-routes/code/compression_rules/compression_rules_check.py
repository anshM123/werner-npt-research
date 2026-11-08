"""One-copy compressions of the test vector (Observation 7.25, addendum of 8 November 2026).

g(y) = lambda_min(P_{3/4}^{(x)n}(yy^T)) / |y|^2, and G_{s,w} = I - ww^T at copy s (the input is re-optimised, since g
already minimises over it).
(a) n = 4, y = p2 - c(p1 + p3), the partner vector of Observation 7.26. Every one-copy compression raises g. By the
    global O(3) invariance of y and the copy permutations that fix it, all copies and directions give one value;
    the script samples 200 directions at every copy (floating point) and recomputes one value in 30-digit arithmetic.
(b) n = 3: inputs at which compressing along a fixed eigenvector of the one-copy marginal raises g at every copy
    (smallest eigenvector; middle eigenvector), with the best of the three eigenvectors and the best direction there."""
import itertools, math
import numpy as np
import mpmath as mp

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
def compress(y, n, s, w):
    T = y.reshape([3] * n); Pw = np.eye(3) - np.outer(w, w)
    return np.moveaxis(np.tensordot(Pw, T, axes=([1], [s])), 0, s).reshape(-1)
def marg(y, n, s):
    T = y.reshape([3] * n); rest = [t for t in range(n) if t != s]
    return np.tensordot(T, T, axes=(rest, rest))

print('(a) four-copy partner vector')
n = 4
def omega_pairs(pairs):
    v = np.zeros([3] * n)
    for idx in itertools.product(range(3), repeat=n):
        if all(idx[s] == idx[t] for s, t in pairs):
            v[idx] = 1
    return v.reshape(-1)
p1, p2, p3 = (omega_pairs(pr) for pr in ([(0, 1), (2, 3)], [(0, 2), (1, 3)], [(0, 3), (1, 2)]))
c = (math.sqrt(6805) - 71) / 36
y = p2 - c * (p1 + p3)
g0 = g(y, n)
print(f'g(y) = {g0 * 256:.12f} x 4^-4, (151 - sqrt 6805)/6 = {(151 - math.sqrt(6805)) / 6:.12f}')
rng = np.random.default_rng(0)
gaps = []
for s in range(n):
    for _ in range(200):
        w = rng.standard_normal(3); w /= np.linalg.norm(w)
        gaps.append(g(compress(y, n, s, w), n) - g0)
gaps = np.array(gaps) * 256
print(f'compressions at 4 copies x 200 directions: gap in [{gaps.min():.10f}, {gaps.max():.10f}] x 4^-4 (all > 0)')
assert gaps.min() > 0.3
mp.mp.dps = 30
z = compress(y, n, 0, np.array([0.0, 0.0, 1.0]))
# recompute both values with 30-digit eigenvalues of the same (double precision) matrices
Wz = mp.matrix((Pn(np.outer(z, z), n) / (z @ z)).tolist())
Wy = mp.matrix((Pn(np.outer(y, y), n) / (y @ y)).tolist())
gz_mp, gy_mp = min(mp.eigsy(Wz)[0]), min(mp.eigsy(Wy)[0])
print(f'30-digit eigenvalues: g(y) = {mp.nstr(gy_mp * 256, 15)}, g(Gy) = {mp.nstr(gz_mp * 256, 15)} (x 4^-4), gap {mp.nstr((gz_mp - gy_mp) * 256, 10)}')

print('(b) three-copy failures of fixed eigenvector rules')
n = 3
for name, f, k in (('smallest eigenvector', 'counterexample_bottom_n3.npy', 0), ('middle eigenvector', 'counterexample_middle_n3.npy', 1)):
    y = np.load(f); y /= np.linalg.norm(y); g0 = g(y, n)
    rule = []; any3 = []
    for s in range(n):
        ev, V = np.linalg.eigh(marg(y, n, s))
        rule.append(g(compress(y, n, s, V[:, k]), n) - g0)
        any3.append(min(g(compress(y, n, s, V[:, j]), n) for j in range(3)) - g0)
    best = []
    for s in range(n):
        vals = []
        for _ in range(2000):
            w = rng.standard_normal(3); w /= np.linalg.norm(w)
            vals.append(g(compress(y, n, s, w), n))
        best.append(min(vals) - g0)
    print(f'{name} rule: g(y) = {g0 * 64:.6f} x 4^-3; rule gaps {np.round(np.array(rule) * 64, 6)}; best of the three '
          f'eigenvectors {np.round(np.array(any3) * 64, 4)}; best sampled direction {np.round(np.array(best) * 64, 4)} (x 4^-3)')
    assert min(rule) > 0 and min(any3) < 0
print('PASS')
