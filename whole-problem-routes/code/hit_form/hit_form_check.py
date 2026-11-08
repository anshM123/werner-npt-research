"""Classical forms of the value of P_{3/4}^{(x)n} (OBSTRUCTIONS.md, "Classical forms of the value").

(a) Hit form. P_{3/4}(X) = int 3<v|X|v> [(I - vv^T) - (1/4) vv^T] dv (v uniform on the real unit sphere), so
    <y|P^{(x)n}(xx^T)|y> = E[(-1/4)^N]: v is drawn from 3^n |<v|x>|^2 dv and N counts the copies at which the test
    {v_s v_s^T, I - v_s v_s^T} on y gives v_s. On the pair (x_s, y_s) of a copy, "hit" is M = (I + F + Omega)/5 and
    "no hit" is I - M, so P(hit set = T) = <x (x) y| M^{(x)T} (x) (I - M)^{(x)T^c} |x (x) y>.
(b) Exact three-copy classical law with conditional hit probabilities in [1/5, 3/5] and a negative value; worst
    exchangeable laws for n <= 16 (the value is linear-fractional in each ratio, so vertices suffice).
(c) Spin-sector form V = (7/12)^n sum_S (-6/7)^{|S|} c_S, c_S = <Pi_2^{(x)S} sigma_S, Pi_2^{(x)S} rho_S>."""
import functools, itertools, math
from fractions import Fraction as Fr
import numpy as np

I9 = np.eye(9); F = np.zeros((9, 9)); Om = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        F[3 * i + j, 3 * j + i] = 1
        Om[3 * i + i, 3 * j + j] = 1
M = (I9 + F + Om) / 5
C = 0.75 * I9 - 0.25 * (F + Om)
assert np.allclose((I9 - M) - 0.25 * M, C)

def pair_vec(x, y, n):
    T = np.multiply.outer(x.reshape([3] * n), y.reshape([3] * n))
    return np.transpose(T, [k for s in range(n) for k in (s, n + s)]).reshape(-1)

def expect(ops, w, n):
    T = w.reshape([9] * n)
    for s, O in enumerate(ops):
        T = np.moveaxis(np.tensordot(O, T, axes=([1], [s])), 0, s)
    return float(w @ T.reshape(-1))

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

def hit_report(name, x, y, n):
    x = x / np.linalg.norm(x); y = y / np.linalg.norm(y)
    w = pair_vec(x, y, n)
    full = {T: expect([M if t else I9 - M for t in T], w, n) for T in itertools.product([0, 1], repeat=n)}
    dist = [sum(p for T, p in full.items() if sum(T) == k) for k in range(n + 1)]
    G = sum(dist[k] * (-0.25) ** k for k in range(n + 1))
    V = float(y @ Pn(np.outer(x, x), n) @ y)
    marg = [sum(p for T, p in full.items() if T[s]) for s in range(n)]
    cov = [sum(p for T, p in full.items() if T[s] and T[t]) - marg[s] * marg[t] for s in range(n) for t in range(s + 1, n)]
    assert abs(V - G) < 1e-12, (V, G)
    print(f'{name} (n={n}): value {V:.8f} = {V * 4 ** n:.4f} x 4^-n = E[(-1/4)^N] (diff {abs(V - G):.1e}); '
          f'P(N=k) {np.round(dist, 5)}; hit probabilities {np.round(marg, 4)}; pair covariances '
          f'{min(cov):+.4f} to {max(cov):+.4f}')

e = np.eye(3)
kron = lambda *vs: functools.reduce(np.kron, vs)
print('(a) hit form')
hit_report('aligned real product', kron(e[0], e[0]), kron(e[0], e[0]), 2)
hit_report('orthogonal real product', kron(e[0], e[0]), kron(e[1], e[1]), 2)
phi = sum(kron(e[i], e[i]) for i in range(3))
hit_report('maximally entangled pair, y = x', phi, phi, 2)
eps = np.zeros(27)
for p in itertools.permutations(range(3)):
    eps[9 * p[0] + 3 * p[1] + p[2]] = round(np.linalg.det(np.eye(3)[list(p)]))
ev, V = np.linalg.eigh(Pn(np.outer(eps, eps) / 6, 3))
hit_report('antisymmetric tensor, bottom test vector', eps, V[:, 0], 3)
n = 4
def omega_pairs(pairs):
    v = np.zeros([3] * n)
    for idx in itertools.product(range(3), repeat=n):
        if all(idx[s] == idx[t] for s, t in pairs):
            v[idx] = 1
    return v.reshape(-1)
p1, p2, p3 = (omega_pairs(pr) for pr in ([(0, 1), (2, 3)], [(0, 2), (1, 3)], [(0, 3), (1, 2)]))
c = (math.sqrt(6805) - 71) / 36
hit_report('singlet pair of Observation 7.26', p3 - p1, p2 - c * (p1 + p3), 4)
rng = np.random.default_rng(1)
for k in range(3):
    hit_report(f'random real pair {k}', rng.standard_normal(27), rng.standard_normal(27), 3)

print('(b) classical relaxation (exact rationals)')
def law(R):
    u = [Fr(1)]
    for r in R:
        u.append(u[-1] * r)
    return u
def value(u, n):
    Z = sum(math.comb(n, k) * u[k] for k in range(n + 1))
    return sum(math.comb(n, k) * u[k] * Fr(-1, 4) ** k for k in range(n + 1)) / Z
R = [Fr(3, 2), Fr(1, 4), Fr(3, 2)]
u = law(R)
val3 = value(u, 3)
print(f'three-copy law, conditional hit probabilities {[str(r / (1 + r)) for r in R]} given 0, 1, 2 hits elsewhere: '
      f'E[(-1/4)^N] = {val3}')
assert val3 == Fr(-13, 1472)
for j in range(2):
    Zb = sum(math.comb(2, i) * u[i + j] for i in range(3))
    vb = sum(math.comb(2, i) * u[i + j] * Fr(-1, 4) ** i for i in range(3)) / Zb
    print(f'  two-copy block given {j} hit(s) at the third copy: {vb}')
    assert vb > 0
for n in range(2, 17):
    best = min((value(law(R), n), R) for R in itertools.product([Fr(1, 4), Fr(3, 2)], repeat=n))
    print(f'  n={n}: worst exchangeable value {float(best[0]):+.6f} '
          f'({"".join("H" if r == Fr(3, 2) else "L" for r in best[1])})')

print('(c) spin-sector form')
def marg(v, n, S):
    T = v.reshape([3] * n)
    rest = [t for t in range(n) if t not in S]
    k = len(S)
    return np.tensordot(T, T, axes=(rest, rest)).reshape(3 ** k, 3 ** k)
def pi2(Mx, k):
    T = Mx.reshape([3] * (2 * k))
    for s in range(k):
        Tm = np.moveaxis(T, [s, k + s], [0, 1])
        tr = np.einsum('ii...->...', Tm)
        out = 0.5 * (Tm + np.swapaxes(Tm, 0, 1))
        for i in range(3):
            out[i, i] -= tr / 3
        T = np.moveaxis(out, [0, 1], [s, k + s])
    return T.reshape(3 ** k, 3 ** k)
def spin_form(x, y, n):
    x = x / np.linalg.norm(x); y = y / np.linalg.norm(y)
    tot, by_size = 0.0, [0.0] * (n + 1)
    for k in range(n + 1):
        for S in itertools.combinations(range(n), k):
            cS = 1.0 if k == 0 else float(np.sum(pi2(marg(x, n, list(S)), k) * pi2(marg(y, n, list(S)), k)))
            by_size[k] += (-6 / 7) ** k * cS
    V = float(y @ Pn(np.outer(x, x), n) @ y)
    return V, (7 / 12) ** n * sum(by_size), by_size
for name, x, y in [('singlet pair', p3 - p1, p2 - c * (p1 + p3)), ('random real pair', rng.standard_normal(81), rng.standard_normal(81))]:
    V, W, by = spin_form(x, y, 4)
    assert abs(V - W) < 1e-12
    print(f'{name} (n=4): value {V:.8f}, spin-sector form {W:.8f}; weighted sums by |S| {np.round(by, 5)}')
ame = np.zeros([3] * 4)
for i in range(3):
    for j in range(3):
        ame[i, j, (i + j) % 3, (i + 2 * j) % 3] = 1
ame = ame.reshape(-1)
ev, Vv = np.linalg.eigh(Pn(np.outer(ame, ame) / 9, 4))
V, W, by = spin_form(ame, Vv[:, 0], 4)
assert abs(V - W) < 1e-12
print(f'AME(4,3) input, bottom test vector: value {V:.8f} = {V * 256:.4f} x 4^-4; weighted sums by |S| {np.round(by, 5)} '
      f'(no terms with |S| <= 2: the input is 2-uniform)')
print('PASS')
