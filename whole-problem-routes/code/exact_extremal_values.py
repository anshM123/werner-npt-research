"""Checks for Observation 7.27(d), Proposition 7.33 and the remark after it.

(a) eps/sqrt6 at three copies: 192 * P^{(x)3}(eps eps^T / 6) has characteristic polynomial (x-21)^3 (x-36)^7 (x-42)^17
    (exact, python-flint), so the spectrum is 7/64 (x3), 12/64 (x7), 14/64 (x17).
(b) the input y = (p (x) e + e (x) p)/sqrt2 on two pairs of copies (p = e0 e0, e = (e1 e1 + e2 e2)/sqrt2): the least output
    eigenvalue is the root (63 - sqrt 2029)/512 of x^2 - 126 x + 1940 (for 512 W, exact).
(c) the balanced input y = (p e e + e p e + e e p)/sqrt3 on three pairs: least output eigenvalue 9/1024 (floating point).
(d) Proposition 7.33 at three copies on random separable inputs (floating point): the ratio to 4^-n 3^{2 sum tau} is >= 1.
(e) the singlet pair of Observation 7.26: lambda_min(P^{(x)4}(x x^T)) = (151 - sqrt 6805)/1536 with a balanced bottom
    eigenvector, below the floor (3/4)^4 3^{-4/3} that every separable input gives at a balanced test vector."""
import itertools
import math
from fractions import Fraction
import numpy as np
import flint


def W_of(y, n, exact=False):
    """P_{3/4}^{(x)n}(y y^T), applied copy by copy"""
    T = np.multiply.outer(y, y).reshape([3] * (2 * n))
    a = Fraction(3, 4) if exact else 0.75
    for s in range(n):
        Tm = np.moveaxis(T, [s, n + s], [0, 1])
        tr = sum(Tm[i, i] for i in range(3))
        new = np.empty_like(Tm)
        for i in range(3):
            for j in range(3):
                val = -(Tm[i, j] + Tm[j, i]) / 4
                new[i, j] = val + a * tr if i == j else val
        T = np.moveaxis(new, [0, 1], [s, n + s])
    return T.reshape(3 ** n, 3 ** n)


def charpoly_factors(W):
    den = 1
    for row in W:
        for v in row:
            den = den * v.denominator // math.gcd(den, v.denominator)
    cp = flint.fmpz_mat([[int(v * den) for v in row] for row in W]).charpoly()
    return den, sorted((str(f), e) for f, e in cp.factor()[1])


# (a)
eps = np.empty(27, dtype=object)
eps[:] = [Fraction(0)] * 27
for (i, j, k), s in [((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)]:
    eps[9 * i + 3 * j + k] = Fraction(s)
den, fac = charpoly_factors(W_of(eps, 3, exact=True) * Fraction(1, 6))
print('(a) eps/sqrt6:', den, '* W has charpoly factors', fac)
assert den == 192 and fac == sorted([('x + (-21)', 3), ('x + (-36)', 7), ('x + (-42)', 17)])

# (b)
yq = np.empty(81, dtype=object)
yq[:] = [Fraction(0)] * 81
for a_ in (4, 8):                         # e = (e1 e1 + e2 e2)/sqrt2 has pair indices 4 and 8
    yq[a_] += Fraction(1, 2)              # p (x) e  (pair 1 = e0 e0 = index 0)
    yq[9 * a_] += Fraction(1, 2)          # e (x) p
den, fac = charpoly_factors(W_of(yq, 4, exact=True))
print('(b) Dicke-over-pairs (2,1):', den, '* W has the factor x^2 - 126x + 1940:', ('x^2 + (-126)*x + 1940', 1) in fac)
assert den == 512 and ('x^2 + (-126)*x + 1940', 1) in fac
lam_b = np.linalg.eigvalsh(np.array(W_of(yq, 4, exact=True), dtype=float))[0]
assert abs(lam_b - (63 - math.sqrt(2029)) / 512) < 1e-14
print(f'    least eigenvalue {lam_b:.15f} = (63 - sqrt 2029)/512')

# (c)
p = np.zeros(9); p[0] = 1
e = np.zeros(9); e[4] = e[8] = 1 / math.sqrt(2)
y = (np.kron(np.kron(p, e), e) + np.kron(np.kron(e, p), e) + np.kron(np.kron(e, e), p)) / math.sqrt(3)
ev = np.linalg.eigvalsh(W_of(y, 6))
print(f'(c) Dicke-over-pairs (3,1): least eigenvalue {ev[0]:.15f}, 9/1024 = {9 / 1024:.15f}')
assert abs(ev[0] - 9 / 1024) < 1e-13
for s in range(6):
    T = y.reshape([3] * 6)
    ax = [t for t in range(6) if t != s]
    assert np.allclose(np.tensordot(T, T, axes=(ax, ax)), np.eye(3) / 3)

# (d)
rng = np.random.default_rng(3)
worst = 9.0
for trial in range(300):
    K = int(rng.integers(1, 6))
    w = rng.random(K)
    w /= w.sum()
    S = np.zeros((27, 27))
    for k in range(K):
        v = np.ones(1)
        for s in range(3):
            u = rng.standard_normal(3)
            v = np.kron(v, u / np.linalg.norm(u))
        S += w[k] * np.outer(v, v)
    evS, U = np.linalg.eigh(S)
    WS = sum(evS[i] * W_of(U[:, i], 3) for i in range(27) if evS[i] > 1e-14)
    T6 = S.reshape([3] * 6)
    taus = [np.linalg.eigvalsh(np.einsum({0: 'aijbij->ab', 1: 'iajibj->ab', 2: 'ijaijb->ab'}[s], T6))[0] for s in range(3)]
    worst = min(worst, np.linalg.eigvalsh(WS)[0] / (4.0 ** -3 * 3 ** (2 * sum(taus))))
print(f'(d) Proposition 7.33 at n = 3, 300 random separable inputs: least ratio {worst:.6f}')
assert worst > 1 - 1e-9

# (e)
def omega_pair(n, pairs):
    v = np.zeros([3] * n)
    for idx in itertools.product(range(3), repeat=n):
        if all(idx[s] == idx[t] for s, t in pairs):
            v[idx] = 1
    return v.reshape(-1)


p1, p2, p3 = (omega_pair(4, pr) for pr in ([(0, 1), (2, 3)], [(0, 2), (1, 3)], [(0, 3), (1, 2)]))
x = (p3 - p1) / np.linalg.norm(p3 - p1)
c = (math.sqrt(6805) - 71) / 36
yb = p2 - c * (p1 + p3)
yb /= np.linalg.norm(yb)
ex, vx = np.linalg.eigh(W_of(x, 4))
floor = 0.75 ** 4 * 3 ** (-4 / 3)
print(f'(e) singlet pair: lambda_min(W_x) = {ex[0]:.15f}, (151 - sqrt 6805)/1536 = {(151 - math.sqrt(6805)) / 1536:.15f}, '
      f'overlap with y {abs(vx[:, 0] @ yb):.12f}; separable floor {floor:.6f}')
assert abs(ex[0] - (151 - math.sqrt(6805)) / 1536) < 1e-14 and abs(abs(vx[:, 0] @ yb) - 1) < 1e-9 and floor > ex[0]
print('PASS')
