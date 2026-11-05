"""Checks for the algebraic proof of Theorem B0 (README, Section 3b).
(a) Exact (Fractions) at random rational points, k = 3..7, real zeta with sum_r x_r zeta_r = 0:
    E - D = G, T <= G, and the sum-of-squares identity (2 + 2 gamma) S - L = line1 + line2 + line3.
    (Both sides are Hermitian forms with real coefficients, so the identity for real zeta gives it for complex zeta.)
(b) The letter-2 termwise inequality in SymPy: (1+e)y/(2e(1+e-2y)) - y/(1-y) - y/(2e) = -y(e-y)(1+e-y)/(e(1-y)(1+e-2y)).
(c) Floating point, complex zeta: L <= (2 + 2 gamma) S and (2 + 2 gamma) <= (12/25)(2 + gamma)^2 at random points."""
import itertools
import numpy as np
import sympy as sp
from fractions import Fraction as Fr


def parts(x, z, absq):
    k = len(x); lam = [1 / (1 - xi) for xi in x]; gam = sum(xi * li for xi, li in zip(x, lam)) / 2
    pairs = list(itertools.combinations(range(k), 2))
    S = sum(x[r] * absq(z[r]) for r in range(k))
    D = sum(x[r] * lam[r] ** 2 * absq(z[r]) for r in range(k))
    T = sum(x[r] * x[t] / (1 - x[r] - x[t]) * absq(lam[r] * z[r] + lam[t] * z[t]) for r, t in pairs) / 4
    E = sum(x[r] * x[t] * (lam[r] + lam[t]) * absq(z[r] - z[t]) for r, t in pairs)
    G = sum(x[r] * x[t] * sum(x[u] * lam[u] ** 2 for u in range(k) if u not in (r, t)) * absq(z[r] - z[t]) for r, t in pairs)
    l1 = sum(x[r] * x[t] * x[u] * lam[u] * absq(z[r] - z[t]) for r, t in pairs for u in range(k) if u not in (r, t))
    l2 = sum(x[r] * x[t] * x[u] * absq(lam[t] * (z[r] - z[u]) - lam[r] * (z[t] - z[u])) for r, t in pairs for u in range(k) if u not in (r, t)) / 4
    l3 = sum(x[r] * x[t] / (1 - x[r] - x[t]) * sum(x[u] * x[v] * (lam[r] + lam[t]) ** 2 * absq(z[u] - z[v])
             for u, v in itertools.combinations([q for q in range(k) if q not in (r, t)], 2)) for r, t in pairs) / 4
    return S, D, T, E, G, gam, l1, l2, l3


ok = True
rng = np.random.default_rng(1)
for trial in range(40):
    k = int(rng.integers(3, 8)); xi = [Fr(int(rng.integers(1, 40))) for _ in range(k)]; tot = sum(xi); x = [v / tot for v in xi]
    z = [Fr(int(rng.integers(-20, 21))) for _ in range(k - 1)]; z.append(-sum(x[r] * z[r] for r in range(k - 1)) / x[k - 1])
    S, D, T, E, G, gam, l1, l2, l3 = parts(x, z, lambda v: v * v)
    ok &= (E - D == G) and (T <= G) and ((2 + 2 * gam) * S - (D + T) == l1 + l2 + l3)
print(f"(a) exact, 40 rational points (k = 3..7): E - D = G, T <= G, sum-of-squares identity: {'PASS' if ok else 'FAIL'}")
e, y = sp.symbols('e y', positive=True)
lhs = (1 + e) * y / (2 * e * (1 + e - 2 * y)) - y / (1 - y) - y / (2 * e)
rhs = -y * (e - y) * (1 + e - y) / (e * (1 - y) * (1 + e - 2 * y))
okb = sp.simplify(lhs - rhs) == 0
a, Gp = sp.symbols('a Gp', positive=True)
okb &= sp.expand((a + 1 + Gp) ** 2 - a ** 2 - a * (Gp + sp.Rational(1, 4)) - (a * (sp.Rational(7, 4) + Gp) + (1 + Gp) ** 2)) == 0
print(f"(b) letter-2 termwise identity and the final determinant bound: {'PASS' if okb else 'FAIL'}")
ok &= okb
worst = 0.0; worst2 = 0.0
for trial in range(3000):
    k = int(rng.integers(3, 9)); x = rng.dirichlet(np.ones(k) * rng.choice([0.2, 1, 5])); x = np.maximum(x, 1e-6); x /= x.sum()
    z = rng.normal(size=k) + 1j * rng.normal(size=k); z -= x @ z
    S, D, T, E, G, gam, l1, l2, l3 = parts(list(x), list(z), lambda v: abs(v) ** 2)
    worst = max(worst, (D + T) / ((2 + 2 * gam) * S)); worst2 = max(worst2, (2 + 2 * gam) / (2 + gam) ** 2 / (12 / 25))
print(f"(c) complex zeta, 3000 points: max L/((2 + 2 gamma) S) = {worst:.6f}; max (2 + 2 gamma)/((12/25)(2 + gamma)^2) = {worst2:.6f}")
ok &= worst <= 1 + 1e-9 and worst2 <= 1 + 1e-12
print("PASS" if ok else "FAIL")
