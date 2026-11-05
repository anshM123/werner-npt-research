"""Exact arithmetic (Fractions) for Proposition Q of NOTES.md: the P- and Q-states of the type I block T = {t} of
(R_inf)_1, written as the letter-1 block of Theorem B0 with an extra copy (the "phantom" copy 0) that enters only through
the constraint and the norm.

Phantom weights.  For |s| = 1, sigma = x_t, put y_c = 2 x_c/(2 - sigma) for c in U and y_0 = sigma/(2 - sigma), so that
y_0 + sum_U y_c = 1.  Put lam_c = 1/(1 - y_c), g_c = y_c/(2(1 - y_c)) = (lam_c - 1)/2, gp = sum_{c in U} g_c and
del_cu = 1 - y_c - y_u (c, u in U).

For zeta indexed by {0} + U with sum_all y_c zeta_c = 0:
  S(zeta)  = sum_all y_c zeta_c^2,
  L(zeta)  = sum_{c in U} y_c lam_c^2 zeta_c^2 + (1/4) sum_{c<u in U} (y_c y_u/del_cu) (lam_c zeta_c + lam_u zeta_u)^2,
and the identity of Proposition Q (proof in NOTES.md, Section 2) reads (2 + 2 gp) S - L = R1 + R2 + R3 + R4 with
  R1 = sum_{r<t in U} 2 y_r y_t (gp - g_r - g_t) d_rt^2 + sum_{r in U} y_r y_0 (1 + 2 gp - 2 g_r) d_r0^2,
  R2 = (1/2) sum_{l in U} y_l lam_l^2 sum_{a in U - l} y_a y_0 d_a0^2,
  R3 = (1/4) sum_{r<t in U} sum_{u not in {r,t}} y_r y_t y_u (lam_t d_ru - lam_r d_tu)^2,
  R4 = (1/4) sum_{r<t in U} (y_r y_t/del_rt) sum_{u<v, u,v not in {r,t}} y_u y_v (lam_r + lam_t)^2 d_uv^2,
where d_ab = zeta_a - zeta_b and u, v range over all copies (the phantom 0 included)."""
from fractions import Fraction as Fr
import itertools


def phantom(x, t):
    """x: squared weights with sum 1 (Fractions), t: the copy with the letter 2.  Returns (U, y) with y a dict on
    {0} + U, key 'ph' for the phantom copy."""
    k = len(x)
    sig = x[t]
    U = [c for c in range(k) if c != t]
    y = {c: 2 * x[c] / (2 - sig) for c in U}
    y['ph'] = sig / (2 - sig)
    assert sum(y.values()) == 1
    return U, y


def forms(U, y, zeta):
    """S, L and the four remainder terms at a vector zeta (dict on 'ph' + U) satisfying the constraint."""
    allc = ['ph'] + list(U)
    assert sum(y[c] * zeta[c] for c in allc) == 0
    lam = {c: 1 / (1 - y[c]) for c in U}
    g = {c: y[c] / (2 * (1 - y[c])) for c in U}
    gp = sum(g.values())
    d = lambda a, b: zeta[a] - zeta[b]
    S = sum(y[c] * zeta[c] ** 2 for c in allc)
    L = sum(y[c] * lam[c] ** 2 * zeta[c] ** 2 for c in U)
    L += Fr(1, 4) * sum(y[c] * y[u] / (1 - y[c] - y[u]) * (lam[c] * zeta[c] + lam[u] * zeta[u]) ** 2
                        for c, u in itertools.combinations(U, 2))
    R1 = sum(2 * y[r] * y[s] * (gp - g[r] - g[s]) * d(r, s) ** 2 for r, s in itertools.combinations(U, 2))
    R1 += sum(y[r] * y['ph'] * (1 + 2 * gp - 2 * g[r]) * d(r, 'ph') ** 2 for r in U)
    R2 = Fr(1, 2) * sum(y[l] * lam[l] ** 2 * sum(y[a] * y['ph'] * d(a, 'ph') ** 2 for a in U if a != l) for l in U)
    R3 = Fr(0)
    R4 = Fr(0)
    for r, s in itertools.combinations(U, 2):
        others = [u for u in allc if u not in (r, s)]
        R3 += Fr(1, 4) * sum(y[r] * y[s] * y[u] * (lam[s] * d(r, u) - lam[r] * d(s, u)) ** 2 for u in others)
        R4 += Fr(1, 4) * (y[r] * y[s] / (1 - y[r] - y[s])) * sum(
            y[u] * y[v] * (lam[r] + lam[s]) ** 2 * d(u, v) ** 2 for u, v in itertools.combinations(others, 2))
    coef_ok = all(gp - g[r] - g[s] >= 0 for r, s in itertools.combinations(U, 2)) and \
        all(1 + 2 * gp - 2 * g[r] > 0 for r in U) and all(1 - y[r] - y[s] > 0 for r, s in itertools.combinations(U, 2))
    return dict(S=S, L=L, gp=gp, R=(R1, R2, R3, R4), coef_ok=coef_ok)
