"""The reduced soft form Phi1 (NOTES.md, Sections 2-3), exact or float.
Unnormalised version (weights w on U, sigma = s_T): Phi1(a, c) = f(a, c) - (sigma/4) [ (12/5) sum_{r!=h} Gm_rh^2
+ sum_h Gm_hh^2 / (3/2 - rho_h) ], with f the soft form (F2) and Gm the coupling matrix (F3) (both as polynomials).
Normalised block formulas (K = 1): Ma, Mc, Mca with entries as in NOTES.md Section 3 (entry formulas), products of
normalised weights passed as rationals (w_h w_k / K)."""
from fractions import Fraction as Fr
from cu_form import inv_mat

BETA = Fr(12, 5)


def phi1_unnorm(w, sig, a, c, F=Fr):
    u = len(w)
    s = [x * x for x in w]
    S = sum(s, F(0)); K = sig + S / 2
    d = [s[h] + S / 2 + sig for h in range(u)]
    E = [[(d[i] if i == j else F(0)) + w[i] * w[j] / 2 for j in range(u)] for i in range(u)]
    Ei = inv_mat(E, F)
    E0Ei = [[(F(1) if i == j else F(0)) - sig * Ei[i][j] for j in range(u)] for i in range(u)]
    z = [sum(Ei[i][j] * w[j] for j in range(u)) for i in range(u)]
    p = sum(ci * wi for ci, wi in zip(c, w))
    f = (sum(a[i] * E0Ei[i][j] * a[j] for i in range(u) for j in range(u))
         + sum(c[i] * E0Ei[i][j] * c[j] for i in range(u) for j in range(u))
         - sig * p * p / (2 * K * K) + p * sum(a_i * z_i for a_i, z_i in zip(a, z))
         + 2 * sum(a[r] * c[r] * w[r] * z[r] for r in range(u)))
    G = [[2 * Ei[i][j] * a[j] * w[j] + 2 * c[i] * w[i] * Ei[i][j] + p * (Ei[i][j] - (F(1) / K if i == j else F(0)))
          for j in range(u)] for i in range(u)]
    rho = [s[h] / d[h] for h in range(u)]
    cpl = (BETA * sum(G[i][j] ** 2 for i in range(u) for j in range(u) if i != j)
           + sum(G[h][h] ** 2 / (Fr(3, 2) - rho[h]) for h in range(u)))
    return f - sig / 4 * cpl


def coupling_matrix(w, sig, a, c, F=Fr):
    """the coupling matrix Gm of (F3) (unnormalised)."""
    u = len(w)
    s = [x * x for x in w]
    S = sum(s, F(0)); K = sig + S / 2
    d = [s[h] + S / 2 + sig for h in range(u)]
    E = [[(d[i] if i == j else F(0)) + w[i] * w[j] / 2 for j in range(u)] for i in range(u)]
    Ei = inv_mat(E, F)
    p = sum(ci * wi for ci, wi in zip(c, w))
    return [[2 * Ei[i][j] * a[j] * w[j] + 2 * c[i] * w[i] * Ei[i][j] + p * (Ei[i][j] - (F(1) / K if i == j else F(0)))
             for j in range(u)] for i in range(u)]


def blocks_norm(x, t, ww, F=Fr):
    """x_h = s_h/K, t = sigma/K, ww[h][k] = w_h w_k / K (normalised weight products; ww[h][h] = x_h)."""
    u = len(x)
    d = [1 + y for y in x]; rho = [y / (1 + y) for y in x]; g = sum(rho, F(0)); G = 2 + g
    eta2 = [y / (1 + y) ** 2 for y in x]; N = sum(eta2, F(0))
    Lam = [Fr(3, 2) - r for r in rho]
    e = [1 / d[h] - eta2[h] / G for h in range(u)]
    m = [rho[h] + eta2[h] / G for h in range(u)]
    b = [BETA * eta2[h] * (N - eta2[h]) / G ** 2 for h in range(u)]
    B = sum(b, F(0)); Pm = sum((m[h] ** 2 / Lam[h] for h in range(u)), F(0))
    Ma = [[F(0)] * u for _ in range(u)]; Mc = [[F(0)] * u for _ in range(u)]; Mca = [[F(0)] * u for _ in range(u)]
    for h in range(u):
        for k in range(u):
            Ma[h][k] = t / G * rho[h] * rho[k] * ww[h][k]
            if h == k:
                Ma[h][h] += (1 - t) + t * rho[h] - t * (BETA * rho[h] ** 2 * (N - eta2[h]) / G ** 2 + x[h] * e[h] ** 2 / Lam[h])
                Mc[h][h] = ((1 - t) + t * rho[h] + t * eta2[h] / G - t * x[h] / 2
                            - t / 4 * x[h] * (8 * b[h] + B + 4 * e[h] * (e[h] - m[h]) / Lam[h] + Pm))
                Mca[h][h] = (2 - x[h]) * rho[h] / G - t * x[h] * e[h] * (e[h] - m[h] / 2) / Lam[h] \
                    - t * BETA / (2 * G * G) * x[h] * eta2[h] * (N - eta2[h])
            else:
                pi = 1 / (G * d[h] * d[k]) - Fr(1, 2) - Fr(1, 4) * (B + 2 * b[h] + 2 * b[k] + Pm - 2 * e[h] * m[h] / Lam[h] - 2 * e[k] * m[k] / Lam[k])
                Mc[h][k] = ww[h][k] * t * pi
                sg = -rho[k] / G + t / 2 * e[k] * m[k] / Lam[k] - t * BETA / G ** 2 * eta2[k] * (eta2[h] + (N - eta2[k]) / 2)
                Mca[h][k] = ww[h][k] * sg
    return Ma, Mc, Mca
