"""Closed form of the letter-1-free type I block of (R_inf)_j (|T| = j + 1 copies carrying the letter 2), every j.

Normalisation |s| = 1, x_t = s_t.  T = (t_0, ..., t_j), a_i = x_{t_i}, x_T = sum a_i, U = [k] minus T,
c = 2 - x_T, t_u = c - 2 x_u (u in U), Phi(d) = sum_{u in U} x_u / (t_u + d).
States X_i = (2@t_i, 2 on T minus t_i), i = 0..j, then X_vac = (vac, 2 on T).  Units 2^{k-j}.
  M_ii = 2,  M_il = (1 + G_i)(1 + G_l)/(2 + G_i + G_l + Phi(a_i + a_l)),  G_i = a_i/(c - a_i + a_l)  (pair (i, l)),
  M_i,vac = -A_i,  A_i = 2(1 + g_i)/(2 + g_i + Phi(a_i)),  g_i = a_i/(c - a_i),
  M_vac,vac = 2/Gamma,  Gamma = 1 + Phi(0)/2.
The A-part (no E^{-1}) of the same block: 2 on the diagonal, 1/2 between X_i and X_l, -1 between X_i and X_vac,
2 at (vac, vac).  bscale multiplies the B-part (A - M) for negative controls.
Works with any field elements supporting + - * / (fractions.Fraction, flint.fmpq, floats)."""
import itertools


def phi(xU, c, d):
    tot = 0 * c
    for xu in xU:
        tot = tot + xu / (c - 2 * xu + d)
    return tot


def block(aT, xU, bscale=1):
    """aT: weights of the copies of T (list), xU: weights of the copies of U; sum of all = 1."""
    n = len(aT)
    xT = sum(aT[1:], aT[0])
    c = 2 - xT
    M = [[0 * c] * (n + 1) for _ in range(n + 1)]
    Aprt = [[0 * c] * (n + 1) for _ in range(n + 1)]
    for i in range(n):
        M[i][i] = 2 * (c / c)
        Aprt[i][i] = M[i][i]
    for i, l in itertools.combinations(range(n), 2):
        a, b = aT[i], aT[l]
        Gi = a / (c - a + b)
        Gl = b / (c - b + a)
        M[i][l] = M[l][i] = (1 + Gi) * (1 + Gl) / (2 + Gi + Gl + phi(xU, c, a + b))
        Aprt[i][l] = Aprt[l][i] = (c / c) / 2
    for i in range(n):
        a = aT[i]
        g = a / (c - a)
        M[i][n] = M[n][i] = -2 * (1 + g) / (2 + g + phi(xU, c, a))
        Aprt[i][n] = Aprt[n][i] = -(c / c)
    M[n][n] = 2 / (1 + phi(xU, c, 0 * c) / 2)
    Aprt[n][n] = 2 * (c / c)
    if bscale != 1:
        M = [[Aprt[i][l] - bscale * (Aprt[i][l] - M[i][l]) for l in range(n + 1)] for i in range(n + 1)]
    return M


def eliminate_vacuum(M):
    """F = C - (Gamma/2) A A^T: the Schur complement of the vacuum entry."""
    n = len(M) - 1
    return [[M[i][l] - M[i][n] * M[n][l] / M[n][n] for l in range(n)] for i in range(n)]


def states(k, T):
    """ylib states of the block, in the order X_0..X_j, X_vac."""
    def st(ds):
        z = [0] * k
        for t in ds:
            z[t] = 2
        return tuple(z)
    out = [(st([T[i]]), st([T[m] for m in range(len(T)) if m != i])) for i in range(len(T))]
    out.append(((0,) * k, st(T)))
    return out
