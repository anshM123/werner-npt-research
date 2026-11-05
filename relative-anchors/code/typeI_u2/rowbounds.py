"""Per-coordinate bounds for the blocks of Phi1 (NOTES.md, Section 4).  Generic in the number type: works with
floats, Fractions and SymPy expressions.  Parameters: S = 2(1 - t) in (0, 2), xi in (0, 1] (x = S xi is the
coordinate's weight s_h/K), theta in [0, 1] (g' = rho(S') + theta (S' - rho(S')), S' = S - x, the rest's mass),
and nu in {0, 1} (N - eta_h^2 = nu g', the two extreme values of a quantity entering linearly).
Returns (all already divided by S/2 = 1 - t):
  LBc   lower bound of  [Mc_hh - Gershgorin row sum] / (1 - t)
  ROWx  upper bound of  [weighted row sum of |Mca|, row h] / (1 - t)   (both signs of Mca_hh: ROWx_plus, ROWx_minus)
  COLx  upper bound of  [weighted column sum of |Mca|, column h] / (1 - t)
"""
beta_num, beta_den = 12, 5


def bounds(S, xi, theta, nu, one=1):
    one = one
    half = one / 2
    beta = one * beta_num / beta_den
    t = 1 - S / 2
    x = S * xi
    d = 1 + x
    rho = x / d
    eta2 = x / (d * d)
    Sp = S - x
    rSp = Sp / (1 + Sp)
    gp = rSp + theta * (Sp - rSp)
    g = rho + gp
    G = 2 + g
    Lam = 3 * half - rho
    e = 1 / d - eta2 / G
    m = rho + eta2 / G
    Nm = nu * gp                      # N - eta_h^2
    N = eta2 + Nm
    # upper bounds of the global moments
    b_h = beta * eta2 * Nm / (G * G)
    B_ub = beta * (eta2 + gp) ** 2 / (G * G)            # B <= beta N^2 / G^2, N <= eta2 + g'
    cm = (1 + 1 / G) ** 2
    LamS = 3 * half - rSp
    Dp = Sp - gp                                          # = sum_{k != h} x_k rho_k >= sum_{k != h} rho_k^2
    Pm_rest = cm * Dp / LamS                              # sum_{k != h} m_k^2 / Lam_k
    Pm_ub = m * m / Lam + Pm_rest
    xb_rest = beta * (eta2 + gp) * Dp / (G * G)           # sum_{k != h} x_k b_k
    absP4_rest = (1 + 1 / G) * (2 + Sp * (1 + 1 / G)) * Dp / LamS   # sum_{k != h} x_k |m_k (2 e_k - m_k)| / Lam_k
    # Mc_hh (lower bound: the moments enter with a minus sign)
    Mc_hh = (1 - t) + t * rho + t * eta2 / G - t * x / 2 - t / 4 * x * (8 * b_h + B_ub + 4 * e * (e - m) / Lam + Pm_ub)
    # Gershgorin row sum (weights w): t * sum_{k != h} x_k (P1 + P3 + |P4|)
    sP1 = Sp / 2 - gp / (d * G)
    sP3 = (Sp * (B_ub + 2 * b_h + Pm_ub - m * m / Lam) + 2 * xb_rest) / 4
    mh = m * (2 * e - m)
    sP4_a = (Sp * mh / Lam + absP4_rest) / 4          # valid if m(2e - m) >= 0
    sP4_b = (-Sp * mh / Lam + absP4_rest) / 4         # valid if m(2e - m) <= 0
    LBc_a = (Mc_hh - t * (sP1 + sP3 + sP4_a)) / (1 - t)
    LBc_b = (Mc_hh - t * (sP1 + sP3 + sP4_b)) / (1 - t)
    # Mca diagonal (exact given N) and off-diagonal sums
    Mca_hh = (2 - x) * rho / G - t * x * e * (e - m / 2) / Lam - t * beta / (2 * G * G) * x * eta2 * Nm
    row_rest = Dp / G + t * beta / (G * G) * (eta2 + (eta2 + gp) / 2) * Dp
    col_rest = Sp * rho / G + t * beta / (G * G) * eta2 * (Dp + Sp * gp / 2)
    ROW_p = (Mca_hh + row_rest) / (1 - t)
    ROW_m = (-Mca_hh + row_rest) / (1 - t)
    COL_p = (Mca_hh + col_rest) / (1 - t)
    COL_m = (-Mca_hh + col_rest) / (1 - t)
    return dict(LBc_a=LBc_a, LBc_b=LBc_b, mh=mh, ROW_p=ROW_p, ROW_m=ROW_m, COL_p=COL_p, COL_m=COL_m)
