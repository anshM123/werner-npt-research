"""Per-coordinate bounds (NOTES.md, Section 4), generic number type (float, Fraction, SymPy).
Parameters: S = 2(1 - t) in (0, 2); xi in (0, 1] (x = S xi); theta in [0, 1]
(g' = rho(S') + theta (S' - rho(S')), S' = S - x); nu in [0, 1] (N - eta_h^2 = nu * Nub', where Nub' is the upper
bound of N - eta_h^2 below).  All outputs are divided by (1 - t) = S/2.
  LB    : lower bound of Mc_hh + t sum_{k != h} x_k pi_hk   (signed Gershgorin row; valid when all pi_hk <= 0)
  ROWp, ROWm, COLp, COLm : upper bounds of the weighted row / column sums of |Mca| for Mca_hh >= 0 / <= 0.
"""


def bounds(S, xi, theta, nu, one=1):
    half = one / 2
    beta = one * 12 / 5
    t = 1 - S / 2
    x = S * xi
    d = 1 + x
    rho = x / d
    eta2 = x / (d * d)
    Sp = S - x
    rSp = Sp / (1 + Sp)
    gp = rSp + theta * (Sp - rSp)
    Dp = Sp - gp                          # = sum_{k != h} x_k rho_k  (exact identity)
    g = rho + gp
    G = 2 + g
    Lam = 3 * half - rho
    e = 1 / d - eta2 / G
    m = rho + eta2 / G
    Nub_p = gp - Dp / (1 + Sp)            # N - eta_h^2 = g' - sum rho_k^2 <= g' - Dp/(1 + S')
    Nm = nu * Nub_p
    N = eta2 + Nm
    Nub = eta2 + Nub_p
    b_h = beta * eta2 * Nm / (G * G)
    B_ub = beta * Nub * Nub / (G * G)
    Pm_ub = m * m / Lam + 2 * (1 + 1 / G) ** 2 * Dp / 3
    xb_rest = beta * Nub * Dp / (G * G)   # sum_{k != h} x_k b_k <= beta N sum rho_k^2 / G^2, sum rho_k^2 <= Dp
    br = S * (B_ub + Pm_ub) + (8 * x + 2 * Sp) * b_h + 4 * x * e * (e - m) / Lam - 2 * Sp * e * m / Lam + 2 * xb_rest
    LB = ((1 - t) ** 2 + t * rho + t * g / (G * d) - t / 4 * br) / (1 - t)
    Mca_hh = (2 - x) * rho / G - t * x * e * (e - m / 2) / Lam - t * beta / (2 * G * G) * x * eta2 * Nm
    row_rest = Dp / G + t * beta / (G * G) * (eta2 + Nub / 2) * Dp
    col_rest = Sp * rho / G + t * beta / (G * G) * eta2 * (Dp + Sp * Nub_p / 2)
    return dict(LB=LB, ROWp=(Mca_hh + row_rest) / (1 - t), ROWm=(-Mca_hh + row_rest) / (1 - t),
                COLp=(Mca_hh + col_rest) / (1 - t), COLm=(-Mca_hh + col_rest) / (1 - t), Mca_hh=Mca_hh)


def pair_Q(xh, xk, g, one=1):
    """Q(x_h, x_k, g) = 1/(G d_h d_k) - 1/2 + (1/4)[m_h(2e_h - m_h)/Lam_h + m_k(2e_k - m_k)/Lam_k]; pi_hk <= Q."""
    G = 2 + g
    out = 1 / (G * (1 + xh) * (1 + xk)) - one / 2
    for x in (xh, xk):
        d = 1 + x; rho = x / d; eta2 = x / (d * d)
        e = 1 / d - eta2 / G; m = rho + eta2 / G; Lam = one * 3 / 2 - rho
        out += m * (2 * e - m) / (4 * Lam)
    return out
