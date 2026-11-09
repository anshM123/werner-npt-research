"""Stars and the stabilizer-domination sums (Section 7 of NOTES.md, "Correction").

Star K_{1,n} (centre 0, unit weights), N = n + 1, shift d = 0 at the centre and on k leaves, d = 1 on m >= 1 leaves.
Claimed closed form, for v = v_c e_c supported at the centre:
    S(d, v_c e_c) = 3^n [ (1-t)^{k+1} + w^{v_c} + w^{2 v_c} (1-t)^m ].
This script
  (1) checks the closed form against the defining sum S(d, v) = sum_a w^{v.a} (1-3t)^{N - w_d(a)};
  (2) scans all v (fast Fourier transform) for the star patterns and prints the largest ratio |S(d,v)|/S(d,0) for
      N = 3..10 at t = 1/2 and t = 2/5 (tau = t/(1-t));
  (3) proves in exact rational arithmetic that the bound |S(d,v)| <= tau S(d,0) fails at t = 2/5, N = 8, k = 3, m = 4;
  (4) checks the sharp form (7.1'), |S(d,v)| <= S(d,0) - 3^N (1-t)^{N-1}(1-2t), on the star family for many (k, m, t),
      and the polynomial g(x,y) = 3x^2y^2 + x + y - xy - 2xy(x+y) >= 0 on a grid of [0,1]^2.
usage: python star_family.py"""
import itertools
from fractions import Fraction
import numpy as np

om = np.exp(2j * np.pi / 3)


def star(N):
    G = np.zeros((N, N), dtype=np.int64); G[0, 1:] = 1; G[1:, 0] = 1
    return G


def S_direct(G, d, v, t):
    N = len(d); lam = 1 - 3 * t
    A = np.array(list(itertools.product(range(3), repeat=N)))
    w = ((A != 0) | ((A @ G.T + d) % 3 != 0)).sum(1)
    return np.sum(om ** ((A @ v) % 3) * lam ** (N - w))


# (1) closed form
err = 0.0
for k in range(4):
    for m in range(1, 4):
        n = k + m; N = n + 1
        for t in (0.3, 0.4, 0.5):
            for vc in range(3):
                d = np.array([0] + [0] * k + [1] * m); v = np.zeros(N, dtype=np.int64); v[0] = vc
                pred = 3 ** n * ((1 - t) ** (k + 1) + om ** vc + om ** (2 * vc) * (1 - t) ** m)
                err = max(err, abs(S_direct(star(N), d, v, t) - pred))
print(f"(1) closed form vs defining sum (k <= 3, 1 <= m <= 3, t in {{0.3, 0.4, 0.5}}, all v_c): max error {err:.2e}")

# (2) scan over all v for the star shift patterns
for t in (0.5, 0.4):
    lam = 1 - 3 * t; tau = t / (1 - t); row = []
    for N in range(3, 11):
        A = np.array(list(itertools.product(range(3), repeat=N)), dtype=np.int8)
        GA = (A.astype(np.int64) @ star(N).T) % 3; nzA = A != 0; pw = lam ** (N - np.arange(N + 1))
        best = 0.0
        for dc in range(3):
            for m in range(N):
                d = np.array([dc] + [0] * (N - 1 - m) + [1] * m)
                w = (nzA | ((GA + d) % 3 != 0)).sum(1)
                S = np.fft.fftn(pw[w].reshape((3,) * N)).ravel()
                best = max(best, np.abs(S[1:]).max() / S[0].real)
        row.append(f"N={N}: {best:.4f} ({best / tau:.4f} tau)")
    print(f"(2) t = {t}: largest ratio over star shift patterns and all v != 0:")
    print("      " + "; ".join(row))

# (3) exact counterexample to the tau-bound
t = Fraction(2, 5); tau = t / (1 - t); k, m = 3, 4; n = k + m
x = (1 - t) ** (k + 1); y = (1 - t) ** m
S0 = 3 ** n * (x + 1 + y)
mod2 = 3 ** (2 * n) * (x * x + 1 + y * y - x - y - x * y)          # |x + w + w^2 y|^2
print(f"(3) t = 2/5, N = 8, k = 3, m = 4: S(d,0) = {S0}, |S(d,e_c)|^2 = {mod2}")
print(f"    |S(d,e_c)|^2 - tau^2 S(d,0)^2 = {mod2 - tau ** 2 * S0 ** 2} > 0 : {mod2 > tau ** 2 * S0 ** 2};"
      f"  |S(d,e_c)| < S(d,0) : {mod2 < S0 ** 2}")

# (4) the sharp form on stars
worst = np.inf
for kk in range(0, 30):
    for mm in range(1, 30):
        for tt in np.linspace(0.0, 0.5, 51):
            xx = (1 - tt) ** (kk + 1); yy = (1 - tt) ** mm; NN = kk + mm + 1
            lhs = abs(xx + om + om ** 2 * yy)
            rhs = xx + 1 + yy - 3 * (1 - tt) ** (NN - 1) * (1 - 2 * tt)
            worst = min(worst, rhs - lhs)
X, Y = np.meshgrid(np.linspace(0, 1, 401), np.linspace(0, 1, 401))
g = 3 * X ** 2 * Y ** 2 + X + Y - X * Y - 2 * X * Y * (X + Y)
print(f"(4) (7.1') on stars, k < 30, 1 <= m < 30, 51 values of t in [0, 1/2]: min of (right side - left side)/3^n = {worst:.3e}")
print(f"    g(x, y) on a 401 x 401 grid of [0,1]^2: min = {g.min():.3e}")

# (5) the merging budget FB(t) on star pairs, from the closed form (centre merged with an unshifted or a shifted leaf):
#     merging replaces x by x/(1-t) (and y by y/(1-t) if the merged leaf is shifted)
worst = np.inf; arg = None
for t in np.linspace(0.34, 0.5, 33):
    gm = (1 - 2 * t) / 2
    for k in range(40):
        for m in range(1, 40):
            N = k + m + 1; x, y = (1 - t) ** (k + 1), (1 - t) ** m
            Q, X = (x + 1 + y) / 3, (x + om + om ** 2 * y) / 3
            cases = ([((1 - t) ** k, y)] if k >= 1 else []) + [((1 - t) ** k, (1 - t) ** (m - 1))]
            for xp, yp in cases:
                Qp, Xp = (xp + 1 + yp) / 3, (xp + om + om ** 2 * yp) / 3
                marg = (Q - gm * Qp - (1 - t) ** (N - 1) / 2) - (abs(X - gm * Xp) - t * (1 - t) ** (N - 2) / 2)
                if marg / Q < worst:
                    worst, arg = marg / Q, (round(float(t), 3), k, m)
print(f"(5) FB on star pairs (closed form), t in [0.34, 0.5], k < 40, 1 <= m < 40: min relative margin {worst:.3e} at (t, k, m) = {arg}")
