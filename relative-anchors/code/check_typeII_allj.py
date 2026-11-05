"""Checks for Proposition T2' (README, Section 3d): type II blocks of every (R_inf)_j, all weights.
(a) One copy, exact: for v = |0>, Y_v = |10><10| + |20><20| + (|12> - |21>)(<12| - <21|) (row (x) column), its kernel is
    (v (x) C^3) + Sym^2(v^perp) (dimension 6), and A^{-1/2} R_v^Gamma A^{-1/2} has spectrum {1, -1/3}.
(b) For p = sum_r w_r |1_r> on k copies and u = |2@t>, the compression of Y_p to u (x) (all columns) is
    >= (4/9) 4^{-k} min{x_t, 1 - x_t} I: direct construction (check_bath_reduction.group_ops), k = 2..5."""
import os, sys, itertools
import numpy as np
import sympy as sp
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', '..', 'single-excitation-range', 'code'))
from check_bath_reduction import group_ops

H = sp.Rational(1, 2)
def m1(a, y, b, x): return (1 if (a == b and y == x) else 0) - (H if (a == y and b == x) else 0)
st = [(a, y) for a in range(3) for y in range(3)]; ix = {s: i for i, s in enumerate(st)}
M = sp.Matrix(9, 9, lambda i, j: m1(st[i][0], st[i][1], st[j][0], st[j][1]))
A = sp.Matrix(9, 9, lambda i, j: m1(st[j][0], st[i][1], st[i][0], st[j][1]))       # A[(a,y'),(b,y)] = M[(b,y'),(a,y)]
V = sp.Matrix(9, 3, lambda i, x: 1 if (st[i][0] == 0 and st[i][1] == x) else 0)     # p = |0>
D = V.T * M * V; R = M * V * D.inv() * V.T * M
B = sp.Matrix(9, 9, lambda i, j: R[ix[(st[j][0], st[i][1])], ix[(st[i][0], st[j][1])]])
Y = A - B
def ket(a, y): v = sp.zeros(9, 1); v[ix[(a, y)]] = 1; return v
d = ket(1, 2) - ket(2, 1)
Yclaim = ket(1, 0) * ket(1, 0).T + ket(2, 0) * ket(2, 0).T + d * d.T
oka = (Y - Yclaim) == sp.zeros(9, 9) and Y.rank() == 3
evT = (A.inv() * B).eigenvals()          # same spectrum as A^{-1/2} B A^{-1/2}
oka &= set(evT.keys()) == {1, sp.Rational(-1, 3)}
print(f"(a) Y_|0> closed form, rank 3 (kernel dimension 6), spec T = {dict(evT)}: {'PASS' if oka else 'FAIL'}")
ok = oka
rng = np.random.default_rng(4); worst = 9.0; cnt = 0
for k in (2, 3, 4, 5):
    for trial in range(5 if k < 5 else 3):
        w = rng.random(k) ** rng.choice([1, 3, 6]) + 1e-3
        if trial == 0: w = np.array([1.0] + [1e-2] * (k - 1))
        pG = {}
        for t in range(k):
            x = [0] * k; x[t] = 1; pG[tuple(x)] = w[t]
        _, Aent, Bent = group_ops(pG, k)
        s = w ** 2; xs = s / s.sum(); cols = list(itertools.product(range(3), repeat=k))
        for t in range(k):
            u = tuple(2 if c == t else 0 for c in range(k))
            Yu = np.array([[Aent(u, yp, u, y) - Bent(u, yp, u, y) for y in cols] for yp in cols])
            lm = np.linalg.eigvalsh((Yu + Yu.T) / 2)[0]
            worst = min(worst, lm / (4 / 9 * 4.0 ** (-k) * min(xs[t], 1 - xs[t]))); cnt += 1
print(f"(b) {cnt} rows u = 2@t, k = 2..5: smallest ratio lambda_min / bound = {worst:.3f} (must be >= 1)")
ok &= worst >= 1
print("PASS" if ok else "FAIL")
