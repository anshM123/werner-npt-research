"""Exact check (python-flint) of the symmetric free energy and its square identity (NOTES.md, Section 2), against the
block of Theorem C built exactly from the closed form (freeblock.py), at rational weights with free squared mass 1.

For symmetric X on the free copies, F(X) = min_q z^T C z (z = (X_U, q)).  With x_i = w_i^2, t = sigma, K = t + 1/2,
f_i = 1/(K + x_i), F = E f, alpha = 1/(2 + F), r_i = x_i f_i, B = E r, D = 2K + 1, the decomposition
X = Y + w_i w_j (l_i + l_j + p), Yw = 0, E l = 0, and
  A0 = E[m_L l^2] + 2t^2 alpha (E l f)^2 + c_p p^2 - 8 t^2 alpha p E l f,  m_L = 2t(x + 1/2 + 2 alpha x) f,
  c_p = t - 2 t^2 alpha F - t^2/(2K^2),
  H+(Y) = ||Y||^2/2 + sum_i K f_i Y_ii^2 + alpha sum_ij r_i r_j Y_ij^2,
  <Y, R> = sum_i A_i Y_ii + alpha <Y, S>,  A_i = 2t r_i l_i + p t^2 r_i/K,  S_ij = a_i a_j (t(l_i + l_j)/K + p t^2/K^2),
  a_i = w_i r_i,
the claims are (1) F(X) = c (A0 + H+(Y) + 2<Y,R>) for one constant c, and (2) the exact square identity
  A0 + H+ + 2<Y,R> - ||Y||^2/2 = sum_i K f_i (Y_ii + A_i/(K f_i))^2 + alpha sum_ij r_i r_j (Y_ij + S_ij/(r_i r_j))^2
                                  + (t/K) E l^2 + (t alpha/K^2) E[r (l + p t/K)^2] + R* p^2,
  R* = t(t^2 + t + 1/2)/(2K^2(t + 1)) + t^2 B/(4 D K^3 (D - B))."""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from freeblock import free_block, q_eliminated
Q = flint.fmpq

# rational unit vectors (sum of squares 1) for the free weights
UNIT = [[Q(3, 5), Q(4, 5)], [Q(5, 13), Q(12, 13)], [Q(1, 3), Q(2, 3), Q(2, 3)], [Q(2, 7), Q(3, 7), Q(6, 7)],
        [Q(1, 2), Q(1, 2), Q(1, 2), Q(1, 2)], [Q(1, 9), Q(4, 9), Q(8, 9)], [Q(2, 11), Q(6, 11), Q(9, 11)],
        [Q(1, 5), Q(2, 5), Q(2, 5), Q(4, 5)], [Q(2, 15), Q(5, 15), Q(14, 15)], [Q(2, 9), Q(4, 9), Q(4, 9), Q(6, 9), Q(3, 9)]]
for v in UNIT:
    assert sum(x * x for x in v) == 1
rng = random.Random(7)
ok = True
consts = set()
for it, wU in enumerate(UNIT):
    for wt in (Q(1, 3), Q(2), Q(7, 5)):
        u = len(wU)
        Fx = q_eliminated(free_block([wt], wU, Q))
        t = wt * wt
        x = [w * w for w in wU]
        K = t + Q(1, 2)
        f = [1 / (K + xi) for xi in x]
        Fm = sum(x[i] * f[i] for i in range(u))
        al = 1 / (2 + Fm)
        r = [x[i] * f[i] for i in range(u)]
        Bm = sum(x[i] * r[i] for i in range(u))
        Dm = 2 * K + 1
        E = lambda vals: sum(x[i] * vals[i] for i in range(u))
        for trial in range(3):
            X = [[Q(0)] * u for _ in range(u)]
            for i in range(u):
                for j in range(i, u):
                    X[i][j] = X[j][i] = Q(rng.randint(-9, 9), rng.randint(1, 5))
            # our q-eliminated form
            vec = [X[i][j] for i in range(u) for j in range(u)]
            ours = sum(vec[i] * Fx[i][j] * vec[j] for i in range(u * u) for j in range(u * u))
            # decomposition
            Xw = [sum(X[i][j] * wU[j] for j in range(u)) for i in range(u)]
            p = sum(wU[i] * Xw[i] for i in range(u))
            l = [Xw[i] / wU[i] - p for i in range(u)]
            assert E(l) == 0
            Y = [[X[i][j] - wU[i] * wU[j] * (l[i] + l[j] + p) for j in range(u)] for i in range(u)]
            assert all(sum(Y[i][j] * wU[j] for j in range(u)) == 0 for i in range(u))
            mL = [2 * t * (x[i] + Q(1, 2) + 2 * al * x[i]) * f[i] for i in range(u)]
            cp = t - 2 * t * t * al * Fm - t * t / (2 * K * K)
            Elf = E([l[i] * f[i] for i in range(u)])
            A0 = E([mL[i] * l[i] * l[i] for i in range(u)]) + 2 * t * t * al * Elf * Elf + cp * p * p - 8 * t * t * al * p * Elf
            Yn2 = sum(Y[i][j] ** 2 for i in range(u) for j in range(u))
            Hp = Yn2 / 2 + sum(K * f[i] * Y[i][i] ** 2 for i in range(u)) + al * sum(r[i] * r[j] * Y[i][j] ** 2 for i in range(u) for j in range(u))
            Acoef = [2 * t * r[i] * l[i] + p * t * t * r[i] / K for i in range(u)]
            a = [wU[i] * r[i] for i in range(u)]
            Sij = [[a[i] * a[j] * (t * (l[i] + l[j]) / K + p * t * t / (K * K)) for j in range(u)] for i in range(u)]
            YR = sum(Acoef[i] * Y[i][i] for i in range(u)) + al * sum(Y[i][j] * Sij[i][j] for i in range(u) for j in range(u))
            formula_value = A0 + Hp + 2 * YR
            if formula_value != 0:
                consts.add(ours / formula_value)
            Rstar = t * (t * t + t + Q(1, 2)) / (2 * K * K * (t + 1)) + t * t * Bm / (4 * Dm * K ** 3 * (Dm - Bm))
            sos = (sum(K * f[i] * (Y[i][i] + Acoef[i] / (K * f[i])) ** 2 for i in range(u))
                   + al * sum(r[i] * r[j] * (Y[i][j] + Sij[i][j] / (r[i] * r[j])) ** 2 for i in range(u) for j in range(u))
                   + (t / K) * E([li * li for li in l]) + (t * al / (K * K)) * E([r[i] * (l[i] + p * t / K) ** 2 for i in range(u)])
                   + Rstar * p * p)
            ident = (formula_value - Yn2 / 2 == sos)
            ok &= ident
            if not ident:
                print('identity FAILS at', it, wt, trial)
print('ratios our form / formula:', sorted(str(c) for c in consts))
print('square identity:', 'holds at all', len(UNIT) * 3 * 3, 'cases' if ok else 'FAILS')
print('PASS' if ok and len(consts) == 1 else 'FAIL')
