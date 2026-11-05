"""Exact check: for a perp w_U, Psi_a = sum_{r,h in U} a_r w_h P_{r,h} has 2<Psi_a, Y Psi_a> = s_T <a, E_0 E_T^{-1} a>
(E_0 = diag(s_h + s_U/2) + w_U w_U^T/2, E_T = s_T I + E_0), Fractions, k = 4, 5, 6."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cf import TypeI, inv_mat
ok = True
for w, T, a in [((Fr(1), Fr(2), Fr(3, 2), Fr(1, 2)), [0, 1], None), ((Fr(2, 3), Fr(1), Fr(3, 2), Fr(1, 2), Fr(5, 4)), [1], None),
                ((Fr(2, 3), Fr(1), Fr(3, 2), Fr(1, 2), Fr(5, 4), Fr(7, 3)), [0, 4], None)]:
    ti = TypeI(w, T); M = ti.build(); U = ti.U; s = ti.s; wU = [ti.w[h] for h in U]
    u = len(U)
    # a perp w_U: a = e_0 * w_1 - e_1 * w_0 (+ e_2 * small) adjusted
    a = [Fr(0)] * u; a[0] = wU[1]; a[1] = -wU[0]
    if u > 2:
        a[2] = Fr(1); a[0] -= wU[2] / wU[0] * 0; a[0] += -wU[2] / wU[0]  # keep a.w = 0
    assert sum(x * y for x, y in zip(a, wU)) == 0
    v = [Fr(0)] * ti.n
    for i, r in enumerate(U):
        for j, h in enumerate(U):
            v[ti.idx[('P', r, h)]] = a[i] * wU[j]
    val = sum(v[i] * M[i][j] * v[j] for i in range(ti.n) for j in range(ti.n))
    sT = sum(s[t] for t in T); sU = sum(s[h] for h in U)
    E0 = [[(s[U[i]] + sU / 2 if i == j else 0) + wU[i] * wU[j] / 2 for j in range(u)] for i in range(u)]
    ET = [[E0[i][j] + (sT if i == j else 0) for j in range(u)] for i in range(u)]
    ETi = inv_mat(ET, Fr)
    E0a = [sum(E0[i][j] * a[j] for j in range(u)) for i in range(u)]
    pred = sT * sum(a[i] * ETi[i][j] * E0a[j] for i in range(u) for j in range(u))
    good = val == pred
    ok &= good
    print(f"k={len(w)} T={T}: 2<Psi_a,Y Psi_a> = {val} ; s_T <a,E_0 E_T^-1 a> = {pred} : {good}")
print("PASS" if ok else "FAIL")
