"""Exact check of Remark 7.15's PPT state.  rho = sum_J p_J P_J / dim(P_J) over patterns J in {Omega-hat, s0, an}^3 of the
three pairs (input s, output s), with the weight of each multiset of letters spread evenly over its arrangements:
OOO 1/27, OOs 10/117, Oaa 4/13, sss 200/351.  Checks rho >= 0, rho^Gamma >= 0 (partial transpose on the outputs, through
the matrix G of Proposition 7.3) and Tr[C_{3/4}^{(x)3} rho] = -3/832.  usage: python ppt_state_n3.py"""
from fractions import Fraction as Fr
import itertools
G = [[Fr(1, 3), Fr(5, 3), Fr(-1)], [Fr(1, 3), Fr(1, 6), Fr(1, 2)], [Fr(-1, 3), Fr(5, 6), Fr(1, 2)]]   # columns: images of P_Om, P_s0, P_an
dims = [1, 5, 3]
c = [Fr(-1, 4), Fr(1, 2), Fr(1)]          # eigenvalues of C_{3/4} on P_Om, P_s0, P_an
w = {'OOO': Fr(1, 27), 'OOs': Fr(10, 117), 'Oaa': Fr(4, 13), 'sss': Fr(200, 351)}
lab = 'Osa'
pats = list(itertools.product(range(3), repeat=3))
key = lambda J: ''.join(sorted(lab[j] for j in J))
count = {k: sum(1 for J in pats if key(J) == k) for k in w}
p = {J: (w[key(J)] / count[key(J)] if key(J) in w else Fr(0)) for J in pats}
r = {J: p[J] / (dims[J[0]] * dims[J[1]] * dims[J[2]]) for J in pats}
rG = {K: sum(G[K[0]][J[0]] * G[K[1]][J[1]] * G[K[2]][J[2]] * r[J] for J in pats) for K in pats}
assert sum(p.values()) == 1
ok_psd = all(v >= 0 for v in r.values()); ok_ppt = all(v >= 0 for v in rG.values())
val = sum(p[J] * c[J[0]] * c[J[1]] * c[J[2]] for J in pats)
print('trace 1:', sum(p.values()) == 1, '; rho >= 0:', ok_psd, '; rho^Gamma >= 0:', ok_ppt)
print('Tr[C_{3/4}^{(x)3} rho] =', val)
print('PASS' if ok_psd and ok_ppt and val == Fr(-3, 832) else 'FAIL')
