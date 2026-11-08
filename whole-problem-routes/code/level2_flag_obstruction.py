"""Exact check of Remark 7.15': the second level of the symmetric-extension hierarchy fails for P_{3/4} at 17 copies.

Seed: the three-copy PPT state rho of Remark 7.15 (Tr[C^{(x)3} rho] = -3/832).  Add k = 14 flag copies.  On each of the two
outputs put V b = b (x) |1^k>, e = |000> (x) |0^k>; fix the input flags to |0^k>; W b = (V b (x) e + e (x) V b)/sqrt 2;
Q = projector onto Sym^2(V C^27).  Then
    X = (I (x) W) rho (I (x) W)^* + alpha I (x) Q + beta I (x) |ee><ee|,   alpha = 1/3499200, beta = 2187/4,
is bosonic on the outputs and positive under every partial transpose.  After transposing one output, the only
off-diagonal block links the sector with no active output (diagonal beta) to the sector with two active outputs
(diagonal alpha (I + Omega)/2 >= alpha/2); that block is T/2 with T T^* = Tr_B rho^2, so alpha*beta >= |T|^2/2 suffices.
With H = C_{3/4}^{(x)17} (x) I on the second output, the three parts of X contribute
    u^k t/2 + s b^k/2,   alpha * 14 (Tr C)^3 u^k,   beta (7/4)^3 b^k,
where u = <0 1|C|0 1> = 3/4 (input flag 0, output flag 1), b = <0 0|C|0 0> = 1/4, s = <000|Tr_A[C^{(x)3}]|000>/27 and
14 = (27 + 1)/2 is the one-output marginal of Q.  Symmetrizing the real indices multiplies the seed terms by 2/3 and
turns Q into (2/3)Q + (1/3)Omega, whose one-output marginal is (27 + 2)/3.  usage: python level2_flag_obstruction.py"""
from fractions import Fraction as Fr
import itertools

# local operator C_{3/4} = (3/4) I - (F + Omega)/4 on C^3 (x) C^3, exact
idx = [(i, j) for i in range(3) for j in range(3)]
def C(p, q):
    (i, j), (k, l) = p, q
    val = Fr(3, 4) * (i == k and j == l)
    val -= Fr(1, 4) * (i == l and j == k)          # swap F
    val -= Fr(1, 4) * (i == j and k == l)          # Omega = |sum_a aa><sum_b bb|
    return val
u = C((0, 1), (0, 1)); b = C((0, 0), (0, 0))
trC = sum(C(p, p) for p in idx)
trA = [[sum(C((a, i), (a, j)) for a in range(3)) for j in range(3)] for i in range(3)]
assert u == Fr(3, 4) and b == Fr(1, 4) and trC == Fr(21, 4)
assert all(trA[i][j] == (Fr(7, 4) if i == j else 0) for i in range(3) for j in range(3))

# seed of Remark 7.15 (same construction as ppt_state_n3.py)
G = [[Fr(1, 3), Fr(5, 3), Fr(-1)], [Fr(1, 3), Fr(1, 6), Fr(1, 2)], [Fr(-1, 3), Fr(5, 6), Fr(1, 2)]]
dims = [1, 5, 3]; cval = [Fr(-1, 4), Fr(1, 2), Fr(1)]; marg = [Fr(1, 3), Fr(5, 3), Fr(1)]   # Tr_B of each local projector / I
w = {'OOO': Fr(1, 27), 'OOs': Fr(10, 117), 'Oaa': Fr(4, 13), 'sss': Fr(200, 351)}
pats = list(itertools.product(range(3), repeat=3))
key = lambda J: ''.join(sorted('Osa'[j] for j in J))
count = {kk: sum(1 for J in pats if key(J) == kk) for kk in w}
p = {J: (w[key(J)] / count[key(J)] if key(J) in w else Fr(0)) for J in pats}
r = {J: p[J] / (dims[J[0]] * dims[J[1]] * dims[J[2]]) for J in pats}
rG = {K: sum(G[K[0]][J[0]] * G[K[1]][J[1]] * G[K[2]][J[2]] * r[J] for J in pats) for K in pats}
prod = lambda f, J: f[J[0]] * f[J[1]] * f[J[2]]
t = sum(p[J] * prod(cval, J) for J in pats)
rhoA = sum(r[J] * prod(marg, J) for J in pats)
trB_rho2 = sum(r[J] ** 2 * prod(marg, J) for J in pats)
assert sum(p.values()) == 1 and all(v >= 0 for v in r.values()) and all(v >= 0 for v in rG.values())
assert t == Fr(-3, 832) and rhoA == Fr(1, 27) and trB_rho2 == Fr(109, 369603)

# feasibility (Schur bound) and value at 17 copies
d, k = 27, 14
alpha, beta = Fr(1, 3499200), Fr(2187, 4)
assert alpha * beta == Fr(1, 6400) and alpha * beta > trB_rho2 / 2
s = Fr(7, 4) ** 3 / d
J = Fr(7, 4) ** 3
value = u ** k * (t / 2 + alpha * Fr(d + 1, 2) * trC ** 3) + b ** k * (s / 2 + beta * J)
real_value = u ** k * (Fr(2, 3) * t / 2 + alpha * Fr(d + 2, 3) * trC ** 3) + b ** k * (Fr(2, 3) * s / 2 + beta * J)
print('seed: t =', t, ', rho_A = I/27:', rhoA == Fr(1, 27), ', Tr_B rho^2 =', trB_rho2, 'I')
print('Schur: alpha*beta =', alpha * beta, '>', trB_rho2 / 2)
print('Tr[(C^{(x)17} (x) I) X] =', value, '=', float(value))
print('real symmetrized value =', real_value, '=', float(real_value))
ok = value == Fr(-52542910001, 4824107266867200) and real_value == Fr(-97742757347, 28944643601203200) and value < 0
print('PASS' if ok else 'FAIL')
