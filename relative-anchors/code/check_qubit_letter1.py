"""Check of the remark after Lemma G of the README (the letter-1 blocks are a qubit problem).  For qubits (d = 2) the
operator Y = (M^{(x)k})^Gamma - (M V D^{-1} V^* M)^Gamma (M = I - |Omega><Omega|/2 on C^2 (x) C^2, V = p (x) I,
D = V^* M^{(x)k} V, Gamma the transpose of the row factor) is built from its definition, compressed to G1_j + G0_{j+1}
for j = 0, ..., k - 1, and compared entry by entry with 2^{-Z} times the letter-1 block of Lemma B (j1lib.py), in exact
rational arithmetic.  python check_qubit_letter1.py K   (K = 3 or 4; a few minutes at K = 4)"""
import itertools, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'middle_j1'))
from fractions import Fraction as Fr
import sympy
from j1lib import J1
k = int(sys.argv[1])
w = [Fr(1), Fr(2, 3), Fr(3, 2), Fr(1, 2)][:k]
strs = list(itertools.product((0, 1), repeat=k))
def M1(a, b, a2, b2):
    return Fr(int(a == a2 and b == b2)) - Fr(1, 2) * int(a == b and a2 == b2)
def Mk(x, y, x2, y2):
    r = Fr(1)
    for t in range(k):
        r *= M1(x[t], y[t], x2[t], y2[t])
        if r == 0:
            return r
    return r
p = {x: (w[x.index(1)] if sum(x) == 1 else Fr(0)) for x in strs}
P1 = [x for x in strs if sum(x) == 1]
# D[y, y'] = sum_{x, x'} p_x p_x' M[(x,y),(x',y')]
D = sympy.Matrix(len(strs), len(strs), lambda i, l: sum(p[x] * p[x2] * Mk(x, strs[i], x2, strs[l]) for x in P1 for x2 in P1))
Di = D.inv()
# (M V)[(x,y), y'] = sum_{x'} M[(x,y),(x',y')] p_x'
def MV(x, y, l):
    return sum(Mk(x, y, x2, strs[l]) * p[x2] for x2 in P1)
idx = {y: i for i, y in enumerate(strs)}
def R(x, y, x2, y2):
    a = [MV(x, y, l) for l in range(len(strs))]
    b = [MV(x2, y2, l) for l in range(len(strs))]
    return sum(a[i] * Di[i, l] * b[l] for i in range(len(strs)) for l in range(len(strs)) if a[i] != 0 and b[l] != 0)
def Y(x, y, x2, y2):   # Gamma transposes the row factor: Y[(x,y),(x2,y2)] = (M - R)[(x2,y),(x,y2)]
    return Mk(x2, y, x, y2) - R(x2, y, x, y2)
J = J1(w)
ok = True
for j in range(0, k):
    Z = k - j
    st = [((1, r), y) for r in range(k) for y in strs if sum(y) == j] + [(None, y) for y in strs if sum(y) == j + 1]
    def full(s):
        a, y = s
        x = tuple(1 if (a is not None and t == a[1]) else 0 for t in range(k))
        return x, y
    Mj = J.matrix(st)
    err = 0
    for i, s in enumerate(st):
        for l, s2 in enumerate(st):
            x, y = full(s); x2, y2 = full(s2)
            v = Y(x, y, x2, y2) * 2 ** Z
            if v != Fr(str(Mj[i][l])):
                err += 1
    print(f'k={k} j={j}: {len(st)} states, {err} entries differ')
    ok &= err == 0
print('PASS' if ok else 'FAIL')
