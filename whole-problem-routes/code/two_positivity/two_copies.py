"""Exact check (rational arithmetic) of Theorem B (two copies):

  for a >= a* = (23 + 6 sqrt 35)/68 = 0.86024...,
  C_a (x) C_a = lam(a) [ M (x) S + S (x) M ] + sum_{j,k} w_jk(a) P_j (x) P_k,
  lam(a) = 12 (1 - a)(a + 1/4) >= 0,  w_jk(a) >= 0,

with C_a = a I - (F + Omega Omega^*)/4, M = I - Omega Omega^*/2, S = P_s0/10 + P_an/6 (the O(3)-twirl of the product
state |e1 e2><e1 e2|), and P_j in {P_Om, P_s0, P_an}.  Here X (x) Y means X on copy 1 (A1 B1) and Y on copy 2 (A2 B2).

What is checked exactly here:
  (1) the projector algebra on C^3 (x) C^3 (P_Om, P_s0, P_an orthogonal projectors summing to I, C_a and M in that basis);
  (2) the three traces <e1 e2|P_j|e1 e2> = (0, 1/2, 1/2), so that the twirl of |e1 e2><e1 e2| is S;
  (3) the 81 x 81 matrix identity above, coefficient by coefficient in a (it is quadratic in a);
  (4) the closed forms of w_jk(a) and their signs on [a*, 1] (polynomial inequalities, by exact root isolation).
The argument that each block is positive on Schmidt-rank-two vectors is in whole-problem-routes/README.md
(Proposition 7.11); (1)-(4) are its only computational inputs.  Prints PASS at the end.
"""
from fractions import Fraction as Fr
import itertools
import sympy as sp

ok = True


def zeros(n):
    return [[Fr(0)] * n for _ in range(n)]


def eye(n):
    Z = zeros(n)
    for i in range(n):
        Z[i][i] = Fr(1)
    return Z


def add(A, B, s=Fr(1)):
    return [[A[i][j] + s * B[i][j] for j in range(len(A))] for i in range(len(A))]


def scal(s, A):
    return [[s * x for x in row] for row in A]


def mul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n) if A[i][k] != 0) for j in range(n)] for i in range(n)]


def kron(A, B):
    n, m = len(A), len(B)
    return [[A[i // m][j // m] * B[i % m][j % m] for j in range(n * m)] for i in range(n * m)]


def eq(A, B):
    return all(A[i][j] == B[i][j] for i in range(len(A)) for j in range(len(A)))


def reorder(W):
    """operator on (A1 B1)(A2 B2), index ((a1 b1) a2 b2), -> index (a1 a2 b1 b2).  (A permutation similarity; the
    identity is checked in the first ordering, which is all that matters.)"""
    return W


# one copy, index 3a + b
I9 = eye(9)
F = zeros(9)
for i in range(3):
    for j in range(3):
        F[3 * i + j][3 * j + i] = Fr(1)
OO = zeros(9)
for i in range(3):
    for j in range(3):
        OO[4 * i][4 * j] = Fr(1)
P_Om = scal(Fr(1, 3), OO)
P_sym = scal(Fr(1, 2), add(I9, F))
P_an = scal(Fr(1, 2), add(I9, F, Fr(-1)))
P_s0 = add(P_sym, P_Om, Fr(-1))
P = [P_Om, P_s0, P_an]
names = ['Om', 's0', 'an']

# (1) projector algebra
c1 = eq(add(add(P_Om, P_s0), P_an), I9)
for X in P:
    c1 &= eq(mul(X, X), X)
for X, Y in itertools.combinations(P, 2):
    c1 &= eq(mul(X, Y), zeros(9))
M = add(I9, OO, Fr(-1, 2))
c1 &= eq(M, add(add(scal(Fr(-1, 2), P_Om), P_s0), P_an))
K0 = scal(Fr(-1, 4), add(F, OO))       # C_a = a I + K0
c1 &= eq(K0, add(add(scal(Fr(-1), P_Om), scal(Fr(-1, 4), P_s0)), scal(Fr(1, 4), P_an)))   # C_a = (a-1)P_Om + (a-1/4)P_s0 + (a+1/4)P_an
print('(1) projector algebra, M = -P_Om/2 + P_s0 + P_an, C_a = (a-1)P_Om + (a-1/4)P_s0 + (a+1/4)P_an:', c1)
ok &= c1

# (2) traces of |e1 e2><e1 e2| (index 3*0 + 1 = 1)
rho = zeros(9)
rho[1][1] = Fr(1)
tr = [sum(mul(rho, X)[i][i] for i in range(9)) for X in P]
c2 = tr == [Fr(0), Fr(1, 2), Fr(1, 2)]
S = add(scal(Fr(1, 10), P_s0), scal(Fr(1, 6), P_an))
c2 &= eq(S, add(add(scal(tr[0] / 1, P_Om), scal(tr[1] / 5, P_s0)), scal(tr[2] / 3, P_an)))
print('(2) <e1e2|P_j|e1e2> =', [str(t) for t in tr], '; twirl = P_s0/10 + P_an/6:', c2)
ok &= c2

# (3) the identity, coefficient by coefficient in a.
a = sp.symbols('a')
lam = 12 * (1 - a) * (a + sp.Rational(1, 4))
cvec = [a - 1, a - sp.Rational(1, 4), a + sp.Rational(1, 4)]
mvec = [sp.Rational(-1, 2), 1, 1]
svec = [0, sp.Rational(1, 10), sp.Rational(1, 6)]
w = [[sp.expand(cvec[j] * cvec[k] - lam * (mvec[j] * svec[k] + svec[j] * mvec[k])) for k in range(3)] for j in range(3)]
# left side: C_a (x) C_a = a^2 I(x)I + a (I(x)K0 + K0(x)I) + K0(x)K0
LHS = [kron(I9, I9), add(kron(I9, K0), kron(K0, I9)), kron(K0, K0)]       # coefficients of a^2, a^1, a^0
MS = add(kron(M, S), kron(S, M))
lam_c = sp.Poly(lam, a).all_coeffs()[::-1]       # [a^0, a^1, a^2]
c3 = True
PP = {(j, k): kron(P[j], P[k]) for j in range(3) for k in range(3)}
for deg in (0, 1, 2):
    R = scal(Fr(int(sp.fraction(sp.Rational(lam_c[deg]))[0]), int(sp.fraction(sp.Rational(lam_c[deg]))[1])), MS)
    for j in range(3):
        for k in range(3):
            cf = sp.Poly(w[j][k], a).all_coeffs()[::-1]
            cf = cf + [0] * (3 - len(cf))
            q = sp.Rational(cf[deg])
            if q != 0:
                R = add(R, PP[(j, k)], Fr(int(q.p), int(q.q)))
    target = LHS[2 - deg]
    c3 &= eq(R, target)
print('(3) 81x81 identity C_a(x)C_a = lam (M(x)S + S(x)M) + sum w_jk P_j(x)P_k, all powers of a:', c3)
ok &= c3

# (4) closed forms and signs
astar = (23 + 6 * sp.sqrt(35)) / 68
closed = {
    (0, 0): (1 - a) ** 2,
    (0, 1): sp.Rational(2, 5) * (1 - a) ** 2,
    (0, 2): 0,
    (1, 1): (sp.Rational(1, 16)) * (272 * a ** 2 - 184 * a - 43) / 5,
    (1, 2): (a + sp.Rational(1, 4)) * ((a - sp.Rational(1, 4)) - sp.Rational(16, 5) * (1 - a)),
    (2, 2): (a + sp.Rational(1, 4)) * (5 * a - sp.Rational(15, 4)),
}
c4 = True
for (j, k), expr in closed.items():
    c4 &= sp.expand(w[j][k] - expr) == 0 and sp.expand(w[k][j] - expr) == 0
c4 &= sp.simplify(sp.Poly(272 * a ** 2 - 184 * a - 43, a).as_expr().subs(a, astar)) == 0
c4 &= bool(sp.nsimplify(astar) > sp.Rational(23, 28)) and bool(astar > sp.Rational(3, 4))
# on [a*, 1]: 272a^2-184a-43 >= 0 (a* is its larger root), a >= 23/28 for w_{s0,an}, a >= 3/4 for w_{an,an}; lam >= 0
roots = sp.solve(272 * a ** 2 - 184 * a - 43, a)
c4 &= max(roots, key=lambda r: float(r)).equals(astar)
c4 &= all(sp.Poly(lam, a).eval(t) >= 0 for t in (astar, 1)) and sp.degree(lam, a) == 2 and sp.LC(sp.Poly(lam, a)) < 0
print('(4) closed forms of w_jk; a* = (23+6 sqrt35)/68 =', float(astar), 'is the larger root of 272a^2-184a-43;',
      'w_jk >= 0 and lam >= 0 on [a*, 1]:', c4)
ok &= c4
for j in range(3):
    print('   w row', names[j], [str(sp.factor(w[j][k])) for k in range(3)])
print('lam(a) =', sp.factor(lam))
print('PASS' if ok else 'FAIL')
