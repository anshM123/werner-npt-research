"""An exact m = 3 operator in case 1a with ratio above 4 in (star), and above 2 in the bound
sqrt(Gamma_X) - sqrt(Sigma_full(X)) <= 2 sigma_2(X).

One free copy (states |0>, |1>), three sector copies with planes K_s = span{|0>, |1>} and e = |2,2,2>, so A is three
qubits. With a_k, b_k in A (Gaussian integers / 1024),
    r_1 = |0> a_1,  l_1 = |1> b_1,  r_2 = |0> a_2 + |1> e,  l_2 = |1> b_2 + |0> e,  C = r_1 l_1^* + r_2 l_2^*.
Then X = |0><1| (x) D with D = a_1 b_1^* + a_2 b_2^*, P = |0><0| (x) a_2, Q = |1><1| (x) b_2^*, Y = |1><0|, so
    q_1(Y) = 1,  q~(P) = |a_2|^2/2,  q~(Q) = |b_2|^2/2,  Gamma_X = 2|D|^2 - |Tr D|^2,
    Sigma_full = sum over the 27 Pauli strings sigma_mu on A of |Tr(sigma_mu D)|^2.
The range and co-range contain r_1 and l_1 in C^3 (x) A (case 1a). All arithmetic is exact (integers and fractions)."""
import itertools
import math
from fractions import Fraction

A_cols = [(324+868j, -1+1j), (177-115j, 2-3j), (345+31j, -4-2j), (18-54j, 4-2j), (119-26j, -2j), (-44-4j, 5j),
          (48-50j, 1-5j), (10-10j, -1j)]
B_cols = [(-807, 41), (-175+317j, -161-12j), (133+279j, -154+52j), (97+35j, -38-439j), (-66+353j, -12+50j),
          (104+63j, 574+187j), (121+12j, -421-270j), (18-39j, -406-136j)]


def gauss(z):
    z = complex(z)
    return int(z.real), int(z.imag)


def mul(x, y):
    return x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0]


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def conj(x):
    return x[0], -x[1]


def n2(x):
    return x[0]**2 + x[1]**2


a1 = [gauss(r[0]) for r in A_cols]
a2 = [gauss(r[1]) for r in A_cols]
b1 = [gauss(r[0]) for r in B_cols]
b2 = [gauss(r[1]) for r in B_cols]
Dint = [[add(mul(a1[i], conj(b1[j])), mul(a2[i], conj(b2[j]))) for j in range(8)] for i in range(8)]  # 2^20 D
scale = Fraction(1, 2**40)
trD = (sum(Dint[i][i][0] for i in range(8)), sum(Dint[i][i][1] for i in range(8)))
Gamma = scale * (2 * sum(n2(Dint[i][j]) for i in range(8) for j in range(8)) - n2(trD))
pauli = {1: [[(0, 0), (1, 0)], [(1, 0), (0, 0)]], 2: [[(0, 0), (0, -1)], [(0, 1), (0, 0)]],
         3: [[(1, 0), (0, 0)], [(0, 0), (-1, 0)]]}
Sigma = 0
for mu in itertools.product((1, 2, 3), repeat=3):
    tr = (0, 0)
    for i in range(8):
        for j in range(8):
            entry = (1, 0)                    # (sigma_mu)_{ji}
            for k in range(3):
                entry = mul(entry, pauli[mu[k]][(j >> (2 - k)) & 1][(i >> (2 - k)) & 1])
            tr = add(tr, mul(entry, Dint[i][j]))
    Sigma += n2(tr)
Sigma = scale * Sigma
y = 1
p = Fraction(sum(map(n2, a2)), 2 * 1024**2)
q = Fraction(sum(map(n2, b2)), 2 * 1024**2)
print(f"Gamma_X = {Gamma.numerator}/2^40,  Sigma_full = {Sigma.numerator}/2^40,  q_1(Y) = {y},  "
      f"q~(P) = {p.numerator}/2^{p.denominator.bit_length() - 1},  q~(Q) = {q.numerator}/2^{q.denominator.bit_length() - 1}")


def exceeds(K):
    """sqrt(y) (sqrt(Gamma) - sqrt(Sigma)) > K sqrt(p q), decided exactly."""
    T = y * (Gamma - Sigma) - K * K * p * q
    return T > 0 and T * T > 4 * K * K * y * Sigma * p * q


ratio = math.sqrt(y) * (math.sqrt(Gamma) - math.sqrt(Sigma)) / math.sqrt(p * q)
print(f"ratio sqrt(q_1(Y)) (sqrt(Gamma) - sqrt(Sigma)) / sqrt(q~(P) q~(Q)) = {ratio:.12f}")
print(f"exceeds 4 (exact): {exceeds(4)};  exceeds 2^3 - 2 = 6 (exact): {exceeds(6)}")
# |a_2||b_2| >= sigma_2(X), since X minus the a_1 b_1^* term has rank one; and 4 |a_2|^2 |b_2|^2 = 16 p q
print(f"so sqrt(Gamma) - sqrt(Sigma) > 2 |a_2| |b_2| >= 2 sigma_2(X): {exceeds(4)};  "
      f"(sqrt(Gamma) - sqrt(Sigma)) / (|a_2| |b_2|) = {ratio / 2:.12f}")
