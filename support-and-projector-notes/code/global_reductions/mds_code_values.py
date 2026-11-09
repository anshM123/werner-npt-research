"""Section 7: q_N on two-dimensional subspaces of qutrit quantum MDS codes [[n,1,(n+1)/2]]_3, from the weight enumerators,
which MDS parameters determine: the normalizer is an additive MDS code [n, (n+1)/2] over F_9 and the stabilizer group an
additive MDS code [n, (n-1)/2]; A_w = C(n,w) sum_j (-1)^j C(w,j) (9^{w-d+1-j} - 1).  Assuming the eight logical cosets have
equal weight enumerators (true for the five-qutrit code), q_N(t, Pi) = (4/3^n)[W_S + (W_N - W_S)/16], W_X = sum (1-3t)^{n - wt}.
Exact rational output.  usage: python mds_code_values.py"""
from fractions import Fraction as Fr
from math import comb


def mds_weights(n, K, Q=9):
    d = n - K + 1
    A = {0: 1}
    for w in range(d, n + 1):
        A[w] = comb(n, w) * sum((-1) ** j * comb(w, j) * (Q ** (w - d + 1 - j) - 1) for j in range(w - d + 1))
    assert sum(A.values()) == Q ** K
    return A


for n in (5, 7, 9):
    AN, AS = mds_weights(n, (n + 1) // 2), mds_weights(n, (n - 1) // 2)
    for t in (Fr(1, 2), Fr(2, 5)):
        lam = 1 - 3 * t
        WN = sum(a * lam ** (n - w) for w, a in AN.items())
        WS = sum(a * lam ** (n - w) for w, a in AS.items())
        q = Fr(4, 3 ** n) * (WS + (WN - WS) / 16)
        print(f"[[{n},1,{(n + 1) // 2}]]_3, t = {t}: q_n(t, Pi)/||Pi||^2 = {q / 2} = {float(q / 2):.6f}")
