"""Search for a Wu--Zou type certificate of positivity on PSD rank-two operators at n copies (n odd).

For orthonormal u, v put z_T = Tr(rho^u_T rho^v_T), e_A(w) = 1 - Tr((rho^w_A)^2).
Cut inequalities (Wu--Zou Thm 1.1 + Lemma 5.1 applied to X_A = Tr_A |u><v| and the split A^c = K + L):
    g_c = z_A + sqrt(e_A(u) e_A(v)) - z_{A u K} - z_{A u L} >= 0.
Target: Z = sum_T (-1)^{n-|T|} 2^{|T|} z_T >= -sqrt(Qh(u) Qh(v)), Qh(w) = 2^n q_n(P_w).
Certificate: m_c >= 0 with
  (z)  gamma_T = (-1)^{n-|T|} 2^{|T|} - sum_{c: A_c = T} m_c + sum_{c: A_c u K_c = T or A_c u L_c = T} m_c >= 0
       for all T != empty, [n];
  (a)  2 sum_{c: |S cap A_c| odd} m_c <= 3^{|S|} - 1 for all even nonempty S.
Then Z >= -1 - sum m_c sqrt(e e') >= -1 - sqrt(X_u X_v) >= -sqrt(Qh_u Qh_v)."""
import itertools, sys
import numpy as np
from scipy.optimize import linprog


def cuts(n):
    U = frozenset(range(n))
    out = []
    for r in range(1, n - 1):
        for A in itertools.combinations(range(n), r):
            A = frozenset(A); rest = sorted(U - A)
            first = rest[0]
            others = rest[1:]
            for k in range(0, len(others)):
                for Kx in itertools.combinations(others, k):
                    K = frozenset((first,) + Kx); L = frozenset(rest) - K
                    if L:
                        out.append((A, K, L))
    return out


def lp(n, verbose=False):
    U = frozenset(range(n))
    C = cuts(n)
    subsets = [frozenset(s) for r in range(1, n) for s in itertools.combinations(range(n), r)]
    idx = {T: i for i, T in enumerate(subsets)}
    nc = len(C)
    # z-condition: gamma_T >= 0  <=>  sum_{A=T} m - sum_{AuK=T or AuL=T} m <= (-1)^{n-|T|} 2^{|T|}
    Aub, bub = [], []
    for T in subsets:
        row = np.zeros(nc)
        for ci, (A, K, L) in enumerate(C):
            if A == T:
                row[ci] += 1
            if A | K == T:
                row[ci] -= 1
            if A | L == T:
                row[ci] -= 1
        Aub.append(row); bub.append((-1) ** (n - len(T)) * 2 ** len(T))
    # a-condition
    for r in range(2, n + 1, 2):
        for S in itertools.combinations(range(n), r):
            S = frozenset(S)
            row = np.zeros(nc)
            for ci, (A, K, L) in enumerate(C):
                if len(S & A) % 2 == 1:
                    row[ci] += 2
            Aub.append(row); bub.append(3 ** len(S) - 1)
    res = linprog(np.zeros(nc), A_ub=np.array(Aub), b_ub=np.array(bub), bounds=[(0, None)] * nc, method='highs')
    return res, C


if __name__ == '__main__':
    for n in [int(a) for a in sys.argv[1:]]:
        res, C = lp(n)
        print(f'n={n}: cuts {len(C)}, LP status {res.status} ({res.message})', flush=True)
