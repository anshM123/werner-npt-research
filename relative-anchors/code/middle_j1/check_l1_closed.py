"""Lemma L1.1 (NOTES.md, Section 4): the vacuum part of the letter-1 block of (R_inf)_1 and its inverse-free form,
exactly against the Lemma B builder j1lib at rational weights with |w| = 1, k = 4..7.  Units 2^{k-1}, unrescaled.
States P_{r,c} = (1@r, 1@c), Q_P = (vac, 1@P), P a pair.
 (i)  Y_QQ = 2 (I + N2/2)^{-1},  N2 = J^- C3^{-1} J^+ on the pairs, C3 = 2(1 - x_a - x_b - x_c) on the triples;
 (ii) for every row r: (I + N1/2) Y_PQ[(r, .), :] = -(I + N1/2) e(r, P) + J^- C2^{-1} T_r |P>, where
      N1 = J^- C2^{-1} J^+ on the singles, C2 = 2(1 - x_a - x_b) on the pairs, e(r, P) = e_{c*} if P = {r, c*} and 0
      otherwise, and T_r |P> = [r in P] (1/2) sum_{u != c*} w_u |{c*, u}> + 2 w_r [r not in P] |P>.
So after the congruence X = X'(I + N1/2) on the P-part (psi_r = (I + N1/2) psi'_r) and the elimination of the vacuum
states with Y_QQ^{-1} = (I + N2/2)/2, the letter-1 block becomes an explicit form in X' with no inverse other than the
Sherman-Morrison inverse of E_1 in its G1-G1 part.  Usage: python check_l1_closed.py"""
import random, itertools
from fractions import Fraction as Fr
from j1lib import J1, block_states, inv

random.seed(9)


def unit_weights(k):
    v = [Fr(random.randint(1, 9), random.randint(1, 9)) for _ in range(k - 1)]
    n2 = sum(a * a for a in v)
    w = [abs(2 * a / (1 + n2)) for a in v] + [abs((1 - n2) / (1 + n2))]
    return w if all(a != 0 for a in w) else unit_weights(k)


ok = True
for k in range(4, 8):
    for trial in range(2):
        w = unit_weights(k)
        x = [a * a for a in w]
        J = J1(w)
        st = block_states(k, ())
        M = J.matrix(st)
        P = [i for i, (a, y) in enumerate(st) if a is not None]
        Q = [i for i, (a, y) in enumerate(st) if a is None]
        pairs = [tuple(c for c in range(k) if st[i][1][c] == 1) for i in Q]
        trip = list(itertools.combinations(range(k), 3))
        # (i)
        N2 = [[Fr(0)] * len(pairs) for _ in pairs]
        for tr in trip:
            c3 = 2 * (1 - sum(x[a] for a in tr))
            for a_, P1 in enumerate(pairs):
                if not set(P1) <= set(tr):
                    continue
                u1 = (set(tr) - set(P1)).pop()
                for b_, P2 in enumerate(pairs):
                    if not set(P2) <= set(tr):
                        continue
                    u2 = (set(tr) - set(P2)).pop()
                    N2[a_][b_] += w[u1] * w[u2] / c3
        A2 = [[(Fr(1) if a_ == b_ else Fr(0)) + N2[a_][b_] / 2 for b_ in range(len(pairs))] for a_ in range(len(pairs))]
        YQQ = [[2 * v for v in row] for row in inv(A2, Fr)]
        good1 = all(M[Q[a_]][Q[b_]] == YQQ[a_][b_] for a_ in range(len(Q)) for b_ in range(len(Q)))
        # (ii)
        N1 = [[(sum(x[u] / (2 * (1 - x[c] - x[u])) for u in range(k) if u != c) if c == cp
                else w[c] * w[cp] / (2 * (1 - x[c] - x[cp]))) for cp in range(k)] for c in range(k)]
        A1 = [[(Fr(1) if c == cp else Fr(0)) + N1[c][cp] / 2 for cp in range(k)] for c in range(k)]
        good2 = True
        for r in range(k):
            rows = [next(i for i in P if st[i][0][1] == r and st[i][1][c] == 1) for c in range(k)]
            for jq, Pq in enumerate(pairs):
                lhs = [sum(A1[c][cp] * M[rows[cp]][Q[jq]] for cp in range(k)) for c in range(k)]
                e = [Fr(0)] * k
                if r in Pq:
                    cstar = Pq[0] if Pq[1] == r else Pq[1]
                    e[cstar] = Fr(1)
                Te = {}
                if r in Pq:
                    for u in range(k):
                        if u != cstar:
                            pr = tuple(sorted((cstar, u)))
                            Te[pr] = Te.get(pr, 0) + w[u] / 2
                else:
                    Te[Pq] = Te.get(Pq, 0) + 2 * w[r]
                v = [Fr(0)] * k
                for pr, cf in Te.items():
                    c2 = 2 * (1 - x[pr[0]] - x[pr[1]])
                    v[pr[1]] += w[pr[0]] * cf / c2
                    v[pr[0]] += w[pr[1]] * cf / c2
                rhs = [-sum(A1[c][cp] * e[cp] for cp in range(k)) + v[c] for c in range(k)]
                good2 &= lhs == rhs
        ok &= good1 and good2
        print(f'k = {k}: (i) vacuum block {good1}, (ii) transformed coupling {good2}', flush=True)
print('PASS' if ok else 'FAIL')
