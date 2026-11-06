"""Exact witness for Step 1 (closed forms through a few symmetric quantities): the diagonal block G_tt is not a
function of (s_t, sigma, S, sum_h s_h^2).  The U-amplitudes (1, 9, 10) and (5, 6, 11) give squared weights with the same
first and second power sums (182 and 16562; a, b, a + b with a^2 + ab + b^2 = 91 in both cases), and the exact G
(Lemma B builder, whole block, gauge P_{r0,h}) differs in its rho-rho entry, for m = 1 and m = 2.
(For u = 3 the multiset of U-weights is fixed by three power sums, so this is all that can be said at u = 3; the
structural reason for general u is the pair structure of Theorem C's form, typeI_G_m1/NOTES.md Section 4.)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import flint
from emlib import lemmab_block, ldl_pivots, qf, Q


def G_rr(wT, wU):
    m = len(wT)
    w = [Fr(x) for x in list(wT) + list(wU)]
    k = len(w)
    T = list(range(m))
    U = list(range(m, k))
    Mb, lab = lemmab_block(w, T)
    n = len(lab)
    idx = {l: i for i, l in enumerate(lab)}
    keep = [i for i in range(n) if not (lab[i][0] == 'P' and lab[i][1] == U[0])]
    tb = []
    for t in T:
        tb.append(idx[('R', t)])
        for h in U:
            tb.append(idx[('Z', t, h)])
    el = [i for i in keep if i not in tb]
    C = flint.fmpq_mat([[qf(Mb[i][j]) for j in el] for i in el])
    B = flint.fmpq_mat([[qf(Mb[i][j]) for j in tb] for i in el])
    A = flint.fmpq_mat([[qf(Mb[i][j]) for j in tb] for i in tb])
    G = A - B.transpose() * C.solve(B)
    return G[0, 0]


U1, U2 = (1, 9, 10), (5, 6, 11)
p = lambda U, e: sum(x ** (2 * e) for x in U)
ok = p(U1, 1) == p(U2, 1) and p(U1, 2) == p(U2, 2) and p(U1, 3) != p(U2, 3)
print('power sums of s_U: %s vs %s' % ([p(U1, e) for e in (1, 2, 3)], [p(U2, e) for e in (1, 2, 3)]))
for wT in ((2,), (2, 3)):
    g1, g2 = G_rr(wT, U1), G_rr(wT, U2)
    print('w_T = %s: G[rho_t, rho_t] = %s (U-amplitudes %s) vs %s (%s); equal: %s' % (wT, g1, U1, g2, U2, g1 == g2))
    ok &= g1 != g2
print('PASS' if ok else 'FAIL')
