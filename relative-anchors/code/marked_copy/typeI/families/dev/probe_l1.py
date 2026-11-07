"""Development probe (not part of any proof): at a member of case MU sym, compare the rows without the T' letter
(O, U, W, G after the split) with the letter-1 two-row form of mlib.py at k' = k - tau', j' = j - tau' (exact, one
value of sigma).  Prints whether X - L0 is >= 0 and the inertia of M' = Q2 - (L0 (+) 0)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import flint
import tlib
import mlib
import tsplit


def inertia(A):
    """exact inertia of a symmetric fmpq matrix (list of lists) by symmetric elimination with 2x2 pivots avoided:
    returns (pos, neg, zero) via LDL with diagonal pivoting on the largest |diagonal|; falls back on rank"""
    n = len(A)
    M = flint.fmpq_mat(n, n, [flint.fmpq(A[a][b].numerator, A[a][b].denominator) for a in range(n) for b in range(n)])
    # characteristic polynomial sign changes (Descartes on charpoly of a symmetric matrix is exact)
    cp = M.charpoly()
    co = [cp[i] for i in range(n + 1)]   # co[i] coefficient of x^i
    zero = 0
    while zero <= n and co[zero] == 0:
        zero += 1
    rest = co[zero:]
    # number of positive roots = sign changes of p(x); negative roots = sign changes of p(-x)
    def changes(c):
        s = [x for x in c if x != 0]
        return sum(1 for a, b in zip(s, s[1:]) if (a > 0) != (b > 0))
    pos = changes(rest)
    neg = changes([c * (-1) ** i for i, c in enumerate(rest)])
    return pos, neg, zero


def main(q0, t, i, tp, sg):
    sg = Fr(sg)
    k = tp + q0 + t + 2 * i + 1
    j, tau = q0 + tp + i, tp
    C = tlib.Calc(k, j, tau, i, 'MU', 'sym', sigma=sg)
    names, _, Q = tlib.reduced(C)
    m = q0 + t
    vecs, two, three, split = tsplit.basis('MU', names, q0, m, i, Fr(1))
    order = two + three
    idx = {nm: a for a, nm in enumerate(names)}
    def ent(p, q):
        s = Fr(0)
        for oa, ca in vecs[p].items():
            for ob, cb in vecs[q].items():
                v = Q[idx[oa]][idx[ob]]
                s += ca * cb * Fr(int(v.p), int(v.q))
        return s
    Xn = [nm for nm in two if nm[0] in 'OUWG' or nm[0] == 'B']
    Yn = [nm for nm in two if nm[0] in 'CD']
    print('member k, j, tau, i =', (k, j, tau, i), 'names', len(names), 'two', two, 'three', three)
    # letter-1 form
    kp, jp = k - tp, j - tp
    Cm = mlib.Calc(kp, jp, i, sigma=sg)
    n1, _, Q1 = mlib.reduced(Cm)
    v1, two1, three1, split1 = tsplit.basis('MU', n1, q0, m, i, Fr(1))
    idx1 = {nm: a for a, nm in enumerate(n1)}
    def ent1(p, q):
        s = Fr(0)
        for oa, ca in v1[p].items():
            for ob, cb in v1[q].items():
                v = Q1[idx1[oa]][idx1[ob]]
                s += ca * cb * Fr(int(v.p), int(v.q))
        return s
    assert set(two1) == set(Xn), (two1, Xn)
    X = [[ent(p, q) for q in Xn] for p in Xn]
    L0 = [[ent1(p, q) for q in Xn] for p in Xn]
    D = [[X[a][b] - L0[a][b] for b in range(len(Xn))] for a in range(len(Xn))]
    print('X - L0 inertia (pos, neg, zero):', inertia(D))
    print('L0 inertia:', inertia(L0), ' X inertia:', inertia(X))
    allv = Xn + Yn
    Q2 = [[ent(p, q) for q in allv] for p in allv]
    Mp = [[Q2[a][b] - (L0[a][b] if (a < len(Xn) and b < len(Xn)) else 0) for b in range(len(allv))]
          for a in range(len(allv))]
    print('Q2 inertia:', inertia(Q2), " M' = Q2 - (L0 + 0) inertia:", inertia(Mp))


if __name__ == '__main__':
    a = sys.argv[1:]
    main(int(a[0]), int(a[1]), int(a[2]), int(a[3]), a[4])
