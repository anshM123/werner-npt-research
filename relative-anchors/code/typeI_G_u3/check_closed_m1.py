"""Exact check of the closed form closed_m1_u{u}.pkl (default u = 3) against direct elimination.

At random rational weights w (so s = w^2 is rational), the type I block cf.TypeI(w, T = {0}) is built over Q
(Fractions), rescaled to the coordinates x_{r,h} = P_{r,h}/(w_r w_h), z_h = zeta_{0,h}/(w_0 w_h), the gauge rows
P_{1,h} are deleted and every other P and Q are eliminated exactly.  The Schur complement must equal the closed form
evaluated at s.  Also checks the gauge independence: deleting the rows P_{2,h} instead gives the same G."""
import sys, os, time, pickle, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cf import TypeI

u = int(sys.argv[1]) if len(sys.argv) > 1 else 3
D = pickle.load(open(f'closed_m1_u{u}.pkl', 'rb'))


import flint
CTX = flint.fmpq_mpoly_ctx.get(tuple(D['names']), 'lex')
NUM = [[CTX.from_dict(n_) for (n_, d_) in row] for row in D['G']]
DEN = [[CTX.from_dict(d_) for (n_, d_) in row] for row in D['G']]


def closed(s):
    q = [flint.fmpq(x.numerator, x.denominator) for x in s]
    out = []
    for i in range(len(NUM)):
        r = []
        for j in range(len(NUM)):
            v = NUM[i][j](*q) / DEN[i][j](*q)
            r.append(Fr(int(v.p), int(v.q)))
        out.append(r)
    return out


def direct(w, gauge_copy):
    T = [0]
    ti = TypeI(w, T)
    M = ti.build()
    U = ti.U
    sc = []
    for lab in ti.labels:
        sc.append(w[lab[1]] * w[lab[2]] if lab[0] in ('P', 'Z') else Fr(1))
    n = ti.n
    Ms = [[M[i][j] * sc[i] * sc[j] for j in range(n)] for i in range(n)]
    tb = [ti.idx[('R', 0)]] + [ti.idx[('Z', 0, h)] for h in U]
    el = [i for i in range(n) if i not in tb and not (ti.labels[i][0] == 'P' and ti.labels[i][1] == gauge_copy)]
    order = el + tb
    A = [[Ms[i][j] for j in order] for i in order]
    ne = len(el)
    N = len(order)
    for p in range(ne):
        pv = A[p][p]
        assert pv != 0
        for i in range(p + 1, N):
            if A[i][p] == 0:
                continue
            f = A[i][p] / pv
            for j in range(p + 1, N):
                A[i][j] -= f * A[p][j]
    return [[A[ne + i][ne + j] for j in range(len(tb))] for i in range(len(tb))]


random.seed(11)
ok = True
t0 = time.time()
cnt = 0
for trial in range(12):
    w = [Fr(random.randint(1, 40), random.randint(1, 40)) for _ in range(u + 1)]
    s = [x * x for x in w]
    Gc = closed(s)
    for gc in (1, 2):
        Gd = direct(w, gc)
        good = all(Gc[i][j] == Gd[i][j] for i in range(u + 1) for j in range(u + 1))
        ok &= good
        cnt += 1
    print(f'w = {[str(x) for x in w]}: closed form = direct elimination (gauge rows P_1h and P_2h): {good}', flush=True)
print(f'{cnt} exact comparisons ({time.time() - t0:.1f} s)')
print('PASS' if ok else 'FAIL')
