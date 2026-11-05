"""Exact check of a closed form closed_u{u}_m{m}.pkl (from closed_gen.py) against direct elimination.

At random rational weights w (so s = w^2 is rational) the type I block cf.TypeI(w, T = {0..m-1}) is built over Q
(Fractions), rescaled to x_{r,h} = P_{r,h}/(w_r w_h), z_{t,h} = zeta_{t,h}/(w_t w_h), the gauge rows P_{g,h} are deleted
for two different choices g in U (one if u = 1), and all other P and Q are eliminated exactly.  Both Schur complements must equal the
closed form evaluated at s.  Usage: python check_closed_gen.py u m [npoints] [dehom] (dehom: the closed form with s0 = 1)."""
import sys, os, time, pickle, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import flint
from cf import TypeI

u, m = int(sys.argv[1]), int(sys.argv[2])
npts = int(sys.argv[3]) if len(sys.argv) > 3 else 6
dehom = len(sys.argv) > 4 and sys.argv[4] == 'dehom'
D = pickle.load(open(f'closed_u{u}_m{m}{"_dehom" if dehom else ""}.pkl', 'rb'))
CTX = flint.fmpq_mpoly_ctx.get(tuple(D['names']), 'lex')
NUM = [[CTX.from_dict(n_) for (n_, d_) in row] for row in D['G']]
DEN = [[CTX.from_dict(d_) for (n_, d_) in row] for row in D['G']]
n = len(NUM)


def closed(s):
    q = [flint.fmpq(x.numerator, x.denominator) for x in s]
    out = []
    for i in range(n):
        r = []
        for j in range(n):
            v = NUM[i][j](*q) / DEN[i][j](*q)
            r.append(Fr(int(v.p), int(v.q)))
        out.append(r)
    return out


def direct(w, gauge_copy):
    T = list(range(m))
    ti = TypeI(w, T)
    M = ti.build()
    U = ti.U
    sc = [w[lab[1]] * w[lab[2]] if lab[0] in ('P', 'Z') else Fr(1) for lab in ti.labels]
    N = ti.n
    Ms = [[M[i][j] * sc[i] * sc[j] for j in range(N)] for i in range(N)]
    tb = []
    for t in T:
        tb.append(ti.idx[('R', t)])
        tb += [ti.idx[('Z', t, h)] for h in U]
    el = [i for i in range(N) if i not in tb and not (ti.labels[i][0] == 'P' and ti.labels[i][1] == gauge_copy)]
    order = el + tb
    A = [[Ms[i][j] for j in order] for i in order]
    ne = len(el)
    for p in range(ne):
        pv = A[p][p]
        assert pv != 0
        for i in range(p + 1, len(order)):
            if A[i][p] == 0:
                continue
            f = A[i][p] / pv
            for j in range(p + 1, len(order)):
                A[i][j] -= f * A[p][j]
    return [[A[ne + i][ne + j] for j in range(len(tb))] for i in range(len(tb))]


random.seed(100 * u + m)
ok = True
t0 = time.time()
cnt = 0
for trial in range(npts):
    w = [Fr(random.randint(1, 40), random.randint(1, 40)) for _ in range(u + m)]
    if dehom:
        w[0] = Fr(1)
    s = [x * x for x in w]
    Gc = closed(s[1:] if dehom else s)
    for g in list(range(m, m + u))[:2]:
        Gd = direct(w, g)
        good = all(Gc[i][j] == Gd[i][j] for i in range(n) for j in range(n))
        ok &= good
        cnt += 1
    print(f'w = {[str(x) for x in w]}: closed form = direct elimination (two gauges): {good}', flush=True)
print(f'u = {u}, m = {m}: {cnt} exact comparisons ({time.time() - t0:.1f} s)')
print('PASS' if ok else 'FAIL')
