"""Steps (0) and (2) of prove_letter2.py for the generic sym family [q0 >= 2, t >= 4, i >= 2, tau >= 1], whose
elimination is done in prove_letter2_split.py: the Newton window is exact on the family (level offsets), and the
symbolic Q, G agree exactly with the numeric orbit calculus Num2 at 3^4 = 81 sample points (k up to 30)."""
import itertools, time
from math import comb
from orbit2 import SymF2, Num2, reduced2
import prove_letter2 as PL
t0 = time.time()
be = SymF2()
SH = PL.shifts(None, None, None, None, 'sym')
names, G, Q = reduced2(be, 'sym')
dmin = min(d for d, dp in be.seen); dpmin = min(dp for d, dp in be.seen)
okw = (dmin + SH['q0'] >= 2) and (dpmin + SH['t'] >= 2) and (SH['q0'] + SH['t'] >= 2) and SH['i'] >= 2
print('level offsets: delta >=', dmin, ', delta\' >=', dpmin, '-> window exact on the family:', okw)
ok = True; npts = 0
for sample in itertools.product(*[(SH[s], SH[s] + 1, SH[s] + 3) for s in be.syms]):
    val = dict(zip(be.syms, sample))
    q0v, tv, iv, tauv = val['q0'], val['t'], val['i'], val['tau']
    m = q0v + tv; k = m + 2 * iv + tauv; j = q0v + iv + tauv
    n2, G2, Q2 = reduced2(Num2(k, j, tauv, iv), 'sym')
    N0 = 2 ** iv * comb(m, q0v)
    npts += 1
    if n2 != names:
        ok = False; print('names differ at', k, j, tauv, iv); continue
    for a in range(len(names)):
        for b in range(len(names)):
            if PL.evalF(Q[a][b], sample) * N0 != Q2[a][b] or PL.evalF(G[a][b], sample) * N0 != G2[a][b]:
                ok = False; print('value mismatch at', k, j, tauv, iv, names[a], names[b])
print(npts, 'sample points;', 'PASS' if (ok and okw) else 'FAIL', f'{time.time()-t0:.0f}s')
