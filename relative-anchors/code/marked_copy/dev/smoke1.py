import sys, time
sys.path.insert(0, '.')
from fractions import Fraction as Fr
import mlib, explicit
bad = 0; tot = 0
t0 = time.time()
for k in (2, 3, 4, 5):
    for w0 in (Fr(2), Fr(1, 3)):
        sig = w0 * w0
        w = [w0] + [Fr(1)] * (k - 1)
        for j in range(-1, k + 1):
            st, M = explicit.block(k, j, w)
            for i in range(0, (k - 1) // 2 + 1):
                C = mlib.Calc(k, j, i, sigma=sig)
                names, norms, Q = mlib.reduced(C)
                st2, vecs = explicit.twisted(k, j, i)
                en = sorted(vecs)
                if sorted(names) != en:
                    print('NAMES', k, j, i, names, en); bad += 1; continue
                for a, na in enumerate(names):
                    if sum(c * c for c in vecs[na].values()) != norms[a]:
                        print('NORM', k, j, i, na); bad += 1
                    for b, nb in enumerate(names):
                        ex = explicit.form(M, vecs[na], vecs[nb])
                        got = Fr(int(Q[a][b].p), int(Q[a][b].q))
                        tot += 1
                        if ex != got:
                            bad += 1
                            if bad < 15: print('DIFF', k, j, i, sig, na, nb, ex, got)
print('compared', tot, 'bad', bad, f'{time.time()-t0:.1f}s')
