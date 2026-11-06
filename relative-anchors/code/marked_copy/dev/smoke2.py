import sys, time
sys.path.insert(0, '.')
from fractions import Fraction as Fr
import mlib
from ufield import RF
t0 = time.time()
bad = 0; tot = 0
for k in (2, 3, 4, 5, 8, 12):
    for j in range(-1, k + 1):
        for i in range(0, (k - 1) // 2 + 1):
            C = mlib.Calc(k, j, i, sigma=None)
            names, norms, Q = mlib.reduced(C)
            # symmetric
            for a in range(len(names)):
                for b in range(len(names)):
                    if Q[a][b] != Q[b][a]: bad += 1; print('ASYM', k, j, i)
            for sig in (Fr(3, 7), Fr(5)):
                C2 = mlib.Calc(k, j, i, sigma=sig)
                n2, nr2, Q2 = mlib.reduced(C2)
                assert n2 == names
                for a in range(len(names)):
                    for b in range(len(names)):
                        tot += 1
                        if Q[a][b](sig) != Fr(int(Q2[a][b].p), int(Q2[a][b].q)):
                            bad += 1
            # kernel vectors
            sigma = RF(mlib.PL([0, 1]))
            for lab, v in mlib.kernel_vectors(names, k, j, i, sigma):
                for a in range(len(names)):
                    s = RF(0)
                    for nm, cf in v.items():
                        s = s + Q[a][names.index(nm)] * cf
                    if not s.is_zero():
                        bad += 1; print('KER', k, j, i, lab, names[a])
print('tot', tot, 'bad', bad, f'{time.time()-t0:.1f}s')
