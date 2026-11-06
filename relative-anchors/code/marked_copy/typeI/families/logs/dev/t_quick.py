import time, sys
from fractions import Fraction as Fr
import tfamlib as T, tprove as P
import cProfile, pstats
fams = [('MT','anti',(4,4,2,2)), ('MT','sym',(0,1,0,0)), ('MT','sym',(4,4,2,2))]
for case, mode, fam in fams:
    t0 = time.time()
    C = T.FamCalc(case, mode, *fam)
    pr = cProfile.Profile(); pr.enable()
    names, Q = T.reduced(C)
    pr.disable()
    t1 = time.time()
    bad, npts = P.sanity(C, names, Q, Fr(1))
    print(case, mode, fam, len(names), f'{t1-t0:.2f}s', bad, npts, len(C.dens), len(C.stats['sectors']), C.stats['dmin'], C.stats['dmax'], flush=True)
    pstats.Stats(pr).sort_stats('cumulative').print_stats(12)
    sys.stdout.flush()
