import sys, time
from fractions import Fraction as Fr
import tprove as P
for s in sys.argv[1:]:
    case, mode, fam = P.parse_fam(s)
    t0 = time.time()
    res, dt, info = P.prove(case, mode, fam, verbose=True)
    print(s, res, f'{dt:.1f}s', info.get('form_time'), info.get('sanity_points'), flush=True)
