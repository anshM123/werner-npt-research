"""Step 1 of the corner proof: det Q7 over QQ(q0, t, nu) at the interpolation nodes tau = 2..16 (exact, python-flint).
Each determinant is saved as (numerator, denominator) coefficient dictionaries in nodes/det_tau<T>.pkl.
Usage: python compute_nodes.py [tau_lo tau_hi]"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
import sys, time, os, pickle
import psutil
from nucalc import SymNu, form_in_basis, det_field, Q7NAMES

lo, hi = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (2, 16)
os.makedirs('nodes', exist_ok=True)
T0 = time.time()
for tau in range(lo, hi + 1):
    fn = f'nodes/det_tau{tau}.pkl'
    if os.path.exists(fn):
        print(f'tau = {tau}: exists', flush=True)
        continue
    t0 = time.time()
    be = SymNu(tau)
    Q7, G7 = form_in_basis(be, Q7NAMES)
    d = det_field(Q7, be.K)
    out = {'n': {e: (int(c.p), int(c.q)) for e, c in d.n.to_dict().items()},
           'd': {e: (int(c.p), int(c.q)) for e, c in d.d.to_dict().items()}, 'vars': be.syms, 'tau': tau}
    with open(fn + '.tmp', 'wb') as f:
        pickle.dump(out, f)
    os.replace(fn + '.tmp', fn)
    rss = psutil.Process(os.getpid()).memory_info().rss / 2 ** 30
    print(f'tau = {tau}: num {len(d.n)} terms (deg {d.n.degrees()}), den {len(d.d)} terms; {time.time()-t0:.1f}s; '
          f'rss {rss:.2f} GB', flush=True)
print(f'done {time.time()-T0:.0f}s')
