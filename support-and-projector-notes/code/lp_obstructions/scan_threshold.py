import sys, time
from exact_lp import lp
n0, n1 = int(sys.argv[1]), int(sys.argv[2])
for n in range(n0, n1 + 1):
    t0 = time.time(); out = []; tstar = None
    t = max(1, n // 6)
    # find least t with positive optimum, scanning upward from a guess and downward if needed
    vals = {}
    def val(t):
        if t not in vals:
            st, q = lp(n, t); vals[t] = (st, q)
        return vals[t]
    t = 1
    while t <= n // 2:
        st, q = val(t)
        if st != 'optimal':
            break
        if q > 0:
            tstar = t; break
        t += 1
    desc = ' '.join(f'{k}:{(float(v[1]) if v[0]=="optimal" else v[0]):+.4g}' if v[0]=='optimal' else f'{k}:{v[0]}' for k, v in sorted(vals.items()))
    print(f'n={n:2d} t*={tstar}  [{desc}]  ({time.time()-t0:.0f}s)', flush=True)
