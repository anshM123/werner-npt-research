import numpy as np, time
from bos import *
for n in [2,3,4,5,6,7]:
    for a in [a_n(n)-0.002, a_n(n)+0.002, 0.75]:
        t0=time.time()
        K = build_K(n, a); N = len(mindices(n))
        xw, yw = witness_sym(n)
        val, x, y = alt_min(K, N, starts=30, iters=400, rng=n, x0=xw)
        wv = 2.0**(-n)*((2*a-0.5)**n-0.5); pv = (a-0.5)**n
        print(f"n={n} a={a:.6f} min={val:.12e} witness={wv:.12e} product={pv:.12e} diff={val-min(wv,pv):.2e} t={time.time()-t0:.1f}s", flush=True)
