# T3 numerics: x in Sym^n, y arbitrary in (R^3)^{(x)n}.  min_y <y|P^{(x)n}(xx^T)|y> = lambda_min(P^{(x)n}(xx^T)).
import numpy as np, sys, time
from bos import *
def alt_T3(n, a, starts, rng, x0=None, iters=200):
    B = sym_basis_full(n); N = B.shape[1]
    best = (np.inf, None)
    for s in range(starts):
        if s == 0 and x0 is not None: xc = x0.copy()
        else: xc = rng.standard_normal(N)
        xc /= np.linalg.norm(xc)
        old = np.inf
        for it in range(iters):
            x = B @ xc
            M = apply_Pn(np.outer(x, x), n, a)
            w, V = np.linalg.eigh(M); y = V[:, 0]
            My = apply_Pn(np.outer(y, y), n, a)     # self-adjointness: <y|P(xx)|y> = <x|P(yy)|x>
            Ms = B.T @ My @ B
            w2, V2 = np.linalg.eigh(Ms); xc = V2[:, 0]
            if abs(old - w2[0]) < 1e-14: break
            old = w2[0]
        if w2[0] < best[0]: best = (w2[0], xc)
    return best
rng = np.random.default_rng(7)
for n in [int(s) for s in sys.argv[1:]]:
    xw, yw = witness_sym(n)
    for a in [a_n(n)-0.002, a_n(n)+0.002, 0.75]:
        t0 = time.time()
        val, xc = alt_T3(n, a, starts=12 if n<6 else 6, rng=rng, x0=xw)
        wv = 2.0**(-n)*((2*a-0.5)**n-0.5); pv = (a-0.5)**n
        print(f"n={n} a={a:.6f}: min over x sym, y arbitrary = {val:.10e}; min(witness, product) = {min(wv,pv):.10e}; t={time.time()-t0:.1f}s", flush=True)
