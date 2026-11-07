# SOS (decomposability) threshold of K_n(a) on Sym x Sym, using the Y-basis span (exact up to the span assumption) and the full O(3) SDP
import numpy as np, cvxpy as cp, sys
from math import factorial
from ybasis import build_Y
from o3red import O3Red
from bos import a_n, build_K
from sdp_o3b import dec_sdp
for n in [int(s) for s in sys.argv[1:]]:
    R = O3Red(n); G = R.gamma_matrix()
    def maxt(a):
        tv, st, Q = dec_sdp(R, G, build_K(n, a), solver='SCS')
        return tv
    lo, hi = a_n(n), 0.75
    # find smallest a with max t >= 0 (SOS-positive)
    for it in range(14):
        mid = (lo+hi)/2
        if maxt(mid) >= -1e-7: hi = mid
        else: lo = mid
    # crossover witness/product
    from scipy.optimize import brentq
    f = lambda a: 2.0**(-n)*((2*a-0.5)**n-0.5) - (a-0.5)**n
    ac = brentq(f, a_n(n), 0.75)
    print(f"n={n}: a_n={a_n(n):.6f}, crossover a_c={ac:.6f}, SOS threshold ~ {hi:.5f}", flush=True)
