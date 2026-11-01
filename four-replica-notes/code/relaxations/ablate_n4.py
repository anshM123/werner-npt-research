"""Which constraints does the letter-resolved relaxation need at n = 4 for projectors (t = +1) and opposite sign (t = -1)?
usage: python ablate_n4.py CONS   (comma list from ORTH, MARG, PT, COV, LE2)"""
import sys, os, time
from alg4 import build_algebras
from herm4 import Herm4
cons = [c for c in sys.argv[1].split(',') if c]
algs = build_algebras(); t0 = time.time()
rel = Herm4(4, algs, cons, verbose=False)
print(f"n=4 cons={'+'.join(cons) or 'none'}: build {time.time() - t0:.0f}s, N = {rel.N}", flush=True)
for t in (1.0, -1.0, 0.5):
    sol = rel.solve(1.0, t, max_iter=400)
    print(f"n=4 cons={'+'.join(cons) or 'none'} t={t:+.1f}: status {sol.status}, value {sol.obj_val:+.6e} (true min {(1 - abs(t)) ** 2 / 16:.6e}); solve {rel.solve_time:.0f}s", flush=True)
