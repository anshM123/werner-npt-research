import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import flint
from cert_row import build_P
t0 = time.time()
ctx, P = build_P()
print('built', time.time() - t0, flush=True)
Q = flint.fmpq
vals = [Q(3, 17), Q(5, 19), Q(2, 23), Q(7, 29), Q(1, 31), Q(4, 37), Q(6, 41), Q(8, 43)]
t1 = time.time()
for z in range(5):
    v = P(*vals)
print('5 evaluations', time.time() - t1, v > 0)
