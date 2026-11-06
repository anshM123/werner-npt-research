"""Floating-point exploration (no claim rests on it): sup of the (A) ratios and inf of the (B) ratios of the weighted
block Gershgorin criterion over two-valued weight vectors (squared weights s and 1 on U and on T), in the frame of
prove_gersh.py.  Prints the extreme configurations (used to place the negative-control witnesses)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ex_numframe import parts

if __name__ == '__main__':
    best, worst = {}, {}
    for u1, u2 in [(1, 2), (1, 3), (1, 5), (2, 2), (2, 3), (3, 3), (1, 8)]:
        for m1, m2 in [(2, 0), (1, 1), (0, 2), (2, 1), (1, 2), (3, 0)]:
            for ls in np.linspace(-6, 6, 25):
                s = float(np.exp(ls))
                sT = [s] * m1 + [1.0] * m2
                o = parts(sT, u1, s, u2, 1.0)
                for k, v in o.items():
                    if not v:
                        continue
                    if k.startswith('A'):
                        if max(v) > best.get(k, (-1,))[0]:
                            best[k] = (max(v), u1, u2, m1, m2, s)
                    else:
                        if min(v) < worst.get(k, (1e9,))[0]:
                            worst[k] = (min(v), u1, u2, m1, m2, s)
    for k in sorted(best):
        print('sup', k, best[k])
    for k in sorted(worst):
        print('inf', k, worst[k])
    print('exploration done')
