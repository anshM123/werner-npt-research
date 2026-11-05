import os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from djlib import block
for j in (1, 2, 3, 5):
    for N in (0, 1, 3):
        for big in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6):
            w = [float(big)] + [1.0 + 0.1 * i for i in range(j)] + [0.5 + 0.2 * i for i in range(N)]
            S = sum(w); w = [v / S for v in w]
            M = np.array(block(w[:j + 1], w[j + 1:], bscale=1.001), float)
            ev = np.linalg.eigvalsh(M).min()
            if ev < 0:
                print(j, N, big, ev); break
        else:
            print(j, N, 'none up to 1e6')
