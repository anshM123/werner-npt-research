"""Comparison of the floating-point builder l1float.py with the exact builder j1lib.py (Lemma B), entry by entry, on
the letter-1 block of (R_inf)_j at general rational weights; also the known kernel p (x) H_j.  python validate.py K J"""
import sys, itertools
import numpy as np
from fractions import Fraction as Fr
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "middle_j1"))
from j1lib import J1
from l1float import L1Block
k, j = int(sys.argv[1]), int(sys.argv[2])
w = [Fr(1), Fr(2, 3), Fr(1, 2), Fr(3, 2), Fr(1), Fr(5, 4), Fr(1, 3)][:k]
B = L1Block([float(x * x) for x in w], j)
Y = B.build()
J = J1(w)
def tup(y):
    return tuple(1 if t in y else 0 for t in range(k))
states = [((1, a), tup(y)) if a is not None else (None, tup(y)) for (a, y) in B.states]
M = J.matrix(states)
err = max(abs(float(M[i][l]) - Y[i, l]) for i in range(len(states)) for l in range(len(states)))
K = B.kernel()
print(f'k={k} j={j}: {len(states)} states, max |exact - float| = {err:.2e}, |Y K| = {np.abs(Y @ K).max():.2e}, '
      f'min eig = {np.linalg.eigvalsh(Y).min():.3e}')
