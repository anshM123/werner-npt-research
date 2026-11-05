"""Does the Schur complement onto the U-variables (X_U, q) depend on the T-weights only through s_T?"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI
np.set_printoptions(linewidth=200, precision=5, suppress=True)
def ured(w, T):
    ti = TypeI(w, T, F=float); M = np.array(ti.build()); U = ti.U
    keep = [ti.idx[('P', r, h)] for r in U for h in U] + [ti.idx[('Q',)]]
    rest = [i for i in range(ti.n) if i not in keep]
    S = M[np.ix_(keep, keep)] - M[np.ix_(keep, rest)] @ np.linalg.solve(M[np.ix_(rest, rest)], M[np.ix_(rest, keep)])
    return S
wU = [1.0, 0.7, 1.3]
for wT in ([0.5], [0.3, 0.4], [0.1, 0.2, np.sqrt(0.25 - 0.05)], [np.sqrt(0.125), np.sqrt(0.125)]):
    w = list(wT) + wU; T = list(range(len(wT)))
    S = ured(w, T)
    print(len(wT), np.sum(np.square(wT)), np.round(np.linalg.eigvalsh(S), 6))
