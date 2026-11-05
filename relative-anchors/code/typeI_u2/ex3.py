import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ex2 import reduced, comparison
np.set_printoptions(linewidth=200, precision=3, suppress=False)
for eps in (1e-2, 1e-4, 1e-6):
    w = np.array([1.0, 0.6, 1.2, 0.8, 1.1]); T = [0, 1]
    w[T] *= np.sqrt(eps)
    S, blocks, labs = reduced(w, T)
    names, C, lmin = comparison(S, blocks)
    print("eps", eps, names, "lmin full", lmin)
    print(C)
    # T-block internal eigenvalues
    for n in names:
        B = np.array(blocks[n]).T
        print(n, np.linalg.eigvalsh(B.T @ S @ B))
