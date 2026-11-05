"""Exploration: Ghat for u >= 2 against the merged u = 1 system (one U-copy of weight sqrt(S))."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numG import Ghat

np.set_printoptions(precision=5, suppress=True, linewidth=160)
rng = np.random.default_rng(1)
for sT, sU in [([0.3, 0.5, 0.2], [0.4, 0.35, 0.25]), ([0.3, 0.5, 0.2], [1/3, 1/3, 1/3]), ([0.3, 0.5, 0.2], [0.9, 0.05, 0.05]),
               ([0.05, 0.02, 0.01], [0.5, 0.3, 0.2]), ([3.0, 1.0, 2.0], [0.5, 0.3, 0.2]), ([0.3, 0.5, 0.2], [0.4, 0.3, 0.2, 0.1])]:
    u = len(sU); d = u + 1; m = len(sT)
    Gh = Ghat(sT, sU)
    G1 = Ghat(sT, [sum(sU)])
    # (rho, z) part
    idx = [t * d + i for t in range(m) for i in (0, 1)]
    A = Gh[np.ix_(idx, idx)]
    print('sT', sT, 'sU', sU)
    print('  max |Ghat(rho,z) - G1hat| =', np.abs(A - G1).max(), ' (G1 max entry', np.abs(G1).max(), ')')
    perp = [t * d + i for t in range(m) for i in range(2, d)]
    if perp:
        B = Gh[np.ix_(perp, perp)]
        C = Gh[np.ix_(idx, perp)]
        print('  perp block eig range', np.linalg.eigvalsh(B).min(), np.linalg.eigvalsh(B).max(), ' max |coupling| =', np.abs(C).max())
        print('  perp diag blocks:', [np.round(np.diag(B)[i], 4) for i in range(len(perp))])
    print('  lambda_min Ghat', np.linalg.eigvalsh(Gh).min(), ' lambda_min G1hat', np.linalg.eigvalsh(G1).min())
    print('  diff matrix (rho,z):\n', A - G1)
