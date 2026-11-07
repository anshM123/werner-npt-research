"""Floating-point check of the margin of Theorem C' (NOTES.md, Section 4): lambda_min of F(X) = min_q (block of
Theorem C) against sigma/(12(sigma + S)) ||X||_F^2, at unnormalised random weights, and of the symmetric half-shell
H_c,+ >= 1/2 (free squared mass 1).  No claim rests on it."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from freeblock import free_block, q_eliminated
rng = np.random.default_rng(3)
worst = 9e9
tdev = 0.0
for it in range(150):
    m = int(rng.integers(1, 4)); u = int(rng.integers(2, 7)); mode = it % 3
    sT = np.exp(rng.normal(0, 2, m)); sU = np.exp(rng.normal(0, 2, u))
    if mode == 1: sT = sT * 1e-3
    if mode == 2: sT = sT * 1e3
    Fx = np.array(q_eliminated(free_block(list(np.sqrt(sT)), list(np.sqrt(sU)), float)))
    Fx = (Fx + Fx.T) / 2
    lam = np.linalg.eigvalsh(Fx)[0]
    X = rng.normal(size=(u, u)); v, vt = X.reshape(-1), X.T.reshape(-1)
    tdev = max(tdev, abs(v @ Fx @ v - vt @ Fx @ vt) / (v @ v))
    sig, S = sT.sum(), sU.sum()
    worst = min(worst, lam / (sig / (12 * (sig + S))))
print(f'150 weight vectors: min lambda_min(F) / (sigma/(12(sigma+S))) = {worst:.4f} (the theorem needs >= 1); '
      f'largest |F(X) - F(X^T)|/|X|^2 = {tdev:.1e}')
print('PASS' if worst >= 1 and tdev < 1e-9 else 'FAIL')
