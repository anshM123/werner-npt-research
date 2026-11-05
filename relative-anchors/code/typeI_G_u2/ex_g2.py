"""Exploration: does G (T-blocks) depend on w_U through a single functional?  Match u=2 split to u=3 equal."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ex_g1 import Gmat
from scipy.optimize import brentq

def aligned(st, sU):
    w = np.sqrt(np.array(list(st) + list(sU)))
    G, U = Gmat(w, list(range(len(st))))
    u = len(U); d = u + 1
    wU = np.sqrt(np.array(sU)); wh = wU / np.linalg.norm(wU)
    P = np.zeros((d, 2)); P[0, 0] = 1; P[1:, 1] = wh
    out = []
    for i in range(len(st)):
        for j in range(len(st)):
            out.append(P.T @ G[i*d:(i+1)*d, j*d:(j+1)*d] @ P)
    return np.array(out)

for st in ([0.3, 0.5, 0.2], [0.05, 0.02, 0.01], [3.0, 1.0, 2.0]):
    ref = aligned(st, [0.8, 0.12, 0.08])
    f = lambda a: aligned(st, [a, 1 - a])[0][0, 0] - ref[0][0, 0]
    try:
        a = brentq(f, 0.55, 0.999999)
    except Exception as e:
        print('no root', e); continue
    A2 = aligned(st, [a, 1 - a])
    print(st, 'matching a', a, ' sum s^2 u2', a*a + (1-a)**2, ' u3', 1/3, ' max diff other entries', np.abs(A2 - ref).max())
