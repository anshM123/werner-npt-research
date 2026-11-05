"""Exploration (floats): the Schur complement G on the T-blocks for u = 2 versus u = 1 at the same S = s_U."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI

def Gmat(w, T):
    ti = TypeI(w, T, F=float)
    M = np.array(ti.build()); n = ti.n
    U = ti.U; r0 = U[0]
    keep_T = [ti.idx[('R', t)] for t in T]
    tb = []
    for t in T:
        tb.append(ti.idx[('R', t)])
        for h in U:
            tb.append(ti.idx[('Z', t, h)])
    el = [i for i in range(n) if i not in tb and not (ti.labels[i][0] == 'P' and ti.labels[i][1] == r0)]
    A11 = M[np.ix_(el, el)]; A12 = M[np.ix_(el, tb)]; A22 = M[np.ix_(tb, tb)]
    G = A22 - A12.T @ np.linalg.solve(A11, A12)
    return (G + G.T) / 2, U

if __name__ == '__main__':
  np.set_printoptions(precision=6, suppress=True, linewidth=150)
  # T = {0,1,2}, weights s0, s1, rest; U weights vary with fixed S
  st = [0.3, 0.5, 0.2]
  for sU in ([1.0], [0.5, 0.5], [0.9, 0.1], [0.99, 0.01], [1/3, 1/3, 1/3]):
      w = np.sqrt(np.array(st + sU))
      G, U = Gmat(w, [0, 1, 2])
      u = len(U); d = u + 1
      wU = np.sqrt(np.array(sU)); wh = wU / np.linalg.norm(wU)
      # coordinates (rho, zeta aligned, zeta perp...)
      Bm = np.linalg.svd(wh[None, :])[2][1:].T if u > 1 else np.zeros((1, 0))
      P = np.zeros((d, d)); P[0, 0] = 1; P[1:, 1] = wh
      if u > 1: P[1:, 2:] = Bm
      Gc = np.zeros_like(G)
      for i in range(3):
          for j in range(3):
              Gc[i*d:(i+1)*d, j*d:(j+1)*d] = P.T @ G[i*d:(i+1)*d, j*d:(j+1)*d] @ P
      print('sU', sU)
      print('G00 aligned 2x2:\n', Gc[0:2, 0:2], '\n G01 aligned 2x2:\n', Gc[0:2, d:d+2])
      if u > 1:
          print(' G00 perp block:\n', Gc[0:d, 2:d].T)
          print(' G01 perp rows:\n', Gc[0:d, d+2:2*d].T)
