"""Exploration (floating point): relax the hard part X' (X'w = 0, w^T X' = 0) of the U-block entrywise AFTER the exact
elimination of q, keep the soft part Omega = a w^T + w c^T, rho, zeta exact.  Compare the resulting Schur complement
on (rho, zeta) with the exact G, in the frame (rho sqrt(1+sig/S), merged sqrt(1+S/sig), perp)."""
import numpy as np, sys
from numlib import block

def frame(st, sU):
    u = len(sU); S = sU.sum(); wU = np.sqrt(sU); wh = wU/np.linalg.norm(wU)
    Pt = np.zeros((u+1,u+1)); Pt[0,0] = np.sqrt(1 + st/S); Pt[1:,1] = np.sqrt(1 + S/st)*wh
    if u > 1: Pt[1:,2:] = np.linalg.svd(wh[None,:])[2][1:].T
    return Pt

def analyse(st, sU, mode='tight'):
    u = len(sU)
    P, lab, ne = block(st, sU, scaled=False)       # order: X (u*u), q, rho, zeta
    nX = u*u
    # eliminate q exactly
    iq = nX
    keep = [i for i in range(P.shape[0]) if i != iq]
    Pq = P[np.ix_(keep, keep)] - np.outer(P[keep, iq], P[iq, keep]) / P[iq, iq]
    # now order: X (nX), rho, zeta (u)
    w = np.sqrt(sU); S = sU.sum()
    # orthonormal basis of hard space {X: Xw=0, w^T X=0}: Pi X Pi
    wh = w/np.linalg.norm(w)
    Bp = np.linalg.svd(wh[None,:])[2][1:].T        # u x (u-1), orthonormal basis of w-perp
    Hb = np.kron(Bp, Bp)                            # vec(X) row-major: X = Bp Theta Bp^T -> vec = kron(Bp,Bp) vec(Theta)
    # soft basis: a w^T with a = Bp alpha, and w c^T
    Sa = np.kron(Bp, w[:, None])                    # vec(a w^T) = kron(a, w)
    Sc = np.kron(w[:, None], np.eye(u))             # vec(w c^T) = kron(w, c)
    nr = Pq.shape[0] - nX
    T = np.zeros((Pq.shape[0], Hb.shape[1] + Sa.shape[1] + Sc.shape[1] + nr))
    T[:nX, :Hb.shape[1]] = Hb
    T[:nX, Hb.shape[1]:Hb.shape[1]+Sa.shape[1]] = Sa
    T[:nX, Hb.shape[1]+Sa.shape[1]:Hb.shape[1]+Sa.shape[1]+Sc.shape[1]] = Sc
    T[nX:, -nr:] = np.eye(nr)
    Pt = T.T @ Pq @ T
    nh = Hb.shape[1]
    Fh = Pt[:nh, :nh]; Cpl = Pt[:nh, nh:]; R = Pt[nh:, nh:]
    # exact G from this representation
    def schur_last(M, k):
        A = M[:-k,:-k]; B = M[:-k,-k:]; D = M[-k:,-k:]
        return D - B.T @ np.linalg.solve(A, B)
    Gex = schur_last(Pt, u+1)
    # relaxation: Fh(theta) >= sum lam_rh X_rh^2 with X = vec^-1(Hb theta)
    d = sU + S/2 + st; rho = sU/d; g = rho.sum()
    lam = 0.5 - np.outer(rho, rho)/(2+g)
    np.fill_diagonal(lam, 1.5 - rho + rho**2/(2+g))
    if mode == 'lemmaH':
        lam = np.full((u,u), 5/12); np.fill_diagonal(lam, 1.5 - rho)
    # coupling representative H = Hb @ Cpl (each column: X-space vector)
    Hx = Hb @ Cpl                                   # nX x (rest)
    Lrel = (Hx / lam.reshape(-1)[:, None]).T @ Hx   # sum H_i H_j / lam
    Rrel = R - Lrel
    # check validity: Fh - Hb^T diag(lam) Hb >= 0
    chk = np.linalg.eigvalsh(Fh - Hb.T @ np.diag(lam.reshape(-1)) @ Hb).min()
    evR = np.linalg.eigvalsh(Rrel).min()
    Grel = schur_last(Rrel, u+1) if evR > -1e300 else None
    Pf = frame(st, sU)
    return (np.linalg.eigvalsh(Pf.T @ Gex @ Pf).min(), np.linalg.eigvalsh(Pf.T @ Grel @ Pf).min(), chk,
            np.linalg.eigvalsh(Rrel[:-(u+1), :-(u+1)]).min())

if __name__ == "__main__":
  rng = np.random.default_rng(4)
  for mode in ('tight', 'lemmaH'):
      for u in (3, 4, 6, 8):
          res = []
          for it in range(40):
              fam = it % 4
              st = float(np.exp(rng.uniform(-8, 8)))
              if fam == 0: sU = np.exp(rng.uniform(-3, 3, u))
              elif fam == 1: sU = np.exp(rng.uniform(-8, 8, u))
              elif fam == 2: sU = np.ones(u); sU[0] = np.exp(rng.uniform(-10, 10))
              else: sU = np.exp(np.linspace(0, rng.uniform(-10, 10), u))
              res.append(analyse(st, sU, mode))
          res = np.array(res)
          print(mode, u, 'min eig exact', res[:,0].min(), ' relaxed', res[:,1].min(), ' validity', res[:,2].min(),
                ' relaxed soft block min eig', res[:,3].min(), flush=True)
