"""Exploration (floating point): the same entrywise relaxation of the hard part X', but BEFORE eliminating q (q kept
as a variable, as in typeI_u2/NOTES.md Section B5 for the whole block).  Valid as a lower bound since
Phi(X', 0, 0, 0) = F(X') + (tr X')^2/2 >= F(X') >= sum lam_rh X'_rh^2 (Lemma H).  Compares the smallest eigenvalue
of the relaxed G in the frame of NOTES.md Section 1 with the exact one, and reports the smallest eigenvalue of the
relaxed form itself (normalised by its largest diagonal entry)."""
import numpy as np
from numlib import block
from relax1 import frame


def analyse(st, sU):
    u = len(sU)
    P, lab, ne = block(st, sU, scaled=False)       # order: X (u*u), q, rho, zeta
    nX = u * u
    w = np.sqrt(sU)
    wh = w / np.linalg.norm(w)
    Bp = np.linalg.svd(wh[None, :])[2][1:].T
    Hb = np.kron(Bp, Bp)
    Sa = np.kron(Bp, w[:, None])
    Sc = np.kron(w[:, None], np.eye(u))
    nr = P.shape[0] - nX                            # q, rho, zeta
    nh = Hb.shape[1]
    T = np.zeros((P.shape[0], nh + Sa.shape[1] + Sc.shape[1] + nr))
    T[:nX, :nh] = Hb
    T[:nX, nh:nh + Sa.shape[1]] = Sa
    T[:nX, nh + Sa.shape[1]:nh + Sa.shape[1] + Sc.shape[1]] = Sc
    T[nX:, -nr:] = np.eye(nr)
    Pt = T.T @ P @ T
    Cpl = Pt[:nh, nh:]
    R = Pt[nh:, nh:]
    d = sU + sU.sum() / 2 + st
    rho = sU / d
    lam = np.full((u, u), 5 / 12)
    np.fill_diagonal(lam, 1.5 - rho)
    Hx = Hb @ Cpl
    Rrel = R - (Hx / lam.reshape(-1)[:, None]).T @ Hx
    k = u + 1

    def schur_last(M, k):
        A = M[:-k, :-k]; B = M[:-k, -k:]; D = M[-k:, -k:]
        return D - B.T @ np.linalg.solve(A, B)
    C = P[:ne, :ne]; B = P[:ne, ne:]; A = P[ne:, ne:]
    Gex = A - B.T @ np.linalg.solve(C, B)
    Pf = frame(st, sU)
    ev_rel_form = np.linalg.eigvalsh(Rrel).min() / np.abs(np.diag(Rrel)).max()
    try:
        Grel = schur_last(Rrel, k)
        evg = np.linalg.eigvalsh(Pf.T @ Grel @ Pf).min()
    except np.linalg.LinAlgError:
        evg = float('nan')
    return np.linalg.eigvalsh(Pf.T @ Gex @ Pf).min(), evg, ev_rel_form


if __name__ == '__main__':
    rng = np.random.default_rng(4)
    for u in (2, 3, 4, 6):
        res = []
        for it in range(24):
            fam = it % 4
            st = float(np.exp(rng.uniform(-8, 8)))
            if fam == 0: sU = np.exp(rng.uniform(-3, 3, u))
            elif fam == 1: sU = np.exp(rng.uniform(-8, 8, u))
            elif fam == 2: sU = np.ones(u); sU[0] = np.exp(rng.uniform(-10, 10))
            else: sU = np.exp(np.linspace(0, rng.uniform(-10, 10), u))
            res.append(analyse(st, sU))
        res = np.array(res)
        print(f'u = {u}: exact {res[:, 0].min():.5f}; relaxed before eliminating q: smallest eigenvalue of G_rel in '
              f'the frame {np.nanmin(res[:, 1]):.5f}, of the relaxed form (normalised) {res[:, 2].min():.3e}', flush=True)
