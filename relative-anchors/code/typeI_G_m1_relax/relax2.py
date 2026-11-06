"""Exploration (floating point): the relaxed form Psi(a, c, rho, zeta) from the explicit formulas, for several choices
of representatives of the hard-part couplings; G_rel = its Schur complement on (rho, zeta) vs the exact G."""
import numpy as np
from numlib import block
from forms import Forms
from relax1 import frame

def exactG(st, sU):
    P, lab, ne = block(st, sU, scaled=False)
    C = P[:ne,:ne]; B = P[:ne,ne:]; A = P[ne:,ne:]
    G = A - B.T @ np.linalg.solve(C, B); return (G+G.T)/2

def Psi_matrix(st, sU, var, lamkind='H'):
    Fm = Forms(st, sU); u = Fm.u; w = Fm.w; S = Fm.S
    Pi = np.eye(u) - np.outer(w, w) / S
    Bp = np.linalg.svd((w/np.linalg.norm(w))[None,:])[2][1:].T
    rho_ = Fm.sU / Fm.d; g = rho_.sum()
    if lamkind == 'H':
        lam = np.full((u,u), 5/12); np.fill_diagonal(lam, 1.5 - rho_)
    else:
        lam = 0.5 - np.outer(rho_, rho_)/(2+g); np.fill_diagonal(lam, 1.5 - rho_ + rho_**2/(2+g))
    n = (u-1) + u + 1 + u
    def unpack(v):
        a = Bp @ v[:u-1]; c = v[u-1:2*u-1]; rho = v[2*u-1]; zeta = v[2*u:]
        return a, c, rho, zeta
    def Q(v):
        a, c, rho, zeta = unpack(v)
        Om = np.outer(a, w) + np.outer(w, c)
        Hf = Fm.Hfull(rho, zeta)
        Gc = Fm.Gcal(a, c)
        if var in ('V1', 'V2'):
            Hr = Pi @ Hf @ Pi
        elif var == 'V3':
            # drop the w-outer-product terms (they vanish on hard X')
            Wd = Fm.W; Y0 = Fm.wt*rho + w@zeta; cc = Y0/2 + w*zeta
            Hr = (rho - (Fm.kappa*rho + Fm.wt*(w@zeta)/(st+S))) * np.eye(u) + Fm.wt * (cc[:,None] * Wd[1:,1:].T)
        else:
            Hr = Hf
        if var == 'V2':
            Gc = Pi @ Gc @ Pi
        R = st * Gc + Hr
        return st * Fm.f(a, c) + np.sum(Hf * Om) + Fm.Acal(rho, zeta) - 0.25 * np.sum(R**2 / lam)
    M = np.zeros((n, n)); E = np.eye(n)
    dg = [Q(E[i]) for i in range(n)]
    for i in range(n):
        M[i,i] = dg[i]
        for j in range(i+1, n):
            M[i,j] = M[j,i] = (Q(E[i]+E[j]) - dg[i] - dg[j]) / 2
    return M

if __name__ == '__main__':
    rng = np.random.default_rng(4)
    for u in (3, 5, 8):
        out = {}
        for it in range(24):
            fam = it % 4
            st = float(np.exp(rng.uniform(-8, 8)))
            if fam == 0: sU = np.exp(rng.uniform(-3, 3, u))
            elif fam == 1: sU = np.exp(rng.uniform(-8, 8, u))
            elif fam == 2: sU = np.ones(u); sU[0] = np.exp(rng.uniform(-10, 10))
            else: sU = np.exp(np.linspace(0, rng.uniform(-10, 10), u))
            Pf = frame(st, sU)
            ge = np.linalg.eigvalsh(Pf.T @ exactG(st, sU) @ Pf).min()
            out.setdefault('exact', []).append(ge)
            for var in ('V1', 'V2', 'V3', 'V4'):
                M = Psi_matrix(st, sU, var)
                k = u + 1
                A = M[:-k,:-k]; B = M[:-k,-k:]; D = M[-k:,-k:]
                evA = np.linalg.eigvalsh(A).min()
                Gr = D - B.T @ np.linalg.solve(A, B)
                out.setdefault(var, []).append(np.linalg.eigvalsh(Pf.T @ Gr @ Pf).min() if evA > 0 else -1)
        print(u, {k: round(min(v), 5) for k, v in out.items()}, flush=True)
