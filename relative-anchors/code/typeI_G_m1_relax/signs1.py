"""Exploration: sign pattern of the relaxed form Psi (variant V2, unconstrained a) in scaled coordinates
a_r = w_r ah_r, c_h = w_h ch_h, zeta_h = w_t w_h z_h, rho.  Is it a Z-matrix after a sign flip per variable class?"""
import numpy as np
from forms import Forms

def Psi_full(st, sU, lamkind='H'):
    Fm = Forms(st, sU); u = Fm.u; w = Fm.w; S = Fm.S
    Pi = np.eye(u) - np.outer(w, w) / S
    rho_ = Fm.sU / Fm.d; g = rho_.sum()
    lam = np.full((u,u), 5/12); np.fill_diagonal(lam, 1.5 - rho_)
    n = 3*u + 1
    def Q(v):
        a = w * v[:u]; c = w * v[u:2*u]; rho = v[2*u]; zeta = Fm.wt * w * v[2*u+1:]
        Om = np.outer(a, w) + np.outer(w, c)
        Hf = Fm.Hfull(rho, zeta); Gc = Fm.Gcal(a, c)
        R = st * (Pi @ Gc @ Pi) + Pi @ Hf @ Pi
        return st * Fm.f(a, c) + np.sum(Hf * Om) + Fm.Acal(rho, zeta) - 0.25 * np.sum(R**2 / lam)
    M = np.zeros((n, n)); E = np.eye(n)
    dg = [Q(E[i]) for i in range(n)]
    for i in range(n):
        M[i,i] = dg[i]
        for j in range(i+1, n):
            M[i,j] = M[j,i] = (Q(E[i]+E[j]) - dg[i] - dg[j]) / 2
    return M

if __name__ == '__main__':
    rng = np.random.default_rng(8)
    u = 4
    SG = []
    mins = []
    for it in range(200):
        st = float(np.exp(rng.uniform(-5, 5))); sU = np.exp(rng.uniform(-4, 4, u))
        M = Psi_full(st, sU)
        mins.append(np.linalg.eigvalsh(M).min() / np.abs(np.diag(M)).max())
        SG.append(np.sign(np.where(np.abs(M) < 1e-13 * np.abs(M).max(), 0, M)))
    SG = np.array(SG)
    mixed = (SG.max(0) > 0) & (SG.min(0) < 0)
    pat = np.where(mixed, 9, SG[0]).astype(int)
    names = [f'a{i}' for i in range(u)] + [f'c{i}' for i in range(u)] + ['rho'] + [f'z{i}' for i in range(u)]
    print('min eig (unconstrained a), normalised:', min(mins))
    print('      ' + ' '.join(f'{n:>4}' for n in names))
    for i in range(3*u+1):
        print(f'{names[i]:>5} ' + ' '.join(f'{x:>4}' for x in pat[i]))
