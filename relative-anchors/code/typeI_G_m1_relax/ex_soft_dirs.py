import numpy as np
from numlib import block
np.set_printoptions(linewidth=200, precision=5, suppress=True)
def Gunscaled(st, sU):
    P, lab, ne = block(st, sU, scaled=False)
    C = P[:ne,:ne]; B = P[:ne,ne:]; A = P[ne:,ne:]
    return A - B.T @ np.linalg.solve(C, B), C, B, A
rng = np.random.default_rng(5)
u = 4
base = np.exp(rng.uniform(-1,1,u)); base /= base.sum()
for lam in [1e-4,1e-3,1e-2,1e-1,1,10,100,1e3,1e4]:
    st = 1.0; sU = base*lam   # S = lam, sigma = 1
    G, C, B, A = Gunscaled(st, sU)
    S = sU.sum(); wU = np.sqrt(sU); wh = wU/np.linalg.norm(wU)
    c = np.sqrt(1 + S/st)
    Pt = np.zeros((u+1,u+1)); Pt[0,0]=1; Pt[1:,1]=c*wh
    Bm = np.linalg.svd(wh[None,:])[2][1:].T; Pt[1:,2:]=Bm
    Gh = Pt.T@G@Pt
    ev, V = np.linalg.eigh(Gh)
    evC = np.linalg.eigvalsh(C)
    print(f'S/sigma={lam:8.0e}  eig(Ghat)={ev}  min eig C={evC.min():.3e}  vec_min={V[:,0]}')
