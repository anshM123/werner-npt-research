import numpy as np, time
from numlib import block
def Gun(st, sU):
    P, lab, ne = block(st, sU, scaled=False)
    C = P[:ne,:ne]; B = P[:ne,ne:]; A = P[ne:,ne:]
    G = A - B.T @ np.linalg.solve(C, B); return (G+G.T)/2
def frame(st, sU):
    u = len(sU); S = sU.sum(); wU = np.sqrt(sU); wh = wU/np.linalg.norm(wU)
    Pt = np.zeros((u+1,u+1)); Pt[0,0] = np.sqrt(1 + st/S); Pt[1:,1] = np.sqrt(1 + S/st)*wh
    if u > 1: Pt[1:,2:] = np.linalg.svd(wh[None,:])[2][1:].T
    return Pt
rng = np.random.default_rng(9)
t0 = time.time()
for u in (3, 4, 6, 8, 12):
    worst = (1e9, None)
    for it in range(60 if u < 12 else 20):
        fam = it % 4
        st = float(np.exp(rng.uniform(-8, 8)))
        if fam == 0: sU = np.exp(rng.uniform(-3, 3, u))
        elif fam == 1: sU = np.exp(rng.uniform(-8, 8, u))
        elif fam == 2: sU = np.ones(u); sU[0] = np.exp(rng.uniform(-10, 10))
        else: sU = np.exp(np.linspace(0, rng.uniform(-10,10), u))
        G = Gun(st, sU); P = frame(st, sU)
        ev = np.linalg.eigvalsh(P.T @ G @ P).min()
        if ev < worst[0]: worst = (ev, (st, sU.sum(), np.round(np.log(sU),1)))
    print(u, 'min eig in frame (rho*sqrt(1+sig/S), merged*sqrt(1+S/sig), perp):', worst[0], worst[1], f'{time.time()-t0:.1f}s', flush=True)
