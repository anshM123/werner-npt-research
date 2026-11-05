import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI
np.set_printoptions(linewidth=200, precision=4, suppress=True)

def analyse(w, T, nshow=4):
    ti = TypeI(w, T, F=float)
    M = np.array(ti.build())
    k = len(w); U = ti.U; lab = ti.labels
    w = np.asarray(w, float)
    # kernel vectors p (x) f_h
    Kv = []
    for h in U:
        v = np.zeros(ti.n)
        for r in range(k):
            v[ti.idx[('P', r, h)]] = w[r]
        Kv.append(v)
    Kv = np.array(Kv).T
    res = np.abs(M @ Kv).max()
    # restrict to orthogonal complement of kernel
    Qk, _ = np.linalg.qr(Kv)
    Pc = np.eye(ti.n) - Qk @ Qk.T
    ev, evec = np.linalg.eigh(Pc @ M @ Pc + 10 * Qk @ Qk.T)
    print(f"w={w}, T={T}, s_T={np.sum(w[T]**2):.3g}, s_U={np.sum(w[U]**2):.3g}; |M kernel|={res:.1e}; lowest eigs {ev[:nshow]}")
    for j in range(min(nshow, 3)):
        v = evec[:, j]
        # describe: X_U, X_T, rho, zeta, q
        XU = np.array([[v[ti.idx[('P', r, h)]] for h in U] for r in U])
        XT = np.array([[v[ti.idx[('P', t, h)]] for h in U] for t in T])
        rho = np.array([v[ti.idx[('R', t)]] for t in T])
        zeta = np.array([[v[ti.idx[('Z', t, h)]] for h in U] for t in T])
        q = v[ti.idx[('Q',)]]
        print(f"  eig {ev[j]:.3e}: |XU|={np.linalg.norm(XU):.3f} |XT|={np.linalg.norm(XT):.3f} |rho|={np.linalg.norm(rho):.3f} |zeta|={np.linalg.norm(zeta):.3f} q={q:.3f}")
        print("   XU=\n", XU)
        print("   XT=", XT, " rho=", rho, "\n   zeta=", zeta)
    return ev

rng = np.random.default_rng(3)
for eps in (1e-1, 1e-2, 1e-3):
    w = np.array([1.0, 0.7, 1.3, 0.9])
    T = [0, 1]
    w2 = w.copy(); w2[T] *= np.sqrt(eps)
    analyse(w2, T)
