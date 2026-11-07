"""Compressed products Z_e = (Pi x Pi)(P_Om^{e1} (x) P_s0^{e2} (x) P_an^{e3})(Pi x Pi) on Sym^n (x) Sym^n."""
import numpy as np
from bos import mindices, lowering

def one_copy_ops(d=3):
    I = np.eye(d*d)
    F = np.zeros((d*d, d*d))
    for i in range(d):
        for j in range(d):
            F[i*d+j, j*d+i] = 1
    om = np.zeros(d*d)
    for i in range(d):
        om[i*d+i] = 1
    Om = np.outer(om, om)
    POm = Om / d
    Psym = (I + F)/2
    Ps0 = Psym - POm
    Pan = (I - F)/2
    return dict(I=I, F=F, Om=Om, POm=POm, Ps0=Ps0, Pan=Pan)

def transfer(Y, m, X, d=3):
    """T_X(Y) = sum_{ijkl} X[(i,j),(k,l)] (L_i (x) L_j)^T Y (L_k (x) L_l), from level m-1 to level m."""
    Nm1 = len(mindices(m-1, d)); N = len(mindices(m, d))
    low = lowering(m, d)
    r = {}; s = {}
    for i in range(d):
        ti, si = low[i]
        for j in range(d):
            tj, sj = low[j]
            r[i, j] = (ti[:, None]*Nm1 + tj[None, :]).reshape(-1)
            s[i, j] = (si[:, None]*sj[None, :]).reshape(-1)
    out = np.zeros((N*N, N*N))
    for i in range(d):
        for j in range(d):
            for k in range(d):
                for l in range(d):
                    c = X[i*d+j, k*d+l]
                    if c != 0:
                        out += c * s[i, j][:, None] * s[k, l][None, :] * Y[np.ix_(r[i, j], r[k, l])]
    return out

def build_Z(n, d=3):
    ops = one_copy_ops(d)
    P = [ops['POm'], ops['Ps0'], ops['Pan']]
    Z = {(0, 0, 0): np.array([[1.0]])}
    for m in range(1, n+1):
        Znew = {}
        for e1 in range(m+1):
            for e2 in range(m+1-e1):
                e3 = m - e1 - e2
                e = (e1, e2, e3)
                # pick first nonzero type
                for tpos in range(3):
                    if e[tpos] > 0:
                        prev = list(e); prev[tpos] -= 1
                        Znew[e] = transfer(Z[tuple(prev)], m, P[tpos], d)
                        break
        Z = Znew
    return Z
