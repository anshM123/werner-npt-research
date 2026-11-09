"""Least output eigenvalue of T_S F_{a,t}^(x)n for the completions F_{a,t} = (P_a + t(id - T)/2)/(3a - 1/2) of P_a
(remark after Proposition M), compared with the product value ((a - 1/2)/(3a - 1/2))^n.
usage: python completion_family_check.py a t n |S| seed starts"""
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[5]))
I3 = np.eye(3)
def lin(a, b, c):
    def f(Tp):
        tr = np.trace(Tp, axis1=-2, axis2=-1)
        return b*Tp + c*np.swapaxes(Tp, -1, -2) + a*tr[..., None, None]*I3
    return f
a, t = float(sys.argv[1]), float(sys.argv[2]); n, k = int(sys.argv[3]), int(sys.argv[4])
c = 3*a - 0.5
F = lin(a/c, (-0.25 + t/2)/c, (-0.25 - t/2)/c); TF = lin(a/c, (-0.25 - t/2)/c, (-0.25 + t/2)/c)
def apply(X, n, maps):
    T = X.reshape([3]*(2*n))
    for s in range(n):
        Tp = np.moveaxis(T, [s, n+s], [2*n-2, 2*n-1]); T = np.moveaxis(maps[s](Tp), [2*n-2, 2*n-1], [s, n+s])
    return T.reshape(3**n, 3**n)
maps = [TF]*k + [F]*(n-k); D = 3**n; best = np.inf
for start in range(int(sys.argv[6])):
    psi = rng.normal(size=D)+1j*rng.normal(size=D)
    if start % 3 == 1: psi = psi.real + 0j
    psi /= np.linalg.norm(psi); old = np.inf
    for it in range(400):
        Y = apply(np.outer(psi, psi.conj()), n, maps); e1, U1 = np.linalg.eigh((Y+Y.conj().T)/2); phi = U1[:, 0]
        Z = apply(np.outer(phi, phi.conj()), n, maps); e2, U2 = np.linalg.eigh((Z+Z.conj().T)/2); psi = U2[:, 0]
        if abs(e2[0]-old) < 1e-15: break
        old = e2[0]
    best = min(best, e2[0])
prod = ((a-0.5)/c)**n
print(f"a = {a}, t = {t}, n = {n}, |S| = {k}: least {best:.6e}, product value {prod:.6e}, ratio {best/prod:.6f}")
