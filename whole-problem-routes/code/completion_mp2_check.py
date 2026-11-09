"""(MP2) of the remark after Proposition M: is T_S F^{(x)m} two-positive?  Least value of
<phi|(id_r (x) T_S F^m)(psi psi^*)|phi> over psi in C^r (x) (C^3)^m and phi (alternating minimisation).
usage: python completion_mp2_check.py m |S| seed r starts   (logs/completion_five_sixths_mp2_m5.log: m = 5, |S| = 5, r = 2)"""
import numpy as np, sys, itertools
rng = np.random.default_rng(int(sys.argv[3]))
I3 = np.eye(3)
def lin(a, b, c):
    def f(Tp):
        tr = np.trace(Tp, axis1=-2, axis2=-1)
        return b*Tp + c*np.swapaxes(Tp, -1, -2) + a*tr[..., None, None]*I3
    return f
F = lin(5/12, -1/12, -2/12); TF = lin(5/12, -2/12, -1/12)
def apply_ref(X, m, maps, r):
    """X on C^r (x) (C^3)^m; apply maps on the m qutrit sites."""
    n = m
    T = X.reshape([r] + [3]*m + [r] + [3]*m)
    for s in range(m):
        a, b = 1+s, 1+m+1+s
        Tp = np.moveaxis(T, [a, b], [-2, -1]); T = np.moveaxis(maps[s](Tp), [-2, -1], [a, b])
    return T.reshape(r*3**m, r*3**m)
m = int(sys.argv[1]); k = int(sys.argv[2]); r = int(sys.argv[4]) if len(sys.argv) > 4 else 2
maps = [TF]*k + [F]*(m-k); D = r*3**m; best = np.inf
for start in range(int(sys.argv[5]) if len(sys.argv) > 5 else 30):
    psi = rng.normal(size=D)+1j*rng.normal(size=D); psi /= np.linalg.norm(psi); old = np.inf
    for it in range(500):
        Y = apply_ref(np.outer(psi, psi.conj()), m, maps, r); e1, U1 = np.linalg.eigh((Y+Y.conj().T)/2); phi = U1[:, 0]
        Z = apply_ref(np.outer(phi, phi.conj()), m, maps, r); e2, U2 = np.linalg.eigh((Z+Z.conj().T)/2); psi = U2[:, 0]
        if abs(e2[0]-old) < 1e-15: break
        old = e2[0]
    best = min(best, e2[0])
print(f"m = {m}, |S| = {k}, reference dim {r}: least value {best:+.6e}  (6^-m = {6.0**-m:.6e})")
