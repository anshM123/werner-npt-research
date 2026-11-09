"""(MP) at n = 7 (whole-problem-routes/README.md, remark after Proposition M): alternating minimisation of
<phi| T_S F^(x)7 (psi psi^*) |phi> with Lanczos least eigenpairs, from GHZ, plane-GHZ, singlet-pair and random starts.
usage: python completion_five_sixths_n7.py |S| seed random_starts   (the log used |S| = 1, 2, 3, seeds 71, 72, 73, 4 starts)"""
import numpy as np, sys, itertools, time
from scipy.sparse.linalg import eigsh
rng = np.random.default_rng(int(sys.argv[2]))
I3 = np.eye(3)
def lin(a, b, c):
    def f(Tp):
        tr = np.trace(Tp, axis1=-2, axis2=-1)
        return b*Tp + c*np.swapaxes(Tp, -1, -2) + a*tr[..., None, None]*I3
    return f
F = lin(5/12, -1/12, -2/12); TF = lin(5/12, -2/12, -1/12)
n = 7; k = int(sys.argv[1]); maps = [TF]*k + [F]*(n-k); D = 3**n
def apply(X):
    T = X.reshape([3]*(2*n))
    for s in range(n):
        Tp = np.moveaxis(T, [s, n+s], [2*n-2, 2*n-1]); T = np.moveaxis(maps[s](Tp), [2*n-2, 2*n-1], [s, n+s])
    return T.reshape(D, D)
def least(Y, v0):
    Y = (Y + Y.conj().T)/2
    w, V = eigsh(Y, k=1, which='SA', v0=v0, tol=1e-12, maxiter=5000)
    return w[0], V[:, 0]
seeds = []
g = np.zeros(D); 
for i in range(3): g[sum(i*3**s for s in range(n))] = 1
seeds.append(g)
p = np.array([1, 1j, 0])/np.sqrt(2); a_, b_ = np.ones(1, complex), np.ones(1, complex)
for _ in range(n): a_, b_ = np.kron(a_, p), np.kron(b_, p.conj())
seeds.append(a_ + b_)
sing = np.kron(I3[0], I3[1]) - np.kron(I3[1], I3[0]); v = np.ones(1)
for _ in range(3): v = np.kron(v, sing)
seeds.append(np.kron(v, I3[2]))
for _ in range(int(sys.argv[3])): seeds.append(rng.normal(size=D) + 1j*rng.normal(size=D))
best = np.inf; t0 = time.time()
for psi in seeds:
    psi = np.asarray(psi, complex) + 1e-3*(rng.normal(size=D)+1j*rng.normal(size=D)); psi /= np.linalg.norm(psi)
    old = np.inf; phi = None
    for it in range(150):
        e1, phi = least(apply(np.outer(psi, psi.conj())), phi)
        e2, psi = least(apply(np.outer(phi, phi.conj())), psi)
        if abs(e2 - old) < 1e-16: break
        old = e2
    best = min(best, e2)
    print(f"  start done: {e2:.6e} (x 6^7 = {e2*6**7:.6f}), elapsed {time.time()-t0:.0f}s", flush=True)
print(f"n = 7, |S| = {k}: least value found {best:.6e}; 6^-7 = {6.0**-7:.6e}; ratio {best*6**7:.6f}")
