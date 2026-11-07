"""Large-n test on translation-invariant inputs (numerics for Observation 7.4).

For periodic real matrix product states x (tensors A^i, D x D) and y (B^k) on n copies,
    V_n = <y| P_a^n (x x^T) |y> = Tr(E^n),
    E = a (sum_i A^i (x) A^i) (x) (sum_k B^k (x) B^k) - 1/4 sum_ij [A^i (x) A^j (x) B^i (x) B^j + A^i (x) A^j (x) B^j (x) B^i].
If V_n >= 0 for every n, the series sum_n V_n z^n has nonnegative coefficients, so by Pringsheim's theorem the spectral
radius of E is an eigenvalue.  The script maximises
    ratio = (largest |lambda| over eigenvalues of E that are not positive reals) / (largest positive eigenvalue)
over A, B (random starts, starts near the plane witness of Theorem 7.1, starts with y close to x; Nelder-Mead then
Powell).  ratio > 1 would mean V_n < 0 for infinitely many n.  ratio = 1 is reached trivially, e.g. by states of period
two, whose transfer operators have spectrum symmetric under lambda -> -lambda.

usage: python ti_transfer_search.py a D starts seed
"""
import numpy as np, sys
from scipy.optimize import minimize
a = float(sys.argv[1]); D = int(sys.argv[2]); starts = int(sys.argv[3]); seed = int(sys.argv[4])
rng = np.random.default_rng(seed)
def Emat(A, B):
    EA = sum(np.kron(A[i], A[i]) for i in range(3)); EB = sum(np.kron(B[k], B[k]) for k in range(3))
    E = a*np.kron(EA, EB)
    for i in range(3):
        for j in range(3):
            AA = np.kron(A[i], A[j])
            E = E - 0.25*(np.kron(AA, np.kron(B[i], B[j])) + np.kron(AA, np.kron(B[j], B[i])))
    return E
def ratio(v):
    A = v[:3*D*D].reshape(3, D, D); B = v[3*D*D:].reshape(3, D, D)
    w = np.linalg.eigvals(Emat(A, B))
    pos = w[(np.abs(w.imag) < 1e-9*max(1, np.abs(w).max())) & (w.real > 0)].real
    bad = w[~((np.abs(w.imag) < 1e-9*max(1, np.abs(w).max())) & (w.real > 0))]
    if len(pos) == 0: return 10.0
    rb = np.abs(bad).max() if len(bad) else 0.0
    return rb/pos.max()
def plane_ghz(D):
    A = np.zeros((3, D, D)); A[0, 0, 0] = A[0, 1, 1] = 1/np.sqrt(2); A[1, 0, 1] = -1/np.sqrt(2); A[1, 1, 0] = 1/np.sqrt(2)
    return A
best = 0; bestv = None
g = plane_ghz(max(D, 2))[:, :D, :D] if D >= 2 else None
if D >= 2:
    v0 = np.concatenate([g.reshape(-1), g.reshape(-1)])
    print(f'plane GHZ (D={D}) at a={a}: ratio {ratio(v0):.6f}', flush=True)
for st in range(starts):
    v0 = rng.standard_normal(2*3*D*D)
    if st % 3 == 1 and D >= 2: v0 = np.concatenate([g.reshape(-1), g.reshape(-1)]) + 0.3*rng.standard_normal(2*3*D*D)
    if st % 3 == 2: v0[3*D*D:] = v0[:3*D*D] + 0.2*rng.standard_normal(3*D*D)
    r = minimize(lambda v: -ratio(v), v0, method='Nelder-Mead', options={'maxiter': 4000, 'maxfev': 4000})
    r = minimize(lambda v: -ratio(v), r.x, method='Powell', options={'maxiter': 4000})
    val = -r.fun
    if val > best: best, bestv = val, r.x
    print(f'start {st}: ratio {val:.6f}  (best {best:.6f})', flush=True)
print(f'a = {a}, D = {D}: largest ratio found {best:.6f}')
