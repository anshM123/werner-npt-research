"""Section 9, purified form: alpha = sum_i |i> (x) a_i, beta = sum_i |i> (x) b_i, C = sum_i a_i b_i^*.
Checks ||Tr_K C||^2 = Tr(rho_{RK} sigma_{RK}), Tr C = <beta, alpha>, ||C||_1 = F(rho_R, sigma_R)."""
import itertools, numpy as np
from scipy.linalg import sqrtm
rng = np.random.default_rng(11)
n = 3; d = 3; D = d ** n
a = rng.standard_normal((2, D)) + 1j * rng.standard_normal((2, D)); b = rng.standard_normal((2, D)) + 1j * rng.standard_normal((2, D))
C = np.outer(a[0], b[0].conj()) + np.outer(a[1], b[1].conj())
alpha = a.reshape([2] + [d] * n); beta = b.reshape([2] + [d] * n)
def marg(psi, keep):          # reduced state of |psi> on axes 'keep' (axis 0 = reference)
    axes = list(range(psi.ndim)); tr = [x for x in axes if x not in keep]
    X = np.tensordot(psi, psi.conj(), axes=(tr, tr)); k = len(keep)
    dim = int(np.prod([psi.shape[x] for x in keep])); return X.reshape(dim, dim)
mx = 0
for r in range(n + 1):
    for K in itertools.combinations(range(n), r):
        T = C.reshape([d] * 2 * n); rem = list(range(n)); X = T
        for s in sorted(K, reverse=True):
            i = rem.index(s); X = np.trace(X, axis1=i, axis2=len(rem) + i); rem.pop(i)
        lhs = np.sum(abs(X) ** 2)
        keep = [0] + [s + 1 for s in K]
        rhs = np.trace(marg(alpha, keep) @ marg(beta, keep)).real
        mx = max(mx, abs(lhs - rhs) / lhs)
print(f"max relative residual of ||Tr_K C||^2 = Tr(rho_RK sigma_RK) over all K: {mx:.1e}")
print(f"Tr C = {np.trace(C):.10f}   <beta,alpha> = {np.vdot(b.reshape(-1), a.reshape(-1)):.10f}")
rR, sR = marg(alpha, [0]), marg(beta, [0])
sq = sqrtm(rR); F = np.trace(sqrtm(sq @ sR @ sq)).real
print(f"||C||_1 = {np.linalg.svd(C, compute_uv=False).sum():.10f}   F(rho_R, sigma_R) = {F:.10f}")
