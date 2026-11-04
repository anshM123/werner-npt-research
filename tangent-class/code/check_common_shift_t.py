"""The shared-shift obstruction away from the endpoint.  For the Werner parameter t, q_t(X) = sum_J (-t)^|J| ||Tr_J X||^2
and the single-copy squeeze forms are built with -t in place of -1/2.  The minimiser argument, run over all rank-two
operators at parameter t, would need: no thick C has all its single-copy pencils PSD at one common lambda in [1 - 3t, 0)
(1 - 3t is the smallest eigenvalue of the n-copy map).  Each stored three-copy operator below is thick (range and
co-range) and has all three pencils positive definite at lambda = 1 - 3t, for t = 0.45, 0.40, 0.35 and 0.335."""
import numpy as np, itertools
from scipy.linalg import sqrtm
from common import Bt, _ptrace_keep, Lmap

def q_t(X, n, t):
    return np.vdot(X, Lmap(X, n, t)).real

def q_t_direct(X, n, t):
    T = X.reshape([3] * (2 * n)); tot = 0.0
    for J in itertools.product([0, 1], repeat=n):
        Y = T; cur = list(range(n))
        for s in sorted([s for s in range(n) if J[s]], reverse=True):
            i = cur.index(s); m = len(cur); Y = np.trace(Y, axis1=i, axis2=m + i); cur.pop(i)
        tot += (-t) ** sum(J) * np.sum(abs(Y) ** 2)
    return tot

def squeeze_forms_t(X, n, s, t):
    K = np.zeros((9, 9)); Kn = np.zeros((9, 9)); others = [u for u in range(n) if u != s]
    for Jbits in itertools.product([0, 1], repeat=n - 1):
        J = [others[i] for i in range(n - 1) if Jbits[i]]; c = (-t) ** len(J)
        T, cur = _ptrace_keep(X, n, J); k = len(cur); i = cur.index(s)
        T = np.moveaxis(T, [i, k + i], [0, 1]).reshape(3, 3, -1)
        G = np.einsum('abr,cdr->acbd', T.conj(), T); val = np.einsum('kac,lfb,acbf->kl', Bt, Bt, G).real
        K += c * val
        if not J:
            Kn += val
        Z = np.einsum('kba,abr->kr', Bt, T); K += (-t * c) * np.real(Z.conj() @ Z.T)
    return (K + K.T) / 2, (Kn + Kn.T) / 2

rng = np.random.default_rng(0); n = 3; D = 27; ok = True
for t, name in ((0.45, "0.45"), (0.40, "0.40"), (0.35, "0.35"), (0.335, "0.335")):
    X = np.load(f"../data/common_shift_t{name}_n3.npy"); X /= np.linalg.norm(X); lam = 1 - 3 * t
    margins, ratios, err = [], [], 0
    for s in range(n):
        K, Kn = squeeze_forms_t(X, n, s, t)
        A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); M = A @ A.conj().T
        ops = [np.eye(3)] * n; ops[s] = sqrtm(M); G = np.kron(np.kron(ops[0], ops[1]), ops[2]); Y = G @ X @ G
        c = np.array([np.trace(B.conj().T @ M).real for B in Bt])
        err = max(err, abs(c @ K @ c - q_t_direct(Y, n, t)) / abs(q_t_direct(Y, n, t)))
        margins.append(np.linalg.eigvalsh(K - lam * Kn)[0] / np.linalg.norm(Kn, 2))
        for Yp in (X @ X.conj().T, X.conj().T @ X):
            T = np.moveaxis(Yp.reshape([3] * (2 * n)), [s, n + s], [0, n]).reshape(3, 9, 3, 9)
            ev = np.linalg.eigvalsh(np.einsum('aibi->ab', T)); ratios.append(ev[0] / ev[-1])
    rk = np.linalg.matrix_rank(X, 1e-10)
    print(f"t = {t}: lambda = 1 - 3t = {lam:+.3f}; lambda_min(Q_s - lambda N_s)/||N_s|| at s = 0, 1, 2: "
          f"{', '.join(f'{m_:+.2e}' for m_ in margins)}; rank {rk}; smallest local-support ratio {min(ratios):.3f}; "
          f"R_t = {q_t(X, n, t):.4f}; forms vs direct {err:.1e}")
    ok &= min(margins) > 0 and rk == 2 and min(ratios) > 0.2 and err < 1e-10
assert ok
print("PASS")
