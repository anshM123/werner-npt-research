"""Conjecture QC (compression-notes) on the star near-zeros of Section 7.
C = Z^d|g><g| + e^{i phi} Z^{d-v}|g><g|Z^v (star K_{1,n}, d = 0 on centre and k leaves, 1 on m leaves, v = e_c),
phi minimising q_N(1/2, C).  QC asks for a copy s and unit w with q_N(C) >= q_N(P C P), P = I - |w><w| at copy s.
We maximise Delta_s(w) = q_N(C) - q_N(C_{w-perp}) over w (random starts + Nelder-Mead) at every copy.
usage: python qc_star.py k m restarts seed   (dense; N = k + m + 1 <= 6 is quick)"""
import sys, itertools
import numpy as np
from scipy.optimize import minimize
k, m, R, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(seed)
n = k + m; N = n + 1; t = 0.5; D = 3 ** N
om = np.exp(2j * np.pi / 3)
X = np.array(list(itertools.product(range(3), repeat=N)), dtype=np.int64)
G = np.zeros((N, N), dtype=np.int64); G[0, 1:] = 1; G[1:, 0] = 1
ph = np.zeros(D, dtype=np.int64)
for i in range(N):
    for j in range(i + 1, N):
        ph += G[i, j] * X[:, i] * X[:, j]
g = om ** (ph % 3) / np.sqrt(D)
d = np.array([0] + [0] * k + [1] * m); v = np.zeros(N, dtype=np.int64); v[0] = 1
def Zvec(u):
    return om ** ((X @ u) % 3)
C1 = np.outer(Zvec(d) * g, g.conj())
C2 = np.outer(Zvec((d - v) % 3) * g, (Zvec(-v % 3) * g).conj())     # Z^{d-v} g g^* Z^{v}
def q(C):
    T = C.reshape((3,) * (2 * N)); tot = 0.0
    for K in itertools.product([0, 1], repeat=N):
        A = T
        # trace over copies in K: contract axis i with axis N+i, from the highest index down
        idx = [i for i in range(N) if K[i]]
        cur_n = N
        for i in sorted(idx, reverse=True):
            A = np.trace(A, axis1=i, axis2=i + cur_n)
            cur_n -= 1
        tot += (-t) ** sum(K) * np.vdot(A, A).real
    return tot
q11 = q(C1); q22 = q(C2); q12 = (q(C1 + C2) - q11 - q22) / 2 - 1j * (q(C1 + 1j * C2) - q11 - q22) / 2
phi = np.pi - np.angle(q12)
C = C1 + np.exp(1j * phi) * C2
qC = q(C)
print(f"star N={N} (k={k}, m={m}): q11 = {q11:.6e}, |q12| = {abs(q12):.6e}, q_N(C)/||C||^2 = {qC/2:.6e}")
def compress(C, s, w):
    w = w / np.linalg.norm(w)
    Ps = np.eye(3) - np.outer(w, w.conj())
    T = C.reshape((3,) * (2 * N))
    T = np.tensordot(Ps, T, axes=([1], [s]))          # rows of copy s
    T = np.moveaxis(T, 0, s)
    T = np.tensordot(T, Ps.T, axes=([N + s], [0]))    # columns of copy s:  C P  (P Hermitian)
    T = np.moveaxis(T, -1, N + s)
    return T.reshape(D, D)
best_all = -np.inf
for s in range(N):
    best = -np.inf
    for r in range(R):
        x0 = rng.normal(size=6)
        f = lambda x: -(qC - q(compress(C, s, x[:3] + 1j * x[3:])))
        res = minimize(f, x0, method='Nelder-Mead', options={'maxiter': 400, 'xatol': 1e-7, 'fatol': 1e-12})
        best = max(best, -res.fun)
    best_all = max(best_all, best)
    print(f"  copy {s}: max_w Delta = {best:+.6e}")
print(f"  max over copies and w of Delta = {best_all:+.6e}  (QC needs >= 0); q_N(C) = {qC:.6e}")
