"""Identity-count distribution of rank-two operators.  For C on (C^3)^{(x)N} let pibar_C(k) = sum_{|K|=k} ||C_K||^2/||C||^2,
C_K = part of C that is the identity exactly on the copies in K.  The endpoint is sum_k pibar(k)(-1/2)^k >= 0.
Support function: s(g) = min over rank-two C of sum_k g(k) pibar_C(k), computed with
  sum_k g(k) pibar_C(k) = sum_K g~(K) ||Tr_K C||^2 / ||C||^2 (Moebius: ||Tr_K C||^2 = 3^{|K|} sum_{L >= K} ||C_L||^2).
usage: python pgf_support.py N directions restarts seed"""
import sys, itertools
import numpy as np
from scipy.optimize import minimize
from math import comb
N, ndir, R, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(seed)
D = 3 ** N
subsets = list(itertools.product([0, 1], repeat=N))
def trace_norms(C):
    """returns array over subsets K of ||Tr_K C||^2"""
    T = C.reshape((3,) * (2 * N)); out = []
    for K in subsets:
        A = T; cur = N
        for i in sorted([i for i in range(N) if K[i]], reverse=True):
            A = np.trace(A, axis1=i, axis2=i + cur); cur -= 1
        out.append(np.vdot(A, A).real)
    return np.array(out)
# ||C_L||^2 by Moebius inversion: mu_K = ||Tr_K C||^2 / 3^{|K|} = sum_{L >= K} e_L  =>  e_L = sum_{K >= L} (-1)^{|K|-|L|} mu_K
sizes = np.array([sum(K) for K in subsets])
sup = np.array([[all(K[i] >= L[i] for i in range(N)) for K in subsets] for L in subsets])   # sup[L][K]: K >= L
mob = np.where(sup, (-1.0) ** (sizes[None, :] - sizes[:, None]), 0.0)
def pibar(C):
    mu = trace_norms(C) / 3.0 ** sizes
    e = mob @ mu
    pb = np.array([e[sizes == k].sum() for k in range(N + 1)])
    return pb / np.vdot(C, C).real
def support(g, restarts):
    best = (np.inf, None)
    for r in range(restarts):
        x0 = rng.normal(size=8 * D)
        def f(x):
            A = (x[:2 * D] + 1j * x[2 * D:4 * D]).reshape(D, 2); B = (x[4 * D:6 * D] + 1j * x[6 * D:]).reshape(D, 2)
            return g @ pibar(A @ B.conj().T)
        res = minimize(f, x0, method='L-BFGS-B', options={'maxiter': 3000})
        if res.fun < best[0]:
            A = (res.x[:2 * D] + 1j * res.x[2 * D:4 * D]).reshape(D, 2); B = (res.x[4 * D:6 * D] + 1j * res.x[6 * D:]).reshape(D, 2)
            best = (res.fun, pibar(A @ B.conj().T))
    return best
# sanity: GHZ-diagonal C = |0^N><0^N| + |1^N><1^N|  -> pibar(k) = 3^-N C(N,k)[2^{N-k} + (-1)^{N-k}]
C = np.zeros((D, D)); C[0, 0] = 1; C[(D - 1) // 2, (D - 1) // 2] = 1     # index of 1...1 in base 3 = (3^N - 1)/2
pb = pibar(C)
pred = np.array([comb(N, k) * (2 ** (N - k) + (-1) ** (N - k)) for k in range(N + 1)]) / 3 ** N
print("GHZ check:", np.abs(pb - pred).max())
print("endpoint functional at GHZ:", pb @ (-0.5) ** np.arange(N + 1))
# directions: the endpoint direction, then random directions; report extreme points
dirs = [(-0.5) ** np.arange(N + 1)] + [rng.normal(size=N + 1) for _ in range(ndir)]
for g in dirs:
    val, p = support(g, R)
    print(f"g = {np.round(g, 3).tolist()}: support {val:+.6f} at pibar = {np.round(p, 5).tolist()}", flush=True)
