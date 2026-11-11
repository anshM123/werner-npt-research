"""Independent check of data/compression_counterexample_n3.npy (numpy and scipy only, no torch).
q_n is computed from partial traces, q_n(X) = sum_K (-1/2)^{|K|} ||Tr_K X||^2, and the leak vector at each copy is
parametrised by four angles, l = (cos a, sin a cos b e^{ic}, sin a sin b e^{ie}).  The minimum of R(P C P) over product
qubits P = (x)_s (I - l_s l_s^*) is searched by BFGS from random angles.
usage: python qubit_compression_check.py FILE starts"""
import sys
import numpy as np
from scipy.optimize import minimize

C = np.load(sys.argv[1]); NS = int(sys.argv[2])
n = round(np.log(C.shape[0]) / np.log(3)); C = C / np.linalg.norm(C)


def q(X):
    T = X.reshape((3,) * (2 * n)); tot = 0.0
    for mask in range(2 ** n):
        K = [s for s in range(n) if mask >> s & 1]; Y = T
        for s in sorted(K, reverse=True):
            Y = np.trace(Y, axis1=s, axis2=Y.ndim // 2 + s)
        tot += (-0.5) ** len(K) * np.sum(np.abs(Y) ** 2)
    return tot


def R(X):
    return q(X) / np.sum(np.abs(X) ** 2)


def f(x):
    P = np.ones((1, 1))
    for s in range(n):
        a, b, c, e = x[4 * s:4 * s + 4]
        l = np.array([np.cos(a), np.sin(a) * np.cos(b) * np.exp(1j * c), np.sin(a) * np.sin(b) * np.exp(1j * e)])
        P = np.kron(P, np.eye(3) - np.outer(l, l.conj()))
    X = P @ C @ P
    return R(X) if np.sum(np.abs(X) ** 2) > 1e-8 else 10.0


rng = np.random.default_rng(2024)
vals = np.sort([minimize(f, rng.uniform(0, 2 * np.pi, 4 * n), method='BFGS', options={'gtol': 1e-10, 'maxiter': 2000}).fun
                for _ in range(NS)])
print(f"{sys.argv[1].replace(chr(92), '/').split('/')[-1]}: singular values {np.round(np.linalg.svd(C, compute_uv=False)[:3], 9).tolist()}")
print(f"R(C) from partial traces = {R(C):.10f}")
print(f"min over {NS} BFGS runs of R(P C P) = {vals[0]:.10f};  gap = {vals[0] - R(C):+.10f};  five smallest {np.round(vals[:5], 9).tolist()}")
