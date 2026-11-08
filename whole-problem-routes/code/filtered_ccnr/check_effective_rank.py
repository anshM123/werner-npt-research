"""q_2(1/2, C) < 0 for an operator with ||C||_1^2 = 2 ||C||_2^2 ("effective rank two"): exact rank two is essential."""
import itertools, os, numpy as np
d, n = 3, 2
C = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'effective_rank_counterexample_n2.npy')).reshape(d**n, d**n)
def ptr(X, J):
    T = X.reshape([d] * (2 * n)); m = n
    for i in sorted([i for i in range(n) if J[i]], reverse=True):
        T = np.trace(T, axis1=i, axis2=i + m); m -= 1
    k = n - sum(J); return T.reshape(d**k, d**k)
q = sum((-0.5) ** sum(J) * np.linalg.norm(ptr(C, J)) ** 2 for J in itertools.product([0, 1], repeat=n))
s = np.linalg.svd(C, compute_uv=False)
print('q_2 / ||C||^2 =', q / np.linalg.norm(C) ** 2)
print('||C||_1^2 / ||C||_2^2 =', s.sum() ** 2 / (s ** 2).sum())
print('singular values / largest:', np.round(s[:4] / s[0], 4))
