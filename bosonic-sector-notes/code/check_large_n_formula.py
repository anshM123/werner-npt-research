"""Independent check of the bosonic q_n formula used by check_large_n.py at n = 5: embed C in (C^3)^(x)5 via the symmetric isometry and
compute q_n = sum_J (-1/2)^|J| ||Tr_J C||^2 directly."""
import itertools, numpy as np, sys
sys.argv = ['x', '5', '1', '0']
exec(open('check_large_n.py').read().split("B = basis(n)")[0])
from math import factorial
B = basis(n); D = len(B)
# isometry Sym^n -> (C^3)^n : |a,b,c> -> sqrt(a! b! c!/n!) sum over words with those counts
V = np.zeros((3 ** n, D))
for k, t in enumerate(B):
    for w in itertools.product(range(3), repeat=n):
        if (w.count(0), w.count(1), w.count(2)) == t:
            V[np.ravel_multi_index(w, [3] * n), k] = 1.0
    V[:, k] /= np.linalg.norm(V[:, k])
rng = np.random.default_rng(3)
for trial in range(3):
    Y = rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2)); Z = rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2))
    C = Y @ Z.conj().T; Cf = V @ C @ V.T
    T = Cf.reshape([3] * (2 * n)); qd = 0
    for r in range(n + 1):
        for J in itertools.combinations(range(n), r):
            X = T
            for j in sorted(J, reverse=True):
                m = X.ndim // 2; X = np.trace(X, axis1=j, axis2=j + m)
            qd += (-0.5) ** r * np.sum(np.abs(X) ** 2)
    qf = q_and_grad(C)[0]
    print(f"trial {trial}: direct {qd:.12e}  formula {qf:.12e}  rel diff {abs(qd - qf) / abs(qd):.1e}")
