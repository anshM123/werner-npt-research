"""Check the four-slot encoding and the single-sandwich condition at n = 1.

Single sandwich (tensor-stable): if for one slot set S both (A+B)^{G_S} >= 0 and (A-B)^{G_S} >= 0,
then (A^{(x)n} - B^{(x)n})^{G_S} >= 0 for all n, hence crossed CS (= EP_n) for all n.
"""
import numpy as np
from scipy.optimize import minimize
from common import AB, ptrans, subsets, rand_vec

rng = np.random.default_rng(1)

for t in [0.5, 0.55]:
    A, B = AB(t)
    # identity check on random vectors
    W = np.eye(9) - t * np.eye(9)[:, [3 * (k % 3) + k // 3 for k in range(9)]]  # I - tF on C3xC3
    a, b, c, dd = (rand_vec(3, rng) for _ in range(4))
    Phi = np.kron(np.kron(np.kron(a, dd), c), b)
    lhsA = np.vdot(Phi, A @ Phi).real
    rhsA = np.vdot(np.kron(a, b), W @ np.kron(a, b)).real * np.vdot(np.kron(c, dd), W @ np.kron(c, dd)).real
    lhsB = np.vdot(Phi, B @ Phi).real
    rhsB = abs(np.vdot(np.kron(a, dd), W @ np.kron(c, b))) ** 2
    print(f"t={t}: encoding check A {lhsA:.6f} vs {rhsA:.6f};  B {lhsB:.6f} vs {rhsB:.6f}")

    # min of crossed-CS form over product vectors (n = 1), local optimisation
    def f(x):
        v = x[:12].reshape(4, 3) + 1j * x[12:].reshape(4, 3)
        v = v / np.linalg.norm(v, axis=1, keepdims=True)
        P = np.kron(np.kron(np.kron(v[0], v[1]), v[2]), v[3])
        return np.vdot(P, (A - B) @ P).real
    best = min(minimize(f, rng.normal(size=24), method='BFGS').fun for _ in range(30))
    print(f"   min over product vectors of <Phi|A-B|Phi> (n=1): {best:.3e}")

    # single sandwich, all 16 slot sets
    res = []
    for S in subsets():
        lp = np.linalg.eigvalsh(ptrans(A + B, S)).min()
        lm = np.linalg.eigvalsh(ptrans(A - B, S)).min()
        res.append((S, lp, lm))
    for S, lp, lm in res:
        print(f"   S={str(S):14s} lam_min[(A+B)^G]={lp:+.4f}  lam_min[(A-B)^G]={lm:+.4f}")
