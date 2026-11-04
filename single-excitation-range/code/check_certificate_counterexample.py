"""The certificates of Proposition 6 do not exist for every POVM: counterexample at n = 7, c = 5 (and n = 8, c = 6).

Family: one heavy copy A_0 = diag(1 - d, 1) and N = n - 1 identical light copies U_t = sqrt(d/N) |1><e_1|, sector
c = n - 2.  The certificate program is convex and the family is invariant under permutations of the light copies, so a
certificate exists iff a symmetric one does.  A symmetric certificate has five unknowns blocks:
  theta_L (Q = all light copies), theta_H (Q containing the heavy copy),
  m_A0, m_AL (patterns containing the heavy copy; block at the heavy / a light copy), m_BL (patterns without it),
and it suffices to impose one pattern inequality and one string inequality per orbit.
(a) for n = 5, 6 the symmetric program reproduces the full program of combined_certificate.py (same maximal slack);
(b) the maximal slack is negative at n = 7, d = 0.2 and at n = 8, d = 0.05: no certificate exists there;
(c) the sector form itself is positive at the n = 7 counterexample (lambda_min of the explicit sector matrix).
"""
import numpy as np, cvxpy as cp, warnings
from combined_certificate import perm_to_front, certificate, sector_matrix
warnings.filterwarnings('ignore')

def heavy(N, d):
    return [np.array([[np.sqrt(1 - d), 0], [0, 1.0]], complex)] + [np.array([[np.sqrt(d / N), 0], [0, 0]], complex) for _ in range(N)]

def symmetric_slack(N, d):
    n = N + 1; c = n - 2; U = heavy(N, d); A = [u.conj().T @ u for u in U]
    ys = [np.conj(np.array([U[t][z, r] for z in range(2) for r in range(2)])) for t in range(n)]
    thL, thH = cp.Variable(nonneg=True), cp.Variable(nonneg=True)
    mA0, mAL, mBL = [cp.Variable((4, 4), hermitian=True) for _ in range(3)]; lam = cp.Variable()
    theta = lambda Q: thL if 0 not in Q else thH
    mvar = lambda S, t: (mA0 if t == 0 else mAL) if 0 in S else mBL
    def pattern(S):
        H = [x for x in range(n) if x not in S]
        budget = np.kron(np.eye(2 ** c), 2 * np.eye(2) + 2 * sum((A[x] for x in H), np.zeros((2, 2)))); tot = 0
        for p, t in enumerate(S):
            Pm = perm_to_front(c, p); tot = tot + Pm.T @ cp.kron(mvar(S, t), np.eye(2 ** (c - 1))) @ Pm
        for u in H:
            Q = tuple(sorted(S + (u,))); tot = tot + 2 * theta(Q) * np.kron(np.eye(2 ** c), sum(A[x] for x in Q))
        return budget - tot - lam * np.eye(2 ** (c + 1)) >> 0
    def string(R):
        free = [t for t in range(n) if t not in R]; k = len(free); blocks = [[None] * k for _ in range(k)]
        for i in range(k):
            for j in range(k):
                if i == j: Nij = 1.0 - sum(theta(tuple(sorted(R + (free[i], free[l])))) for l in range(k) if l != i)
                else: Nij = -1.0 + theta(tuple(sorted(R + (free[i], free[j]))))
                blocks[i][j] = (mvar(tuple(sorted(R + (free[i],))), free[i]) if i == j else 0) - Nij * np.outer(ys[free[i]], ys[free[j]].conj())
        return cp.bmat(blocks) - lam * np.eye(4 * k) >> 0
    L = list(range(1, n))
    cons = [pattern(tuple([0] + L[:c - 1])), pattern(tuple(L[:c])), string(tuple([0] + L[:c - 2])), string(tuple(L[:c - 1]))]
    cp.Problem(cp.Maximize(lam), cons).solve(solver='CLARABEL')
    return lam.value

ok = True
for N, d in [(4, 0.05), (4, 0.2), (5, 0.05), (5, 0.2)]:
    s, (f, _, _) = symmetric_slack(N, d), certificate(heavy(N, d), N + 1, N - 1)
    print(f"(a) n={N + 1} c={N - 1} d={d}: symmetric slack {s:+.6e}, full program {f:+.6e}", flush=True)
    ok &= abs(s - f) < 1e-6 and s > 0
for N, d in [(6, 0.05), (6, 0.2), (7, 0.05)]:
    s = symmetric_slack(N, d); print(f"(b) n={N + 1} c={N - 1} d={d}: maximal slack {s:+.4e}", flush=True)
    if (N, d) in [(6, 0.2), (7, 0.05)]: ok &= s < -1e-5
mu = np.linalg.eigvalsh(sector_matrix(heavy(6, 0.2), 7, 5))[0]
print(f"(c) n=7 c=5 d=0.2: lambda_min of the sector form {mu:+.4e} (the sector itself is positive)")
ok &= mu > 0.1
assert ok
print("PASS")
