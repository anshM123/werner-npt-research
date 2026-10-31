"""Section 11, one-sided separability: is (P_V (x) I) rho^{(x)n} (P_V (x) I) separable (qubit V | Bob's 3^n) for every
plane V of Alice's (C^3)^{(x)n}?  PPT of all these is exactly EP_n.
- Outer test: bosonic K-extension of the qubit side with PPT across (qubit copies | Bob); infeasible => entangled.
- Inner test: eps* = min{eps : sigma + eps I/dim = sum_p |a_p><a_p| (x) B_p, B_p >= 0} over a fixed net of 120 qubit
  states; eps* <= 0 means separable (up to the net).
usage: python check_one_sided.py n K   (n = 2, K = 4: a few minutes; needs cvxpy with Clarabel)
"""
import sys, numpy as np, cvxpy as cp
n = int(sys.argv[1]) if len(sys.argv) > 1 else 2
K = int(sys.argv[2]) if len(sys.argv) > 2 else 4
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2
def power(m):
    R = np.eye(1)
    for _ in range(m): R = np.kron(R, rho)
    perm = [2 * i for i in range(m)] + [2 * i + 1 for i in range(m)]
    return R.reshape([3] * (4 * m)).transpose(perm + [2 * m + q for q in perm]).reshape(9 ** m, 9 ** m)
N = 3 ** n
Rn = power(n)
def sigma_of(V):                      # V: N x 2 with orthonormal columns
    S = np.kron(V.conj().T, np.eye(N)) @ Rn @ np.kron(V, np.eye(N))
    S = (S + S.conj().T) / 2
    return S / np.trace(S).real
def dicke(k):
    from itertools import combinations
    D = np.zeros((2 ** k, k + 1))
    for j in range(k + 1):
        for c in combinations(range(k), j):
            idx = sum(1 << (k - 1 - q) for q in c)
            D[idx, j] = 1
        D[:, j] /= np.linalg.norm(D[:, j])
    return D
def red1(k):
    """R[j][l] = Tr_{2..k} |D_j><D_l| (2 x 2)."""
    D = dicke(k); R = [[None] * (k + 1) for _ in range(k + 1)]
    for j in range(k + 1):
        for l in range(k + 1):
            M = np.outer(D[:, j], D[:, l]).reshape(2, 2 ** (k - 1), 2, 2 ** (k - 1))
            R[j][l] = np.einsum('aibi->ab', M)
    return R
def ptB(X, dA, dB):
    return cp.partial_transpose(X, dims=(dA, dB), axis=1)
def extendible(S, k):
    dB = S.shape[0] // 2; R = red1(k); dA = k + 1
    X = cp.Variable((dA * dB, dA * dB), hermitian=True)
    marg = 0
    for j in range(dA):
        for l in range(dA):
            if np.abs(R[j][l]).max() < 1e-14: continue
            marg = marg + cp.kron(R[j][l], X[j * dB:(j + 1) * dB, l * dB:(l + 1) * dB])
    cons = [X >> 0, ptB(X, dA, dB) >> 0, marg == S]
    prob = cp.Problem(cp.Minimize(0), cons)
    try:
        prob.solve(solver=cp.CLARABEL)
    except Exception as e:
        return f"error {e}"
    return prob.status
def net(P):
    i = np.arange(P) + 0.5; th = np.arccos(1 - 2 * i / P); ph = np.pi * (1 + 5 ** 0.5) * i
    return [np.array([np.cos(t / 2), np.exp(1j * f) * np.sin(t / 2)]) for t, f in zip(th, ph)]
def inner_eps(S, P=120):
    dB = S.shape[0] // 2; A = net(P)
    Bs = [cp.Variable((dB, dB), hermitian=True) for _ in range(P)]
    eps = cp.Variable()
    tot = sum(cp.kron(np.outer(a, a.conj()), B) for a, B in zip(A, Bs))
    cons = [B >> 0 for B in Bs] + [tot == S + eps * np.eye(2 * dB) / (2 * dB)]
    prob = cp.Problem(cp.Minimize(eps), cons)
    prob.solve(solver=cp.CLARABEL)
    return eps.value
def ptmin(S):
    dB = S.shape[0] // 2
    T = S.reshape(2, dB, 2, dB).transpose(0, 3, 2, 1).reshape(2 * dB, 2 * dB)
    return np.linalg.eigvalsh(T)[0]
rng = np.random.default_rng(1)
def orth(M):
    q, _ = np.linalg.qr(M); return q[:, :2]
tests = {}
e = np.eye(N)
if n == 2:
    tests['span{|00>,|11>}'] = orth(np.stack([e[0], e[4]], 1))
    tests['span{(|01>-|10>)/r2,|22>}'] = orth(np.stack([(e[1] - e[3]) / 2 ** .5, e[8]], 1))
    tests['span{Phi3,|01>}'] = orth(np.stack([(e[0] + e[4] + e[8]) / 3 ** .5, e[1]], 1))
for t in range(4):
    tests[f'random {t}'] = orth(rng.standard_normal((N, 2)) + 1j * rng.standard_normal((N, 2)))
for name, V in tests.items():
    S = sigma_of(V)
    print(f"n = {n}, {name}: PT-min {ptmin(S):+.2e}, rank {np.linalg.matrix_rank(S, 1e-10)}, "
          f"bosonic {K}-ext+PPT: {extendible(S, K)}, inner eps {inner_eps(S):+.2e}", flush=True)
