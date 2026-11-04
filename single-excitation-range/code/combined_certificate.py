"""The combined certificate of Conjecture 6 (single-excitation-range/README.md) as a semidefinite program.

Gauge sum_t U_t^* U_t = I, A_t = U_t^* U_t, sector c, h = n - c, co-range vectors f_sigma in C^2 (f-convention):
    Q_c(f) = sum_sigma f^*(2I + 2A_H) f - sum_rho [2 sum_t |tau_{rho,t}|^2 - |sum_t tau_{rho,t}|^2],
    tau_{rho,t} = sum_{z,r} U_t[z,r] f_{rho+z@t}[r],   y_t = conj(vec U_t) in C^4 (index (z,r)).
Combined certificate: theta_Q >= 0 for (c+1)-sets Q, Hermitian 4x4 m_{S,t} on letter(t) (x) reference, with
  (L) for every (c-1)-set R:  (+)_t m_{R+t,t} - Y(N_R) >= 0,  N_R = 2I - J - sum_{s<t} theta_{R+s+t} (e_s-e_t)(e_s-e_t)^T,
      Y(N)_{st} = N_{st} |y_s><y_t|;
  (P) for every c-set S:  sum_{t in S} m_{S,t}^{(t)} + 2 sum_{u not in S} theta_{S+u} I (x) A_{S+u} <= I (x) (2I + 2A_H).
Then Q_c >= 0 (Proposition 6).  certificate() maximises the common slack lam with cvxpy/CLARABEL and re-checks every
LMI of the returned solution with numpy eigenvalues; sector_matrix() is the explicit Hermitian matrix of Q_c.
"""
import numpy as np, cvxpy as cp, itertools, warnings
warnings.filterwarnings('ignore')

def gauge(X):
    S = sum(x.conj().T @ x for x in X); w, V = np.linalg.eigh(S)
    return [x @ (V @ np.diag(w ** -0.5) @ V.conj().T) for x in X]

def sector_matrix(U, n, c):
    S = [s for s in itertools.product(range(3), repeat=n) if sum(x != 0 for x in s) == c]; pos = {s: i for i, s in enumerate(S)}; m = len(S)
    A = [u.conj().T @ u for u in U]; M = np.zeros((2 * m, 2 * m), complex)
    for j, s in enumerate(S):
        M[2 * j:2 * j + 2, 2 * j:2 * j + 2] += 2 * np.eye(2) + 2 * sum((A[t] for t in range(n) if s[t] == 0), np.zeros((2, 2)))
    for rho in itertools.product(range(3), repeat=n):
        if sum(x != 0 for x in rho) != c - 1: continue
        free = [t for t in range(n) if rho[t] == 0]; Vm = []
        for t in free:
            v = np.zeros(2 * m, complex)
            for z in (1, 2):
                sp = list(rho); sp[t] = z; v[2 * pos[tuple(sp)]:2 * pos[tuple(sp)] + 2] += U[t][z - 1, :]
            Vm.append(v)
        Vm = np.array(Vm); k = len(free); M -= Vm.T @ (2 * np.eye(k) - np.ones((k, k))) @ Vm.conj()
    return (M + M.conj().T) / 2

def perm_to_front(c, pos):
    D = 2 ** (c + 1); idx = np.arange(D).reshape([2] * (c + 1)); order = [pos, c] + [i for i in range(c) if i != pos]
    new = np.transpose(idx, order).reshape(D); P = np.zeros((D, D)); P[np.arange(D), new] = 1.0; return P

def certificate(U, n, c, local_only=False):
    A = [u.conj().T @ u for u in U]; Qs = list(itertools.combinations(range(n), c + 1))
    th = {Q: (0.0 if local_only else cp.Variable(nonneg=True)) for Q in Qs}; m = {}; cons = []; Pm = {p: perm_to_front(c, p) for p in range(c)}
    ys = [np.conj(np.array([U[t][z, r] for z in range(2) for r in range(2)])) for t in range(n)]
    lam = cp.Variable(); P_exprs = {}; L_exprs = {}
    for S in itertools.combinations(range(n), c):
        H = [x for x in range(n) if x not in S]; budget = np.kron(np.eye(2 ** c), 2 * np.eye(2) + 2 * sum((A[x] for x in H), np.zeros((2, 2))))
        tot = 0
        for p, t in enumerate(S):
            m[S, t] = cp.Variable((4, 4), hermitian=True); tot = tot + Pm[p].T @ cp.kron(m[S, t], np.eye(2 ** (c - 1))) @ Pm[p]
        for u in H:
            Q = tuple(sorted(S + (u,))); tot = tot + 2 * th[Q] * np.kron(np.eye(2 ** c), sum(A[x] for x in Q))
        P_exprs[S] = budget - tot; cons.append(P_exprs[S] - lam * np.eye(2 ** (c + 1)) >> 0)
    for R in itertools.combinations(range(n), c - 1):
        free = [t for t in range(n) if t not in R]; k = len(free); blocks = [[None] * k for _ in range(k)]
        for i in range(k):
            for j in range(k):
                if i == j: Nij = 1.0 - sum(th[tuple(sorted(R + (free[i], free[l])))] for l in range(k) if l != i)
                else: Nij = -1.0 + th[tuple(sorted(R + (free[i], free[j])))]
                blocks[i][j] = (m[tuple(sorted(R + (free[i],))), free[i]] if i == j else 0) - Nij * np.outer(ys[free[i]], ys[free[j]].conj())
        L_exprs[R] = cp.bmat(blocks); cons.append(L_exprs[R] - lam * np.eye(4 * k) >> 0)
    prob = cp.Problem(cp.Maximize(lam), cons); prob.solve(solver='CLARABEL')
    # (b) independent re-check of every LMI with the returned values
    worst = min([np.linalg.eigvalsh((E.value + E.value.conj().T) / 2)[0] for E in P_exprs.values()] +
                [np.linalg.eigvalsh((E.value + E.value.conj().T) / 2)[0] for E in L_exprs.values()])
    thmin = min([0.0] if local_only else [th[Q].value for Q in Qs])
    return lam.value, worst, thmin

