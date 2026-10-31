"""Section 11: no functional W with W (x) rho decomposable (W (x) rho = P + Q^Gamma, P, Q >= 0, which implies W in
K_1(3)^*) separates the Tiles state from K_1(3).  SDP: minimise <W, T> over Tr W = 1 and W (x) rho decomposable, for T
the Tiles state and Tiles with 5% noise; the optimum is 0 (attained by a positive W).
usage: python check_tiles_decomposable.py   (needs cvxpy with Clarabel; a few minutes)
"""
import numpy as np, cvxpy as cp
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2
e0, e1, e2 = np.eye(3)
upb = [np.kron(e0, (e0 - e1) / 2 ** .5), np.kron(e2, (e1 - e2) / 2 ** .5), np.kron((e0 - e1) / 2 ** .5, e2),
       np.kron((e1 - e2) / 2 ** .5, e0), np.kron((e0 + e1 + e2) / 3 ** .5, (e0 + e1 + e2) / 3 ** .5)]
Tiles = (np.eye(9) - sum(np.outer(x, x) for x in upb)) / 4
def arrange_np(ops):
    """ops[k] on pair k (each 3 x 3 local dims): operator on (A_0 A_1 ..)(B_0 B_1 ..)."""
    m = len(ops); R = np.eye(1)
    for o in ops: R = np.kron(R, o)
    perm = [2 * i for i in range(m)] + [2 * i + 1 for i in range(m)]
    R = R.reshape([3] * (4 * m)).transpose(perm + [2 * m + q for q in perm])
    return R.reshape(9 ** m, 9 ** m)
def arrange_cp_W(W, extra):
    """W (cvxpy 9x9 on (Aw Bw)) tensored with fixed operator 'extra' on the remaining pairs, reordered."""
    m = 1 + int(round(np.log(extra.shape[0]) / np.log(9)))
    # build linear map: sum over basis of W
    D = 9 ** m; out = 0
    for i in range(9):
        for j in range(9):
            E = np.zeros((9, 9)); E[i, j] = 1
            B = arrange_np_from_pairop(np.kron(E, extra), m)
            out = out + W[i, j] * B
    return out
def arrange_np_from_pairop(R, m):
    perm = [2 * i for i in range(m)] + [2 * i + 1 for i in range(m)]
    R = R.reshape([3] * (4 * m)).transpose(perm + [2 * m + q for q in perm])
    return R.reshape(9 ** m, 9 ** m)
def ptB(X, d):
    return cp.partial_transpose(X, dims=(d, d), axis=1)
def best_witness(target):
    W = cp.Variable((9, 9), hermitian=True)
    WR = arrange_cp_W(W, rho)                      # 81 x 81 with Alice = (Aw, A1), Bob = (Bw, B1)
    P = cp.Variable((81, 81), hermitian=True); Q = cp.Variable((81, 81), hermitian=True)
    cons = [cp.real(cp.trace(W)) == 1, P >> 0, Q >> 0, WR == P + ptB(Q, 9)]
    prob = cp.Problem(cp.Minimize(cp.real(cp.trace(W @ target))), cons)
    prob.solve(solver=cp.CLARABEL)
    return prob.value, (W.value + W.value.conj().T) / 2
for label, target in [("Tiles", Tiles), ("Tiles + 5% noise", 0.95 * Tiles + 0.05 * np.eye(9) / 9)]:
    val, Wv = best_witness(target)
    print(f"[{label}] min <W, target> over W with W (x) rho decomposable: {val:+.5f}", flush=True)
    print(f"   spectrum of W: {np.round(np.linalg.eigvalsh(Wv), 4)}", flush=True)
