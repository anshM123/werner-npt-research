"""Floating-point exploration behind the remarks of NOTES.md, Section 5 (no claim rests on it).  Double precision,
closed forms of nlib.py, copy t = 0, |s| = 1.

1. Ratio X22 : Pi on w-perp (Proposition Q gives >= (2 + 2gp + gp^2)/(2 + gp)^2 >= 1/2).
2. H-matrix tests (comparison matrix positive definite?) for Phi, Phi' and X (gauge-fixed, rescaled).
3. Block norm-comparison test for X with 2 x 2 blocks (eta_c, pi_c).
4. Decomposition over pairs of the state groups R, P, Q (instead of pairs of copies): SDP margin.
5. Smallest normalised eigenvalue of the (R, Q) block C and where it occurs.
Usage: python explore_structure.py"""
import itertools
import numpy as np
import cvxpy as cp
from nlib import blocks, phi, phip, full

rng = np.random.default_rng(2026)


def sample(n):
    out = []
    for i in range(n):
        k = int(rng.integers(3, 9))
        m = i % 4
        if m == 0:
            x = rng.uniform(0, 1, k)
        elif m == 1:
            x = rng.dirichlet(0.2 * np.ones(k))
        elif m == 2:
            x = 3.0 ** rng.integers(-6, 3, k)
        else:
            x = rng.dirichlet(0.3 * np.ones(k)); x[0] = 10 ** rng.uniform(-4, 0)
        x = np.maximum(x, 1e-6)
        out.append(x / x.sum())
    return out


def nmin(A):
    d = np.sqrt(np.abs(np.diag(A)))
    return np.linalg.eigvalsh(A / d[:, None] / d[None, :]).min()


def comp(A):
    C = -np.abs(A)
    np.fill_diagonal(C, np.abs(np.diag(A)))
    return C


def Xg(b):
    k = b['k']; w = b['w']; A1 = b['A1']; A2 = b['A2']; H = b['H']; N = b['N']
    X11 = A1 @ np.diag(N) @ A1 - 0.5 * H @ A2 @ H.T
    X12 = A1 @ b['YRP'] + 0.5 * H @ A2 @ b['YPQ'].T
    X22 = b['Pi'] - 0.5 * b['YPQ'] @ A2 @ b['YPQ'].T
    X = np.block([[X11, X12], [X12.T, X22]])
    D = np.concatenate([w, w])
    Xt = X * D[:, None] * D[None, :]
    keep = [i for i in range(2 * k) if i != k]          # delete pi_t
    return Xt[np.ix_(keep, keep)], X22


pts = sample(1500)
r1 = 9
cnt = {'Phi': 0, 'Phip': 0, 'X': 0, 'Xblock': 0}
for x in pts:
    b = blocks(x); k = b['k']; w = b['w']
    G, X22 = Xg(b)
    Q, _ = np.linalg.qr(np.column_stack([w, rng.standard_normal((k, k - 1))]))
    B = Q[:, 1:]
    L = np.linalg.cholesky(B.T @ b['Pi'] @ B); Li = np.linalg.inv(L)
    r1 = min(r1, np.linalg.eigvalsh(Li @ (B.T @ X22 @ B) @ Li.T).min())
    cnt['Phi'] += nmin(comp(phi(b))) <= 0
    cnt['Phip'] += nmin(comp(phip(b))) <= 0
    cnt['X'] += nmin(comp(G)) <= 0
    # block norm comparison, blocks {eta_t}, {eta_c, pi_c}
    blk = [[0]] + [[c, k - 1 + c] for c in range(1, k)]
    W = []
    for g in blk:
        e, V = np.linalg.eigh(G[np.ix_(g, g)])
        W.append(V @ np.diag(1 / np.sqrt(e)) @ V.T)
    Nm = np.eye(k)
    for i in range(k):
        for j in range(k):
            if i != j:
                Nm[i, j] = -np.linalg.norm(W[i] @ G[np.ix_(blk[i], blk[j])] @ W[j], 2)
    cnt['Xblock'] += np.linalg.eigvalsh(Nm).min() <= 0
print(f'1. smallest ratio X22 : Pi on w-perp over {len(pts)} samples: {r1:.4f}')
print(f'2. comparison matrix not positive definite: Phi {cnt["Phi"]}, Phi\' {cnt["Phip"]}, X (gauge-fixed) {cnt["X"]} '
      f'of {len(pts)}')
print(f'3. block norm-comparison matrix of X not positive definite: {cnt["Xblock"]} of {len(pts)}')


def cone(Gm, groups):
    n = Gm.shape[0]
    tot = 0
    for a, c in itertools.combinations(range(len(groups)), 2):
        idx = groups[a] + groups[c]
        Z = cp.Variable((len(idx), len(idx)), PSD=True)
        E = np.zeros((len(idx), n))
        for i, j in enumerate(idx):
            E[i, j] = 1
        tot = tot + E.T @ Z @ E
    eps = cp.Variable()
    cp.Problem(cp.Maximize(eps), [tot + eps * np.eye(n) == Gm, eps <= 1]).solve(solver=cp.SCS, eps=1e-9, max_iters=50000)
    return float(eps.value)


neg = 0
tests = pts[:10] + [np.array([0.01, 0.495, 0.495]), np.array([0.99, 0.0098, 0.0002]), np.array([0.3, 0.3, 0.2, 0.2])]
for x in tests:
    b = blocks(x); k = b['k']
    M = full(b)
    idx = [i for i in range(3 * k - 1) if i != k]
    Gm = M[np.ix_(idx, idx)]
    d = np.sqrt(np.diag(Gm)); Gm = Gm / d[:, None] / d[None, :]
    groups = [list(range(0, k)), list(range(k, 2 * k - 1)), list(range(2 * k - 1, 3 * k - 2))]
    neg += cone(Gm, groups) < 0
print(f'4. decomposition over pairs of the state groups R, P, Q: SDP margin negative at {neg} of {len(tests)} points')
worst = []
for x in pts:
    b = blocks(x)
    C = np.block([[b['DR'], b['YRQ']], [b['YRQ'].T, b['YQQ']]])
    worst.append((nmin(C), x))
worst.sort(key=lambda r: r[0])
print('5. smallest normalised eigenvalues of C:', '; '.join(f'{v:.2e} at {np.round(x, 4)}' for v, x in worst[:3]))
# 6. the undeleted block (group of t = {R_t, P_t}) and the inverse-free congruent form diag(A1, I, A2)^T M diag(A1, I, A2)
#    (both singular): SDP margin, which is 0 up to solver accuracy exactly when the matrix is in the closed cone
mx = 0
negs = []
for x in pts[10:22]:
    b = blocks(x); k = b['k']; w = b['w']
    for cong in (False, True):
        M = full(b)
        if cong:
            T = np.zeros((3 * k - 1, 3 * k - 1)); T[:k, :k] = b['A1']; T[k:2 * k, k:2 * k] = np.eye(k); T[2 * k:, 2 * k:] = b['A2']
            M = T.T @ M @ T
        mu = np.concatenate([w, w, w[1:]]); M = M * mu[:, None] * mu[None, :]
        d = np.sqrt(np.diag(M)); M = M / d[:, None] / d[None, :]
        groups = [[0, k]] + [[c, k + c, 2 * k + c - 1] for c in range(1, k)]
        e = cone(M, groups)
        if cong:
            negs.append(e)
        else:
            mx = max(mx, abs(e))
print(f'6. undeleted block, 12 points: largest |SDP margin| {mx:.1e}; inverse-free congruent form, same points: margin '
      f'below -1e-6 at {sum(1 for e in negs if e < -1e-6)} of 12 (smallest {min(negs):.2e})')
print('numerical exploration done')
