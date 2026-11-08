"""Adversarial searches over MIXED real states rho on n qutrit copies (floating point, L-BFGS on rho = G G^T / Tr).

Write Phi_k = P_{3/4}^{(x)k}, tau_s = least eigenvalue of the one-copy reduced state of rho at copy s, T_k = prod (1/4 + tau/2)
over the copies concerned, and vertex_s(rho) = (P_{1/2} at copy s (x) P_{3/4} elsewhere)(rho).  Modes (the printed number is
the ratio left side / right side; a value below 1 is a violation):

  tc        lambda_min(Phi_n(rho)) >= T_n                                   (Conjecture TC for mixed inputs)
  peel1     lambda_min(Phi_n(rho)) >= (1/4 + tau_1/2) lambda_min(Phi_{n-1}(Tr_1 rho))
  peelbest  the same for SOME copy s (ratio against the smallest right side over s)
  vtc1      lambda_min(vertex_1(rho)) >= (tau_1/2) T_{n-1}(Tr_1 rho)
  vtcbest   the same for SOME copy s (ratio = largest left/right over s)

Remark 7.31: mixed TC at n follows from mixed TC at n - 1 and vtcbest at n.
usage: python peel_search.py n rank mode starts seed"""
import sys
import numpy as np
import torch
from scipy.optimize import minimize

torch.set_default_dtype(torch.float64)
n, rank, mode, starts, seed = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5])


def pmap(X, m, avec):
    """(P_{a_1} (x) ... (x) P_{a_m})(X) with P_a(X) = a Tr(X) I - (X + X^T)/4"""
    if m == 0:
        return X.reshape(1, 1)
    T = X.reshape([3] * (2 * m))
    for s in range(m):
        T = torch.movedim(T, (s, m + s), (0, 1))
        sh = T.shape
        Mi = T.reshape(3, 3, -1)
        tr = torch.einsum('iir->r', Mi)
        T = torch.movedim((-0.25 * (Mi + Mi.transpose(0, 1)) + avec[s] * torch.eye(3).unsqueeze(-1) * tr).reshape(sh),
                          (0, 1), (s, m + s))
    return T.reshape(3 ** m, 3 ** m)


def ptrace(rho, m, s):
    T = rho.reshape([3] * (2 * m))
    return torch.diagonal(T, dim1=s, dim2=m + s).sum(-1).reshape(3 ** (m - 1), 3 ** (m - 1))


def marg(rho, m, s):
    T = rho.reshape([3] * (2 * m))
    letters = 'abcdefghijklmnop'
    idx = [None] * (2 * m)
    for i in range(m):
        if i == s:
            idx[i], idx[m + i] = 'X', 'Y'
        else:
            idx[i] = idx[m + i] = letters[i]
    return torch.einsum(''.join(idx) + '->XY', T)


def permute_first(rho, m, s):
    perm = [s] + [t for t in range(m) if t != s]
    return rho.reshape([3] * (2 * m)).permute(perm + [m + p for p in perm]).reshape(3 ** m, 3 ** m)


def lmin(M):
    return torch.linalg.eigvalsh((M + M.T) / 2)[0]


A34 = [0.75] * n


def ratio(rho):
    taus = torch.stack([torch.linalg.eigvalsh(marg(rho, n, s))[0] for s in range(n)])
    lam = lmin(pmap(rho, n, A34))
    if mode == 'tc':
        return lam / torch.prod(0.25 + 0.5 * taus)
    if mode in ('peel1', 'peelbest'):
        cand = [(0.25 + 0.5 * taus[s]) * lmin(pmap(ptrace(rho, n, s), n - 1, A34[1:])) for s in range(n)]
        return lam / (cand[0] if mode == 'peel1' else torch.min(torch.stack(cand)))
    rs = []
    for s in range(n):
        R = permute_first(rho, n, s)
        v = lmin(pmap(R, n, [0.5] + A34[1:]))
        b = 0.5 * taus[s] * torch.prod(torch.stack([0.25 + 0.5 * taus[t] for t in range(n) if t != s]))
        rs.append(v / b)
    return rs[0] if mode == 'vtc1' else torch.max(torch.stack(rs))


def obj(z):
    zt = torch.tensor(z, requires_grad=True)
    G = zt.reshape(3 ** n, rank)
    rho = G @ G.T
    rho = rho / torch.trace(rho)
    val = ratio(rho)
    val.backward()
    return val.item(), zt.grad.numpy()


rng = np.random.default_rng(seed)
best = np.inf
for st in range(starts):
    kind = st % 4
    z0 = rng.standard_normal(3 ** n * rank)
    if kind == 1:      # near a pure product
        e = np.zeros(3 ** n)
        e[0] = 1
        z0 = (np.outer(e, np.ones(rank)) + 0.2 * rng.standard_normal((3 ** n, rank))).reshape(-1)
    elif kind == 2:    # near a planar maximally entangled pair on copies 1, 2
        v = np.zeros((3,) * n)
        v[(0, 0) + (0,) * (n - 2)] = 1
        v[(1, 1) + (0,) * (n - 2)] = 1
        z0 = (np.outer(v.reshape(-1), np.ones(rank)) + 0.2 * rng.standard_normal((3 ** n, rank))).reshape(-1)
    elif kind == 3:    # nearly thin at the first copy
        z0 = z0.reshape(3, 3 ** (n - 1), rank)
        z0[2] *= 0.05
        z0 = z0.reshape(-1)
    try:
        r = minimize(obj, z0, jac=True, method='L-BFGS-B', options={'maxiter': 3000})
    except Exception as ex:
        print(f'start {st}: failed ({ex})', flush=True)
        continue
    print(f'start {st} (kind {kind}): ratio {r.fun:.10f}', flush=True)
    best = min(best, r.fun)
print(f'n = {n}, rank = {rank}, mode = {mode}: smallest ratio found {best:.10f}', flush=True)
