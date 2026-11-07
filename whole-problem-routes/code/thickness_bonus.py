"""Checks for Proposition 7.28 (the thickness bonus at two copies) and the one-sided form of Conjecture TC.

Proposition 7.28: for a real unit y on two copies whose one-copy reduced state has least eigenvalue tau,
    lambda_min(P_{3/4}^{(x)2}(y y^T)) >= min{1/16 + (3/4) tau, 1/8 + (7/16) tau} >= (1/4 + tau/2)^2.
1. With y = sum_i s_i e_i (x) e_i (local rotations), P^{(x)2}(yy^T) is block diagonal (Theorem 7.2): the block
   diag(a^2 - (a - 1/8) s_k^2) + (1/8) s s^T on span{|kk>}, and for i < j a 2x2 block with eigenvalues
   a^2 - (a/2)(s_i^2 + s_j^2) +- (1/8) s_i s_j.  On a fine grid of the Schmidt simplex: the block formula agrees with
   the least eigenvalue of the explicit 9 x 9 matrix, and both inequalities hold (least gaps reported).
2. The largest kappa with lambda_min >= (1/4 + kappa tau)^2 at two copies (it is attained at the maximally entangled y).
3. One-sided form at three and four copies: adversarial minimisation of lambda_min(P^{(x)n}(yy^T)) - prod_s (1/4 + tau_s/2)
   over real unit y (random and near-product starts); it should not go below zero.
Ends with PASS or FAIL."""
import numpy as np
import torch
from scipy.optimize import minimize

a = 0.75
ok = True


def pmap_np(X, m):
    T = X.reshape([3] * (2 * m))
    for s in range(m):
        T = np.moveaxis(T, [s, m + s], [0, 1])
        sh = T.shape
        Mi = T.reshape(3, 3, -1)
        tr = np.einsum('iir->r', Mi)
        out = -0.25 * (Mi + Mi.transpose(1, 0, 2)) + a * np.eye(3)[:, :, None] * tr
        T = np.moveaxis(out.reshape(sh), [0, 1], [s, m + s])
    return T.reshape(3**m, 3**m)


def block_min(s):
    MA = np.diag(a**2 - (a - 1 / 8) * s**2) + np.outer(s, s) / 8
    vals = [np.linalg.eigvalsh(MA)[0]]
    for i in range(3):
        for j in range(i + 1, 3):
            vals.append(a**2 - (a / 2) * (s[i]**2 + s[j]**2) - s[i] * s[j] / 8)
    return min(vals)


# ------------------------------------------------------------------ parts 1 and 2
worst_formula = 0.0
gap1 = np.inf
gap2 = np.inf
kap = np.inf
for u in np.linspace(0, np.pi / 2, 241):
    for v in np.linspace(0, np.pi / 2, 241):
        s = np.array([np.sin(u) * np.cos(v), np.sin(u) * np.sin(v), np.cos(u)])
        tau = float(np.min(s**2))
        bm = block_min(s)
        if (int(u * 1000) + int(v * 1000)) % 7 == 0:          # spot-check the block formula on part of the grid
            y = np.zeros((3, 3))
            np.fill_diagonal(y, s)
            y = y.reshape(-1)
            W = pmap_np(np.outer(y, y), 2)
            worst_formula = max(worst_formula, abs(np.linalg.eigvalsh((W + W.T) / 2)[0] - bm))
        gap1 = min(gap1, bm - min(1 / 16 + 0.75 * tau, 1 / 8 + 7 * tau / 16))
        gap2 = min(gap2, min(1 / 16 + 0.75 * tau, 1 / 8 + 7 * tau / 16) - (0.25 + tau / 2)**2)
        if tau > 1e-4:
            kap = min(kap, (np.sqrt(bm) - 0.25) / tau)
p1 = worst_formula < 1e-12 and gap1 > -1e-12 and gap2 > -1e-12
print(f"Part 1: block formula vs 9x9 eigenvalues: max deviation {worst_formula:.1e}; least of lambda_min - min(...) "
      f"{gap1:+.2e}; least of min(...) - (1/4 + tau/2)^2 {gap2:+.2e}{'' if p1 else '  FAIL'}")
ok &= bool(p1)
s_me = np.ones(3) / np.sqrt(3)
k_me = (np.sqrt(block_min(s_me)) - 0.25) / (1 / 3)
p2 = abs(kap - k_me) < 1e-3
print(f"Part 2: largest kappa with lambda_min >= (1/4 + kappa tau)^2 at two copies: {kap:.4f} (maximally entangled y: "
      f"{k_me:.4f}, lambda_min = {block_min(s_me):.6f} = 13/48){'' if p2 else '  FAIL'}")
ok &= bool(p2)

# ------------------------------------------------------------------ part 3: one-sided form at n = 3, 4
torch.set_default_dtype(torch.float64)


def pmap(X, m):
    T = X.reshape([3] * (2 * m))
    for s in range(m):
        T = torch.movedim(T, (s, m + s), (0, 1))
        sh = T.shape
        Mi = T.reshape(3, 3, -1)
        tr = torch.einsum('iir->r', Mi)
        out = -0.25 * (Mi + Mi.transpose(0, 1)) + a * torch.eye(3).unsqueeze(-1) * tr
        T = torch.movedim(out.reshape(sh), (0, 1), (s, m + s))
    return T.reshape(3**m, 3**m)


def margmin(v, n):
    rows = []
    for s in range(n):
        V = torch.movedim(v.reshape([3] * n), s, 0).reshape(3, -1)
        rows.append(torch.linalg.eigvalsh(V @ V.T)[0])
    return torch.stack(rows)


rng = np.random.default_rng(728)
for n, R in ((3, 12), (4, 9)):
    D = 3**n
    e1 = np.zeros(3)
    e1[0] = 1
    prod = e1
    for _ in range(n - 1):
        prod = np.kron(prod, e1)
    worst = np.inf
    for st in range(R):
        y0 = rng.standard_normal(D) if st % 3 else prod + 0.15 * rng.standard_normal(D)

        def obj(z):
            zt = torch.tensor(z, requires_grad=True)
            y = zt / zt.norm()
            W = pmap(torch.outer(y, y), n)
            lm = torch.linalg.eigvalsh((W + W.T) / 2)[0]
            tau = margmin(y, n).clamp(min=0)
            val = lm - torch.prod(0.25 + 0.5 * tau)
            val.backward()
            return val.item(), zt.grad.numpy()
        worst = min(worst, minimize(obj, y0, jac=True, method='L-BFGS-B', options={'maxiter': 3000}).fun)
    print(f"Part 3: n = {n}: least lambda_min(P(yy^T)) - prod(1/4 + tau_s/2) found {worst:+.2e}", flush=True)
    ok &= bool(worst > -1e-9)
print("PASS" if ok else "FAIL")
