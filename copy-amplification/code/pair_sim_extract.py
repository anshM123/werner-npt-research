"""Exact certificate that two qutrit Werner copies are a separable-operation image of one 9-dim endpoint core.

Claim to certify at a rational beta:  rho_beta (x) rho_beta  (two qutrit sites, Alice A1A2 | Bob B1B2) lies in
    conv{ twirl[(A (x) B) E_9 (A (x) B)^T] / Tr ,  twirl[x x^T (x) y y^T] / Tr }      (E_9 = I - F/2 on C^9 (x) C^9)
with REAL filters A, B : C^9 -> C^9 and real product vectors x, y in C^9.  Twirled operators on the two sites are
determined by their four sector weights (ss, sa, as, aa), so membership is a 4-dimensional convex-hull question.
Consequence: rho_beta^{(x)4} is a nonnegative combination of separable-operation images of E_9 (x) E_9,
E_9 (x) separable and separable operators, hence one-copy undistillable by the two-copy theorem at d = 9.

Pipeline: column generation (float) -> basic LP solution (<= 4 generators) -> round filters to rationals ->
exact sector weights (Fractions) -> exact barycentric weights -> check all >= 0.
usage: pair_sim_extract.py beta_num beta_den starts steps rounds denom
"""
import itertools
import sys
from fractions import Fraction as Fr

import numpy as np
import torch
from scipy.optimize import linprog

torch.set_default_dtype(torch.float64)
bn, bd = int(sys.argv[1]), int(sys.argv[2])
beta = Fr(bn, bd)
starts, steps, rounds, denom = (int(x) for x in sys.argv[3:7])
d = 3
N = d * d            # Alice's two qutrit sites
gen = torch.Generator().manual_seed(23)


def swap_np(m):
    F = np.zeros((m * m, m * m))
    for i in range(m):
        for j in range(m):
            F[i * m + j, j * m + i] = 1
    return F


# sector projectors on (A1 A2) (x) (B1 B2) ordered as a1 a2 b1 b2
def sector_ops_np():
    Ps = (np.eye(d * d) + swap_np(d)) / 2; Pa = (np.eye(d * d) - swap_np(d)) / 2
    ops = []
    for P1 in (Ps, Pa):
        for P2 in (Ps, Pa):
            # P1 acts on (a1, b1), P2 on (a2, b2); reorder to (a1 a2 b1 b2)
            O = np.kron(P1, P2).reshape(d, d, d, d, d, d, d, d).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(N * N, N * N)
            ops.append(O)
    return ops


OPS = sector_ops_np()
OPS_t = [torch.tensor(O) for O in OPS]
E9 = np.eye(N * N) - swap_np(N) / 2
E9_t = torch.tensor(E9)
bf = float(beta)
rho = np.eye(d * d) - bf * swap_np(d)
T = np.kron(rho, rho).reshape(d, d, d, d, d, d, d, d).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(N * N, N * N)
target = np.array([np.trace(T @ O) for O in OPS]); target /= target.sum()


def prof_t(X):
    v = torch.stack([torch.trace(X @ O) for O in OPS_t])
    return v / v.sum()


def oracle_core(w):
    best = (-np.inf, None, None)
    wt = torch.tensor(w)
    for s in range(starts):
        A = torch.randn(N, N, generator=gen).requires_grad_(); B = torch.randn(N, N, generator=gen).requires_grad_()
        opt = torch.optim.Adam([A, B], lr=0.05)
        for it in range(steps):
            opt.zero_grad()
            AB = torch.kron(A, B)
            X = AB @ E9_t @ AB.T
            val = (wt * prof_t(X)).sum()
            (-val).backward(); opt.step()
        v = float(val.detach())
        if v > best[0]:
            best = (v, ('core', A.detach().numpy().copy(), B.detach().numpy().copy()), prof_t(X).detach().numpy())
    return best


def oracle_prod(w):
    best = (-np.inf, None, None)
    wt = torch.tensor(w)
    for s in range(starts):
        x = torch.randn(N, generator=gen).requires_grad_(); y = torch.randn(N, generator=gen).requires_grad_()
        opt = torch.optim.Adam([x, y], lr=0.05)
        for it in range(steps):
            opt.zero_grad()
            X = torch.kron(torch.outer(x, x), torch.outer(y, y))
            val = (wt * prof_t(X)).sum()
            (-val).backward(); opt.step()
        v = float(val.detach())
        if v > best[0]:
            best = (v, ('prod', x.detach().numpy().copy(), y.detach().numpy().copy()), prof_t(X).detach().numpy())
    return best


cols = []
for w in [np.eye(4)[k] for k in range(4)] + [np.array([-1, -1, -1, 1.0]), np.array([-1, 0, 0, 1.0]),
                                             np.array([0, 1, 1, -1.0]), np.array([-1, 0.9, -1, 1.0]),
                                             np.array([-1, -1, 0.9, 1.0])]:
    for orc in (oracle_core, oracle_prod):
        v, g, p = orc(w)
        cols.append((g, p))
print(f'beta = {beta}: target = {np.round(target, 5)}', flush=True)
for r in range(rounds):
    P = np.array([c[1] for c in cols])
    res = linprog(-np.concatenate([target, [-1.0]]), A_ub=np.hstack([P, -np.ones((len(P), 1))]),
                  b_ub=np.zeros(len(P)), bounds=[(-1, 1)] * 4 + [(None, None)], method='highs')
    w, c = res.x[:4], res.x[4]
    gap = float(target @ w - c)
    print(f'  round {r}: gap {gap:+.6f}', flush=True)
    if gap <= 1e-9:
        break
    cands = [oracle_core(w), oracle_prod(w)]
    v, g, p = max(cands, key=lambda t: t[0])
    if v <= c + 1e-9:
        print('  oracle cannot improve: NOT in hull (numerical)'); sys.exit(0)
    cols.append((g, p))

# basic feasible solution
P = np.array([c[1] for c in cols])
res = linprog(np.zeros(len(P)), A_eq=np.vstack([P.T, np.ones(len(P))]), b_eq=np.concatenate([target, [1.0]]),
              bounds=[(0, None)] * len(P), method='highs-ds')
lam = res.x
supp = [i for i in range(len(P)) if lam[i] > 1e-10]
import pickle
pickle.dump({"beta": (bn, bd), "support": [cols[i] for i in supp], "lam": lam[supp]}, open(f"pair_sim_support_{bn}_{bd}.pkl", "wb"))
print('float support:', [(cols[i][0][0], round(lam[i], 6)) for i in supp], flush=True)


# ---- exact verification ----
def to_fr(a):
    return np.vectorize(lambda v: Fr(int(round(v * denom)), denom), otypes=[object])(a)


def fr_kron(a, b):
    m, n = a.shape; p, q = b.shape
    out = np.empty((m * p, n * q), dtype=object)
    for i in range(m):
        for j in range(n):
            out[i * p:(i + 1) * p, j * q:(j + 1) * q] = a[i, j] * b
    return out


OPS_fr = [np.vectorize(Fr, otypes=[object])(np.round(O * 4).astype(int)) / 4 for O in OPS]  # entries in {0, +-1/4, +-1/2, 1}
for O, Of in zip(OPS, OPS_fr):
    assert np.allclose(O, np.vectorize(float)(Of))
E9_fr = np.vectorize(Fr, otypes=[object])(np.round(E9 * 2).astype(int)) / 2


def exact_profile(g):
    kind, a, b = g
    if kind == 'core':
        A, B = to_fr(a), to_fr(b)
        AB = fr_kron(A, B)
        X = AB.dot(E9_fr).dot(AB.T)
    else:
        x, y = to_fr(a).reshape(-1, 1), to_fr(b).reshape(-1, 1)
        X = fr_kron(x.dot(x.T), y.dot(y.T))
    v = [sum(X[i, j] * O[j, i] for i in range(X.shape[0]) for j in range(X.shape[0]) if O[j, i] != 0) for O in OPS_fr]
    s = sum(v)
    return [vi / s for vi in v]


rho_fr = np.vectorize(Fr, otypes=[object])(np.eye(d * d).astype(int)) - beta * np.vectorize(Fr, otypes=[object])(swap_np(d).astype(int))
T_fr = fr_kron(rho_fr, rho_fr)
# reorder T_fr from (a1 b1 a2 b2) to (a1 a2 b1 b2)
idx = np.arange(N * N).reshape(d, d, d, d).transpose(0, 2, 1, 3).reshape(-1)
T_fr = T_fr[np.ix_(idx, idx)]
tv = [sum(T_fr[i, j] * O[j, i] for i in range(N * N) for j in range(N * N) if O[j, i] != 0) for O in OPS_fr]
ts = sum(tv)
target_fr = [x / ts for x in tv]
print('exact target:', [str(x) for x in target_fr], flush=True)
assert np.allclose([float(x) for x in target_fr], target), 'exact target mismatch'
profs = [exact_profile(cols[i][0]) for i in supp]
k = len(profs)
# solve sum lam_i prof_i = target, sum lam_i = 1 (least-squares exact via normal equations on chosen rows)
import sympy as sp
M = sp.Matrix([[profs[i][r] for i in range(k)] for r in range(4)] + [[1] * k])
rhs = sp.Matrix([target_fr[r] for r in range(4)] + [1])
sol = M.gauss_jordan_solve(rhs)
lam_exact = sol[0]
free = sol[1]
if free.shape[0] > 0:
    lam_exact = lam_exact.subs({f: 0 for f in free})
ok = all(l >= 0 for l in lam_exact) and (M * lam_exact - rhs).is_zero_matrix
print('exact weights:', [str(l) for l in lam_exact])
print('EXACT CERTIFICATE VALID:', ok)
