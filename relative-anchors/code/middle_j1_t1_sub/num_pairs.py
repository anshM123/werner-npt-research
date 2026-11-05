"""Floating-point survey for Conjecture TC (NOTES.md, Section 4); no claim rests on it.

For random and degenerate weight vectors (k = 3, ..., 10) the gauge-fixed block G (the state P_t deleted, diagonally
normalised; closed forms of nlib.py, double precision) is tested for the pair-cone property: the largest eps such that
G - eps I is a sum of positive semidefinite matrices each supported on the states of two copies (groups {R_t} and
{R_c, P_c, Q_c}).  Reported: eps, the smallest eigenvalue lambda_min of G, and their ratio (eps <= lambda_min always).
Controls: random positive definite matrices of the same sizes, normalised the same way (Wishart with n degrees of
freedom plus 0.02 I), are not in the cone (eps < 0).
Weights below 1e-6 are raised to 1e-6: at x_t of order 1e-7 the smallest eigenvalue is about 6e-8, below the accuracy
of the solver, and the survey once returned eps = -5e-8 there; at the same shape with x_t = 1e-4, ..., 1e-6 the ratio
eps/lambda_min is 1.000 (and cert_pairs.py certifies a rational point of that shape exactly).
PASS: eps > 0 at every weight vector and every control outside the cone.
Usage: python num_pairs.py [n_random]"""
import sys
import numpy as np
import cvxpy as cp
import itertools
from nlib import blocks, full

STATUS = []      # solver status of every solve (SCS reports optimal_inaccurate when it stops at reduced accuracy)


def nrm(A):
    d = np.sqrt(np.diag(A))
    return A / d[:, None] / d[None, :]


def gauge(x):
    b = blocks(x)
    k = b['k']
    M = full(b)
    idx = [i for i in range(3 * k - 1) if i != k]
    G = nrm(M[np.ix_(idx, idx)])
    groups = [[0]] + [[c, k - 1 + c, 2 * k - 2 + c] for c in range(1, k)]
    return G, groups


def cone_eps(G, groups):
    n = G.shape[0]
    tot = 0
    for a, b in itertools.combinations(range(len(groups)), 2):
        idx = groups[a] + groups[b]
        Z = cp.Variable((len(idx), len(idx)), PSD=True)
        E = np.zeros((len(idx), n))
        for i, j in enumerate(idx):
            E[i, j] = 1
        tot = tot + E.T @ Z @ E
    eps = cp.Variable()
    prob = cp.Problem(cp.Maximize(eps), [tot + eps * np.eye(n) == G, eps <= 1])
    prob.solve(solver=cp.SCS, eps=1e-9, max_iters=100000)
    STATUS.append(prob.status)
    return float(eps.value) if eps.value is not None else float('nan')


if __name__ == '__main__':
    nr = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    rng = np.random.default_rng(17)
    tests = []
    for i in range(nr):
        k = int(rng.integers(3, 11))
        m = i % 5
        if m == 0:
            x = rng.uniform(0, 1, k)
        elif m == 1:
            x = rng.dirichlet(0.2 * np.ones(k))
        elif m == 2:
            x = 3.0 ** rng.integers(-6, 3, k)
        elif m == 3:
            x = rng.dirichlet(0.3 * np.ones(k)); x[0] = 10 ** rng.uniform(-4, -1)
        else:
            x = rng.dirichlet(0.3 * np.ones(k)); x[0] = 1 - 10 ** rng.uniform(-3, -0.5)
            x[1:] *= (1 - x[0]) / x[1:].sum()
        x = np.maximum(x, 1e-6)
        tests.append(x / x.sum())
    for s in (1e-2, 1e-3):
        tests.append(np.array([s, (1 - s) / 2, (1 - s) / 2]))
        tests.append(np.array([s, 0.45, 0.45] + [(0.1 - s) / 5] * 5))
    for s in (0.99, 0.999):
        tests.append(np.array([s, (1 - s) * 0.98, (1 - s) * 0.02]))
        tests.append(np.array([s] + [(1 - s) / 8] * 8))
    worst = 9
    bad = 0
    for x in tests:
        G, groups = gauge(x)
        e = cone_eps(G, groups)
        lmin = np.linalg.eigvalsh(G).min()
        worst = min(worst, e / lmin)
        bad += (not e > 0)
        print(f'k = {len(x)}, x = {np.round(x, 4)}: eps = {e:.3e}, lambda_min = {lmin:.3e}, ratio {e / lmin:.4f} '
              f'[{STATUS[-1]}]', flush=True)
    print(f'{len(tests)} weight vectors: eps <= 0 at {bad}; smallest ratio eps/lambda_min = {worst:.4f}; solver status '
          f'optimal_inaccurate at {sum(1 for st in STATUS if st != "optimal")}', flush=True)
    cbad = 0
    for k in (4, 6, 8):
        n = 3 * (k - 1) + 1
        groups = [[0]] + [[1 + 3 * (c - 1), 2 + 3 * (c - 1), 3 + 3 * (c - 1)] for c in range(1, k)]
        for trial in range(3):
            A = rng.standard_normal((n, n))
            G = nrm(A @ A.T / n + 0.02 * np.eye(n))
            e = cone_eps(G, groups)
            cbad += (e < 0)
            print(f'  control k = {k}: random positive definite matrix, eps = {e:.3e}, lambda_min = '
                  f'{np.linalg.eigvalsh(G).min():.3e}', flush=True)
    print(f'controls outside the cone: {cbad} of 9', flush=True)
    print('PASS' if (bad == 0 and cbad == 9) else 'FAIL')
