"""Check of the Bethe ansatz for E_c on letters 0, 1 (exploration): fit Bethe roots to each eigenvector, then test
the Bethe equations  sum_t s_t/(s_t - x_a) + 2 sum_{b != a} x_b/(x_a - x_b) = 4  and the eigenvalue 2|s| - 2 sum x_a."""
import itertools
import numpy as np
from scipy.optimize import least_squares
from l1float import L1Block
rng = np.random.default_rng(7)
for (k, c) in [(6, 2), (7, 3)]:
    s = np.exp(rng.normal(0, 0.7, size=k))
    B = L1Block(list(s), c)
    E = B.E(c)
    ev, V = np.linalg.eigh(E)
    H = [tuple(sorted(y)) for y in B.H[c]]
    w = np.sqrt(s)
    perms = list(itertools.permutations(range(c)))
    def bethe(x):
        xs = [x[2 * m] + 1j * x[2 * m + 1] for m in range(c)]
        F = [w / (s - xm) for xm in xs]
        return np.array([sum(np.prod([F[p[m]][y[m]] for m in range(c)]) for p in perms) for y in H]), xs
    worst_eq, worst_ev, nfit = 0.0, 0.0, 0
    for n in range(len(ev)):
        v = V[:, n]; best = (1.0, None)
        for trial in range(80):
            x0 = rng.normal(0, 3, 2 * c)
            def r(x):
                b, _ = bethe(x)
                cc = np.vdot(b, v) / np.vdot(b, b)
                d = v - cc * b
                return np.concatenate([d.real, d.imag])
            sol = least_squares(r, x0, max_nfev=800, xtol=1e-15, ftol=1e-15, gtol=1e-15)
            res = np.linalg.norm(sol.fun)
            if res < best[0]:
                best = (res, sol.x)
            if res < 1e-11:
                break
        if best[0] > 1e-8:
            continue
        nfit += 1
        _, xs = bethe(best[1])
        for a in range(c):
            lhs = sum(s / (s - xs[a])) + 2 * sum(xs[b] / (xs[a] - xs[b]) for b in range(c) if b != a)
            worst_eq = max(worst_eq, abs(lhs - 4))
        worst_ev = max(worst_ev, abs(2 * s.sum() - 2 * sum(xs) - ev[n]))
    print(f'k={k} c={c}: {nfit} of {len(ev)} eigenvectors fitted; max |Bethe equation - 4| = {worst_eq:.2e}; '
          f'max |2|s| - 2 sum x - eigenvalue| = {worst_ev:.2e}')
