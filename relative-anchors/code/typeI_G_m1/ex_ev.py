"""Floating-point exploration (no claim rests on it): where does the normalised margin of the rho-row of (D),
f(s) = (G_z v)_rho / (G_z[rho,rho] v_0), attain its minimum over the U-weights with sigma and S fixed?  (A test of the
"equal variable" route of NOTES.md Section 4.)  u = 4 and 5, S = 1, sigma in {0.05, 1, 20}; Nelder-Mead on softmax
coordinates from 12 random starts; prints the best minimiser (sorted U-weights) and f there, and the best values of f
on the two-level families (a, b, ..., b) and (a, a, b, ..., b).  Then a concavity test of the unnormalised
g(s) = (G_z v)_rho along random segments of the simplex (it fails: g is not concave)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize
from ex_dep import Gz

rng = np.random.default_rng(11)
S = 1.0


def margins(sig, sU):
    u = len(sU)
    G = Gz(sig, list(sU))
    v = np.array([4 * sig * (sig + S)] + [3 * (2 * sig + S)] * u)
    return G @ v, G, v


def f(sig, sU):
    Gv, G, v = margins(sig, sU)
    return Gv[0] / (G[0, 0] * v[0])


def from_logits(z):
    e = np.exp(z - z.max())
    return S * e / e.sum()


for u in (4, 5):
    for sig in (0.05, 1.0, 20.0):
        best = (np.inf, None)
        for it in range(12):
            z0 = rng.normal(0, 1.5, u)
            r = minimize(lambda z: f(sig, from_logits(z)), z0, method='Nelder-Mead',
                         options=dict(xatol=1e-7, fatol=1e-12, maxiter=4000))
            if r.fun < best[0]:
                best = (r.fun, np.sort(from_logits(r.x)))
        ts = np.linspace(0.001, 0.999, 300)
        one = min(f(sig, [t] + [(S - t) / (u - 1)] * (u - 1)) for t in ts)
        two = min(f(sig, [t / 2, t / 2] + [(S - t) / (u - 2)] * (u - 2)) for t in ts)
        eq = f(sig, [S / u] * u)
        print(f'u = {u}, sigma = {sig}: min f (numerical) = {best[0]:.6f} at sorted s_U = {np.round(best[1], 4)}; '
              f'f at equal weights {eq:.6f}; best on (a, b, ..., b) {one:.6f}; best on (a, a, b, ..., b) {two:.6f}',
              flush=True)

u = 4
for sig in (0.05, 1.0, 20.0):
    worst = np.inf
    for it in range(300):
        a = rng.dirichlet(np.ones(u) * 0.7) * S
        b = rng.dirichlet(np.ones(u) * 0.7) * S
        a = np.maximum(a, 1e-6); b = np.maximum(b, 1e-6)
        a *= S / a.sum(); b *= S / b.sum()
        g = lambda x: margins(sig, x)[0][0]
        gap = g((a + b) / 2) - (g(a) + g(b)) / 2
        worst = min(worst, gap / max(abs(g(a)), abs(g(b))))
    print(f'u = 4, sigma = {sig}: min over 300 random segments of [g(mid) - mean g(ends)]/max|g(ends)| = {worst:.3e} '
          f'({"concave on all samples" if worst >= -1e-12 else "not concave"})', flush=True)
print('PASS')
