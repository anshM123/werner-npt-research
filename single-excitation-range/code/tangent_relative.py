"""Evidence for the tangent form of Conjecture 12 (README, remark after Proposition 14).
For p = sum_t w_t |1_t> in H_1, the compression of Y_p = (M - R_p)^Gamma to (H_0 (x) H_k) + (H_1 (x) H_{k-1}) should be
>= 0 for every k; by Lemma 11' this covers every plane with range in span{|0...0>} + H_1 that contains p.
Usage: python tangent_relative.py n   (random, structured and Nelder-Mead weights; evidence, not proof)."""
import numpy as np, itertools, sys
from scipy.optimize import minimize
from relative_anchor import Mel, strings
def build(n, k):
    vac = (0,) * n; H1 = strings(n, 1)
    Hk = strings(n, k) if 0 <= k <= n else []; Hk1 = strings(n, k - 1) if 1 <= k <= n + 1 else []
    Sk = [(vac, y) for y in Hk] + [(a, y) for a in H1 for y in Hk1]
    allcols = list(itertools.product(range(3), repeat=n))
    U = sorted(set([(b, yp) for (a, yp) in Sk for (b, y) in Sk] + [(a, y) for (a, yp) in Sk for (b, y) in Sk]))
    ui = {kk: i for i, kk in enumerate(U)}
    Gt = np.array([[[Mel(r, y, b, x) for b in H1] for x in allcols] for (r, y) in U])
    Dt = np.array([[[[Mel(a, x, b, xp) for b in H1] for a in H1] for xp in allcols] for x in allcols])
    Mu = np.array([[Mel(r, y, rp, yp) for (rp, yp) in U] for (r, y) in U])
    P1 = np.array([[ui[(b, yp)] for (b, y) in Sk] for (a, yp) in Sk]); P2 = np.array([[ui[(a, y)] for (b, y) in Sk] for (a, yp) in Sk])
    return dict(H1=H1, Gt=Gt, Dt=Dt, Mu=Mu, P1=P1, P2=P2, n=n, dim=len(Sk))
def lmin(w, B):
    p = np.array([w[a.index(1)] if 1 in a else 0.0 for a in B['H1']])
    D = np.einsum('ijab,a,b->ij', B['Dt'], p, p); G = B['Gt'] @ p
    W = B['Mu'] - G @ np.linalg.solve(D, G.T)
    Y = W[B['P1'], B['P2']]; return np.linalg.eigvalsh((Y + Y.T) / 2)[0]
if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3; rng = np.random.default_rng(1); worst = 9
    for k in range(0, n + 2):
        B = build(n, k)
        vals = [lmin(rng.random(n) + 0.05, B) for _ in range(4)] + [lmin(np.array([1.0] + [0.05] * (n - 1)), B), lmin(np.ones(n), B)]
        r = minimize(lambda x: lmin(np.exp(x), B), rng.normal(size=n), method='Nelder-Mead', options=dict(maxiter=100 * n))
        worst = min(worst, min(vals), r.fun)
        print(f"n={n} k={k} (dim {B['dim']}): random/structured min {min(vals):+.3e}, Nelder-Mead {r.fun:+.3e}", flush=True)
    print("PASS" if worst > -1e-12 else "FAIL")
