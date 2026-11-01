"""Random and structured scan of BM targets with the one-atom fit (fit_robust) and the PPT margin mu*.
families: 'nu' near-unitary C, D (unitary + 0.3 Gaussian), weight w in [0.2, 5];  'uni' rho + G(I, W), W Haar unitary;
'ill' C or D with condition number 1e-1..1e-3;  'comm' C = D diagonal positive;  'one' C = I, D Gaussian;
'gauss' C, D Gaussian.
usage: python scan_bm.py family seed count"""
import sys, time, numpy as np
from bm_core import G, RHO, ptmin, rand_c, D, N
from bm_fit import fit_robust, norm_t, candidates
from ppt_margin import maximize_mu

fam, seed, count = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
g = np.random.default_rng(seed)


def haar():
    Q, R = np.linalg.qr(rand_c(g, D, D)); return Q * (np.diag(R) / abs(np.diag(R)))


def make():
    if fam == 'nu':
        C = haar() + 0.3 * rand_c(g, D, D); Dm = haar() + 0.3 * rand_c(g, D, D)
    elif fam == 'uni':
        C = np.eye(D); Dm = haar()
    elif fam == 'ill':
        C = rand_c(g, D, D); Dm = rand_c(g, D, D)
        for M in ((C,) if g.random() < 0.5 else (C, Dm)):
            U, s, Vh = np.linalg.svd(M); s[-1] = s[0] * 10 ** g.uniform(-3, -1); M[:] = U @ np.diag(s) @ Vh
    elif fam == 'comm':
        C = np.diag(np.exp(g.uniform(-1.5, 1.5, D))).astype(complex); Dm = C.copy()
    elif fam == 'one':
        C = np.eye(D); Dm = rand_c(g, D, D)
    else:
        C, Dm = rand_c(g, D, D), rand_c(g, D, D)
    Gc = G(C, Dm); w = np.exp(g.uniform(np.log(0.2), np.log(5)))
    X = norm_t(RHO / np.trace(RHO).real + w * Gc / np.trace(Gc).real)
    return X, C, Dm, w


worst_fit = 0; worst_mu = np.inf; t0 = time.time()
for k in range(count):
    X, C, Dm, w = make()
    out = fit_robust(X, P=60, restarts=2, rounds=8, seed=int(g.integers(1 << 20)), extra=[(C, Dm), (np.eye(D), np.eye(D))])
    m, E, Fm, lam = maximize_mu(X, [(out['E'], out['F'])] + candidates(X) + [(C, Dm), (np.eye(D), np.eye(D))], iters=200)
    worst_fit = max(worst_fit, out['vp']); worst_mu = min(worst_mu, m)
    print(f"{fam} #{k}: w {w:.2f}, PT-min {ptmin(X):+.4f}, lambda_min {np.linalg.eigvalsh(X)[0]:.4f}; one-atom fit {out['vp']:.2e} "
          f"(atom wt {out['weight']:.3f}); mu* >= {m:+.5f} (lam {lam:.3f})  [{time.time() - t0:.0f}s]", flush=True)
print(f"{fam} seed {seed}: {count} targets; worst one-atom residual {worst_fit:.2e}; smallest mu* lower bound {worst_mu:+.5f}", flush=True)
