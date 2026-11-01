"""normal-form ('vdm') rule on the saved F3 targets of compression-notes (3 x 3 filterings of rho^{(x)M}) and on the
commuting boundary target; prints the rule's PPT margin, the best margin found over all atoms, and the product fit of the
rule's remainder.  usage: python rule_f3.py"""
import sys, glob, itertools, numpy as np
from bm_core import G, RHO, ptmin, rand_c, D, N, pt, vecr, SepDict, nnls_fit, best_products, unvecr
from bm_fit import norm_t, candidates, vdm_frame, polar_u
from ppt_margin import maximize_mu
g = np.random.default_rng(3)
CODE = "../"
def filt(A, B, m):
    subsets = [S for r in range(m + 1) for S in itertools.combinations(range(m), r)]
    At, Bt = A.reshape([D] + [3] * m), B.reshape([D] + [3] * m); Y = 0
    for S in subsets:
        Sc = [s for s in range(m) if s not in S]; pa = [0] + [1 + s for s in S] + [1 + s for s in Sc]
        Ap = At.transpose(pa).reshape(D, 3 ** len(S), 3 ** len(Sc)); Bp = Bt.transpose(pa).reshape(D, 3 ** len(S), 3 ** len(Sc))
        M1 = np.einsum('axc,eyc->aexy', Ap, Ap.conj()); M2 = np.einsum('byc,fxc->bfyx', Bp, Bp.conj())
        Y = Y + (-0.5) ** len(S) * np.einsum('aexy,bfyx->abef', M1, M2).reshape(N, N)
    return Y
def best_lam(X, Gh):
    lams = np.linspace(0, 1.2, 2401)
    vals = [min(np.linalg.eigvalsh(X - l * Gh)[0], np.linalg.eigvalsh(pt(X - l * Gh))[0]) for l in lams]
    i = int(np.argmax(vals)); return lams[i], vals[i]
def sep_fit(S, rounds=60):
    nS = np.linalg.norm(vecr(S)); dic = SepDict(g, 1500); x = vecr(S); rr = 1
    for rd in range(rounds):
        coef, _ = nnls_fit(x, dic.V); r = x - dic.V @ coef; rr = np.linalg.norm(r) / nS
        if rr < 1e-13: break
        new = best_products(unvecr(r), g, nstart=24, keep=12)
        xs = [p for v, p, q in new if v > 0]; ys = [q for v, p, q in new if v > 0]
        if not xs: break
        dic.add_many(xs, ys)
    return rr
def vdm_atom(X):
    a, b, Y = vdm_frame(X)
    w, U = np.linalg.eigh(pt(Y)); Wp = polar_u(U[:, 0].reshape(D, D))
    Gh = G(np.linalg.inv(a) @ Wp, np.linalg.inv(b)); return Gh / np.trace(Gh).real, (np.linalg.inv(a) @ Wp, np.linalg.inv(b))
targets = []
for f in sorted(glob.glob(CODE + "adv_f3_filters_M*_s*.npz")):
    Fl = np.load(f); M = int(f.split('_M')[1].split('_')[0])
    targets.append((f.split('/')[-1], norm_t(filt(Fl['A'], Fl['B'], M))))
Z = np.load("comm_s21_0.npz"); targets.append(("commuting boundary target (scan comm #0)", norm_t(Z['X'])))
for tag, X in targets:
    Gh, (E, F) = vdm_atom(X); lam, mu = best_lam(X, Gh)
    ms, _, _, _ = maximize_mu(X, [(E, F)] + candidates(X) + [(rand_c(g, D, D), rand_c(g, D, D)) for _ in range(4)], iters=300)
    sf = sep_fit(X - lam * Gh) if mu > -1e-9 else np.nan
    print(f"{tag}: PT-min {ptmin(X):+.4f}; rule margin {mu:+.6f} (lam {lam:.3f}); best margin over atoms >= {ms:+.6f}; rule remainder product-fit {sf:.1e}", flush=True)
