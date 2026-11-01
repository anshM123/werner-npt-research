"""Candidate rules for the merged atom (E, F) of X = rho + G(C, D), tested on random NPT targets.
For a rule's (E, F): lam = argmax of the PPT margin mu(lam) = min(lambda_min(S), lambda_min(S^Gamma)),
S = X - lam G^(E, F); then S is fitted by products only (NNLS + column generation).  A rule works on X if mu >= 0 and
the product fit of S is exact (separable remainder).
rules: 'neg'  E F^* = mat(negative PT eigenvector of X), E = that matrix, F = I (canonical rule of NOTES Sec. 11);
       'vdm'  normal-form frame: (a (x) b) X (a (x) b)^* has maximally mixed marginals; W = polar part of the negative
              PT eigenvector there; E = a^-1, F = b^-1 W^*  (and the other polar orientation);
       'mid'  principal square roots E = C^{1/2}, F = D^{1/2} (on rho + G(I, W) this is the geodesic midpoint);
       'mu'   the maximiser of the PPT margin over all atoms (local search from many starts; empirically unique).
usage: python rules_test.py seed count"""
import sys, time, numpy as np, scipy.linalg as sla
from bm_core import G, RHO, ptmin, rand_c, D, N, pt, vecr, SepDict, nnls_fit, best_products, unvecr
from bm_fit import norm_t, candidates, vdm_frame, polar_u
from ppt_margin import maximize_mu, mu_hard

seed, count = int(sys.argv[1]), int(sys.argv[2])
g = np.random.default_rng(seed)


def best_lam(X, E, F):
    Gh = G(E, F); Gh = Gh / np.trace(Gh).real
    lams = np.linspace(0, 1.2, 1201)
    vals = np.array([min(np.linalg.eigvalsh(X - l * Gh)[0], np.linalg.eigvalsh(pt(X - l * Gh))[0]) for l in lams])
    i = int(np.argmax(vals)); return lams[i], vals[i], Gh


def sep_fit(S, rounds=40):
    """relative residual of the best product fit of S (NNLS + column generation)"""
    nS = np.linalg.norm(vecr(S))
    if nS == 0:
        return 0.0
    dic = SepDict(g, 1500); x = vecr(S)
    for rd in range(rounds):
        coef, _ = nnls_fit(x, dic.V); r = x - dic.V @ coef
        rr = np.linalg.norm(r) / nS
        if rr < 1e-13:
            break
        new = best_products(unvecr(r), g, nstart=24, keep=12)
        xs = [a for v, a, b in new if v > 0]; ys = [b for v, a, b in new if v > 0]
        if not xs:
            break
        dic.add_many(xs, ys)
    return rr


def rule_atoms(X, C, Dm):
    out = {}
    w, U = np.linalg.eigh(pt(X)); M = U[:, 0].reshape(D, D)
    out['neg'] = [(M, np.eye(D)), (np.eye(D), M.conj().T)]
    a, b, Y = vdm_frame(X)
    w, U = np.linalg.eigh(pt(Y)); Mn = U[:, 0].reshape(D, D); Wp = polar_u(Mn)
    out['vdm'] = [(np.linalg.inv(a), np.linalg.inv(b) @ Wp.conj().T), (np.linalg.inv(a) @ Wp, np.linalg.inv(b))]
    out['mid'] = [(sla.sqrtm(C), sla.sqrtm(Dm))]
    return out


t0 = time.time()
summary = {k: [] for k in ('neg', 'vdm', 'mid', 'mu')}
done = 0
while done < count:
    C = np.linalg.qr(rand_c(g, D, D))[0] + 0.35 * rand_c(g, D, D); Dm = np.linalg.qr(rand_c(g, D, D))[0] + 0.35 * rand_c(g, D, D)
    Gc = G(C, Dm); wgt = np.exp(g.uniform(np.log(0.3), np.log(3)))
    s = np.sqrt(wgt * np.trace(RHO).real / np.trace(Gc).real); C = C * s
    X = norm_t(RHO + G(C, Dm))
    if ptmin(X) > -0.005:
        continue
    done += 1
    line = f"#{done}: PT-min {ptmin(X):+.4f}"
    ra = rule_atoms(X, C, Dm)
    for k in ('neg', 'vdm', 'mid'):
        best = None
        for (E, F) in ra[k]:
            lam, mu, Gh = best_lam(X, E, F)
            if best is None or mu > best[1]:
                best = (lam, mu, Gh)
        lam, mu, Gh = best
        res = sep_fit(X - lam * Gh) if mu >= -1e-12 else np.nan
        summary[k].append((mu, res))
        line += f" | {k}: mu {mu:+.5f}, sep-fit {res:.1e}"
    m, E, F, lam = maximize_mu(X, candidates(X) + [(np.eye(D), np.eye(D)), (C, Dm)] + [(rand_c(g, D, D), rand_c(g, D, D)) for _ in range(3)], iters=300)
    Gh = G(E, F); Gh /= np.trace(Gh).real
    res = sep_fit(X - lam * Gh)
    summary['mu'].append((m, res))
    line += f" | mu-max: mu {m:+.5f}, sep-fit {res:.1e}  [{time.time() - t0:.0f}s]"
    print(line, flush=True)
for k, v in summary.items():
    v = np.array(v)
    ok = np.sum((v[:, 0] >= -1e-12) & (v[:, 1] < 1e-10))
    print(f"rule {k}: PPT-feasible on {int(np.sum(v[:, 0] >= -1e-12))}/{len(v)}, separable remainder (fit < 1e-10) on {ok}/{len(v)}; "
          f"worst mu {v[:, 0].min():+.5f}", flush=True)
