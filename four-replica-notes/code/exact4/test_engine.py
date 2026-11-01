"""Validation of engine4 against the certified exact data of ../relaxations/exact_coreN.py (n = 2, 3).
At a random rational point of the common variable space (built in exact_coreN's basis, converted to ours by the exact
projection), compare exactly: objective, all moment functionals, raw images of every identity map (all sorted target
patterns), and every raw covariance block.  Also compares the exact hats of a true point (our integer Fourier
transform vs exact_coreN.fourier_from_moments).
usage: python test_engine.py n  (n = 4: objective and the covariance blocks in TEST4_TARGETS only)"""
import sys, time, itertools, numpy as np, flint
from gmpy2 import mpq
from core4 import *
from engine4 import *
import exact_coreN as XN

n = int(sys.argv[1]); t0 = time.time()
A, B1, A3 = site_algebras()
vs = VarSpace(A, n)
rel = XN.ExactRelaxN(n, A, B1, A3, log=lambda s: None)
FULL = n <= 3
if FULL: rel.build(log=lambda s: None)
else:
    w01_ = PERM[(1, 0, 2, 3)] * 2 - np.eye(81, dtype=np.int64); w12_ = PERM[(0, 2, 1, 3)] * 2 - np.eye(81, dtype=np.int64)
    rel.c = (rel.functional('uuuu', [w01_] * n) + rel.functional('vvvv', [w01_] * n) + rel.functional('uuvv', [w12_] * n) * 2) * mpq(1, 2 ** n)
LIN = len(sys.argv) > 2 and sys.argv[2] == 'more'
TARGETS = [('uu', (0, 0, 1, 2)), ('uv', (0, 1, 2, 3))]
if LIN:
    TARGETS = [('uu', (0, 0, 1, 2)), ('uv', (0, 1, 2, 3)), ('uu', (1, 1, 1, 1)), ('vv', (3, 3, 3, 3)), ('uv', (2, 3, 3, 3)),
               ('uv', (1, 1, 2, 2)), ('uu', (0, 1, 3, 3)), ('vv', (1, 2, 3, 3)), ('uv', (0, 0, 0, 0)), ('uu', (2, 2, 3, 3))]
print(f"built both: ours N={vs.N}, exact_coreN N={rel.N}; {time.time() - t0:.1f}s", flush=True)
assert vs.N == rel.N
nb = len(A.blocks)
def m2f(q): return fq(int(q.numerator), int(q.denominator))
def obj2fm(H): return flint.fmpq_mat(H.shape[0], H.shape[1], [m2f(mpq(v)) for v in H.ravel()])

rng = np.random.default_rng(5)
xr = np.array([mpq(int(rng.integers(-5, 6)), int(rng.integers(1, 4))) for _ in range(rel.N)], dtype=object)
hats = {c: {} for c in FAMS}
for c in FAMS:
    for lam in sorted_patterns(nb, n):
        D = int(np.prod([A.blocks[l]['f'] for l in lam]))
        if (c, lam) in rel.var:
            vd = rel.var[(c, lam)]; H = XN.zeros_o((D, D))
            for k, i in enumerate(vd['idx']): H = H + vd['Xh'][k] * xr[i]
        else:
            H = XN.zeros_o((D, D))
        hats[c][lam] = obj2fm(H)
xo, ok = vs.x_from_hats(hats)
print(f"conversion to our basis exact: {ok}", flush=True)
xo_m = np.array([mpq(int(v.p), int(v.q)) for v in xo], dtype=object)

def images_of_point(Ti, dT, tg_nb, c, x, mus):
    """sum_i x_i * (exact image of parameter i) for every sorted target mu: dict mu -> object (Dm, Dm)"""
    out = {}
    for key in vs.order:
        if key[0] != c: continue
        g = vs.groups[key]
        Xp = vs.pair_float(key)
        level = tree(Xp, key[1], Ti, tg_nb)
        bound = g['Xmax'] * rowsum_bound(Ti) ** n * len(g['orbit'])
        for mu in mus:
            acc = assemble(level, g['orbit'], mu)
            if acc is None: continue
            check_exact(acc, bound)
            gs = tuple(int(round(np.sqrt(sz))) for sz in acc.shape[1:])
            Y = pair_to_mat(acc, gs).astype(np.int64).astype(object)
            sc = mpq(1, dT ** n * g['dX'])
            tot = sum((Y[k] * x[i] for k, i in enumerate(g['idx'])), XN.zeros_o(Y.shape[1:]))
            out[mu] = out[mu] + tot * sc if mu in out else tot * sc
    return out

# objective and moments
fails = 0
w01 = PERM[(1, 0, 2, 3)] * 2 - np.eye(81, dtype=np.int64); w12 = PERM[(0, 2, 1, 3)] * 2 - np.eye(81, dtype=np.int64)
def functional_value(c, w, x):
    Ti, d = functional_map_int(A, w)
    tot = mpq(0)
    for key in vs.order:
        if key[0] != c: continue
        g = vs.groups[key]
        level = tree(vs.pair_float(key), key[1], Ti, 1)
        acc = assemble(level, g['orbit'], (0,) * n)
        if acc is None: continue
        check_exact(acc, g['Xmax'] * rowsum_bound(Ti) ** n * len(g['orbit']))
        vals = acc.reshape(g['K']).astype(np.int64)
        for k, i in enumerate(g['idx']): tot += int(vals[k]) * x[i] * mpq(g['mlam'], d ** n * g['dX'])
    return tot
obj_o = (functional_value('uuuu', w01, xo_m) + functional_value('vvvv', w01, xo_m) + 2 * functional_value('uuvv', w12, xo_m)) * mpq(1, 2 ** n)
obj_r = rel.c.dot(xr)
print(f"objective: ours {float(obj_o):.12f}, exact_coreN {float(obj_r):.12f}, equal {obj_o == obj_r}", flush=True)
fails += obj_o != obj_r
k = 0; bad = 0
if FULL: Mm, em = rel.mom
elif LIN: Mm = [rel.functional(c, [PERM[tau]] * n) for c in FAMS for tau in S4]
for c in (FAMS if (FULL or LIN) else []):
    for tau in S4:
        if functional_value(c, PERM[tau], xo_m) != Mm[k].dot(xr): bad += 1
        k += 1
print(f"moment functionals: {k} compared, {bad} mismatches; {time.time() - t0:.1f}s", flush=True)
fails += bad
# identities: raw images
bad = 0; cnt = 0
for name, terms in (identity_defs() if (FULL or LIN) else []):
    for mu in sorted_patterns(len(A3.blocks), n):
        ours = None
        for coef, c, phi in terms:
            Ti, dT = site_map_int(A, A3, phi)
            im = images_of_point(Ti, dT, len(A3.blocks), c, xo_m, [mu])
            if mu in im: ours = im[mu] * coef if ours is None else ours + im[mu] * coef
        theirs = None
        for coef, c, phi in terms:
            smo = rel.sitemap_obj(SiteMap(A, A3, phi))
            acc = rel.target_image(A3, mu, [(mpq(coef), c, smo)])
            for i, Y in acc.items():
                theirs = Y * xr[i] if theirs is None else theirs + Y * xr[i]
        cnt += 1
        if (ours is None) != (theirs is None) or (ours is not None and not all(a == b for a, b in zip(ours.ravel(), theirs.ravel()))):
            bad += 1
print(f"identity images: {cnt} (map, target) pairs compared, {bad} mismatches; {time.time() - t0:.1f}s", flush=True)
fails += bad
# covariance raw blocks
bad = 0; cnt = 0
covmaps = {}
for (xy, r, s_, c, P) in cov_entries():
    Pm = PERM[P]
    covmaps[(xy, r, s_)] = (c,) + site_map_int(A, B1, lambda X, Pm=Pm: realign1(Pm @ X @ Pm.T))
for xy in ['uu', 'uv', 'vv']:
    mlist = sorted_patterns(len(B1.blocks), n) if FULL else [mu for (x_, mu) in TARGETS if x_ == xy]
    if not mlist: continue
    ims = {(r, s_): images_of_point(covmaps[(xy, r, s_)][1], covmaps[(xy, r, s_)][2], len(B1.blocks), covmaps[(xy, r, s_)][0], xo_m,
                                    mlist) for r in range(2) for s_ in range(2)}
    for mu in mlist:
        g = [B1.blocks[q]['f'] for q in mu]; Dm = int(np.prod(g))
        Gt = XN.okron([XN.fm2np(B1.blocks[q]['G']) for q in mu])
        Bk = XN.zeros_o((2 * Dm, 2 * Dm)); nz = False
        for (r, s_), im in ims.items():
            if mu in im: Bk[r * Dm:(r + 1) * Dm, s_ * Dm:(s_ + 1) * Dm] = Gt.dot(im[mu]); nz = True
        Bk = (Bk + Bk.T) * mpq(1, 2)
        if FULL:
            blk = [b for b in rel.psd if b[0] == f"COV{xy}{mu}"]
            theirs = rel.block_at(blk[0], xr) if blk else None
        else:
            theirs = XN.zeros_o((2 * Dm, 2 * Dm)); anyd = False
            for r, lr in enumerate('uv'):
                for s_, ls in enumerate('uv'):
                    c_, sp_ = locate(xy + lr + ls); Pm_ = PERM[sp_]
                    smo = rel.sitemap_obj(SiteMap(A, B1, lambda X, P=Pm_: realign1(P @ X @ P.T)))
                    acc = rel.target_image(B1, mu, [(mpq(1), c_, smo)])
                    for i, Y in acc.items():
                        theirs[r * Dm:(r + 1) * Dm, s_ * Dm:(s_ + 1) * Dm] += Gt.dot(Y) * xr[i]; anyd = True
            theirs = (theirs + theirs.T) * mpq(1, 2)
            print(f"   COV{xy}{mu}: exact_coreN block computed; {time.time() - t0:.1f}s", flush=True)
        cnt += 1
        if theirs is None:
            if any(v != 0 for v in Bk.ravel()): bad += 1
        elif not all(a == b for a, b in zip(Bk.ravel(), theirs.ravel())): bad += 1
print(f"covariance raw blocks: {cnt} compared, {bad} mismatches; {time.time() - t0:.1f}s", flush=True)
fails += bad
print(f"n={n}: ENGINE {'AGREES' if fails == 0 else 'DISAGREES'} with exact_coreN ({fails} failures)", flush=True)
