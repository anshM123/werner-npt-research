"""Three-copy thickness certificate, step 1: exact data and the floating-point SDP.

Exact data: per-copy algebras (exact_alg.py), embeddings, target, orbits, forced kernels and complement bases R_beta
(integer), relation span.  Then a floating-point SDP in the parametrization Gamma_beta = R_beta Y_beta R_beta^T
(Gamma = Gram blocks V^T X V in the integer spin bases), maximising the least eigenvalue of all Y_beta at
kappa = 98/27.  Saves exact data (pickle) and the float solution (npz).

Certificate shape (copy 0 is the special copy; everything is averaged over copy permutations):
   target = Z0 (x) I_w  +  Z1  +  G(L1)  (+ operator relations),  modulo relabeling null forms,
   Z0: 3 strands (x, y, y') per copy;  Z1: 4 strands (x, y, y', w) at copy 0, 3 strands (x, y, y') at copies 1, 2;
   L1: 3 strands (idx = y, x, y') at copy 0, 2 strands (x, y') at copies 1, 2 (identity on w everywhere else)."""
import itertools
import pickle
import sys
import time
from fractions import Fraction
import numpy as np
import flint
import cvxpy as cp
from exact_alg import build, canon, nullspace_int, fz

T0 = time.time()


def log(*a):
    print(*a, f'[{time.time() - T0:.1f}s]', flush=True)


n = 3
A2, A3, A4 = build(2), build(3), build(4)
N2, N3, N4 = len(A2['diags']), len(A3['diags']), len(A4['diags'])
idx4 = {m: i for i, m in enumerate(A4['diags'])}
ta = Fraction(3, 4)
tt = Fraction(1, 64)
kap = Fraction(98, 27)


def dg4(pairs):
    return idx4[canon(pairs)]


def embed(m_small, sigma):
    k = len(sigma)
    pts = {}
    for i in range(k):
        pts[i] = sigma[i]
        pts[k + i] = 4 + sigma[i]
    m = [(pts[p], pts[q]) for p, q in m_small] + [(r, 4 + r) for r in range(4) if r not in sigma]
    return idx4[canon(m)]


def trace_first(m3):
    pairs = [tuple(p) for p in m3]
    if (0, 3) in pairs:
        return 3, [p for p in pairs if p != (0, 3)]
    pt = [q for p in pairs for q in p if 0 in p and q != 0][0]
    pb = [q for p in pairs for q in p if 3 in p and q != 3][0]
    rest = [p for p in pairs if 0 not in p and 3 not in p]
    return 1, rest + [(pt, pb)]


# integer embedding matrices (N4 x N3) and the 2-strand embedding (N4 x N2)
E_xyy = np.zeros((N4, N3), dtype=np.int64)
E_s1 = np.zeros((N4, N3), dtype=np.int64)
E_s2 = np.zeros((N4, N3), dtype=np.int64)
for ai, m3 in enumerate(A3['diags']):
    E_xyy[embed(m3, [0, 1, 2]), ai] += 1
    E_s1[embed(m3, [1, 0, 2]), ai] += 1
    fac, rest = trace_first(m3)
    mp = {1: 0, 2: 2, 4: 4, 5: 6}
    E_s2[idx4[canon([(mp[p], mp[q]) for p, q in rest] + [(1, 3), (5, 7)])], ai] += fac
E_2 = np.zeros((N4, N2), dtype=np.int64)
for ai, m2 in enumerate(A2['diags']):
    E_2[embed(m2, [0, 2]), ai] += 1
assert all(E_xyy[:, a].sum() == 1 for a in range(N3))

I4 = dg4([(0, 4), (1, 5), (2, 6), (3, 7)])
F01 = dg4([(0, 5), (1, 4), (2, 6), (3, 7)])
O01 = dg4([(0, 1), (4, 5), (2, 6), (3, 7)])
O13 = dg4([(1, 3), (5, 7), (0, 4), (2, 6)])
cvec = {I4: ta, F01: Fraction(-1, 4), O01: Fraction(-1, 4)}
ivec = {I4: Fraction(1)}
omyw = {O13: Fraction(1)}

# ---------------------------------------------------------------- relabeling group and orbits
def perm_on_diagrams(pp):
    return np.array([idx4[canon([(pp[p], pp[q]) for p, q in m])] for m in A4['diags']])


glob = []
for ex in (0, 1):
    for py in itertools.permutations([1, 2, 5, 6]):
        pp = list(range(8))
        if ex:
            pp[0], pp[4] = 4, 0
        for src, dst in zip([1, 2, 5, 6], py):
            pp[src] = dst
        glob.append(perm_on_diagrams(pp))
wsw = list(range(8))
wsw[3], wsw[7] = 7, 3
wperm = perm_on_diagrams(wsw)
w_id = sorted(set(int(E_xyy[:, a].argmax()) for a in range(N3)))      # the 15 diagrams that are identity on w
w_id_set = set(w_id)
rows = [tup for tup in itertools.product(w_id, repeat=n)]
for s in range(n):
    for tup in itertools.product(*[range(N4) if r == s else w_id for r in range(n)]):
        if tup[s] not in w_id_set:
            rows.append(tup)
tuples = np.array(rows, dtype=np.int64)
NT = len(tuples)
weights = N4 ** np.arange(n - 1, -1, -1)
best = None
for g in glob:
    tg = g[tuples]
    for ew in itertools.product((0, 1), repeat=n):
        tw = tg.copy()
        for s in range(n):
            if ew[s]:
                tw[:, s] = wperm[tw[:, s]]
        for perm in itertools.permutations(range(n)):
            key = tw[:, list(perm)] @ weights
            best = key if best is None else np.minimum(best, key)
reps, orb = np.unique(best, return_inverse=True)
NO = len(reps)
orbit_of = {int(f): int(o) for f, o in zip(tuples @ weights, orb)}
log('tuples', NT, 'orbits', NO)


def orbit_sums_exact(coeffs):
    """coeffs: dict tuple -> Fraction (4-strand diagram tuples) -> list of NO Fractions"""
    out = [Fraction(0)] * NO
    for tup, c in coeffs.items():
        if c == 0:
            continue
        o = orbit_of[int(np.dot(tup, weights))]
        out[o] += c
    return out


def kron_dicts(ds):
    out = {}
    for combo in itertools.product(*[list(d.items()) for d in ds]):
        tup = tuple(k for k, _ in combo)
        v = Fraction(1)
        for _, c in combo:
            v *= c
        out[tup] = out.get(tup, Fraction(0)) + v
    return out


lhs_c = orbit_sums_exact(kron_dicts([cvec] * n))
lhs_i = orbit_sums_exact(kron_dicts([ivec] * n))
lhs_k = [Fraction(0)] * NO
for s in range(n):
    part = orbit_sums_exact(kron_dicts([omyw if r == s else ivec for r in range(n)]))
    lhs_k = [x + y for x, y in zip(lhs_k, part)]
# lhs_k is the orbit sum of +Omega_{yw}; the target is f - t - t kap sum <w|rho|w>
target = [c - tt * i - tt * kap * k for c, i, k in zip(lhs_c, lhs_i, lhs_k)]
log('exact target built')

# ---------------------------------------------------------------- forced kernels (exact, in V_j coordinates)
def perm_matrix(k, perm):
    """permutation of the k strands of (R^3)^{(x)k}: (P v)[i_perm] = v[i]"""
    dim = 3 ** k
    P = np.zeros((dim, dim), dtype=np.int64)
    for ii in itertools.product(range(3), repeat=k):
        jj = tuple(ii[perm[s]] for s in range(k))
        r = int(np.dot(ii, 3 ** np.arange(k - 1, -1, -1)))
        c = int(np.dot(jj, 3 ** np.arange(k - 1, -1, -1)))
        P[r, c] = 1
    return P


def sym_constraint(k, strands):
    """integer matrix whose kernel is the set of vectors symmetric in the given strands"""
    tot = np.zeros((3 ** k, 3 ** k), dtype=np.int64)
    cnt = 0
    for p in itertools.permutations(strands):
        perm = list(range(k))
        for a_, b_ in zip(strands, p):
            perm[a_] = b_
        tot += perm_matrix(k, perm)
        cnt += 1
    return cnt * np.eye(3 ** k, dtype=np.int64) - tot


def kernel_coords(V, constraint):
    """integer coordinates u (columns) with constraint @ V u = 0"""
    return nullspace_int((constraint @ V).astype(object))


K3 = {j: kernel_coords(A3['V'][j], sym_constraint(3, [0, 1, 2])) for j in A3['spins']}
# 4 strands: Sym^3 on strands 0,1,2 tensor R^3, orthogonal to T_A(e,e,e,w) = (e^T A e)(e.w) for symmetric A
TA = []
for (i0, j0) in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]:
    A = np.zeros((3, 3), dtype=np.int64)
    A[i0, j0] = 1
    A[j0, i0] = 1
    T = np.zeros((3, 3, 3, 3), dtype=np.int64)
    for p in itertools.permutations([0, 1, 2]):
        # symmetrise A_{ab} delta_{cl} over the first three slots
        for a_, b_, c_, l_ in itertools.product(range(3), repeat=4):
            idx = [a_, b_, c_]
            T[idx[p[0]], idx[p[1]], idx[p[2]], l_] += A[a_, b_] * (1 if c_ == l_ else 0)
    TA.append(T.reshape(-1))
TA = np.array(TA, dtype=np.int64)
cons4 = np.vstack([sym_constraint(4, [0, 1, 2]), TA])
K4 = {j: kernel_coords(A4['V'][j], cons4) for j in A4['spins']}
K2 = {0: np.array([[1]], dtype=object), 1: np.zeros((1, 0), dtype=object), 2: np.array([[1]], dtype=object)}
log('forced kernel dims 3-strand', {j: K3[j].shape[1] for j in K3}, '4-strand', {j: K4[j].shape[1] for j in K4})

# numerical cross-check of K4 against the sampled span of e e e w with w perp e
rng = np.random.default_rng(1)
vs = []
for _ in range(80):
    e = rng.standard_normal(3)
    w = rng.standard_normal(3)
    w -= (w @ e) / (e @ e) * e
    vs.append(np.kron(np.kron(np.kron(e, e), e), w))
Ks = np.linalg.svd(np.array(vs).T, full_matrices=False)[0][:, :24]
for j in A4['spins']:
    Vj = A4['V'][j].astype(float)
    if K4[j].shape[1]:
        vec = Vj @ K4[j].astype(float)
        assert np.abs(vec - Ks @ (Ks.T @ vec)).max() < 1e-8 * np.abs(vec).max()
    # dimension check: intersection of span(V_j) with the sampled span
    Q = np.linalg.qr(Vj)[0]
    M = Q.T @ Ks
    sv = np.linalg.svd(M, compute_uv=False)
    assert int((sv > 1 - 1e-8).sum()) == K4[j].shape[1], (j, sv)
log('K4 agrees with the sampled span')

alg = {'A2': A2, 'A3': A3, 'A4': A4}
kers = {'A2': K2, 'A3': K3, 'A4': K4}
terms = {
    'Z0': ['A3', 'A3', 'A3'],
    'Z1': ['A4', 'A3', 'A3'],
    'L1': ['A3', 'A2', 'A2'],
}


def block_list(names):
    out = []
    for beta in itertools.product(*[alg[nm]['spins'] for nm in names]):
        ms = [alg[nm]['mj'][j] for nm, j in zip(names, beta)]
        b = int(np.prod(ms))
        Kb = np.ones((1, 1), dtype=object)
        for nm, j in zip(names, beta):
            Kb = np.kron(Kb, kers[nm][j]) if kers[nm][j].shape[1] else np.zeros((Kb.shape[0] * kers[nm][j].shape[0], 0), dtype=object)
        if Kb.shape[1] == 0:
            R = np.eye(b, dtype=np.int64).astype(object)
        else:
            R = nullspace_int(Kb.T)
        out.append(dict(beta=beta, ms=ms, b=b, K=Kb, R=R))
    return out


blocks = {nm: block_list(names) for nm, names in terms.items()}
log('blocks', {nm: (len(bl), sum(B['R'].shape[1] for B in bl)) for nm, bl in blocks.items()})


# ---------------------------------------------------------------- float maps: reduced block entries -> orbit sums
def coord_index(nm_list, beta, rr, cc):
    f = 0
    for nm, j, r, c in zip(nm_list, beta, rr, cc):
        al = alg[nm]
        f = f * len(al['coords']) + al['coords'].index((j, r, c))
    return f


Psi_f = {nm: np.array([[float(x) for x in row] for row in alg[nm]['Psi']]) for nm in alg}
copy_embed = {
    'Z0': [E_xyy, E_xyy, E_xyy],
    'Z1': [np.eye(N4, dtype=np.int64), E_xyy, E_xyy],
    'L1': [E_s1 - E_s2, E_2, E_2],
}
pos = {int(f): i for i, f in enumerate(tuples @ weights)}


def orbit_sum_float(factors):
    cs = [F.shape[1] for F in factors]
    out = np.zeros((NO, int(np.prod(cs))))
    supp = [np.nonzero(np.abs(F).sum(axis=1) > 1e-14)[0] for F in factors]
    sel = []
    for tup in itertools.product(*supp):
        sel.append(pos[int(np.dot(tup, weights))])
    sel = np.array(sel)
    tt_ = tuples[sel]
    oo = orb[sel]
    order = np.argsort(oo)
    tt_, oo = tt_[order], oo[order]
    bounds = np.searchsorted(oo, np.arange(NO + 1))
    for o in range(NO):
        lo, hi = bounds[o], bounds[o + 1]
        if lo < hi:
            out[o] = np.einsum('ka,kb,kc->abc', *[factors[s][tt_[lo:hi, s]] for s in range(n)], optimize=True).reshape(-1)
    return out


maps = {}
for nm, names in terms.items():
    factors = [copy_embed[nm][s].astype(float) @ Psi_f[names[s]] for s in range(n)]
    maps[nm] = orbit_sum_float(factors)
    log('map', nm, maps[nm].shape)

# relations: relation (copy 0, 4 strands) x identity-on-w 3-strand diagrams at copies 1, 2
rel4 = np.array(A4['relations'], dtype=np.int64).T                     # N4 x 14
SRf = orbit_sum_float([rel4.astype(float), E_xyy.astype(float), E_xyy.astype(float)])
uR, sR, _ = np.linalg.svd(SRf, full_matrices=False)
rr_ = int((sR > 1e-9 * sR[0]).sum())
# exact independent subset of relation columns
import scipy.linalg as sla
_, _, pv = sla.qr(SRf, pivoting=True)
rel_cols = sorted(pv[:rr_].tolist())
log('relation span', rr_, 'columns', rel_cols[:5], '...')

# ---------------------------------------------------------------- float SDP
target_f = np.array([float(x) for x in target])
smin = cp.Variable()
cons = []
expr = 0
Yvars = []
for nm in terms:
    M = maps[nm]
    for B in blocks[nm]:
        R = B['R'].astype(float)
        br = R.shape[1]
        if br == 0:
            continue
        b = B['b']
        rows_ = list(itertools.product(*[range(x) for x in B['ms']]))
        fl = np.array([[coord_index(terms[nm], B['beta'], r1, c1) for c1 in rows_] for r1 in rows_])
        Mb = M[:, fl.reshape(-1)].reshape(NO, b, b)
        MB = np.einsum('oij,ia,jb->oab', Mb, R, R).reshape(NO, -1)
        # normalise the variable scale: Y' = D Y D with D = 1/|R col| keeps conditioning sane
        Y = cp.Variable((br, br), symmetric=True)
        cons.append(Y >> smin * np.eye(br))
        Yvars.append((nm, B['beta'], Y, MB))
        expr = expr + MB @ cp.reshape(Y, (br * br,), order='C')
rho = cp.Variable(len(rel_cols))
expr = expr + SRf[:, rel_cols] @ rho
cons.append(64 * expr == 64 * target_f)
cons.append(smin <= 1)
prob = cp.Problem(cp.Maximize(smin), cons)
log('solving float SDP', len(Yvars), 'blocks')
prob.solve(solver=cp.CLARABEL, verbose=False)
log('status', prob.status, 'least eigenvalue', smin.value, 'residual', np.abs(expr.value - target_f).max())
sol = {f'{nm}_{"".join(map(str, beta))}': (Y.value + Y.value.T) / 2 for nm, beta, Y, _ in Yvars}
np.savez('tc3_float.npz', rho=rho.value, **sol)
pickle.dump(dict(target=target, NO=NO, tuples=tuples, orb=orb, weights=weights, blocks=blocks, terms=terms,
                 copy_embed=copy_embed, rel_cols=rel_cols, kers=kers), open('tc3_setup.pkl', 'wb'))
log('saved')
