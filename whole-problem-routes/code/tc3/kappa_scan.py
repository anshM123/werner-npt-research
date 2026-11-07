"""Best bonus constant of the certificate family of Theorem 7.29 (floating point).

Same certificate shape as the exact proof (strands x, y, y', w per copy; w non-trivial at one copy per term; least-
eigenvector constraint with a PSD multiplier; identity on orbit sums of relabelings, with operator relations), with
the forced kernels removed, maximising the constant kappa in f >= 4^-n (1 + kappa sum_s tau_s).
usage: python kappa_scan.py n restricted none sn CLARABEL        (n = 2 gives 5.0000, n = 3 gives 4.1253)"""
import itertools
import sys
import time
import numpy as np
import scipy.sparse as sps
import cvxpy as cp

T0 = time.time()


def log(*a):
    print(*a, f'[{time.time() - T0:.1f}s]', flush=True)


n = int(sys.argv[1])
mode = sys.argv[2]
USE_SN = (len(sys.argv) <= 4 or sys.argv[4] != "nosn")
kap_fixed = float(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != "none" else None
a = 0.75
t = 4.0 ** -n
D = 3


def matchings(points):
    if not points:
        yield []
        return
    p = points[0]
    for i in range(1, len(points)):
        rest = points[1:i] + points[i + 1:]
        for m in matchings(rest):
            yield [(p, points[i])] + m


def canon(m):
    return tuple(sorted(tuple(sorted(p)) for p in m))


def algebra(K):
    ident = canon([(s, s + K) for s in range(K)])
    diag = [ident] + [m for m in sorted(set(canon(m) for m in matchings(list(range(2 * K))))) if m != ident]
    idx = {m: i for i, m in enumerate(diag)}
    mats = []
    for m in diag:
        M = np.zeros((D ** K, D ** K))
        for ii in itertools.product(range(D), repeat=2 * K):
            if all(ii[p] == ii[q] for p, q in m):
                r = c = 0
                for s in range(K):
                    r = r * D + ii[s]
                    c = c * D + ii[K + s]
                M[r, c] = 1
        mats.append(M)
    eps = np.zeros((3, 3, 3))
    for (i, j, k), s in [((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)]:
        eps[i, j, k] = s
    L = [-eps[k] for k in range(3)]
    tot = []
    for Lk in L:
        out = np.zeros((D ** K, D ** K))
        for s in range(K):
            f = [np.eye(3)] * K
            f[s] = Lk
            mm = f[0]
            for x in f[1:]:
                mm = np.kron(mm, x)
            out += mm
        tot.append(out)
    _, sv, vt = np.linalg.svd(tot[2])
    N0 = vt[sv < 1e-9].T
    Cas = sum(T_ @ T_ for T_ in tot)
    kk = N0.T @ Cas @ N0
    ev, evec = np.linalg.eigh((kk + kk.T) / 2)
    U = {}
    for j in range(K + 1):
        sel = np.abs(ev + j * (j + 1)) < 1e-6
        if sel.any():
            U[j] = N0 @ evec[:, sel]
    spins = sorted(U)
    mj = {j: U[j].shape[1] for j in spins}
    coords = [(j, r, c) for j in spins for r in range(mj[j]) for c in range(mj[j])]
    Phi = np.zeros((len(coords), len(diag)))
    for ai, M in enumerate(mats):
        for i, (j, r, c) in enumerate(coords):
            Phi[i, ai] = U[j][:, r] @ M @ U[j][:, c]
    Phinv = np.linalg.pinv(Phi)
    assert np.allclose(Phi @ Phinv, np.eye(len(coords)))
    u_, s_, vt_ = np.linalg.svd(Phi)
    Ker = vt_[int((s_ > 1e-9).sum()):].T
    return dict(K=K, diag=diag, idx=idx, mats=mats, spins=spins, mj=mj, coords=coords,
                cidx={c: i for i, c in enumerate(coords)}, Phi=Phi, Phinv=Phinv, Ker=Ker)


A4 = algebra(4)
A3 = algebra(3)
N4, N3 = len(A4['diag']), len(A3['diag'])
log('4 strands:', N4, 'diagrams', len(A4['coords']), 'block coords', A4['Ker'].shape[1], 'relations')
for i in range(A4['Ker'].shape[1]):
    assert np.abs(sum(A4['Ker'][k, i] * A4['mats'][k] for k in range(N4))).max() < 1e-9


def dg4(pairs):
    return A4['idx'][canon(pairs)]


def v4(d):
    v = np.zeros(N4)
    for k, c in d.items():
        v[k] += c
    return v


def embed(m_small, sigma):
    """matching on len(sigma) strands placed on the 4-strand strands sigma; identity elsewhere"""
    k = len(sigma)
    pts = {}
    for i in range(k):
        pts[i] = sigma[i]
        pts[k + i] = 4 + sigma[i]
    m = [(pts[p], pts[q]) for p, q in m_small] + [(r, 4 + r) for r in range(4) if r not in sigma]
    return A4['idx'][canon(m)]


def trace_first(m3):
    """partial trace of a 3-strand matching over strand 0 -> (factor, pairs on points {1,2,4,5})"""
    top, bot = 0, 3
    pairs = [tuple(p) for p in m3]
    if (0, 3) in pairs:
        return 3, [p for p in pairs if p != (0, 3)]
    pt = [q for p in pairs for q in p if top in p and q != top][0]
    pb = [q for p in pairs for q in p if bot in p and q != bot][0]
    rest = [p for p in pairs if top not in p and bot not in p]
    return 1, rest + [(pt, pb)]


I4 = dg4([(0, 4), (1, 5), (2, 6), (3, 7)])
cvec = v4({I4: a, dg4([(0, 5), (1, 4), (2, 6), (3, 7)]): -0.25, dg4([(0, 1), (4, 5), (2, 6), (3, 7)]): -0.25})
ivec = v4({I4: 1})
omyw = v4({dg4([(1, 3), (5, 7), (0, 4), (2, 6)]): 1})

# embedding matrices (N4 x N3): 3-strand algebra on (x, y, y') with identity on w   [for Z0 and Z_s others]
E_xyy = np.zeros((N4, N3))
# multiplier at the special copy: 3-strand (idx, x, y') -> (1, 0, 2), identity on w; minus trace term
E_s1 = np.zeros((N4, N3))
E_s2 = np.zeros((N4, N3))
# multiplier at other copies, full mode: 3-strand (x, y', w) -> (0, 2, 3)
E_xyw = np.zeros((N4, N3))
for ai, m3 in enumerate(A3['diag']):
    E_xyy[embed(m3, [0, 1, 2]), ai] += 1
    E_s1[embed(m3, [1, 0, 2]), ai] += 1
    E_xyw[embed(m3, [0, 2, 3]), ai] += 1
    fac, rest = trace_first(m3)
    mp = {1: 0, 2: 2, 4: 4, 5: 6}
    pairs = [(mp[p], mp[q]) for p, q in rest] + [(1, 3), (5, 7)]
    E_s2[A4['idx'][canon(pairs)], ai] += fac
# 2-strand projectors on (x, y') with identity on y, w   [restricted-mode multiplier at other copies]
F02 = dg4([(0, 6), (2, 4), (1, 5), (3, 7)])
O02 = dg4([(0, 2), (4, 6), (1, 5), (3, 7)])
P2 = [v4({O02: 1 / 3}), v4({I4: 0.5, F02: -0.5}), v4({I4: 0.5, F02: 0.5, O02: -1 / 3})]

# ---------------------------------------------------------------- relabeling group on 8 points
def perm_on_diagrams(pp):
    return np.array([A4['idx'][canon([(pp[p], pp[q]) for p, q in m])] for m in A4['diag']])


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
idperm = np.arange(N4)
w_id_set = set(np.where(np.array([(3, 7) in m for m in A4['diag']]))[0].tolist())   # identity on w


# ---------------------------------------------------------------- tuple set (restricted: w non-trivial at <= 1 copy)
wid = sorted(w_id_set)
rows = [tup for tup in itertools.product(wid, repeat=n)]
for s in range(n):
    for tup in itertools.product(*[range(N4) if r == s else wid for r in range(n)]):
        if tup[s] not in w_id_set:
            rows.append(tup)
tuples = np.array(rows, dtype=np.int64)
NT = len(tuples)
weights = N4 ** np.arange(n - 1, -1, -1)
flat_all = tuples @ weights
log('tuples', NT)
best = None
for g in glob:
    tg = g[tuples]
    for ew in itertools.product((0, 1), repeat=n):
        tw = tg.copy()
        for s in range(n):
            if ew[s]:
                tw[:, s] = wperm[tw[:, s]]
        for perm in (itertools.permutations(range(n)) if USE_SN else [tuple(range(n))]):
            key = tw[:, list(perm)] @ weights
            best = key if best is None else np.minimum(best, key)
reps, orb = np.unique(best, return_inverse=True)
NO = len(reps)
log('orbits', NO)
pos = {f: i for i, f in enumerate(flat_all.tolist())}


def orbit_sum(factors):
    """S @ (kron of per-copy N4 x c_s matrices restricted to the tuple set) as a dense (NO x prod c_s) array"""
    cs = [F.shape[1] for F in factors]
    out = np.zeros((NO, int(np.prod(cs))))
    supp = [np.nonzero(np.abs(F).sum(axis=1) > 1e-14)[0] for F in factors]
    sel = []
    for tup in itertools.product(*supp):
        f = int(np.dot(tup, weights))
        i = pos.get(f)
        if i is None:
            raise AssertionError('term leaves the tuple set')
        sel.append(i)
    sel = np.array(sel)
    tt = tuples[sel]
    oo = orb[sel]
    order = np.argsort(oo)
    tt, oo = tt[order], oo[order]
    bounds = np.searchsorted(oo, np.arange(NO + 1))
    letters = 'abcdefgh'
    ein = ','.join('k' + letters[s] for s in range(n)) + '->' + letters[:n]
    for o in range(NO):
        lo, hi = bounds[o], bounds[o + 1]
        if lo == hi:
            continue
        ops = [factors[s][tt[lo:hi, s]] for s in range(n)]
        out[o] = np.einsum(ein, *ops, optimize=True).reshape(-1)
    return out


def orbit_vec(factors):
    return orbit_sum([f.reshape(-1, 1) for f in factors])[:, 0]


lhs_c = orbit_vec([cvec] * n) - t * orbit_vec([ivec] * n)
lhs_k = np.zeros(NO)
for s in range(n):
    lhs_k += -t * orbit_vec([omyw if r == s else ivec for r in range(n)])
log('target built')


def blocks_of(alg_list):
    out = []
    for beta in itertools.product(*[al['spins'] for al in alg_list]):
        ms = [al['mj'][j] for al, j in zip(alg_list, beta)]
        rws = list(itertools.product(*[range(x) for x in ms]))
        b = len(rws)
        fl = np.zeros((b, b), dtype=np.int64)
        for R, rr in enumerate(rws):
            for C_, cc in enumerate(rws):
                f = 0
                for s, al in enumerate(alg_list):
                    f = f * len(al['coords']) + al['cidx'][(beta[s], rr[s], cc[s])]
                fl[R, C_] = f
        out.append((beta, b, fl))
    return out


A2 = dict(spins=[0, 1, 2], mj={0: 1, 1: 1, 2: 1}, coords=[(0, 0, 0), (1, 0, 0), (2, 0, 0)],
          cidx={(0, 0, 0): 0, (1, 0, 0): 1, (2, 0, 0): 2})
P2mat = np.array(P2).T
Q3 = E_xyy @ A3['Phinv']
terms = []
terms.append((orbit_sum([Q3] * n), blocks_of([A3] * n), 'Z0'))
log('Z0 map', terms[-1][0].shape)
terms.append((orbit_sum([A4['Phinv']] + [Q3] * (n - 1)), blocks_of([A4] + [A3] * (n - 1)), 'Z1'))
log('Z1 map', terms[-1][0].shape)
terms.append((orbit_sum([(E_s1 - E_s2) @ A3['Phinv']] + [P2mat] * (n - 1)), blocks_of([A3] + [A2] * (n - 1)), 'L1'))
SR = orbit_sum([A4['Ker']] + [E_xyy] * (n - 1))
log('maps built; relations', SR.shape)
np.savez(f'maps_n{n}.npz', lhs_c=lhs_c, lhs_k=lhs_k, SR=SR, **{nm: M for M, _, nm in terms})


# ---------------------------------------------------------------- forced kernels (facial reduction)
rngk = np.random.default_rng(5)


def span_proj(vecs, tol=1e-9):
    M = np.array(vecs).T
    u_, s_, _ = np.linalg.svd(M, full_matrices=False)
    r = int((s_ > s_[0] * tol).sum())
    return u_[:, :r]


def spin_bases(alg):
    K = alg['K']
    U = {}
    eps = np.zeros((3, 3, 3))
    for (i, j, k), s in [((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)]:
        eps[i, j, k] = s
    L = [-eps[k] for k in range(3)]
    tot = []
    for Lk in L:
        out = np.zeros((D ** K, D ** K))
        for s in range(K):
            f = [np.eye(3)] * K
            f[s] = Lk
            mm = f[0]
            for x in f[1:]:
                mm = np.kron(mm, x)
            out += mm
        tot.append(out)
    _, sv, vt = np.linalg.svd(tot[2])
    N0 = vt[sv < 1e-9].T
    Cas = sum(T_ @ T_ for T_ in tot)
    kk = N0.T @ Cas @ N0
    ev, evec = np.linalg.eigh((kk + kk.T) / 2)
    for j in range(K + 1):
        sel = np.abs(ev + j * (j + 1)) < 1e-6
        if sel.any():
            U[j] = N0 @ evec[:, sel]
    for ai, Mx in enumerate(alg['mats'][:3]):
        for (j, r, c) in alg['coords']:
            assert abs(U[j][:, r] @ Mx @ U[j][:, c] - alg['Phi'][alg['cidx'][(j, r, c)], ai]) < 1e-9
    return U


U3 = spin_bases(A3)
U4 = spin_bases(A4)


def kernel_in_blocks(U, Kbasis):
    out = {}
    P = Kbasis @ Kbasis.T
    for j in U:
        M = U[j].T @ P @ U[j]
        w_, v_ = np.linalg.eigh((M + M.T) / 2)
        out[j] = v_[:, w_ > 0.5]
    return out


vs3, vs4 = [], []
for _ in range(80):
    e = rngk.standard_normal(3)
    e /= np.linalg.norm(e)
    vs3.append(np.kron(np.kron(e, e), e))
    wv = rngk.standard_normal(3)
    wv -= (wv @ e) * e
    vs4.append(np.kron(np.kron(np.kron(e, e), e), wv))
K3 = span_proj(vs3)
K4 = span_proj(vs4)
log('forced kernel dims: 3 strands', K3.shape[1], '4 strands (e e e w, w perp e)', K4.shape[1])
kz3 = kernel_in_blocks(U3, K3)
kz4 = kernel_in_blocks(U4, K4)
log('per-spin kernel dims 3-strand', {j: v.shape[1] for j, v in kz3.items()}, '4-strand', {j: v.shape[1] for j, v in kz4.items()})
kz2 = {0: np.ones((1, 1)), 1: np.zeros((1, 0)), 2: np.ones((1, 1))}
for j in U3:
    if j not in kz3:
        kz3[j] = np.zeros((U3[j].shape[1], 0))
for j in U4:
    if j not in kz4:
        kz4[j] = np.zeros((U4[j].shape[1], 0))


def reduced_basis(beta, kers):
    Kb = np.ones((1, 1))
    dims = []
    for s, j in enumerate(beta):
        kj = kers[s][j]
        dims.append(kj.shape[0])
        Kb = np.kron(Kb, kj)
    b = int(np.prod(dims))
    if Kb.shape[1] == 0:
        return np.eye(b), 0
    P = np.eye(b) - Kb @ Kb.T
    w_, v_ = np.linalg.eigh(P)
    return v_[:, w_ > 0.5], Kb.shape[1]


kerlist = {'Z0': [kz3] * n, 'Z1': [kz4] + [kz3] * (n - 1), 'L1': [kz3] + [kz2] * (n - 1)}
kap = kap_fixed if kap_fixed is not None else cp.Variable()
smin = cp.Variable() if kap_fixed is not None else 0.0
cons = []
expr = 0
Xvars = []
nker = 0
for M, blocks, nm in terms:
    for beta, b, fl in blocks:
        Bb, kd = reduced_basis(beta, kerlist[nm])
        nker += kd
        br = Bb.shape[1]
        if br == 0:
            continue
        Y = cp.Variable((br, br), symmetric=True)
        cons.append(Y >> smin * np.eye(br))
        Xvars.append((nm, beta, fl, Bb, Y))
        Mb = M[:, fl.reshape(-1)]
        MB = np.einsum('oij,ia,jb->oab', Mb.reshape(-1, b, b), Bb, Bb).reshape(Mb.shape[0], -1)
        expr = expr + MB @ cp.reshape(Y, (br * br,), order='C')
log('forced kernel vectors removed:', nker)
uR, sR, _ = np.linalg.svd(SR, full_matrices=False)
SRb = uR[:, sR > 1e-9 * sR[0]]
log('relations span', SRb.shape[1], 'orbit directions')
rho = cp.Variable(SRb.shape[1])
expr = expr + SRb @ rho
scale = 4.0 ** n
cons.append(scale * expr == scale * (lhs_c + lhs_k * kap))
if kap_fixed is not None:
    cons.append(smin <= 1.0)
    prob = cp.Problem(cp.Maximize(smin), cons)
else:
    prob = cp.Problem(cp.Maximize(kap), cons)
solver = sys.argv[5] if len(sys.argv) > 5 else 'CLARABEL'
log('solving; constraints', NO, 'solver', solver)
if solver == 'SCS':
    prob.solve(solver=cp.SCS, verbose=False, max_iters=200000, eps_abs=1e-10, eps_rel=1e-10)
else:
    prob.solve(solver=cp.CLARABEL, verbose=(len(sys.argv) > 6), max_iter=500)
log('status', prob.status, 'value', prob.value)
kv = kap if kap_fixed is not None else kap.value
worst = 1e9
res = {}
for nm, beta, fl, Bb, Y in Xvars:
    Yv = (Y.value + Y.value.T) / 2
    worst = min(worst, np.linalg.eigvalsh(Yv)[0])
    res[f'{nm}_{"".join(map(str, beta))}'] = Bb @ Yv @ Bb.T
resid = np.abs(expr.value - (lhs_c + lhs_k * kv)).max()
log('kappa', kv, 'least reduced-block eigenvalue', worst, 'equality residual', resid)
np.savez(f'solfr_n{n}_{kv:.6f}.npz', rho=rho.value, kap=kv, **res)
