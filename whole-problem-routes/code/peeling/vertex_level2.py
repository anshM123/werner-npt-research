"""Real second-level certificate for block-positivity of Choi(Sym_2) (x) C_a^{(x)m} on real product vectors.

Real 2-positivity of P_a^{(x)m} (a = 3/4: the vertex of Conjecture QT_m, the core of the map route) is positivity of
Sym_2 (x) P_a^{(x)m} on real states, with Sym_2(X) = (X + X^T)/2 on a rebit (Choi (F + Omega)/2).  Copy 0 is the rebit
(d = 2), copies 1..m are qutrits (d = 3).  Per copy the strands are (x, y, y') = (0, 1, 2) (points 0..2 top, 3..5
bottom).  Certificate:  target (x) I_{y'} - t I == Z + R  (mod orbit sums),  Z >= 0 (O(2) x O(3)^m invariant, block
form),  R a combination of the 5 linear relations among the 15 three-strand Brauer diagrams at d = 2 (rebit copy).
Orbits: exchanging the two x legs and permuting the four y legs, identically at every copy; permuting qutrit copies.
usage: python vertex_level2.py m a_num a_den [rebit|qutrit]   (qutrit: copy 0 is a qutrit at a = 1/2, the full vertex)"""
import itertools
import sys
import time
import numpy as np
import cvxpy as cp

T0 = time.time()


def log(*a):
    print(*a, f'[{time.time() - T0:.1f}s]', flush=True)


m = int(sys.argv[1])
a = int(sys.argv[2]) / int(sys.argv[3])
kind0 = sys.argv[4] if len(sys.argv) > 4 else 'rebit'
n = m + 1
K = 3


def matchings(points):
    if not points:
        yield []
        return
    p = points[0]
    for i in range(1, len(points)):
        rest = points[1:i] + points[i + 1:]
        for mm in matchings(rest):
            yield [(p, points[i])] + mm


def canon(mm):
    return tuple(sorted(tuple(sorted(p)) for p in mm))


IDENT = canon([(s, s + K) for s in range(K)])
DIAG = [IDENT] + [mm for mm in sorted(set(canon(mm) for mm in matchings(list(range(2 * K))))) if mm != IDENT]
ND = len(DIAG)
IDX = {mm: i for i, mm in enumerate(DIAG)}


def dmat(mm, D):
    M = np.zeros((D ** K, D ** K))
    for idx in itertools.product(range(D), repeat=2 * K):
        if all(idx[p] == idx[q] for p, q in mm):
            r = c = 0
            for s in range(K):
                r = r * D + idx[s]
                c = c * D + idx[K + s]
            M[r, c] = 1
    return M


def act(mm, perm):
    return IDX[canon([(perm[p], perm[q]) for p, q in mm])]


elems = []
ylegs = [1, 2, 4, 5]
for ex in (0, 1):
    for yp in itertools.permutations(ylegs):
        perm = list(range(6))
        if ex:
            perm[0], perm[3] = 3, 0
        for src, dst in zip(ylegs, yp):
            perm[src] = dst
        elems.append(np.array([act(mm, perm) for mm in DIAG]))


def tot(Lk, D):
    out = np.zeros((D ** K, D ** K))
    for s in range(K):
        f = [np.eye(D)] * K
        f[s] = Lk
        mat = f[0]
        for x in f[1:]:
            mat = np.kron(mat, x)
        out += mat
    return out


def copy_data(D):
    """block vectors U (one m=0-type vector per irrep copy), Phi (block coords x diagrams), Phinv (right inverse),
    kernel of the diagram -> operator map"""
    DM = [dmat(mm, D) for mm in DIAG]
    Vmat = np.array([M.reshape(-1) for M in DM]).T
    _, sv, vt = np.linalg.svd(Vmat)
    rank = int((sv > 1e-9).sum())
    kern = vt[rank:].T
    U = {}
    if D == 3:
        eps = np.zeros((3, 3, 3))
        for (i, j, k), s in [((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)]:
            eps[i, j, k] = s
        Lt = [tot(-eps[k], 3) for k in range(3)]
        _, s2, v2 = np.linalg.svd(Lt[2])
        N0 = v2[s2 < 1e-9].T
        Cas = sum(L @ L for L in Lt)
        kk = N0.T @ Cas @ N0
        ev, evec = np.linalg.eigh((kk + kk.T) / 2)
        for j in range(K + 1):
            sel = np.abs(ev + j * (j + 1)) < 1e-6
            if sel.any():
                U[j] = N0 @ evec[:, sel]
    else:
        J = np.array([[0., -1.], [1., 0.]])
        L = tot(J, 2)
        L2 = L @ L
        Rf = np.diag([1., -1.])
        Rt = Rf
        for _ in range(K - 1):
            Rt = np.kron(Rt, Rf)
        ev, evec = np.linalg.eigh((L2 + L2.T) / 2)
        for k in range(K + 1):
            sel = np.abs(ev + k * k) < 1e-6
            if not sel.any():
                continue
            Vk = evec[:, sel]
            r_ev, r_vec = np.linalg.eigh(Vk.T @ Rt @ Vk)
            U[k] = Vk @ r_vec[:, r_ev > 0.5]
    SP = sorted(U)
    mj = {j: U[j].shape[1] for j in SP}
    coords = [(j, r, c) for j in SP for r in range(mj[j]) for c in range(mj[j])]
    Phi = np.zeros((len(coords), ND))
    for ai, M in enumerate(DM):
        for i, (j, r, c) in enumerate(coords):
            Phi[i, ai] = U[j][:, r] @ M @ U[j][:, c]
    Phinv = np.linalg.pinv(Phi)
    assert np.allclose(Phi @ Phinv, np.eye(len(coords)), atol=1e-9)
    # sanity: invariant operators are determined by their block coordinates (rank check)
    assert np.linalg.matrix_rank(Phi, tol=1e-9) == rank, (np.linalg.matrix_rank(Phi, tol=1e-9), rank)
    return dict(D=D, U=U, SP=SP, mj=mj, coords=coords, cidx={c: i for i, c in enumerate(coords)}, Phi=Phi,
                Phinv=Phinv, kern=kern, rank=rank)


CD = {2: copy_data(2), 3: copy_data(3)}
dims = [2 if kind0 == 'rebit' else 3] + [3] * m
log('rebit kernel dim', CD[2]['kern'].shape[1], 'rank', CD[2]['rank'], 'sectors', CD[2]['mj'], '| qutrit sectors', CD[3]['mj'])

# orbits of diagram tuples
NT = ND ** n
tuples = np.array(list(itertools.product(range(ND), repeat=n)), dtype=np.int64)
weights = ND ** np.arange(n - 1, -1, -1)
best = None
for g in elems:
    tg = g[tuples]
    for perm in itertools.permutations(range(1, n)):
        cols = [0] + list(perm)
        key = tg[:, cols] @ weights
        best = key if best is None else np.minimum(best, key)
reps, orb = np.unique(best, return_inverse=True)
NO = len(reps)
log('coordinates', NT, 'orbits', NO)

Pinvs = [CD[d]['Phinv'] for d in dims]


def modeprod_list(mats, X):
    T = X.reshape([ND] * n)
    for s in range(n):
        T = np.moveaxis(np.tensordot(mats[s], T, axes=([1], [s])), 0, s)
    return T.reshape(-1)


blocksizes = [CD[d]['Phinv'].shape[1] for d in dims]
NB = int(np.prod(blocksizes))
Morb = np.zeros((NO, NB))
order = np.argsort(orb)
bounds = np.searchsorted(orb[order], np.arange(NO + 1))
for o in range(NO):
    ind = np.zeros(NT)
    ind[order[bounds[o]:bounds[o + 1]]] = 1.0
    T = ind.reshape([ND] * n)
    for s in range(n):
        T = np.moveaxis(np.tensordot(Pinvs[s].T, T, axes=([1], [s])), 0, s)
    Morb[o] = T.reshape(-1)
log('orbit map built', Morb.shape)


def vec(dct):
    v = np.zeros(ND)
    for kk, c in dct.items():
        v[IDX[canon(kk)]] += c
    return v


F_AB = ((0, 4), (1, 3), (2, 5))
OM_AB = ((0, 1), (3, 4), (2, 5))
ID = ((0, 3), (1, 4), (2, 5))
cq = vec({ID: a, F_AB: -0.25, OM_AB: -0.25})
crebit = vec({F_AB: 0.5, OM_AB: 0.5})
cvertex = vec({ID: 0.5, F_AB: -0.25, OM_AB: -0.25})
ivec = vec({ID: 1.0})


def kron_all(vs):
    out = np.array([1.0])
    for v in vs:
        out = np.kron(out, v)
    return out


first = crebit if kind0 == 'rebit' else cvertex
Sc = np.bincount(orb, weights=kron_all([first] + [cq] * m), minlength=NO)
Si = np.bincount(orb, weights=kron_all([ivec] * n), minlength=NO)

# block variables: one PSD block per sector tuple
blocks = []
cd_list = [CD[d] for d in dims]
for beta in itertools.product(*[cd['SP'] for cd in cd_list]):
    ms = [cd['mj'][j] for cd, j in zip(cd_list, beta)]
    rows = list(itertools.product(*[range(x) for x in ms]))
    b = len(rows)
    flat = np.zeros((b, b), dtype=np.int64)
    for R_, rr in enumerate(rows):
        for C_, cc in enumerate(rows):
            f = 0
            for s in range(n):
                f = f * blocksizes[s] + cd_list[s]['cidx'][(beta[s], rr[s], cc[s])]
            flat[R_, C_] = f
    blocks.append((beta, b, flat))

t = cp.Variable()
expr = 0
cons = []
for beta, b, flat in blocks:
    Zb = cp.Variable((b, b), symmetric=True)
    cons.append(Zb >> 0)
    expr = expr + Morb[:, flat.reshape(-1)] @ cp.reshape(Zb, (b * b,), order='C')
# relation slack for the rebit copy: kern (x) anything, pushed to orbit sums
if kind0 == 'rebit':
    kern = CD[2]['kern']
    nk = kern.shape[1]
    # orbit-sum matrix in diagram coordinates restricted to (kern at copy 0) (x) (all diagrams elsewhere)
    Odiag = np.zeros((NO, NT))
    Odiag[orb, np.arange(NT)] = 1.0
    Odiag = Odiag.reshape(NO, ND, ND ** (n - 1))
    Mrel = np.einsum('oar,ak->okr', Odiag, kern).reshape(NO, nk * ND ** (n - 1))
    _, sr, vr = np.linalg.svd(Mrel, full_matrices=False)
    keep = sr > 1e-9 * sr[0]
    Brel = (Mrel @ vr[keep].T)
    w = cp.Variable(int(keep.sum()))
    expr = expr + Brel @ w
    log('relation slack dims', int(keep.sum()))
cons.append(expr == Sc - t * Si)
prob = cp.Problem(cp.Maximize(t), cons)
try:
    prob.solve(solver=cp.CLARABEL)
except Exception as ex:
    log('CLARABEL failed:', ex, '-> SCS')
    prob.solve(solver=cp.SCS, eps=1e-8, max_iters=200000)
log(f'm = {m} qutrit copies + {kind0} copy, a = {a:.6f}: status {prob.status}, best t = {t.value:.6e}')
