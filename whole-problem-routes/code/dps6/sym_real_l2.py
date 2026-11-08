"""Real second-level certificate for positivity of P_a^{(x)n}, reduced by ALL symmetries: relabelings (x legs, S4 on
the y legs) acting identically on every copy, permutations of the copies (S_n), and O(3) at every copy.

Same relaxation as real_l2.py:  C_a^{(x)n} (x) I_{y'} - t I == Z (orbit sums), Z >= 0, maximise t.
Z is taken S_n-symmetric and O(3)^n-invariant: one block per spin type (k_0, k_1, k_2, k_3) (k_j copies with spin j,
multiplicities m = 1, 3, 2, 1) and per pair of partitions (lam_1 of k_1, at most 3 rows; lam_2 of k_2, at most 2 rows),
of size d_{lam_1}(GL3) d_{lam_2}(GL2) (Gelfand-Tsetlin bases, gtbasis.py).
Constraint coefficients: G[O, type, lam] = sum over tuples tau in orbit O of (x)_j B^{(j)}_{lam_j}(tau restricted to the
copies of spin j), B^{(j)}_lam(nu) = sum_T W_T^T ((x)_i M^{(j)}_{nu_i}) W_T,  M^{(j)}_d[r, c] = Phinv[d, (j, r, c)].
usage: python sym_real_l2.py n a_num a_den [solver]"""
import itertools
import math
import sys
import time
from collections import defaultdict
import numpy as np
import cvxpy as cp
from gtbasis import gt_bases, partitions

T0 = time.time()


def log(*a):
    print(*a, f'[{time.time() - T0:.1f}s]', flush=True)


n = int(sys.argv[1])
a = int(sys.argv[2]) / int(sys.argv[3])
solver = sys.argv[4] if len(sys.argv) > 4 else 'CLARABEL'
K = 3
D = 3


# ------------------------------------------------------------------ per copy (as in real_l2.py)
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


def dmat(mm):
    M = np.zeros((D ** K, D ** K))
    for idx in itertools.product(range(D), repeat=2 * K):
        if all(idx[p] == idx[q] for p, q in mm):
            r = c = 0
            for s in range(K):
                r = r * D + idx[s]
                c = c * D + idx[K + s]
            M[r, c] = 1
    return M


DM = [dmat(mm) for mm in DIAG]
ylegs = [1, 2, 4, 5]
elems = []
for ex in (0, 1):
    for yp in itertools.permutations(ylegs):
        perm = list(range(6))
        if ex:
            perm[0], perm[3] = 3, 0
        for src, dst in zip(ylegs, yp):
            perm[src] = dst
        elems.append(np.array([IDX[canon([(perm[p], perm[q]) for p, q in mm])] for mm in DIAG], dtype=np.int64))


def tot(Lk):
    out = np.zeros((D ** K, D ** K))
    for s in range(K):
        f = [np.eye(3)] * K
        f[s] = Lk
        mat = f[0]
        for x in f[1:]:
            mat = np.kron(mat, x)
        out += mat
    return out


eps3 = np.zeros((3, 3, 3))
for (i, j, k), s in [((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)]:
    eps3[i, j, k] = s
Lt = [tot(-eps3[k]) for k in range(3)]
_, sv, vt = np.linalg.svd(Lt[2])
N0 = vt[sv < 1e-9].T
Cas = sum(L @ L for L in Lt)
kk = N0.T @ Cas @ N0
ev, evec = np.linalg.eigh((kk + kk.T) / 2)
U = {}
for j in range(K + 1):
    sel = np.abs(ev + j * (j + 1)) < 1e-6
    U[j] = N0 @ evec[:, sel]
SP = [0, 1, 2, 3]
mj = {j: U[j].shape[1] for j in SP}
assert [mj[j] for j in SP] == [1, 3, 2, 1]
coords = [(j, r, c) for j in SP for r in range(mj[j]) for c in range(mj[j])]
cidx = {c: i for i, c in enumerate(coords)}
Phi = np.zeros((ND, ND))
for ai, M in enumerate(DM):
    for i, (j, r, c) in enumerate(coords):
        Phi[i, ai] = U[j][:, r] @ M @ U[j][:, c]
Phinv = np.linalg.inv(Phi)
Mj = {j: np.array([[[Phinv[d, cidx[(j, r, c)]] for c in range(mj[j])] for r in range(mj[j])] for d in range(ND)])
      for j in SP}                                               # Mj[j][d] is mj x mj


def vec(dct):
    v = np.zeros(ND)
    for k_, c in dct.items():
        v[IDX[canon(k_)]] += c
    return v


cvec = vec({((0, 3), (1, 4), (2, 5)): a, ((0, 4), (1, 3), (2, 5)): -0.25, ((0, 1), (3, 4), (2, 5)): -0.25})
ivec = vec({((0, 3), (1, 4), (2, 5)): 1.0})

# ------------------------------------------------------------------ orbits of multisets of size n
def multisets(k):
    if k == 0:
        return np.zeros((1, 0), dtype=np.int64)
    return np.array(list(itertools.combinations_with_replacement(range(ND), k)), dtype=np.int64)


def ordcount(ms):
    """number of distinct orderings of each multiset (rows of a sorted int array)"""
    out = np.full(len(ms), math.factorial(ms.shape[1]), dtype=np.float64)
    for d in range(ND):
        cnt = (ms == d).sum(axis=1)
        out /= np.array([math.factorial(int(c)) for c in cnt])
    return out


W_ENC = (ND ** np.arange(n - 1, -1, -1)).astype(np.int64)
MS = multisets(n)
mskey = MS @ W_ENC
best = None
for g in elems:
    kk_ = np.sort(g[MS], axis=1) @ W_ENC
    best = kk_ if best is None else np.minimum(best, kk_)
reps, orb_of_ms = np.unique(best, return_inverse=True)
NO = len(reps)
key_order = np.argsort(mskey)
sorted_keys = mskey[key_order]


def orbit_of_rows(rows):
    """orbit index of multisets given as (unsorted) rows of diagram indices"""
    keys = np.sort(rows, axis=1) @ W_ENC
    pos = np.searchsorted(sorted_keys, keys)
    return orb_of_ms[key_order[pos]]


oc = ordcount(MS)
Sc = np.bincount(orb_of_ms, weights=oc * np.prod(cvec[MS], axis=1), minlength=NO)
Si = np.bincount(orb_of_ms, weights=oc * np.prod(ivec[MS], axis=1), minlength=NO)
log(f'n = {n}, a = {a}: multisets {len(MS)}, orbits {NO}')

# ------------------------------------------------------------------ B^{(j)}_lam(nu) for spins 1 (m=3) and 2 (m=2)
GT = {}


def gt(m, k):
    if (m, k) not in GT:
        GT[(m, k)] = gt_bases(m, k)
    return GT[(m, k)]


def B_blocks(j, k, ms):
    """for multisets ms (rows, size k) of diagrams: dict lam -> array (len(ms), d, d)"""
    m = mj[j]
    if k == 0:
        return {(): np.ones((len(ms), 1, 1))}
    bases = gt(m, k)
    lams = list(bases)
    Wall = np.hstack([bases[l]['W'][T] for l in lams for T in bases[l]['tabs']])      # (m^k, m^k)
    out = {l: np.zeros((len(ms), bases[l]['d'], bases[l]['d'])) for l in lams}
    batch = max(1, int(2e7 // (m ** (2 * k))))
    for st in range(0, len(ms), batch):
        sub = ms[st:st + batch]
        nb = len(sub)
        Y = np.broadcast_to(Wall.reshape([m] * k + [m ** k]), (nb,) + tuple([m] * k) + (m ** k,)).copy()
        for leg in range(k):
            Mleg = Mj[j][sub[:, leg]]                                  # (nb, m, m)
            Y = np.moveaxis(np.einsum('bij,bj...->bi...', Mleg, np.moveaxis(Y, leg + 1, 1)), 1, leg + 1)
        Y = Y.reshape(nb, m ** k, m ** k)
        col = 0
        for l in lams:
            d = bases[l]['d']
            for T in bases[l]['tabs']:
                WT = bases[l]['W'][T]
                out[l][st:st + nb] += np.einsum('ip,bij->bpj', WT, Y[:, :, col:col + d], optimize=True)
                col += d
    return out


# ------------------------------------------------------------------ assemble the constraints
types = [tk for tk in itertools.product(range(n + 1), repeat=4) if sum(tk) == n]
blocks = []           # (type, lam1, lam2, d, G array (NO, d, d), multiplicity factor)
for tk in types:
    k0, k1, k2, k3 = tk
    lam1s = partitions(k1, 3) if k1 else [()]
    lam2s = partitions(k2, 2) if k2 else [()]
    ms0, ms1, ms2, ms3 = multisets(k0), multisets(k1), multisets(k2), multisets(k3)
    s0 = np.prod(Mj[0][ms0][:, :, 0, 0], axis=1)                     # spin 0 and 3: scalar factors
    s3 = np.prod(Mj[3][ms3][:, :, 0, 0], axis=1)
    B1 = B_blocks(1, k1, ms1)
    B2 = B_blocks(2, k2, ms2)
    o0, o1, o2, o3 = ordcount(ms0), ordcount(ms1), ordcount(ms2), ordcount(ms3)
    G = {(l1, l2): np.zeros((NO, B1[l1].shape[1] * B2[l2].shape[1], B1[l1].shape[1] * B2[l2].shape[1]))
         for l1 in lam1s for l2 in lam2s}
    sizes = (len(ms0), len(ms1), len(ms2), len(ms3))
    total = int(np.prod(sizes))
    CH = 400000
    for st in range(0, total, CH):
        lin = np.arange(st, min(total, st + CH))
        I0, I1, I2, I3 = np.unravel_index(lin, sizes)
        rows = np.hstack([ms0[I0], ms1[I1], ms2[I2], ms3[I3]])
        orbs = orbit_of_rows(rows)
        wgt = o0[I0] * o1[I1] * o2[I2] * o3[I3] * s0[I0] * s3[I3]
        for l1 in lam1s:
            for l2 in lam2s:
                d = G[(l1, l2)].shape[1]
                acc = G[(l1, l2)].reshape(NO, d * d)
                sub = max(1, int(2.5e7 // (d * d)))          # at most ~200 MB per contribution array
                for s2 in range(0, len(lin), sub):
                    sl = slice(s2, s2 + sub)
                    contrib = np.einsum('b,bij,bkl->bikjl', wgt[sl], B1[l1][I1[sl]], B2[l2][I2[sl]]).reshape(-1, d * d)
                    for e in range(d * d):
                        acc[:, e] += np.bincount(orbs[sl], weights=contrib[:, e], minlength=NO)
    mult = math.factorial(n) / (math.factorial(k0) * math.factorial(k1) * math.factorial(k2) * math.factorial(k3))
    for (l1, l2), Garr in G.items():
        blocks.append((tk, l1, l2, Garr.shape[1], Garr, mult))
    log('type', tk, 'done; blocks so far', len(blocks))
log('blocks', len(blocks), 'max size', max(b[3] for b in blocks), 'variables', sum(b[3] * (b[3] + 1) // 2 for b in blocks))

# ------------------------------------------------------------------ save the a-independent data
if len(sys.argv) > 5 and sys.argv[5] == 'save':
    np.savez_compressed(f'symdata_n{n}.npz', MS=MS, orb_of_ms=orb_of_ms, oc=oc, Si=Si, NO=NO,
                        types=np.array([b[0] for b in blocks]), mults=np.array([b[5] for b in blocks]),
                        dims=np.array([b[3] for b in blocks]), lam1=np.array([str(b[1]) for b in blocks]),
                        lam2=np.array([str(b[2]) for b in blocks]),
                        **{f'G{i}': b[4] for i, b in enumerate(blocks)})
    log('saved', f'symdata_n{n}.npz')
    sys.exit(0)
# ------------------------------------------------------------------ SDP
t = cp.Variable()
expr = 0
cons = []
for tk, l1, l2, d, Garr, mult in blocks:
    Zb = cp.Variable((d, d), symmetric=True)
    cons.append(Zb >> 0)
    expr = expr + mult * (Garr.reshape(NO, d * d) @ cp.reshape(Zb, (d * d,), order='C'))
cons.append(expr == Sc - t * Si)
prob = cp.Problem(cp.Maximize(t), cons)
if solver == 'SCS':
    prob.solve(solver=cp.SCS, eps=1e-9, max_iters=200000)
else:
    prob.solve(solver=cp.CLARABEL)
an = (1 + 2 ** (1 - 1 / n)) / 4
log(f'n = {n}, a = {a:.6f}: status {prob.status}, best t = {t.value:.6e}  (4^-n = {4.0 ** -n:.6e}, a_n = {an:.6f})')
