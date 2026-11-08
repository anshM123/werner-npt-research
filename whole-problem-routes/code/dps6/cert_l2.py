"""Exact certificate for positivity of P_a^{(x)n} from a solution of the reduced level-2 SDP (sym_real_l2.py, sym_solve.py).

  python cert_l2.py n anum aden extract MODE      diagram coefficients z (one per multiset of diagrams) from
                                                 symsol_n{n}_{anum}_{aden}_{MODE}.npz; checks the orbit identity in floats
  python cert_l2.py n anum aden round BITS [t]   rational z: rounded to 2^-BITS, then the orbit residual is put on one
                                                 multiset per orbit, so that the orbit identity holds exactly (target t)
  python cert_l2.py n anum aden check            exact PSD check of Z = sum_tau z_tau D_tau: one block per spin type,
                                                 integer bases V_j (exact_alg.py), Young-symmetrizer images for S_k1 x S_k2,
                                                 exact LDL^T of every block

Z (x) is an operator on (x_i, y_i, y'_i)_{i<=n}; the identity  (x(x)y(x)y)^T Z (x(x)y(x)y) = F(x,y)|y|^2 - t|x|^2|y|^4
with F(x,y) = <y|P_a^{(x)n}(x x^T)|y> holds iff the orbit sums agree (relabelings x-legs / S_4 on y-legs acting on every
copy at once, and the copy permutations)."""
import ast
import itertools
import math
import os
import pickle
import sys
import time
from fractions import Fraction
import numpy as np

T0 = time.time()


def log(*a):
    print(*a, f'[{time.time() - T0:.1f}s]', flush=True)


n, anum, aden = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
stage = sys.argv[4]
a_fr = Fraction(anum, aden)
K = 3
D = 3


# ------------------------------------------------------------------ per copy, verbatim from sym_real_l2.py
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
      for j in SP}


def cvec_exact(a):
    v = [Fraction(0)] * ND
    v[IDX[canon(((0, 3), (1, 4), (2, 5)))]] += a
    v[IDX[canon(((0, 4), (1, 3), (2, 5)))]] -= Fraction(1, 4)
    v[IDX[canon(((0, 1), (3, 4), (2, 5)))]] -= Fraction(1, 4)
    return v


# ------------------------------------------------------------------ multisets and orbits (as in sym_real_l2.py)
def multisets(k):
    if k == 0:
        return np.zeros((1, 0), dtype=np.int64)
    return np.array(list(itertools.combinations_with_replacement(range(ND), k)), dtype=np.int64)


W_ENC = (ND ** np.arange(n - 1, -1, -1)).astype(np.int64)
MS = multisets(n)
NMS = len(MS)
mskey = MS @ W_ENC
best = None
for g in elems:
    kk_ = np.sort(g[MS], axis=1) @ W_ENC
    best = kk_ if best is None else np.minimum(best, kk_)
reps, orb_of_ms = np.unique(best, return_inverse=True)
NO = len(reps)
OC = [math.factorial(n) // math.prod(math.factorial(int((row == d).sum())) for d in range(ND)) for row in MS]
oc = np.array(OC, dtype=np.float64)
key_order = np.argsort(mskey)
sorted_keys = mskey[key_order]


def full_index():
    """multiset index of every ordered tuple tau in [ND]^n (flattened, C order)"""
    grids = np.indices((ND,) * n).reshape(n, -1).T
    keys = np.sort(grids, axis=1) @ W_ENC
    return key_order[np.searchsorted(sorted_keys, keys)].astype(np.int32)


types = [tk for tk in itertools.product(range(n + 1), repeat=4) if sum(tk) == n]


def mult_of(tk):
    return math.factorial(n) // math.prod(math.factorial(x) for x in tk)


def target_exact(t):
    cv = cvec_exact(a_fr)
    Sc = [Fraction(0)] * NO
    Si = [Fraction(0)] * NO
    for i, row in enumerate(MS):
        p = Fraction(1)
        for d in row:
            p *= cv[d]
        Sc[orb_of_ms[i]] += OC[i] * p
        if all(d == 0 for d in row):
            Si[orb_of_ms[i]] += OC[i]
    return [Sc[o] - t * Si[o] for o in range(NO)]


# ------------------------------------------------------------------ stage: extract
if stage == 'extract':
    from gtbasis import gt_bases
    mode = sys.argv[5]
    sol = np.load(f'symsol_n{n}_{anum}_{aden}_{mode}.npz')
    Dd = np.load(f'symdata_n{n}.npz')
    btypes = [tuple(int(x) for x in r) for r in Dd['types']]
    lam1 = [ast.literal_eval(s) for s in Dd['lam1']]
    lam2 = [ast.literal_eval(s) for s in Dd['lam2']]
    nb = len(btypes)
    tval = float(sol['t'])
    GT = {}

    def gt(m, k):
        if (m, k) not in GT:
            GT[(m, k)] = gt_bases(m, k)
        return GT[(m, k)]

    fidx = full_index()
    zacc = np.zeros(NMS)
    for tk in types:
        k0, k1, k2, k3 = tk
        b1, b2 = gt(3, k1), gt(2, k2)
        dimx = 3 ** k1 * 2 ** k2
        X = np.zeros((dimx, dimx))
        for i in range(nb):
            if btypes[i] != tk:
                continue
            Zb = sol[f'Z{i}']
            l1, l2 = lam1[i], lam2[i]
            for T1 in b1[l1]['tabs']:
                for T2 in b2[l2]['tabs']:
                    W = np.kron(b1[l1]['W'][T1], b2[l2]['W'][T2])
                    X += W @ Zb @ W.T
        # c(tau) = < (x)_i M^{(j_i)}_{tau_i}, X >, copies ordered spin 0 | 1 | 2 | 3
        kq = k1 + k2
        ms_ = [3] * k1 + [2] * k2
        Tn = X.reshape(ms_ + ms_)
        Tn = Tn.transpose([x for i in range(kq) for x in (i, i + kq)]) if kq else Tn
        Tn = Tn.reshape([m * m for m in ms_]) if kq else Tn.reshape(())
        for i in range(kq):
            Mfl = Mj[1 if i < k1 else 2].reshape(ND, -1)                 # (15, m^2)
            Tn = np.moveaxis(np.tensordot(Tn, Mfl, axes=([i], [1])), -1, i)
        s0 = Mj[0][:, 0, 0]
        s3 = Mj[3][:, 0, 0]
        full = np.ones(())
        for _ in range(k0):
            full = np.multiply.outer(full, s0)
        full = np.multiply.outer(full, Tn)
        for _ in range(k3):
            full = np.multiply.outer(full, s3)
        zacc += mult_of(tk) * np.bincount(fidx, weights=full.reshape(-1), minlength=NMS)
    z = zacc / oc
    # orbit identity in floats
    cv = np.array([float(x) for x in cvec_exact(a_fr)])
    Sc = np.bincount(orb_of_ms, weights=oc * np.prod(cv[MS], axis=1), minlength=NO)
    Si = np.bincount(orb_of_ms, weights=oc * np.all(MS == 0, axis=1), minlength=NO)
    lhs = np.bincount(orb_of_ms, weights=oc * z, minlength=NO)
    log(f'n = {n}, a = {anum}/{aden}, mode {mode}: t = {tval:.6e}; |z| max {np.abs(z).max():.3e}; '
        f'orbit identity residual {np.abs(lhs - (Sc - tval * Si)).max():.3e} (max |Sc| {np.abs(Sc).max():.2e})')
    np.save(f'zfloat_n{n}_{anum}_{aden}_{mode}.npy', z)
    sys.exit(0)

# ------------------------------------------------------------------ stage: round
if stage == 'round':
    bits = int(sys.argv[5])
    tt = Fraction(sys.argv[6]) if len(sys.argv) > 6 else Fraction(0)
    mode = sys.argv[7] if len(sys.argv) > 7 else 'center'
    if mode == 'mix':        # theta * (maxt solution) + (1 - theta) * (center solution at t = 0), theta = t / t_max
        tmax = float(np.load(f'symsol_n{n}_{anum}_{aden}_maxt.npz')['t'])
        assert float(np.load(f'symsol_n{n}_{anum}_{aden}_center.npz')['t']) == 0.0, 'mix needs the center solution at t = 0'
        theta = float(tt) / tmax
        z = theta * np.load(f'zfloat_n{n}_{anum}_{aden}_maxt.npy') + (1 - theta) * np.load(f'zfloat_n{n}_{anum}_{aden}_center.npy')
        log(f'mix: t_max = {tmax:.9e}, theta = {theta:.6f}')
    else:
        z = np.load(f'zfloat_n{n}_{anum}_{aden}_{mode}.npy')
    den0 = 2 ** bits
    zq = [Fraction(int(round(float(x) * den0)), den0) for x in z]
    tgt = target_exact(tt)
    lhs = [Fraction(0)] * NO
    for i in range(NMS):
        lhs[orb_of_ms[i]] += OC[i] * zq[i]
    # representative of each orbit: the multiset with the largest number of orderings
    rep = {}
    for i in range(NMS):
        o = orb_of_ms[i]
        if o not in rep or OC[i] > OC[rep[o]]:
            rep[o] = i
    maxres = 0.0
    for o in range(NO):
        r = tgt[o] - lhs[o]
        maxres = max(maxres, abs(float(r)))
        zq[rep[o]] += r / OC[rep[o]]
    # symmetrize under transposition of every diagram (Z -> (Z + Z^T)/2; same orbit, same number of orderings)
    tp = [3, 4, 5, 0, 1, 2]
    TR = np.array([IDX[canon([(tp[p], tp[q]) for p, q in mm])] for mm in DIAG], dtype=np.int64)
    keysT = np.sort(TR[MS], axis=1) @ W_ENC
    idxT = key_order[np.searchsorted(sorted_keys, keysT)]
    assert all(orb_of_ms[i] == orb_of_ms[idxT[i]] and OC[i] == OC[idxT[i]] for i in range(NMS))
    zq = [(zq[i] + zq[idxT[i]]) / 2 for i in range(NMS)]
    # recheck exactly
    lhs = [Fraction(0)] * NO
    for i in range(NMS):
        lhs[orb_of_ms[i]] += OC[i] * zq[i]
    assert all(lhs[o] == tgt[o] for o in range(NO))
    L = 1
    for x in zq:
        L = L * x.denominator // math.gcd(L, x.denominator)
    num = [int(x * L) for x in zq]
    pickle.dump(dict(n=n, a=(anum, aden), t=(tt.numerator, tt.denominator), den=L, num=num), open(f'zexact_n{n}_{anum}_{aden}.pkl', 'wb'))
    log(f'rounded to 2^-{bits}; max orbit residual before correction {maxres:.3e}; common denominator {L} '
        f'(bits {L.bit_length()}); orbit identity exact for t = {tt}')
    sys.exit(0)


# ------------------------------------------------------------------ stage: export (self-describing JSON)
if stage == 'export':
    import json
    Z = pickle.load(open(f'zexact_n{n}_{anum}_{aden}.pkl', 'rb'))
    out = dict(n=n, a=[anum, aden], t=list(Z['t']),
               points='top x, y, yp = 0, 1, 2; bottom x, y, yp = 3, 4, 5',
               diagrams=[[list(p) for p in mm] for mm in DIAG], den=str(Z['den']),
               coefficients=[[list(map(int, MS[i])), str(Z['num'][i])] for i in range(NMS)])
    fn = sys.argv[5] if len(sys.argv) > 5 else f'cert_n{n}_{anum}_{aden}.json'
    if fn.endswith('.xz'):
        import lzma
        with lzma.open(fn, 'wt') as fh:
            json.dump(out, fh)
    else:
        json.dump(out, open(fn, 'w'))
    log('wrote', fn)
    sys.exit(0)


# ------------------------------------------------------------------ stage: check (exact)
def young_image(m, lam):
    """integer basis (columns) of the image of the Young symmetrizer of the row-reading tableau of shape lam on
    (R^m)^{(x)k}: a copy V_lam (x) u of the GL_m irreducible; one column per semistandard filling."""
    k = sum(lam)
    if k == 0:
        return np.ones((1, 1), dtype=object)
    T0 = []
    e = 0
    for r, l in enumerate(lam):
        T0.append(list(range(e, e + l)))
        e += l
    cols = [[T0[r][c] for r in range(len(lam)) if c < lam[r]] for c in range(lam[0])]
    rowperms = [list(itertools.permutations(row)) for row in T0]
    colperms = [list(itertools.permutations(col)) for col in cols]

    def sign(src, dst):
        p = [src.index(x) for x in dst]
        s = 1
        seen = [False] * len(p)
        for i in range(len(p)):
            if not seen[i]:
                j, L = i, 0
                while not seen[j]:
                    seen[j] = True
                    j = p[j]
                    L += 1
                s *= (-1) ** (L - 1)
        return s

    def apply_rows(vec):            # vec: dict multi-index -> coeff; symmetrize over row group
        out = {}
        for combo in itertools.product(*rowperms):
            perm = list(range(k))
            for row, img in zip(T0, combo):
                for s_, d_ in zip(row, img):
                    perm[s_] = d_
            for I, c in vec.items():
                J = [0] * k
                for s_ in range(k):
                    J[perm[s_]] = I[s_]
                J = tuple(J)
                out[J] = out.get(J, 0) + c
        return out

    def apply_cols(vec):
        out = {}
        for combo in itertools.product(*colperms):
            perm = list(range(k))
            sg = 1
            for col, img in zip(cols, combo):
                for s_, d_ in zip(col, img):
                    perm[s_] = d_
                sg *= sign(list(col), list(img))
            for I, c in vec.items():
                J = [0] * k
                for s_ in range(k):
                    J[perm[s_]] = I[s_]
                J = tuple(J)
                out[J] = out.get(J, 0) + sg * c
        return out

    # semistandard fillings with entries < m
    shape_cells = [(r, c) for r in range(len(lam)) for c in range(lam[r])]
    ssyt = []

    def rec(i, fill):
        if i == len(shape_cells):
            ssyt.append(dict(fill))
            return
        r, c = shape_cells[i]
        lo = 0
        if c > 0:
            lo = max(lo, fill[(r, c - 1)])
        if r > 0:
            lo = max(lo, fill[(r - 1, c)] + 1)
        for v in range(lo, m):
            fill[(r, c)] = v
            rec(i + 1, fill)
            del fill[(r, c)]
    rec(0, {})
    vecs = []
    for S in ssyt:
        I = [0] * k
        for (r, c), v in S.items():
            I[T0[r][c]] = v
        vec = apply_cols(apply_rows({tuple(I): 1}))
        col = np.zeros(m ** k, dtype=object)
        col[:] = 0
        for J, cf in vec.items():
            if cf:
                col[int(np.ravel_multi_index(J, (m,) * k))] += cf
        vecs.append(col)
    W = np.array(vecs, dtype=object).T if vecs else np.zeros((m ** k, 0), dtype=object)
    return W


def ldl_pivots(A):
    """exact LDL^T (no pivoting) of a symmetric rational matrix; returns the list of pivots (stops at a non-positive)"""
    d = len(A)
    M = [[Fraction(A[i][j]) for j in range(d)] for i in range(d)]
    piv = []
    for k in range(d):
        p = M[k][k]
        piv.append(p)
        if p <= 0:
            return piv
        for i in range(k + 1, d):
            f = M[i][k] / p
            if f:
                Mi, Mk = M[i], M[k]
                for j in range(k + 1, i + 1):
                    Mi[j] -= f * Mk[j]
        for i in range(k + 1, d):
            for j in range(k + 1, i):
                M[j][i] = M[i][j]
    return piv


if stage == 'check':
    import flint
    from gtbasis import partitions
    here = os.path.dirname(os.path.abspath(__file__))
    cands = [os.path.join(here, 'exact_alg_k3.pkl'), os.path.join(here, '..', 'tc3', 'exact_alg_k3.pkl'), 'exact_alg_k3.pkl']
    pk = next((c for c in cands if os.path.exists(c)), None)
    if pk is None:                               # build the integer spin bases with ../tc3/exact_alg.py
        sys.path.insert(0, os.path.join(here, '..', 'tc3'))
        import exact_alg
        A = exact_alg.build(3)
    else:
        A = pickle.load(open(pk, 'rb'))
    assert [tuple(m) for m in A['diags']] == [tuple(m) for m in DIAG]
    V = A['V']
    PhiE = {j: np.array([(V[j].T.astype(object) @ np.array(A['mats'][d], dtype=object) @ V[j].astype(object))
                         for d in range(ND)], dtype=object) for j in SP}          # PhiE[j][d] is m_j x m_j, integer
    if len(sys.argv) > 5:                       # read the exported certificate (JSON, optionally .xz)
        import json
        if sys.argv[5].endswith('.xz'):
            import lzma
            Cj = json.load(lzma.open(sys.argv[5], 'rt'))
        else:
            Cj = json.load(open(sys.argv[5]))
        assert Cj['n'] == n and Fraction(*Cj['a']) == a_fr
        assert [tuple(tuple(p) for p in mm) for mm in Cj['diagrams']] == [tuple(mm) for mm in DIAG]
        assert [tuple(r) for r, _ in Cj['coefficients']] == [tuple(int(v) for v in row) for row in MS]
        Z = dict(t=tuple(Cj['t']), den=int(Cj['den']), num=[int(v) for _, v in Cj['coefficients']])
    else:
        Z = pickle.load(open(f'zexact_n{n}_{anum}_{aden}.pkl', 'rb'))
    # the orbit identity, exactly
    tgt = target_exact(Fraction(*Z['t']))
    lhs = [Fraction(0)] * NO
    for i in range(NMS):
        lhs[orb_of_ms[i]] += OC[i] * Fraction(Z['num'][i], Z['den'])
    assert all(lhs[o] == tgt[o] for o in range(NO)), 'orbit identity'
    log(f'orbit identity exact on all {NO} orbits (t = {Fraction(*Z["t"])})')
    num = np.empty(NMS, dtype=object)
    num[:] = Z['num']
    fidx = full_index()
    zfull = num[fidx].reshape((ND,) * n)
    log(f'exact z: {NMS} multisets, denominator bits {Z["den"].bit_length()}, max numerator bits '
        f'{max(abs(x).bit_length() for x in Z["num"])}')
    worst = None
    nblocks = 0
    YI = {}
    for tk in types:
        k0, k1, k2, k3 = tk
        js = [0] * k0 + [1] * k1 + [2] * k2 + [3] * k3
        # contract the scalar copies first, then spin 2, then spin 1 (smallest growth first)
        order = sorted(range(n), key=lambda i: {0: 0, 3: 0, 2: 1, 1: 2}[js[i]])
        Tn = zfull
        axes = list(range(n))                      # current position of each original copy axis
        for i in order:
            j = js[i]
            Pf = PhiE[j].reshape(ND, -1)            # (15, m_j^2)
            pos = axes.index(i)
            Tn = np.tensordot(Tn, Pf, axes=([pos], [0]))         # new axis at the end
            axes = axes[:pos] + axes[pos + 1:] + [i]
        # bring axes back to copy order, drop the scalar copies, split (r, c)
        perm = [axes.index(i) for i in range(n)]
        Tn = Tn.transpose(perm)
        keep = [i for i in range(n) if js[i] in (1, 2)]
        Tn = Tn.reshape([mj[js[i]] ** 2 for i in range(n) if js[i] in (1, 2)] or [1])
        msz = [mj[js[i]] for i in keep]
        if keep:
            Tn = Tn.reshape([x for m_ in msz for x in (m_, m_)])
            q = len(keep)
            Tn = Tn.transpose([2 * i for i in range(q)] + [2 * i + 1 for i in range(q)])
        dimx = int(np.prod(msz)) if keep else 1
        X = Tn.reshape(dimx, dimx)
        assert all(X[r, c] == X[c, r] for r in range(dimx) for c in range(r))
        for l1 in (partitions(k1, 3) if k1 else [()]):
            for l2 in (partitions(k2, 2) if k2 else [()]):
                if (3, l1) not in YI:
                    YI[(3, l1)] = young_image(3, l1)
                if (2, l2) not in YI:
                    YI[(2, l2)] = young_image(2, l2)
                W1, W2 = YI[(3, l1)], YI[(2, l2)]
                Wk = np.array([[W1[i1, c1] * W2[i2, c2] for c1 in range(W1.shape[1]) for c2 in range(W2.shape[1])]
                               for i1 in range(W1.shape[0]) for i2 in range(W2.shape[0])], dtype=object)
                Xb = Wk.T @ X @ Wk
                db = Xb.shape[0]
                rk = flint.fmpz_mat([[int(v) for v in row] for row in Wk.tolist()]).rank() if db else 0
                assert rk == db, (tk, l1, l2, rk, db)
                piv = ldl_pivots(Xb.tolist())
                ok = len(piv) == db and all(p > 0 for p in piv)
                # relative size of the smallest pivot (diagnostic only)
                rel = min(float(p) for p in piv) / max(float(abs(Xb[i, i])) for i in range(db))
                nblocks += 1
                if worst is None or rel < worst[0]:
                    worst = (rel, tk, l1, l2)
                if not ok:
                    log('NOT POSITIVE DEFINITE:', tk, l1, l2, 'pivots', [float(p) for p in piv[-3:]])
                    sys.exit(1)
        log('type', tk, 'ok')
    log(f'ALL {nblocks} blocks positive definite (exact LDL^T). Smallest relative pivot {worst[0]:.3e} at {worst[1:]}. '
        f'Hence (x(x)y(x)y)^T Z (x(x)y(x)y) = <y|P_a^(x)n(xx^T)|y>|y|^2 - t|x|^2|y|^4 with Z > 0, a = {anum}/{aden}, '
        f't = {Z["t"][0]}/{Z["t"][1]}.')
