"""Independent check of a level-2 positivity certificate for P_a^{(x)n},  P_a(X) = a Tr(X) I - (X + X^T)/4 on 3 x 3 matrices.

Certificate (JSON, written by cert_l2.py export): n, a, t, the 15 Brauer diagrams on three strands (x, y, y') given as
perfect matchings of the points 0..5 (top x, y, y' = 0, 1, 2; bottom = 3, 4, 5), a common denominator and one integer
numerator per multiset of n diagrams.  Z = sum over ordered tuples tau of z_{sorted tau} D_{tau_1} (x) ... (x) D_{tau_n}.

Claim: (x (x) y (x) y)^T Z (x (x) y (x) y) = <y|P_a^{(x)n}(x x^T)|y> |y|^2 - t |x|^2 |y|^4 for all real x, y, and Z >= 0.
Hence <y|P_a^{(x)n}(x x^T)|y> >= t |x|^2 |y|^2 for real x, y (and complex positivity follows, see the README).

Checked here without using cert_l2.py or sym_real_l2.py:
 1. the diagrams are the 15 distinct perfect matchings; the diagram matrices on (R^3)^{(x)3} are linearly independent;
 2. the identity, as equality of the sums over the orbits of the relabelling group S_2 (x points {0,3}) x S_4 (y points
    {1,2,4,5}) acting on all copies at once: two O(3)^n-invariant tensors give the same polynomial in (x, y) iff their
    S_2 x S_4 symmetrisations agree, and product diagrams form a basis.  The coefficients of C_a (x) I_{y'} in the diagram
    basis are recomputed from the 27 x 27 matrix;
 3. Z >= 0: Z commutes with O(3) at every copy, so it is the direct sum of its blocks on the m = 0 vectors of each spin
    (integer bases recomputed here); copies are grouped by spin, and each block is reduced by S_{k1} x S_{k2} using
    Jucys-Murphy eigenvectors (exact null spaces); every reduced block has an exact LDL^T with positive pivots;
 4. sanity (floating point): the identity at random x, y of tensor rank 2.
usage: python verify_l2.py cert.json"""
import itertools
import json
import math
import sys
import time
from fractions import Fraction
import numpy as np
import flint

T0 = time.time()


def log(*a):
    print(*a, f'[{time.time() - T0:.1f}s]', flush=True)


if sys.argv[1].endswith('.xz'):
    import lzma
    C = json.load(lzma.open(sys.argv[1], 'rt'))
else:
    C = json.load(open(sys.argv[1]))
n = C['n']
a = Fraction(*C['a'])
t = Fraction(*C['t'])
den = int(C['den'])
DG = [tuple(sorted(tuple(sorted(p)) for p in mm)) for mm in C['diagrams']]
ND = len(DG)

# ---------------------------------------------------------------- 1. diagrams
allm = set()
for perm in itertools.permutations(range(6)):
    allm.add(tuple(sorted(tuple(sorted((perm[2 * i], perm[2 * i + 1]))) for i in range(3))))
assert len(allm) == 15 and set(DG) == allm and len(set(DG)) == 15, 'diagram list'
POS = {d: i for i, d in enumerate(DG)}


def dmatrix(mm):
    """27 x 27 integer matrix; rows = top indices (points 0,1,2), columns = bottom indices (points 3,4,5)"""
    M = np.zeros((27, 27), dtype=np.int64)
    for top in itertools.product(range(3), repeat=3):
        for bot in itertools.product(range(3), repeat=3):
            v = top + bot
            if all(v[p] == v[q] for p, q in mm):
                M[top[0] * 9 + top[1] * 3 + top[2], bot[0] * 9 + bot[1] * 3 + bot[2]] = 1
    return M


DMAT = [dmatrix(mm) for mm in DG]
Amat = flint.fmpz_mat([[int(DMAT[d].reshape(-1)[e]) for d in range(ND)] for e in range(729)])
assert Amat.rank() == 15, 'diagram matrices independent'

# coefficients of C_a (x) I_{y'} (legs x, y, y'), C_a = sum_ij E_ij (x) P_a(E_ij)
Cmat = [[Fraction(0)] * 27 for _ in range(27)]
for i in range(3):
    for j in range(3):
        E = [[Fraction(int(r == i and c == j)) for c in range(3)] for r in range(3)]
        P = [[(a * (1 if i == j else 0) * (1 if r == c else 0)) - Fraction(1, 4) * (E[r][c] + E[c][r]) for c in range(3)]
             for r in range(3)]
        for r in range(3):
            for c in range(3):
                for k in range(3):
                    Cmat[i * 9 + r * 3 + k][j * 9 + c * 3 + k] += P[r][c]
rhs = flint.fmpq_mat(729, 1, [flint.fmpq(Cmat[e // 27][e % 27].numerator, Cmat[e // 27][e % 27].denominator) for e in range(729)])
Aq = flint.fmpq_mat(Amat)
AtA = Aq.transpose() * Aq
cvec = AtA.solve(Aq.transpose() * rhs)
assert Aq * cvec == rhs, 'C_a (x) I is in the diagram span'
cv = [Fraction(int(cvec[d, 0].p), int(cvec[d, 0].q)) for d in range(ND)]
ident = POS[((0, 3), (1, 4), (2, 5))]
log(f'n = {n}, a = {a}, t = {t}; per-copy target coefficients:',
    {DG[d]: str(cv[d]) for d in range(ND) if cv[d]})

# ---------------------------------------------------------------- 2. orbit identity
coef = {}
for ms, num in C['coefficients']:
    ms = tuple(ms)
    assert len(ms) == n and list(ms) == sorted(ms) and ms not in coef
    coef[ms] = int(num)
allms = list(itertools.combinations_with_replacement(range(ND), n))
assert len(coef) == len(allms) and all(ms in coef for ms in allms), 'one coefficient per multiset'
group = []
for perm in itertools.permutations(range(6)):
    if {perm[0], perm[3]} == {0, 3}:
        group.append([POS[tuple(sorted(tuple(sorted((perm[p], perm[q]))) for p, q in mm))] for mm in DG])
assert len(group) == 48
orbit = {}
for ms in allms:
    orbit[ms] = min(tuple(sorted(g[d] for d in ms)) for g in group)


def ordcount(ms):
    out = math.factorial(n)
    for d in set(ms):
        out //= math.factorial(ms.count(d))
    return out


lhs, tgt = {}, {}
for ms in allms:
    o = orbit[ms]
    oc = ordcount(ms)
    lhs[o] = lhs.get(o, 0) + oc * Fraction(coef[ms], den)
    p = Fraction(1)
    for d in ms:
        p *= cv[d]
    tg = oc * p - (t if all(d == ident for d in ms) else 0)
    tgt[o] = tgt.get(o, 0) + tg
bad = [o for o in lhs if lhs[o] != tgt[o]]
log(f'orbit identity: {len(lhs)} orbits, {len(bad)} mismatches')
assert not bad

# ---------------------------------------------------------------- 3. Z >= 0
# real rotation generators on R^3 and their total action on (R^3)^{(x)3}
Lg = []
for q in range(3):
    L = np.zeros((3, 3), dtype=np.int64)
    a1, a2 = [(1, 2), (2, 0), (0, 1)][q]
    L[a2, a1], L[a1, a2] = 1, -1
    tot = np.zeros((27, 27), dtype=np.int64)
    for s in range(3):
        f = [np.eye(3, dtype=np.int64)] * 3
        f[s] = L
        tot += np.kron(np.kron(f[0], f[1]), f[2])
    Lg.append(tot)
Cas = sum(L @ L for L in Lg)


def int_nullspace(M):
    X, k = flint.fmpz_mat([[int(v) for v in row] for row in M]).nullspace()
    cols = []
    for c in range(k):
        col = [int(X[r, c]) for r in range(X.nrows())]
        g = 0
        for v in col:
            g = math.gcd(g, abs(v))
        cols.append([v // g for v in col])
    return np.array(cols, dtype=object).T if cols else np.zeros((M.shape[1], 0), dtype=object)


Vsp = {}
for j in range(4):
    Vsp[j] = int_nullspace(np.vstack([Lg[2], Cas + j * (j + 1) * np.eye(27, dtype=np.int64)]))
msp = {j: Vsp[j].shape[1] for j in range(4)}
assert [msp[j] for j in range(4)] == [1, 3, 2, 1] and sum(msp[j] * (2 * j + 1) for j in range(4)) == 27
PH = {j: [Vsp[j].T @ DMAT[d].astype(object) @ Vsp[j] for d in range(ND)] for j in range(4)}


def jm_block_basis(m, lam):
    """exact integer basis of the joint eigenspace of the Jucys-Murphy elements on (R^m)^{(x)k} for the row-reading
    tableau of shape lam (a copy of the GL_m irreducible V_lam)"""
    k = sum(lam)
    if k == 0:
        return np.ones((1, 1), dtype=object)
    content = []
    for r, l in enumerate(lam):
        content += [c - r for c in range(l)]
    Dm = m ** k
    idx = list(itertools.product(range(m), repeat=k))
    pos = {v: i for i, v in enumerate(idx)}
    rows = []
    for i in range(1, k):
        J = np.zeros((Dm, Dm), dtype=np.int64)
        for l in range(i):
            for v in idx:
                w = list(v)
                w[l], w[i] = w[i], w[l]
                J[pos[tuple(w)], pos[v]] += 1
        rows.append(J - content[i] * np.eye(Dm, dtype=np.int64))
    if not rows:                                    # k = 1: no conditions
        return np.eye(Dm, dtype=np.int64).astype(object)
    return int_nullspace(np.vstack(rows))


def parts(k, maxrows):
    out = []

    def rec(rem, mx, cur):
        if rem == 0:
            out.append(tuple(cur))
            return
        if len(cur) == maxrows:
            return
        for p in range(min(rem, mx), 0, -1):
            rec(rem - p, p, cur + [p])
    rec(k, k, [])
    return out


def hooks(lam):
    conj = [sum(1 for l in lam if l > c) for c in range(lam[0])] if lam else []
    return [lam[r] - c + conj[c] - r - 1 for r in range(len(lam)) for c in range(lam[r])]


def n_syt(lam):
    return math.factorial(sum(lam)) // math.prod(hooks(lam)) if lam else 1


def dim_gl(m, lam):
    num = math.prod(m + c - r for r in range(len(lam)) for c in range(lam[r]))
    return num // math.prod(hooks(lam)) if lam else 1


def ldl_ok(A):
    d = len(A)
    M = [[Fraction(A[i][j]) for j in range(d)] for i in range(d)]
    for k in range(d):
        p = M[k][k]
        if p <= 0:
            return False
        for i in range(k + 1, d):
            f = M[i][k] / p
            if f:
                for j in range(k + 1, d):
                    M[i][j] -= f * M[k][j]
    return True


num_full = np.empty((ND,) * n, dtype=object)
for tau in itertools.product(range(ND), repeat=n):
    num_full[tau] = coef[tuple(sorted(tau))]
JB = {}
nblocks = 0
for tk in [x for x in itertools.product(range(n + 1), repeat=4) if sum(x) == n]:
    js = sum(([j] * tk[j] for j in range(4)), [])
    T = num_full
    for i in range(n):                              # contract copy 0, 1, ..., n-1 in turn (axis 0 each time)
        P = np.array([PH[js[i]][d].reshape(-1) for d in range(ND)], dtype=object)       # (15, m^2)
        T = np.tensordot(T, P, axes=([0], [0]))                                         # new axis at the end
    keep = [i for i in range(n) if msp[js[i]] > 1]
    msz = [msp[js[i]] for i in keep]
    T = T.reshape([x for m_ in msz for x in (m_, m_)] or [1, 1])
    q = len(keep)
    if q:
        T = T.transpose([2 * i for i in range(q)] + [2 * i + 1 for i in range(q)])
    dimx = int(np.prod(msz)) if q else 1
    X = T.reshape(dimx, dimx)
    k1, k2 = tk[1], tk[2]
    tot_dim = 0
    for l1 in (parts(k1, 3) if k1 else [()]):
        for l2 in (parts(k2, 2) if k2 else [()]):
            for key, m_ in (((3, l1), 3), ((2, l2), 2)):
                if key not in JB:
                    JB[key] = jm_block_basis(m_, key[1])
                    assert JB[key].shape[1] == dim_gl(m_, key[1]), key
            W1, W2 = JB[(3, l1)], JB[(2, l2)]
            tot_dim += n_syt(l1) * n_syt(l2) * W1.shape[1] * W2.shape[1]
            Wk = np.empty((W1.shape[0] * W2.shape[0], W1.shape[1] * W2.shape[1]), dtype=object)
            for i1 in range(W1.shape[0]):
                for i2 in range(W2.shape[0]):
                    for c1 in range(W1.shape[1]):
                        for c2 in range(W2.shape[1]):
                            Wk[i1 * W2.shape[0] + i2, c1 * W2.shape[1] + c2] = W1[i1, c1] * W2[i2, c2]
            B = Wk.T @ X @ Wk
            assert all(B[r, c] == B[c, r] for r in range(len(B)) for c in range(r))
            assert ldl_ok(B.tolist()), ('block not positive definite', tk, l1, l2)
            nblocks += 1
    assert tot_dim == dimx, (tk, tot_dim, dimx)                   # the isotypic pieces fill the whole block
log(f'Z >= 0: {nblocks} reduced blocks, all positive definite (exact LDL^T)')

# ---------------------------------------------------------------- 4. sanity at random points (floating point)
rng = np.random.default_rng(7)
zf = np.vectorize(lambda v: float(Fraction(v, den)))(num_full).astype(float)
worst = 0.0
for trial in range(3):
    xs = [[rng.standard_normal(3) for _ in range(n)] for _ in range(2)]
    ys = [[rng.standard_normal(3) for _ in range(n)] for _ in range(2)]
    DMf = [M.astype(float) for M in DMAT]
    lhs_v = 0.0
    for s, sb, u1, u2, v1, v2 in itertools.product(range(2), repeat=6):
        Tz = zf
        for i in range(n):
            ket = np.kron(np.kron(xs[s][i], ys[u1][i]), ys[u2][i])
            bra = np.kron(np.kron(xs[sb][i], ys[v1][i]), ys[v2][i])
            vec = np.array([bra @ DMf[d] @ ket for d in range(ND)])
            Tz = np.tensordot(Tz, vec, axes=([0], [0]))
        lhs_v += float(Tz)

    def Pa(X):
        return float(a) * np.trace(X) * np.eye(3) - (X + X.T) / 4
    F = 0.0
    for s, sb, u, ub in itertools.product(range(2), repeat=4):
        F += np.prod([ys[ub][i] @ Pa(np.outer(xs[s][i], xs[sb][i])) @ ys[u][i] for i in range(n)])
    ny = sum(np.prod([ys[u][i] @ ys[ub][i] for i in range(n)]) for u, ub in itertools.product(range(2), repeat=2))
    nx = sum(np.prod([xs[s][i] @ xs[sb][i] for i in range(n)]) for s, sb in itertools.product(range(2), repeat=2))
    rhs_v = F * ny - float(t) * nx * ny * ny
    worst = max(worst, abs(lhs_v - rhs_v) / max(1.0, abs(rhs_v), F * ny))
log(f'random rank-2 points: max relative mismatch {worst:.2e}')
assert worst < 1e-9
log(f'PASS: P_a^(x)n is positive for n = {n}, a = {a}; <y|P(xx^T)|y> >= {t} |x|^2 |y|^2 on real vectors')
