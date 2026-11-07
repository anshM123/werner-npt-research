"""Independent exact verifier for the three-copy thickness certificate (cert/tc3_certificate.json.xz).

Written separately from the generating scripts: it rebuilds diagrams, the target, the embeddings, the relabeling
group and its orbits, the spin bases and the forced kernels from scratch.  From the certificate it only decodes the
operators (blocks R Y R^T pushed through the stored decoders Psi to diagram coefficients) and reads the integer
relation vectors, which it re-checks to be zero operators.  Nothing in the decoding needs to be trusted: positivity
and the identity are checked on the decoded operators.  Needs NumPy and python-flint.

Claim: for real x, y in (R^3)^{(x)3} and real w_0, w_1, w_2 in R^3,
   P(x, y, w) = f(x,y)|y|^2 W - (1/64)|x|^2|y|^4 W - (1/64)(98/27)|x|^2|y|^2 sum_s <w_s|rho^y_s|w_s> W/|w_s|^2
equals Avg_{sigma in S_3} [ <Psi|Z0 (x) I_w + Z1|Psi> + Tr(Lambda(x, y') G_0(y, w_0)) |w_1|^2 |w_2|^2 ] o sigma,
W = prod |w_s|^2, Psi = x (x) y (x) y (x) w (strands per copy: x, y, y', w), Z0, Z1, L >= 0.
Checks: (A) every relation used is the zero operator; (B) all orbit sums of target - terms vanish exactly;
(C) Z0, Z1, L are positive semidefinite (Gram blocks on m = 0 spin vectors; forced kernel K checked exactly
with Gamma K = 0 and Gamma + K K^T > 0 by exact LDL^T)."""
import itertools
import pickle
import sys
import time
from fractions import Fraction
import numpy as np
import flint

T0 = time.time()
ok = True


def log(*a):
    print(*a, f'[{time.time() - T0:.1f}s]', flush=True)


# ---------------------------------------------------------------- diagrams (independent enumeration)
def all_matchings(pts):
    pts = list(pts)
    if not pts:
        return [frozenset()]
    p, rest = pts[0], pts[1:]
    out = []
    for i, q in enumerate(rest):
        for m in all_matchings(rest[:i] + rest[i + 1:]):
            out.append(m | {frozenset((p, q))})
    return out


def key_of(m):
    return tuple(sorted(tuple(sorted(p)) for p in m))


def ordered(k):
    ident = key_of([(s, s + k) for s in range(k)])
    keys = sorted(set(key_of(m) for m in all_matchings(range(2 * k))))
    keys.remove(ident)
    return [ident] + keys


DG = {k: ordered(k) for k in (2, 3, 4)}
IX = {k: {m: i for i, m in enumerate(DG[k])} for k in (2, 3, 4)}
assert [len(DG[k]) for k in (2, 3, 4)] == [3, 15, 105]


def dmat(m, k):
    M = np.zeros((3 ** k, 3 ** k), dtype=np.int64)
    for ii in itertools.product(range(3), repeat=2 * k):
        if all(ii[p] == ii[q] for p, q in m):
            M[int(''.join(map(str, ii[:k])), 3), int(''.join(map(str, ii[k:])), 3)] = 1
    return M


MAT = {k: [dmat(m, k) for m in DG[k]] for k in (2, 3, 4)}


def lift(m_small, strands, k_small):
    """matching on k_small strands placed on 4-strand strands `strands`, identity on the other strands"""
    tr = {}
    for i, s in enumerate(strands):
        tr[i] = s
        tr[k_small + i] = 4 + s
    pairs = [(tr[p], tr[q]) for p, q in m_small] + [(s, 4 + s) for s in range(4) if s not in strands]
    return key_of(pairs)


# ---------------------------------------------------------------- read and decode the certificate
import json
import lzma
CERT = json.load(lzma.open(sys.argv[1] if len(sys.argv) > 1 else 'cert/tc3_certificate.json.xz', 'rt', encoding='utf-8'))
for k in (2, 3, 4):
    theirs = [tuple(tuple(p) for p in m) for m in CERT['diagrams'][str(k)]]
    if theirs != DG[k]:
        ok = False
        log('diagram order differs for k =', k)


def fr(pq):
    return Fraction(pq[0], pq[1])


def q_of(sx):
    p, q_ = sx.split('/')
    return flint.fmpq(int(p), int(q_))


def decode(term):
    """diagram coefficients {(a0, a1, a2): Fraction} of the operator of a term, from the blocks R Y R^T"""
    ks = CERT['terms'][term]
    algs = [CERT['algebras'][str(k)] for k in ks]
    cidx = [{tuple(c): i for i, c in enumerate(al['coords'])} for al in algs]
    ncoord = [len(al['coords']) for al in algs]
    z = {}
    for B in CERT['blocks']:
        if B['term'] != term:
            continue
        beta = B['beta']
        R = flint.fmpq_mat([[flint.fmpq(int(x)) for x in row] for row in B['R']])
        Y = flint.fmpq_mat([[q_of(x) for x in row] for row in B['Y']])
        G = R * Y * R.transpose()
        ms = [len(al['V'][str(j)][0]) for al, j in zip(algs, beta)]
        rows_ = list(itertools.product(*[range(m) for m in ms]))
        for Ri, rr in enumerate(rows_):
            for Ci, cc in enumerate(rows_):
                v = G[Ri, Ci]
                if v != 0:
                    key = tuple(cidx[s][(beta[s], rr[s], cc[s])] for s in range(3))
                    z[key] = z.get(key, flint.fmpq(0)) + v
    # c = (x)_s Psi_s z, mode by mode with flint matrices
    Psis = [flint.fmpq_mat([[q_of(x) for x in row] for row in al['Psi']]) for al in algs]
    cur = z
    for s in range(3):
        # group by the other two coordinates
        others = sorted(set((k_[:s] + k_[s + 1:]) for k_ in cur))
        oidx = {o: i for i, o in enumerate(others)}
        M = flint.fmpq_mat(ncoord[s], len(others))
        for k_, v in cur.items():
            M[k_[s], oidx[k_[:s] + k_[s + 1:]]] = v
        P = Psis[s] * M
        nxt = {}
        for a in range(P.nrows()):
            for o, i in oidx.items():
                v = P[a, i]
                if v != 0:
                    nxt[o[:s] + (a,) + o[s:]] = v
        cur = nxt
    return {k_: Fraction(int(v.p), int(v.q)) for k_, v in cur.items()}


CZ = {'Z0': decode('Z0'), 'Z1': decode('Z1'), 'L1': decode('L')}
rel4 = CERT['relations']
log('decoded: Z0', len(CZ['Z0']), 'Z1', len(CZ['Z1']), 'L', len(CZ['L1']), 'coefficients')

# (A) relations are zero operators
for r in rel4:
    if np.abs(sum(c * MAT[4][a] for a, c in enumerate(r))).max() != 0:
        ok = False
log('(A) relations are zero operators:', ok)

# ---------------------------------------------------------------- all terms as coefficients over 4-strand tuples
terms = {}


def add(tup, c):
    terms[tup] = terms.get(tup, Fraction(0)) + c


I4 = key_of([(0, 4), (1, 5), (2, 6), (3, 7)])
Fxy = key_of([(0, 5), (1, 4), (2, 6), (3, 7)])
Oxy = key_of([(0, 1), (4, 5), (2, 6), (3, 7)])
Oyw = key_of([(1, 3), (5, 7), (0, 4), (2, 6)])
t = Fraction(1, 64)
kap = Fraction(98, 27)
cdict = {I4: Fraction(3, 4), Fxy: Fraction(-1, 4), Oxy: Fraction(-1, 4)}
for combo in itertools.product(cdict.items(), repeat=3):
    add(tuple(k for k, _ in combo), combo[0][1] * combo[1][1] * combo[2][1])
add((I4, I4, I4), -t)
for s in range(3):
    add(tuple(Oyw if r == s else I4 for r in range(3)), -t * kap)
# subtract the certificate terms
for (a0, a1, a2), c in CZ['Z0'].items():
    add(tuple(lift(DG[3][a], [0, 1, 2], 3) for a in (a0, a1, a2)), -c)
for (a0, a1, a2), c in CZ['Z1'].items():
    add((DG[4][a0], lift(DG[3][a1], [0, 1, 2], 3), lift(DG[3][a2], [0, 1, 2], 3)), -c)


def mult_images(m3):
    """G term at the special copy for a 3-strand diagram on (idx, x, y'): list of (4-strand matching, coefficient)"""
    first = lift(m3, [1, 0, 2], 3)                       # |w|^2 Tr(Lambda rho^y): idx -> y strand, identity on w
    pairs = [tuple(p) for p in m3]
    if (0, 3) in pairs:
        fac = 3
        rest = [p for p in pairs if p != (0, 3)]
    else:
        fac = 1
        a_ = [q for p in pairs for q in p if 0 in p and q != 0][0]
        b_ = [q for p in pairs for q in p if 3 in p and q != 3][0]
        rest = [p for p in pairs if 0 not in p and 3 not in p] + [(a_, b_)]
    mp = {1: 0, 2: 2, 4: 4, 5: 6}
    second = key_of([(mp[p], mp[q]) for p, q in rest] + [(1, 3), (5, 7)])   # <w|rho^y|w> Tr Lambda
    return [(first, 1), (second, -fac)]


for (a0, a1, a2), c in CZ['L1'].items():
    for mm, e in mult_images(DG[3][a0]):
        add((mm, lift(DG[2][a1], [0, 2], 2), lift(DG[2][a2], [0, 2], 2)), -c * e)
for r, a1, a2, sx in CERT['relation_terms']:
    c = Fraction(*map(int, sx.split('/')))
    for a, rc in enumerate(rel4[r]):
        if rc:
            add((DG[4][a], lift(DG[3][a1], [0, 1, 2], 3), lift(DG[3][a2], [0, 1, 2], 3)), -c * rc)
log('collected', len(terms), 'diagram tuples')

# ---------------------------------------------------------------- relabeling group and orbit sums
def act(m, perm):
    return key_of([(perm[p], perm[q]) for p, q in m])


glob = []
for sx in (False, True):
    for py in itertools.permutations([1, 2, 5, 6]):
        perm = list(range(8))
        if sx:
            perm[0], perm[4] = 4, 0
        for a_, b_ in zip([1, 2, 5, 6], py):
            perm[a_] = b_
        glob.append(perm)
wperm = [0, 1, 2, 7, 4, 5, 6, 3]
cache = {}


def images(m):
    if m not in cache:
        cache[m] = [(act(m, g), act(act(m, g), wperm)) for g in glob]
    return cache[m]


orbits = {}
for tup, c in terms.items():
    if c == 0:
        continue
    ims = [images(m) for m in tup]
    best = None
    for gi in range(len(glob)):
        opts = [ims[s][gi] for s in range(3)]
        for ws in itertools.product((0, 1), repeat=3):
            cand = [opts[s][ws[s]] for s in range(3)]
            for p in itertools.permutations(range(3)):
                tt = (cand[p[0]], cand[p[1]], cand[p[2]])
                if best is None or tt < best:
                    best = tt
    orbits[best] = orbits.get(best, Fraction(0)) + c
nz = [k for k, v in orbits.items() if v != 0]
log('(B) orbits touched', len(orbits), '; non-zero orbit sums', len(nz))
ok &= (len(nz) == 0)

# ---------------------------------------------------------------- PSD of Z0, Z1, L (independent spin bases)
def nullspace_q(M):
    """exact null space (columns) of a rational/integer matrix via flint rref"""
    Mq = flint.fmpq_mat([[flint.fmpq(int(x)) if not isinstance(x, (Fraction, flint.fmpq)) else (flint.fmpq(x.numerator, x.denominator) if isinstance(x, Fraction) else x) for x in row] for row in M])
    R, rk = Mq.rref()
    ncol = Mq.ncols()
    piv = []
    r = 0
    for c in range(ncol):
        if r < rk and R[r, c] != 0:
            piv.append(c)
            r += 1
    free = [c for c in range(ncol) if c not in piv]
    out = []
    for f in free:
        v = [flint.fmpq(0)] * ncol
        v[f] = flint.fmpq(1)
        for i, pc in enumerate(piv):
            v[pc] = -R[i, f]
        out.append(v)
    return out


def so3(k):
    eps = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}
    gens = []
    for q in range(3):
        L = np.zeros((3, 3), dtype=np.int64)
        for (i, j, l), s in eps.items():
            if i == q:
                L[j, l] = -s
        tot = np.zeros((3 ** k, 3 ** k), dtype=np.int64)
        for s in range(k):
            f = [np.eye(3, dtype=np.int64)] * k
            f[s] = L
            mm = f[0]
            for x in f[1:]:
                mm = np.kron(mm, x)
            tot += mm
        gens.append(tot)
    return gens


def spin_bases(k):
    L = so3(k)
    Cas = L[0] @ L[0] + L[1] @ L[1] + L[2] @ L[2]
    out = {}
    for j in range(k + 1):
        stacked = np.vstack([L[2], Cas + j * (j + 1) * np.eye(3 ** k, dtype=np.int64)])
        ns = nullspace_q(stacked)
        if ns:
            out[j] = flint.fmpq_mat([[v[i] for v in ns] for i in range(3 ** k)])
    return out


def span_basis(vecs):
    M = flint.fmpq_mat([[flint.fmpq(int(x)) for x in v] for v in vecs]).transpose()
    R, rk = M.transpose().rref()
    return [[R[i, c] for c in range(M.nrows())] for i in range(rk)], rk


SB = {k: spin_bases(k) for k in (2, 3, 4)}
log('spin multiplicities', {k: {j: SB[k][j].ncols() for j in SB[k]} for k in SB})
# forced kernel subspaces from exact rational samples
es = [(1, 2, 3), (1, -1, 2), (2, 1, -1), (0, 1, 1), (1, 0, -2), (3, -2, 1), (1, 1, 1), (2, -3, 1), (1, 4, -1), (0, 0, 1),
      (1, 0, 0), (0, 1, 0), (2, 2, -1), (1, -2, -3)]
S3v = [np.kron(np.kron(e, e), e) for e in map(np.array, es)]
S4v = []
for e in map(np.array, es):
    for w in (np.array([e[1], -e[0], 0]), np.array([e[2], 0, -e[0]]), np.array([0, e[2], -e[1]])):
        if np.any(w):
            S4v.append(np.kron(np.kron(np.kron(e, e), e), w))
S2v = [np.kron(e, e) for e in map(np.array, es)]


def kernel_coords(k, j, samples):
    """coordinates u with V_j u in span(samples), as columns"""
    V = SB[k][j]
    basis, rk = span_basis(samples)
    Bm = flint.fmpq_mat(basis).transpose() if rk else None
    # solve V u - B v = 0
    cols = V.ncols()
    big = [[V[i, c] for c in range(cols)] + ([-Bm[i, c] for c in range(Bm.ncols())] if Bm is not None else []) for i in range(V.nrows())]
    ns = nullspace_q(big)
    us = [v[:cols] for v in ns]
    if not us:
        return flint.fmpq_mat(cols, 0)
    Mu = flint.fmpq_mat(us)          # rows
    R, rk2 = Mu.rref()
    return flint.fmpq_mat([[R[i, c] for c in range(cols)] for i in range(rk2)]).transpose() if rk2 else flint.fmpq_mat(cols, 0)


KER = {3: {j: kernel_coords(3, j, S3v) for j in SB[3]}, 4: {j: kernel_coords(4, j, S4v) for j in SB[4]},
       2: {j: kernel_coords(2, j, S2v) for j in SB[2]}}
log('forced kernels', {k: {j: KER[k][j].ncols() for j in KER[k]} for k in KER})


def gram_tables(k):
    """G[j][a] = V_j^T D_a V_j (fmpq_mat)"""
    out = {}
    for j, V in SB[k].items():
        Vt = V.transpose()
        out[j] = [Vt * flint.fmpq_mat([[int(x) for x in row] for row in MAT[k][a]]) * V for a in range(len(DG[k]))]
    return out


GT = {k: gram_tables(k) for k in (2, 3, 4)}


def kron_q(A, B):
    m, n_ = A.nrows(), A.ncols()
    p, q_ = B.nrows(), B.ncols()
    out = flint.fmpq_mat(m * p, n_ * q_)
    for i in range(m):
        for j in range(n_):
            a = A[i, j]
            if a != 0:
                for r in range(p):
                    for c in range(q_):
                        out[i * p + r, j * q_ + c] = a * B[r, c]
    return out


def ldl_pos(M):
    m = M.nrows()
    A_ = [[M[i, j] for j in range(m)] for i in range(m)]
    for k in range(m):
        if A_[k][k] <= 0:
            return False
        for i in range(k + 1, m):
            if A_[i][k] != 0:
                f = A_[i][k] / A_[k][k]
                for j in range(k + 1, m):
                    A_[i][j] -= f * A_[k][j]
    return True


def check_psd(name, coeffs, ks):
    coeffs = {k_: (v.numerator, v.denominator) for k_, v in coeffs.items()}
    good = True
    # group coefficients by (a0, a1) for a mode-wise contraction
    for beta in itertools.product(*[sorted(SB[k]) for k in ks]):
        Gam = None
        acc = {}
        for (a0, a1, a2), pq in coeffs.items():
            acc.setdefault((a0, a1), []).append((a2, fr(pq)))
        # sum_{a0,a1} G0[a0] (x) G1[a1] (x) (sum_{a2} c G2[a2])
        for (a0, a1), lst in acc.items():
            inner = None
            for a2, c in lst:
                term = GT[ks[2]][beta[2]][a2] * flint.fmpq(c.numerator, c.denominator)
                inner = term if inner is None else inner + term
            blk = kron_q(kron_q(GT[ks[0]][beta[0]][a0], GT[ks[1]][beta[1]][a1]), inner)
            Gam = blk if Gam is None else Gam + blk
        K = kron_q(kron_q(KER[ks[0]][beta[0]], KER[ks[1]][beta[1]]), KER[ks[2]][beta[2]]) if all(KER[ks[s]][beta[s]].ncols() for s in range(3)) else None
        if K is not None and K.ncols():
            if not (Gam * K).entries() or any(x != 0 for x in (Gam * K).entries()):
                good = False
                log('  forced kernel not annihilated in', name, beta)
                continue
            M = Gam + K * K.transpose()
        else:
            M = Gam
        if not (M == M.transpose()) or not ldl_pos(M):
            good = False
            log('  block not positive:', name, beta)
    log(f'(C) {name} positive semidefinite:', good)
    return good


ok &= check_psd('Z0', CZ['Z0'], (3, 3, 3))
ok &= check_psd('L', CZ['L1'], (3, 2, 2))
ok &= check_psd('Z1', CZ['Z1'], (4, 3, 3))
print('PASS' if ok else 'FAIL')
