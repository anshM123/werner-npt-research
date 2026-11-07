"""Independent exact check of a level-2 symmetric-extension certificate for C_a^{(x)n}.

Statement.  P_a(X) = a Tr(X) I - (X + X^T)/4 on M_3(C), Choi matrix C_a = sum_ij |i><j| (x) P_a(|i><j|) on
A (x) B = C^3 (x) C^3.  For n copies, W = C_a^{(x)n} on A_1..A_n B_1..B_n.  The certificate claims
    <x (x) y| W |x (x) y>  >=  t |x|^2 |y|^2      for all x, y in C^{3^n},
i.e. <y|P_a^{(x)n}(|conj x><conj x|)|y> >= t |x|^2 |y|^2; for t >= 0 this says P_a^{(x)n} is positive.

Certificate.  Operators QA, QB on A B B' (B' a second copy of B), as rational combinations of products of per-copy
Brauer diagrams (perfect matchings of the six points: rows of strands A_s, B_s, B'_s = points 0, 1, 2, columns = 3, 4, 5).
With S = SWAP(B, B') and P = (I + S)/2:   QA' = P QA P,   Z = P [W (x) I_B' - t I - Gamma_A(QA') - Gamma_B(QB)] P.
If QA', QB, Z >= 0 then for psi = x (x) y (x) y (P psi = psi):
    <xy|W|xy> |y|^2 - t |x|^2 |y|^4 = <psi|Z|psi> + <conj(x) y y|QA'|conj(x) y y> + <x conj(y) y|QB|x conj(y) y> >= 0.

This script shares no code with the scripts that produced the certificates.  It rebuilds every ingredient from the
matchings in the file and from the definition of P_a, computes QA', QB and Z exactly as integer coefficient tensors,
and decides positive semidefiniteness exactly, block by block:
  every product of diagrams commutes with SO(3) acting on the three strands of each copy, so X >= 0 iff each
  compression G_beta = B_beta^T X B_beta >= 0, where the columns of B_beta = (x)_s B_{j_s} span the m = 0 vectors of the
  spin-(j_1..j_n) isotypic component (j in {0,1,2,3}, multiplicities 1, 3, 2, 1 per copy).
Two independent exact tests are used:
  (L) fraction-free symmetric elimination with diagonal pivoting (exact integers), for blocks up to a size limit;
  (C) congruence test: the forced kernel vector of the block (the real product minimisers, if any) is checked to be an
      exact null vector and one coordinate on its support is deleted; for the rest, an integer matrix M is built from
      floating eigenvectors and H = M^T G M is computed exactly; strict (scaled) diagonal dominance of H with positive
      diagonal proves H > 0, hence M invertible and G > 0.
usage: python verify_independent.py CERT [ldl_max] [t=p/q]     (CERT may be a .txt or a .txt.xz file; t=p/q replaces the
constant of the certificate, for negative controls)"""
import sys, time, itertools, math, lzma
from fractions import Fraction
import numpy as np
import sympy as sp
import flint

T0 = time.time()


def log(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- reading
def read_cert(path):
    f = {}
    opener = lzma.open if path.endswith('.xz') else open
    with opener(path, 'rt') as fh:
        for line in fh:
            p = line.split()
            if p:
                f[p[0]] = p[1:]
    n = int(f['n'][0])
    a = Fraction(int(f['a'][0]), int(f['a'][1]))
    t = Fraction(int(f['t'][0]), int(f['t'][1]))
    L = int(f['L'][0])
    match = [tuple(int(v) for v in f['D%d' % k]) for k in range(15)]
    if 'QA' in f:
        return n, a, t, L, match, 'full', [int(v) for v in f['QA']], [int(v) for v in f['QB']]
    return n, a, t, L, match, 'orbit', [int(v) for v in f['QAO']], [int(v) for v in f['QBO']]


# ---------------------------------------------------------------- one copy
def diagram_matrix(m):
    assert sorted(m) == list(range(6)), m
    pairs = [(m[0], m[1]), (m[2], m[3]), (m[4], m[5])]
    D = np.zeros((27, 27), dtype=np.int64)
    for idx in itertools.product(range(3), repeat=6):
        if all(idx[p] == idx[q] for p, q in pairs):
            D[idx[0] * 9 + idx[1] * 3 + idx[2], idx[3] * 9 + idx[4] * 3 + idx[5]] = 1
    return D


def pt(M, strand):            # partial transpose of a 27x27 matrix on strand 0 (A), 1 (B) or 2 (B')
    T = M.reshape(3, 3, 3, 3, 3, 3)
    ax = list(range(6))
    ax[strand], ax[strand + 3] = ax[strand + 3], ax[strand]
    return T.transpose(ax).reshape(27, 27)


def swap_BBp():
    S = np.zeros((27, 27), dtype=np.int64)
    for i, j, k in itertools.product(range(3), repeat=3):
        S[i * 9 + k * 3 + j, i * 9 + j * 3 + k] = 1
    return S


def choi_times_identity(a):   # (C_a (x) I_B') on strands (A, B, B'), built from the definition of P_a
    C = [[Fraction(0)] * 27 for _ in range(27)]
    for i, j in itertools.product(range(3), repeat=2):
        for k, l in itertools.product(range(3), repeat=2):
            # <k| P_a(|i><j|) |l> = a d_ij d_kl - (d_ki d_lj + d_kj d_li)/4
            v = (a if (i == j and k == l) else Fraction(0)) - Fraction(int(k == i and l == j) + int(k == j and l == i), 4)
            for m in range(3):
                C[i * 9 + k * 3 + m][j * 9 + l * 3 + m] = v
    return C


def so3_generators():
    Lx = np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]])
    Ly = np.array([[0, 0, 1], [0, 0, 0], [-1, 0, 0]])
    Lz = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
    I3 = np.eye(3, dtype=np.int64)
    out = []
    for Lk in (Lx, Ly, Lz):
        out.append(np.kron(np.kron(Lk, I3), I3) + np.kron(np.kron(I3, Lk), I3) + np.kron(np.kron(I3, I3), Lk))
    return out


def int_nullspace(M):
    """Exact integer basis (columns) of the null space of an integer matrix."""
    ns = sp.Matrix(M.tolist()).nullspace()
    cols = []
    for v in ns:
        den = math.lcm(*[int(sp.fraction(x)[1]) for x in v])
        w = [int(x * den) for x in v]
        g = math.gcd(*w)
        cols.append([x // g for x in w])
    return np.array(cols, dtype=object).T if cols else np.zeros((M.shape[1], 0), dtype=object)


def setup_copy(match, a):
    D = [diagram_matrix(m) for m in match]
    flat = np.array([d.ravel() for d in D])
    assert np.linalg.matrix_rank(flat.astype(float)) == 15, 'diagrams not independent'
    assert sp.Matrix(flat.tolist()).rank() == 15
    log('15 diagram matrices built from the matchings; exactly linearly independent')

    def ident(M):
        hits = [k for k in range(15) if np.array_equal(D[k], M)]
        assert len(hits) == 1, 'not a single diagram'
        return hits[0]
    S = swap_BBp()
    perms = {
        'gA': [ident(pt(D[k], 0)) for k in range(15)],
        'gB': [ident(pt(D[k], 1)) for k in range(15)],
        'sL': [ident(S @ D[k]) for k in range(15)],
        'sR': [ident(D[k] @ S) for k in range(15)],
        'tr': [ident(D[k].T) for k in range(15)],
    }
    for name, p in perms.items():
        assert sorted(p) == list(range(15)), name
    log('partial transposes on A and on B, multiplication by SWAP(B,B\') on the left and on the right, and the full '
        'transpose map diagrams to diagrams bijectively (checked on explicit matrices)')
    # C_a (x) I in the diagram basis (exact solve, then exact check)
    C = choi_times_identity(a)
    Mx = sp.Matrix(flat.T.tolist())
    rhs = sp.Matrix([C[r][c] for r in range(27) for c in range(27)])
    sol = (Mx.T * Mx).LUsolve(Mx.T * rhs)
    assert Mx * sol == rhs, 'C_a (x) I is not in the span of the diagrams'
    c = [Fraction(int(sp.fraction(x)[0]), int(sp.fraction(x)[1])) for x in sol]
    log('C_a (x) I_B\' from the definition of P_a, in diagram coordinates:',
        {k: str(v) for k, v in enumerate(c) if v != 0})
    idx_id = ident(np.eye(27, dtype=np.int64))
    # so(3) invariance and spin bases
    G = so3_generators()
    for d in D:
        for g in G:
            assert np.array_equal(d @ g, g @ d)
    log('every diagram commutes exactly with the three so(3) generators of the copy')
    Cas = sum(g @ g for g in G)
    bases, signs, mult = {}, {}, {}
    for j in range(4):
        cols, sg = [], []
        for s in (1, -1):
            Bjs = int_nullspace(np.vstack([G[2], Cas + j * (j + 1) * np.eye(27, dtype=np.int64), S - s * np.eye(27, dtype=np.int64)]))
            for q in range(Bjs.shape[1]):
                cols.append(Bjs[:, q]); sg.append(s)
        Bj = np.array(cols, dtype=object).T
        bases[j], signs[j], mult[j] = Bj, sg, Bj.shape[1]
        assert sp.Matrix(Bj.tolist()).rank() == mult[j]
    assert int_nullspace(G[2]).shape[1] == 7 and [mult[j] for j in range(4)] == [1, 3, 2, 1], mult
    log('m = 0 spin bases (adapted to SWAP(B,B\')): multiplicities', mult, 'signs', signs)
    dcomp = {}
    for j in range(4):
        Bj = bases[j]
        dcomp[j] = np.array([(Bj.T.dot(D[k].astype(object))).dot(Bj) for k in range(15)], dtype=object)
    rows = []
    for k in range(15):
        rows.append([x for j in range(4) for x in dcomp[j][k].ravel()])
    assert sp.Matrix(rows).rank() == 15
    log('compression of the copy algebra to the four spin blocks is a bijection (exact rank 15)')
    # forced kernel: m = 0 part of Sym^3(C^3) inside each spin block
    Psym = np.zeros((27, 27), dtype=object)
    for perm in itertools.permutations(range(3)):
        for idx in itertools.product(range(3), repeat=3):
            jdx = tuple(idx[perm[q]] for q in range(3))
            Psym[jdx[0] * 9 + jdx[1] * 3 + jdx[2], idx[0] * 9 + idx[1] * 3 + idx[2]] += 1
    forced = {}
    for j in range(4):
        Bj = bases[j]
        Mj = (6 * np.eye(27, dtype=object) - Psym).dot(Bj)          # (I - Pi_sym) B_j c = 0
        ns = int_nullspace(np.array(Mj.tolist(), dtype=object))
        forced[j] = None if ns.shape[1] == 0 else [int(x) for x in ns[:, 0]]
        assert ns.shape[1] == (1 if j in (1, 3) else 0)
    log('forced (Sym^3) vectors in the spin blocks:', forced)
    return perms, c, idx_id, dcomp, signs, forced


# ---------------------------------------------------------------- tensors
def apply_perm(X, perm):      # Y[perm[a_1], ..., perm[a_n]] = X[a_1, ..., a_n]
    inv = np.argsort(np.array(perm))
    Y = X
    for ax in range(X.ndim):
        Y = np.take(Y, inv, axis=ax)
    return Y


def expand_orbits(vals, n, tr):
    keys = {}
    for tup in itertools.product(range(15), repeat=n):
        k1 = tuple(sorted(tup)); k2 = tuple(sorted(tr[x] for x in tup))
        keys[tup] = min(k1, k2)
    order = sorted(set(keys.values()))
    assert len(order) == len(vals), (len(order), len(vals))
    pos = {k: i for i, k in enumerate(order)}
    X = np.empty((15,) * n, dtype=object)
    for tup, k in keys.items():
        X[tup] = vals[pos[k]]
    return X


def compress(X, beta, dcomp):
    T = X
    for j in beta:
        T = np.tensordot(T, dcomp[j], axes=([0], [0]))
    n = len(beta)
    T = T.transpose([2 * s for s in range(n)] + [2 * s + 1 for s in range(n)])
    dim = 1
    for j in beta:
        dim *= dcomp[j].shape[1]
    return T.reshape(dim, dim)


# ---------------------------------------------------------------- exact PSD tests
def psd_bareiss(G):
    """Fraction-free symmetric elimination with diagonal pivoting.  Returns (psd, kernel_dim)."""
    m = G.shape[0]
    A = [[int(G[i, j]) for j in range(m)] for i in range(m)]
    act = list(range(m)); prev = 1
    while act:
        p = max(act, key=lambda i: A[i][i])
        d = A[p][p]
        if d < 0:
            return False, None
        if d == 0:
            for i in act:
                for j in act:
                    if A[i][j] != 0:
                        return False, None
            return True, len(act)
        act.remove(p)
        rp = A[p]
        for i in act:
            Ai = A[i]; aip = Ai[p]
            for j in act:
                q, r = divmod(d * Ai[j] - aip * rp[j], prev)
                assert r == 0
                Ai[j] = q
        prev = d
    return True, 0


def to_fmpz_mat(A):
    return flint.fmpz_mat([[int(x) for x in row] for row in A])


def psd_congruence(G, f):
    """G symmetric integer (object) matrix.  f: forced null vector (list of ints) or None.  Returns (psd, kernel_dim)."""
    m = G.shape[0]
    keep = list(range(m)); kdim = 0
    if f is not None:
        Gf = G.dot(np.array(f, dtype=object))
        if all(v == 0 for v in Gf):          # exact null vector: delete one coordinate on its support
            d = max(range(m), key=lambda i: abs(f[i]))
            keep.remove(d); kdim = 1
    Gk = G[np.ix_(keep, keep)]
    k = len(keep)
    if k == 0 or all(v == 0 for v in Gk.ravel()):
        return True, kdim + k
    diag = [int(Gk[i, i]) for i in range(k)]
    if min(diag) <= 0:
        return False, None
    e = [(x.bit_length() - 1) // 2 for x in diag]                       # 2^{2e} ~ diagonal
    F = np.empty((k, k))
    for i in range(k):
        for j in range(k):
            F[i, j] = float(Fraction(int(Gk[i, j]), 1 << (e[i] + e[j])))
    lam, V = np.linalg.eigh(F)
    if lam[0] <= 0:
        return False, None
    Mf = (V / np.sqrt(lam))                                              # columns: S V Lambda^{-1/2} up to the scaling S
    E = max(e)
    Mint = []
    for i in range(k):
        Mint.append([int(round(Mf[i, c] * 2.0 ** 52)) << (E - e[i]) for c in range(k)])   # row i scaled by 2^{-e_i} (times 2^E)
    Mm = to_fmpz_mat(Mint)
    H = Mm.transpose() * to_fmpz_mat(Gk.tolist()) * Mm
    for i in range(k):
        hii = int(H[i, i])
        if hii <= 0:
            return False, None
        off = sum(abs(int(H[i, j])) for j in range(k) if j != i)
        if not off < hii:
            return False, None
    return True, kdim


# ---------------------------------------------------------------- main
def main():
    path = sys.argv[1]
    ldl_max = int(sys.argv[2]) if len(sys.argv) > 2 else 81
    n, a, t, L, match, fmt, qa, qb = read_cert(path)
    log(f'certificate {path}: n = {n}, a = {a}, t = {t}, common denominator L ({L.bit_length()} bits), format {fmt}')
    if len(sys.argv) > 3 and sys.argv[3].startswith('t='):
        t = Fraction(sys.argv[3][2:])
        log(f'NEGATIVE CONTROL: the constant is replaced by t = {t}; this run is expected to FAIL')
    perms, c, idx_id, dcomp, signs, forced = setup_copy(match, a)
    if fmt == 'full':
        QA = np.array(qa, dtype=object).reshape((15,) * n)
        QB = np.array(qb, dtype=object).reshape((15,) * n)
    else:
        QA = expand_orbits(qa, n, perms['tr'])
        QB = expand_orbits(qb, n, perms['tr'])
    log(f'QA, QB read ({time.time() - T0:.1f} s)')
    # W (x) I - t I, scaled to integers
    cvec = np.array(c, dtype=object)
    Wf = cvec
    for _ in range(n - 1):
        Wf = np.multiply.outer(Wf, cvec)
    Wf = Wf.copy()
    Wf[(idx_id,) * n] -= t
    lam1 = math.lcm(*[x.denominator for x in Wf.ravel()])
    assert all((x * lam1).denominator == 1 for x in Wf.ravel())
    Wint = np.vectorize(lambda x: int(x * lam1), otypes=[object])(Wf)
    sL, sR, gA, gB = perms['sL'], perms['sR'], perms['gA'], perms['gB']
    QAp = QA + apply_perm(QA, sL) + apply_perm(QA, sR) + apply_perm(apply_perm(QA, sL), sR)      # = 4 L QA'
    del QA
    Bint = 4 * L * Wint - lam1 * apply_perm(QAp, gA) - 4 * lam1 * apply_perm(QB, gB)
    del Wint
    Zint = Bint + apply_perm(Bint, sL) + apply_perm(Bint, sR) + apply_perm(apply_perm(Bint, sL), sR)  # = 16 L lam1 Z
    del Bint
    log(f'QA\' = P QA P and Z = P[W(x)I - tI - Gamma_A(QA\') - Gamma_B(QB)]P computed exactly ({time.time() - T0:.1f} s)')
    ops = {'QB': QB, "QA'": QAp, 'Z': Zint}
    sym = True
    for name, X in ops.items():
        assert np.array_equal(apply_perm(X, perms['tr']), X), name + ' not symmetric'
        for ax in range(1, n):
            p = list(range(n)); p[0], p[ax] = p[ax], p[0]
            if not np.array_equal(X.transpose(p), X):
                sym = False
    log('QB, QA\', Z are exactly symmetric (transpose-invariant coefficient tensors); copy-permutation invariant:', sym)
    betas = [b for b in itertools.product(range(4), repeat=n) if (not sym) or list(b) == sorted(b)]
    log(f'blocks checked per operator: {len(betas)}' + (' (sorted spin labels; the others are permutation conjugates)' if sym else ''))
    allok = True
    for name, X in ops.items():
        if all(v == 0 for v in X.ravel()):
            log(f'{name}: identically zero (PSD)')
            continue
        dims = nc = nl = ncert = 0; refuted = False; kern = 0; tb = time.time()
        for beta in betas:
            G = compress(X, beta, dcomp)
            assert np.array_equal(G, G.T)
            sg = np.array([int(np.prod([signs[j][i] for j, i in zip(beta, ii)]))
                           for ii in itertools.product(*[range(len(signs[j])) for j in beta])])
            f = None
            if all(j in (1, 3) for j in beta):
                f = forced[beta[0]]
                for j in beta[1:]:
                    f = [u * v for u in f for v in forced[j]]
            if name != 'QB':                                   # P X P: the S = -1 rows must vanish exactly
                neg = np.where(sg < 0)[0]; pos = np.where(sg > 0)[0]
                assert all(v == 0 for v in G[neg, :].ravel()), (name, beta)
                G = G[np.ix_(pos, pos)]
                if f is not None:
                    assert all(f[i] == 0 for i in neg)
                    f = [f[i] for i in pos]
            dims += G.shape[0]
            okc, kc = psd_congruence(G, f)
            nc += okc
            okl = None
            if G.shape[0] <= ldl_max:
                okl, kl = psd_bareiss(G); nl += 1
                if okl is False:
                    refuted = True; log(f'  {name} block {beta}: exact elimination finds it NOT positive semidefinite')
                elif (f is None and kl != 0) or (f is not None and kl != 1):
                    log(f'  {name} block {beta}: kernel dimension {kl} (forced: {0 if f is None else 1})')
                if okl:
                    kern += kl
                if okc and okl and kc != kl:
                    log(f'  {name} block {beta}: kernel dimensions differ between the tests ({kc} vs {kl})')
            elif okc:
                kern += kc
            if okc or okl:
                ncert += 1
            else:
                log(f'  {name} block {beta} (size {G.shape[0]}): not certified by either test')
        ok = (ncert == len(betas)) and not refuted
        log(f'{name}: {"PSD" if ok else "NOT CERTIFIED"}: {ncert}/{len(betas)} blocks certified (dimension {dims}, '
            f'exact kernel dimension {kern}); congruence test passed on {nc}, exact elimination run on {nl} blocks of size '
            f'<= {ldl_max}{", one refuted" if refuted else ", none refuted"} ({time.time() - tb:.1f} s)')
        allok = allok and ok
    if allok:
        log(f'PASS: for all x, y in C^{3 ** n}: <x(x)y| C_{a}^(x){n} |x(x)y> >= {t} |x|^2 |y|^2')
    else:
        log('FAIL')
    log(f'total time {time.time() - T0:.1f} s')


if __name__ == '__main__':
    main()
