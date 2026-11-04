"""Equal-weight (R_inf)_j, letter-1 torus block, reduced by S_k symmetry.

Everything follows Lemma B of relative-anchors/README.md (unrescaled, scaled by 2^Z, Z = k - j), at w_t = 1.
Column strings are Boolean (letters 0/1) bitmasks on k sites.  Row a = (letter 1 at site a) or vac.

Reduction: H_i = (S_2 wr S_i) x S_m acting on pairs (2l, 2l+1), l < i, and rest R = {2i..k-1}, m = k - 2i;
chi = sign on each pair swap, trivial on pair permutations and on S_R.  The (H_i, chi)-twisted invariants of the
letter-1 block are exactly the multiplicity spaces of the S_k-irreps (k-i, i) and (k-i-1, i, 1) (Pieri), so
the block is PSD iff the Gram form Q_i on these (at most 7) vectors is PSD for every i.
"""
import flint
from fractions import Fraction

fq = flint.fmpq
HALF = fq(1, 2)


def popc(x):
    return bin(x).count('1')


class Level:
    """operators on Boolean strings of k sites"""

    def __init__(self, k):
        self.k = k

    # vectors: dict mask -> fmpq
    @staticmethod
    def add(acc, v, c=1):
        for x, val in v.items():
            nv = acc.get(x, 0) + c * val
            if nv == 0:
                acc.pop(x, None)
            else:
                acc[x] = nv
        return acc

    def U(self, v):
        out = {}
        k = self.k
        for x, val in v.items():
            for t in range(k):
                if not (x >> t) & 1:
                    y = x | (1 << t)
                    out[y] = out.get(y, 0) + val
        return {a: b for a, b in out.items() if b != 0}

    def D(self, v):
        out = {}
        k = self.k
        for x, val in v.items():
            for t in range(k):
                if (x >> t) & 1:
                    y = x & ~(1 << t)
                    out[y] = out.get(y, 0) + val
        return {a: b for a, b in out.items() if b != 0}

    @staticmethod
    def sm(a, v):  # sigma_a^- = |0><1|_a
        return {x & ~(1 << a): val for x, val in v.items() if (x >> a) & 1}

    @staticmethod
    def sp(a, v):  # sigma_a^+ = |1><0|_a
        return {x | (1 << a): val for x, val in v.items() if not (x >> a) & 1}

    @staticmethod
    def nop(a, v):  # n_a = |1><1|_a
        return {x: val for x, val in v.items() if (x >> a) & 1}

    @staticmethod
    def Nop(a, v):  # N_a = 3/2 |0><0| + 1/2 |1><1|
        return {x: (HALF if (x >> a) & 1 else fq(3, 2)) * val for x, val in v.items()}

    def T(self, a, v):  # calT_a = 1/2 sigma_a^- U + N_a   (Lemma N, equal weights)
        out = {x: HALF * val for x, val in self.sm(a, self.U(v)).items()}
        return self.add(out, self.Nop(a, v))

    def Einv(self, c, v):
        """E_c^{-1} v, E_c = (3k-2c)/2 + 1/2 DU on level c; Lagrange interpolation in DU."""
        if not v:
            return {}
        k = self.k
        d = min(c, k - c)
        lam = [fq((k - c - ip) * (c + 1 - ip)) for ip in range(d + 1)]
        alpha = fq(3 * k - 2 * c, 2)
        f = [1 / (alpha + HALF * l) for l in lam]
        # Newton divided differences of f at nodes lam
        n = len(lam)
        coef = list(f)
        for lev in range(1, n):
            for idx in range(n - 1, lev - 1, -1):
                coef[idx] = (coef[idx] - coef[idx - 1]) / (lam[idx] - lam[idx - lev])
        # p(x) = coef0 + coef1 (x-l0) + coef2 (x-l0)(x-l1) + ...
        out = {}
        w = dict(v)
        for idx in range(n):
            self.add(out, w, coef[idx])
            if idx < n - 1:
                DUw = self.D(self.U(w))
                w = self.add(DUw, w, -lam[idx])
        return out


def dot(u, v):
    if len(u) > len(v):
        u, v = v, u
    s = fq(0)
    for x, val in u.items():
        w = v.get(x)
        if w is not None:
            s += val * w
    return s


def subsets(sites, r):
    from itertools import combinations
    for c in combinations(sites, r):
        m = 0
        for t in c:
            m |= 1 << t
        yield m


def singlet_prod(pairs):
    """dict mask -> sign for prod over pairs (|1 at first> - |1 at second>)"""
    out = {0: 1}
    for (p, q) in pairs:
        new = {}
        for x, s in out.items():
            new[x | (1 << p)] = s
            new[x | (1 << q)] = -s
        out = new
    return out


def tensor(base, sites, cnt):
    """base (dict mask->coef) times Dicke(sites, cnt)"""
    out = {}
    if cnt < 0 or cnt > len(sites):
        return out
    for y in subsets(sites, cnt):
        for x, s in base.items():
            out[x | y] = out.get(x | y, 0) + s
    return out


def twisted_vectors(k, j, i):
    """returns list of (name, psi: dict site -> level-j vector, phi: level-(j+1) vector)"""
    pairs = [(2 * l, 2 * l + 1) for l in range(i)]
    R = list(range(2 * i, k))
    vecs = []
    sP = singlet_prod(pairs)
    # A_x : marked site in R
    for x in (0, 1):
        psi = {}
        for a in R:
            base = {(m | (x << a)): s for m, s in sP.items()}
            v = tensor(base, [t for t in R if t != a], j - i - x)
            if v:
                psi[a] = {y: fq(c) for y, c in v.items()}
        if psi:
            vecs.append((f'A{x}', psi, {}))
    # B_xy : marked site in a pair
    for x in (0, 1):
        for y in (0, 1):
            psi = {}
            for l, (p, q) in enumerate(pairs):
                others = [pp for ll, pp in enumerate(pairs) if ll != l]
                sO = singlet_prod(others)
                for a, ap, sg in ((p, q, 1), (q, p, -1)):
                    base = {(m | (x << a) | (y << ap)): s * sg for m, s in sO.items()}
                    v = tensor(base, R, j - i + 1 - x - y)
                    if v:
                        psi[a] = {yy: fq(c) for yy, c in v.items()}
            if psi:
                vecs.append((f'B{x}{y}', psi, {}))
    # G0
    if j + 1 <= k:
        v = tensor(sP, R, j + 1 - i)
        if v:
            vecs.append(('G0', {}, {yy: fq(c) for yy, c in v.items()}))
    return vecs


def orbit_reps(k, i):
    """H_i-orbits on ordered pairs (a, b) of sites: list of (a, b, size)"""
    m = k - 2 * i
    R = list(range(2 * i, k))
    reps = []
    if m >= 1:
        reps.append((R[0], R[0], m))
    if m >= 2:
        reps.append((R[0], R[1], m * (m - 1)))
    if m >= 1 and i >= 1:
        reps.append((R[0], 0, 2 * i * m))
        reps.append((0, R[0], 2 * i * m))
    if i >= 1:
        reps.append((0, 0, 2 * i))
        reps.append((0, 1, 2 * i))
    if i >= 2:
        reps.append((0, 2, 2 * i * (2 * i - 2)))
    return reps


def site_reps(k, i):
    m = k - 2 * i
    reps = []
    if m >= 1:
        reps.append((2 * i, m))
    if i >= 1:
        reps.append((0, 2 * i))
    return reps


def reduced(k, j, i, check_full=False):
    """Gram matrix G and form matrix Q (2^Z-scaled Lemma B form) on the twisted vectors."""
    L = Level(k)
    vecs = twisted_vectors(k, j, i)
    nv = len(vecs)
    pr = orbit_reps(k, i)
    sr = site_reps(k, i)
    G = [[fq(0)] * nv for _ in range(nv)]
    Q = [[fq(0)] * nv for _ in range(nv)]
    # precompute E_j^{-1} T_a psi^beta_b for orbit reps
    cacheE = {}
    for be, (_, psib, phib) in enumerate(vecs):
        for (a, b, sz) in pr:
            pb = psib.get(b)
            if pb:
                cacheE[(be, a, b)] = L.Einv(j, L.T(a, pb))
    # coupling vectors: (-sigma_a^- + D E_{j+1}^{-1} T_a) phi^beta
    cacheC = {}
    for be, (_, psib, phib) in enumerate(vecs):
        if phib:
            for (a, sz) in sr:
                v = {x: -c for x, c in L.sm(a, phib).items()}
                L.add(v, L.D(L.Einv(j + 1, L.T(a, phib))))
                cacheC[(be, a)] = v
    for al in range(nv):
        _, psia, phia = vecs[al]
        for be in range(nv):
            _, psib, phib = vecs[be]
            g = fq(0)
            q = fq(0)
            for (a, sz) in sr:
                pa, pb = psia.get(a), psib.get(a)
                if pa and pb:
                    g += sz * dot(pa, pb)
                    # (2 - n_a) term
                    q += sz * (2 * dot(pa, pb) - dot(pa, L.nop(a, pb)))
            if phia and phib:
                g += dot(phia, phib)
            for (a, b, sz) in pr:
                pa = psia.get(a)
                pb = psib.get(b)
                if pa and pb:
                    coef = HALF - (1 if a == b else 0)
                    q += sz * coef * dot(L.sm(b, pa), L.sm(a, pb))
                    q -= sz * dot(L.T(b, pa), cacheE[(be, a, b)])
            for (a, sz) in sr:
                pa = psia.get(a)
                if pa and phib:
                    q += sz * dot(pa, cacheC[(be, a)])
                pb = psib.get(a)
                if pb and phia:
                    q += sz * dot(pb, cacheC[(al, a)])
            if phia and phib:
                q += 2 * dot(phia, phib)
                if j + 2 <= k:
                    q -= dot(L.U(phia), L.Einv(j + 2, L.U(phib)))
            G[al][be] = g
            Q[al][be] = q
    return [v[0] for v in vecs], G, Q, vecs


def ldl_psd(Mx):
    """exact symmetric elimination with diagonal pivoting; returns (is_psd, rank, pivots)."""
    n = len(Mx)
    M = [row[:] for row in Mx]
    act = list(range(n))
    pivs = []
    while act:
        # choose a nonzero diagonal
        cand = [t for t in act if M[t][t] != 0]
        if not cand:
            # all remaining diagonals zero: PSD iff remaining block zero
            ok = all(M[a][b] == 0 for a in act for b in act)
            return ok, len(pivs), pivs
        t = cand[0]
        piv = M[t][t]
        pivs.append(piv)
        if piv < 0:
            return False, len(pivs), pivs
        act.remove(t)
        for a in act:
            f = M[a][t] / piv
            if f != 0:
                for b in act:
                    M[a][b] -= f * M[t][b]
    return True, len(pivs), pivs
