"""Isotypic reduction at fixed level sizes by Young symmetrizers.

G = S_{u0} x S_{u1} permutes the copies of each level (level 0 = copies 0..u0-1, level 1 = copies u0..k-1).  For an
irrep lam = (lam0, lam1) take the tableaux filled row by row with the copies of each level, R = row group, C = column
group, and the Young symmetrizer e = (sum_{s in R} s)(sum_{t in C} sgn(t) t).  On any permutation module M, e is zero on
the other isotypic components and has rank one on each copy of the Specht module, so e.M is a space of dimension
mult(lam, M) whose vectors are s0 (x) c (s0 a fixed vector of the irrep).  For a G-invariant symmetric Y and any basis
B of e.M, B^T Y B = |s0|^2 C^T Q_lam C with C invertible, where Y on the lam-component is <,> (x) Q_lam.

States: P-states (r, X) (X a pair) and Q-states T (triples), as in l2block.  Vectors are dicts or dense lists over the
P-states (or the Q-states).  The row symmetrization is done with orbit sums: sum_{s in R} e_{s y} = |Stab_R(y)| 1_{R y},
and R-orbits are determined by the row (block) of each index and the coincidences."""
import itertools
from fractions import Fraction as Fr
from math import factorial


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield (p,) + rest


def hook_dim(lam):
    n = sum(lam)
    if n == 0:
        return 1
    conj = [sum(1 for x in lam if x > j) for j in range(lam[0])] if lam else []
    h = 1
    for i, li in enumerate(lam):
        for j in range(li):
            h *= (li - j - 1) + (conj[j] - i - 1) + 1
    return factorial(n) // h


def off(lam):
    return sum(lam) - (lam[0] if lam else 0)


class Young:
    def __init__(self, k, u0, lam0, lam1):
        self.k, self.u0 = k, u0
        levels = [list(range(u0)), list(range(u0, k))]
        self.rows, self.cols = [], []
        for copies, lam in zip(levels, (lam0, lam1)):
            assert sum(lam) == len(copies)
            pos = 0
            rws = []
            for li in lam:
                rws.append(copies[pos:pos + li])
                pos += li
            self.rows += rws
            for j in range(lam[0] if lam else 0):
                col = [rw[j] for rw in rws if len(rw) > j]
                if len(col) > 1:
                    self.cols.append(col)
        self.block = {}
        for bi, rw in enumerate(self.rows):
            for c in rw:
                self.block[c] = bi
        # column group: product over columns of all permutations of the column, with signs
        self.C = [({}, 1)]
        for col in self.cols:
            new = []
            for perm in itertools.permutations(col):
                sg = perm_sign(col, perm)
                mp = dict(zip(col, perm))
                for (m0, s0) in self.C:
                    m = dict(m0)
                    m.update(mp)
                    new.append((m, s0 * sg))
            self.C = new
        self.Pst = [(r, X) for r in range(k) for X in itertools.combinations(range(k), 2)]
        self.Qst = list(itertools.combinations(range(k), 3))
        self.okeys_P = [self.key_P(x) for x in self.Pst]
        self.okeys_Q = [self.key_Q(x) for x in self.Qst]
        self.osize_P, self.osize_Q = {}, {}
        for kk in self.okeys_P:
            self.osize_P[kk] = self.osize_P.get(kk, 0) + 1
        for kk in self.okeys_Q:
            self.osize_Q[kk] = self.osize_Q.get(kk, 0) + 1
        self.olist_P = sorted(self.osize_P)
        self.olist_Q = sorted(self.osize_Q)

    def key_P(self, st):
        r, X = st
        b = self.block
        if r in X:
            c = X[0] if X[1] == r else X[1]
            return ('in', b[r], b[c])
        return ('out', b[r], tuple(sorted((b[X[0]], b[X[1]]))))

    def key_Q(self, T):
        return tuple(sorted(self.block[t] for t in T))

    def act_P(self, m, st):
        r, X = st
        g = lambda c: m.get(c, c)
        return (g(r), tuple(sorted((g(X[0]), g(X[1])))))

    def act_Q(self, m, T):
        return tuple(sorted(m.get(c, c) for c in T))

    def sym_P(self, v):
        """e.v for a dict v over P-states, in orbit coordinates (dict orbit key -> coefficient of 1_orbit), up to the
        common factor |R|"""
        out = {}
        for st, c in v.items():
            for m, sg in self.C:
                y = self.act_P(m, st)
                kk = self.key_P(y)
                out[kk] = out.get(kk, 0) + Fr(sg) * c / self.osize_P[kk]
        return {kk: c for kk, c in out.items() if c != 0}

    def sym_Q(self, v):
        out = {}
        for st, c in v.items():
            for m, sg in self.C:
                y = self.act_Q(m, st)
                kk = self.key_Q(y)
                out[kk] = out.get(kk, 0) + Fr(sg) * c / self.osize_Q[kk]
        return {kk: c for kk, c in out.items() if c != 0}

    def dense_P(self, ov):
        return [ov.get(kk, 0) for kk in self.okeys_P]

    def dense_Q(self, ov):
        return [ov.get(kk, 0) for kk in self.okeys_Q]


def perm_sign(seq, perm):
    idx = {x: i for i, x in enumerate(seq)}
    p = [idx[x] for x in perm]
    sg = 1
    seen = [False] * len(p)
    for i in range(len(p)):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            if L % 2 == 0:
                sg = -sg
    return sg


def rref_basis(vectors, keys):
    """basis (reduced row echelon form, integer-scaled rows) of the span of the dict vectors over the given keys"""
    rows = []
    pivots = []
    for v in vectors:
        r = [Fr(v.get(kk, 0)) for kk in keys]
        for (pc, pr) in zip(pivots, rows):
            if r[pc] != 0:
                f = r[pc]
                r = [a - f * b for a, b in zip(r, pr)]
        nz = next((i for i, a in enumerate(r) if a != 0), None)
        if nz is None:
            continue
        f = r[nz]
        r = [a / f for a in r]
        for i, pr in enumerate(rows):
            if pr[nz] != 0:
                g = pr[nz]
                rows[i] = [a - g * b for a, b in zip(pr, r)]
        rows.append(r)
        pivots.append(nz)
    order = sorted(range(len(rows)), key=lambda i: pivots[i])
    return [rows[i] for i in order], [pivots[i] for i in order]


def extend_basis(base, cands):
    """greedily extend the list of dense vectors base by candidates (dense) that raise the rank; returns the added ones"""
    rows, piv = [], []

    def reduce(r):
        r = list(r)
        for (pc, pr) in zip(piv, rows):
            if r[pc] != 0:
                f = r[pc] / pr[pc]
                r = [a - f * b for a, b in zip(r, pr)]
        return r

    for v in base:
        r = reduce([Fr(x) for x in v])
        nz = next((i for i, a in enumerate(r) if a != 0), None)
        assert nz is not None, 'base not independent'
        rows.append(r)
        piv.append(nz)
    added = []
    for v in cands:
        r = reduce([Fr(x) for x in v])
        nz = next((i for i, a in enumerate(r) if a != 0), None)
        if nz is None:
            continue
        rows.append(r)
        piv.append(nz)
        added.append(v)
    return added


def irreps(u0, u1, maxoff=3):
    """the irreps (lam0, lam1) of S_u0 x S_u1 with at most maxoff boxes outside the first rows (the only ones that can
    occur in modules of states with at most three indices); the dimension count in the callers checks this."""
    out = []
    for l0 in partitions(u0):
        for l1 in partitions(u1):
            if off(l0) + off(l1) <= maxoff:
                out.append((l0, l1))
    return out


def part_data(k, u0, lam0, lam1):
    """test vectors for the irrep (lam0, lam1): bases of e.M_P (P-states), e.K (kernel span) and a complement of e.K in
    e.M_P built from e.W, W = span{(r, X) : r not in X} + span{(a, X) - (b, X) : X = {a, b}} (a G-invariant complement of
    K in M_P), and a basis of e.M_Q.  Vectors are dense lists over the P-states (resp. Q-states)."""
    Y = Young(k, u0, lam0, lam1)
    keysP = Y.olist_P
    keysQ = Y.olist_Q
    eM = [Y.sym_P({x: 1}) for x in Y.Pst]
    basisM, _ = rref_basis(eM, keysP)
    mP = len(basisM)
    pairs = list(itertools.combinations(range(k), 2))
    eK = [Y.sym_P({(r, X): 1 for r in range(k)}) for X in pairs]
    basisK, _ = rref_basis(eK, keysP)
    d = len(basisK)
    Wgen = []
    for X in pairs:
        for r in range(k):
            if r not in X:
                Wgen.append({(r, X): 1})
        Wgen.append({(X[0], X): 1, (X[1], X): -1})
    eW = [Y.sym_P(v) for v in Wgen]
    basisW, _ = rref_basis(eW, keysP)
    eQ = [Y.sym_Q({T: 1}) for T in Y.Qst]
    basisQ, _ = rref_basis(eQ, keysQ)

    def todense_P(row):
        ov = dict(zip(keysP, row))
        return [ov[kk] for kk in Y.okeys_P]

    def todense_Q(row):
        ov = dict(zip(keysQ, row))
        return [ov[kk] for kk in Y.okeys_Q]

    # integer scaling of each vector (positive factor; spans unchanged)
    def intscale(v):
        from math import lcm
        den = 1
        for x in v:
            den = lcm(den, Fr(x).denominator)
        return [int(Fr(x) * den) for x in v]

    KP = [intscale(todense_P(r)) for r in basisK]
    WP = [intscale(todense_P(r)) for r in basisW]
    QQ = [intscale(todense_Q(r)) for r in basisQ]
    return dict(Y=Y, mP=mP, d=d, K=KP, W=WP, Q=QQ, dim=hook_dim(lam0) * hook_dim(lam1))
