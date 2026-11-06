"""The letter-1 block of (R_inf)_2 at two weight levels, with the level sizes symbolic.

Notation (NOTES.md, Sections 1-3).  Level 0: u0 copies of squared weight 1; level 1: u1 copies of squared weight q.
Rescaled frame and units 2^{k-2} as in middle_j2_l1 (Lemma A).  States: P = (r, X) (X a pair), Q = T (a triple).

Operators (all G-invariant, G = S_u0 x S_u1), written on joint index patterns (slots with levels):
  Y_PP   Lemma A item 1, with Ehat_2^{-1} by Woodbury through the singles (NOTES Lemma 1);
  A2     the pair operator Ahat_2 (NOTES Lemma 2), used as I (x) A2 on the P-states;
  K      the inverse-free P-Q operator with Y_PQ = (I (x) A2^{-1}) K (NOTES Lemma 2);
  W3     What_3 of Lemma V (Y_QQ = 2 W3^{-1}).
The congruence by I (x) A2 turns the Schur complement S = Y_PP - 1/2 Y_PQ W3 Y_QP into
  S' = (I (x) A2) Y_PP (I (x) A2) - 1/2 K W3 K^T,
which has the same kernel 1 (x) M_pairs (NOTES Lemma 3).

Reduction: for an irrep lambda of G, a tableau with named copies (the copies in columns of length >= 2), C_T its column
group, R_T its row group; test vectors C_T 1_O for R_T-orbits O of states (they span e'_lambda M, e' = C_T R_T).  Every
reduced matrix is a pattern sum: sum over joint index patterns relative to the named copies, with realisation counts
prod_l ff(u_l - n_l, m_l), of (coefficient) x (operator value on the pattern).  All arithmetic is generic over a field
of numbers (fractions, or flint rational functions of q and the level sizes)."""
import itertools
from fractions import Fraction as Fr
import numpy as np


# ---------------------------------------------------------------- structured matrices on the singles
class SM:
    """M[c, c'] = [c = c'] a[l(c)] + b[l(c)][l(c')]  (a two-level invariant operator on the singles)."""
    __slots__ = ('a', 'b')

    def __init__(self, a, b):
        self.a, self.b = a, b


def sm_mul(X, Y, u):
    a = [X.a[l] * Y.a[l] for l in (0, 1)]
    b = [[X.a[l] * Y.b[l][m] + X.b[l][m] * Y.a[m] + X.b[l][0] * u[0] * Y.b[0][m] + X.b[l][1] * u[1] * Y.b[1][m]
          for m in (0, 1)] for l in (0, 1)]
    return SM(a, b)


def sm_inv(X, u, one):
    M = [[(X.a[l] if l == m else one * 0) + X.b[l][m] * u[m] for m in (0, 1)] for l in (0, 1)]
    det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    Mi = [[M[1][1] / det, -M[0][1] / det], [-M[1][0] / det, M[0][0] / det]]
    ia = [one / X.a[0], one / X.a[1]]
    bt = [[-(Mi[l][0] * X.b[0][m] + Mi[l][1] * X.b[1][m]) * ia[m] for m in (0, 1)] for l in (0, 1)]
    return SM(ia, bt)


# ---------------------------------------------------------------- the operators on index patterns
class Ctx:
    """one: field unit; q: squared weight of level 1; u0, u1: level sizes (field elements); bscale: factor on every
    term of Lemma B that contains an E^{-1} (only the numeric negative controls use bscale != 1)."""

    def __init__(self, one, q, u0, u1):
        self.one, self.zero = one, one * 0
        self.half = one / 2
        self.s = [one, q]
        self.u = [u0, u1]
        L = self.L = u0 + q * u1
        s, u = self.s, self.u
        self.D2l = [[2 * (L - s[l] - s[m]) / (s[l] * s[m]) for m in (0, 1)] for l in (0, 1)]
        bN = [[one / self.D2l[l][m] for m in (0, 1)] for l in (0, 1)]
        diagN = [(u[0] - (1 if l == 0 else 0)) * bN[l][0] + (u[1] - (1 if l == 1 else 0)) * bN[l][1] for l in (0, 1)]
        aN = [diagN[l] - bN[l][l] for l in (0, 1)]
        self.N = SM(aN, bN)
        self.Psi = SM([s[l] + self.half * aN[l] for l in (0, 1)], [[self.half * bN[l][m] for m in (0, 1)] for l in (0, 1)])
        self.Pi = sm_inv(self.Psi, u, one)
        self.PiN = sm_mul(self.Pi, self.N, u)
        self.NPiN = sm_mul(self.N, self.PiN, u)

    # helpers on sets of slots (lev: slot -> level)
    def se(self, X, c, d, lev):
        v = X.b[lev[c]][lev[d]]
        return v + X.a[lev[c]] if c == d else v

    def pi(self, Z, lev):
        v = self.one
        for z in Z:
            v = v * self.s[lev[z]]
        return v

    def sig(self, Z, lev):
        v = self.zero
        for z in Z:
            v = v + self.s[lev[z]]
        return v

    def d2(self, Z, lev):
        return 2 * (self.L - self.sig(Z, lev)) / self.pi(Z, lev)

    def dsum(self, Z, lev):
        """sum_{t not in Z} s_t / (2(L - sig_Z - s_t))"""
        sg = self.sig(Z, lev)
        v = self.zero
        for l in (0, 1):
            n = sum(1 for z in Z if lev[z] == l)
            c = self.u[l] - n
            if c == 0:
                continue                 # no copy of this level outside Z
            v = v + c * self.s[l] / (2 * (self.L - sg - self.s[l]))
        return v

    # Lemma A item 1 (row state (r, Xp), column state (rp, X)), Ehat_2^{-1} by Woodbury (NOTES Lemma 1)
    def pp(self, r, Xp, rp, X, lev, b=1):
        s = self.s
        sX, sXp = set(X), set(Xp)
        val = self.zero
        if r in sX:
            rest = sX - {r}
            if rp not in rest and (rest | {rp}) == sXp:
                val = val + (self.half - (1 if r == rp else 0)) * s[lev[rp]] * self.pi(X, lev)
        if r == rp and sX == sXp:
            val = val + (2 - (1 if r in sX else 0)) * s[lev[r]] * self.pi(X, lev)
        rin, rpin = r in sX, rp in sXp
        if rin and rpin:
            c = next(iter(sX - {r}))
            cp = next(iter(sXp - {rp}))
            B = self.se(self.N, cp, c, lev) / 4 - self.se(self.NPiN, cp, c, lev) / 8
        elif rin:
            c = next(iter(sX - {r}))
            t = self.zero
            for y in Xp:
                t = t + self.se(self.PiN, y, c, lev)
            B = ((1 if c in sXp else 0) * self.one - self.half * t) / self.d2(Xp, lev)
        elif rpin:
            cp = next(iter(sXp - {rp}))
            t = self.zero
            for x in X:
                t = t + self.se(self.PiN, x, cp, lev)
            B = ((1 if cp in sX else 0) * self.one - self.half * t) / self.d2(X, lev)
        else:
            t = self.zero
            for x in X:
                for y in Xp:
                    t = t + self.se(self.Pi, y, x, lev)
            dX = self.d2(X, lev)
            B = ((4 if sX == sXp else 0) * self.one) / dX - 2 * t / (dX * self.d2(Xp, lev))
        return val - b * s[lev[r]] * s[lev[rp]] * B

    # the pair operator Ahat_2 (and What_3 on triples): same pattern, m = |Z|
    def adj(self, Zp, Z, lev):
        sZ, sZp = set(Z), set(Zp)
        if sZ == sZp:
            return (self.one + self.half * self.dsum(sZ, lev)) / self.pi(sZ, lev)
        inter = sZ & sZp
        if len(inter) == len(sZ) - 1:
            return self.one / (4 * self.pi(inter, lev) * (self.L - self.sig(sZ | sZp, lev)))
        return self.zero

    def kk(self, r, Z, T, lev):
        """K[(r, Z), T]: Y_PQ = (I (x) Ahat_2^{-1}) K"""
        sT, sZ = set(T), set(Z)
        if r in sT:
            return -self.s[lev[r]] * self.one if sZ == sT - {r} else self.zero
        if sZ <= sT:
            t = next(iter(sT - sZ))
            return self.s[lev[r]] * self.s[lev[t]] / (self.L - self.sig(sT, lev))
        return self.zero


# ---------------------------------------------------------------- tableaux, test vectors, pattern sums
def perms_with_sign(elems):
    out = []
    for p in itertools.permutations(range(len(elems))):
        sg = 1
        for i in range(len(p)):
            for j in range(i + 1, len(p)):
                if p[i] > p[j]:
                    sg = -sg
        out.append(({elems[i]: elems[p[i]] for i in range(len(p))}, sg))
    return out


class Family:
    """An irrep (lambda0, lambda1) of S_u0 x S_u1 with at most three boxes outside the first rows.
    named: list of (level, row) of the named copies; columns: lists of named indices (column heads first, row 1)."""

    def __init__(self, name, named, columns, shape):
        self.name, self.named, self.columns, self.shape = name, named, columns, shape
        self.N = len(named)
        self.nlev = [sum(1 for (l, _) in named if l == lv) for lv in (0, 1)]
        kap = [({}, 1)]
        for col in columns:
            new = []
            for (d, sg) in kap:
                for (e, sg2) in perms_with_sign(col):
                    dd = dict(d)
                    dd.update(e)
                    new.append((dd, sg * sg2))
            kap = new
        # store the inverse permutations
        self.kinv = []
        for (d, sg) in kap:
            inv = {v: k for k, v in d.items()}
            self.kinv.append((inv, sg))

    def rc(self, tok):
        if isinstance(tok, int):
            return self.named[tok]
        return (tok[1], 1)

    def act(self, inv, p):
        return tuple(inv.get(t, t) if isinstance(t, int) else t for t in p)


def local_patterns(N, nposs, distinct):
    """tuples of tokens: named index (int) or ('A', level, i) anonymous (i by first appearance); distinct: list of
    position pairs that must differ."""
    out = []

    def rec(pos, cur, nanon, levs):
        if pos == nposs:
            ok = all(cur[i] != cur[j] for (i, j) in distinct)
            if ok:
                out.append(tuple(cur))
            return
        for j in range(N):
            rec(pos + 1, cur + [j], nanon, levs)
        for i in range(nanon):
            rec(pos + 1, cur + [('A', levs[i], i)], nanon, levs)
        for l in (0, 1):
            rec(pos + 1, cur + [('A', l, nanon)], nanon + 1, levs + [l])

    rec(0, [], 0, [])
    return out


P_DIST = [(1, 2)]
Q_DIST = [(0, 1), (0, 2), (1, 2)]


def pkey(F, p):
    r, x1, x2 = (F.rc(t) for t in p)
    if p[0] == p[1]:
        return ('in', r, x2)
    if p[0] == p[2]:
        return ('in', r, x1)
    return ('out', r, tuple(sorted((x1, x2))))


def candidates(F):
    """candidate test vectors (dicts pattern -> int) for the K-part, the W-part (P-states) and the Q-states."""
    Ppat = local_patterns(F.N, 3, P_DIST)
    Qpat = local_patterns(F.N, 3, Q_DIST)
    coefP, coefX = {}, {}
    for p in Ppat:
        cp, cx = {}, {}
        for (inv, sg) in F.kinv:
            pp_ = F.act(inv, p)
            O = pkey(F, pp_)
            cp[O] = cp.get(O, 0) + sg
            Ox = tuple(sorted((F.rc(pp_[1]), F.rc(pp_[2]))))
            cx[Ox] = cx.get(Ox, 0) + sg
        coefP[p] = cp
        coefX[p] = cx
    outs = sorted({O for c in coefP.values() for O in c if O[0] == 'out'})
    ins = sorted({O for c in coefP.values() for O in c if O[0] == 'in'})
    pairs = sorted({O for c in coefX.values() for O in c})
    Wc = []
    for O in outs:
        Wc.append((('W', O), {p: c.get(O, 0) for p, c in coefP.items() if c.get(O, 0)}))
    for O in ins:
        O2 = ('in', O[2], O[1])
        if O[1] < O[2]:
            v = {p: c.get(O, 0) - c.get(O2, 0) for p, c in coefP.items() if c.get(O, 0) - c.get(O2, 0)}
            Wc.append((('W', O, '-'), v))
    Kc = []
    for O in pairs:
        Kc.append((('K', O), {p: c.get(O, 0) for p, c in coefX.items() if c.get(O, 0)}))
    coefQ = {}
    for T in Qpat:
        cq = {}
        for (inv, sg) in F.kinv:
            O = tuple(sorted(F.rc(t) for t in F.act(inv, T)))
            cq[O] = cq.get(O, 0) + sg
        coefQ[T] = cq
    trip = sorted({O for c in coefQ.values() for O in c})
    Qc = [(('Q', O), {T: c.get(O, 0) for T, c in coefQ.items() if c.get(O, 0)}) for O in trip]
    Kc = [c for c in Kc if c[1]]
    Wc = [c for c in Wc if c[1]]
    Qc = [c for c in Qc if c[1]]
    return Kc, Wc, Qc


def ff(x, m, one):
    v = one
    for i in range(m):
        v = v * (x - i)
    return v


def nanon(p):
    a = [0, 0]
    seen = set()
    for t in p:
        if not isinstance(t, int) and t[2] not in seen:
            seen.add(t[2])
            a[t[1]] += 1
    return a


def gram_numeric(F, vecs, u, factor):
    """exact Gram matrix at integer sizes u = (u0, u1): factor * sum_p count(p) v_i(p) v_j(p)"""
    pats = sorted({p for (_, v) in vecs for p in v}, key=str)
    n = len(vecs)
    G = [[Fr(0)] * n for _ in range(n)]
    for p in pats:
        a = nanon(p)
        cnt = ff(u[0] - F.nlev[0], a[0], 1) * ff(u[1] - F.nlev[1], a[1], 1)
        if cnt == 0:
            continue
        vals = [v.get(p, 0) for (_, v) in vecs]
        for i in range(n):
            if vals[i]:
                for j in range(n):
                    if vals[j]:
                        G[i][j] += factor * cnt * vals[i] * vals[j]
    return G


def rank_fr(M):
    A = [list(r) for r in M]
    n = len(A)
    m = len(A[0]) if n else 0
    rk = 0
    for c in range(m):
        piv = next((i for i in range(rk, n) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        for i in range(n):
            if i != rk and A[i][c] != 0:
                f = A[i][c] / A[rk][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[rk])]
        rk += 1
    return rk


def select(F, cands, u, factor, base=()):
    """greedy maximal independent subset (exact Gram ranks at the integer sizes u), after the vectors in base"""
    chosen = list(base)
    for c in cands:
        trial = chosen + [c]
        G = gram_numeric(F, trial, u, factor)
        if rank_fr(G) == len(trial):
            chosen = trial
    return chosen[len(base):]


def canon(slots, lev):
    m = {}
    out = []
    levs = []
    for s in slots:
        if s not in m:
            m[s] = len(m)
            levs.append(lev[s])
        out.append(m[s])
    return tuple(out), tuple(levs)


def joint_sum(F, vecs1, vecs2, factor):
    """dict (orbital key, m0, m1) -> integer matrix sum over joint patterns of v_i(p1) w_j(p2), times nothing;
    factor is returned unchanged for the caller.  The key is the canonical form of the 6-tuple p1 + p2."""
    pats1 = sorted({p for (_, v) in vecs1 for p in v}, key=str)
    pats2 = sorted({p for (_, v) in vecs2 for p in v}, key=str)
    V1 = {p: np.array([v.get(p, 0) for (_, v) in vecs1], dtype=np.int64) for p in pats1}
    V2 = {p: np.array([v.get(p, 0) for (_, v) in vecs2], dtype=np.int64) for p in pats2}
    N = F.N
    acc = {}
    info2 = {}
    for p2 in pats2:
        ids, levs = [], []
        for t in p2:
            if not isinstance(t, int) and t[2] not in ids:
                ids.append(t[2])
                levs.append(t[1])
        info2[p2] = (ids, levs)
    for p1 in pats1:
        a1ids, a1lev = [], []
        for t in p1:
            if not isinstance(t, int) and t[2] not in a1ids:
                a1ids.append(t[2])
                a1lev.append(t[1])
        lev = {j: F.named[j][0] for j in range(N)}
        s1 = []
        for t in p1:
            if isinstance(t, int):
                s1.append(t)
            else:
                sl = N + a1ids.index(t[2])
                lev[sl] = t[1]
                s1.append(sl)
        n1 = len(a1ids)
        for p2 in pats2:
            ids2, levs2 = info2[p2]
            # matchings: each anon of p2 -> None (new) or an unused anon of p1 of the same level
            choices = []
            for l2 in levs2:
                choices.append([None] + [i for i in range(n1) if a1lev[i] == l2])
            for mt in itertools.product(*choices):
                used = [x for x in mt if x is not None]
                if len(used) != len(set(used)):
                    continue
                lev2 = dict(lev)
                mapping = {}
                nxt = N + n1
                for idx, x in enumerate(mt):
                    if x is None:
                        mapping[ids2[idx]] = nxt
                        lev2[nxt] = levs2[idx]
                        nxt += 1
                    else:
                        mapping[ids2[idx]] = N + x
                s2 = [t if isinstance(t, int) else mapping[t[2]] for t in p2]
                key = canon(s1 + s2, lev2)
                m = [0, 0]
                for sl in range(N, nxt):
                    m[lev2[sl]] += 1
                kk = (key, m[0], m[1])
                o = np.outer(V1[p1], V2[p2])
                if kk in acc:
                    acc[kk] += o
                else:
                    acc[kk] = o.copy()
    return acc


def assemble(F, ctx, acc, oracle, factor, cache):
    """sum over keys of oracle(key) * count * matrix, as field elements"""
    one = ctx.one
    # combine the counts per orbital key first
    bykey = {}
    ffc = {}
    for (key, m0, m1), M in acc.items():
        if not M.any():
            continue
        if (m0, m1) not in ffc:
            ffc[(m0, m1)] = ff(ctx.u[0] - F.nlev[0], m0, one) * ff(ctx.u[1] - F.nlev[1], m1, one)
        c = ffc[(m0, m1)]
        if key not in bykey:
            bykey[key] = {}
        d = bykey[key]
        nz = np.argwhere(M)
        for (i, j) in nz:
            v = int(M[i, j]) * c
            d[(i, j)] = d[(i, j)] + v if (i, j) in d else v
    n1 = max(i for (_, _, _), M in acc.items() for i in [M.shape[0]]) if acc else 0
    n2 = max(M.shape[1] for M in acc.values()) if acc else 0
    R = [[ctx.zero] * n2 for _ in range(n1)]
    for key, d in bykey.items():
        d = {ij: c for ij, c in d.items() if not (c == 0)}
        if not d:
            continue                     # pattern not realisable at these sizes (count 0)
        if key not in cache:
            cache[key] = oracle(key)
        val = cache[key]
        if val == 0:
            continue
        for (i, j), c in d.items():
            R[i][j] = R[i][j] + val * c
    return [[x * factor for x in row] for row in R]


def oracle_pp(ctx, b=1):
    def f(key):
        sl, levs = key
        lev = dict(enumerate(levs))
        return ctx.pp(sl[0], (sl[1], sl[2]), sl[3], (sl[4], sl[5]), lev, b)
    return f


def oracle_a2(ctx):
    def f(key):
        sl, levs = key
        lev = dict(enumerate(levs))
        if sl[0] != sl[3]:
            return ctx.zero
        return ctx.adj((sl[1], sl[2]), (sl[4], sl[5]), lev)
    return f


def oracle_k(ctx):
    def f(key):
        sl, levs = key
        lev = dict(enumerate(levs))
        return ctx.kk(sl[0], (sl[1], sl[2]), (sl[3], sl[4], sl[5]), lev)
    return f


def oracle_w3(ctx):
    def f(key):
        sl, levs = key
        lev = dict(enumerate(levs))
        return ctx.adj((sl[0], sl[1], sl[2]), (sl[3], sl[4], sl[5]), lev)
    return f


def oracle_id(ctx, kind):
    def f(key):
        sl, levs = key
        if kind == 'P':
            return ctx.one if (sl[0] == sl[3] and {sl[1], sl[2]} == {sl[4], sl[5]}) else ctx.zero
        return ctx.one if {sl[0], sl[1], sl[2]} == {sl[3], sl[4], sl[5]} else ctx.zero
    return f


# ---------------------------------------------------------------- the families
def families():
    """the 16 families of irreps (lambda0, lambda1); 'level' l means the named copies sit on level l."""
    out = []
    for l in (0, 1):
        o = 1 - l
        out.append(Family('std%d' % l, [(l, 1), (l, 2)], [[0, 1]], {l: (1,), o: ()}))
        out.append(Family('two%d' % l, [(l, 1), (l, 2), (l, 1), (l, 2)], [[0, 1], [2, 3]], {l: (2,), o: ()}))
        out.append(Family('anti%d' % l, [(l, 1), (l, 2), (l, 3)], [[0, 1, 2]], {l: (1, 1), o: ()}))
        out.append(Family('three%d' % l, [(l, 1), (l, 2), (l, 1), (l, 2), (l, 1), (l, 2)], [[0, 1], [2, 3], [4, 5]],
                          {l: (3,), o: ()}))
        out.append(Family('tto%d' % l, [(l, 1), (l, 2), (l, 3), (l, 1), (l, 2)], [[0, 1, 2], [3, 4]],
                          {l: (2, 1), o: ()}))
        out.append(Family('twostd%d' % l, [(l, 1), (l, 2), (l, 1), (l, 2), (o, 1), (o, 2)], [[0, 1], [2, 3], [4, 5]],
                          {l: (2,), o: (1,)}))
        out.append(Family('antistd%d' % l, [(l, 1), (l, 2), (l, 3), (o, 1), (o, 2)], [[0, 1, 2], [3, 4]],
                          {l: (1, 1), o: (1,)}))
    out.insert(0, Family('inv', [], [], {0: (), 1: ()}))
    out.insert(1, Family('mix', [(0, 1), (0, 2), (1, 1), (1, 2)], [[0, 1], [2, 3]], {0: (1,), 1: (1,)}))
    return out


def fam_by_name(name):
    return next(F for F in families() if F.name == name)


def hook_dim(u, tail):
    """dimension of the irrep (u - |tail|, tail...) of S_u (0 if not a partition)"""
    from math import factorial
    shape = [u - sum(tail)] + list(tail)
    if any(shape[i] < shape[i + 1] for i in range(len(shape) - 1)) or shape[0] < 0:
        return 0
    if u == 0:
        return 1
    hooks = 1
    for i, row in enumerate(shape):
        for j in range(row):
            arm = row - j - 1
            leg = sum(1 for i2 in range(i + 1, len(shape)) if shape[i2] > j)
            hooks *= arm + leg + 1
    return factorial(u) // hooks


# ---------------------------------------------------------------- the reduced forms of one family
def ginv(M, one):
    n = len(M)
    A = [list(r) + [one if i == j else one * 0 for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        piv = next(i for i in range(c, n) if not (A[i][c] == 0))
        A[c], A[piv] = A[piv], A[c]
        pinv = one / A[c][c]
        A[c] = [x * pinv for x in A[c]]
        for i in range(n):
            if i != c and not (A[i][c] == 0):
                f = A[i][c]
                rc = A[c]
                A[i] = [x - f * y for x, y in zip(A[i], rc)]
    return [r[n:] for r in A]


def mm(A, B, zero):
    n, m, p = len(A), len(B), len(B[0]) if B else 0
    out = [[zero] * p for _ in range(n)]
    for i in range(n):
        for l in range(m):
            a = A[i][l]
            if a == 0:
                continue
            Bl = B[l]
            row = out[i]
            for j in range(p):
                if not (Bl[j] == 0):
                    row[j] = row[j] + a * Bl[j]
    return out


def tr(A):
    return [list(r) for r in zip(*A)] if A else []


def gram_sym(F, ctx, vecs, factor):
    pats = sorted({p for (_, v) in vecs for p in v}, key=str)
    n = len(vecs)
    G = [[ctx.zero] * n for _ in range(n)]
    one = ctx.one
    for p in pats:
        a = nanon(p)
        cnt = ff(ctx.u[0] - F.nlev[0], a[0], one) * ff(ctx.u[1] - F.nlev[1], a[1], one)
        vals = [v.get(p, 0) for (_, v) in vecs]
        for i in range(n):
            if vals[i]:
                for j in range(n):
                    if vals[j]:
                        G[i][j] = G[i][j] + cnt * (vals[i] * vals[j])
    return [[x * factor for x in r] for r in G]


def choose(F, sel_u):
    """the bases: K-vectors, W-vectors (independent of K), Q-vectors, selected at the integer sizes sel_u"""
    if any(sel_u[l] < F.nlev[l] for l in (0, 1)):
        return [], [], []          # the tableau does not fit: the irrep does not occur at these sizes
    Kc, Wc, Qc = candidates(F)
    K = select(F, Kc, sel_u, Fr(1, 2))
    W = select(F, Wc, sel_u, Fr(1, 2), base=K)
    Q = select(F, Qc, sel_u, Fr(1, 6))
    return K, W, Q


def reduce_family(F, ctx, bases, b=1, full=True, log=None):
    """S' restricted to the W-vectors (and, if full, the whole reduced S' on K + W for the kernel check).
    Returns dict with GP, GQ, S (W-block), Sfull (if full), names."""
    import time
    t0 = time.time()
    K, W, Q = bases
    P = K + W
    nK, nW, nQ = len(K), len(W), len(Q)
    out = dict(names=([n for n, _ in K], [n for n, _ in W], [n for n, _ in Q]))
    if nW == 0:
        out['S'] = []
        return out
    h = ctx.one / 4
    GP = gram_sym(F, ctx, P, ctx.one / 2)
    out['GP'] = GP
    cache = {}
    accPP = joint_sum(F, P, P, None)
    RPP = assemble(F, ctx, accPP, oracle_pp(ctx, b), h, {})
    RA = assemble(F, ctx, accPP, oracle_a2(ctx), h, {})
    if log:
        log(f'    R_PP, R_A ({len(accPP)} joint classes) [{time.time() - t0:.1f}s]')
    GPi = ginv(GP, ctx.one)
    T1 = mm(RA if full else RA[nK:], GPi, ctx.zero)
    S1 = mm(mm(T1, RPP, ctx.zero), tr(T1), ctx.zero)
    if nQ:
        GQ = gram_sym(F, ctx, Q, ctx.one / 6)
        out['GQ'] = GQ
        accPQ = joint_sum(F, P, Q, None)
        RK = assemble(F, ctx, accPQ, oracle_k(ctx), ctx.one / 12, {})
        if b != 1:
            # with the B-part scaled, Y_PQ = A-part + b (Y_PQ - A-part); handled only numerically (see neg.py)
            raise ValueError('bscale != 1 is not supported in the congruence form')
        accQQ = joint_sum(F, Q, Q, None)
        RW = assemble(F, ctx, accQQ, oracle_w3(ctx), ctx.one / 36, {})
        GQi = ginv(GQ, ctx.one)
        T2 = mm(RK if full else RK[nK:], GQi, ctx.zero)
        S2 = mm(mm(T2, RW, ctx.zero), tr(T2), ctx.zero)
        S = [[S1[i][j] - ctx.half * S2[i][j] for j in range(len(S1))] for i in range(len(S1))]
    else:
        S = S1
    if log:
        log(f'    products done [{time.time() - t0:.1f}s]')
    if full:
        out['Sfull'] = S
        out['S'] = [row[nK:] for row in S[nK:]]
    else:
        out['S'] = S
    return out


# ---------------------------------------------------------------- the B-scaled block (negative controls, numeric only)
def oracle_apq(ctx):
    """the A-part of Y_PQ: -[r in T][X' = T - r] pi_T"""
    def f(key):
        sl, levs = key
        lev = dict(enumerate(levs))
        r, Z, T = sl[0], {sl[1], sl[2]}, {sl[3], sl[4], sl[5]}
        if r in T and Z == T - {r}:
            return -ctx.pi(T, lev)
        return ctx.zero
    return f


def oracle_dq(ctx):
    def f(key):
        sl, levs = key
        lev = dict(enumerate(levs))
        T1, T2 = {sl[0], sl[1], sl[2]}, {sl[3], sl[4], sl[5]}
        return ctx.pi(T1, lev) if T1 == T2 else ctx.zero
    return f


def reduce_scaled(F, ctx, bases, b):
    """reduced Schur complement S_b = Y_PP(b) - Y_PQ(b) Y_QQ(b)^{-1} Y_QP(b) on the full basis K + W, with every term
    of Lemma B that contains an E^{-1} multiplied by b (no congruence; products and inverses through the reduced
    bases).  Numeric fields only (it inverts the reduced A2 on the whole P-basis)."""
    K, W, Q = bases
    P = K + W
    one, zero = ctx.one, ctx.zero
    GP = gram_sym(F, ctx, P, one / 2)
    accPP = joint_sum(F, P, P, None)
    RPP = assemble(F, ctx, accPP, oracle_pp(ctx, b), one / 4, {})
    if not Q:
        return RPP
    RA = assemble(F, ctx, accPP, oracle_a2(ctx), one / 4, {})
    GQ = gram_sym(F, ctx, Q, one / 6)
    accPQ = joint_sum(F, P, Q, None)
    RK = assemble(F, ctx, accPQ, oracle_k(ctx), one / 12, {})
    RAPQ = assemble(F, ctx, accPQ, oracle_apq(ctx), one / 12, {})
    accQQ = joint_sum(F, Q, Q, None)
    RW = assemble(F, ctx, accQQ, oracle_w3(ctx), one / 36, {})
    RD = assemble(F, ctx, accQQ, oracle_dq(ctx), one / 36, {})
    RPQ1 = mm(GP, mm(ginv(RA, one), RK, zero), zero)
    RPQ = [[RAPQ[i][j] + b * (RPQ1[i][j] - RAPQ[i][j]) for j in range(len(Q))] for i in range(len(P))]
    RQQ1 = mm(mm(GQ, ginv(RW, one), zero), GQ, zero)
    RQQ = [[2 * RD[i][j] + b * (2 * RQQ1[i][j] - 2 * RD[i][j]) for j in range(len(Q))] for i in range(len(Q))]
    C = mm(mm(RPQ, ginv(RQQ, one), zero), tr(RPQ), zero)
    return [[RPP[i][j] - C[i][j] for j in range(len(P))] for i in range(len(P))]


def pd_exact(A):
    """True iff the symmetric matrix A (exact numbers) is positive definite (all pivots > 0)."""
    A = [list(r) for r in A]
    n = len(A)
    for c in range(n):
        p = A[c][c]
        if not p > 0:
            return False
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / p
                for j in range(c + 1, n):
                    A[i][j] -= f * A[c][j]
    return True
