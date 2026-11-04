"""Equal-weight (R_inf)_j, LETTER-2 torus blocks, reduced by S_T x S_U symmetry (explicit vectors, numeric k).
Lemmas B and N of relative-anchors/README.md at w = 1, letters 0, 1, 2; scaled by 2^(k-j).
Case I pattern: T = copies with one letter 2 (last tau copies), U = the others (pairs + rest R, as in eqlib).
  T-sym vectors (S_T-invariant): A_x, B_xy (row 1 in U), G0, C1 (row 1 in T), D_z (row 2 at t in T, column z at t).
  T-anti vectors (sign on (t0, t1), symmetric on the rest of T): C1a, D_za.
Case II pattern: copy t* with letter 2 in row and column, T1 = copies with a column 2: one vector per i.
"""
import flint
from itertools import combinations

fq = flint.fmpq
HALF = fq(1, 2)


def addv(acc, v, c=1):
    for x, val in v.items():
        nv = acc.get(x, 0) + c * val
        if nv == 0:
            acc.pop(x, None)
        else:
            acc[x] = nv
    return acc


def setl(s, r, l):
    return s[:r] + (l,) + s[r + 1:]


def U(v):
    out = {}
    for s, c in v.items():
        for t, l in enumerate(s):
            if l == 0:
                y = setl(s, t, 1)
                out[y] = out.get(y, 0) + c
    return {a: b for a, b in out.items() if b != 0}


def D(v):
    out = {}
    for s, c in v.items():
        for t, l in enumerate(s):
            if l == 1:
                y = setl(s, t, 0)
                out[y] = out.get(y, 0) + c
    return {a: b for a, b in out.items() if b != 0}


def sm(r, l, v):
    return {setl(s, r, 0): c for s, c in v.items() if s[r] == l}


def nop(r, v):
    return {s: c for s, c in v.items() if s[r] != 0}


def Nop(r, v):
    f = {0: fq(3, 2), 1: HALF, 2: fq(1)}
    return {s: f[s[r]] * c for s, c in v.items()}


def T(a, v):
    r, l = a
    out = {s: HALF * c for s, c in sm(r, l, U(v)).items()}
    if l == 1:
        addv(out, Nop(r, v))
    else:
        addv(out, {setl(s, r, 1): -HALF * c for s, c in v.items() if s[r] == 2})
    return out


def Einv(v):
    """E = calN + 1/2 DU, block diagonal in (2-pattern, number of ones)."""
    groups = {}
    for s, c in v.items():
        key = (tuple(x == 2 for x in s), sum(x == 1 for x in s))
        groups.setdefault(key, {})[s] = c
    out = {}
    for (pat, c1), w0 in groups.items():
        kB = sum(1 for x in pat if not x)
        n2 = len(pat) - kB
        n0 = kB - c1
        alpha = fq(3, 2) * n0 + HALF * c1 + n2
        lam = [fq((kB - c1 - ip) * (c1 + 1 - ip)) for ip in range(0, min(c1, kB - c1) + 1)]
        w = dict(w0)
        den = fq(1)
        for n in range(len(lam)):
            den = den * (alpha + lam[n] / 2)
            addv(out, w, fq((-1) ** n, 2 ** n) / den)
            if n < len(lam) - 1:
                w = addv(D(U(w)), w, -lam[n])
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


def boolpart(k, pairs, R, fixed, q):
    """dict tuple->coef over sites of U: singlets on pairs not in fixed, letters in fixed, Dicke(q) on R - fixed"""
    if q < 0:
        return {}
    base = {(): 1}
    out = {}
    freeR = [t for t in R if t not in fixed]
    if q > len(freeR):
        return {}
    sing = [p for p in pairs if p[0] not in fixed]
    for ones in combinations(freeR, q):
        for bits in range(2 ** len(sing)):
            s = [0] * k
            c = 1
            for t, l in fixed.items():
                s[t] = l
            for n, (p, qq) in enumerate(sing):
                if (bits >> n) & 1:
                    s[qq] = 1
                    c = -c
                else:
                    s[p] = 1
            for t in ones:
                s[t] = 1
            out[tuple(s)] = out.get(tuple(s), 0) + c
    return out


def withT(v, Tset, special=None):
    """put letter 2 on Tset (except site special=(t, z))"""
    out = {}
    for s, c in v.items():
        s = list(s)
        for t in Tset:
            s[t] = 2
        if special is not None:
            s[special[0]] = special[1]
        out[tuple(s)] = fq(c)
    return out


def caseI_vectors(k, j, tau, i, mode):
    u = k - tau
    Tset = list(range(u, k))
    pairs = [(2 * l, 2 * l + 1) for l in range(i)]
    R = list(range(2 * i, u))
    vecs = []
    if mode == 'sym':
        for x in (0, 1):
            psi = {}
            for a in R:
                v = withT(boolpart(k, pairs, R, {a: x}, j - tau - i - x), Tset)
                if v:
                    psi[(a, 1)] = v
            if psi:
                vecs.append((f'A{x}', psi, {}))
        for x in (0, 1):
            for y in (0, 1):
                psi = {}
                for l, (p, qq) in enumerate(pairs):
                    for a, ap, sg in ((p, qq, 1), (qq, p, -1)):
                        others = [pp for pp in pairs if pp != (p, qq)]
                        v = boolpart(k, others, R, {a: x, ap: y}, j - tau - (i - 1) - x - y)
                        v = withT({s: sg * c for s, c in v.items()}, Tset)
                        if v:
                            psi[(a, 1)] = v
                if psi:
                    vecs.append((f'B{x}{y}', psi, {}))
        phi = withT(boolpart(k, pairs, R, {}, j + 1 - tau - i), Tset)
        if phi:
            vecs.append(('G0', {}, phi))
        psi = {}
        for r in Tset:
            v = withT(boolpart(k, pairs, R, {}, j - tau - i), Tset)
            if v:
                psi[(r, 1)] = v
        if psi:
            vecs.append(('C1', psi, {}))
        for z in (0, 1):
            psi = {}
            for t in Tset:
                v = withT(boolpart(k, pairs, R, {}, j - tau + 1 - z - i), Tset, special=(t, z))
                if v:
                    psi[(t, 2)] = v
            if psi:
                vecs.append((f'D{z}', psi, {}))
    else:
        t0, t1 = Tset[0], Tset[1]
        psi = {}
        for r, sg in ((t0, 1), (t1, -1)):
            v = withT(boolpart(k, pairs, R, {}, j - tau - i), Tset)
            if v:
                psi[(r, 1)] = {s: sg * c for s, c in v.items()}
        if psi:
            vecs.append(('C1a', psi, {}))
        for z in (0, 1):
            psi = {}
            for t, sg in ((t0, 1), (t1, -1)):
                v = withT(boolpart(k, pairs, R, {}, j - tau + 1 - z - i), Tset, special=(t, z))
                if v:
                    psi[(t, 2)] = {s: sg * c for s, c in v.items()}
            if psi:
                vecs.append((f'D{z}a', psi, {}))
    return vecs


def caseII_vectors(k, j, tau1, i):
    """t* = k-1 (row 2, column 2), T1 = the tau1 copies before it, U = the rest"""
    ts = k - 1
    T1 = list(range(k - 1 - tau1, k - 1))
    u = k - 1 - tau1
    pairs = [(2 * l, 2 * l + 1) for l in range(i)]
    R = list(range(2 * i, u))
    v = withT(boolpart(k, pairs, R, {}, j - 1 - tau1 - i), T1 + [ts])
    return [('II', {(ts, 2): v}, {})] if v else []


def row_reps(k, tau, i, mode):
    """H-orbit representatives of rows (r, l) with orbit sizes (H = S_T-part x (S_2 wr S_i) x S_R)"""
    u = k - tau
    reps = []
    m = u - 2 * i
    Tset = list(range(u, k))
    for l in (1, 2):
        if m >= 1:
            reps.append(((2 * i, l), m))
        if i >= 1:
            reps.append(((0, l), 2 * i))
        if tau >= 1:
            if mode == 'anti':
                reps.append(((Tset[0], l), 2))
                if tau >= 3:
                    reps.append(((Tset[2], l), tau - 2))
            else:
                reps.append(((Tset[0], l), tau))
    return reps


def form(k, j, vecs, reps=None):
    rows = [(r, l) for r in range(k) for l in (1, 2)]
    if reps is None:
        reps = [(a, 1) for a in rows]
    n = len(vecs)
    G = [[fq(0)] * n for _ in range(n)]
    Q = [[fq(0)] * n for _ in range(n)]
    Ecache = {}
    for be_, (_, psib, phib) in enumerate(vecs):
        for a, _w in reps:
            for b, pb in psib.items():
                Ecache[(be_, a, b)] = Einv(T(a, pb))
    Ccache = {}
    for be_, (_, psib, phib) in enumerate(vecs):
        if phib:
            for a, _w in reps:
                v = {s: -c for s, c in sm(a[0], a[1], phib).items()}
                addv(v, D(Einv(T(a, phib))))
                Ccache[(be_, a)] = v
    for al in range(n):
        _, psia, phia = vecs[al]
        for be_ in range(n):
            _, psib, phib = vecs[be_]
            g = fq(0)
            q = fq(0)
            for a, wa in reps:
                pa = psia.get(a)
                if pa:
                    pb = psib.get(a)
                    if pb:
                        g += wa * dot(pa, pb)
                        q += wa * (2 * dot(pa, pb) - dot(pa, nop(a[0], pb)))
                    for b, pb2 in psib.items():
                        coef = HALF - (1 if a[0] == b[0] else 0)
                        q += wa * coef * dot(sm(b[0], b[1], pa), sm(a[0], a[1], pb2))
                        q -= wa * dot(T(b, pa), Ecache[(be_, a, b)])
                    if phib:
                        q += wa * dot(pa, Ccache[(be_, a)])
                pb = psib.get(a)
                if pb and phia:
                    q += wa * dot(pb, Ccache[(al, a)])
            if phia and phib:
                g += dot(phia, phib)
                q += 2 * dot(phia, phib)
                Ua, Ub = U(phia), U(phib)
                if Ua and Ub:
                    q -= dot(Ua, Einv(Ub))
            G[al][be_] = g
            Q[al][be_] = q
    return G, Q
