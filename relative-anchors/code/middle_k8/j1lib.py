"""Direct implementation of Lemma B (README of relative-anchors, Section 2) for the blocks of (R_inf)_1.
(J1.entry is Lemma B for any j: E^{-1} is applied on whatever sector the vector lies in; check_t3j2.py uses it for
j = 2, with units 2^{k-2}.)

Strings are tuples over {0, 1, 2}; vectors are dicts string -> number.  Weights w (Fractions or floats), s = w^2.
Lemma B, units 2^Z (Z = k - j):
  G1-G1:  <y'| (1/2 - [r = r']) s_b^+ s_a^- + [a = b](2 - n_r) - T_b^T E_j^{-1} T_a |y>
  G1-G0:  -<y'| s_a^- - J^- E_{j+1}^{-1} T_a |y>
  G0-G0:  2 <y'| I - 1/2 J^- E_{j+2}^{-1} J^+ |y>
States are pairs (a, y): a = None for the vacuum row, else a = (letter, copy).
block_states(k, T) lists the states of the torus block of (R_inf)_1 labelled by T (T = () is the letter-1 block).
"""
import itertools
from fractions import Fraction as Fr


def add(d, key, v):
    if v != 0:
        d[key] = d.get(key, 0) + v


def setl(y, t, l):
    z = list(y)
    z[t] = l
    return tuple(z)


class J1:
    def __init__(self, w, F=Fr, bscale=1):
        self.F = F
        self.bscale = F(bscale)   # multiplies the B-part (every term with E^{-1}); 1 for the true block
        self.w = [F(x) for x in w]
        self.k = len(w)
        self.s = [x * x for x in self.w]
        self.S = sum(self.s)
        self._Einv = {}

    # ---- operators on dict vectors ----
    def sig(self, a, v):
        l, r = a
        out = {}
        for y, c in v.items():
            if y[r] == l:
                add(out, setl(y, r, 0), c)
        return out

    def sigp(self, a, v):
        l, r = a
        out = {}
        for y, c in v.items():
            if y[r] == 0:
                add(out, setl(y, r, l), c)
        return out

    def Jp(self, v):
        out = {}
        for y, c in v.items():
            for t in range(self.k):
                if y[t] == 0:
                    add(out, setl(y, t, 1), self.w[t] * c)
        return out

    def Jm(self, v):
        out = {}
        for y, c in v.items():
            for t in range(self.k):
                if y[t] == 1:
                    add(out, setl(y, t, 0), self.w[t] * c)
        return out

    def Tm(self, a, v):
        l, r = a
        out = {}
        for y, c in self.Jp(self.sig(a, v)).items():
            add(out, y, c / 2)
        for y, c in v.items():
            if l == 1:
                if y[r] == 0:
                    add(out, y, 2 * self.w[r] * c)
                elif y[r] == 2:
                    add(out, y, self.w[r] * c)
            else:
                if y[r] == 2:
                    add(out, setl(y, r, 1), -self.w[r] * c)
        return out

    def Cd(self, y):
        return sum(self.s[t] * (2 if y[t] == 0 else (1 if y[t] == 2 else 0)) for t in range(self.k))

    def Eop(self, v):
        out = {}
        for y, c in v.items():
            add(out, y, self.Cd(y) * c)
        for y, c in self.Jp(self.Jm(v)).items():
            add(out, y, c / 2)
        return out

    def sector(self, y):
        """all strings with the same number of nonzero letters, same letter-2 positions, same number of 1s."""
        k = self.k
        twos = tuple(t for t in range(k) if y[t] == 2)
        n1 = sum(1 for t in range(k) if y[t] == 1)
        free = [t for t in range(k) if y[t] != 2]
        out = []
        for P in itertools.combinations(free, n1):
            z = [0] * k
            for t in twos:
                z[t] = 2
            for t in P:
                z[t] = 1
            out.append(tuple(z))
        return out

    def Einv(self, v):
        """E^{-1} v, sector by sector (E preserves the sector)."""
        bysec = {}
        for y, c in v.items():
            key = tuple(sorted(self.sector(y)))
            bysec.setdefault(key, {})[y] = c
        out = {}
        for key, vv in bysec.items():
            if key not in self._Einv:
                strs = list(key)
                idx = {y: i for i, y in enumerate(strs)}
                n = len(strs)
                M = [[self.F(0)] * n for _ in range(n)]
                for j, y in enumerate(strs):
                    for z, c in self.Eop({y: self.F(1)}).items():
                        M[idx[z]][j] += c
                self._Einv[key] = (strs, idx, inv(M, self.F))
            strs, idx, Mi = self._Einv[key]
            for i, y in enumerate(strs):
                val = sum((Mi[i][idx[z]] * c for z, c in vv.items()), self.F(0))
                add(out, y, val)
        return out

    # ---- entries of 2^Z Y ----
    def entry(self, X1, X2):
        (a, yp), (b, y) = X1, X2
        F = self.F
        if a is not None and b is not None:
            val = F(0)
            ra, rb = a[1], b[1]
            coef = F(1) / 2 - (1 if ra == rb else 0)
            u1 = self.sig(b, {yp: F(1)})
            u2 = self.sig(a, {y: F(1)})
            val += coef * dot(u1, u2)
            if a == b and yp == y:
                val += 2 - (1 if y[ra] != 0 else 0)
            t1 = self.Tm(b, {yp: F(1)})
            t2 = self.Einv(self.Tm(a, {y: F(1)}))
            val -= self.bscale * dot(t1, t2)
            return val
        if a is not None and b is None:
            # -<yp| s_a^- - J^- E^{-1} T_a |y>
            v1 = self.sig(a, {y: F(1)})
            v2 = self.Jm(self.Einv(self.Tm(a, {y: F(1)})))
            return -(v1.get(yp, 0) - self.bscale * v2.get(yp, 0))
        if a is None and b is not None:
            return self.entry(X2, X1)
        # both vacuum
        v = self.Jm(self.Einv(self.Jp({y: F(1)})))
        return 2 * ((1 if yp == y else 0) - self.bscale * v.get(yp, 0) / 2)

    def matrix(self, states):
        n = len(states)
        M = [[None] * n for _ in range(n)]
        for i in range(n):
            for j in range(i, n):
                M[i][j] = M[j][i] = self.entry(states[i], states[j])
        return M


def dot(u, v):
    return sum((c * v[y] for y, c in u.items() if y in v), 0)


def inv(M, F):
    n = len(M)
    A = [list(row) + [F(1) if i == j else F(0) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        r = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[r] = A[r], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * yv for x, yv in zip(A[i], A[c])]
    return [row[n:] for row in A]


def onehot(k, pos):
    """string with letters given by dict pos {copy: letter}"""
    z = [0] * k
    for t, l in pos.items():
        z[t] = l
    return tuple(z)


def block_states(k, T):
    """torus block of (R_inf)_1 labelled by T: T = () letter-1 block; T = (t,) ; T = (t, t2)."""
    T = tuple(T)
    if len(T) == 0:
        st = [((1, r), onehot(k, {c: 1})) for r in range(k) for c in range(k)]
        st += [(None, onehot(k, {c: 1, c2: 1})) for c, c2 in itertools.combinations(range(k), 2)]
        return st
    if len(T) == 1:
        t = T[0]
        st = [((2, t), onehot(k, {c: 1})) for c in range(k)]
        st += [((1, r), onehot(k, {t: 2})) for r in range(k)]
        st += [(None, onehot(k, {t: 2, c: 1})) for c in range(k) if c != t]
        return st
    if len(T) == 2:
        t, t2 = T
        return [((2, t), onehot(k, {t2: 2})), ((2, t2), onehot(k, {t: 2})), (None, onehot(k, {t: 2, t2: 2}))]
    raise ValueError


def mu(state, w):
    """rescaling factor of ylib: prod_t w_t^{N_t}, N_t = number of letters 1 at copy t in row and column."""
    a, y = state
    k = len(y)
    m = 1
    for t in range(k):
        n = (1 if (a is not None and a[0] == 1 and a[1] == t) else 0) + (1 if y[t] == 1 else 0)
        m *= w[t] ** n
    return m


def to_ylib_state(state, k):
    a, y = state
    if a is None:
        return ((0,) * k, y)
    return (onehot(k, {a[1]: a[0]}), y)
