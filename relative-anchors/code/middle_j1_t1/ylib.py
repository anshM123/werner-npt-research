"""Builder of the (R_inf)_j blocks in the squared-weight rescaled basis (squared weights s_t = w_t^2, rescaled by the diagonal factors mu).
States X = (a, y): row string a (vacuum or one nonzero letter), column string y.  Letters 0, 1, 2.
Rescaled entries (mu_X = prod_t w_t^{N_t(X)}, N_t = number of letters 1 at copy t in row and column of X):
  A'[(a,y'),(b,y)] = M[(b,y'),(a,y)] * s^{(N(a)+N(y')+N(b)+N(y))/2}
  B'[(a,y'),(b,y)] = R'[(b,y'),(a,y)],   R' = MV' D'^{-1} MV'^T
  MV'[(b,y), z]    = sum_t M[(b,y),(e_t,z)] s^{(N(b)+N(y)+N(z)+e_t)/2}
  D'[y,y']         = sum_{t,t'} M[(e_t,y),(e_t',y')] s^{(N(y)+N(y')+e_t+e_t')/2}
Y' = A' - B'.  The ring is either the flint rational-function field (symbolic s) or flint.fmpq (rational point).
"""
import itertools
import flint
import sys
from flintfield import Field

HALF = flint.fmpq(1, 2)


def m1(b, yp, a, y):
    return (1 if (a == b and y == yp) else 0) - (HALF if (b == yp and a == y) else 0)


def Mel(b, yp, a, y):
    v = flint.fmpq(1)
    for s in range(len(b)):
        v *= m1(b[s], yp[s], a[s], y[s])
        if v == 0:
            return v
    return v


class Ring:
    """s: list of ring elements (the squared weights)."""

    def __init__(self, s, zero, one):
        self.s, self.zero, self.one = s, zero, one
        self.k = len(s)
        self._mon = {}

    def mono(self, ex):
        if ex not in self._mon:
            v = self.one
            for u, x in enumerate(ex):
                for _ in range(x):
                    v = v * self.s[u]
            self._mon[ex] = v
        return self._mon[ex]


def symbolic_ring(k, names=None):
    names = names or [f's{t}' for t in range(k)]
    K = Field(names)
    return Ring(K.gens, K(0), K(1)), K


def point_ring(svals):
    return Ring([flint.fmpq(x) if not isinstance(x, flint.fmpq) else x for x in svals], flint.fmpq(0), flint.fmpq(1))


def Nv(*strs):
    k = len(strs[0])
    return tuple(sum(int(st[u] == 1) for st in strs) for u in range(k))


def halfex(*vecs):
    k = len(vecs[0])
    tot = [sum(v[u] for v in vecs) for u in range(k)]
    assert all(x % 2 == 0 for x in tot), tot
    return tuple(x // 2 for x in tot)


def inv(Mx, zero, one):
    n = len(Mx)
    A = [row[:] + [one if i == j else zero for j in range(n)] for i, row in enumerate(Mx)]
    for c in range(n):
        r = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[r] = A[r], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]


def wcol(y):
    return (sum(v == 1 for v in y),) + tuple(int(v == 2) for v in y)


class Builder:
    def __init__(self, ring, letters=(0, 1, 2)):
        self.R = ring
        k = self.k = ring.k
        self.e = [tuple(1 if r == t else 0 for r in range(k)) for t in range(k)]
        self.cols = list(itertools.product(letters, repeat=k))
        self._D = None
        self._MV = {}
        self._T = {}
        self._Rc = {}

    def Dinv(self):
        if self._D is None:
            R, k, e = self.R, self.k, self.e
            blocks = {}
            for y in self.cols:
                blocks.setdefault(wcol(y), []).append(y)
            Dinv = {}
            for key, ys in blocks.items():
                Bm = []
                for y in ys:
                    row = []
                    for yp in ys:
                        v = R.zero
                        for t in range(k):
                            for tp in range(k):
                                mv = Mel(e[t], y, e[tp], yp)
                                if mv:
                                    v += R.one * mv * R.mono(halfex(Nv(y), Nv(yp), e[t], e[tp]))
                        row.append(v)
                    Bm.append(row)
                Bi = inv(Bm, R.zero, R.one)
                for i, y in enumerate(ys):
                    for j, yp in enumerate(ys):
                        if Bi[i][j] != 0:
                            Dinv.setdefault(y, {})[yp] = Bi[i][j]
            self._D = Dinv
        return self._D

    def MV(self, b, y):
        key = (b, y)
        if key not in self._MV:
            R, k, e = self.R, self.k, self.e
            out = {}
            for z in self.cols:
                v = R.zero
                for t in range(k):
                    mv = Mel(b, y, e[t], z)
                    if mv:
                        v += R.one * mv * R.mono(halfex(Nv(b, y), Nv(z), e[t]))
                if v != 0:
                    out[z] = v
            self._MV[key] = out
        return self._MV[key]

    def T(self, b, y):
        key = (b, y)
        if key not in self._T:
            Dinv = self.Dinv()
            acc = {}
            for z, v in self.MV(b, y).items():
                for zp, dv in Dinv.get(z, {}).items():
                    acc[zp] = acc.get(zp, self.R.zero) + v * dv
            self._T[key] = acc
        return self._T[key]

    def Rel(self, b, yp, a, y):
        key = (b, yp, a, y)
        if key not in self._Rc:
            r1 = self.T(b, yp)
            r2 = self.MV(a, y)
            s = self.R.zero
            for z, v in r1.items():
                v2 = r2.get(z)
                if v2 is not None:
                    s += v * v2
            self._Rc[key] = s
        return self._Rc[key]

    def A(self, a, yp, b, y):
        mv = Mel(b, yp, a, y)
        return self.R.one * mv * self.R.mono(halfex(Nv(a, yp), Nv(b, y))) if mv else self.R.zero

    def B(self, a, yp, b, y):
        return self.Rel(b, yp, a, y)

    def Y(self, X1, X2):
        (a, yp), (b, y) = X1, X2
        return self.A(a, yp, b, y) - self.B(a, yp, b, y)

    def matrix(self, states):
        return [[self.Y(X1, X2) for X2 in states] for X1 in states]


def q(st):
    return sum(x != 0 for x in st)


def rinf_states(k, j, letters=(0, 1, 2)):
    allst = list(itertools.product(letters, repeat=k))
    H = {i: [st for st in allst if q(st) == i] for i in range(k + 1)}
    vac = (0,) * k
    S = ([(a, y) for a in H[1] for y in H[j]] if 0 <= j <= k else []) + \
        ([(vac, y) for y in H[j + 1]] if j + 1 <= k else [])
    return S


def wstate(a, y):
    k = len(a)
    return (sum(v == 1 for v in a) + sum(v == 1 for v in y),) + tuple(int(a[s] == 2) + int(y[s] == 2) for s in range(k))


def ldl(Mx, order=None, zero=None):
    """symmetric elimination without pivoting in the given order; returns list of (index, pivot) and the final
    Schur complement on the remaining (zero-pivot) indices.  Stops pivoting on an index with zero diagonal
    (it is skipped and kept)."""
    n = len(Mx)
    M = [row[:] for row in Mx]
    order = list(range(n)) if order is None else list(order)
    pivs = []
    act = list(order)
    skipped = []
    while act:
        i = act[0]
        if M[i][i] == 0:
            skipped.append(i)
            act = act[1:]
            continue
        piv = M[i][i]
        pivs.append((i, piv))
        rest = act[1:] + []
        allrest = rest + skipped
        for jj, j in enumerate(allrest):
            if M[j][i] == 0:
                continue
            f = M[j][i] / piv
            for l in allrest[jj:]:
                if M[i][l] != 0:
                    M[j][l] = M[j][l] - f * M[i][l]
                    if l != j:
                        M[l][j] = M[j][l]
        act = rest
    rem = [[M[i][l] for l in skipped] for i in skipped]
    return pivs, skipped, rem
