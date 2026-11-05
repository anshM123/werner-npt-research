"""Exact (Fraction) evaluation of Lemma B of relative-anchors/README.md on the type I block of (R_inf)_{k-1}.

Lemma B (j = k - 1, Z = 1), unrescaled:
  2 Y[(a,y'),(b,y)] = <y'| (1/2 - [r=r']) s_b^+ s_a^- + [a=b](2 - n_r) |y> - <T_b y', E^{-1} T_a y>     (G1-G1)
  2 Y[(a,y'),(vac,y)] = -<y'| s_a^- |y> + <y'| J^- E^{-1} T_a |y>                                          (G1-G0)
  2 Y[(vac,y'),(vac,y)] = 2 [y = y']                                                                     (G0-G0)
with T_a = 1/2 J^+ s_a^- + [l=1] w_r (2|0><0|_r + |2><2|_r) - [l=2] w_r |1><2|_r and E = C + 1/2 J^+ J^- on the
column strings (block diagonal over 2-patterns and number of nonzero letters).
Returns Y (not 2Y) in the state order of fastdirect.typeI_states."""
from fractions import Fraction as Fr
import itertools
from fastdirect import typeI_states


def solve(Mx, b):
    """Exact Gaussian elimination, Mx square list of lists, b list.  Returns x with Mx x = b."""
    n = len(Mx)
    A = [list(Mx[i]) + [b[i]] for i in range(n)]
    for c in range(n):
        r = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[r] = A[r], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [A[i][n] for i in range(n)]


class LemmaB:
    def __init__(self, w):
        self.w = [Fr(x) for x in w]
        self.s = [x * x for x in self.w]
        self.k = len(w)

    # operators on column strings, as dicts {string: coeff}
    def Jp(self, vec):
        out = {}
        for y, c in vec.items():
            for t in range(self.k):
                if y[t] == 0:
                    z = y[:t] + (1,) + y[t + 1:]
                    out[z] = out.get(z, 0) + c * self.w[t]
        return out

    def Jm(self, vec):
        out = {}
        for y, c in vec.items():
            for t in range(self.k):
                if y[t] == 1:
                    z = y[:t] + (0,) + y[t + 1:]
                    out[z] = out.get(z, 0) + c * self.w[t]
        return out

    def sm(self, a, vec):   # sigma_a^- = |0><l|_r
        r, l = a
        out = {}
        for y, c in vec.items():
            if y[r] == l:
                z = y[:r] + (0,) + y[r + 1:]
                out[z] = out.get(z, 0) + c
        return out

    def sp(self, a, vec):   # sigma_a^+ = |l><0|_r
        r, l = a
        out = {}
        for y, c in vec.items():
            if y[r] == 0:
                z = y[:r] + (l,) + y[r + 1:]
                out[z] = out.get(z, 0) + c
        return out

    def Top(self, a, vec):
        r, l = a
        out = {}
        for z, c in self.Jp(self.sm(a, vec)).items():
            out[z] = out.get(z, 0) + c / 2
        for y, c in vec.items():
            if l == 1:
                f = 2 if y[r] == 0 else (1 if y[r] == 2 else 0)
                if f:
                    out[y] = out.get(y, 0) + c * self.w[r] * f
            else:
                if y[r] == 2:
                    z = y[:r] + (1,) + y[r + 1:]
                    out[z] = out.get(z, 0) - c * self.w[r]
        return {y: c for y, c in out.items() if c != 0}

    def Eapply_inv(self, vec):
        """E^{-1} vec, E = C + 1/2 J^+ J^-, solved on the span of all strings with the same 2-pattern and the same
        number of nonzero letters as those in vec (E preserves both)."""
        out = {}
        groups = {}
        for y, c in vec.items():
            key = (tuple(int(v == 2) for v in y), sum(v != 0 for v in y))
            groups.setdefault(key, {})[y] = c
        for (pat, nz), sub in groups.items():
            free = [t for t in range(self.k) if not pat[t]]
            n1 = nz - sum(pat)
            basis = []
            for ones in itertools.combinations(free, n1):
                basis.append(tuple(2 if pat[t] else (1 if t in ones else 0) for t in range(self.k)))
            idx = {y: i for i, y in enumerate(basis)}
            n = len(basis)
            Em = [[Fr(0)] * n for _ in range(n)]
            for j, y in enumerate(basis):
                col = {y: Fr(sum(self.s[t] * (2 if y[t] == 0 else (1 if y[t] == 2 else 0)) for t in range(self.k)))}
                jj = self.Jp(self.Jm({y: Fr(1)}))
                for z, c in jj.items():
                    col[z] = col.get(z, 0) + c / 2
                for z, c in col.items():
                    Em[idx[z]][j] += c
            rhs = [sub.get(y, Fr(0)) for y in basis]
            x = solve(Em, rhs)
            for y, c in zip(basis, x):
                if c != 0:
                    out[y] = out.get(y, 0) + c
        return out

    @staticmethod
    def dot(u, v):
        return sum((c * v[y] for y, c in u.items() if y in v), Fr(0))

    def entry2(self, X1, X2):
        """2 Y[X1, X2] by Lemma B."""
        (a, yp), (b, y) = X1, X2
        k = self.k
        ra = next((t for t in range(k) if a[t]), None)
        rb = next((t for t in range(k) if b[t]), None)
        if ra is None and rb is None:
            return Fr(2) if y == yp else Fr(0)
        if ra is None:   # symmetric
            return self.entry2(X2, X1)
        A_ = (ra, a[ra])
        if rb is None:
            v = -self.sm(A_, {y: Fr(1)}).get(yp, Fr(0))
            v += self.Jm(self.Eapply_inv(self.Top(A_, {y: Fr(1)}))).get(yp, Fr(0))
            return v
        B_ = (rb, b[rb])
        val = Fr(0)
        coef = Fr(1, 2) - (1 if ra == rb else 0)
        val += coef * self.sp(B_, self.sm(A_, {y: Fr(1)})).get(yp, Fr(0))
        if A_ == B_ and y == yp:
            val += 2 - (1 if y[ra] != 0 else 0)
        val -= self.dot(self.Top(B_, {yp: Fr(1)}), self.Eapply_inv(self.Top(A_, {y: Fr(1)})))
        return val

    def typeI(self, T):
        st, lab = typeI_states(self.k, T)
        Y = [[self.entry2(X1, X2) / 2 for X2 in st] for X1 in st]
        return Y, lab
