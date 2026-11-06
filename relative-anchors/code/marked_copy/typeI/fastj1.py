"""A fast exact builder of the blocks of Lemma B (README of relative-anchors, Section 2), for any list of states of one
(R_inf)_j (any letters), with the same output as j1lib.J1.matrix (2^Z Y, Z = k - j).

j1lib.J1.entry computes every entry separately, applying E^{-1} each time.  Here the E^{-1} terms are grouped:
  G1-G1   <y'| T_b^T E^{-1} T_a |y>  = (T_b^T E^{-1} T_a)[y', y]     for X1 = (a, y'), X2 = (b, y),
  G1-G0   <y'| J^- E^{-1} T_a |y>    = ((J^+)^T E^{-1} T_a)[y', y]   for X1 = (a, y'), X2 = (vac, y),
  G0-G0   <y'| J^- E^{-1} J^+ |y>    = ((J^+)^T E^{-1} J^+)[y', y],
(J^- is the transpose of J^+ and E is symmetric), with E^{-1} applied once per row letter a to the matrix of the
vectors T_a|y>, sector by sector (E preserves the positions of the letters 2 and the number of letters 1), in exact
rational arithmetic (python-flint).  The operators sig, Jp, Tm and Eop are those of j1lib.J1.
"""
from fractions import Fraction as Fr
import flint
from j1lib import J1


def fq(x):
    x = Fr(x)
    return flint.fmpq(x.numerator, x.denominator)


def tofr(x):
    return Fr(int(x.p), int(x.q))


class FastJ1:
    def __init__(self, w, bscale=1):
        self.J = J1([Fr(x) for x in w], bscale=bscale)
        self.k = len(w)
        self.bscale = Fr(bscale)
        self._E = {}

    def seckey(self, y):
        return (tuple(t for t in range(self.k) if y[t] == 2), sum(1 for t in range(self.k) if y[t] == 1))

    def sector(self, key):
        if key not in self._E:
            any_y = None
            twos, n1 = key
            y = [0] * self.k
            for t in twos:
                y[t] = 2
            free = [t for t in range(self.k) if t not in twos]
            for t in free[:n1]:
                y[t] = 1
            strs = sorted(self.J.sector(tuple(y)))
            idx = {s: i for i, s in enumerate(strs)}
            n = len(strs)
            ent = [[flint.fmpq(0)] * n for _ in range(n)]
            for jj, s in enumerate(strs):
                for z, c in self.J.Eop({s: Fr(1)}).items():
                    ent[idx[z]][jj] += fq(c)
            self._E[key] = (strs, idx, flint.fmpq_mat(n, n, [x for row in ent for x in row]))
        return self._E[key]

    def gram(self, P, Q):
        """dict (p, q) -> <P_p, E^{-1} Q_q> for dicts of named vectors P, Q (sector by sector)"""
        bysec = {}
        for side, vecs in ((0, P), (1, Q)):
            for name, v in vecs.items():
                for z, c in v.items():
                    if c != 0:
                        bysec.setdefault(self.seckey(z), ({}, {}))[side].setdefault(name, {})[z] = c
        out = {}
        for key, (Pp, Qq) in bysec.items():
            if not Pp or not Qq:
                continue
            strs, idx, E = self.sector(key)
            n = len(strs)
            pn, qn = list(Pp), list(Qq)
            Pm = flint.fmpq_mat(n, len(pn))
            for c, name in enumerate(pn):
                for z, x in Pp[name].items():
                    Pm[idx[z], c] = fq(x)
            Qm = flint.fmpq_mat(n, len(qn))
            for c, name in enumerate(qn):
                for z, x in Qq[name].items():
                    Qm[idx[z], c] = fq(x)
            R = Pm.transpose() * E.solve(Qm)
            for a, pa in enumerate(pn):
                for b, qb in enumerate(qn):
                    v = R[a, b]
                    if v != 0:
                        out[(pa, qb)] = out.get((pa, qb), 0) + tofr(v)
        return out

    def matrix(self, states):
        J = self.J
        n = len(states)
        one = Fr(1)
        rows = sorted({a for a, _ in states if a is not None})
        cols1 = sorted({y for a, y in states if a is not None})
        cols0 = sorted({y for a, y in states if a is None})
        # vectors
        T1 = {(a, y): J.Tm(a, {y: one}) for a in rows for y in cols1}          # T_a |y>, y a column of a G1 state
        T0 = {(a, y): J.Tm(a, {y: one}) for a in rows for y in cols0}          # T_a |y>, y a column of a G0 state
        P1 = {y: J.Jp({y: one}) for y in cols1}                               # J^+ |y'>, y' a column of a G1 state
        P0 = {y: J.Jp({y: one}) for y in cols0}
        G11 = self.gram(T1, T1)
        G10 = self.gram(P1, T0)
        G00 = self.gram(P0, P0)
        bs = self.bscale
        M = [[None] * n for _ in range(n)]
        for p in range(n):
            a, yp = states[p]
            for q in range(p, n):
                b, y = states[q]
                if a is not None and b is not None:
                    ra, rb = a[1], b[1]
                    coef = Fr(1, 2) - (1 if ra == rb else 0)
                    val = Fr(0)
                    if yp[rb] == b[0] and y[ra] == a[0]:
                        z1 = list(yp); z1[rb] = 0
                        z2 = list(y); z2[ra] = 0
                        if z1 == z2:
                            val += coef
                    if a == b and yp == y:
                        val += 2 - (1 if y[ra] != 0 else 0)
                    val -= bs * G11.get(((b, yp), (a, y)), 0)
                elif a is not None and b is None:
                    sg = 0
                    if y[a[1]] == a[0]:
                        z = list(y); z[a[1]] = 0
                        if tuple(z) == yp:
                            sg = 1
                    val = -(sg - bs * G10.get((yp, (a, y)), 0))
                elif a is None and b is not None:
                    sg = 0
                    if yp[b[1]] == b[0]:
                        z = list(yp); z[b[1]] = 0
                        if tuple(z) == y:
                            sg = 1
                    val = -(sg - bs * G10.get((y, (b, yp)), 0))
                else:
                    val = 2 * (1 if yp == y else 0) - bs * G00.get((yp, y), 0)
                M[p][q] = M[q][p] = Fr(val)
        return M
