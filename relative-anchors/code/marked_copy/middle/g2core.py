"""Slices of the determinant of the gauge-fixed generic two-row form (g2lib.build): the 4-variable entries are
specialised in two of the variables and the determinant is computed exactly over the rational functions in the
other two (python-flint fmpq_mpoly with reduced fractions, flintfield)."""
import flint
from flintfield import Field, F
import g2lib

VAR = {'q0': 0, 't': 1, 'i': 2, 's': 3}


class Slicer:
    def __init__(self, check=False, log=print):
        self.C, self.keep, self.M = g2lib.build(check=check, log=log)
        self.n = len(self.keep)
        facs = {}
        for a in self.keep:
            for b in self.keep:
                e = self.M[a, b]
                if e == 0:
                    continue
                c, fl = e.d.factor()
                for f, m in fl:
                    facs[str(f)] = f
        self.factors = list(facs.values())
        self.fields = {}

    def field(self, free):
        if free not in self.fields:
            self.fields[free] = Field(list(free))
        return self.fields[free]

    def spec(self, P, free, vals):
        """substitute the fixed variables (vals: name -> int) into the 4-variable polynomial P; the free variables
        (a pair of names, in this order) become the generators of a 2-variable context"""
        K = self.field(free)
        gens = dict(zip(free, K.ctx.gens()))
        args = [gens[v] if v in gens else K.ctx.from_dict({(0,) * len(free): vals[v]}) for v in ('q0', 't', 'i', 's')]
        return P.compose(*args, ctx=K.ctx)

    def det(self, free, vals):
        """det of the gauge-fixed form with the variables in vals fixed, as a reduced fraction over Q(free)"""
        K = self.field(free)
        A = [[(F(K, self.spec(self.M[a, b].n, free, vals), self.spec(self.M[a, b].d, free, vals))
               if self.M[a, b] != 0 else K(0)) for b in self.keep] for a in self.keep]
        n = self.n
        det = K(1)
        for c in range(n):
            p = A[c][c]
            if p == 0:
                # exchange with a later row/column (symmetric permutation keeps the determinant)
                r = next(r for r in range(c + 1, n) if A[r][r] != 0)
                A[c], A[r] = A[r], A[c]
                for row in A:
                    row[c], row[r] = row[r], row[c]
                p = A[c][c]
            det = det * p
            for r in range(c + 1, n):
                if A[r][c] != 0:
                    f = A[r][c] / p
                    for cc in range(c + 1, n):
                        A[r][cc] = A[r][cc] - f * A[c][cc]
        return det


def multiplicity(f, P):
    """the largest e with f^e | P (f non-constant)"""
    e = 0
    while True:
        q, r = divmod(P, f)
        if r != 0:
            return e
        P = q
        e += 1
