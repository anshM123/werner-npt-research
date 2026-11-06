"""Exact univariate rational functions over Q, on top of python-flint fmpq_poly (fast gcd).
Elements num/den, reduced, with monic den.  Also: evaluation at a rational point, Sturm counts of real roots of a
polynomial in an open interval (exact), and a sign test for 'positive for every q > 0'."""
from fractions import Fraction as Fr
import flint

P = flint.fmpq_poly


def _const(c):
    if isinstance(c, Fr):
        return P([flint.fmpq(c.numerator, c.denominator)])
    if isinstance(c, flint.fmpq):
        return P([c])
    return P([c])


class RF:
    __slots__ = ('n', 'd')

    def __init__(self, n, d=None, reduce=True):
        if not isinstance(n, P):
            n = _const(n)
        if d is None:
            d = P([1])
            reduce = False
        if n.is_zero():
            self.n, self.d = n, P([1])
            return
        if reduce:
            g = n.gcd(d)
            if not g.is_one():
                n = n // g
                d = d // g
        lc = d.leading_coefficient()
        if lc != 1:
            n = n / lc
            d = d / lc
        self.n, self.d = n, d

    @staticmethod
    def _c(o):
        return o if isinstance(o, RF) else RF(o)

    def __add__(self, o):
        o = RF._c(o)
        if self.n.is_zero():
            return o
        if o.n.is_zero():
            return self
        if self.d == o.d:
            return RF(self.n + o.n, self.d)
        g = self.d.gcd(o.d)
        if g.is_one():
            return RF(self.n * o.d + o.n * self.d, self.d * o.d, reduce=False)
        d1, d2 = self.d // g, o.d // g
        return RF(self.n * d2 + o.n * d1, self.d * d2)

    __radd__ = __add__

    def __neg__(self):
        r = RF.__new__(RF)
        r.n, r.d = -self.n, self.d
        return r

    def __sub__(self, o):
        return self + (-RF._c(o))

    def __rsub__(self, o):
        return RF._c(o) + (-self)

    def __mul__(self, o):
        if not isinstance(o, RF):
            if o == 0 or self.n.is_zero():
                return RF(0)
            r = RF.__new__(RF)
            r.n, r.d = self.n * _const(o), self.d
            return r
        if self.n.is_zero() or o.n.is_zero():
            return RF(0)
        g1 = self.n.gcd(o.d)
        g2 = o.n.gcd(self.d)
        n = (self.n // g1) * (o.n // g2)
        d = (self.d // g2) * (o.d // g1)
        return RF(n, d, reduce=False)

    __rmul__ = __mul__

    def inv(self):
        if self.n.is_zero():
            raise ZeroDivisionError
        return RF(self.d, self.n, reduce=False)

    def __truediv__(self, o):
        if not isinstance(o, RF):
            if o == 0:
                raise ZeroDivisionError
            r = RF.__new__(RF)
            r.n, r.d = self.n / _const(o)[0], self.d
            return r
        return self * o.inv()

    def __rtruediv__(self, o):
        return RF._c(o) * self.inv()

    def __eq__(self, o):
        o = RF._c(o)
        return self.n == o.n and self.d == o.d

    def __ne__(self, o):
        return not self.__eq__(o)

    def __hash__(self):
        return hash((str(self.n), str(self.d)))

    def is_zero(self):
        return self.n.is_zero()

    def __call__(self, x):
        """value at a rational x (Fraction)"""
        xf = flint.fmpq(x.numerator, x.denominator) if isinstance(x, Fr) else flint.fmpq(x)
        v = self.n(xf) / self.d(xf)
        return Fr(int(v.p), int(v.q))

    def __repr__(self):
        return f'({self.n})/({self.d})'


Q = RF(P([0, 1]), P([1]), reduce=False)   # the variable q


def sturm_count(p, lo, hi=None):
    """number of distinct real roots of the polynomial p (fmpq_poly) in the open interval (lo, hi) (hi = None: +inf),
    by Sturm's theorem, exact.  lo is rational and must not be a root."""
    if p.degree() <= 0:
        return 0
    g = p.gcd(p.derivative())
    if g.degree() > 0:
        p = p // g          # squarefree part: same distinct roots
    seq = [p, p.derivative()]
    while seq[-1].degree() > 0:
        r = -(seq[-2] % seq[-1])
        if r.is_zero():
            break
        seq.append(r)

    def changes_at(x):
        vals = [s(x) for s in seq]
        vals = [v for v in vals if v != 0]
        return sum(1 for a, b in zip(vals, vals[1:]) if (a > 0) != (b > 0))

    def changes_inf():
        vals = [s.leading_coefficient() for s in seq if not s.is_zero()]
        return sum(1 for a, b in zip(vals, vals[1:]) if (a > 0) != (b > 0))

    xl = flint.fmpq(lo.numerator, lo.denominator)
    assert p(xl) != 0
    if hi is None:
        return changes_at(xl) - changes_inf()
    xh = flint.fmpq(hi.numerator, hi.denominator)
    assert p(xh) != 0
    return changes_at(xl) - changes_inf() if False else changes_at(xl) - changes_at(xh)


def poly_positive_on_open_halfline(p):
    """(True, how) if the polynomial p (fmpq_poly) is > 0 for every q > 0; (False, why) otherwise.  Exact: either all
    coefficients >= 0 with some > 0 (then p(q) > 0 for q > 0), or no root in (0, inf) by Sturm and p(1) > 0."""
    cs = [c for c in p.coeffs()]
    if p.is_zero():
        return False, 'zero'
    if all(c >= 0 for c in cs):
        return True, 'nonnegative coefficients'
    # roots at q = 0 do not matter on (0, inf); strip the power of q
    k = 0
    while cs[k] == 0:
        k += 1
    pp = P(cs[k:])
    if pp(flint.fmpq(1)) <= 0:
        return False, f'value at q = 1 is {pp(flint.fmpq(1))}'
    # Sturm on (0, inf): 0 is not a root of pp
    n = sturm_count(pp, Fr(0), None)
    if n == 0:
        return True, 'Sturm: no root in (0, inf), positive at q = 1'
    return False, f'Sturm: {n} root(s) in (0, inf)'


def rf_positive(f):
    """f = n/d (RF): True if f > 0 for every q > 0 (both n and d positive there, or both negative), with a reason."""
    okn, hn = poly_positive_on_open_halfline(f.n)
    okd, hd = poly_positive_on_open_halfline(f.d)
    if okn and okd:
        return True, f'num ({f.n.degree()}, {hn}); den ({f.d.degree()}, {hd})'
    okn2, hn2 = poly_positive_on_open_halfline(-f.n)
    okd2, hd2 = poly_positive_on_open_halfline(-f.d)
    if okn2 and okd2:
        return True, f'num and den both negative on (0, inf)'
    return False, f'num: {hn}; den: {hd}'
