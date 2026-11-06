"""Minimal exact rational-function field QQ(x_1..x_r) on top of python-flint fmpq_mpoly (fast multivariate gcd).
Elements are reduced fractions num/den with den monic-positive (leading coefficient > 0).  Used by exact_group_flint.py
as a drop-in replacement for sympy's field (same interface subset: + - * /, == 0, .numer, .denom)."""
import flint


class Poly:
    """thin wrapper so that .coeffs() and .terms() behave like sympy PolyElement for the certificate code."""
    __slots__ = ('p',)

    def __init__(self, p):
        self.p = p

    def coeffs(self):
        return list(self.p.coeffs())

    def terms(self):
        return list(range(len(self.p)))


class Field:
    def __init__(self, names):
        self.ctx = flint.fmpq_mpoly_ctx.get(tuple(names), 'lex')
        self.gens = [F(self, g, self.ctx.from_dict({(0,) * len(names): 1})) for g in self.ctx.gens()]

    def __call__(self, c):
        if isinstance(c, F):
            return c
        return F(self, self.ctx.from_dict({(0,) * self.ctx.nvars(): flint.fmpq(c)}) if c != 0 else self.ctx.from_dict({}),
                 self.ctx.from_dict({(0,) * self.ctx.nvars(): 1}), reduce=False)


class F:
    __slots__ = ('K', 'n', 'd')

    def __init__(self, K, n, d, reduce=True):
        self.K = K
        if reduce and not n.is_zero():
            g = n.gcd(d)
            if not g.is_one():
                n = n / g
                d = d / g
        if n.is_zero():
            d = K.ctx.from_dict({(0,) * K.ctx.nvars(): 1})
        lc = d.leading_coefficient()
        if lc != 1:
            n = n / lc
            d = d / lc
        self.n, self.d = n, d

    def _c(self, o):
        return o if isinstance(o, F) else self.K(o)

    def __add__(self, o):
        o = self._c(o)
        if self.n.is_zero():
            return o
        if o.n.is_zero():
            return self
        if self.d == o.d:
            return F(self.K, self.n + o.n, self.d)
        g = self.d.gcd(o.d)
        d1, d2 = self.d / g, o.d / g
        return F(self.K, self.n * d2 + o.n * d1, self.d * d2)

    __radd__ = __add__

    def __neg__(self):
        return F(self.K, -self.n, self.d, reduce=False)

    def __sub__(self, o):
        return self + (-self._c(o))

    def __rsub__(self, o):
        return self._c(o) - self

    def __mul__(self, o):
        o = self._c(o)
        if self.n.is_zero() or o.n.is_zero():
            return self.K(0)
        g1 = self.n.gcd(o.d)
        g2 = o.n.gcd(self.d)
        n = (self.n / g1) * (o.n / g2)
        d = (self.d / g2) * (o.d / g1)
        return F(self.K, n, d, reduce=False)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._c(o)
        if o.n.is_zero():
            raise ZeroDivisionError
        return self * F(self.K, o.d, o.n, reduce=False)

    def __eq__(self, o):
        if isinstance(o, (int,)) and o == 0:
            return self.n.is_zero()
        o = self._c(o)
        return (self - o).n.is_zero()

    def __ne__(self, o):
        return not self.__eq__(o)

    def __hash__(self):
        return hash((str(self.n), str(self.d)))

    @property
    def numer(self):
        return Poly(self.n)

    @property
    def denom(self):
        return Poly(self.d)
