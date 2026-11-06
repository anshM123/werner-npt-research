"""Symbolic derivation (SymPy) of the reduced finite system of the relaxed form Psi (NOTES.md, Section 2):

    Psi = sum_h  v_h^T D(s_h) v_h  +  y^T C y,      v_h = (ah_h, ch_h, zh_h),   y = (rho, F_1, ..., F_J),

in the scaled coordinates a_h = w_h ah_h, c_h = w_h ch_h, zeta_h = w_t w_h zh_h (rho unscaled).  Every F_j is a global
linear functional F_j = sum_h phi_j(s_h) x_h with x in {ah, ch, zh} and phi_j a fixed basis function of s (a monomial
s^k or 1/(s + alpha)^k, from partial fractions in s over Q(sigma, S)); so the columns of U are the vectors
(phi_j(s_h))_h placed in the block of x.  The entries of D(s) are rational in (s, sigma, S) and polynomial in IG =
1/(2 + g), IGT = 1/(2 + gt) and moment symbols; the entries of C are polynomials in IG, IGT and the moment symbols
M_k = sum_h phi_k(s_h), with coefficients rational in (sigma, S).

Generic index: the symbol s is the squared weight of the index h being summed over.  sum1(E) returns
sum_h E(s_h, ah_h, ch_h, zh_h), classified into a local quadratic part (dict) and a global part (SymPy expression in
rho, the F_j and the constants)."""
import sympy as sp

s, sig, S = sp.symbols('s sigma S', positive=True)
ah, ch, zh = sp.symbols('ah ch zh')
rho = sp.Symbol('rho')
IG, IGT = sp.symbols('IG IGT')
USYM = sp.Symbol('u')
LOC = (ah, ch, zh)

d = s + S / 2 + sig            # d_h
d0 = s + (sig + S) / 2         # d0_h
dt0 = (3 * sig + S) / 2        # d0_t
K = sig + S / 2
kap = (2 * sig + S) / (sig + S)
dL = s + 3 * S / 2 + 3 * sig   # Lambda_h = dL/(2 d)
DOM = sp.QQ.frac_field(s, sig, S)


class Ctx:
    def __init__(self):
        self.F = {}        # (x, key) -> Symbol
        self.Finfo = {}    # Symbol -> (x, key)
        self.M = {}        # key -> Symbol
        self.Minfo = {}
        self.cc = {}

    def gvars(self):
        return [rho] + list(self.Finfo)

    def consts(self):
        return [IG, IGT, USYM] + list(self.Minfo)

    def canon(self, phi):
        """phi rational in s (parameters sigma, S) -> list of (coefficient in Q(sigma,S), basis key)."""
        phi = sp.cancel(phi)
        if phi in self.cc:
            return self.cc[phi]
        if phi == 0:
            return []
        pf = sp.apart(phi, s)
        out = {}
        for t in sp.Add.make_args(pf):
            c, b = t.as_independent(s, as_Add=False)
            num, den = sp.fraction(sp.together(b))
            lcd = sp.Poly(den, s).LC()
            lcn = sp.Poly(num, s).LC()
            key = sp.cancel(b * lcd / lcn)
            out[key] = out.get(key, 0) + c * lcn / lcd
        res = [(sp.cancel(c), k) for k, c in out.items() if sp.cancel(c) != 0]
        self.cc[phi] = res
        return res

    def fsym(self, x, key):
        if (x, key) not in self.F:
            sym = sp.Symbol(f'F_{x}_{len(self.F)}')
            self.F[(x, key)] = sym
            self.Finfo[sym] = (x, key)
        return self.F[(x, key)]

    def msym(self, key):
        if key == 1:
            return USYM
        if key == s:
            return S
        if key not in self.M:
            sym = sp.Symbol(f'M_{len(self.M)}')
            self.M[key] = sym
            self.Minfo[sym] = key
        return self.M[key]


class Qty:
    """a summed scalar: glob (SymPy expression, polynomial in global variables and constants) + loc (dict
    (x1, x2) -> coefficient in s; meaning sum_h coef(s_h) x1_h x2_h)."""

    def __init__(self, glob=0, loc=None):
        self.glob = sp.sympify(glob)
        self.loc = dict(loc or {})

    def __add__(self, o):
        o = o if isinstance(o, Qty) else Qty(o)
        loc = dict(self.loc)
        for k, v in o.loc.items():
            loc[k] = loc.get(k, 0) + v
        return Qty(self.glob + o.glob, loc)

    __radd__ = __add__

    def scale(self, c):
        return Qty(self.glob * c, {k: v * c for k, v in self.loc.items()})

    def __mul__(self, o):
        o = o if isinstance(o, Qty) else Qty(o)
        if self.loc and o.loc:
            raise ValueError('product of two local quadratic parts')
        if self.loc:
            return Qty(self.glob * o.glob, {k: v * o.glob for k, v in self.loc.items()})
        if o.loc:
            return Qty(self.glob * o.glob, {k: v * self.glob for k, v in o.loc.items()})
        return Qty(self.glob * o.glob)

    __rmul__ = __mul__

    def __neg__(self):
        return self.scale(-1)

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Qty) else Qty(o)))


def sum1(ctx, E):
    """sum over the generic index of E (polynomial of degree <= 2 in the local and global variables)."""
    E = sp.sympify(E)
    if E == 0:
        return Qty(0)
    gv = ctx.gvars()
    cs = ctx.consts()
    gens = list(LOC) + gv + cs
    P = sp.Poly(E, *gens, domain=DOM)
    nl = len(LOC)
    ng = len(gv)
    glob = 0
    loc = {}
    for mon, cf in P.terms():
        cfe = P.domain.to_sympy(cf)
        el = mon[:nl]
        eg = mon[nl:nl + ng]
        ec = mon[nl + ng:]
        cmon = sp.Mul(*[c ** e for c, e in zip(cs, ec)])
        gmon = sp.Mul(*[g ** e for g, e in zip(gv, eg)])
        dl = sum(el)
        if dl == 2:
            if sum(eg) != 0:
                raise ValueError('degree > 2')
            idx = [i for i in range(nl) for _ in range(el[i])]
            key = (LOC[idx[0]], LOC[idx[1]])
            loc[key] = loc.get(key, 0) + cfe * cmon
        elif dl == 1:
            x = LOC[el.index(1)]
            for c, b in ctx.canon(cfe):
                glob += c * cmon * gmon * ctx.fsym(str(x), b)
        else:
            for c, b in ctx.canon(cfe):
                glob += c * cmon * gmon * ctx.msym(b)
    return Qty(glob, loc)


class Mat:
    """scaled U x U matrix: diag(dg) + sum_i outer(x_i, y_i); entries are SymPy expressions at the generic index."""

    def __init__(self, dg=0, outers=None):
        self.dg = sp.sympify(dg)
        self.outers = list(outers or [])

    def __add__(self, o):
        return Mat(self.dg + o.dg, self.outers + o.outers)


def center(ctx, Mm):
    """(I - 1 p^T) M (I - p 1^T), p = s/S."""
    out = Mat(Mm.dg, [])
    D = Mm.dg
    if D != 0:
        mD = sum1(ctx, s * s * D / S ** 2)
        out.outers += [(-D * s / S, sp.Integer(1)), (sp.Integer(1), -D * s / S), (mD.glob, sp.Integer(1))]
        assert not mD.loc
    for x, y in Mm.outers:
        px = sum1(ctx, s * x / S)
        py = sum1(ctx, s * y / S)
        assert not px.loc and not py.loc
        out.outers.append((x - px.glob, y - py.glob))
    return out


def build(verbose=False):
    ctx = Ctx()
    # global functionals used by f, Acal, the coupling
    A1 = sum1(ctx, s / d * ah).glob          # sum rho_h ah_h
    C1 = sum1(ctx, s / d * ch).glob          # sum rho_h ch_h
    P = sum1(ctx, s * ch).glob               # p = w.c
    Z1 = sum1(ctx, s * zh).glob              # sum s_h z_h
    Z0r = sum1(ctx, s / d0 * zh).glob        # sum rho0_h z_h
    sum1(ctx, s * ah)                         # the constraint functional sum s_h ah_h (a orthogonal to w)
    # ---- sig f(a, c)
    fl = sum1(ctx, (s - sig * s / d) * (ah ** 2 + ch ** 2) + 4 * IG * s * (s / d) * ah * ch)
    fg = sig * IG * (A1 ** 2 + C1 ** 2) - sig / (2 * K ** 2) * P ** 2 + 2 * IG * P * A1
    Psi = (fl + Qty(fg)).scale(sig)
    # ---- Acal
    m = kap * rho + sig * Z1 / (sig + S)
    Psi = Psi + sum1(ctx, sig * s * zh ** 2) + Qty(2 * rho ** 2 - m ** 2 / 2)
    # ---- coupling <Hfull, Omega> = sum_{r,h} s_r s_h Hhat_rh (ah_r + ch_h)
    Y = rho + Z1
    Wyt = (rho - sig * (rho / dt0 + Z0r) * IGT) / dt0
    Hf = Mat((rho - m) / s + sig * (Y / 2 + s * zh) / (s * d0),
             [(m / (2 * K), sp.Integer(1)),
              (sig * IGT * (Y / 2 + s * zh) / (dt0 * d0), sp.Integer(1)),
              (-sig * IGT * (Y / 2 + s * zh) / d0, 1 / d0),
              (-Wyt / 2, sp.Integer(1)),
              (sp.Integer(1), sig * (zh - (rho / dt0 + Z0r) * IGT) / (2 * d0))])
    cpl = sum1(ctx, s * s * Hf.dg * (ah + ch))
    for x, y in Hf.outers:
        cpl = cpl + sum1(ctx, s * x * ah) * sum1(ctx, s * y) + sum1(ctx, s * x) * sum1(ctx, s * y * ch)
    Psi = Psi + cpl
    # ---- relaxation: R = Pi (sig Gcal + Hhat) Pi, scaled; T = (12/5) sum s_r s_h R^2 + sum (2d/dL - 12/5) s^2 R_hh^2
    Mh = Mat(sig * (2 * (ah + ch) / d - P / (K * d) - Y / ((sig + S) * d0) + zh / d0),
             [(1 / d, -2 * sig * IG * s * ah / d),
              (-2 * sig * IG * s * ch / d, 1 / d),
              (-sig * IG * P / d, 1 / d),
              (-sig * IGT * (Y / 2 + s * zh) / d0, 1 / d0)])
    R = center(ctx, Mh)
    Dc = R.dg
    T = sum1(ctx, s * s * Dc * Dc)
    for x, y in R.outers:
        T = T + sum1(ctx, 2 * s * s * Dc * x * y)
    n = len(R.outers)
    for i in range(n):
        xi, yi = R.outers[i]
        for j in range(n):
            xj, yj = R.outers[j]
            T = T + sum1(ctx, s * xi * xj) * sum1(ctx, s * yi * yj)
    T = T.scale(sp.Rational(12, 5))
    Rhh = Dc + sum((x * y for x, y in R.outers), sp.Integer(0))
    T = T + sum1(ctx, (2 * d / dL - sp.Rational(12, 5)) * s * s * Rhh * Rhh)
    Psi = Psi - T.scale(sp.Rational(1, 4))
    return ctx, Psi


def ints(pe):
    return [(tuple(m), int(sp.Rational(c).p), int(sp.Rational(c).q)) for m, c in pe.terms()]


def ratterms(expr, vars_):
    num, den = sp.fraction(sp.together(expr))
    tn = [(tuple(m), int(sp.Rational(c).p), int(sp.Rational(c).q)) for m, c in sp.Poly(num, *vars_).terms()]
    td = [(tuple(m), int(sp.Rational(c).p), int(sp.Rational(c).q)) for m, c in sp.Poly(den, *vars_).terms()]
    return tn, td


if __name__ == '__main__':
    import time, pickle
    t0 = time.time()
    ctx, Psi = build()
    print(f'built in {time.time() - t0:.1f} s; {len(ctx.Finfo)} functionals, {len(ctx.Minfo)} moment symbols', flush=True)
    for sym, (x, key) in ctx.Finfo.items():
        print(' ', sym, x, key)
    for sym, key in ctx.Minfo.items():
        print(' ', sym, key)
    gens = ctx.gvars()
    consts = ctx.consts()
    ng = len(gens)
    dom2 = sp.QQ.frac_field(sig, S)
    acc = {}
    terms = sp.Add.make_args(Psi.glob)
    for k, t in enumerate(terms):
        P = sp.Poly(t, *gens, *consts, domain=dom2)
        for mon, cf in P.as_dict(native=True).items():
            acc[mon] = acc[mon] + cf if mon in acc else cf
    C = []
    for mon, cf in acc.items():
        if cf == 0:
            continue
        eg = mon[:ng]
        assert sum(eg) == 2
        idx = [i for i in range(ng) for _ in range(eg[i])]
        C.append((idx[0], idx[1], tuple(mon[ng:]), ints(cf.numer), ints(cf.denom)))
    dom3 = sp.QQ.frac_field(s, sig, S)
    loc = []
    for (a, b), v in Psi.loc.items():
        Pl = sp.Poly(v, *consts, domain=dom3)
        for mon, cf in Pl.as_dict(native=True).items():
            loc.append((str(a), str(b), tuple(mon), ints(cf.numer), ints(cf.denom)))
    F = [(x,) + ratterms(key, (s, sig, S)) for sym, (x, key) in ctx.Finfo.items()]
    Mb = [ratterms(key, (s, sig, S)) for sym, key in ctx.Minfo.items()]
    j0 = [1 + j for j, (sym, (x, key)) in enumerate(ctx.Finfo.items()) if x == 'ah' and sp.simplify(key - s) == 0]
    assert len(j0) == 1
    pickle.dump(dict(gens=[str(g) for g in gens], consts=[str(c) for c in consts], C=C, loc=loc, F=F, M=Mb, j0=j0[0],
                     Fkeys=[str(key) for sym, (x, key) in ctx.Finfo.items()],
                     Mkeys=[str(key) for sym, key in ctx.Minfo.items()]),
                open('reduced_compact.pkl', 'wb'))
    print(f'saved reduced_compact.pkl: {len(C)} C-terms, {len(loc)} D-terms ({time.time() - t0:.1f} s)')
