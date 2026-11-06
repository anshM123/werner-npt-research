"""Exact evaluation of the reduced finite system (reduced.pkl, written by reduce_cas.py) at rational weights.

  Psi = blockdiag_h D(s_h) + U C U^T   in the coordinates (ah_1..ah_u, ch_1..ch_u, zh_1..zh_u, rho),
  a_h = w_h ah_h, c_h = w_h ch_h, zeta_h = w_t w_h zh_h.
Columns of U: column 0 is e_rho; column j >= 1 is the functional F_j = sum_h phi_j(s_h) x_h (x in {ah, ch, zh}),
i.e. the vector (phi_j(s_h))_h in the block of x.  Constants: IG = 1/(2 + g), g = sum_h s_h/d_h; IGT = 1/(2 + gt),
gt = sigma/d0_t + sum_h s_h/d0_h; u = |U|; M_k = sum_h psi_k(s_h).  The constraint (a orthogonal to w) is
sum_h s_h ah_h = 0, i.e. y_{j0} = 0 for the functional (ah, s)."""
import pickle
from fractions import Fraction as Fr
import sympy as sp
import flint

Q = flint.fmpq
s, sig, S = sp.symbols('s sigma S', positive=True)
IG, IGT, USYM = sp.symbols('IG IGT u')
rho = sp.Symbol('rho')


def eval_srepr(x):
    return sp.sympify(eval(x, {k: getattr(sp, k) for k in dir(sp) if not k.startswith('_')}))


def _poly_terms(pe):
    """PolyElement -> list of (exponent tuple, fmpq)."""
    return [(m, Q(int(sp.Rational(c).p), int(sp.Rational(c).q))) for m, c in pe.terms()]


def _ev(terms, vals):
    tot = Q(0)
    for m, c in terms:
        t = c
        for v, e in zip(vals, m):
            if e:
                t *= v ** e
        tot += t
    return tot


class Reduced:
    def __init__(self, fn='reduced.pkl'):
        D = pickle.load(open(fn, 'rb'))
        self.Finfo = [(sp.Symbol(k), v[0], eval_srepr(v[1])) for k, v in D['F'].items()]
        self.Minfo = [(sp.Symbol(k), eval_srepr(v)) for k, v in D['M'].items()]
        glob = eval_srepr(D['glob'])
        self.Fsyms = [f for f, _, _ in self.Finfo]
        self.gens = [rho] + self.Fsyms
        self.consts = [IG, IGT, USYM] + [m for m, _ in self.Minfo]
        dom2 = sp.QQ.frac_field(sig, S)
        P = sp.Poly(glob, *self.gens, *self.consts, domain=dom2)
        ng = len(self.gens)
        self.Cterms = []          # (i, j, const exponents, numer terms, denom terms)
        for mon, cf in P.as_dict(native=True).items():
            eg = mon[:ng]
            assert sum(eg) == 2, eg
            idx = [i for i in range(ng) for _ in range(eg[i])]
            self.Cterms.append((idx[0], idx[1], mon[ng:], _poly_terms(cf.numer), _poly_terms(cf.denom)))
        dom3 = sp.QQ.frac_field(s, sig, S)
        self.loc = []
        for (a, b), v in D['loc'].items():
            Pl = sp.Poly(eval_srepr(v), *self.consts, domain=dom3)
            for mon, cf in Pl.as_dict(native=True).items():
                self.loc.append((a, b, mon, _poly_terms(cf.numer), _poly_terms(cf.denom)))
        # basis functions as (numer terms, denom terms) in (s, sigma, S)
        def ratterms(expr):
            num, den = sp.fraction(sp.together(expr))
            tn = [(m, Q(int(sp.Rational(c).p), int(sp.Rational(c).q))) for m, c in sp.Poly(num, s, sig, S).terms()]
            td = [(m, Q(int(sp.Rational(c).p), int(sp.Rational(c).q))) for m, c in sp.Poly(den, s, sig, S).terms()]
            return (tn, td)
        self.Fbasis = [(x, ratterms(key)) for _, x, key in self.Finfo]
        self.Mbasis = [ratterms(key) for _, key in self.Minfo]
        self.j0 = None
        for j, (f, x, key) in enumerate(self.Finfo):
            if x == 'ah' and sp.simplify(key - s) == 0:
                self.j0 = 1 + j
        assert self.j0 is not None

    def numeric(self, wt, wU):
        wt = to_q(wt); wU = [to_q(x) for x in wU]
        u = len(wU)
        sv = [x * x for x in wU]
        sg = wt * wt
        SS = sum(sv, Q(0))
        dd = [x + SS / 2 + sg for x in sv]
        d0 = [x + (sg + SS) / 2 for x in sv]
        dt0 = (3 * sg + SS) / 2
        g = sum((sv[i] / dd[i] for i in range(u)), Q(0))
        gt = sg / dt0 + sum((sv[i] / d0[i] for i in range(u)), Q(0))

        def fun(terms, sh):
            return _ev(terms[0], (sh, sg, SS)) / _ev(terms[1], (sh, sg, SS))
        cv = [1 / (2 + g), 1 / (2 + gt), Q(u)] + [sum((fun(t, x) for x in sv), Q(0)) for t in self.Mbasis]

        def cmon(ex):
            t = Q(1)
            for v, e in zip(cv, ex):
                if e:
                    t *= v ** e
            return t
        J = len(self.gens)
        C = [[Q(0)] * J for _ in range(J)]
        for i, j, ex, nt, dt in self.Cterms:
            val = cmon(ex) * _ev(nt, (sg, SS)) / _ev(dt, (sg, SS))
            if i == j:
                C[i][i] += val
            else:
                C[i][j] += val / 2
                C[j][i] += val / 2
        names = ['ah', 'ch', 'zh']
        Db = []
        for h in range(u):
            B3 = [[Q(0)] * 3 for _ in range(3)]
            for a, b, ex, nt, dt in self.loc:
                val = cmon(ex) * _ev(nt, (sv[h], sg, SS)) / _ev(dt, (sv[h], sg, SS))
                ia, ib = names.index(a), names.index(b)
                if ia == ib:
                    B3[ia][ia] += val
                else:
                    B3[ia][ib] += val / 2
                    B3[ib][ia] += val / 2
            Db.append(B3)
        Ucols = [(names.index(x), [fun(t, sv[h]) for h in range(u)]) for x, t in self.Fbasis]
        return dict(D=Db, C=C, U=Ucols, u=u, sv=sv, sg=sg, w=wU, wt=wt, consts=cv)


def to_q(x):
    if isinstance(x, flint.fmpq):
        return x
    if isinstance(x, int):
        return Q(x)
    return Q(x.numerator, x.denominator)


def assemble(num):
    """Psi = blockdiag D + U C U^T as an fmpq_mat in the order (ah (u), ch (u), zh (u), rho)."""
    u = num['u']
    n = 3 * u + 1
    Mx = flint.fmpq_mat(n, n)
    for i in range(u):
        for a in range(3):
            for b in range(3):
                Mx[a * u + i, b * u + i] = to_q(num['D'][i][a][b])
    J = len(num['C'])
    Um = flint.fmpq_mat(n, J)
    Um[n - 1, 0] = Q(1)
    for j, (blk, vals) in enumerate(num['U']):
        for i in range(u):
            Um[blk * u + i, 1 + j] = to_q(vals[i])
    Cm = flint.fmpq_mat(J, J, [to_q(num['C'][i][j]) for i in range(J) for j in range(J)])
    return Mx + Um * Cm * Um.transpose(), Mx, Um, Cm
