"""Exact evaluation (python-flint fmpq) of the reduced finite system from reduced_compact.pkl (written by
reduce_cas.py; no SymPy needed).  Same conventions as reduced_eval.py:

  Psi = blockdiag_h D(s_h) + U C U^T  in the coordinates (ah_1..ah_u, ch_1..ch_u, zh_1..zh_u, rho),
  y = (rho, F_1, ..., F_30), F_j = sum_h phi_j(s_h) x_h;  constraint y_j0 = 0 (a orthogonal to w)."""
import pickle
import flint

Q = flint.fmpq


def to_q(x):
    if isinstance(x, flint.fmpq):
        return x
    if isinstance(x, int):
        return Q(x)
    return Q(x.numerator, x.denominator)


def _conv(terms):
    return [(m, Q(p, q)) for m, p, q in terms]


def _ev(terms, vals):
    tot = Q(0)
    for m, c in terms:
        t = c
        for v, e in zip(vals, m):
            if e:
                t *= v ** e
        tot += t
    return tot


class ReducedC:
    def __init__(self, fn='reduced_compact.pkl'):
        D = pickle.load(open(fn, 'rb'))
        self.gens = D['gens']
        self.consts = D['consts']
        self.j0 = D['j0']
        self.Cterms = [(i, j, ex, _conv(n), _conv(d)) for i, j, ex, n, d in D['C']]
        self.loc = [(a, b, ex, _conv(n), _conv(d)) for a, b, ex, n, d in D['loc']]
        self.Fbasis = [(x, (_conv(n), _conv(d))) for x, n, d in D['F']]
        self.Mbasis = [(_conv(n), _conv(d)) for n, d in D['M']]
        self.Fkeys = D['Fkeys']
        self.Mkeys = D['Mkeys']

    def numeric(self, wt, wU):
        wt = to_q(wt)
        wU = [to_q(x) for x in wU]
        u = len(wU)
        sv = [x * x for x in wU]
        sg = wt * wt
        SS = sum(sv, Q(0))
        dd = [x + SS / 2 + sg for x in sv]
        d0 = [x + (sg + SS) / 2 for x in sv]
        dt0 = (3 * sg + SS) / 2
        g = sum((sv[i] / dd[i] for i in range(u)), Q(0))
        gt = sg / dt0 + sum((sv[i] / d0[i] for i in range(u)), Q(0))

        def fun(t, sh):
            return _ev(t[0], (sh, sg, SS)) / _ev(t[1], (sh, sg, SS))
        cv = [1 / (2 + g), 1 / (2 + gt), Q(u)] + [sum((fun(t, x) for x in sv), Q(0)) for t in self.Mbasis]
        cache = {}

        def cmon(ex):
            if ex not in cache:
                t = Q(1)
                for v, e in zip(cv, ex):
                    if e:
                        t *= v ** e
                cache[ex] = t
            return cache[ex]
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


def assemble(num):
    u = num['u']
    n = 3 * u + 1
    Mx = flint.fmpq_mat(n, n)
    for i in range(u):
        for a in range(3):
            for b in range(3):
                Mx[a * u + i, b * u + i] = num['D'][i][a][b]
    J = len(num['C'])
    Um = flint.fmpq_mat(n, J)
    Um[n - 1, 0] = Q(1)
    for j, (blk, vals) in enumerate(num['U']):
        for i in range(u):
            Um[blk * u + i, 1 + j] = vals[i]
    Cm = flint.fmpq_mat(J, J, [num['C'][i][j] for i in range(J) for j in range(J)])
    return Mx + Um * Cm * Um.transpose(), Mx, Um, Cm
