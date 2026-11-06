"""The relaxed form Psi on the perpendicular subspace p = {sum s_h ah_h = sum s_h ch_h = sum s_h zh_h = 0, rho = 0},
with natural functionals (no partial fractions): same symbolic program as reduce_cas.py, but every functional keeps its
natural weight function phi(s) (normalised to leading coefficient 1 in numerator and denominator), so the global
part on p is a quadratic form in sums such as sum_h s_h^2 ah_h/d_h^2 or sum_h s_h ch_h/d0_h.
Output: perp_compact.pkl (same layout as reduced_compact.pkl, restricted to p)."""
import time, pickle
import sympy as sp
import reduce_cas as rc
from reduce_cas import s, sig, S, ah, ch, zh, rho, IG, IGT, USYM, ints, ratterms


class NatCtx(rc.Ctx):
    def canon(self, phi):
        phi = sp.cancel(phi)
        if phi in self.cc:
            return self.cc[phi]
        if phi == 0:
            return []
        num, den = sp.fraction(phi)
        lcn = sp.Poly(num, s).LC()
        lcd = sp.Poly(den, s).LC()
        c = sp.cancel(lcn / lcd)
        key = sp.cancel(phi / c)
        res = [(c, key)]
        self.cc[phi] = res
        return res


if __name__ == '__main__':
    t0 = time.time()
    rc.Ctx = NatCtx
    ctx, Psi = rc.build()
    print(f'built in {time.time() - t0:.1f} s; {len(ctx.Finfo)} natural functionals, {len(ctx.Minfo)} moments',
          flush=True)
    zero = {rho: 0}
    for sym, (x, key) in ctx.Finfo.items():
        if sp.simplify(key - s) == 0:
            zero[sym] = 0
            print('  set to zero on p:', sym, x, key)
    glob = Psi.glob.subs(zero)
    gens = [g for g in ctx.gvars() if g not in zero]
    consts = ctx.consts()
    ng = len(gens)
    dom2 = sp.QQ.frac_field(sig, S)
    acc = {}
    for t in sp.Add.make_args(sp.expand(glob) if False else glob):
        P = sp.Poly(t, *gens, *consts, domain=dom2)
        for mon, cf in P.as_dict(native=True).items():
            acc[mon] = acc[mon] + cf if mon in acc else cf
    C = []
    for mon, cf in acc.items():
        if cf == 0:
            continue
        eg = mon[:ng]
        if sum(eg) == 0:
            continue
        assert sum(eg) == 2, eg
        idx = [i for i in range(ng) for _ in range(eg[i])]
        C.append((idx[0], idx[1], tuple(mon[ng:]), ints(cf.numer), ints(cf.denom)))
    dom3 = sp.QQ.frac_field(s, sig, S)
    loc = []
    for (a, b), v in Psi.loc.items():
        Pl = sp.Poly(v, *consts, domain=dom3)
        for mon, cf in Pl.as_dict(native=True).items():
            loc.append((str(a), str(b), tuple(mon), ints(cf.numer), ints(cf.denom)))
    Fl = [(ctx.Finfo[g][0],) + ratterms(ctx.Finfo[g][1], (s, sig, S)) for g in gens]
    Mb = [ratterms(key, (s, sig, S)) for sym, key in ctx.Minfo.items()]
    pickle.dump(dict(gens=[str(g) for g in gens], consts=[str(c) for c in consts], C=C, loc=loc, F=Fl, M=Mb,
                     Fkeys=[str(ctx.Finfo[g][1]) for g in gens], Mkeys=[str(k) for k in ctx.Minfo.values()]),
                open('perp_compact.pkl', 'wb'))
    print(f'saved perp_compact.pkl: {len(gens)} functionals on p, {len(C)} C-terms, {len(loc)} D-terms '
          f'({time.time() - t0:.1f} s)')
    for g in gens:
        print('  ', g, ctx.Finfo[g][0], ctx.Finfo[g][1])
