"""Exploration: factor the common denominators of G00, G01 (u = 2) and check the dependence claims of (S4)."""
import pickle, flint, time
ctx = flint.fmpq_mpoly_ctx.get(('s0', 's1', 's2', 's3', 's4'), 'lex')
s0, s1, s2, s3, s4 = ctx.gens()
D = pickle.load(open('closed_u2.pkl', 'rb'))
P = lambda d: ctx.from_dict(d)
for name in ('G00', 'G01'):
    dens = [P(D[name][i][j][1]) for i in range(3) for j in range(3)]
    print(name, 'all dens equal:', all(d == dens[0] for d in dens))
    c, fac = dens[0].factor()
    for f_, e_ in fac:
        print('   ', e_, f_.total_degree(), str(f_)[:150])
# S4 check for G00: substitute s1 -> sig - s0 - s2 with a new variable and see s2 disappears
ctx2 = flint.fmpq_mpoly_ctx.get(('s0', 'sg', 's2', 's3', 's4'), 'lex')
a0, sg, a2, a3, a4 = ctx2.gens()
for i in range(3):
    for j in range(3):
        n_, d_ = D['G00'][i][j]
        nn = P(n_).compose(a0, sg - a0 - a2, a2, a3, a4, ctx=ctx2)
        dd = P(d_).compose(a0, sg - a0 - a2, a2, a3, a4, ctx=ctx2)
        g = nn.gcd(dd)
        nn, dd = nn / g, dd / g
        print('G00', i, j, 'degree in s2 after s1 -> sg - s0 - s2:', nn.degrees()[2], dd.degrees()[2])
