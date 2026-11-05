"""Exploration: factor the common denominator of the u = 3, m = 1 closed form."""
import pickle, flint, sys
u = int(sys.argv[1]) if len(sys.argv) > 1 else 3
D = pickle.load(open(f'closed_m1_u{u}.pkl', 'rb'))
names = D['names']
ctx = flint.fmpq_mpoly_ctx.get(tuple(names), 'lex')
G = D['G']
d = len(G)
dens = [ctx.from_dict(G[i][j][1]) for i in range(d) for j in range(d)]
print('all dens equal:', all(x == dens[0] for x in dens))
c, fac = dens[0].factor()
print('den const', c)
for f, e in fac:
    print('  factor deg', f.total_degree(), 'degrees', f.degrees(), 'terms', len(f), '^', e, ':', str(f)[:200])
