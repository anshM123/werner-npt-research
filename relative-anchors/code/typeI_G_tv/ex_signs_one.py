"""Exploration (no claim rests on it beyond the structural remarks of NOTES.md, Sections 3 and 8): monomial content (powers of a = s_A, b = s_B) and coefficient signs of the u1 = 1 closed forms."""
import pickle, time, flint, sys
D = pickle.load(open('tl_inv_one.pkl', 'rb'))
names = D['names']; keep = D['keep']
ctx = flint.fmpq_mpoly_ctx.get(tuple(names), 'lex')
Q = flint.fmpq
def P(d): return ctx.from_dict({m: Q(p, q) for m, (p, q) in d.items()})
G = D['G']
for (i, j), (n, d) in sorted(G.items()):
    N = P(n); Dd = P(d)
    mins = [min(m[k] for m in N.monoms()) for k in range(len(names))]
    neg = sum(1 for c in N.coeffs() if c < 0); pos = len(N) - neg
    negd = sum(1 for c in Dd.coeffs() if c < 0)
    print('G[%s,%s]: monomial content %s, num coeffs +%d -%d; den negative coeffs %d' % (keep[i], keep[j], dict(zip(names, mins)), pos, neg, negd), flush=True)
print('exploration done')
