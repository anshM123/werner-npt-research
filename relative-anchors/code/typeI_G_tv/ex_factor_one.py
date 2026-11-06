"""Exploration (no claim rests on it beyond the structural remarks of NOTES.md, Sections 3 and 8): factor structure of the denominators and numerators of the u1 = 1 closed forms with arbitrary T-weights."""
import pickle, time, flint, sys
D = pickle.load(open('tl_inv_one.pkl', 'rb'))
names = D['names']
ctx = flint.fmpq_mpoly_ctx.get(tuple(names), 'lex')
Q = flint.fmpq
def P(d): return ctx.from_dict({m: Q(p, q) for m, (p, q) in d.items()})
G = D['G']
for key in [(0, 0), (0, 3)]:
    n, d = G[key]
    t0 = time.time()
    c, fl = P(d).factor()
    print(key, 'den factors:', [(len(f), e, f.degrees()) for f, e in fl], 'const', c, '%.1fs' % (time.time() - t0), flush=True)
    for f, e in fl:
        if len(f) < 40:
            print('   ', f, '^', e)
    t0 = time.time()
    c, fl = P(n).factor()
    print(key, 'num factors:', [(len(f), e, f.degrees()) for f, e in fl], '%.1fs' % (time.time() - t0), flush=True)
    for f, e in fl:
        if len(f) < 40:
            print('   ', f, '^', e)
print(names)
print('exploration done')
