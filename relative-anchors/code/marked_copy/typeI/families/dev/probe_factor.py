"""Development probe (no proof): does the slice polynomial D_c det of a D2 family factor?  (A factorisation would
point to a further exact decomposition of the first part.)"""
import sys, pickle, time
import flint
fn = sys.argv[1]
d = pickle.load(open(fn, 'rb'))
ctx = flint.fmpq_mpoly_ctx.get(('nu', 's'), 'lex')
P = ctx.from_dict({e: flint.fmpq(p, q) for e, (p, q) in d.items()})
t0 = time.time()
c, fl = P.factor()
print(fn, ':', len(d), 'terms, degrees', P.degrees(), '; factors (terms, degrees, multiplicity):',
      [(len(f), f.degrees(), m) for f, m in fl], f'[{time.time()-t0:.1f}s]')
