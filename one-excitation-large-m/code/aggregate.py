"""Aggregate per-job JSON results into tables (max over seeds)."""
import json, glob, collections, sys, os

rows = collections.defaultdict(list)
for fn in sorted(glob.glob(os.path.join('results', '*.json'))):
    base = os.path.basename(fn)
    if base.endswith('_s1.json') or (base.startswith('s45') and base.endswith('_s2.json')):
        continue          # skip the initial timing tests
    d = json.load(open(fn))
    if 'lemma' in d:
        key = ('lemma' if base.startswith('lemma_') else 'conc', d['lemma'], d['D'])
        rows[key].append(d)
    else:
        rows[('s45', d['kind'], d['param'])].append(d)
out = []
for key in sorted(rows):
    L = rows[key]
    best = max(L, key=lambda d: d['sup_found'])
    if key[0] == 's45':
        lim = best.get('claimed_bound', best.get('needed'))
        inits = sorted(set(d['init'] for d in L))
        out.append(f"{key[1]:8s} param={key[2]}  sup_found={best['sup_found']:.6f}  limit={lim}  runs={len(L)} "
                   f"(restarts {sum(d['restarts'] for d in L)} x {L[0]['steps']} steps) inits={inits} best_init={best['init']} "
                   f"Y_svals={['%.3g' % x for x in best['Y_svals'][:2]]}")
    else:
        ok = all(d['holds'] for d in L)
        out.append(f"{key[0]:5s} {key[1]:9s} D={key[2]}  sup_found={best['sup_found']:.6g}  bound={best['bound']:.6g}  holds={ok}  "
                   f"restarts={sum(d['restarts'] for d in L)} x steps={L[0]['steps']}")
print('\n'.join(out))
if len(sys.argv) > 1:
    open(sys.argv[1], 'w').write('\n'.join(out) + '\n')
