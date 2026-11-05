"""Every principal minor of every row group (states of a gauge-fixed block with one fixed row string: the vacuum, a
letter 1 at r, or a letter 2 at t), exact, with the method of allminors.py (all subsets).  Proposition R predicts that
these minors are subtraction-free for every k for the vacuum and the letter-1 rows.
Usage: python rowgroup_minors.py k"""
import sys, time
from psdtools import load_blocks
from allminors import run_block

k = int(sys.argv[1])
t0 = time.time()
ok = True
tot = dict(n=0, pos=0, mixed=0, bad=0)
for blk in load_blocks(k):
    groups = {}
    for i, st in enumerate(blk['states']):
        groups.setdefault(st[0], []).append(i)
    for a, G in sorted(groups.items()):
        sub = dict(blk)
        sub['Mg'] = [[blk['Mg'][p][q] for q in G] for p in G]
        sub['labels'] = [blk['labels'][p] for p in G]
        sub['states'] = [blk['states'][p] for p in G]
        kind = 'vacuum' if all(c == 0 for c in a) else ('letter 1' if 1 in a else 'letter 2')
        sub['name'] = blk['name'] + ' row ' + ''.join(map(str, a)) + f' ({kind})'
        st, fp = run_block(sub, None, quiet=True)
        print(f'{sub["name"]}: {len(G)} states, {st["n"]} minors, positive cores {st["pos"]}, mixed {st["mixed"]}, '
              f'point-check failures {st["badcheck"]}', flush=True)
        tot['n'] += st['n']
        tot['pos'] += st['pos']
        tot['mixed'] += st['mixed']
        tot['bad'] += st['badcheck']
        ok &= st['badcheck'] == 0 and fp
print(f'TOTAL: {tot["n"]} principal minors of row groups, positive cores {tot["pos"]}, mixed {tot["mixed"]}, '
      f'point-check failures {tot["bad"]}  [{time.time() - t0:.1f}s]')
print('PASS' if ok else 'FAIL')
