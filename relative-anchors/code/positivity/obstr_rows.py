"""Where do the obstructed pairs of Lemma 2 sit?  For every block, the ordered pairs (a, b) with an exact certificate
(analyze.obstruction, directions {-3..3}^r) are split into pairs inside one row group (same row string) and pairs
between two row groups.  Usage: python obstr_rows.py k"""
import sys, time
from psdtools import load_blocks
from analyze import obstruction

k = int(sys.argv[1])
t0 = time.time()
tot_in = tot_x = 0
for blk in load_blocks(k):
    ob = obstruction(blk['Mg'])
    if not ob or not ob[0]:
        continue
    rows = [st[0] for st in blk['states']]
    same = sum(1 for (a, b, l) in ob[0] if rows[a] == rows[b])
    cross = len(ob[0]) - same
    tot_in += same
    tot_x += cross
    print(f'{blk["name"]} ({blk["kind"]}): {len(ob[0])} obstructed ordered pairs, {same} inside one row group, '
          f'{cross} between row groups  [{time.time() - t0:.1f}s]', flush=True)
print(f'RESULT k = {k}: inside row groups {tot_in}, between row groups {tot_x}')
print('PASS')
