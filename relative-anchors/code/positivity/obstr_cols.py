"""Like obstr_rows.py, with a three-way split of the obstructed ordered pairs (a, b): same row string, same column
string, neither.  Usage: python obstr_cols.py k"""
import sys, time
from psdtools import load_blocks
from analyze import obstruction

k = int(sys.argv[1])
t0 = time.time()
tot = [0, 0, 0]
for blk in load_blocks(k):
    ob = obstruction(blk['Mg'])
    if not ob or not ob[0]:
        continue
    st = blk['states']
    cnt = [0, 0, 0]
    for (a, b, l) in ob[0]:
        if st[a][0] == st[b][0]:
            cnt[0] += 1
        elif st[a][1] == st[b][1]:
            cnt[1] += 1
        else:
            cnt[2] += 1
    tot = [x + y for x, y in zip(tot, cnt)]
    print(f'{blk["name"]} ({blk["kind"]}): same row {cnt[0]}, same column {cnt[1]}, neither {cnt[2]}  '
          f'[{time.time() - t0:.1f}s]', flush=True)
print(f'RESULT k = {k}: same row {tot[0]}, same column {tot[1]}, neither {tot[2]}')
print('PASS')
