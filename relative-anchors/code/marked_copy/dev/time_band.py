import sys, time
sys.path.insert(0, '.')
import prove_fam
for fam in [(2, None, 4), (None, 2, 4), (5, None, 3), (None, 5, 3), (8, None, 0), (None, 8, 1)]:
    t0 = time.time()
    res, dt, info = prove_fam.prove(*fam)
    print(fam, res, f'{dt:.1f}s', info.get('maxterms'), flush=True)
