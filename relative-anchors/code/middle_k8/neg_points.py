"""Negative controls at explicit points, from the Lemma B builder (j1lib.py, bscale), independent of the certificates:
for every family and every configuration with two levels, a rational square q0 at which the gauge-fixed block with
the B-part multiplied by 1001/1000 is not positive definite (exact symmetric elimination), while the true gauge-fixed
block is positive definite there.
python neg_points.py FAMILIES [k=5] [reps=1] [first=N]
Prints PASS (every two-level configuration has such a point) or FAIL on the last line."""
import sys, time
from fractions import Fraction as Fr
import k5lib as K

t0 = time.time()
args = [a for a in sys.argv[1:] if '=' not in a]
opts = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
k = int(opts.get('k', 5))
fams = args[0].split(',')
ok = True
cands = [Fr(1, 10 ** 6), Fr(10 ** 4), Fr(10 ** 6), Fr(1, 10 ** 8), Fr(10 ** 8), Fr(1, 10 ** 10), Fr(10 ** 10)]
for f in fams:
    j, m = K.fam(f)
    T = tuple(range(m))
    confs = K.rep_configs(k, m) if opts.get('reps') == '1' else K.configs(k, m)
    confs = [c for c in confs if c not in ((0, 0), (m, k - m))]
    if 'first' in opts:   # only the first N two-level configurations of each family (k = 8: time)
        confs = confs[:int(opts['first'])]
    for (a, b) in confs:
        if (a, b) in ((0, 0), (m, k - m)):
            continue          # one level: no ratio
        found = None
        for q0 in cands:
            s = K.config_s(k, m, a, b, q0, Fr(1))
            w = [K.rsqrt(x) for x in s]
            st, Mb = K.build_j1(k, j, T, w, bscale=Fr(1001, 1000))
            vecs, dels = K.known_kernel(k, j, T, st, s, Fr(0), Fr(1))
            keep = [i for i in range(len(st)) if i not in dels]
            pdb, piv = K.ldl_pd(K.principal(Mb, keep))
            if not pdb:
                st, M1 = K.build_j1(k, j, T, w)
                pd1, _ = K.ldl_pd(K.principal(M1, keep))
                if pd1:
                    found = (q0, piv)
                    break
        ok &= found is not None
        if found:
            print(f'{f} {(a, b)}: at q = {found[0]} the scaled gauge-fixed block is not PD (pivot {found[1] + 1} <= 0), '
                  f'the true one is PD', flush=True)
        else:
            print(f'{f} {(a, b)}: NO failing point found among {[str(c) for c in cands]}', flush=True)
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
