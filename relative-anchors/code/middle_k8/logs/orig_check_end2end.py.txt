"""End-to-end exact checks at two-value rational points, from the Lemma B builder (j1lib.py), independent of the
certificates: for every family and configuration and several rational squares q0 (including near-degenerate ones),
 - every known kernel vector is annihilated,
 - the nullity of the whole block equals the number of known kernel vectors (exact rank over Q),
 - the gauge-fixed block is positive definite (exact symmetric elimination), in the gauge of the certificates and in
   a second gauge (the states (1@0, y) deleted instead of (1@k-1, y)).
python check_end2end.py FAMILIES [k=5]
Prints PASS or FAIL on the last line."""
import sys, time
from fractions import Fraction as Fr
import k5lib as K

t0 = time.time()
args = [a for a in sys.argv[1:] if '=' not in a]
opts = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
k = int(opts.get('k', 5))
fams = args[0].split(',')
ok = True
qs = [Fr(1, 4), Fr(9, 4), Fr(1, 10 ** 6), Fr(10 ** 6), Fr(121, 100)]
for f in fams:
    j, m = K.fam(f)
    T = tuple(range(m))
    cnt = 0
    for (a, b) in (K.rep_configs(k, m) if opts.get('reps') == '1' else K.configs(k, m)):
        for q0 in qs:
            s = K.config_s(k, m, a, b, q0, Fr(1))
            w = [K.rsqrt(x) for x in s]
            st, M = K.build_j1(k, j, T, w)
            n = len(st)
            vecs, dels = K.known_kernel(k, j, T, st, s, Fr(0), Fr(1))
            ann = all(all(x == 0 for x in K.matvec(M, v)) for v in vecs)
            nul = n - K.rank_q(M)
            pd1, _ = K.ldl_pd(K.principal(M, [i for i in range(n) if i not in dels]))
            vecs2, dels2 = K.known_kernel(k, j, T, st, s, Fr(0), Fr(1), r0=0)
            pd2, _ = K.ldl_pd(K.principal(M, [i for i in range(n) if i not in dels2]))
            good = ann and nul == len(vecs) and pd1 and pd2
            ok &= good
            cnt += 1
            if not good:
                print(f'  {f} {(a, b)} q = {q0}: annihilated {ann}, nullity {nul} (known {len(vecs)}), PD {pd1} {pd2}')
    print(f'{f}: {cnt} (configuration, q) points, known kernel exactly and gauge-fixed PD in two gauges: {ok} '
          f'[{time.time() - t0:.1f}s]', flush=True)
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
