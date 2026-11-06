"""Exact kernel check at two-value rational points (exploration): for every family and a few configurations, build
the block from Lemma B (j1lib.py) at rational weights, and report the nullity, whether the known kernel vectors are
annihilated, and whether the gauge-fixed block is positive definite (exact LDL).
python explore_kernel.py [families]"""
import sys, time
from fractions import Fraction as Fr
import k5lib as K

t0 = time.time()
k = int(sys.argv[2]) if len(sys.argv) > 2 else 5
fams = sys.argv[1].split(',') if len(sys.argv) > 1 else ['T33', 'T32', 'T42', 'T31', 'L3']
qs = [Fr(1, 4), Fr(9, 4), Fr(1, 100)]
for f in fams:
    j, m = K.fam(f)
    T = tuple(range(m))
    for (a, b) in K.rep_configs(k, m):
        for qv in qs:
            s = K.config_s(k, m, a, b, qv, Fr(1))
            w = [K.rsqrt(x) for x in s]
            assert all(x * x == y for x, y in zip(w, s)), (w, s)
            states, M = K.build_j1(k, j, T, w)
            vecs, dels = K.known_kernel(k, j, T, states, s, Fr(0), Fr(1))
            ann = all(all(x == 0 for x in K.matvec(M, v)) for v in vecs)
            nul = len(states) - K.rank_q(M)
            keep = [i for i in range(len(states)) if i not in dels]
            pd, bad = K.ldl_pd(K.principal(M, keep))
            print(f'{f} j={j} T={T} config {(a, b)} q={qv}: states {len(states)}, known {len(vecs)}, annihilated {ann}, '
                  f'nullity {nul}, gauge-fixed PD {pd} [{time.time() - t0:.1f}s]', flush=True)
