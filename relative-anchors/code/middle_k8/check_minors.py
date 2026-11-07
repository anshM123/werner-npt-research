"""The certified leading minors (saved by cert5.py) against direct minors: for every family and configuration, at
rational squares q0, every saved leading principal minor of the gauge-fixed block equals the determinant of the
corresponding leading principal submatrix of the gauge-fixed block built by j1lib.py (Lemma B) at q0 (exact, flint).
python check_minors.py FAMILIES [k=5]
Prints PASS or FAIL on the last line."""
import sys, time
from fractions import Fraction as Fr
import k5lib as K
from saved import load, minors_at, cert_path

t0 = time.time()
args = [a for a in sys.argv[1:] if '=' not in a]
opts = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
k = int(opts.get('k', 5))
fams = args[0].split(',')
ok = True
for f in fams:
    j, m = K.fam(f)
    T = tuple(range(m))
    cnt = 0
    for (a, b) in (K.rep_configs(k, m) if opts.get('reps') == '1' else K.configs(k, m)):
        data = load(cert_path(f, a, b, k))
        keep = data['keep']
        for q0 in (Fr(4, 9), Fr(49, 16)):
            s = K.config_s(k, m, a, b, q0, Fr(1))
            w = [K.rsqrt(x) for x in s]
            st, Mj = K.build_j1(k, j, T, w)
            G = K.principal(Mj, keep)
            mins = minors_at(data, q0)
            if len(mins) != len(G):
                ok = False
                print(f'  {f} {(a, b)}: {len(mins)} saved minors for a gauge-fixed block of size {len(G)}')
                continue
            for r in range(1, len(G) + 1):
                dv = K.to_fmpq_mat([row[:r] for row in G[:r]]).det()
                eq = Fr(int(dv.p), int(dv.q)) == mins[r - 1]
                ok &= eq
                cnt += 1
                if not eq:
                    print(f'  MISMATCH {f} {(a, b)} q = {q0} minor {r}')
    print(f'{f}: {cnt} saved minors equal direct minors of the Lemma B block: {ok} [{time.time() - t0:.1f}s]', flush=True)
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
