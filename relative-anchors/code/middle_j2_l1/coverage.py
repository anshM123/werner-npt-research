"""Coverage of the fixed-k statement by the certificate runs of run_all.sh (bookkeeping).

A weight vector with at most two distinct weights on k copies has m copies of squared weight alpha and k - m of beta
(0 <= m <= k).  Permuting the copies permutes the states (the block is permutation-equivalent), and dividing all
squared weights by a positive c multiplies the rescaled block by c^{-3} (homogeneity, checked in check_builder.py), which
keeps inertia and kernel.  So with u0 = max(m, k - m), u1 = k - u0, the block is that of the certified pair (u0, u1) at
q = (squared weight of the smaller level)/(squared weight of the larger level) > 0 (any q if alpha = beta, in which case
u1 = 0 also works).  For every k in [5, KMAX] and every m this script finds the certified pair and checks that its log
logs/run_cert_k{k}_{u0}_{u1}.log ends with PASS and reports every part positive.  Usage: python coverage.py KMAX"""
import sys, os

kmax = int(sys.argv[1])
ok = True
here = os.path.dirname(os.path.abspath(__file__))
for k in range(5, kmax + 1):
    pairs = set()
    for m in range(0, k + 1):
        u0 = max(m, k - m)
        pairs.add((u0, k - u0))
    for (u0, u1) in sorted(pairs, reverse=True):
        fn = os.path.join(here, 'logs', f'run_cert_k{k}_{u0}_{u1}.log')
        if not os.path.exists(fn):
            ok = False
            print(f'k={k} ({u0}, {u1}): log missing')
            continue
        lines = open(fn).read().splitlines()
        nparts = sum(1 for ln in lines if ln.strip().startswith('part ') and 'gauge-fixed S of size' in ln)
        allpos = all(ln.rstrip().split('every leading minor > 0 on (0, inf): ')[1].startswith('True')
                     for ln in lines if 'every leading minor > 0 on (0, inf): ' in ln)
        last = lines[-1].strip() if lines else ''
        good = last == 'PASS' and allpos and nparts > 0
        ok &= good
        print(f'k={k} (u0, u1) = ({u0}, {u1}): {nparts} parts, all positive: {allpos}, log ends with {last}: '
              f'{"covered" if good else "NOT covered"}')
    print(f'k={k}: every m in 0..{k} maps to one of {sorted(pairs, reverse=True)}')
print('PASS' if ok else 'FAIL')
