"""Exact PSD certificate of every LETTER-2 torus block of (R_inf)_j at equal weights (j = 0..k), via the
S_T x S_U reduction of eqlib2 (orbit-reduced sums).  Usage: python run2.py kmin kmax [check_full]"""
import sys, time
from eqlib2 import caseI_vectors, caseII_vectors, form, row_reps
from eqlib import ldl_psd
kmin, kmax = int(sys.argv[1]), int(sys.argv[2])
full = len(sys.argv) > 3
ok = True
T0 = time.time()
for k in range(kmin, kmax + 1):
    nblk = 0; nul = {}
    for j in range(0, k + 1):
        jobs = []
        for tau in range(1, min(k, j + 1) + 1):
            for i in range(0, (k - tau) // 2 + 1):
                jobs.append(('I', tau, i, 'sym'))
                if tau >= 2:
                    jobs.append(('I', tau, i, 'anti'))
        for tau1 in range(0, k):
            for i in range(0, (k - 1 - tau1) // 2 + 1):
                jobs.append(('II', tau1, i, None))
        for (cs, tau, i, mode) in jobs:
            vecs = caseI_vectors(k, j, tau, i, mode) if cs == 'I' else caseII_vectors(k, j, tau, i)
            if not vecs:
                continue
            reps = row_reps(k, tau, i, mode) if cs == 'I' else None
            G, Q = form(k, j, vecs, None if full else reps)
            if full is False and cs == 'I' and k <= 5:
                G2, Q2 = form(k, j, vecs, None)
                if G2 != G or Q2 != Q:
                    ok = False; print('REDUCTION MISMATCH', k, j, cs, tau, i, mode)
            psd, rk, piv = ldl_psd(Q)
            nblk += 1
            nn = len(vecs) - rk
            nul[nn] = nul.get(nn, 0) + 1
            if not psd:
                ok = False; print('NOT PSD', k, j, cs, tau, i, mode, [v[0] for v in vecs])
    print(f'k={k}: {nblk} reduced letter-2 forms, all PSD so far: {ok}; nullity histogram {nul}  ({time.time()-T0:.0f}s)', flush=True)
print('PASS' if ok else 'FAIL')
