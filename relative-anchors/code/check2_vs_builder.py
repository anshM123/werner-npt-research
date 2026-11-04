"""eqlib2 (letter-2 blocks, equal weights) vs the independent builder ylib.py: k = 3, 4, every j,
every case I pattern size tau (sym/anti), every case II, every i:  Q == 2^(k-j) v^T Y v and G == v^T v."""
import sys, time
import ylib
from eqlib2 import caseI_vectors, caseII_vectors, form, fq
ok = True; ncmp = 0
for k in (3, 4):
    ring = ylib.point_ring([1] * k)
    Bd = ylib.Builder(ring)
    for j in range(0, k + 1):
        states = ylib.rinf_states(k, j)
        idx = {}
        for n_, (a, y) in enumerate(states):
            if sum(a) == 0:
                idx[('vac', y)] = n_
            else:
                r = next(t for t in range(k) if a[t] != 0)
                idx[((r, a[r]), y)] = n_
        cache = {}
        def Yent(r, c):
            if (r, c) not in cache:
                cache[(r, c)] = Bd.Y(states[r], states[c])
            return cache[(r, c)]
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
            G, Q = form(k, j, vecs)
            full = []
            for (_, psi, phi) in vecs:
                v = {}
                for a, pv in psi.items():
                    for y, c in pv.items():
                        v[idx[(a, y)]] = v.get(idx[(a, y)], 0) + c
                for y, c in phi.items():
                    v[idx[('vac', y)]] = v.get(idx[('vac', y)], 0) + c
                full.append(v)
            for x in range(len(vecs)):
                for z in range(len(vecs)):
                    qb = sum((full[x][r] * Yent(r, c) * full[z][c] for r in full[x] for c in full[z]), fq(0)) * 2 ** (k - j)
                    gb = sum((full[x][r] * full[z].get(r, 0) for r in full[x]), fq(0))
                    ncmp += 1
                    if qb != Q[x][z] or gb != G[x][z]:
                        ok = False
                        print('MISMATCH', k, j, cs, tau, i, mode, vecs[x][0], vecs[z][0], qb, Q[x][z], gb, G[x][z])
        print(f'k={k} j={j} done', flush=True)
print(ncmp, 'entries compared;', 'PASS' if ok else 'FAIL')
