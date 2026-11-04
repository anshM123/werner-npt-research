"""Check the reduced matrices of eqlib against the independent exact builder (ylib.py, which is
checked against a dense floating-point construction of M, V_p, D, R): for k = 3, 4, 5, every j, every i,
  Q_mine[al][be] == 2^(k-j) * v_al^T Y_builder v_be   and   G_mine == v_al^T v_be,
and the letter-1 block of the builder is PSD iff all Q_i are, with matching kernel dimension:
  dim ker(block) = sum_i [ nullity of Q_i restricted to irreps ] * dim S^lambda  (checked via ranks below).
"""
import sys, time
import flint
import ylib
from eqlib import reduced, fq, ldl_psd

ok_all = True
for k in (3, 4, 5):
    ring = ylib.point_ring([1] * k)
    B = ylib.Builder(ring, letters=(0, 1))
    for j in range(-1, k + 1):
        states = ylib.rinf_states(k, j, letters=(0, 1))
        if not states:
            continue
        t0 = time.time()
        Ymat = B.matrix(states)
        idx = {}
        for n_, (a, y) in enumerate(states):
            ym = sum(1 << t for t in range(k) if y[t] == 1)
            if sum(a) == 0:
                idx[('vac', ym)] = n_
            else:
                r = a.index(1)
                idx[(r, ym)] = n_
        Mfl = flint.fmpq_mat([[x for x in row] for row in Ymat])
        rank_full = Mfl.rank()
        nul_full = len(states) - rank_full
        for i in range(0, k // 2 + 1):
            if j < 0:
                continue
            names, G, Q, vecs = reduced(k, j, i)
            # full vectors
            full = []
            for (_, psi, phi) in vecs:
                v = [fq(0)] * len(states)
                for a, pv in psi.items():
                    for ym, c in pv.items():
                        v[idx[(a, ym)]] += c
                for ym, c in phi.items():
                    v[idx[('vac', ym)]] += c
                full.append(v)
            n = len(full)
            for al in range(n):
                Yv = [sum((Ymat[r][c] * full[al][c] for c in range(len(states)) if full[al][c] != 0), fq(0)) for r in range(len(states))]
                for be in range(n):
                    qb = sum((Yv[r] * full[be][r] for r in range(len(states))), fq(0)) * 2 ** (k - j)
                    gb = sum((full[al][r] * full[be][r] for r in range(len(states))), fq(0))
                    if qb != Q[al][be] or gb != G[al][be]:
                        ok_all = False
                        print('MISMATCH', k, j, i, names[al], names[be], qb, Q[al][be], gb, G[al][be])
            psd, rk, piv = ldl_psd(Q)
            print(f'k={k} j={j} i={i} vecs={names} Q psd={psd} rank={rk}/{n}')
        print(f'  k={k} j={j}: full letter-1 block size {len(states)}, nullity {nul_full}  ({time.time()-t0:.1f}s)')
print('PASS' if ok_all else 'FAIL')
