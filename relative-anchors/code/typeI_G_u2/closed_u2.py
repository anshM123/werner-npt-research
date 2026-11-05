"""Closed forms of the T-blocks of G for u = 2 (exact, python-flint).

The type I block of (R_inf)_{k-1} with T = {0, 1, 2}, U = {3, 4} (cf.TypeI, units of 2Y) is built over the rational
function field Q(w0, ..., w4).  All P_{r,h} except the gauge row r0 = 3 (the kernel p (x) f_h has P_{3,h}-coordinate
w_3 != 0) and Q are eliminated by exact Gaussian elimination; the result G is the Schur complement on
(rho_t, zeta_{t,3}, zeta_{t,4}), t = 0, 1, 2.  In the coordinates zeta_{t,h} = w_t w_h z_{t,h} every entry is even in
every w_i (checked), hence a rational function of s_i = w_i^2.  Output (closed_u2.pkl): the 3 x 3 blocks G_00 and G_01
as (numerator, denominator) dicts over (s0, s1, s2, s3, s4)."""
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from flintfield import Field
from cf import TypeI


def halve(p):
    d = p.to_dict()
    out = {}
    for m, c in d.items():
        assert all(e % 2 == 0 for e in m), 'odd exponent'
        out[tuple(e // 2 for e in m)] = c
    return out


def schur_T(K, w, T, U_first_gauge=True):
    ti = TypeI(w, T, F=K)
    M = ti.build()
    U = ti.U; r0 = U[0]
    tb = []; scale = []
    for t in ti.T:
        tb.append(ti.idx[('R', t)]); scale.append(K(1))
        for h in U:
            tb.append(ti.idx[('Z', t, h)]); scale.append(w[t] * w[h])
    el = [i for i in range(ti.n) if i not in tb and not (ti.labels[i][0] == 'P' and ti.labels[i][1] == r0)]
    order = el + tb
    A = [[M[i][j] for j in order] for i in order]
    ne = len(el)
    for p in range(ne):
        pv = A[p][p]
        assert pv != 0
        for i in range(p + 1, len(order)):
            if A[i][p] == 0:
                continue
            f = A[i][p] / pv
            for j in range(p + 1, len(order)):
                if A[p][j] != 0:
                    A[i][j] = A[i][j] - f * A[p][j]
    nt = len(tb)
    return [[A[ne + i][ne + j] * scale[i] * scale[j] for j in range(nt)] for i in range(nt)]


if __name__ == '__main__':
    t0 = time.time()
    K = Field(['w0', 'w1', 'w2', 'w3', 'w4'])
    w = list(K.gens)
    G = schur_T(K, w, [0, 1, 2])
    print(f'Schur complement done ({time.time() - t0:.1f} s)', flush=True)
    out = {}
    for name, (bi, bj) in (('G00', (0, 0)), ('G01', (0, 1))):
        blk = []
        for i in range(3):
            row = []
            for j in range(3):
                x = G[3 * bi + i][3 * bj + j]
                row.append((halve(x.n), halve(x.d)))
            blk.append(row)
        out[name] = blk
    # sanity: symmetry of G00 and the transpose relation G10 = G01^T
    ok = all(G[i][j] == G[j][i] for i in range(9) for j in range(9))
    print('G symmetric:', ok, flush=True)
    pickle.dump(out, open('closed_u2.pkl', 'wb'))
    for name in ('G00', 'G01'):
        for i in range(3):
            for j in range(3):
                n_, d_ = out[name][i][j]
                print(name, i, j, 'num terms', len(n_), 'den terms', len(d_),
                      'num deg', max(sum(m) for m in n_) if n_ else None, 'den deg', max(sum(m) for m in d_), flush=True)
    print(f'total {time.time() - t0:.1f} s')
    print('PASS' if ok else 'FAIL')
