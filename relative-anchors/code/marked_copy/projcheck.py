"""The reduced forms of the marked-copy reduction (mlib.py, numeric mode) against the projections of the true block,
at sizes beyond those of the original checks.  The block is built by the fast builder of markcheck.py (checked
entry by entry against j1lib.J1.matrix, i.e. Lemma B, for k <= 6), rescaled as in explicit.block (D M D with
D = diag(mu), mu from j1lib), and projected onto the explicit twisted vectors of explicit.twisted (sums of states with
signs +-1 built from the states).  Every entry of Q_i must equal v_a^T Y~ v_b, and the norms must agree.

  python projcheck.py K W0,W0,...      every j = -1..K and every i, at sigma = W0^2 for each W0
"""
import sys, os, time
from fractions import Fraction as Fr
import flint

import mlib, explicit
from markcheck import Fast, fq
from j1lib import mu

T0 = time.time()


def main():
    k = int(sys.argv[1])
    ws = [Fr(x) for x in sys.argv[2].split(',')]
    ok = True
    cnt = 0
    for w0 in ws:
        sig = w0 * w0
        w = [w0] + [Fr(1)] * (k - 1)
        fb = Fast(w)
        for j in range(-1, k + 1):
            t0 = time.time()
            st, M = fb.block(j)
            n = len(st)
            mus = [mu(s, w) for s in st]
            Y = flint.fmpq_mat(n, n, [fq(M[a][b] * mus[a] * mus[b]) for a in range(n) for b in range(n)])
            del M
            pos = {s: p for p, s in enumerate(st)}
            stE = explicit.states(k, j)
            assert len(stE) == n and set(stE) == set(st)
            for i in range(0, (k - 1) // 2 + 1):
                names, norms, Q = mlib.reduced(mlib.Calc(k, j, i, sigma=sig))
                st2, vecs = explicit.twisted(k, j, i)
                if sorted(names) != sorted(vecs):
                    ok = False
                    print('NAMES DIFFER', k, j, i, names, sorted(vecs), flush=True)
                    continue
                m = len(names)
                if m == 0:
                    continue
                V = flint.fmpq_mat(n, m)
                for c, nm in enumerate(names):
                    for x, cf in vecs[nm].items():
                        V[pos[st2[x]], c] = fq(cf)
                P = V.transpose() * (Y * V)
                for a, na in enumerate(names):
                    if fq(sum(cf * cf for cf in vecs[na].values())) != flint.fmpq(norms[a]):
                        ok = False
                        print('NORM DIFFERS', k, j, i, na, flush=True)
                    for b in range(m):
                        cnt += 1
                        if P[a, b] != Q[a][b]:
                            ok = False
                            print('ENTRY DIFFERS', k, j, i, sig, na, names[b], flush=True)
            print(f'k={k} sigma={sig} j={j}: {n} states, all i compared [{time.time() - t0:.1f}s]', flush=True)
    print(f'{cnt} entries of the reduced forms compared with projections of the Lemma B block', flush=True)
    print(f'total {time.time() - T0:.0f}s', flush=True)
    print('PASS' if ok else 'FAIL', flush=True)


if __name__ == '__main__':
    main()
