"""The reduced forms of tlib.py against the Lemma B block (fastj1.FastJ1, exact), and the completeness of the reduction.

python check_explicit.py KMIN KMAX JMIN [W0S] [JLIST]   (JLIST: only these j, comma-separated)
For every k in [KMIN, KMAX], every j in [JMIN, k - 1], every tau = 1..j and both cases (marked copy in T or in U):
 (a) for each weight w0 in W0S (sigma = w0^2) and B scaled by 1 and by 1001/1000: every entry of the reduced form of
     tlib.reduced (numeric mode, every mode and i) equals v_x^T Y~ v_y with Y~ the rescaled block of FastJ1 and v_x the
     explicit twisted vectors of explicit_t.py; the norms agree; no explicit twisted vector is missing in tlib;
 (b) completeness (first w0 only): the class sum K_U of the transpositions of U' maps each twisted space V_(mode,i)
     into itself, with eigenvalues only the content sums of (u' - i, i) and (u' - i - 1, i, 1) (multiplicities m by
     exact ranks); the class sum K_T of T' acts on V_(sym,i) by the content sum of (tau') and on V_(anti,i) by that of
     (tau' - 1, 1); and sum over modes, i and the two irreps of dim(T'-irrep) dim(U'-irrep) m equals the number of
     states of the block;
 (c) the known kernel vectors p (x) y lie in the kernel of the explicit block (first w0, B unscaled).
"""
import sys, time
from fractions import Fraction as Fr
import flint
import tlib
from explicit_t import (block, block_states, twisted, form, class_sum, content_sum, dim_irrep, layout, known_kernel)


def tofr(x):
    return Fr(int(x.p), int(x.q))


def rank_fr(rows):
    if not rows:
        return 0
    n, m = len(rows), len(rows[0])
    M = flint.fmpq_mat(n, m, [flint.fmpq(x.numerator, x.denominator) for r in rows for x in r])
    return M.rank()


def main():
    kmin, kmax, jmin = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    w0s = [Fr(x) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 else [Fr(2), Fr(1, 3)]
    jlist = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else None
    T0 = time.time()
    ok = True
    nent = nforms = nblocks = 0
    for k in range(kmin, kmax + 1):
        for j in (range(jmin, k) if jlist is None else [x for x in jlist if jmin <= x < k]):
            for tau in range(1, j + 1):
                for case in ('MT', 'MU'):
                    if case == 'MU' and tau > k - 1:
                        continue
                    T, Tp, Up, pairs_, R = layout(k, tau, case, 0)
                    st = block_states(k, j, T)
                    if not st:
                        continue
                    nblocks += 1
                    total = 0
                    for wn, w0 in enumerate(w0s):
                        for bs in (Fr(1), Fr(1001, 1000)):
                            st2, Y = block(k, j, T, w0, bscale=bs)
                            assert st2 == st
                            if wn == 0 and bs == 1:
                                for kv in known_kernel(k, j, T, st, w0):
                                    for x in range(len(st)):
                                        if sum(Y[x][c] * v for c, v in kv.items()) != 0:
                                            ok = False
                                            print('KERNEL VECTOR NOT ANNIHILATED', k, j, tau, case, flush=True)
                                            break
                            for mode in ('sym', 'anti'):
                                for i in tlib.i_range(k, tau, case):
                                    C = tlib.Calc(k, j, tau, i, case, mode, sigma=w0 * w0, bscale=bs)
                                    if not C.valid():
                                        continue
                                    names, norms, Q = tlib.reduced(C)
                                    vecs = twisted(k, j, tau, case, mode, i, st)
                                    if set(vecs) != set(names):
                                        ok = False
                                        print('NAME SETS DIFFER', k, j, tau, case, mode, i, sorted(vecs), names,
                                              flush=True)
                                        continue
                                    if wn == 0 and bs == 1:
                                        nforms += 1
                                    for x, na in enumerate(names):
                                        if Fr(norms[x]) != sum(c * c for c in vecs[na].values()):
                                            ok = False
                                            print('NORM DIFFERS', k, j, tau, case, mode, i, na, flush=True)
                                        for y, nb in enumerate(names):
                                            ex = form(Y, vecs[na], vecs[nb])
                                            nent += 1
                                            if tofr(Q[x][y]) != ex:
                                                ok = False
                                                print('ENTRY DIFFERS', k, j, tau, case, mode, i, na, nb, str(w0), str(bs),
                                                      tofr(Q[x][y]), ex, flush=True)
                                    # (b) completeness, first point only
                                    if wn == 0 and bs == 1 and names:
                                        up = len(Up)
                                        tp = len(Tp)
                                        nmv = [vecs[nm] for nm in names]
                                        nrm = [sum(c * c for c in v.values()) for v in nmv]
                                        lamU = [(up - i, i), (up - i - 1, i, 1)]
                                        lamT = (tp,) if mode == 'sym' else (tp - 1, 1)
                                        Kmat = []
                                        layout_i = layout(k, tau, case, i)
                                        for v in nmv:
                                            Kv = class_sum(st, layout_i[2], v)
                                            col = [Fr(sum(c * Kv.get(xx, 0) for xx, c in w.items())) / nn
                                                   for w, nn in zip(nmv, nrm)]
                                            rec = {}
                                            for cc, w in zip(col, nmv):
                                                for xx, c in w.items():
                                                    rec[xx] = rec.get(xx, 0) + cc * c
                                            rec = {xx: c for xx, c in rec.items() if c != 0}
                                            if rec != Kv:
                                                ok = False
                                                print('K_U DOES NOT PRESERVE V', k, j, tau, case, mode, i, flush=True)
                                            Kmat.append(col)
                                            KT = class_sum(st, Tp, v)
                                            cT = content_sum(lamT)
                                            if KT != {xx: c * cT for xx, c in v.items() if c * cT != 0}:
                                                ok = False
                                                print('K_T NOT SCALAR', k, j, tau, case, mode, i, flush=True)
                                        nV = len(names)
                                        # Kmat[x] = column x (coefficients of K v_x)
                                        A = [[Kmat[x][y] for x in range(nV)] for y in range(nV)]
                                        ms = []
                                        for lam in lamU:
                                            if min(lam) < 0 or (len(lam) == 3 and lam[1] < 1) or lam[0] < lam[1]:
                                                ms.append(0)
                                                continue
                                            c = content_sum(lam)
                                            Mc = [[A[a][b] - (c if a == b else 0) for b in range(nV)] for a in range(nV)]
                                            ms.append(nV - rank_fr(Mc))
                                        if sum(ms) != nV:
                                            ok = False
                                            print('EIGENVALUES OUTSIDE THE TWO IRREPS', k, j, tau, case, mode, i, ms,
                                                  nV, flush=True)
                                        for lam, mm in zip(lamU, ms):
                                            if mm:
                                                total += dim_irrep(lamT) * dim_irrep(lam) * mm
                    if total != len(st):
                        ok = False
                        print('DIMENSION COUNT FAILS', k, j, tau, case, total, len(st), flush=True)
            print(f'k={k} j={j} done [{time.time() - T0:.1f}s]', flush=True)
    print(f'{nblocks} blocks, {nforms} reduced forms, {nent} entries compared')
    print('PASS' if ok else 'FAIL', f'{time.time() - T0:.0f}s')


if __name__ == '__main__':
    main()
