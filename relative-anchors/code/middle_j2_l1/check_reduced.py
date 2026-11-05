"""The certified objects against the block built by j1lib.py (Lemma B, j = 2), exact, at rational points.

For every size pair (u0, u1), u0 >= u1, u0 + u1 = k in the given list: the gauge-fixed Schur complements S certified
by cert_fixed.py (loaded from the files logs/cert_k*_*_*.pkl written by its runs in run_all.sh) are evaluated at rational q0 = w^2 and compared with the
projections W^T (Y_PP - Y_PQ Y_QQ^{-1} Y_QP) W of the rescaled j1lib block at w (exact inverse of the vacuum block,
no Woodbury); the kernel test vectors e.K are annihilated by that Schur complement; every pivot is positive.
Spectral check (floating point, numpy): at one point per k the spectrum of the Schur complement onto the P-states
equals the union over the parts of the eigenvalues of G^{-1} R (R, G on the basis e.K + e.W, each with multiplicity the
irrep dimension).  Usage: python check_reduced.py K1,K2,...  [spectral=K]"""
import sys, time
from fractions import Fraction as Fr
import numpy as np
import flint
from j1lib import J1, mu
from r2lib import states_L1
from fixedred import leading_pivots
from young import part_data, irreps
from saved import load_saved

t0 = time.time()
ks = [int(x) for x in sys.argv[1].split(',')]
opts = dict(a.split('=', 1) for a in sys.argv[2:] if '=' in a)
spec_ks = [int(x) for x in opts.get('spectral', '').split(',') if x]
ok = True
WPTS = [Fr(1, 3), Fr(2), Fr(7, 5)]


def j1lib_schur(u0, u1, w1):
    """the Schur complement Y_PP - Y_PQ Y_QQ^{-1} Y_QP of the rescaled j1lib block at w = (1,..,1, w1,..,w1), exact
    (flint fmpq_mat)"""
    k = u0 + u1
    w = [Fr(1)] * u0 + [w1] * u1
    st = states_L1(k)
    M = J1(w).matrix(st)
    mus = [mu(x, w) for x in st]
    n = len(st)
    nP = k * (k * (k - 1) // 2)
    nQ = n - nP
    f = lambda x: flint.fmpq(x.numerator, x.denominator)
    Y = [[M[i][j] * mus[i] * mus[j] for j in range(n)] for i in range(n)]
    YPP = flint.fmpq_mat(nP, nP, [f(Y[i][j]) for i in range(nP) for j in range(nP)])
    YPQ = flint.fmpq_mat(nP, nQ, [f(Y[i][j]) for i in range(nP) for j in range(nP, n)])
    YQQ = flint.fmpq_mat(nQ, nQ, [f(Y[i][j]) for i in range(nP, n) for j in range(nP, n)])
    return YPP - YPQ * YQQ.solve(YPQ.transpose())


def proj(S, A, B):
    """A^T S B for the fmpq_mat S and lists A, B of integer vectors, as lists of Fractions"""
    n = S.nrows()
    Am = flint.fmpq_mat(len(A), n, [flint.fmpq(x) for a in A for x in a])
    Bm = flint.fmpq_mat(n, len(B), [flint.fmpq(b[i]) for i in range(n) for b in B])
    R = Am * S * Bm
    return [[Fr(int(R[i, j].p), int(R[i, j].q)) for j in range(len(B))] for i in range(len(A))]


for k in ks:
    for u1 in range(0, k // 2 + 1):
        u0 = k - u1
        saved = load_saved(k, u0, u1)
        pts = WPTS if k <= 6 else WPTS[:1]
        if u1 == 0:
            pts = [Fr(1)]
        for w1 in pts:
            q0 = w1 * w1
            Sref = j1lib_schur(u0, u1, w1)
            good = True
            nparts = 0
            for (l0, l1) in irreps(u0, u1):
                pd = part_data(k, u0, l0, l1)
                if pd['mP'] == 0 and not pd['Q']:
                    continue
                nparts += 1
                Ssym = [[f(q0) for f in row] for row in saved.get((l0, l1), [])]
                Sp = proj(Sref, pd['W'], pd['W']) if pd['W'] else []
                good &= Sp == Ssym
                if pd['K']:
                    KK = proj(Sref, pd['K'] + pd['W'], pd['K'])
                    good &= all(x == 0 for row in KK for x in row)
                piv = leading_pivots(Ssym, Fr(1)) if Ssym else []
                good &= all(p is not None and p > 0 for p in piv)
            good &= nparts == len(saved)
            ok &= good
            print(f'k={k} (u0, u1) = ({u0}, {u1}), q = {q0}: the certified S of all {nparts} parts = projections of the '
                  f'j1lib Schur complement, e.K annihilated, pivots > 0: {good} [{time.time() - t0:.1f}s]', flush=True)
        if k in spec_ks and u1 > 0:
            w1 = Fr(1, 2)
            Sref = j1lib_schur(u0, u1, w1)
            ev = np.sort(np.linalg.eigvalsh(np.array([[int(Sref[i, j].p) / int(Sref[i, j].q) for j in range(Sref.ncols())]
                                                   for i in range(Sref.nrows())])))
            union = []
            for (l0, l1) in irreps(u0, u1):
                pd = part_data(k, u0, l0, l1)
                Bv = pd['K'] + pd['W']
                if not Bv:
                    continue
                R = np.array([[float(x) for x in row] for row in proj(Sref, Bv, Bv)])
                G = np.array([[float(sum(a * b for a, b in zip(x, y))) for y in Bv] for x in Bv])
                e = np.linalg.eigvals(np.linalg.solve(G, R)).real
                union += list(e) * pd['dim']
            union = np.sort(np.array(union))
            same = len(union) == len(ev) and np.max(np.abs(union - ev)) < 1e-9 * max(1.0, np.max(np.abs(ev)))
            ok &= bool(same)
            print(f'k={k} (u0, u1) = ({u0}, {u1}), q = {w1 * w1}: spectrum of the Schur complement ({len(ev)} '
                  f'eigenvalues, {int(np.sum(np.abs(ev) < 1e-9 * np.max(np.abs(ev))))} zero) = union of the parts with '
                  f'multiplicities: {same} [{time.time() - t0:.1f}s]', flush=True)
print('PASS' if ok else 'FAIL')
