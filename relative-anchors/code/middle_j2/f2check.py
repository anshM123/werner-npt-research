"""Checks of the determinant route of f2det.py (not needed by the proof).
 points: the interpolant N~ (recomputed from the saved slices, as in f2det.py interp) at points (q, x0, x1) with x1
         outside the grid 0..B equals det S' * D_c computed directly at that point (every variable specialised, an
         11 x 11 rational determinant);
 neg:    with B replaced by (1 + 1/1000) B the reduced form of the part changes, the kernel vectors are no longer
         annihilated, and det S' on the slice x1 = 0 has negative coefficients (so the sign step cannot pass).
python f2check.py points | neg"""
import sys, os, pickle, time
from fractions import Fraction as Fr
import flint
import f2det
from f2det import Field, FF

mode = sys.argv[1]
t0 = time.time()
if mode == 'points':
    with open(os.path.join(f2det.OUT, 'setup.pkl'), 'rb') as fh:
        st = pickle.load(fh)
    B, dR = st['bounds']['x1']
    L = B - dR
    vals, mons = [], set()
    for c in range(B + 1):
        with open(os.path.join(f2det.OUT, f's_{c}.pkl'), 'rb') as fh:
            d = pickle.load(fh)
        vals.append(d)
        mons.update(d.keys())
    mons = sorted(mons)
    V = flint.fmpq_mat(B + 1, B + 1, [flint.fmpq(c) ** l for c in range(B + 1) for l in range(B + 1)]).inv()
    Y = flint.fmpq_mat(B + 1, len(mons))
    for c, d in enumerate(vals):
        for j, mo in enumerate(mons):
            if mo in d:
                Y[c, j] = flint.fmpq(*d[mo])
    Cf = V * Y
    K, S, ok, facs, estar = f2det.setup(log=lambda m: None)
    Dc = K.ctx.from_dict({(0, 0, 0): 1})
    for f, e in zip(facs, st['dexp']):
        Dc = Dc * f ** e
    good = True
    for (qv, x0v, x1v) in ((Fr(2, 7), Fr(0), Fr(300)), (Fr(5), Fr(40), Fr(1000)), (Fr(1, 1000), Fr(3), Fr(517, 2)),
                           (Fr(1000), Fr(250), Fr(10 ** 4))):
        pt = [flint.fmpq(v.numerator, v.denominator) for v in (qv, x0v, x1v)]
        nt = flint.fmpq(0)
        for j, mo in enumerate(mons):
            a, b = mo
            for l in range(L + 1):
                c = Cf[l, j]
                if c != 0:
                    nt += c * pt[0] ** a * pt[1] ** b * pt[2] ** l
        n = len(S)
        A = flint.fmpq_mat(n, n, [(S[i][j].n(*pt) / S[i][j].d(*pt)) if not S[i][j].n.is_zero() else flint.fmpq(0)
                                  for i in range(n) for j in range(n)])
        direct = A.det() * Dc(*pt)
        eq = direct == nt
        good &= eq
        mag = lambda v: f'{"-" if v < 0 else ""}10^{len(str(abs(int(v.p)))) - len(str(int(v.q)))} (approx.)'
        print(f'(q, x0, x1) = ({qv}, {x0v}, {x1v}): det S\' * D_c = {mag(direct)}, N~ = {mag(nt)}, '
              f'equal: {eq} [{time.time()-t0:.0f}s]', flush=True)
    print('PASS' if good else 'FAIL')
else:
    # the same construction as f2det.build with B scaled
    from f2sym import Model2, reduce_part2, assemble, part_spec2, keep_for2, pair_sums2
    from l1sym import det, minv, mm
    part, o0, o1 = 'inv', 3, 3
    K = Field(['q', 'x0', 'x1'])
    q, x0, x1 = K.gens
    one, zero = K(1), K(0)
    M = Model2(one, one, [one, q], [x0 + o0, x1 + o1], Fr(1001, 1000))
    red = reduce_part2(M, part, o0, o1)
    F, iR, iP, iQ = assemble(red)
    n = len(F)
    ann = True
    for kv in red['ker']:
        vec = [zero] * n
        for j, c in enumerate(kv):
            vec[iP[j]] = K(c)
        ann &= all(sum((F[i][j] * vec[j] for j in range(n)), zero) == 0 for i in range(n))
    print(f'kernel vectors annihilated with B scaled: {ann}', flush=True)
    # the gauge-fixed form with the vacuum states eliminated by the plain Schur complement (no identity for R_QQ)
    iPg = [iP[j] for j in range(len(iP)) if j not in red['gauge']]
    X = iPg + iR
    FQQ = [[F[a][b] for b in iQ] for a in iQ]
    FQQi = minv(FQQ, one)
    FXQ = [[F[a][b] for b in iQ] for a in X]
    T = mm(FXQ, FQQi, zero)
    S = [[F[X[a]][X[b]] - sum((T[a][l] * F[iQ[l]][X[b]] for l in range(len(iQ))), zero) for b in range(len(X))]
         for a in range(len(X))]
    K2 = Field(['q', 'x0'])
    dt = f2det.det_slice(K2, S, 2, 0)
    co = list(dt.n.to_dict().values())
    npos, nneg = sum(1 for v in co if v > 0), sum(1 for v in co if v < 0)
    print(f'det S\' at x1 = 0 with B scaled: {npos} positive, {nneg} negative numerator coefficients '
          f'[{time.time()-t0:.0f}s]', flush=True)
    print('PASS' if (ann and nneg == 0) else 'FAIL')
