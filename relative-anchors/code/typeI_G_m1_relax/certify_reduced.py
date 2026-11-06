"""Exact pointwise certificates through the reduced finite system (NOTES.md, Section 8), for larger u.

Lemma F (reduction): write Psi = sum_h v_h^T D_h v_h + y^T C y with y = (rho, F_1..F_30), F = sum_h U_h^T v_h, and
the constraint y_j0 = 0 (a orthogonal to w).  If every D_h > 0 and B_F = sum_h U_h^T D_h^{-1} U_h > 0, then
Psi > 0 on the constraint subspace iff the matrix N = C + (0 (+) B_F^{-1}) with row and column j0 deleted is > 0.

usage: python certify_reduced.py u_list n_list seed [mode]
mode 'check' (default): conditions (i), (ii); mode 'both' also certifies Psi directly (psi_exact) and compares;
mode 'controlB': replaces B_F by 1001/1000 B_F in condition (ii) (a negative control)."""
import sys, time, random, math
import flint
from fractions import Fraction as Fr
from reduced_c import ReducedC as Reduced, to_q
from certify_points import leading_minors_bareiss, family

Q = flint.fmpq


def conditions(R, wt, wU, scaleB=Fr(1)):
    num = R.numeric(Fr(wt.p, wt.q) if isinstance(wt, flint.fmpq) else wt,
                    [Fr(x.p, x.q) if isinstance(x, flint.fmpq) else x for x in wU])
    u = num['u']
    J = len(num['C'])
    # (i) every D_h > 0
    ok_i = True
    Dinv = []
    for h in range(u):
        Dm = flint.fmpq_mat(3, 3, [to_q(num['D'][h][a][b]) for a in range(3) for b in range(3)])
        m1 = Dm[0, 0]
        m2 = Dm[0, 0] * Dm[1, 1] - Dm[0, 1] * Dm[1, 0]
        m3 = Dm.det()
        if not (m1 > 0 and m2 > 0 and m3 > 0):
            ok_i = False
            return dict(ok_i=False, ok_B=None, ok_ii=None)
        Dinv.append(Dm.inv())
    # B_F
    nF = J - 1
    Uh = [flint.fmpq_mat(3, nF) for _ in range(u)]
    for j, (blk, vals) in enumerate(num['U']):
        for h in range(u):
            Uh[h][blk, j] = to_q(vals[h])
    BF = flint.fmpq_mat(nF, nF)
    for h in range(u):
        BF += Uh[h].transpose() * Dinv[h] * Uh[h]
    BF = BF * to_q(scaleB)
    mB, _ = leading_minors_bareiss(BF)
    ok_B = len(mB) == nF and all(x > 0 for x in mB)
    if not ok_B:
        return dict(ok_i=True, ok_B=False, ok_ii=None)
    BFi = BF.inv()
    Cm = flint.fmpq_mat(J, J, [to_q(num['C'][i][j]) for i in range(J) for j in range(J)])
    for i in range(nF):
        for j in range(nF):
            Cm[1 + i, 1 + j] += BFi[i, j]
    keep = [k for k in range(J) if k != R.j0]
    Nk = flint.fmpq_mat(len(keep), len(keep), [Cm[i, j] for i in keep for j in keep])
    mN, _ = leading_minors_bareiss(Nk)
    ok_ii = len(mN) == len(keep) and all(x > 0 for x in mN)
    return dict(ok_i=True, ok_B=True, ok_ii=ok_ii)


if __name__ == '__main__':
    us = [int(x) for x in sys.argv[1].split(',')]
    nl = [int(x) for x in sys.argv[2].split(',')]
    nl = nl * len(us) if len(nl) == 1 else nl
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    mode = sys.argv[4] if len(sys.argv) > 4 else 'check'
    rng = random.Random(seed)
    t0 = time.time()
    R = Reduced()
    print(f'reduced system loaded ({time.time() - t0:.1f} s)', flush=True)
    allok = True
    for u, npts in zip(us, nl):
        t1 = time.time()
        cnt = dict(i=0, B=0, ii=0, agree=0)
        for k in range(npts):
            wt, wU = family(u, k, rng)
            res = conditions(R, Fr(int(wt.p), int(wt.q)), [Fr(int(x.p), int(x.q)) for x in wU],
                             scaleB=Fr(1001, 1000) if mode == 'controlB' else Fr(1))
            cnt['i'] += bool(res['ok_i'])
            cnt['B'] += bool(res['ok_B'])
            cnt['ii'] += bool(res['ok_ii'])
            if mode == 'both':
                from psi_exact import PsiExact
                Pm, iy = PsiExact(wt, wU).matrix()
                mins, _ = leading_minors_bareiss(Pm)
                direct = len(mins) == Pm.nrows() and all(x > 0 for x in mins)
                cnt['agree'] += (direct == bool(res['ok_i'] and res['ok_ii'])) if res['ok_B'] else 0
            if mode == 'controlB' and not res['ok_ii']:
                print(f'  control: condition (ii) fails with B -> 1.001 B at family {k % 8}: w_t = {wt}, '
                      f'w_U = {wU}', flush=True)
            if res['ok_B'] is False:
                print(f'  B_F singular at family {k % 8} (route not applicable there)', flush=True)
            else:
                allok &= bool(res['ok_i'] and res['ok_ii'])
            if mode == 'both' and res['ok_B']:
                allok &= (direct == bool(res['ok_i'] and res['ok_ii']))
        extra = f", agreement with the direct certificate {cnt['agree']}/{npts}" if mode == 'both' else ''
        print(f"u = {u}: (i) D_h > 0 for all h: {cnt['i']}/{npts}; B_F > 0: {cnt['B']}/{npts}; (ii) N > 0: "
              f"{cnt['ii']}/{npts}{extra}; {time.time() - t1:.1f} s", flush=True)
    if mode == 'controlB':
        print('control finished (failures expected)')
    else:
        print('PASS' if allok else 'FAIL')
