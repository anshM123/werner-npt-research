"""Negative controls (exact).
 1. The certificate of prove_gersh.py (u1 = 1) fails with the constant moved too low (H = H_LO) and too high
    (H = H_HI): some condition then has coefficients of both signs.
 2. Explicit rational witnesses, computed from the exact G of the whole block (Lemma B builder, no closed form), that
    the inequalities themselves are false there: an (A) inequality at H = H_LO and a (B) inequality at H = H_HI, at
    two-valued weight vectors; at H = 1/2 every quantity there is positive.
 3. B -> (1 + 1/1000) B: at explicit two-valued configurations (among them the marked configuration of Theorem M1),
    the gauge-reduced block (P_{r0,h} removed) has a negative LDL^T pivot, while the true block has all pivots
    positive there; the scaled builder differs from the closed form (so the comparison of check_e2e2.py is not vacuous).
 4. Theorem M1: subtracting eps I (eps = 1/1000) from D - O and D + (sigma - 1) O, the certificate of prove_eqT.py
    fails (no uniform margin), and the perturbed matrix is indefinite at an explicit point.
python controls2.py"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import flint
from cf import TypeI
from check_e2e2 import direct, quantities
from emlib import qf, ldl_pivots

Q = flint.fmpq
H_LO, H_HI = Q(1, 10), Q(4, 5)


def cert_all(h):
    import prove_gersh as PG
    conds = PG.conditions('one', h)
    PG.COLLECT = None
    bad = [name for name, val, free in conds if not PG.cert(name, val, free)]
    return bad


def quant_at(wT, wU, u1, u2, h):
    import check_e2e2 as CE
    old = CE.H
    CE.H = Q(h)
    Gz, piv, Mb, lab, w, T, U = direct(wT, wU)
    s = qf(wU[0]) ** 2
    sT = [qf(x) ** 2 for x in wT]
    qn = quantities(Gz, sT, s, u1, u2)
    CE.H = old
    return qn


def minval(qn, kind):
    if kind == 'A':
        vals = [v for (_, _, v) in qn['A_inv'] + qn['A_std1'] + qn['A_std2']]
    else:
        vals = [min(p) for p in qn['B_inv']] + qn['B_std1'] + qn['B_std2']
    return min(vals)


def main():
    ok = True
    t0 = time.time()
    # 1. certificate fails for H too small / too large
    for h in (H_LO, H_HI):
        bad = cert_all(h)
        print('certificate (u1 = 1) with H = %s: %d conditions not certified %s' % (h, len(bad), bad[:4]), flush=True)
        ok &= len(bad) > 0
    # 2. witnesses (exact, whole block)
    wit = [('A', H_LO, [Fr(20), Fr(20)], [Fr(20), Fr(1), Fr(1)], 1, 2),
           ('B', H_HI, [Fr(1, 10), Fr(1)], [Fr(1, 10), Fr(1), Fr(1)], 1, 2)]
    for kind, h, wT, wU, u1, u2 in wit:
        qn = quant_at(wT, wU, u1, u2, h)
        mv = minval(qn, kind)
        qh = quant_at(wT, wU, u1, u2, Q(1, 2))
        good_half = minval(qh, 'A') > 0 and minval(qh, 'B') > 0
        print('witness (%s) at H = %s, w_T = %s, w_U = %s: smallest value %.4g -> %s; at H = 1/2 all positive: %s'
              % (kind, h, [str(x) for x in wT], [str(x) for x in wU], float(mv),
                 'FALSE there (as it must be)' if mv < 0 else 'not a witness', good_half), flush=True)
        ok &= mv < 0 and good_half
    # 3. B -> (1 + 1/1000) B
    found = 0
    for wT, wU in (([Fr(1, 100), Fr(1, 100)], [Fr(1, 100), Fr(1), Fr(1)]),
                   ([Fr(1), Fr(1)], [Fr(1, 100), Fr(1), Fr(1)]),
                   ([Fr(1), Fr(1), Fr(1)], [Fr(1, 30), Fr(1), Fr(1), Fr(1)]),
                   ([Fr(1, 20), Fr(1)], [Fr(1, 20), Fr(1, 20), Fr(1), Fr(1)])):
        _, piv1, M1, lab, w, T, U = direct(wT, wU)
        _, pivs, Ms, _, _, _, _ = direct(wT, wU, bscale=Fr(1001, 1000))
        Mc = TypeI(w, T).build()
        n = len(lab)
        differs = any(Ms[i][j] != Mc[i][j] for i in range(n) for j in range(n))
        negs = sum(1 for p in pivs if p < 0)
        print('B -> (1001/1000) B at w_T = %s, w_U = %s: negative pivots %d (true block: %d nonpositive, min pivot %.3e); '
              'scaled builder differs from closed form: %s' % ([str(x) for x in wT], [str(x) for x in wU], negs,
                                                                sum(1 for p in piv1 if p <= 0), float(min(piv1)), differs),
              flush=True)
        ok &= all(p > 0 for p in piv1) and differs
        found += negs > 0
    ok &= found > 0
    # 4. Theorem M1 with eps I subtracted: the certificate fails, and the perturbed matrix is indefinite somewhere
    import prove_eqT as PE
    PE.EPS = Q(1, 1000)
    res = PE.main()
    print('prove_eqT with eps = 1/1000: all certified %s (must be False)' % res, flush=True)
    ok &= not res
    Dd, names, K, G = PE.load('tl_inv_one_T22.pkl')
    indef = None
    for s1 in (Q(1, 10 ** 2), Q(1, 10 ** 4), Q(1, 10 ** 6), Q(10 ** 4)):
        for u2 in (2, 5):
            for sg in (2, 3):
                vals = {'u2': Q(u2), 'sig': Q(sg), 's1': s1}
                args = [vals[nm] for nm in names]
                ev = lambda f: f.n(*args) / f.d(*args)
                Dm = [[ev(G[(min(i, j), max(i, j))]) for j in range(3)] for i in range(3)]
                Om = [[ev(G[(i, 3 + j)]) for j in range(3)] for i in range(3)]
                for lab, M in (('D - O', [[Dm[i][j] - Om[i][j] for j in range(3)] for i in range(3)]),
                               ('D + (sigma-1) O', [[Dm[i][j] + (sg - 1) * Om[i][j] for j in range(3)] for i in range(3)])):
                    true_ok = all(p > 0 for p in ldl_pivots(flint.fmpq_mat(M)))
                    Me = [[M[i][j] - (Q(1, 1000) if i == j else 0) for j in range(3)] for i in range(3)]
                    if true_ok and any(p < 0 for p in ldl_pivots(flint.fmpq_mat(Me))) and indef is None:
                        indef = (lab, u2, sg, s1)
    print('perturbed M1 matrix indefinite at (matrix, u2, sigma, s1) = %s; unperturbed positive definite there' % (indef,),
          flush=True)
    ok &= indef is not None
    print('time %.1f s' % (time.time() - t0))
    return ok


if __name__ == '__main__':
    print('PASS' if main() else 'FAIL')
