"""Polynomial data for inequalities (A) and (B), normalised S_T = 1, derived EXACTLY (SymPy) from gstruct.pkl.

Closed forms (scaled coordinates, Gt = D G D, D_t = diag(1, w_t omega)), S = 1:
  off-diagonal:  Gt_off = [[-s s2 p00, -s s2 sig p01], [-s s2 sig p10, s s2 sig p11]] / (2 Den)
  diagonal:      Gt_diag = [[pd00, s sig pd01], [s sig pd01, s sig pd11]] / (2 Dd)
with
  Den = (2 - s + sig)(2 - s2 + sig)(3 + 2 sig + sig^2) det2(s) det2(s2) E3(s, s2, sig),
  Dd  = (2 - s + sig)^2 (3 + 2 sig + sig^2) det2(s)^2,
  det2(x) = 2 + x + 5 sig - x^2 + 3 x sig + 2 sig^2   (= 2 det_t),
  E3 = the cubic factor printed by gstruct (proportional to det E_{tt'}).
All p's are polynomials with integer coefficients (checked below by exact cancellation).

Rescaled 2x2 blocks used in the proof (c^2 = 1 + sig):
  Ghat_tt' = diag(1, c) G_tt' diag(1, c),   N = Ghat_tt' / sqrt(s s2),
  ||N||_F^2 = [ s s2 p00^2 + (1 + sig) sig (s p01^2 + s2 p10^2) + (1 + sig)^2 p11^2 ] / (4 Den^2)
  Ghat_tt = [[pd00, c sqrt(s sig) pd01], [c sqrt(s sig) pd01, (1 + sig) pd11]] / (2 Dd)
Inequalities (H rational):
  FA  = 4 H^2 Den^2 - [ ... ] >= 0                         <=>  ||N||_F <= H
  FB1 = (1 + sig) pd11 - 2 Dd h > 0,  h = H (1 - s)
  FB2 = (pd00 - 2 Dd h)((1 + sig) pd11 - 2 Dd h) - (1 + sig) s sig pd01^2 > 0
  (FB1 > 0 and FB2 > 0  <=>  Ghat_tt - h I > 0  <=>  lambda_min(Ghat_tt) > H (1 - s)).
"""
import pickle, sympy as sp

s, s2, S, sg = sp.symbols('s s2 S sigma', positive=True)


def build(H=sp.Rational(2, 5), verbose=True):
    Go, Gd = pickle.load(open('gstruct.pkl', 'rb'))
    sub = {S: 1}
    Go = Go.subs(sub)
    Gd = Gd.subs(sub)
    det2 = lambda x: 2 + x + 5 * sg - x ** 2 + 3 * x * sg + 2 * sg ** 2
    E3 = E3_expr
    Den = sp.expand((2 - s + sg) * (2 - s2 + sg) * (3 + 2 * sg + sg ** 2) * det2(s) * det2(s2) * E3)
    Dd = sp.expand((2 - s + sg) ** 2 * (3 + 2 * sg + sg ** 2) * det2(s) ** 2)
    pre_off = [[-s * s2, -s * s2 * sg], [-s * s2 * sg, s * s2 * sg]]
    pre_dg = [[1, s * sg], [s * sg, s * sg]]
    P = [[None, None], [None, None]]
    Pd = [[None, None], [None, None]]
    for i in range(2):
        for j in range(2):
            q = sp.cancel(Go[i, j] * 2 * Den / pre_off[i][j])
            P[i][j] = sp.Poly(q, s, s2, sg)          # raises if not a polynomial
            assert all(c.is_integer for c in P[i][j].coeffs())
            q = sp.cancel(Gd[i, j] * 2 * Dd / pre_dg[i][j])
            Pd[i][j] = sp.Poly(q, s, sg)
            assert all(c.is_integer for c in Pd[i][j].coeffs())
    assert sp.expand(Pd[0][1].as_expr() - Pd[1][0].as_expr()) == 0
    p00, p01, p10, p11 = [P[i][j].as_expr() for i, j in ((0, 0), (0, 1), (1, 0), (1, 1))]
    pd00, pd01, pd11 = Pd[0][0].as_expr(), Pd[0][1].as_expr(), Pd[1][1].as_expr()
    # exact re-check of the representation
    for i in range(2):
        for j in range(2):
            assert sp.cancel(Go[i, j] - pre_off[i][j] * P[i][j].as_expr() / (2 * Den)) == 0
            assert sp.cancel(Gd[i, j] - pre_dg[i][j] * Pd[i][j].as_expr() / (2 * Dd)) == 0
    br = s * s2 * p00 ** 2 + (1 + sg) * sg * (s * p01 ** 2 + s2 * p10 ** 2) + (1 + sg) ** 2 * p11 ** 2
    FA = sp.expand(4 * H ** 2 * Den ** 2 - br)
    h = H * (1 - s)
    FB1 = sp.expand((1 + sg) * pd11 - 2 * Dd * h)
    FB2 = sp.expand((pd00 - 2 * Dd * h) * ((1 + sg) * pd11 - 2 * Dd * h) - (1 + sg) * s * sg * pd01 ** 2)
    out = dict(H=H, Den=Den, Dd=Dd, p=(p00, p01, p10, p11), pd=(pd00, pd01, pd11), FA=FA, FB1=FB1, FB2=FB2, E3=E3)
    if verbose:
        for k in ('FA', 'FB1', 'FB2'):
            pk = sp.Poly(out[k], *( (s, s2, sg) if k == 'FA' else (s, sg) ))
            print(k, 'degrees', pk.degree_list(), 'terms', len(pk.terms()), flush=True)
    return out


def check_homogeneity():
    """Numerators and denominators of the closed forms (general S) are homogeneous in (s, s2, S, sigma); the degree
    of Gt_off[i][j] (= deg num - deg den) is 0, 1, 1, 2 for (0,0), (0,1), (1,0), (1,1), and 0, 1, 2 for the diagonal
    (0,0), (0,1), (1,1).  Hence N = Ghat_tt' S_T / sqrt(s_t s_t') and Ghat_tt are invariant under
    (s_t, s_t', S_T, sigma) -> lambda (s_t, s_t', S_T, sigma): one may normalise S_T = 1."""
    Go, Gd = pickle.load(open('gstruct.pkl', 'rb'))
    want_off = [[0, 1], [1, 2]]
    want_dg = [[0, 1], [1, 2]]
    ok = True
    for i in range(2):
        for j in range(2):
            for M, gens, want in ((Go, (s, s2, S, sg), want_off), (Gd, (s, S, sg), want_dg)):
                num, den = sp.fraction(sp.together(M[i, j]))
                pn, pd = sp.Poly(num, *gens), sp.Poly(den, *gens)
                ok &= pn.is_homogeneous and pd.is_homogeneous and (pn.total_degree() - pd.total_degree() == want[i][j])
    return ok


def check_E3():
    """E3 = 4 det E_{tt'} at S = 1 (E_{tt'} = diag(s_x + K - (s + s2)/2) + w_F w_F^T / 2 on F = {0, t, t'})."""
    om, w, w2 = sp.sqrt(sg), sp.sqrt(s), sp.sqrt(s2)
    K = 1 + sg / 2
    base = K - s / 2 - s2 / 2
    wF = [om, w, w2]; sF = [sg, s, s2]
    E = sp.Matrix(3, 3, lambda i, j: (sF[i] + base if i == j else 0) + wF[i] * wF[j] / 2)
    return sp.expand(4 * E.det() - E3_expr) == 0


E3_expr = (4 + 12 * sg - 3 * s ** 2 + 2 * s * s2 + 2 * s * sg - 3 * s2 ** 2 + 2 * s2 * sg + 9 * sg ** 2 + s ** 3
           - s ** 2 * s2 - 4 * s ** 2 * sg - s * s2 ** 2 + 6 * s * s2 * sg + s * sg ** 2 + s2 ** 3 - 4 * s2 ** 2 * sg
           + s2 * sg ** 2 + 2 * sg ** 3)


if __name__ == '__main__':
    import time
    t0 = time.time()
    out = build()
    pickle.dump(out, open('polys.pkl', 'wb'))
    hom = check_homogeneity()
    print('closed forms homogeneous with the stated degrees (normalisation S_T = 1 legitimate):', hom, flush=True)
    e3 = check_E3()
    print('E3 = 4 det E_tt\' (S = 1):', e3, flush=True)
    with open('polys.txt', 'w') as f:
        f.write('# Polynomials of polys.py (S_T = 1, H = %s); variables s, s2, sigma\n' % out['H'])
        for k in ('Den', 'Dd', 'E3', 'FB1', 'FB2', 'FA'):
            f.write(f'\n{k} =\n{out[k]}\n')
        for nm, p in zip(('p00', 'p01', 'p10', 'p11'), out['p']):
            f.write(f'\n{nm} =\n{sp.expand(p)}\n')
        for nm, p in zip(('pd00', 'pd01', 'pd11'), out['pd']):
            f.write(f'\n{nm} =\n{sp.expand(p)}\n')
    print(f'built and checked ({time.time() - t0:.1f} s)')
    print('PASS' if (hom and e3) else 'FAIL')
