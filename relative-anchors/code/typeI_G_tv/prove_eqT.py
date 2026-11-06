"""Theorem M1 (one marked copy in U, equal weights elsewhere): sign certificates.

Setting: U = {h0} + U2, h0 of squared weight s1 > 0, the u2 >= 2 copies of U2 and all m >= 2 copies of T of squared
weight 1 (= s2).  Then sigma = m and G commutes with the permutations of T and of U2.  With D = G_tt and O = G_tt'
(t != t', the same for all pairs), restricted to the invariant part (3 x 3 in (rho, z^1, z^2)) and to the standard part
of U2 (1 x 1; u2 >= 2):
    G > 0  <=>  D - O > 0  and  D + (m - 1) O > 0     on each part
(the S_m-isotypic decomposition: x_t = xbar + xtilde_t with sum_t xtilde_t = 0).  Since sigma = m, D + (m - 1) O =
D + (sigma - 1) O.  We prove both for every real sigma >= 2, u2 >= 2 and s1 > 0 by Sylvester's criterion (invariant
part: the three leading principal minors in the order rho, z^1, z^2) and by coefficient signs: after u2 = 2 + x
(u2 = 3 + x for the standard part with u2 >= 3; u2 = 2 is a separate form) and sigma = 2 + y, every numerator has only
positive coefficients and a monomial free of x and y, and every irreducible factor of the denominator has
coefficients of one sign and a monomial free of x and y; the overall sign is +.
Usage: python prove_eqT.py [eps]   (eps: optional negative control, subtract eps * I from both matrices)"""
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from flintfield import Field, F as FE

Q = flint.fmpq
EPS = Q(sys.argv[1]) if len(sys.argv) > 1 else Q(0)


def load(fn):
    D = pickle.load(open(fn, 'rb'))
    names = D['names']
    K = Field(names)
    P = lambda d: K.ctx.from_dict({m: Q(p, q) for m, (p, q) in d.items()})
    G = {k: FE(K, P(n), P(d)) for k, (n, d) in D['G'].items()}
    return D, names, K, G


VN = ('x', 'y', 's1')
KX = Field(list(VN))
cx = KX.ctx
X, Y, S1 = cx.gens()
ONE = cx.from_dict({(0, 0, 0): 1})


def comp(f, names, u2shift):
    """rational function in names (subset of u2, sig, s1) -> (x, y, s1) with u2 = u2shift + x, sig = 2 + y."""
    m = {'u2': u2shift * ONE + X, 'sig': 2 * ONE + Y, 's1': S1}
    args = [m[nm] for nm in names]
    return FE(KX, f.n.compose(*args, ctx=cx), f.d.compose(*args, ctx=cx))


def cert(name, val, free=('x', 'y')):
    t0 = time.time()
    num, den = val.n, val.d
    idx = [VN.index(v) for v in free]
    cn = list(num.coeffs())
    neg = sum(1 for c in cn if c < 0)
    pos = len(cn) - neg
    sgn_num = 1 if neg == 0 else (-1 if pos == 0 else 0)
    freem = any(all(m[i] == 0 for i in idx) for m in num.monoms())
    c, fl = den.factor()
    sgn = 1 if c > 0 else -1
    den_ok = True
    for f, e in fl:
        fc = list(f.coeffs())
        if not any(all(m[i] == 0 for i in idx) for m in f.monoms()):
            den_ok = False
        if all(q > 0 for q in fc):
            continue
        if all(q < 0 for q in fc):
            if e % 2:
                sgn = -sgn
            continue
        den_ok = False
    ok = (sgn_num * sgn == 1) and freem and den_ok and sgn_num != 0
    print('%-28s numerator %6d terms (+%d, -%d), monomial free of %s: %s; denominator: %d factors, one-signed %s, '
          'sign %+d -> %s [%.1fs]' % (name, len(cn), pos, neg, free, freem, len(fl), den_ok, sgn,
                                      'CERTIFIED' if ok else 'not certified', time.time() - t0), flush=True)
    return ok


def minors3(M):
    m1 = M[0][0]
    m2 = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    m3 = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
          + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    return m1, m2, m3


def main():
    ok = True
    t0 = time.time()
    # ---- invariant part, u2 >= 2
    Dd, names, K, G = load('tl_inv_one_T22.pkl')
    sig = K.gens[names.index('sig')]
    Dm = [[G[(min(i, j), max(i, j))] for j in range(3)] for i in range(3)]
    Om = [[G[(i, 3 + j)] for j in range(3)] for i in range(3)]
    symO = all(Om[i][j] == Om[j][i] for i in range(3) for j in range(3))
    print('O = G_tt\' symmetric: %s' % symO, flush=True)
    ok &= symO
    for label, M in (('D - O', [[Dm[i][j] - Om[i][j] for j in range(3)] for i in range(3)]),
                     ('D + (sigma - 1) O', [[Dm[i][j] + (sig - 1) * Om[i][j] for j in range(3)] for i in range(3)])):
        if EPS != 0:
            M = [[M[i][j] - (EPS if i == j else 0) for j in range(3)] for i in range(3)]
        mins = minors3(M)
        for k, mk in enumerate(mins):
            ok &= cert('inv %s minor %d' % (label, k + 1), comp(mk, names, 2))
    # ---- standard part of U2: u2 >= 3 (form 'one') and u2 = 2 (form 'onetwo')
    for fn, shift in (('tl_std2_one_T22.pkl', 3), ('tl_std2_onetwo_T22.pkl', None)):
        Dd, names, K, G = load(fn)
        sig = K.gens[names.index('sig')]
        d, o = G[(0, 0)], G[(0, 1)]
        for label, v in (('d - o', d - o), ('d + (sigma - 1) o', d + (sig - 1) * o)):
            if EPS != 0:
                v = v - EPS
            ok &= cert('std2 (%s) %s' % ('u2 = 3 + x' if shift else 'u2 = 2', label),
                       comp(v, names, shift if shift else 0), free=('x', 'y') if shift else ('y',))
    print('all certified: %s [%.1fs]' % (ok, time.time() - t0))
    return ok


if __name__ == '__main__':
    ok = main()
    if EPS != 0:
        print('negative control (eps = %s): %s' % (EPS, 'certificate fails as expected' if not ok else 'UNEXPECTED: still certified'))
        print('PASS' if not ok else 'FAIL')
    else:
        print('PASS' if ok else 'FAIL')
