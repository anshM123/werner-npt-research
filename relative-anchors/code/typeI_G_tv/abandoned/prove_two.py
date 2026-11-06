"""Two-valued weight vectors: G > 0 by the exact decomposition under S_{u1} x S_{u2} x S_{m1} x S_{m2}.

Weights: squared weight s (= s1) or 1 (= s2, homogeneity).  U = U1 + U2 (u1 copies of weight s, u2 of weight 1),
T = T1 + T2 (m1 copies of weight s, m2 of weight 1), sigma = m1 s + m2, m = m1 + m2 >= 2.  With the closed forms of
tl_m.py at (a, b) = (s, s) [T11], (s, 1) [T12], (1, 1) [T22]:
    D1 = G_tt (s_t = s), D2 = G_tt (s_t = 1), O11 = G_tt' (s_t = s_t' = s), O12 = G_tt' (s_t = s, s_t' = 1), O22,
on each part of the U-decomposition (invariant: 3 x 3 in (rho, z^1, z^2); standard part of a level: 1 x 1).
G restricted to such a part is the m x m block matrix with these blocks; under S_{m1} x S_{m2} it splits into
  (S1) D1 - O11            (multiplicity m1 - 1; needed if m1 >= 2),
  (S2) D2 - O22            (multiplicity m2 - 1; needed if m2 >= 2),
  (I)  the T-invariant part: [[m1 (D1 + (m1-1) O11), m1 m2 O12], [m1 m2 O12^T, m2 (D2 + (m2-1) O22)]]
       (only the first block if m2 = 0, only the second if m1 = 0).
Each condition is proved positive definite on the whole parameter range by LDL^T pivots (Sylvester), every pivot a
rational function of orthant variables with one-signed numerator coefficients (and a monomial free of the shift
variables) and one-signed denominator factors:
  (S1) sigma = 2 s + r (r >= 0; covers m1 >= 2, any m2);   (S2) sigma = 2 + r;
  (I) m2 = 0: m1 = 2 + p, sigma = m1 s;   m1 = 0: m2 = 2 + p, sigma = m2;   m1, m2 >= 1: m1 = 1 + p, m2 = 1 + q.
Usage: python prove_two.py one            (u1 = 1, u2 >= 2)
       python prove_two.py gen            (u1, u2 >= 2)
       python prove_two.py one eps        (negative control: subtract eps*I; the certificate must fail)"""
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from flintfield import Field, F as FE

Q = flint.fmpq
CASE = sys.argv[1]
EPS = Q(sys.argv[2]) if len(sys.argv) > 2 else Q(0)


def load(fn):
    D = pickle.load(open(fn, 'rb'))
    names = D['names']
    ctx = flint.fmpq_mpoly_ctx.get(tuple(names), 'lex')
    P = lambda d: ctx.from_dict({m: Q(p, q) for m, (p, q) in d.items()})
    return names, D['keep'], {k: (P(n), P(d)) for k, (n, d) in D['G'].items()}


VN = ('x1', 'x2', 'p', 'q', 'r', 's')
KX = Field(list(VN))
cx = KX.ctx
GX = dict(zip(VN, cx.gens()))
ONE = cx.from_dict({(0,) * len(VN): 1})


def comp(nd, names, mp):
    n, d = nd
    args = [mp[nm] for nm in names]
    return FE(KX, n.compose(*args, ctx=cx), d.compose(*args, ctx=cx))


def cert(name, val, free):
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
        if all(z > 0 for z in fc):
            continue
        if all(z < 0 for z in fc):
            if e % 2:
                sgn = -sgn
            continue
        den_ok = False
    ok = (sgn_num * sgn == 1) and freem and den_ok and sgn_num != 0
    print('    %-44s num %7d terms (+%d -%d), free monomial %s, den %d factors one-signed %s -> %s [%.1fs]'
          % (name, len(cn), pos, neg, freem, len(fl), den_ok, 'CERTIFIED' if ok else 'NOT certified',
             time.time() - t0), flush=True)
    return ok


def pivots_cert(M, label, free):
    """LDL^T pivots of the symmetric matrix M (list of lists of FE); certify each pivot > 0."""
    n = len(M)
    A = [row[:] for row in M]
    ok = True
    for p in range(n):
        pv = A[p][p]
        if EPS != 0:
            pass
        ok &= cert('%s pivot %d' % (label, p + 1), pv, free)
        if not ok and EPS != 0:
            return False
        for i in range(p + 1, n):
            if A[i][p] == 0:
                continue
            f = A[i][p] / pv
            for j in range(p + 1, n):
                A[i][j] = A[i][j] - f * A[p][j]
    return ok


def blocks(prefix, half, mp):
    """the blocks D1, D2, O11, O12, O22 (half x half) composed with the substitution mp."""
    out = {}
    for mode in ('T11', 'T12', 'T22'):
        names, keep, G = load('%s_%s.pkl' % (prefix, mode))
        A = [[comp(G[(min(i, j), max(i, j))], names, mp) for j in range(half)] for i in range(half)]
        B = [[comp(G[(i, half + j)], names, mp) for j in range(half)] for i in range(half)]
        if mode == 'T11':
            out['D1'], out['O11'] = A, B
        elif mode == 'T12':
            out['O12'] = B
        else:
            out['D2'], out['O22'] = A, B
    return out


def sub(M, N, c=1):
    return [[M[i][j] - c * N[i][j] for j in range(len(M))] for i in range(len(M))]


def addm(M, N, c):
    return [[M[i][j] + c * N[i][j] for j in range(len(M))] for i in range(len(M))]


def scal(M, c):
    return [[c * M[i][j] for j in range(len(M))] for i in range(len(M))]


def eps_shift(M):
    if EPS == 0:
        return M
    return [[M[i][j] - (EPS if i == j else 0) for j in range(len(M))] for i in range(len(M))]


def run_part(prefix, half, umap, free_u):
    """all conditions (S1), (S2), (I) for one part (prefix = closed-form file prefix)."""
    ok = True
    s, r, p, q = GX['s'], GX['r'], GX['p'], GX['q']
    regimes = []
    # (S1): sigma = 2 s + r
    regimes.append(('S1 (m1 >= 2): D1 - O11, sigma = 2s + r', dict(umap, sig=2 * s + r, s1=s), 'S1', ('r',)))
    # (S2): sigma = 2 + r
    regimes.append(('S2 (m2 >= 2): D2 - O22, sigma = 2 + r', dict(umap, sig=2 * ONE + r, s1=s), 'S2', ('r',)))
    # (I) m2 = 0: m1 = 2 + p, sigma = m1 s
    regimes.append(('I (m2 = 0, m1 = 2 + p)', dict(umap, sig=(2 * ONE + p) * s, s1=s), 'I10', ('p',)))
    # (I) m1 = 0: m2 = 2 + p, sigma = m2
    regimes.append(('I (m1 = 0, m2 = 2 + p)', dict(umap, sig=2 * ONE + p, s1=s), 'I01', ('p',)))
    # (I) m1, m2 >= 1
    regimes.append(('I (m1 = 1 + p, m2 = 1 + q)', dict(umap, sig=(ONE + p) * s + ONE + q, s1=s), 'I11', ('p', 'q')))
    for label, mp, kind, fr in regimes:
        t0 = time.time()
        Bk = blocks(prefix, half, mp)
        print('  %s: blocks composed [%.1fs]' % (label, time.time() - t0), flush=True)
        free = tuple(free_u) + fr
        if kind == 'S1':
            M = sub(Bk['D1'], Bk['O11'])
        elif kind == 'S2':
            M = sub(Bk['D2'], Bk['O22'])
        elif kind == 'I10':
            M = addm(Bk["D1"], Bk["O11"], FE(KX, ONE + p, ONE))
        elif kind == 'I01':
            M = addm(Bk['D2'], Bk['O22'], FE(KX, ONE + p, ONE))
        else:
            m1 = FE(KX, ONE + p, ONE)
            m2 = FE(KX, ONE + q, ONE)
            P1 = scal(addm(Bk['D1'], Bk['O11'], FE(KX, p, ONE)), m1)
            P2 = scal(addm(Bk['D2'], Bk['O22'], FE(KX, q, ONE)), m2)
            C = scal(Bk['O12'], m1 * m2)
            n = half
            M = [[None] * (2 * n) for _ in range(2 * n)]
            for i in range(n):
                for j in range(n):
                    M[i][j] = P1[i][j]
                    M[n + i][n + j] = P2[i][j]
                    M[i][n + j] = C[i][j]
                    M[n + j][i] = C[i][j]
        M = eps_shift(M)
        res = pivots_cert(M, prefix + ' ' + kind, free)
        ok &= res
        if not res and EPS != 0:
            return False
    return ok


def main():
    t0 = time.time()
    ok = True
    x1, x2 = GX['x1'], GX['x2']
    if CASE == 'one':
        parts = [('tl_inv_one', 3, dict(u2=2 * ONE + x2), ('x2',)),
                 ('tl_std2_one', 1, dict(u2=3 * ONE + x2), ('x2',)),
                 ('tl_std2_onetwo', 1, dict(), ())]
    else:
        parts = [('tl_inv_gen', 3, dict(u1=2 * ONE + x1, u2=2 * ONE + x2), ('x1', 'x2')),
                 ('tl_std1_gen', 1, dict(u1=3 * ONE + x1, u2=ONE + x2), ('x1', 'x2')),
                 ('tl_std1_two', 1, dict(u2=ONE + x2), ('x2',))]
    for prefix, half, umap, free_u in parts:
        print('part %s' % prefix, flush=True)
        res = run_part(prefix, half, umap, free_u)
        ok &= res
        if not res and EPS != 0:
            ok = False
            break
    print('all certified: %s [%.1fs]' % (ok, time.time() - t0))
    return ok


if __name__ == '__main__':
    ok = main()
    if EPS != 0:
        print('negative control (eps = %s): %s' % (EPS, 'certificate fails as expected' if not ok else 'UNEXPECTED'))
        print('PASS' if not ok else 'FAIL')
    else:
        print('PASS' if ok else 'FAIL')
