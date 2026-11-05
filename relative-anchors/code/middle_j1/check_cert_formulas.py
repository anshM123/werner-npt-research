"""The certificates certify the formulas they contain; this script ties those formulas to the true blocks.
It executes the part of cert_t2.py and cert_t3j2.py that builds the vertex functions (their source up to the first
line 'ok = True', so the very same symbolic objects are used), evaluates them exactly at rational points, and compares
with the blocks computed from Lemma B by j1lib (independently of the closed forms):
(a) Theorem D1, k = 3 (one other copy, so lam = 1 and all eight vertices coincide): every vertex value equals
    f = det M / M[2][2] of the true 3 x 3 block (det M = (2/G3) f, Lemma D1.2), exactly, at 30 weight vectors.
(b) Theorem D1, k = 4..8: at the configuration's lam, every vertex value is positive and the smallest is <= the true f.
(c) Theorem D2, k = 4 (lam = 1): for all 64 vertices, 2 - G_T a^T C^{-1} a (exact linear algebra) equals the true
    value of the 4 x 4 block, and the certificate's expression P/det C (adjugate formula) agrees, at 20 weight vectors.
(d) Theorem D2, k = 5..8: at the configuration's lam, the smallest vertex value is positive and <= the true value.
Usage: python check_cert_formulas.py"""
import random, itertools, io, contextlib
from fractions import Fraction as Fr
import flint
from j1lib import J1, block_states, onehot

random.seed(17)


def prefix_namespace(path):
    src = open(path, encoding='utf-8').read()
    cut = src.index('\nok = True')
    ns = {'__name__': 'cert_prefix'}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(src[:cut], path, 'exec'), ns)
    return ns


def ev(v, pt):
    """exact value of a flintfield element at a point (list of fmpq)."""
    num = v.n(*pt)
    den = v.d(*pt)
    q = num / den
    return Fr(int(q.p), int(q.q))


def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


def solve(A, b):
    n = len(A)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for c in range(n):
        p = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[p] = M[p], M[c]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c] / M[c][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


def rand_w(k, mode):
    if mode == 0:
        return [Fr(random.randint(1, 20), random.randint(1, 20)) for _ in range(k)]
    if mode == 1:
        return [Fr(2) ** random.randint(-4, 4) for _ in range(k)]
    w = [Fr(random.randint(1, 9), random.randint(1, 9)) for _ in range(k)]
    w[random.randrange(k)] = Fr(1, 40)
    return w


def fq(x):
    return flint.fmpq(x.numerator, x.denominator)


ok = True
# ---------------- Theorem D1
N2 = prefix_namespace('cert_t2.py')
Frng, kap, G0, ap, bp, K1, K2, G3 = (N2[n] for n in ('Frng', 'kap', 'G0', 'ap', 'bp', 'K1', 'K2', 'G3'))


def t2_vertices(pt):
    vals = []
    for i0, i1, i2 in itertools.product((0, 1), repeat=3):
        c = ev(kap, pt) / (2 * (ev(G0, pt) + ev(Frng[0][i0], pt) / 2))
        A = ev(ap, pt) / (ev(K1, pt) + ev(Frng[1][i1], pt) / 2)
        B = ev(bp, pt) / (ev(K2, pt) + ev(Frng[2][i2], pt) / 2)
        g3 = ev(G3, pt)
        vals.append(4 - c * c - g3 * (A * A + B * B - c * A * B))
    return vals


def t2_true(w):
    k = len(w)
    M = J1(w).matrix(block_states(k, (0, 1)))
    return det3(M) / M[2][2]


def t2_point(w):
    s = [v * v for v in w]; S = sum(s); x = [v / S for v in s]
    al, be = x[0], x[1]; e = 1 - al - be
    F3 = sum(x[u] / (2 - 2 * x[u] - al - be) for u in range(2, len(w)))
    Z = F3 / e; Z0 = 1 / (al + be + 2 * e); Z1 = 1 / (al + be)
    lam = (Z - Z0) / (Z1 - Z0)
    return [fq(al), fq(be), fq(e), fq(lam)], lam


gooda = True
for i in range(30):
    w = rand_w(3, i % 3)
    pt, lam = t2_point(w)
    ft = t2_true(w)
    vals = t2_vertices(pt)
    gooda &= lam == 1 and all(v == ft for v in vals)
print(f'(a) D1, k = 3, lam = 1: all 8 vertex values equal the true f, 30 weight vectors: {gooda}', flush=True)
goodb = True
for i in range(30):
    k = random.randint(4, 8)
    w = rand_w(k, i % 3)
    pt, lam = t2_point(w)
    ft = t2_true(w)
    vals = t2_vertices(pt)
    goodb &= 0 < lam < 1 and min(vals) > 0 and min(vals) <= ft
print(f'(b) D1, k = 4..8: vertex values positive and min <= true f, 30 weight vectors: {goodb}', flush=True)
ok &= gooda and goodb
# ---------------- Theorem D2
N3 = prefix_namespace('cert_t3j2.py')
crng, arng, GT = N3['crng'], N3['arng'], N3['GT']


def states_j2T3(k, T):
    st = []
    for t in T:
        rest = [u for u in T if u != t]
        st.append(((2, t), onehot(k, {rest[0]: 2, rest[1]: 2})))
    st.append((None, onehot(k, {u: 2 for u in T})))
    return st


def t3_true(w):
    k = len(w)
    J = J1(w)
    st = states_j2T3(k, (0, 1, 2))
    M = [[J.entry(a, b) for b in st] for a in st]
    C = [row[:3] for row in M[:3]]
    a = [-M[i][3] for i in range(3)]
    gT = 2 / M[3][3]
    y = solve(C, a)
    return 2 - gT * sum(a[i] * y[i] for i in range(3))


def t3_point(w):
    s = [v * v for v in w]; S = sum(s); x = [v / S for v in s]
    a1, a2, a3 = x[0], x[1], x[2]; xT = a1 + a2 + a3; e = 1 - xT
    FT = sum(x[u] / (2 - 2 * x[u] - xT) for u in range(3, len(w)))
    Z = FT / e; Z0 = 1 / (xT + 2 * e); Z1 = 1 / xT
    lam = (Z - Z0) / (Z1 - Z0)
    return [fq(a1), fq(a2), fq(a3), fq(e), fq(lam)], lam


def t3_vertices(pt):
    cv_ = [[ev(v, pt) for v in pair] for pair in crng]
    av_ = [[ev(v, pt) for v in pair] for pair in arng]
    g = ev(GT, pt)
    out = []
    for cv in itertools.product((0, 1), repeat=3):
        c12, c13, c23 = cv_[0][cv[0]], cv_[1][cv[1]], cv_[2][cv[2]]
        C = [[Fr(2), c12, c13], [c12, Fr(2), c23], [c13, c23, Fr(2)]]
        det = 8 - 2 * (c12 * c12 + c13 * c13 + c23 * c23) + 2 * c12 * c13 * c23
        adj = [[4 - c23 * c23, c13 * c23 - 2 * c12, c12 * c23 - 2 * c13],
               [c13 * c23 - 2 * c12, 4 - c13 * c13, c12 * c13 - 2 * c23],
               [c12 * c23 - 2 * c13, c12 * c13 - 2 * c23, 4 - c12 * c12]]
        for av in itertools.product((0, 1), repeat=3):
            a = [av_[i][av[i]] for i in range(3)]
            y = solve(C, a)
            f = 2 - g * sum(a[i] * y[i] for i in range(3))
            P = 2 * det - g * sum(a[i] * adj[i][j] * a[j] for i in range(3) for j in range(3))
            out.append((f, P / det, det))
    return out


goodc = True
for i in range(20):
    w = rand_w(4, i % 3)
    pt, lam = t3_point(w)
    ft = t3_true(w)
    vals = t3_vertices(pt)
    goodc &= lam == 1 and all(f == ft and q == f and d > 0 for f, q, d in vals)
print(f'(c) D2, k = 4, lam = 1: all 64 vertex values equal the true value (and P/det C agrees), 20 weight vectors: '
      f'{goodc}', flush=True)
goodd = True
for i in range(20):
    k = random.randint(5, 8)
    w = rand_w(k, i % 3)
    pt, lam = t3_point(w)
    ft = t3_true(w)
    vals = t3_vertices(pt)
    fs = [f for f, q, d in vals]
    goodd &= 0 < lam < 1 and min(fs) > 0 and min(fs) <= ft and all(q == f for f, q, d in vals)
print(f'(d) D2, k = 5..8: vertex values positive and min <= true value, 20 weight vectors: {goodd}', flush=True)
ok &= goodc and goodd
print('PASS' if ok else 'FAIL')
