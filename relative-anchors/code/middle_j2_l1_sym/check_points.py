"""The reduction of r2v.py against the explicit block built by the Lemma B builder (fastj1.FastJ1 = j1lib.J1.matrix),
exact, at two-level rational points.

python check_points.py U0 U1 Q [Q ...] [sample=S0,S1]   (Q a rational square, e.g. 1/9, 4, 49/25; with sample, the
bases are chosen at the sizes S0, S1, as cert.py does for a regime)

At each point: the rescaled block Y from FastJ1 (units 2^{k-2}, times mu(row) mu(column)); the identities
  (I (x) A2) Y_PQ = K   and   Y_QQ W3 = 2 I        (NOTES Lemma 2, Lemma V),
with A2, K, W3 built from their formulas on the concrete copies; S_ref = Y_PP - Y_PQ Y_QQ^{-1} Y_QP (exact inverse) and
S'_ref = (I (x) A2) S_ref (I (x) A2); then, for every family, the test vectors realised on the concrete copies (named
copies first in each level): the Gram matrices, the projections B^T S'_ref B on the whole basis K + W, and the kernel
identity (the rows of the K-vectors vanish) are compared with reduce_family at the same numbers.  The dimension
counts sum_lambda dim(lambda) * (number of vectors) are checked against the numbers of states."""
import sys, time, itertools
from fractions import Fraction as Fr
from math import comb
import flint
from fastj1 import FastJ1
from j1lib import mu
from r2lib import states_L1
import r2v

t0 = time.time()


def fq(x):
    return flint.fmpq(x.numerator, x.denominator)


def tofr(x):
    return Fr(int(x.p), int(x.q))


def sqrt_fr(q):
    from math import isqrt
    a, b = isqrt(q.numerator), isqrt(q.denominator)
    assert a * a == q.numerator and b * b == q.denominator, 'q must be a rational square'
    return Fr(a, b)


def run(u0, u1, q, sample=None):
    k = u0 + u1
    w = [Fr(1)] * u0 + [sqrt_fr(q)] * u1
    s = [x * x for x in w]
    lev = [0] * u0 + [1] * u1
    st = states_L1(k)
    M = FastJ1(w).matrix(st)
    n = len(st)
    mus = [mu(x, w) for x in st]
    Y = flint.fmpq_mat(n, n, [fq(M[i][j] * mus[i] * mus[j]) for i in range(n) for j in range(n)])
    Pst = [(a[1], tuple(t for t in range(k) if y[t] == 1)) for (a, y) in st if a is not None]
    Qst = [tuple(t for t in range(k) if y[t] == 1) for (a, y) in st if a is None]
    nP, nQ = len(Pst), len(Qst)
    pairs = list(itertools.combinations(range(k), 2))
    pidx = {X: i for i, X in enumerate(pairs)}
    YPP = flint.fmpq_mat(nP, nP, [Y[i, j] for i in range(nP) for j in range(nP)])
    YPQ = flint.fmpq_mat(nP, nQ, [Y[i, nP + j] for i in range(nP) for j in range(nQ)])
    YQQ = flint.fmpq_mat(nQ, nQ, [Y[nP + i, nP + j] for i in range(nQ) for j in range(nQ)])
    ctx = r2v.Ctx(Fr(1), q, Fr(u0), Fr(u1))
    levd = dict(enumerate(lev))
    A2 = flint.fmpq_mat(len(pairs), len(pairs), [fq(ctx.adj(Z1, Z2, levd)) for Z1 in pairs for Z2 in pairs])
    IA = flint.fmpq_mat(nP, nP)
    for i, (r, X) in enumerate(Pst):
        for j, (r2, X2) in enumerate(Pst):
            if r == r2:
                IA[i, j] = A2[pidx[X], pidx[X2]]
    Kf = flint.fmpq_mat(nP, nQ, [fq(ctx.kk(r, X, T, levd)) for (r, X) in Pst for T in Qst])
    W3 = flint.fmpq_mat(nQ, nQ, [fq(ctx.adj(T1, T2, levd)) for T1 in Qst for T2 in Qst])
    ok = True
    e1 = (IA * YPQ == Kf)
    e2 = (YQQ * W3 == 2 * flint.fmpq_mat(nQ, nQ, [1 if i == j else 0 for i in range(nQ) for j in range(nQ)]))
    # Lemma A item 1 with the Woodbury form: the PP oracle on the concrete copies
    e3 = all(YPP[i, j] == fq(ctx.pp(r, X, r2, X2, levd)) for i, (r, X) in enumerate(Pst)
             for j, (r2, X2) in enumerate(Pst))
    ok &= e1 and e2 and e3
    print(f'({u0}, {u1}) q = {q}: (I x A2) Y_PQ = K: {e1}; Y_QQ W3 = 2I: {e2}; Y_PP = PP oracle entrywise: {e3} '
          f'[{time.time() - t0:.1f}s]', flush=True)
    S = YPP - YPQ * YQQ.inv() * YPQ.transpose()
    Sp = IA * S * IA
    dimP = dimK = dimQ = 0
    for F in r2v.families():
        bases = r2v.choose(F, sample or (u0, u1))
        K, W, Qv = bases
        if not (K or W or Qv):
            continue
        d = r2v.hook_dim(u0, F.shape[0]) * r2v.hook_dim(u1, F.shape[1])
        dimP += d * (len(K) + len(W))
        dimK += d * len(K)
        dimQ += d * len(Qv)
        if not W:
            continue
        # realise the named copies: level-0 named on copies 0, 1, ..; level-1 named on u0, u0 + 1, ..
        cnt = [0, 0]
        place = {}
        for j, (l, _) in enumerate(F.named):
            place[j] = (0 if l == 0 else u0) + cnt[l]
            cnt[l] += 1
        inv = {c: j for j, c in place.items()}

        def pattern(tup):
            out, an = [], {}
            for c in tup:
                if c in inv:
                    out.append(inv[c])
                else:
                    if c not in an:
                        an[c] = len(an)
                    out.append(('A', lev[c], an[c]))
            return tuple(out)

        P = K + W
        BP = flint.fmpq_mat(nP, len(P))
        for i, (r, X) in enumerate(Pst):
            p = pattern((r,) + X)
            for c, (_, v) in enumerate(P):
                x = v.get(p, 0)
                if x:
                    BP[i, c] = x
        red = r2v.reduce_family(F, ctx, bases)
        GP = BP.transpose() * BP
        g1 = all(GP[i, j] == fq(red['GP'][i][j]) for i in range(len(P)) for j in range(len(P)))
        R = BP.transpose() * Sp * BP
        Sf = red['Sfull']
        g2 = all(R[i, j] == fq(Sf[i][j]) for i in range(len(P)) for j in range(len(P)))
        nK = len(K)
        g3 = all(R[i, j] == 0 for i in range(nK) for j in range(len(P)))
        g4 = True
        if Qv:
            BQ = flint.fmpq_mat(nQ, len(Qv))
            for i, T in enumerate(Qst):
                p = pattern(T)
                for c, (_, v) in enumerate(Qv):
                    x = v.get(p, 0)
                    if x:
                        BQ[i, c] = x
            GQ = BQ.transpose() * BQ
            g4 = all(GQ[i, j] == fq(red['GQ'][i][j]) for i in range(len(Qv)) for j in range(len(Qv)))
        good = g1 and g2 and g3 and g4
        ok &= good
        print(f'  {F.name}: |K| {len(K)}, |W| {len(W)}, |Q| {len(Qv)}: Gram P {g1}, Gram Q {g4}, '
              f'reduced S\' = projection {g2}, K-rows vanish {g3} [{time.time() - t0:.1f}s]', flush=True)
    dims = (dimP == k * comb(k, 2), dimK == comb(k, 2), dimQ == comb(k, 3))
    ok &= all(dims)
    print(f'  dimension counts (P, K, Q): {dims}', flush=True)
    return ok


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if '=' not in a]
    opts = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
    sample = tuple(int(v) for v in opts['sample'].split(',')) if 'sample' in opts else None
    u0, u1 = int(args[0]), int(args[1])
    if sample:
        print(f'bases chosen at the sample size {sample} (as in cert.py for the regime containing ({u0}, {u1}))')
    ok = True
    for qs in args[2:]:
        ok &= run(u0, u1, Fr(qs), sample)
    print(f'time {time.time() - t0:.1f}s')
    print('PASS' if ok else 'FAIL')
