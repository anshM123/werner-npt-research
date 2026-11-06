"""Negative controls, exact (NOTES Section 6).  B -> (1 + 1/1000) B: every term of Lemma B that contains an E^{-1} is
multiplied by b = 1001/1000.

(a) Explicit blocks, k = 5, 6: at (u0, u1) = (3, 2), (4, 2) the Schur complement S_b of the vacuum states of the block
    built by fastj1 (bscale b) is projected onto the bases of every family; these projections equal r2v.reduce_scaled
    at the same point (b = 1 and b = 1001/1000), all W-blocks are positive definite for b = 1, and for b = 1001/1000
    some W-block is not positive definite at an explicit weight ratio.
(b) Inside the certified regimes (k >= 11): with r2v.reduce_scaled at (u0, u1) = (9, 2), (8, 3), (7, 4), (6, 5),
    (12, 9), (21, 20), every W-block is positive definite for b = 1 and some W-block is not for b = 1001/1000, at an
    explicit q.  So the inequality certified in these regimes does not survive a 0.1% increase of the B-part.
PASS iff all of this holds."""
import time
from fractions import Fraction as Fr
import flint
from fastj1 import FastJ1
from j1lib import mu
from r2lib import states_L1
import r2v

t0 = time.time()
B1 = Fr(1001, 1000)


def fq(x):
    return flint.fmpq(x.numerator, x.denominator)


def tofr(x):
    return Fr(int(x.p), int(x.q))


def wblocks(u0, u1, q, b):
    """dict family -> W-block of the reduced S_b (r2v.reduce_scaled)"""
    out = {}
    ctx = r2v.Ctx(Fr(1), q, Fr(u0), Fr(u1))
    for F in r2v.families():
        bases = r2v.choose(F, (u0, u1))
        if not bases[1]:
            continue
        S = r2v.reduce_scaled(F, ctx, bases, b)
        nK = len(bases[0])
        out[F.name] = (S, [row[nK:] for row in S[nK:]])
    return out


def explicit(u0, u1, w1, b):
    k = u0 + u1
    w = [Fr(1)] * u0 + [w1] * u1
    lev = [0] * u0 + [1] * u1
    st = states_L1(k)
    M = FastJ1(w, bscale=b).matrix(st)
    mus = [mu(x, w) for x in st]
    n = len(st)
    Y = flint.fmpq_mat(n, n, [fq(M[i][j] * mus[i] * mus[j]) for i in range(n) for j in range(n)])
    Pst = [(a[1], tuple(t for t in range(k) if y[t] == 1)) for (a, y) in st if a is not None]
    nP, nQ = len(Pst), n - len(Pst)
    YPP = flint.fmpq_mat(nP, nP, [Y[i, j] for i in range(nP) for j in range(nP)])
    YPQ = flint.fmpq_mat(nP, nQ, [Y[i, nP + j] for i in range(nP) for j in range(nQ)])
    YQQ = flint.fmpq_mat(nQ, nQ, [Y[nP + i, nP + j] for i in range(nQ) for j in range(nQ)])
    S = YPP - YPQ * YQQ.inv() * YPQ.transpose()
    out = {}
    for F in r2v.families():
        bases = r2v.choose(F, (u0, u1))
        K, W, Qv = bases
        if not W:
            continue
        cnt, place = [0, 0], {}
        for j, (l, _) in enumerate(F.named):
            place[j] = (0 if l == 0 else u0) + cnt[l]
            cnt[l] += 1
        inv = {c: j for j, c in place.items()}

        def pattern(tup):
            o, an = [], {}
            for c in tup:
                if c in inv:
                    o.append(inv[c])
                else:
                    an.setdefault(c, len(an))
                    o.append(('A', lev[c], an[c]))
            return tuple(o)

        P = K + W
        BP = flint.fmpq_mat(nP, len(P))
        for i, (r, X) in enumerate(Pst):
            p = pattern((r,) + X)
            for c, (_, v) in enumerate(P):
                x = v.get(p, 0)
                if x:
                    BP[i, c] = x
        R = BP.transpose() * S * BP
        out[F.name] = [[tofr(R[i, j]) for j in range(len(P))] for i in range(len(P))]
    return out


ok = True
# (a) explicit blocks at small k
for (u0, u1) in ((3, 2), (4, 2)):
    found = None
    for w1 in (Fr(1, 30), Fr(1, 10), Fr(10), Fr(30)):
        q = w1 * w1
        good = True
        verdict = {}
        for b in (Fr(1), B1):
            E = explicit(u0, u1, w1, b)
            Rr = wblocks(u0, u1, q, b)
            same = all(E[nm] == Rr[nm][0] for nm in Rr) and set(E) == set(Rr)
            good &= same
            verdict[b] = [nm for nm in Rr if not r2v.pd_exact(Rr[nm][1])]
        ok &= good
        if good and not verdict[Fr(1)] and verdict[B1]:
            found = (w1, verdict[B1])
            break
    ok &= found is not None
    print(f'(a) ({u0}, {u1}): explicit projections = reduced S_b for b = 1, 1001/1000: {good}; '
          f'failing point w = {found[0] if found else None} (q = {found[0] ** 2 if found else None}), '
          f'families with a W-block not PD at b = 1001/1000: {found[1] if found else None} [{time.time() - t0:.1f}s]',
          flush=True)
# (b) inside the regimes
for (u0, u1) in ((9, 2), (8, 3), (7, 4), (6, 5), (12, 9), (21, 20)):
    found = None
    for q in (Fr(1, 10 ** 4), Fr(10 ** 4), Fr(1, 100), Fr(100), Fr(1, 10 ** 6), Fr(10 ** 6)):
        t = wblocks(u0, u1, q, Fr(1))
        if not all(r2v.pd_exact(v[1]) for v in t.values()):
            ok = False
            print(f'(b) ({u0}, {u1}) q = {q}: a W-block is NOT PD for the true B (contradiction)')
            break
        s = wblocks(u0, u1, q, B1)
        bad = [nm for nm, v in s.items() if not r2v.pd_exact(v[1])]
        if bad:
            found = (q, bad)
            break
    ok &= found is not None
    print(f'(b) ({u0}, {u1}): all W-blocks PD at b = 1; at b = 1001/1000, q = {found[0] if found else None}: '
          f'not PD for {found[1] if found else None} [{time.time() - t0:.1f}s]', flush=True)
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
