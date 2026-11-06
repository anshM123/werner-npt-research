"""Exact checks of the irrep split (splitlib.py) against explicit objects at small k.
For k = 2..KMAX, every j = -1..k, every i with V_i nonzero and a split block, and sigma in {4, 1/9} (squared weight of
the marked copy), with B scaled by 1 and by 1001/1000:
 (a) the new basis vectors U, W, S, built as vectors over the states from the explicit twisted vectors of explicit.py,
     give v_a^T Y~ v_b (Y~ the rescaled Lemma B block of j1lib.py) equal to the entry of T^T Q T computed from the
     numeric reduced form of mlib.py; in particular the off-diagonal block between the two parts vanishes;
 (b) (B unscaled) the class sum C2 of the transpositions of S_u acts on every S vector as the content sum of
     (u - i - 1, i, 1) and on every vector of the first part as that of (u - i, i), exactly;
 (c) (B unscaled) the known kernel vectors, written in the new basis by splitlib.to_new, are annihilated by T^T Q T.
Usage: python check_split.py KMAX"""
import sys, time
from fractions import Fraction as Fr
import flint
import mlib, explicit, splitlib
from prove_fam import fam_kernel

KMAX = int(sys.argv[1])
fq = flint.fmpq
T0 = time.time()
ok = True


def F(x):
    return Fr(int(x.p), int(x.q))


def content(lam):
    return sum(c - r for r, row in enumerate(lam) for c in range(row))


def act(state, a, b):
    row, y = state
    y = list(y)
    y[a], y[b] = y[b], y[a]
    if row is not None:
        r = row[1]
        r = b if r == a else a if r == b else r
        row = (row[0], r)
    return (row, tuple(y))


class NumK:
    """stand-in for the kernel builder of prove_fam (which expects a family calculator)"""
    def __init__(self, k, j, i, sig):
        self.ft, self.fi = k - 1 - j - i, i
        self.u = fq(k - 1)
        self.sig = fq(sig.numerator, sig.denominator)
        self.one = fq(1)
        self.K = lambda c: fq(c)


nforms = nent = 0
for k in range(2, KMAX + 1):
    u = k - 1
    for sig in (Fr(4), Fr(1, 9)):
        w0 = Fr(2) if sig == 4 else Fr(1, 3)
        w = [w0] + [Fr(1)] * (k - 1)
        for bs in (Fr(1), Fr(1001, 1000)):
            for j in range(-1, k + 1):
                st, M = explicit.block(k, j, w, bscale=bs)
                idx = {s: n for n, s in enumerate(st)}
                for i in range(0, u // 2 + 1):
                    names, norms, Q = mlib.reduced(mlib.Calc(k, j, i, sigma=sig, bscale=bs))
                    if not names:
                        continue
                    one = fq(1)
                    vecs, split = splitlib.basis(names, fq(j - i), fq(u - 2 * i), fq(i), one)
                    if not split:
                        continue
                    two, three, Mn = splitlib.congruence(names, Q, vecs, fq(0))
                    st2, ex = explicit.twisted(k, j, i)
                    sv = {}
                    for nn, comb in vecs.items():
                        v = {}
                        for on, cf in comb.items():
                            for x, c in ex[on].items():
                                v[x] = v.get(x, 0) + F(cf) * c
                        sv[nn] = {x: c for x, c in v.items() if c != 0}
                    new = two + three
                    for a, na in enumerate(new):
                        for nb in new[a:]:
                            val = explicit.form(M, sv[na], sv[nb])
                            nent += 1
                            if val != F(Mn[na, nb]):
                                ok = False
                                print('SPLIT FORM DIFFERS', k, j, i, sig, bs, na, nb, val, F(Mn[na, nb]))
                            if (na in three) != (nb in three) and val != 0:
                                ok = False
                                print('OFF-DIAGONAL BLOCK NONZERO', k, j, i, sig, bs, na, nb)
                    if bs == 1:
                        lam2, lam3 = (u - i, i), (u - i - 1, i, 1)
                        for nn in new:
                            c = content(lam3 if nn in three else lam2)
                            img = {}
                            for x, cx in sv[nn].items():
                                for a in range(1, k):
                                    for b in range(a + 1, k):
                                        y = idx[act(st[x], a, b)]
                                        img[y] = img.get(y, 0) + cx
                            diff = {x: img.get(x, 0) - c * sv[nn].get(x, 0) for x in set(img) | set(sv[nn])}
                            if any(v != 0 for v in diff.values()):
                                ok = False
                                print('CLASS SUM EIGENVALUE WRONG', k, j, i, nn)
                        for lab, v, dl in fam_kernel(names, NumK(k, j, i, sig)):
                            wv = splitlib.to_new(v, split)
                            if wv is None or any(nm not in two for nm in wv):
                                ok = False
                                print('KERNEL VECTOR NOT IN THE FIRST PART', k, j, i, lab)
                                continue
                            for a in two:
                                s = fq(0)
                                for nm, cf in wv.items():
                                    s += Mn[a, nm] * cf
                                if s != 0:
                                    ok = False
                                    print('KERNEL VECTOR NOT ANNIHILATED', k, j, i, lab)
                    nforms += 1
    print(f'k={k}: done ({time.time() - T0:.0f}s)', flush=True)
print(f'{nforms} split forms, {nent} entries against the Lemma B block; class sums and kernel vectors checked')
print('PASS' if ok else 'FAIL', f'{time.time() - T0:.0f}s')
