"""Exact symbolic certificate for the type I block T = {t} of (R_inf)_1 at a given k (default 5), all positive weights.

Rescaled basis of ylib (entries homogeneous of degree 1 in the squared weights s), built with a lazy copy of
ylib.Builder (D^{-1} only on the column sectors that occur) over Q(s_0, ..., s_{k-1}) with s_1 = 1 (t = 0).
Certificate of Conjecture L at this k (Section 3 of NOTES.md, with Proposition T1.2):
 (i)  the all-ones vector on the P-states (1@r, 2@t) is annihilated identically;
 (ii) the (R, Q) block (states (2@t, 1@c) and (vac, 2@t 1@c)) is positive definite: all pivots of its symmetric
      elimination over the rational function field (ratios of consecutive leading principal minors) are ratios of
      polynomials with only positive coefficients;
 (iii) the Schur complement S_P onto the P-states is a Laplacian: S_P 1 = 0 identically, and for r != r' the
      numerator of S_P[r, r'] has only negative coefficients and its denominator only positive ones.
Then for every positive weight vector the block is >= 0 with kernel exactly p (x) |2@t>.
Usage: python sym_t1.py [k]"""
import sys, time, pickle
import flint
from flintfield import Field
from ylib import Ring
from lazyb import LazyBuilder
from j1lib import block_states, to_ylib_state

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
t0 = time.time()
names = ['s0'] + [f's{i}' for i in range(2, k)]
K = Field(names)
g = K.gens
s = [g[0], K(1)] + list(g[1:])
R = Ring(s, K(0), K(1))


B = LazyBuilder(R)
st = block_states(k, (0,))
ys = [to_ylib_state(z, k) for z in st]
n = len(st)
M = [[None] * n for _ in range(n)]
for i in range(n):
    for j in range(i, n):
        M[i][j] = M[j][i] = B.Y(ys[i], ys[j])
print(f'k = {k}: block of size {n} built over Q({",".join(names)}) with s1 = 1 [{time.time() - t0:.1f}s]', flush=True)
Ridx = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 2]
Pidx = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 1]
Qidx = [i for i, (a, y) in enumerate(st) if a is None]
ok = True
# (i)
oki = all(sum((M[i][j] for j in Pidx), K(0)) == 0 for i in range(n))
print(f'(i) all-ones on P annihilated identically: {oki}', flush=True)
ok &= oki


def allpos(p):
    return all(c > 0 for c in p.coeffs())


def allneg(p):
    return all(c < 0 for c in p.coeffs())


# (ii) + (iii): symmetric elimination of E = R + Q, in this order, on the whole matrix (E first, then P kept)
E = Ridx + Qidx
order = E + Pidx
A = [[M[i][j] for j in order] for i in order]
nE = len(E)
okii = True
pivinfo = []
for c in range(nE):
    p = A[c][c]
    okp = allpos(p.n) and allpos(p.d)
    if not okp:
        # pivot is a ratio; check sign-definiteness via factors
        cn, fn = p.n.factor(); cd, fd = p.d.factor()
        okp = (cn > 0 and all(allpos(f) for f, e in fn) and cd > 0 and all(allpos(f) for f, e in fd))
    pivinfo.append((p.n.total_degree(), len(p.n), okp))
    okii &= okp
    for i in range(c + 1, len(order)):
        if A[i][c] != 0:
            f = A[i][c] / p
            for j in range(c + 1, len(order)):
                if A[c][j] != 0:
                    A[i][j] = A[i][j] - f * A[c][j]
    print(f'  pivot {c + 1}/{nE}: numerator degree {pivinfo[-1][0]}, {pivinfo[-1][1]} terms, positive: {okp} '
          f'[{time.time() - t0:.1f}s]', flush=True)
print(f'(ii) all {nE} pivots of the (R, Q) block are ratios of positive-coefficient polynomials: {okii}')
ok &= okii
SP = [[A[nE + a][nE + b] for b in range(k)] for a in range(k)]
okiii = all(sum(SP[a], K(0)) == 0 for a in range(k))
offinfo = []
for a in range(k):
    for b in range(a + 1, k):
        v = SP[a][b]
        good = (allneg(v.n) and allpos(v.d))
        if not good:
            cn, fn = v.n.factor(); cd, fd = v.d.factor()
            good = (cn < 0 and all(allpos(f) for f, e in fn) and cd > 0 and all(allpos(f) for f, e in fd))
        offinfo.append((a, b, v.n.total_degree(), len(v.n), good))
        okiii &= good
print(f'(iii) S_P 1 = 0 identically and every off-diagonal entry negative (one-signed coefficients): {okiii}')
for x in offinfo:
    print('   ', x)
ok &= okiii
pickle.dump(dict(k=k, names=names, SP=[[(v.n.to_dict(), v.d.to_dict()) for v in row] for row in SP]),
            open(f'logs/sym_t1_k{k}.pkl', 'wb'))
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
