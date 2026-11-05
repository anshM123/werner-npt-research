"""Checks for Section 3 of NOTES.md: the type I blocks T = {t} of (R_inf)_1 (size 3k - 1).
States R_c = (2@t, 1@c) (c in [k]), P_r = (1@r, 2@t) (r in [k]), Q_c = (vac, 2@t 1@c) (c != t).
Rescaled coordinates (ylib): every state divided by mu (R_c, Q_c: w_c; P_r: w_r).  Units 2^{k-1}, |s| arbitrary.
(a) Structure (Lemma T1.1), exactly at rational weights: the R-block is diagonal with entries 2 s_c (c != t) and s_t
    (c = t); the P-block is Pi = diag(N s) - (N s)(N s)^T/(2|s| - s_t) with N_r = 2 (r != t), 1 (r = t); the
    all-ones vector on P is in the kernel of the block.
(b) Conjecture L, exact checks: the (R, Q) block is positive definite, and the Schur complement S_P onto the
    P-states is a Laplacian: S_P 1 = 0 and S_P[r, r'] < 0 for r != r'.  Hence the block is >= 0 with kernel exactly
    p (x) |2@t> at that weight vector.  Rational weights, k = 3..10, six families (random, spread over powers of 3,
    one tiny weight, two tiny weights, one dominant weight, all weights equal but one), t at random.
(c) Negative control: with the B-part scaled by 1 + 1/1000 the gauge-fixed block (one P-state deleted) is not positive
    definite at an explicit rational point where the true one is.
Usage: python check_t1.py [kmax] [points_per_k]"""
import sys, random
from fractions import Fraction as Fr
from j1lib import J1, block_states, mu

kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 10
npts = int(sys.argv[2]) if len(sys.argv) > 2 else 12
random.seed(7)


def rescaled(J, st):
    M = J.matrix(st)
    m = [mu(z, J.w) for z in st]
    return [[m[i] * M[i][j] * m[j] for j in range(len(st))] for i in range(len(st))]


def ldl_pd(A):
    A = [row[:] for row in A]
    n = len(A)
    for c in range(n):
        if A[c][c] <= 0:
            return False
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / A[c][c]
                for j in range(c + 1, n):
                    A[i][j] -= f * A[c][j]
    return True


def schur(M, keep, elim):
    """exact Schur complement M[keep,keep] - M[keep,elim] M[elim,elim]^{-1} M[elim,keep] (Gauss-Jordan)."""
    n = len(elim)
    C = [[M[i][j] for j in elim] for i in elim]
    Bt = [[M[i][j] for j in keep] for i in elim]
    A = [C[i] + Bt[i] for i in range(n)]
    for c in range(n):
        p = A[c][c]
        A[c] = [v / p for v in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[c])]
    X = [row[n:] for row in A]   # C^{-1} B^T
    S = [[M[keep[a]][keep[b]] - sum(M[keep[a]][elim[l]] * X[l][b] for l in range(n)) for b in range(len(keep))]
         for a in range(len(keep))]
    return S


def weights(k, mode):
    if mode == 0:
        return [Fr(random.randint(1, 40), random.randint(1, 40)) for _ in range(k)]
    if mode == 1:
        return [Fr(3) ** random.randint(-4, 4) * Fr(random.randint(1, 4)) for _ in range(k)]
    if mode == 2:   # one tiny
        w = [Fr(random.randint(1, 9)) for _ in range(k)]; w[random.randrange(k)] = Fr(1, 500); return w
    if mode == 3:   # two tiny
        w = [Fr(random.randint(1, 9)) for _ in range(k)]; w[0] = Fr(1, 300); w[-1] = Fr(1, 700); return w
    if mode == 4:   # one dominant
        w = [Fr(1, random.randint(1, 9)) for _ in range(k)]; w[random.randrange(k)] = Fr(60); return w
    return [Fr(1)] * (k - 1) + [Fr(random.randint(1, 5), 7)]


ok = True
cnt = 0; minratio = None
for k in range(3, kmax + 1):
    for i in range(npts):
        w = weights(k, i % 6)
        t = random.randrange(k)
        J = J1(w)
        st = block_states(k, (t,))
        M = rescaled(J, st)
        R = [i2 for i2, (a, y) in enumerate(st) if a is not None and a[0] == 2]
        P = [i2 for i2, (a, y) in enumerate(st) if a is not None and a[0] == 1]
        Q = [i2 for i2, (a, y) in enumerate(st) if a is None]
        s = J.s; S = J.S
        # (a)
        oka = all(M[R[a]][R[b]] == (0 if a != b else (s[a] if a == t else 2 * s[a])) for a in range(k) for b in range(k))
        Ns = [(1 if r == t else 2) * s[r] for r in range(k)]
        oka &= all(M[P[a]][P[b]] == (Ns[a] if a == b else 0) - Ns[a] * Ns[b] / (2 * S - s[t])
                   for a in range(k) for b in range(k))
        oka &= all(sum(M[i2][j] for j in P) == 0 for i2 in range(len(st)))
        # (b)
        E = R + Q
        okb = ldl_pd([[M[a][b] for b in E] for a in E])
        SP = schur(M, P, E)
        okb &= all(sum(row) == 0 for row in SP)
        offs = [SP[a][b] for a in range(k) for b in range(k) if a != b]
        okb &= all(v < 0 for v in offs)
        rat = min(-v / SP[a][a] for a in range(k) for b in range(k) if a != b for v in [SP[a][b]])
        minratio = rat if minratio is None else min(minratio, rat)
        if not (oka and okb):
            print(f'FAIL at k={k}, t={t}, w={[str(v) for v in w]}: structure {oka}, Laplacian {okb}', flush=True)
        ok &= oka and okb
        cnt += 1
    print(f'k = {k}: {npts} weight vectors done (structure and Laplacian property: {ok})', flush=True)
print(f'(a), (b): {cnt} blocks, all pass: {ok}; smallest |S_P[r,r\']|/S_P[r,r] = {float(minratio):.3e}')
# (c) negative control
found = None
for a in range(2, 80):
    for (k, t) in ((3, 0), (4, 0), (4, 1)):
        w = [Fr(1)] + [Fr(1, a * (u + 2)) for u in range(k - 1)] if t == 0 else [Fr(1, a * 3), Fr(1)] + [Fr(1, a * (u + 2)) for u in range(k - 2)]
        J = J1(w, bscale=Fr(1001, 1000))
        st = block_states(k, (t,))
        M = rescaled(J, st)
        keep = [i2 for i2, (aa, y) in enumerate(st) if not (aa is not None and aa[0] == 1 and aa[1] == (1 if t != 1 else 0))]
        if not ldl_pd([[M[i2][j] for j in keep] for i2 in keep]):
            # confirm: true block is PD in the same gauge
            J0 = J1(w); M0 = rescaled(J0, st)
            if ldl_pd([[M0[i2][j] for j in keep] for i2 in keep]):
                found = (k, t, w)
                break
    if found:
        break
print(f'(c) negative control: B -> (1 + 1/1000) B fails (gauge-fixed block not PD) at {found}')
ok &= found is not None
print('PASS' if ok else 'FAIL')
