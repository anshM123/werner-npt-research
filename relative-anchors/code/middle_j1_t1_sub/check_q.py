"""Exact checks of Proposition Q (NOTES.md, Section 2) and of Lemma K (Section 3).

(a) Reduction.  At rational weight vectors with |w| = 1 the Lemma B builder j1lib.py gives the block; on its P- and
    Q-states, with X22 = Pi - Y_PQ Y_QQ^{-1} Y_PQ^T, the matrices of the quadratic forms satisfy, exactly,
        Pi = (2 - sigma) Z^T diag(y) Z,    X22 = (2 - sigma) Z^T (diag(y) - tp^2/4 Lmat) Z,
    where zeta = Z pi is zeta_c = pi_c/w_c - m (m = sum_r N_r w_r pi_r/(2 - sigma); the copy t is the phantom copy),
    y are the phantom weights, Lmat the matrix of L and tp = 1/(1 + gp/2).
(b) Identity.  (2 + 2 gp) S - L = R1 + R2 + R3 + R4 on the constraint subspace, as a quadratic form (checked by
    polarisation on a basis), at random rational phantom weights with |U| = 1, ..., 8; every coefficient >= 0.
    For |U| = 2, 3 also symbolically (sympy, phantom weights as symbols).
(c) Conclusion.  X22 - c Pi with c = (2 + 2 gp + gp^2)/(2 + gp)^2 is >= 0 with kernel exactly span(w) (the matrix with
    the row and column of P_t deleted is positive definite, exact LDL, and (X22 - c Pi) w = 0).
(d) The block of j1lib on the P- and Q-states, with the state P_t deleted, is positive definite (exact LDL), and
    (w, 0) is in the kernel of the undeleted block, at rational weight vectors (not normalised), k = 2, ..., 8.
(e) Lemma K.  Y_PQ = Y_RP^T K exactly with the explicit K of NOTES.md, Section 3, against j1lib.
(f) Negative controls with the B-part scaled by 1001/1000 (j1lib bscale): the (P, Q) block is negative on (w, 0) at
    every point; at explicit rational points with a small weight on t the scaled block with P_t deleted is not
    positive definite, while the true one is.
Usage: python check_q.py"""
import random
import itertools
from fractions import Fraction as Fr
from j1lib import J1, block_states
from qlib import phantom, forms
from tlib import inv

random.seed(5)


def unit_weights(k):
    v = [Fr(random.randint(1, 9), random.randint(1, 9)) for _ in range(k - 1)]
    n2 = sum(a * a for a in v)
    w = [abs(2 * a / (1 + n2)) for a in v] + [abs((1 - n2) / (1 + n2))]
    if any(a == 0 for a in w):
        return unit_weights(k)
    return w


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


def mm(A, B):
    return [[sum(A[i][l] * B[l][j] for l in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def tr(A):
    return [list(r) for r in zip(*A)]


def pq_blocks(w, t, bscale=1):
    k = len(w)
    st = block_states(k, (t,))
    M = J1(w, bscale=bscale).matrix(st)
    P = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 1]
    P = sorted(P, key=lambda i: st[i][0][1])                      # P_r ordered by r
    Q = [i for i, (a, y) in enumerate(st) if a is None]
    Q = sorted(Q, key=lambda i: [c for c in range(k) if st[i][1][c] == 1][0])
    R = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 2]
    R = sorted(R, key=lambda i: [c for c in range(k) if st[i][1][c] == 1][0])
    sub = lambda I, J: [[M[i][j] for j in J] for i in I]
    return dict(Pi=sub(P, P), YPQ=sub(P, Q), YQQ=sub(Q, Q), YRP=sub(R, P), M=M, P=P, Q=Q, R=R)


def quad_matrix(f, n):
    """matrix of a quadratic form f on Q^n, by polarisation."""
    e = lambda i: [Fr(1) if j == i else Fr(0) for j in range(n)]
    d = [f(e(i)) for i in range(n)]
    A = [[Fr(0)] * n for _ in range(n)]
    for i in range(n):
        A[i][i] = d[i]
        for j in range(i + 1, n):
            v = [a + b for a, b in zip(e(i), e(j))]
            A[i][j] = A[j][i] = (f(v) - d[i] - d[j]) / 2
    return A


ok_all = True


def report(name, ok):
    global ok_all
    ok_all &= ok
    print(f'{name}: {"PASS" if ok else "FAIL"}', flush=True)


# ---------------- (a) reduction, (c) conclusion, (e) Lemma K ----------------
cnt = 0
oka = okc = oke = True
for k in range(3, 9):
    for trial in range(3):
        w = unit_weights(k)
        t = (trial * 2 + k) % k
        x = [a * a for a in w]
        sig = x[t]
        U, y = phantom(x, t)
        b = pq_blocks(w, t)
        Pi, YPQ, YQQ = b['Pi'], b['YPQ'], b['YQQ']
        X22 = [[Pi[i][j] - v for j, v in enumerate(row)] for i, row in enumerate(mm(mm(YPQ, inv(YQQ, Fr(1))), tr(YPQ)))]
        N = [Fr(1) if r == t else Fr(2) for r in range(k)]
        gp = sum(x[u] / (2 - sig - 2 * x[u]) for u in U)
        tp = 1 / (1 + gp / 2)

        def zeta_of(pi):
            m = sum(N[r] * w[r] * pi[r] for r in range(k)) / (2 - sig)
            z = {c: pi[c] / w[c] - m for c in U}
            z['ph'] = pi[t] / w[t] - m
            return z
        fS = lambda pi: forms(U, y, zeta_of(pi))['S']
        fL = lambda pi: forms(U, y, zeta_of(pi))['L']
        MS = quad_matrix(fS, k)
        ML = quad_matrix(fL, k)
        okP = all(Pi[i][j] == (2 - sig) * MS[i][j] for i in range(k) for j in range(k))
        okX = all(X22[i][j] == (2 - sig) * (MS[i][j] - tp * tp / 4 * ML[i][j]) for i in range(k) for j in range(k))
        oka &= okP and okX
        # gp of the closed form equals the phantom gp
        oka &= gp == sum(y[c] / (2 * (1 - y[c])) for c in U)
        # (c)
        c = (2 + 2 * gp + gp * gp) / (2 + gp) ** 2
        D = [[X22[i][j] - c * Pi[i][j] for j in range(k)] for i in range(k)]
        kern = all(sum(D[i][j] * w[j] for j in range(k)) == 0 for i in range(k))
        keep = [i for i in range(k) if i != t]
        okc &= kern and ldl_pd([[D[i][j] for j in keep] for i in keep])
        # (e) Lemma K
        g = [xc / (2 * (1 - xc)) for xc in x]
        tau = (1 + g[t]) / (1 + sum(g) / 2)
        K = [[Fr(0)] * len(U) for _ in range(k)]
        for j, cc in enumerate(U):
            K[cc][j] = -(2 * tp / tau) * (2 - sig) * (1 - x[cc]) / (2 - sig - 2 * x[cc])
            K[t][j] = 2 * w[t] * w[cc] * tp * (1 - sig) / ((2 - tau) * (2 - sig - 2 * x[cc]))
        YPR = tr(b['YRP'])
        oke &= mm(YPR, K) == YPQ
        cnt += 1
print(f'{cnt} weight vectors with |w| = 1, k = 3..8')
report('(a) Pi and X22 equal the phantom forms S and S - tp^2 L/4 (times 2 - sigma), exactly', oka)
report('(c) X22 - c Pi >= 0 with kernel exactly span(w), c = (2 + 2gp + gp^2)/(2 + gp)^2, exactly', okc)
report('(e) Lemma K: Y_PQ = Y_RP^T K with the explicit K, exactly', oke)

# ---------------- (b) the identity ----------------
okb = True
npts = 0
for n in range(1, 9):
    for trial in range(4):
        raw = [Fr(random.randint(1, 30), random.randint(1, 30)) for _ in range(n + 1)]
        if trial == 3:
            raw[0] = Fr(1, 10 ** 4)              # small phantom weight
        if trial == 2:
            raw[1] = Fr(10 ** 3)                 # one dominant copy of U
        tot = sum(raw)
        U = list(range(n))
        y = {c: raw[c + 1] / tot for c in U}
        y['ph'] = raw[0] / tot

        def F(v):
            z = {c: v[c] for c in U}
            z['ph'] = -sum(y[c] * v[c] for c in U) / y['ph']
            f = forms(U, y, z)
            return (2 + 2 * f['gp']) * f['S'] - f['L'] - sum(f['R'])
        A = quad_matrix(F, n)
        f0 = forms(U, y, {**{c: Fr(0) for c in U}, 'ph': Fr(0)})
        okb &= all(v == 0 for row in A for v in row) and f0['coef_ok']
        npts += 1
report(f'(b) identity (2 + 2gp) S - L = R1 + R2 + R3 + R4 as a quadratic form on the constraint subspace, {npts} '
       f'rational points, |U| = 1..8; all coefficients >= 0', okb)

try:
    import sympy as sp
    oks = True
    for n in (2, 3):
        ys = sp.symbols(f'y1:{n + 1}', positive=True)
        zs = sp.symbols(f'z1:{n + 1}')
        y0 = 1 - sum(ys)
        z0 = -sum(a * b for a, b in zip(ys, zs)) / y0
        Y = dict(zip(range(n), ys)); Y['ph'] = y0
        Z = dict(zip(range(n), zs)); Z['ph'] = z0
        lam = {c: 1 / (1 - Y[c]) for c in range(n)}
        g = {c: Y[c] / (2 * (1 - Y[c])) for c in range(n)}
        gp = sum(g.values())
        allc = ['ph'] + list(range(n))
        d = lambda a, b: Z[a] - Z[b]
        S = sum(Y[c] * Z[c] ** 2 for c in allc)
        L = sum(Y[c] * lam[c] ** 2 * Z[c] ** 2 for c in range(n)) + sp.Rational(1, 4) * sum(
            Y[c] * Y[u] / (1 - Y[c] - Y[u]) * (lam[c] * Z[c] + lam[u] * Z[u]) ** 2 for c, u in itertools.combinations(range(n), 2))
        R1 = sum(2 * Y[r] * Y[s] * (gp - g[r] - g[s]) * d(r, s) ** 2 for r, s in itertools.combinations(range(n), 2)) + \
            sum(Y[r] * y0 * (1 + 2 * gp - 2 * g[r]) * d(r, 'ph') ** 2 for r in range(n))
        R2 = sp.Rational(1, 2) * sum(Y[l] * lam[l] ** 2 * sum(Y[a] * y0 * d(a, 'ph') ** 2 for a in range(n) if a != l) for l in range(n))
        R3 = 0
        R4 = 0
        for r, s in itertools.combinations(range(n), 2):
            others = [u for u in allc if u not in (r, s)]
            R3 += sp.Rational(1, 4) * sum(Y[r] * Y[s] * Y[u] * (lam[s] * d(r, u) - lam[r] * d(s, u)) ** 2 for u in others)
            R4 += sp.Rational(1, 4) * Y[r] * Y[s] / (1 - Y[r] - Y[s]) * sum(
                Y[u] * Y[v] * (lam[r] + lam[s]) ** 2 * d(u, v) ** 2 for u, v in itertools.combinations(others, 2))
        expr = sp.together((2 + 2 * gp) * S - L - R1 - R2 - R3 - R4)
        num = sp.expand(sp.numer(expr))
        oks &= num == 0
        print(f'  symbolic, |U| = {n}: numerator of the difference is {num}', flush=True)
    report('(b) identity symbolically in the phantom weights, |U| = 2, 3', oks)
except ImportError:
    report('(b) symbolic check skipped (no sympy)', False)

# ---------------- (d) the (P, Q) block itself, (f) negative controls ----------------
okd = True
nd = 0
for k in range(2, 9):
    for trial in range(3):
        w = [Fr(random.randint(1, 12), random.randint(1, 12)) for _ in range(k)]
        t = random.randrange(k)
        b = pq_blocks(w, t)
        idx = b['P'] + b['Q']
        M = b['M']
        sub = [[M[i][j] for j in idx] for i in idx]
        vec = [w[r] for r in range(k)] + [Fr(0)] * len(b['Q'])
        kern = all(sum(sub[i][j] * vec[j] for j in range(len(idx))) == 0 for i in range(len(idx)))
        keep = [i for i in range(len(idx)) if i != t]
        okd &= kern and ldl_pd([[sub[i][j] for j in keep] for i in keep])
        nd += 1
report(f'(d) (P, Q) block of j1lib: (w, 0) in the kernel and positive definite with P_t deleted, {nd} weight vectors, k = 2..8', okd)

okf = True
pts = [([Fr(1, 300), Fr(1), Fr(1)], 0), ([Fr(1, 200), Fr(1), Fr(1), Fr(1, 3)], 0), ([Fr(1, 2), Fr(1, 500), Fr(1), Fr(2)], 1),
       ([Fr(1), Fr(2), Fr(1, 400), Fr(1), Fr(1, 2)], 2)]
for w, t in pts:
    k = len(w)
    vals = []
    for bs in (1, Fr(1001, 1000)):
        b = pq_blocks(w, t, bscale=bs)
        idx = b['P'] + b['Q']
        M = b['M']
        sub = [[M[i][j] for j in idx] for i in idx]
        vec = [w[r] for r in range(k)] + [Fr(0)] * len(b['Q'])
        val = sum(vec[i] * sub[i][j] * vec[j] for i in range(len(idx)) for j in range(len(idx)))
        keep = [i for i in range(len(idx)) if i != t]
        pd = ldl_pd([[sub[i][j] for j in keep] for i in keep])
        vals.append((val, pd))
    (v1, pd1), (v2, pd2) = vals
    good = v1 == 0 and pd1 and v2 < 0 and not pd2
    print(f'  w = {[str(a) for a in w]}, t = {t}: true block value on (w,0) = {v1}, gauge-fixed PD = {pd1}; '
          f'scaled: value {float(v2):.3e}, gauge-fixed PD = {pd2}', flush=True)
    okf &= good
report('(f) negative controls (B scaled by 1001/1000 fails, true block passes) at 4 explicit points', okf)

print('ALL PASS' if ok_all else 'SOME FAIL')
