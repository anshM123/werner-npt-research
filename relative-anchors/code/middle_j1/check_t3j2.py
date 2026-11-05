"""Checks for Theorem D2 (NOTES.md, Section 5): the blocks of (R_inf)_2 with |T| = 3 (no letter 1).
(a) The closed form (Lemma D2.1) equals the block computed from Lemma B (j1lib.entry, which is Lemma B for every j),
    exactly at rational weights, k = 4..8, several T; for k <= 6 that block also equals ylib's (the rescaling is
    trivial, since there is no letter 1).
(b) At random rational configurations (k = 4..9): the moment bounds e*l_d(Z) <= F <= e*h_d(Z) for the six sums,
    C > 0, and f = 2 - G_T a^T C^{-1} a >= min over the 64 vertices > 0 (exact).
(c) Negative control: with B scaled by 1 + 1/1000 the block is indefinite at an explicit rational point (exact).
(d) The block is positive definite (exact LDL) at 40 rational weight vectors, k = 4..9.
Usage: python check_t3j2.py"""
import itertools, random
from fractions import Fraction as Fr
import flint
from j1lib import J1, onehot, to_ylib_state
from ylib import Builder, point_ring

random.seed(31)


def states(k, T):
    st = []
    for t in T:
        rest = [u for u in T if u != t]
        st.append(((2, t), onehot(k, {rest[0]: 2, rest[1]: 2})))
    st.append((None, onehot(k, {u: 2 for u in T})))
    return st


def parts(x, T):
    k = len(x)
    U = [u for u in range(k) if u not in T]
    xs = [x[t] for t in T]
    xT = sum(xs)
    e = 1 - xT
    out = dict(xT=xT, e=e, U=U)
    c = {}
    for i, j in itertools.combinations(range(3), 2):
        kk = 3 - i - j
        g = lambda u: u / (2 - 2 * u - xs[kk])
        F = sum(g(x[u]) for u in U)
        G0 = 1 + (g(xs[i]) + g(xs[j])) / 2
        c[(i, j)] = (1 + g(xs[i])) * (1 + g(xs[j])) / (2 * (G0 + F / 2))
        out[('F', i, j)] = (F, xT - xs[kk], (1 + g(xs[i])) * (1 + g(xs[j])), G0)
    A = []
    for i in range(3):
        sh = xT - xs[i]
        gp = xs[i] / (2 - 2 * xs[i] - sh)
        F = sum(x[u] / (2 - 2 * x[u] - sh) for u in U)
        A.append((1 + gp) / (1 + (gp + F) / 2))
        out[('Fp', i)] = (F, xs[i], gp)
    FT = sum(x[u] / (2 - 2 * x[u] - xT) for u in U)
    GT = 1 + FT / 2
    M = [[Fr(2), c[(0, 1)], c[(0, 2)], -A[0]], [c[(0, 1)], Fr(2), c[(1, 2)], -A[1]],
         [c[(0, 2)], c[(1, 2)], Fr(2), -A[2]], [-A[0], -A[1], -A[2], 2 / GT]]
    out.update(M=M, c=c, A=A, FT=FT, GT=GT)
    return out


def det3(C):
    return (C[0][0] * (C[1][1] * C[2][2] - C[1][2] ** 2) - C[0][1] * (C[0][1] * C[2][2] - C[1][2] * C[0][2])
            + C[0][2] * (C[0][1] * C[1][2] - C[1][1] * C[0][2]))


def fval(c12, c13, c23, a, GT):
    C = [[Fr(2), c12, c13], [c12, Fr(2), c23], [c13, c23, Fr(2)]]
    adj = [[4 - c23 * c23, c13 * c23 - 2 * c12, c12 * c23 - 2 * c13],
           [c13 * c23 - 2 * c12, 4 - c13 * c13, c12 * c13 - 2 * c23],
           [c12 * c23 - 2 * c13, c12 * c13 - 2 * c23, 4 - c12 * c12]]
    d = det3(C)
    q = sum(a[i] * adj[i][j] * a[j] for i in range(3) for j in range(3))
    return d, 2 - GT * q / d


def pd_exact(M):
    A = [row[:] for row in M]
    n = len(A)
    for c in range(n):
        if A[c][c] <= 0:
            return False
        for i in range(c + 1, n):
            f = A[i][c] / A[c][c]
            for j in range(c + 1, n):
                A[i][j] -= f * A[c][j]
    return True


def rand_w(k, mode):
    if mode == 0:
        return [Fr(random.randint(1, 30), random.randint(1, 30)) for _ in range(k)]
    if mode == 1:
        return [Fr(2) ** random.randint(-6, 6) for _ in range(k)]
    if mode == 2:
        w = [Fr(1, random.randint(20, 200)) for _ in range(k)]; w[0] = Fr(1); return w
    w = [Fr(1, random.randint(30, 300)) for _ in range(k)]
    w[0] = Fr(random.randint(1, 3)); w[1] = Fr(random.randint(1, 3)); w[2] = Fr(random.randint(1, 3))
    return w


ok = True
# (a)
good = True; cnt = 0
for k in range(4, 9):
    w = rand_w(k, k % 4)
    J = J1(w)
    x = [s / J.S for s in J.s]
    Bd = Builder(point_ring([flint.fmpq(v.numerator ** 2, v.denominator ** 2) for v in w])) if k <= 6 else None
    for T in [(0, 1, 2), (1, 2, 3), (0, 2, k - 1)]:
        st = states(k, T)
        M = [[J.entry(a, b) for b in st] for a in st]
        P = parts(x, T)
        good &= all(M[i][j] == P['M'][i][j] for i in range(4) for j in range(4))
        if Bd is not None:
            ys = [to_ylib_state(z, k) for z in st]
            for i in range(4):
                for j in range(4):
                    yv = Bd.Y(ys[i], ys[j])
                    good &= Fr(int(yv.p), int(yv.q)) == M[i][j] / 2 ** (k - 2)
        cnt += 1
print(f"(a) closed form = Lemma B block = ylib block, exactly, {cnt} blocks (k = 4..8): {good}")
ok &= good


def h(d, z):
    return z / (1 + d * z)


# (b)
goodb = True; n = 0
for trial in range(200):
    k = random.randint(4, 9)
    w = rand_w(k, trial % 4)
    s = [v * v for v in w]; S = sum(s); x = [v / S for v in s]
    P = parts(x, (0, 1, 2))
    xT, e = P['xT'], P['e']
    Z = P['FT'] / e
    Z0, Z1 = 1 / (xT + 2 * e), 1 / xT
    lam = (Z - Z0) / (Z1 - Z0)
    goodb &= 0 <= lam <= 1
    crng = {}
    for (i, j) in [(0, 1), (0, 2), (1, 2)]:
        F, d, num, G0 = P[('F', i, j)]
        lo, hi = e * ((1 - lam) * h(d, Z0) + lam * h(d, Z1)), e * h(d, Z)
        goodb &= lo <= F <= hi
        crng[(i, j)] = [num / (2 * (G0 + FF / 2)) for FF in (lo, hi)]
    arng = []
    for i in range(3):
        F, d, gp = P[('Fp', i)]
        lo, hi = e * ((1 - lam) * h(d, Z0) + lam * h(d, Z1)), e * h(d, Z)
        goodb &= lo <= F <= hi
        arng.append([(1 + gp) / (1 + (gp + FF) / 2) for FF in (lo, hi)])
    dtrue, ftrue = fval(P['c'][(0, 1)], P['c'][(0, 2)], P['c'][(1, 2)], P['A'], P['GT'])
    vals = []
    for cv in itertools.product((0, 1), repeat=3):
        for av in itertools.product((0, 1), repeat=3):
            d, f = fval(crng[(0, 1)][cv[0]], crng[(0, 2)][cv[1]], crng[(1, 2)][cv[2]], [arng[i][av[i]] for i in range(3)], P['GT'])
            goodb &= d > 0
            vals.append(f)
    goodb &= dtrue > 0 and ftrue >= min(vals) > 0
    n += 1
print(f"(b) moment bounds, C > 0, and f >= min of the 64 vertex values > 0 at {n} rational configurations: {goodb}")
ok &= goodb
# (c)
found = None
cands = [w for a in (5, 10, 20, 40, 80, 160, 320, 640, 1280)
         for w in ([Fr(1), Fr(1, a), Fr(1, a + 1), Fr(1, a + 2)], [Fr(1), Fr(1, a), Fr(1, a), Fr(1, a)],
                   [Fr(1), Fr(1), Fr(1, a), Fr(1, a)])]
for w in cands:
    st = states(4, (0, 1, 2))
    Js = J1(w, bscale=Fr(1001, 1000)); J0 = J1(w)
    Ms = [[Js.entry(p, q) for q in st] for p in st]; M0 = [[J0.entry(p, q) for q in st] for p in st]
    if not pd_exact(Ms) and pd_exact(M0):
        found = [str(v) for v in w]
        break
print(f"(c) negative control: B -> (1 + 1/1000) B indefinite at k = 4, T = (0, 1, 2), w = {found}")
ok &= found is not None
# (d)
goodd = True
for trial in range(40):
    k = random.randint(4, 9)
    w = rand_w(k, trial % 4)
    J = J1(w)
    st = states(k, (0, 1, 2))
    goodd &= pd_exact([[J.entry(p, q) for q in st] for p in st])
print(f"(d) positive definite at 40 rational weight vectors (k = 4..9): {goodd}")
ok &= goodd
print("PASS" if ok else "FAIL")
