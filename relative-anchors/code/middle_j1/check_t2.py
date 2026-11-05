"""Checks for Theorem D1 (NOTES.md, Section 2): the type I blocks with |T| = 2 of (R_inf)_1.
(a) The closed form (Lemma D1.1) equals the block computed from Lemma B (j1lib), exactly, at rational weights,
    k = 3..8, several pairs T.
(b) The moment bounds (Lemma D1.3): e*chord_i(Z) <= F_i <= e*h_i(Z), Z = F3/e, exactly at random rational
    configurations (k = 3..9), including spread and near-degenerate ones.
(c) Lemmas D1.2 and D1.4: at the same configurations, c < 2, det M = 2f/G3, and f(true) >= min over the 8 vertex
    values at the true lam, exactly.
(d) Negative control: with the B-part scaled by 1 + 1/1000 the block is indefinite at an explicit rational weight
    vector (exact; the scaled block is also compared with ylib's A - (1001/1000) B at that point).
(e) The block is positive definite (exact LDL) at 60 rational weight vectors, k = 3..9, including extreme ones.
Usage: python check_t2.py"""
import random, itertools
from fractions import Fraction as Fr
import flint
from j1lib import J1, block_states, to_ylib_state, mu
from ylib import Builder, point_ring

random.seed(2026)
ok = True


def closed(x, t, tp):
    """closed form, x = s/|s| (Fractions, sum 1); order (2@t,2@t'), (2@t',2@t), (vac,2@t 2@t')."""
    k = len(x)
    al, be = x[t], x[tp]
    oth = [u for u in range(k) if u not in (t, tp)]
    ga = al / (2 * (1 - al)); gb = be / (2 * (1 - be))
    F0 = sum(x[u] / (2 - 2 * x[u]) for u in oth)
    F1 = sum(x[u] / (2 - 2 * x[u] - be) for u in oth)
    F2 = sum(x[u] / (2 - 2 * x[u] - al) for u in oth)
    F3 = sum(x[u] / (2 - 2 * x[u] - al - be) for u in oth)
    kap = (1 + ga) * (1 + gb); Gam = 1 + (ga + gb + F0) / 2
    ap = 1 + al / (2 - 2 * al - be); bp = 1 + be / (2 - 2 * be - al)
    Gp = (1 + ap) / 2 + F1 / 2; Gpp = (1 + bp) / 2 + F2 / 2; G3 = 1 + F3 / 2
    c = kap / (2 * Gam); A = ap / Gp; B = bp / Gpp
    M = [[Fr(2), c, -A], [c, Fr(2), -B], [-A, -B, 2 / G3]]
    return M, dict(al=al, be=be, F=(F0, F1, F2, F3), c=c, A=A, B=B, G3=G3, kap=kap, ap=ap, bp=bp,
                   G0=1 + (ga + gb) / 2)


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
    if mode == 1:   # spread
        return [Fr(2) ** random.randint(-6, 6) * Fr(random.randint(1, 5)) for _ in range(k)]
    if mode == 2:   # one dominant T-copy
        w = [Fr(1, random.randint(20, 200)) for _ in range(k)]
        w[1] = Fr(1)
        return w
    w = [Fr(1, random.randint(30, 300)) for _ in range(k)]   # both T copies dominant
    w[0] = Fr(random.randint(1, 3)); w[1] = Fr(random.randint(1, 3))
    return w


# (a)
cnt = 0; good = True
for k in range(3, 9):
    for trial in range(2 if k <= 6 else 1):
        w = rand_w(k, trial % 4)
        J = J1(w)
        x = [s / J.S for s in J.s]
        for (t, tp) in [(0, 1), (1, 2), (0, k - 1)]:
            st = block_states(k, (t, tp))
            M = J.matrix(st)
            C, _ = closed(x, t, tp)
            good &= all(M[i][j] == C[i][j] for i in range(3) for j in range(3))
            cnt += 1
print(f"(a) closed form = Lemma B block, exactly, {cnt} blocks (k = 3..8): {good}")
ok &= good

# (b), (c)
def h(d, z):
    return z / (1 + d * z)


nb = 0; goodb = True; goodc = True
for trial in range(300):
    k = random.randint(3, 9)
    w = rand_w(k, trial % 4)
    s = [v * v for v in w]
    S = sum(s)
    x = [v / S for v in s]
    M, P = closed(x, 0, 1)
    al, be = P['al'], P['be']
    e = 1 - al - be
    F0, F1, F2, F3 = P['F']
    Z = F3 / e
    Z0, Z1 = 1 / (al + be + 2 * e), 1 / (al + be)
    lam = (Z - Z0) / (Z1 - Z0)
    goodb &= (0 <= lam <= 1)
    rngF = []
    for Fi, d in zip((F0, F1, F2), (al + be, al, be)):
        lo = e * ((1 - lam) * h(d, Z0) + lam * h(d, Z1)); hi = e * h(d, Z)
        goodb &= (lo <= Fi <= hi)
        rngF.append((lo, hi))
    nb += 1
    # (c)
    G3 = 1 + F3 / 2
    goodc &= P['c'] < 2
    ftrue = 4 - P['c'] ** 2 - G3 * (P['A'] ** 2 + P['B'] ** 2 - P['c'] * P['A'] * P['B'])
    vals = []
    for i0, i1, i2 in itertools.product((0, 1), repeat=3):
        c = P['kap'] / (2 * (P['G0'] + rngF[0][i0] / 2))
        A = P['ap'] / ((1 + P['ap']) / 2 + rngF[1][i1] / 2)
        B = P['bp'] / ((1 + P['bp']) / 2 + rngF[2][i2] / 2)
        vals.append(4 - c * c - G3 * (A * A + B * B - c * A * B))
    goodc &= ftrue >= min(vals) and min(vals) > 0
    # determinant identity: det M = (2/G3) * ftrue... check det M = 2 f / G3
    det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] ** 2) - M[0][1] * (M[0][1] * M[2][2] - M[1][2] * M[0][2])
           + M[0][2] * (M[0][1] * M[1][2] - M[1][1] * M[0][2]))
    goodc &= det == 2 * ftrue / G3
print(f"(b) moment bounds at {nb} rational configurations: {goodb}")
print(f"(c) c < 2, det = 2f/G3, and f >= min of the 8 vertex values > 0 at the same configurations: {goodc}")
ok &= goodb and goodc

# (d) negative control
found = None
for k in (3, 4):
    for a in range(2, 60):
        w = [Fr(1, a * 7), Fr(1), Fr(1, a * 5)] + [Fr(1, a * 9)] * (k - 3)
        J = J1(w, bscale=Fr(1001, 1000))
        st = block_states(k, (0, 1))
        M = J.matrix(st)
        if not pd_exact(M):
            found = (k, w, M, st)
            break
    if found:
        break
goodd = found is not None
if found:
    k, w, M, st = found
    # compare with ylib A - 1.001 B (rescaled basis; mu = 1 here since there is no letter 1)
    Bd = Builder(point_ring([flint.fmpq(v.numerator ** 2, v.denominator ** 2) for v in w]))
    ys = [to_ylib_state(z, k) for z in st]
    same = True
    for i in range(3):
        for j in range(3):
            av = Bd.A(ys[i][0], ys[i][1], ys[j][0], ys[j][1]); bv = Bd.B(ys[i][0], ys[i][1], ys[j][0], ys[j][1])
            yv = av - flint.fmpq(1001, 1000) * bv
            same &= Fr(int(yv.p), int(yv.q)) == M[i][j] / 2 ** (k - 1)
    J0 = J1(w)
    goodd &= same and pd_exact(J0.matrix(st))
    print(f"(d) negative control: B -> (1 + 1/1000) B is indefinite at k = {k}, w = {[str(v) for v in w]} "
          f"(the true block is PD there; scaled block equals ylib A - 1.001 B: {same})")
else:
    print("(d) negative control: no witness found")
ok &= goodd

# (e)
gooda = True
for trial in range(60):
    k = random.randint(3, 9)
    w = rand_w(k, trial % 4)
    J = J1(w)
    st = block_states(k, (0, 1))
    gooda &= pd_exact(J.matrix(st))
print(f"(e) block positive definite at 60 rational weight vectors (k = 3..9): {gooda}")
ok &= gooda
print("PASS" if ok else "FAIL")
