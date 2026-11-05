"""Exact checks (python-flint) for Theorem 12 and Proposition 11 of NOTES.md, Section 5.
(1) g_m = sum_a |a..a> (m = 2, 3): A - B - P/64 >= 0 (kernel exactly g (x) C^{3^m}) and B + (215/216) A > 0,
    with A = (M^{(x)m})^Gamma, B = R_g^Gamma, P = (I - |g><g|/3) (x) I.  Exact LDL on connected components.
(2) b = |10> + |01> (unnormalised): D_b spectrum {1/2 x3, 3/4 x5, 1 x1} (b normalised); generalized eigenvector v with
    B_b v = -(11/9) A_2 v, |v|^2 = 28, <v, A_2 v> = 9.
(3) E_0 = span{|01>-|10>, |02>-|20>} (normalised rows): compression of M_2 to E_0 (x) C^9 has spectrum
    {1/4 x3, 1/2 x2, 3/4 x9, 1 x4}."""
import itertools, flint
fq = flint.fmpq
H = fq(1, 2)
def m1(a, y, b, x): return (1 if (a == b and y == x) else 0) - (H if (a == y and b == x) else 0)
def Mel(a, y, b, x):
    v = fq(1)
    for s in range(len(a)):
        v *= m1(a[s], y[s], b[s], x[s])
        if v == 0: return v
    return v
def build(m, rows):
    """rows: dict string -> coefficient (the row vector p).  Returns states, A, B as flint fmpq_mat."""
    cols = list(itertools.product(range(3), repeat=m)); N = len(cols)
    st = [(a, y) for a in cols for y in cols]; idx = {s: i for i, s in enumerate(st)}
    D = flint.fmpq_mat(N, N)
    for i, x in enumerate(cols):
        for j, xp in enumerate(cols):
            D[i, j] = sum((c1 * c2 * Mel(r1, x, r2, xp) for r1, c1 in rows.items() for r2, c2 in rows.items()), fq(0))
    Dinv = D.inv()
    G = flint.fmpq_mat(len(st), N)
    for i, (a, y) in enumerate(st):
        for j, x in enumerate(cols):
            G[i, j] = sum((c * Mel(a, y, r, x) for r, c in rows.items()), fq(0))
    R = G * Dinv * G.transpose()
    A = flint.fmpq_mat(len(st), len(st)); B = flint.fmpq_mat(len(st), len(st))
    for i, (a, yp) in enumerate(st):
        for j, (b, y) in enumerate(st):
            A[i, j] = Mel(b, yp, a, y)
            B[i, j] = R[idx[(b, yp)], idx[(a, y)]]
    return st, A, B, D
def components(Mx):
    n = Mx.nrows(); seen = [False] * n; comps = []
    for s in range(n):
        if seen[s]: continue
        stack = [s]; seen[s] = True; c = []
        while stack:
            i = stack.pop(); c.append(i)
            for j in range(n):
                if not seen[j] and Mx[i, j] != 0: seen[j] = True; stack.append(j)
        comps.append(sorted(c))
    return comps
def ldl_psd(Mx):
    """exact symmetric elimination; returns (psd, rank)."""
    n = Mx.nrows(); S = [[Mx[i, j] for j in range(n)] for i in range(n)]; rank = 0
    for k in range(n):
        p = S[k][k]
        if p < 0: return False, rank
        if p == 0:
            if any(S[k][j] != 0 for j in range(k + 1, n)): return False, rank
            continue
        rank += 1
        for i in range(k + 1, n):
            if S[i][k] == 0: continue
            f = S[i][k] / p
            for j in range(k + 1, n): S[i][j] -= f * S[k][j]
    return True, rank
def psd_by_components(Mx):
    tot = 0; ok = True; sizes = {}
    for c in components(Mx):
        sub = flint.fmpq_mat(len(c), len(c))
        for a, i in enumerate(c):
            for b, j in enumerate(c): sub[a, b] = Mx[i, j]
        good, r = ldl_psd(sub); ok &= good; tot += r; sizes[len(c)] = sizes.get(len(c), 0) + 1
    return ok, tot, sizes
ok_all = True
for m in (2, 3):
    rows = {tuple([a] * m): fq(1) for a in range(3)}
    st, A, B, D = build(m, rows)
    n = len(st); P = flint.fmpq_mat(n, n)
    for i, (a, yp) in enumerate(st):
        for j, (b, y) in enumerate(st):
            if yp == y:
                P[i, j] = (1 if a == b else 0) - (fq(1, 3) if (a in rows and b in rows) else 0)
    gap = A - B - P * fq(1, 64); low = B + A * fq(215, 216)
    ok1, r1, sz = psd_by_components(gap); ok2, r2, _ = psd_by_components(low)
    print(f"m={m}: A - B - P/64 PSD {ok1}, rank {r1} (nullity {n - r1}); B + 215/216 A PSD {ok2}, rank {r2} of {n}; component sizes {sz}", flush=True)
    ok_all &= ok1 and ok2 and n - r1 == 3 ** m and r2 == n
# (2) b
rows = {(1, 0): fq(1), (0, 1): fq(1)}
st, A2, Bb, Db = build(2, rows)
Dn = Db * fq(1, 2)    # b normalised: |b|^2 = 2
cp = Dn.charpoly()
print("D_b (normalised) charpoly:", cp.factor())
idx = {s: i for i, s in enumerate(st)}
v = flint.fmpq_mat(81, 1)
for (r, c, val) in [((0, 0), (1, 1), 2), ((0, 1), (0, 1), 3), ((0, 1), (1, 0), 1), ((1, 0), (0, 1), 1), ((1, 0), (1, 0), 3), ((1, 1), (0, 0), 2)]:
    v[idx[(r, c)], 0] = val
lhs = Bb * v; rhs = A2 * v * fq(-11, 9)
nv = (v.transpose() * v)[0, 0]; av = (v.transpose() * A2 * v)[0, 0]
print(f"B_b v = -(11/9) A_2 v: {lhs == rhs};  |v|^2 = {nv}, <v, A_2 v> = {av}")
ok_all &= lhs == rhs and nv == 28 and av == 9
# (3) E_0 compression of M_2
cols = list(itertools.product(range(3), repeat=2))
e1 = {(0, 1): fq(1), (1, 0): fq(-1)}; e2 = {(0, 2): fq(1), (2, 0): fq(-1)}   # norms^2 = 2
Q = flint.fmpq_mat(18, 18)
for i in range(18):
    for j in range(18):
        u, y = (e1 if i < 9 else e2), cols[i % 9]; up, yp = (e1 if j < 9 else e2), cols[j % 9]
        Q[i, j] = sum((c1 * c2 * Mel(r1, y, r2, yp) for r1, c1 in u.items() for r2, c2 in up.items()), fq(0)) / 2
print("E_0 compression charpoly:", Q.charpoly().factor())
print("PASS" if ok_all else "FAIL")
