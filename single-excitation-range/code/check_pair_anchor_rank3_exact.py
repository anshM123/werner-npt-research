"""Exact proof (README, Proposition 10 (a)) that every two-copy vector phi = |00> + b|11> + c|22> (b, c > 0) is a universal anchor with
-3 <= T_phi <= I, i.e. A - R^Gamma >= 0 and R^Gamma + 3A >= 0 on the 81-dimensional row (x) column space, for all
b, c > 0 (open quadrant; boundary by continuity).  Exact arithmetic over Q(b, c) (python-flint via flintfield.py).
Blocks: connected components of the exact nonzero pattern of A and R^Gamma.  PSD certificate: symmetric Gaussian
elimination with pivots that are ratios of polynomials with one-sign coefficients (positive on the open quadrant)."""
import itertools, sys, time
from flintfield import Field
t0 = time.time()
K = Field(['b', 'c']); bb, cc = K.gens; ONE, ZERO = K(1), K(0); half = K(1) / K(2)
bscale = K(1001) / K(1000) if 'sanity' in sys.argv else ONE
def m1(a, y, ap, yp): return (ONE if (a == ap and y == yp) else ZERO) - (half if (a == y and ap == yp) else ZERO)
S2 = list(itertools.product(range(3), repeat=2))
def Mel(row, col, rowp, colp): return m1(row[0], col[0], rowp[0], colp[0]) * m1(row[1], col[1], rowp[1], colp[1])
phi = {(0, 0): ONE, (1, 1): bb, (2, 2): cc}
full = [(r, c) for r in S2 for c in S2]
MV = {(r, c): {z: sum((Mel(r, c, a, z) * pa for a, pa in phi.items()), ZERO) for z in S2} for (r, c) in full}
D = [[sum((pa * pb * Mel(a, y, b_, yp) for a, pa in phi.items() for b_, pb in phi.items()), ZERO) for yp in S2] for y in S2]
def inv(Mx):
    n = len(Mx); A = [row[:] + [ONE if i == j else ZERO for j in range(n)] for i, row in enumerate(Mx)]
    for col in range(n):
        r = next(i for i in range(col, n) if A[i][col] != 0); A[col], A[r] = A[r], A[col]
        pv = A[col][col]; A[col] = [x / pv for x in A[col]]
        for i in range(n):
            if i != col and A[i][col] != 0:
                f = A[i][col]; A[i] = [x - f * y for x, y in zip(A[i], A[col])]
    return [row[n:] for row in A]
Di = inv(D)
def Rel(b_, yp, a, y):   # R[(b,y'),(a,y)] = sum_{z,z'} MV[(b,y'),z] Dinv[z,z'] MV[(a,y),z']  (real parameters)
    r1, r2 = MV[(b_, yp)], MV[(a, y)]
    return sum((r1[S2[i]] * Di[i][j] * r2[S2[j]] for i in range(9) for j in range(9) if r1[S2[i]] != 0 and r2[S2[j]] != 0), ZERO)
idx = full
A = [[Mel(b_, yp, a, y) for (b_, y) in idx] for (a, yp) in idx]           # A[(a,y'),(b,y)] = M[(b,y'),(a,y)]
B = [[Rel(b_, yp, a, y) for (b_, y) in idx] for (a, yp) in idx]           # B = R^Gamma
print(f"built exact A, R^Gamma (81 x 81)  [{time.time()-t0:.1f}s]", flush=True)
def onesign(poly):
    cs = poly.coeffs(); return all(x > 0 for x in cs) or all(x < 0 for x in cs)
def positive(f):
    if f == 0: return False
    return onesign(f.numer) and onesign(f.denom) and ((f.numer.coeffs()[0] > 0) == (f.denom.coeffs()[0] > 0))
def psd_cert(Mx):
    M = [row[:] for row in Mx]; act = list(range(len(M))); npiv = 0
    while act:
        for i in act:
            if M[i][i] == 0 and any(M[i][j] != 0 for j in act): return False, None
        cand = [i for i in act if M[i][i] != 0]
        if not cand: return True, npiv
        cand.sort(key=lambda i: len(M[i][i].numer.terms()) + len(M[i][i].denom.terms()))
        i = next((c_ for c_ in cand if positive(M[c_][c_])), None)
        if i is None: return False, None
        npiv += 1; rest = [j for j in act if j != i]
        for j in rest:
            if M[j][i] == 0: continue
            f = M[j][i] / M[i][i]
            for l in rest:
                if M[i][l] != 0: M[j][l] -= f * M[i][l]
        act = rest
    return True, npiv
ok = True
for name, coef in (("A - R^Gamma", -bscale), ("R^Gamma + 3A", None)):
    Mt = [[(A[i][j] + coef * B[i][j]) if coef is not None else (B[i][j] + 3 * A[i][j]) for j in range(81)] for i in range(81)]
    sym = all(Mt[i][j] == Mt[j][i] for i in range(81) for j in range(81))
    # connected components of the nonzero pattern
    seen = [False] * 81; comps = []
    for s in range(81):
        if seen[s]: continue
        st = [s]; seen[s] = True; comp = []
        while st:
            i = st.pop(); comp.append(i)
            for j in range(81):
                if not seen[j] and (Mt[i][j] != 0): seen[j] = True; st.append(j)
        comps.append(sorted(comp))
    allok = sym; rk = 0
    for comp in comps:
        g, r = psd_cert([[Mt[i][j] for j in comp] for i in comp]); allok &= g; rk += r or 0
        if not g: print("   not certified: component of size", len(comp))
    print(f"{name}: symmetric {sym}, {len(comps)} components (largest {max(map(len, comps))}), certified PSD on the open quadrant: {allok}, rank {rk}  [{time.time()-t0:.1f}s]", flush=True)
    ok &= allok
print("PASS" if ok else "FAIL")
