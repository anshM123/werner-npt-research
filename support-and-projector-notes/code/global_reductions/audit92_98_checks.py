"""Checks of statements from an independent audit (rounds 92-98), referenced from OBSTRUCTIONS.md,
bosonic-sector-notes/NOTES.md and whole-problem-routes/README.md.  Every printed line ends in ok or FAIL.
 (92) sharp correlation bound ||T||_* <= (k+1)/2 - 1/d; the k-positive ball Phi_{r,O} and its absolute value;
      the separable identity for E_s; the integer cutoff 18*5^16 < 6^16; GHZ outputs of L_3^(x)n.
 (93) the tetrahedral matrix-product planes J_n(p): transfer identity, strongly connected components, the return
      recurrence, odd-n positivity (and what happens at even n); the four-triple transfer identity.
 (94) the full-support two-copy family against product-spectrum majorization; positivity margins; p_*; b(p).
 (95) the covariant completions F_u, the balanced-cut eigenvalue, the 6-copy value at u = 2, the integer
      certificates, and the a = 5/6 completion (Choi spectra, Kraus form, EB decompositions, 2-EB, 2-positivity).
 (96) the regularised local form R_alpha: spectra, inverse, three-copy values, flag lift at r = 1, 2 and r = 6.
 (97) closed bosonic spectrum, its minimum -1/(2n), hook bound, the twirl examples and the full-group family.
 (98) the rank-two calibrations, polarisation identity and the spectator-encoding coefficients.
usage: python audit92_98_checks.py"""
import itertools
import math
from fractions import Fraction as Fr
import numpy as np
from scipy.linalg import sqrtm
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

rng = np.random.default_rng(9298)
I3 = np.eye(3)
RESULTS = []


def report(label, cond, extra=""):
    RESULTS.append(bool(cond))
    print(f"{label}: {extra} {'ok' if cond else 'FAIL'}")


def site_apply(X, n, f, sites=None):
    T = np.asarray(X, dtype=complex).reshape([3] * (2 * n))
    for s in (range(n) if sites is None else sites):
        Tp = np.moveaxis(T, [s, n + s], [2 * n - 2, 2 * n - 1])
        T = np.moveaxis(f(Tp), [2 * n - 2, 2 * n - 1], [s, n + s])
    return T.reshape(3 ** n, 3 ** n)


def lin(a=0.0, b=0.0, c=0.0, z=0.0):
    """X -> a Tr(X) I + b X + c X^T + z diag(X), acting on the last two axes."""
    def f(Tp):
        tr = np.trace(Tp, axis1=-2, axis2=-1)
        return b * Tp + c * np.swapaxes(Tp, -1, -2) + z * (Tp * I3) + a * tr[..., None, None] * I3
    return f


RW = lin(a=-0.5, b=1.0)                      # Werner endpoint R = id - D/2


def qform(C, n, f=RW, C2=None):
    return np.vdot(C, site_apply(C if C2 is None else C2, n, f))


def ptrans(X, n, S):
    T = X.reshape([3] * (2 * n)); perm = list(range(2 * n))
    for s in S: perm[s], perm[n + s] = n + s, s
    return np.transpose(T, perm).reshape(3 ** n, 3 ** n)


def permute_sites(X, n, pi):
    return np.transpose(X.reshape([3] * (2 * n)), list(pi) + [n + p for p in pi]).reshape(3 ** n, 3 ** n)


def twirl(X, n, group=None):
    group = list(itertools.permutations(range(n))) if group is None else group
    return sum(permute_sites(X, n, pi) for pi in group) / len(group)


def superop(f):
    """9 x 9 matrix of a map on 3 x 3 matrices (column-stacking of row-major vec)."""
    M = np.zeros((9, 9), dtype=complex)
    for a in range(3):
        for b in range(3):
            E = np.zeros((1, 1, 3, 3)); E[0, 0, a, b] = 1
            M[:, 3 * a + b] = f(E.astype(complex)).reshape(9)
    return M


def choi(f):
    out = np.zeros((9, 9), dtype=complex)
    for a in range(3):
        for b in range(3):
            E = np.zeros((1, 1, 3, 3)); E[0, 0, a, b] = 1
            out += np.kron(np.outer(I3[a], I3[b]), f(E.astype(complex)).reshape(3, 3))
    return out


F9 = np.zeros((9, 9)); Om = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        F9[3 * i + j, 3 * j + i] = 1; Om[4 * i, 4 * j] = 1

# ============================================================================ 92
print("== round 92")


def sym_basis(d):
    B = []
    for i in range(d):
        for j in range(i + 1, d):
            M = np.zeros((d, d)); M[i, j] = M[j, i] = 1 / np.sqrt(2); B.append(M)
    Q = []
    for i in range(d - 1):
        M = np.zeros((d, d)); M[i, i] = 1; M[d - 1, d - 1] = -1
        for N in Q: M = M - np.sum(M * N) * N
        Q.append(M / np.linalg.norm(M))
    return B + Q


for d, k in ((3, 2), (4, 2), (4, 3)):
    S = sym_basis(d); bound = (k + 1) / 2 - 1 / d; worst = 0
    for _ in range(3000):
        U = np.linalg.qr(rng.normal(size=(d, k)) + 1j * rng.normal(size=(d, k)))[0]
        V = np.linalg.qr(rng.normal(size=(d, k)) + 1j * rng.normal(size=(d, k)))[0]
        s = np.abs(rng.normal(size=k)); s /= np.linalg.norm(s)
        Psi = U @ np.diag(s) @ V.T
        T = np.array([[np.trace(Psi.conj().T @ Si @ Psi @ Sj).real for Sj in S] for Si in S])
        worst = max(worst, np.linalg.svd(T, compute_uv=False).sum())
    Psi = np.diag([1.0] * k + [0.0] * (d - k)) / np.sqrt(k)
    T = np.array([[np.trace(Psi.T @ Si @ Psi @ Sj) for Sj in S] for Si in S])
    eq = np.linalg.svd(T, compute_uv=False).sum()
    report(f"(92.1) d = {d}, k = {k}: largest nuclear norm over random Schmidt-rank-k vectors {worst:.4f} <= {bound:.4f};"
           f" real maximally entangled vector gives {eq:.6f}", worst <= bound + 1e-9 and abs(eq - bound) < 1e-12)
P = lin(a=3 / 7, b=-1 / 7, c=-1 / 7)
Pi_s0 = lambda Tp: (Tp + np.swapaxes(Tp, -1, -2)) / 2 - np.trace(Tp, axis1=-2, axis2=-1)[..., None, None] * I3 / 3
Phi = lambda Tp: np.trace(Tp, axis1=-2, axis2=-1)[..., None, None] * I3 / 3 - (2 / 7) * Pi_s0(Tp)
Es = lin(a=1 / 5, b=1 / 5, c=1 / 5); E0 = lin(a=1 / 3)
absPhi = sqrtm(superop(Phi).conj().T @ superop(Phi))
report("(92.2) P = Phi_{2/7,-id} and |Phi| = (5/7) E_s + (2/7) D/3",
       np.allclose(superop(Phi), superop(P)) and np.allclose(absPhi, 5 / 7 * superop(Es) + 2 / 7 * superop(E0)))
lhs = 8 * (np.eye(9) + F9 + Om)
rhs = 16 * sum(np.outer(np.kron(I3[i], I3[i]), np.kron(I3[i], I3[i])) for i in range(3))
rhs = rhs + sum(np.outer(np.kron(e, e), np.kron(e, e)) for e in itertools.product((1, -1), repeat=3) for e in [np.array(e, float)])
report("(92.2) separable identity 2^d (I + F + Omega) = 2^{d+1} sum |ii><ii| + sum_eps |eps eps><eps eps| (d = 3)", np.allclose(lhs, rhs))
CP = choi(P); vals = []
for _ in range(4000):
    U = np.linalg.qr(rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2)))[0]
    V = np.linalg.qr(rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2)))[0]
    s = np.abs(rng.normal(size=2)); s /= np.linalg.norm(s)
    psi = (U @ np.diag(s) @ V.T).reshape(9)
    vals.append(np.vdot(psi, CP @ psi).real)
wbell = (np.kron(I3[0], I3[0]) + np.kron(I3[1], I3[1])) / np.sqrt(2)
report(f"(92.2) P two-positive (min over random Schmidt-rank-two vectors {min(vals):.4f}), witness value at the ball edge",
       min(vals) > -1e-12 and abs(wbell @ CP.real @ wbell) < 1e-12, f"{wbell @ CP.real @ wbell:.1e}")
report("(92.3) integer cutoff 18*5^16 < 6^16 (and 18*5^15 > 6^15)", 18 * 5 ** 16 < 6 ** 16 and 18 * 5 ** 15 > 6 ** 15)
roots = np.roots([3, 8, 40, 68, 29]); roots = roots[np.abs(roots.imag) < 1e-9].real; capok = True
for t in roots:
    c = 3 * t * t + 4 * t + 8
    for _ in range(200):
        v = rng.normal(size=3); v /= np.linalg.norm(v)
        Bv = [np.outer(I3[i], v) + np.outer(v, I3[i]) + t * v[i] * I3 for i in range(3)]
        Kv = sum(B @ B for B in Bv) / c
        capok &= abs(np.trace(Kv) - 1) < 1e-12 and np.trace(Kv @ Kv) <= 5 / 6 + 1e-12
report("(92.3) one-copy purity cap Tr K(vv^T)^2 <= 5/6 for the audit-89 pencils at both roots (random v)", capok)
L3 = lin(a=3, b=-1, c=-1); ghz_ok = True
for n in (2, 3, 4):
    v = np.zeros(3 ** n)
    for i in range(3): v[sum(i * 3 ** s for s in range(n))] = 1
    X = site_apply(np.outer(v, v), n, L3).real
    words = list(itertools.product(range(3), repeat=n)); Y = np.zeros_like(X)
    for a, w1 in enumerate(words):
        Y[a, a] = sum(3 ** (n - w1.count(i)) for i in range(3))
        for i, j in itertools.combinations(range(3), 2):
            if set(w1) <= {i, j}:
                w2 = tuple(j if x == i else i for x in w1)
                Y[a, words.index(w2)] = 2 * (-1) ** n
    ghz_ok &= np.allclose(X, Y) and np.linalg.eigvalsh(X).min() > 0
report("(92.4) L_3^(x)n(vv^T), v = 0^n + 1^n + 2^n: block formula and positive definiteness, n = 2, 3, 4", ghz_ok)

# ============================================================================ 93
print("== round 93")


def Amats(p):
    A = [None] * 3
    for (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        M = np.zeros((3, 3)); M[j, k] = 1; M[k, j] = p; A[i] = M
    return A


def mps_vec(A, n):
    prods = [np.eye(A[0].shape[0])]
    for _ in range(n):
        prods = [Pm @ Ai for Pm in prods for Ai in A]
    return np.array([np.trace(Pm) for Pm in prods])


def kron4(a, b, c, d):
    return np.kron(np.kron(np.kron(a, b), c), d)


def transfers(p):
    A = Amats(p); B = [M.T for M in A]
    S = sum(kron4(A[i], B[i], A[j], B[j]) for i in range(3) for j in range(3))
    T = 2 * sum(kron4(A[i], B[j], A[i], B[j]) for i in range(3) for j in range(3)) - S
    U = 2 * sum(kron4(A[i], B[j], A[j], B[i]) for i in range(3) for j in range(3)) - S
    return T, U


idok = True
for p in (2.0, -2.0, 0.7):
    T, U = transfers(p); A = Amats(p); B = [M.T for M in A]
    for n in (2, 3, 4, 5):
        z, w = mps_vec(A, n), mps_vec(B, n)
        J = np.outer(z, w) - np.outer(w, z)
        lhs = 2 ** n * qform(J, n).real
        rhs = 2 * (np.trace(np.linalg.matrix_power(T, n)) - np.trace(np.linalg.matrix_power(U, n)))
        idok &= abs(lhs - rhs) <= 1e-8 * max(1, abs(lhs))
        if p == 2.0 and n in (3, 5): print(f"   p = 2, n = {n}: 2^n q_n(J_n) = {lhs:.1f}")
report("(93.1) 2^n q_n(J_n) = 2[Tr T^n - Tr U^n] for n = 2..5, p = 2, -2, 0.7", idok)
T, U = transfers(2.0)
for name, M in (("T", T), ("U", U)):
    G = csr_matrix((np.abs(M) > 1e-12).astype(int))
    ncomp, lab = connected_components(G, directed=True, connection='strong')
    sizes = sorted(np.bincount(lab).tolist(), reverse=True)
    rows, cols = G.nonzero(); periods = {}
    for c_ in range(ncomp):
        nodes = np.where(lab == c_)[0]
        if len(nodes) == 1 and G[nodes[0], nodes[0]] == 0: continue
        level = {nodes[0]: 0}; queue = [nodes[0]]; g = 0
        for u in queue:
            for v in G[u].indices:
                if lab[v] != c_: continue
                if v not in level: level[v] = level[u] + 1; queue.append(v)
                else: g = math.gcd(g, level[u] + 1 - level[v])
        periods[len(nodes)] = periods.get(len(nodes), []) + [g]
    odd_only_15 = all((g % 2 == 1) == (sz == 15) for sz, gs in periods.items() for g in gs)
    report(f"(93.2) {name}(2): strongly connected component sizes {[s for s in sizes if s > 1]},"
           f" only the 15-dimensional one has odd cycles", odd_only_15 and sorted([s for s in sizes if s > 1]) == sorted([15, 16, 16, 16, 2, 4, 4, 2, 4, 2]))
one = np.ones((3, 3)); B6 = np.block([[one - I3, -np.sqrt(2) * I3], [-np.sqrt(2) * I3, 2 * I3 - one]])
Gk = [I3, one - I3]
for _ in range(60): Gk.append(Gk[-1] + 4 * Gk[-2])
report("(93.2) B_6^2 - B_6 - 4I = 0 and G_k (G_{k+2} = G_{k+1} + 4G_k) entrywise nonnegative, k <= 61",
       np.allclose(B6 @ B6 - B6 - 4 * np.eye(6), 0) and all((G >= 0).all() for G in Gk))
odd_min, even_min = np.inf, np.inf
for p in (0.2, 0.5, 0.8, 1.3, 2.0, 3.0, -0.6, -1.7):
    T, U = transfers(p); ev_T = np.linalg.eigvals(T); ev_U = np.linalg.eigvals(U)
    scale = max(np.abs(ev_T).max(), np.abs(ev_U).max())
    for n in range(1, 42):
        val = (np.sum((ev_T / scale) ** n) - np.sum((ev_U / scale) ** n)).real
        if n % 2: odd_min = min(odd_min, val)
        else: even_min = min(even_min, val)
report(f"(93.3) Tr T^n - Tr U^n (scaled by rho^n) over n <= 41, eight p: smallest odd value {odd_min:.2e}",
       odd_min > -1e-9, f"(smallest even value {even_min:.2e}, even n are not covered by the round-93 proof)")
idok = True
for n in (2, 3, 4):
    Ts = [[rng.normal(size=(2, 2)) for _ in range(3)] for _ in range(4)]
    xs = [mps_vec(M, n) for M in Ts]
    C0, C1 = np.outer(xs[0], xs[1]), np.outer(xs[2], xs[3])
    A_, B_, C_, D_ = Ts
    Tm = 2 * sum(kron4(A_[i], B_[j], C_[i], D_[j]) for i in range(3) for j in range(3)) \
        - sum(kron4(A_[i], B_[i], C_[j], D_[j]) for i in range(3) for j in range(3))
    idok &= abs(2 ** n * qform(C0, n, C2=C1).real - np.trace(np.linalg.matrix_power(Tm, n))) < 1e-8 * (1 + abs(np.trace(np.linalg.matrix_power(Tm, n))))
report("(93.4) four-triple transfer 2^n q_n(C_0, C_1) = Tr T_01^n (random bond-two triples, n = 2, 3, 4)", idok)

# ============================================================================ 94
print("== round 94")
P7 = lin(a=3, b=-1, c=-1)
Q = np.zeros((9, 9))
for word in ((1, 1), (1, 2), (2, 1), (2, 2)): Q[3 * word[0] + word[1], 3 * word[0] + word[1]] = 1
sv = np.zeros(9); sv[1] = sv[3] = 1 / np.sqrt(2); Q += np.outer(sv, sv)
ok941 = True; marg = True
for m in (0, 1, 5, 17, 18, 20, 100):
    x = np.zeros(9); x[0] = m; x[4] = 1; x[8] = 1
    Z = site_apply(np.outer(x, x), 2, P7).real          # = 49 (m^2 + 2) P^(x)2(rho_m)
    ok941 &= abs(np.trace(Q @ Z) - (39 * m * m + 2 * m + 44)) < 1e-8 and abs(np.trace(Z) - 49 * (m * m + 2)) < 1e-8
    marg &= np.linalg.eigvalsh(Z).min() > 0
spec = sorted(np.round(np.linalg.eigvalsh(site_apply(np.outer(np.kron(I3[0], I3[0]), np.kron(I3[0], I3[0])), 2, P7).real), 8), reverse=True)
report("(94.1) Tr[Q 49(m^2+2) P^(x)2(rho_m)] = 39m^2 + 2m + 44, trace 49(m^2+2); product spectrum top five sum to 39",
       ok941 and abs(sum(spec[:5]) - 39) < 1e-8, f"(m = 20: {39*400+40+44} vs {39*402} over {49*402})")
report("(94.2) the outputs P^(x)2(rho_m) are positive definite (m = 0, 1, 5, 17, 18, 20, 100)", marg)
pstar = 1 / (0.5 + math.log(7 / math.sqrt(19)) / math.log(3))
report(f"(94.4) p_* = {pstar:.16f}", abs(pstar - 1.0739156609098173) < 1e-13)
bok = True
for _ in range(200):
    v = rng.normal(size=3) + 1j * rng.normal(size=3); v /= np.linalg.norm(v); t = abs(v @ v)
    ev = np.sort(np.linalg.eigvalsh(site_apply(np.outer(v, v.conj()), 1, lin(a=3 / 7, b=-1 / 7, c=-1 / 7))))
    bok &= np.allclose(ev, np.sort(np.array([3, 2 + t, 2 - t]) / 7))
report("(94.4) one-site output spectrum (3, 2 + t, 2 - t)/7 with t = |v^T v|", bok)
csok = True
for n in (1, 2, 3, 4):
    e0 = np.zeros(3 ** n); e0[0] = 1; e1 = np.zeros(3 ** n); e1[sum(3 ** s for s in range(n))] = 1
    C0, C1 = np.outer(e0, e0), np.outer(e1, e1)
    h = [[2 ** n * qform(a, n, C2=b).real for b in (C0, C1)] for a in (C0, C1)]
    crossed = 2 ** n * qform(np.outer(e0, e1), n).real
    csok &= np.allclose(h, [[1, (-1) ** n], [(-1) ** n, 1]]) and abs(crossed - 2 ** n) < 1e-9
M = np.full((3, 3), -3 / 5); np.fill_diagonal(M, 1)
report("(94.5) product planes: h = [[1, (-1)^n], [(-1)^n, 1]], crossed diagonal 2^n; signed vs unsigned compression",
       csok and np.allclose(np.linalg.eigvalsh(M), [-0.2, 1.6, 1.6]) and np.allclose(np.linalg.eigvalsh(np.abs(M)), [0.4, 0.4, 2.2]))

# ============================================================================ 95
print("== round 95")
okc = True
for u in (1.0, 2.0, 3.3, 4.0):
    Fu = lin(a=6 / 14, b=(u - 2) / 14, c=-(u + 2) / 14)
    ev = np.sort(np.linalg.eigvalsh(choi(Fu)))
    okc &= np.allclose(ev, np.sort([2 * (u - 1) / 14] + [(4 - u) / 14] * 5 + [(8 + u) / 14] * 3))
    X = rng.normal(size=(1, 1, 3, 3)) + 1j * rng.normal(size=(1, 1, 3, 3)); Y = Fu(X)
    okc &= np.allclose((Y + np.swapaxes(Y, -1, -2)) / 2, lin(a=3 / 7, b=-1 / 7, c=-1 / 7)(X))
report("(95.3) F_u = (6D + (u-2)id - (u+2)T)/14: Choi eigenvalues (2(u-1), 4-u, 8+u)/14 and J F_u = P", okc)
plus = np.array([1, 1j, 0]) / np.sqrt(2); minus = np.array([1, -1j, 0]) / np.sqrt(2)


def ghz_pm(n):
    a, b = np.array([1.0 + 0j]), np.array([1.0 + 0j])
    for _ in range(n): a, b = np.kron(a, plus), np.kron(b, minus)
    g = (a + b) / np.sqrt(2); return np.outer(g, g.conj())


lamok = True
for m in (1, 2, 3):
    G = ghz_pm(2 * m)
    for u in (1.0, 2.0, 3.5):
        Fu = lin(a=6 / 14, b=(u - 2) / 14, c=-(u + 2) / 14)
        ev = np.linalg.eigvalsh(ptrans(site_apply(G, 2 * m, Fu), 2 * m, range(m)))
        lam = ((16 - u * u) ** m - 2 ** (4 * m - 1)) / 14 ** (2 * m)
        lamok &= np.abs(ev - lam).min() < 1e-12
report("(95.4) lambda_{u,m} = ((16-u^2)^m - 2^{4m-1})/14^{2m} is an eigenvalue of the balanced partial transpose, m = 1, 2, 3", lamok)
G = ghz_pm(6); F2 = lin(a=6 / 14, c=-4 / 14); Y = site_apply(G, 6, F2); tn = 0; negs = {}; least = 0
for k in range(7):
    for S in itertools.combinations(range(6), k):
        ev = np.linalg.eigvalsh(ptrans(Y, 6, S)); tn += np.abs(ev).sum()
        if ev.min() < -1e-13: negs[k] = negs.get(k, 0) + 1
        if k == 3: least = min(least, ev.min())
tn /= 64
report(f"(95.4) u = 2, six copies: least balanced eigenvalue {least:.6e} vs -5/117649 = {-5/117649:.6e},"
       f" negative cuts by size {negs}, average trace norm - 1 = {tn - 1:.6e} vs 25/941192 = {25/941192:.6e}",
       abs(least + 5 / 117649) < 1e-12 and negs == {3: 20} and abs(tn - 1 - 25 / 941192) < 1e-10)
eps = Fr(2 ** 43 - 15 ** 11, 14 ** 22); delta = Fr(math.comb(22, 11), 2 ** 21) * eps; K = 10454529937283015932099844
report("(95.4/95.5) 2*15^11 < 16^11, delta = C(22,11) eps / 2^21 exact, K delta > 3^22 - 1",
       2 * 15 ** 11 < 16 ** 11 and delta == Fr(1843409240207301, 614127672205764052568277778432) and K * delta > 3 ** 22 - 1)
ghzpos = True
for n in (2, 3, 4, 5):
    ev = np.linalg.eigvalsh(site_apply(ghz_pm(n), n, lin(a=3 / 7, b=-1 / 7, c=-1 / 7)))
    ghzpos &= abs(ev.min() - 0.5 * (2 / 7) ** n) < 1e-13
report("(95.5) P^(x)n(g_n g_n^*) has least eigenvalue (1/2)(2/7)^n, n = 2..5", ghzpos)
F56 = lin(a=5 / 12, b=-1 / 12, c=-2 / 12); TF56 = lin(a=5 / 12, b=-2 / 12, c=-1 / 12)
spF = np.round(12 * np.linalg.eigvalsh(choi(F56)), 9); spT = np.round(12 * np.linalg.eigvalsh(choi(TF56)), 9)
X = rng.normal(size=(1, 1, 3, 3)) + 1j * rng.normal(size=(1, 1, 3, 3))
Pa = lambda a: lin(a=a, b=-0.25, c=-0.25)
avg_ok = np.allclose((F56(X) + TF56(X)) / 2, Pa(5 / 6)(X) / 2)
kraus = 0
for i, j in itertools.combinations(range(3), 2):
    Sij = np.outer(I3[i], I3[j]) + np.outer(I3[j], I3[i]); Dij = np.outer(I3[i], I3[i]) - np.outer(I3[j], I3[j])
    Aij = np.outer(I3[i], I3[j]) - np.outer(I3[j], I3[i])
    kraus = kraus + (3 * Sij @ X[0, 0] @ Sij.T + 2 * Dij @ X[0, 0] @ Dij.T + 7 * Aij @ X[0, 0] @ Aij.T) / 24
Ec = lin(a=1 / 4, b=1 / 4); Sf = superop(F56)
absF = sqrtm(Sf.conj().T @ Sf)
dec_ok = (np.allclose(absF, superop(lin(a=3 / 12, b=2 / 12, c=1 / 12)))
          and np.allclose(absF, 5 / 12 * superop(Es) + 1 / 3 * superop(Ec) + 1 / 4 * superop(E0))
          and np.allclose(Sf @ Sf, superop(lin(a=45 / 144, b=5 / 144, c=4 / 144)))
          and np.allclose(Sf @ Sf, 5 / 36 * superop(Es) + 1 / 36 * superop(Ec) + 5 / 6 * superop(E0))
          and np.linalg.eigvalsh(choi(Ec)).min() > -1e-12)
report("(95.7) a = 5/6: Choi spectra x12 of F and TF, (F + TF)/2 = P_{5/6}/2, Kraus form (95.8), |F| and F^2 decompositions",
       sorted(set(spF)) == [0, 3, 7] and sorted(set(spT)) == [-2, 4, 6] and avg_ok and np.allclose(kraus, F56(X)[0, 0]) and dec_ok,
       f"(F: {sorted(set(spF))}, TF: {sorted(set(spT))})")
worst_ppt = 0
for _ in range(3000):
    psi = rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6)); rho = psi @ psi.conj().T; rho /= np.trace(rho)
    R6 = rho.reshape(2, 3, 2, 3)
    out = np.zeros((2, 3, 2, 3), dtype=complex)
    for a in range(2):
        for b in range(2):
            out[a, :, b, :] = F56(R6[a, :, b, :].reshape(1, 1, 3, 3))[0, 0]
    worst_ppt = min(worst_ppt, np.linalg.eigvalsh(out.transpose(2, 1, 0, 3).reshape(6, 6)).min())
CT = choi(TF56); vals = []
for _ in range(4000):
    U = np.linalg.qr(rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2)))[0]
    V = np.linalg.qr(rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2)))[0]
    s = np.abs(rng.normal(size=2)); s /= np.linalg.norm(s)
    psi = (U @ np.diag(s) @ V.T).reshape(9); vals.append(np.vdot(psi, CT @ psi).real)
report(f"(95.7) F is 2-entanglement breaking (least eigenvalue of (id_2 x F)(rho)^Gamma over random rho {worst_ppt:.2e});"
       f" TF two-positive (least Schmidt-rank-two value {min(vals):.2e}, zero at (|00> + |11>)/sqrt 2)",
       worst_ppt > -1e-12 and min(vals) > -1e-12 and abs(wbell @ CT.real @ wbell) < 1e-12)
rng_ok = True
for a in (0.8, 0.9, 1.0):
    for t in np.linspace(0, 1.6, 33):
        Fat = lin(a=a, b=-0.25 + t / 2, c=-0.25 - t / 2)
        cp = np.linalg.eigvalsh(choi(Fat)).min() > -1e-12
        rng_ok &= cp == (1 - a - 1e-12 <= t <= 2 * a - 0.5 + 1e-12)
report("(95.6) CP range of F_{a,t} (unnormalised) is exactly 1 - a <= t <= 2a - 1/2 (a = 0.8, 0.9, 1)", rng_ok)

# ============================================================================ 96
print("== round 96")
al = 1001 / 1000
Ra = lin(a=-1, b=al, z=5 / 3)
WaG = al * np.eye(9) + 5 / 3 * np.diag([1 if i % 4 == 0 else 0 for i in range(9)]) - F9
s_ = al + 5 / 3; A_ = 1 / (s_ * (3 - s_)); C_ = 1000 / 1001; B_ = C_ - 3000 / 8003
Psi_a = lin(a=A_, b=-C_, z=B_)
X = rng.normal(size=(1, 1, 3, 3)) + 1j * rng.normal(size=(1, 1, 3, 3))
PsiG = A_ * np.eye(9) + B_ * np.diag([1 if i % 4 == 0 else 0 for i in range(9)]) - C_ * F9
report("(96.2-96.5) spectra of W_alpha^Gamma and of the Choi partial transpose of Psi_alpha; R_alpha Psi_alpha = -id",
       np.allclose(sorted(set(np.round(np.linalg.eigvalsh(WaG), 9))), sorted({5003 / 3000, 1 / 1000, 2001 / 1000}))
       and np.allclose(Ra(Psi_a(X)), -X)
       and np.allclose(sorted(set(np.round(np.linalg.eigvalsh(PsiG), 9))),
                       sorted({round(6009000 / 7978991, 9), round(1030009000 / 7986969991, 9), round(16987991000 / 7986969991, 9)})))
x = np.zeros(27); y = np.zeros(27)
for wd, sg in (('011', 1), ('102', -1), ('220', -1)): x[int(wd, 3)] = sg / np.sqrt(3)
for wd, sg in (('000', 1), ('121', 1), ('212', -1)): y[int(wd, 3)] = sg / np.sqrt(3)
Lm, Hm = np.outer(y, y) + np.outer(x, x), np.outer(y, y) - np.outer(x, x)
qL, qH = qform(Lm, 3, Ra).real, qform(Hm, 3, Ra).real
report(f"(96.6/96.7) q(L) = {qL:.10f}, q(H) = {qH:.10f} under R_alpha^(x)3; Werner value q(H) = {qform(Hm, 3).real:.10f}",
       abs(qL - 260468297081 / 40500000000) < 1e-12 and abs(qH + 9693702919 / 40500000000) < 1e-12 and abs(qform(Hm, 3).real - 1) < 1e-12)
f1 = np.ones(3) / np.sqrt(3); g1 = np.array([1, -1, 0]) / np.sqrt(2)
b1 = lambda A, B: np.vdot(A, Ra(B.reshape(1, 1, 3, 3).astype(complex))[0, 0]).real
tab = [b1(np.outer(f1, f1), np.outer(f1, f1)), b1(np.outer(g1, g1), np.outer(g1, g1)), b1(np.outer(f1, f1), np.outer(g1, g1)),
       b1(np.outer(f1, g1), np.outer(f1, g1)), b1(np.outer(f1, g1), np.outer(g1, f1))]
mixed = max(abs(b1(Aa, Bb)) for Aa in (np.outer(f1, f1), np.outer(g1, g1)) for Bb in (np.outer(f1, g1), np.outer(g1, f1)))
report("(96.9) flag table 9000 x (a, e, b, c, d) and vanishing mixed pairings",
       np.allclose(np.array(tab) * 9000, [5009, 7509, -4000, 14009, 5000]) and mixed < 1e-14, f"{np.round(np.array(tab) * 9000, 6).tolist()}")
liftok = True
for r in (1, 2):
    Fr_, Gr_ = f1, g1
    for _ in range(r - 1): Fr_, Gr_ = np.kron(Fr_, f1), np.kron(Gr_, g1)
    hp, hm = (Fr_ + Gr_) / np.sqrt(2), (Fr_ - Gr_) / np.sqrt(2)
    Pi = np.kron(np.outer(y, y), np.outer(hp, hp)) + np.kron(np.outer(x, x), np.outer(hm, hm))
    a_, e_, b_, c_, d_ = tab
    formula = 0.25 * ((a_ ** r + e_ ** r + 2 * b_ ** r) * qL + 2 * (c_ ** r + d_ ** r) * qH)
    liftok &= abs(qform(Pi, 3 + r, Ra).real - formula) < 1e-12 and np.allclose(Pi @ Pi, Pi) and abs(np.trace(Pi) - 2) < 1e-12
tabF = [Fr(5009, 9000), Fr(7509, 9000), Fr(-4000, 9000), Fr(14009, 9000), Fr(5000, 9000)]
qLF, qHF = Fr(260468297081, 40500000000), Fr(-9693702919, 40500000000)
a_, e_, b_, c_, d_ = tabF
q9 = Fr(1, 4) * ((a_ ** 6 + e_ ** 6 + 2 * b_ ** 6) * qLF + 2 * (c_ ** 6 + d_ ** 6) * qHF)
report("(96.10/96.11) flag lift: identity at r = 1, 2 (direct, 4 and 5 copies), exact nine-copy value",
       liftok and q9 == Fr(-23476199337924225624801226904672029, 21523360500000000000000000000000000), f"{float(q9):.11f}")
R1 = lin(a=-1, b=1, z=5 / 3)
b1 = lambda A, B: np.vdot(A, R1(B.reshape(1, 1, 3, 3).astype(complex))[0, 0]).real
tab1 = [Fr(b1(A, B)).limit_denominator(1000) for A, B in ((np.outer(f1, f1), np.outer(f1, f1)), (np.outer(g1, g1), np.outer(g1, g1)),
        (np.outer(f1, f1), np.outer(g1, g1)), (np.outer(f1, g1), np.outer(f1, g1)), (np.outer(f1, g1), np.outer(g1, f1)))]
qL1, qH1 = qform(Lm, 3, R1).real, qform(Hm, 3, R1).real
a_, e_, b_, c_, d_ = tab1
q91 = Fr(1, 4) * ((a_ ** 6 + e_ ** 6 + 2 * b_ ** 6) * Fr(520, 81) + 2 * (c_ ** 6 + d_ ** 6) * Fr(-20, 81))
report("(96.6) alpha = 1: q(L) = 520/81, q(H) = -20/81, nine-copy value -1574982175/1377495072",
       abs(qL1 - 520 / 81) < 1e-12 and abs(qH1 + 20 / 81) < 1e-12 and q91 == Fr(-1574982175, 1377495072))

# ============================================================================ 97
print("== round 97")


def gamma_sum(n, k):
    return sum(Fr(-1, 2) ** j * Fr(math.comb(n - k, j) * math.comb(n + k + 2, j), math.comb(n, j)) for j in range(n - k + 1))


def gamma_closed(n, k):
    m = n - k
    if m % 2 == 0:
        l = m // 2
        return Fr(math.factorial(2 * l) * math.factorial(n - l + 1), 4 ** l * math.factorial(l) * math.factorial(n) * (n - 2 * l + 1))
    l = (m - 1) // 2
    return -Fr(math.factorial(2 * l + 1) * math.factorial(n - l), 2 * 4 ** l * math.factorial(l) * math.factorial(n) * (n - 2 * l))


report("(97.4-97.7) closed forms equal the finite sums, every k, n <= 60",
       all(gamma_sum(n, k) == gamma_closed(n, k) for n in range(1, 61) for k in range(n + 1)))
report("(97.9) min_k gamma_{n,k} = -1/(2n) for 4 <= n <= 250, and -1/2, -1/4, -1/4 at n = 1, 2, 3",
       all(min(gamma_closed(n, k) for k in range(n + 1)) == Fr(-1, 2 * n) for n in range(4, 251))
       and [min(gamma_closed(n, k) for k in range(n + 1)) for n in (1, 2, 3)] == [Fr(-1, 2), Fr(-1, 4), Fr(-1, 4)])


def hook_dim(lam):
    n = sum(lam); cols = [sum(1 for r in lam if r > c) for c in range(lam[0])]; prod = 1
    for i, r in enumerate(lam):
        for c in range(r): prod *= (r - c - 1) + (cols[c] - i - 1) + 1
    return math.factorial(n) // prod


hook_ok = all(hook_dim(lam) >= n - 1 for n in range(5, 41) for a in range(n, 0, -1) for b in range(min(a, n - a), -1, -1)
              for c in [n - a - b] if 0 <= c <= b for lam in [tuple(v for v in (a, b, c) if v > 0)] if lam != (n,))
s_t, a_t, b_t = rng.uniform(1, 5, 50), rng.uniform(0, 3, 50), rng.uniform(0, 3, 50)
ident = np.allclose(s_t ** 2 * (2 * s_t - 1) * ((a_t + b_t) ** 2 / (2 * s_t - 1) - 2 * a_t * b_t / s_t - b_t ** 2 / s_t ** 2),
                    (s_t * a_t - (s_t - 1) * b_t) ** 2)
report("(97.11/97.12) f_lambda >= n - 1 for nontrivial lambda with <= 3 rows (5 <= n <= 40); the polynomial square identity", hook_ok and ident)
P0 = np.outer(I3[0], I3[0]); P1 = np.outer(I3[1], I3[1]); P2 = np.outer(I3[2], I3[2])
Pex = np.kron(P1 + P2, P0); Eex = twirl(Pex, 2)
ex2 = [qform(Pex, 2).real, qform(Eex, 2).real, qform(Pex - Eex, 2).real]
w3 = lambda s: np.outer(*(2 * [np.eye(27)[int(s, 3)]]))
Pex3 = w3('012') + w3('120')
cyc = [(0, 1, 2), (1, 2, 0), (2, 0, 1)]
ex3 = [qform(Pex3, 3).real, qform(twirl(Pex3, 3), 3).real, qform(twirl(Pex3, 3, cyc), 3).real]
report(f"(97.7) two-copy example {np.round(ex2, 12).tolist()} (expect 0, 1/2, -1/2); three-copy example {np.round(ex3, 12).tolist()}"
       " (expect 0, 1/6, -1/6)", np.allclose(ex2, [0, 0.5, -0.5]) and np.allclose(ex3, [0, 1 / 6, -1 / 6]))
V97 = np.array([[28,-633],[3,-590],[-14,-1411],[3,-590],[-48,-244],[-163,-297],[-14,-1411],[-163,-297],[-533,264],[3,-590],[-48,-244],[-163,-297],[-48,-244],[-29,-81],[-55,-56],[-163,-297],[-55,-56],[-39,80],[-14,-1411],[-163,-297],[-533,264],[-163,-297],[-55,-56],[-39,80],[-533,264],[-39,80],[344,31],[-1715,1277],[-1561,1156],[-3695,2730],[-1561,1156],[-566,405],[-522,340],[-3695,2730],[-522,340],[1550,-1302],[-1561,1156],[-566,405],[-522,340],[-566,405],[-166,114],[-59,28],[-522,340],[-59,28],[272,-213],[-3695,2730],[-522,340],[1550,-1302],[-522,340],[-59,28],[272,-213],[1550,-1302],[272,-213],[-471,448],[-1324,-1640],[-1197,-1485],[-2824,-3508],[-1197,-1485],[-416,-522],[-340,-441],[-2824,-3508],[-340,-441],[1384,1655],[-1197,-1485],[-416,-522],[-340,-441],[-416,-522],[-117,-148],[-25,-37],[-340,-441],[-25,-37],[223,273],[-2824,-3508],[-340,-441],[1384,1655],[-340,-441],[-25,-37],[223,273],[1384,1655],[223,273],[-488,-565]], float)
P4 = V97 @ np.linalg.inv(V97.T @ V97) @ V97.T; E4 = twirl(P4, 4)
q4, q4f = qform(P4, 4).real, qform(P4 - E4, 4).real
report(f"(97.17) four-copy rational projector: q_4(P) = {q4:.5f}, q_4(P - E_4 P) = {q4f:.11f} < -1/4",
       abs(q4f + 0.27193698951) < 1e-10 and abs(q4 - 0.03824) < 1e-4)
fam_ok = True
f3 = np.ones(3) / np.sqrt(3); Kp = np.eye(3) - np.outer(f3, f3)
for n in (3, 4, 5):
    g = np.zeros(3 ** (n - 1))
    for i in range(3): g[sum(i * 3 ** s for s in range(n - 1))] = 1 / np.sqrt(3)
    Pn = np.kron(Kp, np.outer(g, g)); En = twirl(Pn, n)
    pred = (n - 1) / (27 * n) * (2 + 2.0 ** (-(n - 2)) + 2 * (-0.5) ** (n - 2))
    fam_ok &= abs(qform(Pn, n).real) < 1e-12 and abs(qform(En, n).real - pred) < 1e-12 and abs(qform(Pn - En, n).real + pred) < 1e-12
report("(97.18-97.20) P_n = K x P_g: q_n(P_n) = 0, q_n(E_n P_n) = (n-1)/(27n)[2 + 2^{-(n-2)} + 2(-1/2)^{n-2}], n = 3, 4, 5", fam_ok)

# ============================================================================ 98
print("== round 98")
cal_ok = True
for m in (2, 3):
    tail = np.outer(I3[0], I3[1])
    for _ in range(m - 2): tail = np.kron(tail, np.outer(I3[0], I3[1]))
    A = [np.kron(np.outer(I3[i], I3[i]), tail) for i in range(3)]
    Z = sum(A)
    cal_ok &= all(abs(qform(A[i] + A[j], m).real) < 1e-12 and abs(qform(A[i] - A[j], m).real - 2) < 1e-12
                  for i, j in itertools.combinations(range(3), 2)) and abs(qform(Z, m).real + 1.5) < 1e-12
for m in (1, 3):
    Bm = []
    for i in range(3):
        e = np.zeros(3 ** m); e[sum(i * 3 ** s for s in range(m))] = 1; Bm.append(np.outer(e, e))
    cal_ok &= all(abs(qform(Bm[i], m).real - 2.0 ** -m) < 1e-12 for i in range(3))
    cal_ok &= all(abs(qform(Bm[i] + Bm[j], m).real) < 1e-12 for i, j in itertools.combinations(range(3), 2))
    cal_ok &= abs(qform(sum(Bm), m).real + 3 * 2.0 ** -m) < 1e-12
Gm = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); Gm = Gm + Gm.conj().T
fq = lambda c: np.vdot(c, Gm @ c).real; e3 = np.eye(3)
pol = abs(fq(np.ones(3)) - (0.75 * sum(fq(e3[i] + e3[j]) for i, j in itertools.combinations(range(3), 2))
                            - 0.25 * sum(fq(e3[i] - e3[j]) for i, j in itertools.combinations(range(3), 2)))) < 1e-12
pol &= abs(fq(np.ones(3)) - (sum(fq(e3[i] + e3[j]) for i, j in itertools.combinations(range(3), 2)) - sum(fq(e3[i]) for i in range(3)))) < 1e-12
report("(98.3/98.4/98.8) rank-two calibrations (m = 2, 3 and odd m = 1, 3) and both polarisation identities", cal_ok and pol)


def occ_basis(N):
    return [m for m in itertools.product(range(N + 1), repeat=3) if sum(m) == N]


def annihilators(N):
    top, bot = occ_basis(N), occ_basis(N - 1); idx = {m: i for i, m in enumerate(bot)}
    ops = []
    for i in range(3):
        a = np.zeros((len(bot), len(top)))
        for c, m in enumerate(top):
            if m[i] > 0:
                mm = list(m); mm[i] -= 1; a[idx[tuple(mm)], c] = np.sqrt(m[i])
        ops.append(a)
    return ops


def q_bos(C, N, t):
    tot, X = 0.0, C
    for j in range(N + 1):
        tot += (-t) ** j * math.factorial(N - j) / (math.factorial(j) * math.factorial(N)) * np.vdot(X, X).real
        if j < N:
            ops = annihilators(N - j); X = sum(a @ X @ a.T for a in ops)
    return tot


# validate the bosonic formula against direct partial traces at N = 3
N = 3; basis = occ_basis(N)
emb = np.zeros((27, len(basis)))
for c, m in enumerate(basis):
    words = [w for w in itertools.product(range(3), repeat=N) if tuple(w.count(i) for i in range(3)) == m]
    for w in words: emb[int(''.join(map(str, w)), 3), c] = 1 / np.sqrt(len(words))
Cb = rng.normal(size=(len(basis), 2)) @ rng.normal(size=(2, len(basis)))
bos_ok = abs(q_bos(Cb, N, 0.5) - qform(emb @ Cb @ emb.T, N).real) < 1e-10
spec_ok = True
for N in range(2, 8):
    basis = occ_basis(N); idx = {m: i for i, m in enumerate(basis)}
    Vn = np.zeros((len(basis), 2)); Vn[idx[(1, 0, N - 1)], 0] = 1; Vn[idx[(0, 1, N - 1)], 1] = 1
    for t in (0.3, 0.5, 0.6):
        C2 = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        aN = (1 - t) ** (N - 2) * (1 - t / N); bN = t * (1 - t) ** (N - 2) * (t - 1 / N)
        spec_ok &= abs(q_bos(Vn @ C2 @ Vn.T, N, t) - (aN * np.vdot(C2, C2).real + bN * abs(np.trace(C2)) ** 2)) < 1e-10
basis = occ_basis(2); idx = {m: i for i, m in enumerate(basis)}
Vn = np.zeros((len(basis), 2)); Vn[idx[(1, 0, 1)], 0] = 1; Vn[idx[(0, 1, 1)], 1] = 1
diag = q_bos(Vn @ np.eye(2) @ Vn.T, 2, 0.6)
report(f"(98.11) bosonic formula = partial traces (N = 3); spectator coefficients a_N, b_N (N = 2..7); I_2 at t = 3/5, N = 2: {diag:.6f}",
       bos_ok and spec_ok and abs(diag - 41 / 25) < 1e-12)

print(f"\n{sum(RESULTS)} of {len(RESULTS)} checks passed")
