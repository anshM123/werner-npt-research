"""MODULAR exact verification (all global tensors formed only modulo primes; CRT on the compressed blocks) of a level-2 DPS certificate for block-positivity of W = C_a^{(x)n}.

Self-contained: rebuilds every per-copy object from scratch (does not import brauer.py or blocks.py).

Statement checked.  Let A = A_1..A_n (inputs), B = B_1..B_n (outputs), B' = a copy of B, all copies C^3.
P = I_A (x) (I + SWAP_{B,B'})/2.  Given rational coefficient vectors of two operators Q_A, Q_B in the span of
products of Brauer diagrams (one diagram per copy, acting on A_s B_s B'_s), define
    Q_A  := P Q_A_raw P,
    Z    := P [ W (x) I_B'  -  t I  -  Gamma_A(Q_A)  -  Gamma_B(Q_B) ] P,
and check EXACTLY that Q_A >= 0, Q_B >= 0, Z >= 0.  Then for all x in C^{3^n}, y in C^{3^n}, with
psi = x (x) y (x) y (P psi = psi):
    <x y|W|x y> |y|^2 = <psi|W (x) I|psi> = t|psi|^2 + <psi|Z|psi> + <x~ y y|Q_A|x~ y y> + <x y~ y|Q_B|x y~ y> >= t |x|^2 |y|^4
(x~ = complex conjugate), i.e. <x (x) y|W|x (x) y> >= t |x|^2 |y|^2.

Positive semidefiniteness of an operator Q in the span of Brauer-diagram products is checked block by block:
Q commutes with SO(3)^n (one rotation per copy, acting on all three strands of that copy), so Q >= 0 iff for every
spin label (j_1..j_n) the Gram compression B_beta^T Q B_beta >= 0, where the columns of
B_beta = B_{j_1} (x) .. (x) B_{j_n} span the m = 0 subspace of that isotypic component (B_j exact integer matrices).

usage: python exact_verify.py cert.txt
"""
import sys
import math
import time
import itertools
from fractions import Fraction
import numpy as np
import sympy as sp

T0 = time.time()
d = 3
NS = 3
DIM = d ** NS


def log(*a):
    print(*a, flush=True)


# ------------------------------------------------------------------ read certificate
cert = {}
with open(sys.argv[1]) as f:
    for line in f:
        parts = line.split()
        cert[parts[0]] = parts[1:]
n = int(cert['n'][0])
a = Fraction(int(cert['a'][0]), int(cert['a'][1]))
t = Fraction(int(cert['t'][0]), int(cert['t'][1]))
L = int(cert['L'][0])
diagrams = []
for k in range(15):
    s = cert[f'D{k}']
    pts = list(map(int, s))
    diagrams.append(tuple((pts[2 * i], pts[2 * i + 1]) for i in range(3)))
if 'QA' in cert:
    QA_rawI = np.array([int(v) for v in cert['QA']], dtype=object).reshape([15] * n)     # = L * Q_A_raw
    QBI = np.array([int(v) for v in cert['QB']], dtype=object).reshape([15] * n)          # = L * Q_B
log(f'certificate: n = {n}, a = {a}, t = {t}, L = {L}')

# ------------------------------------------------------------------ per-copy objects from scratch
# basis |i_A i_B i_B'>, index 9 i_A + 3 i_B + i_B'; diagram point p<3: row index of strand p, p>=3: column index of strand p-3
assert len(set(frozenset(frozenset(p) for p in m) for m in diagrams)) == 15
for m in diagrams:
    assert sorted(x for p in m for x in p) == list(range(6))


def dmat(m):
    M = np.zeros((DIM, DIM), dtype=np.int64)
    for idx in itertools.product(range(d), repeat=6):
        if all(idx[p] == idx[q] for p, q in m):
            M[9 * idx[0] + 3 * idx[1] + idx[2], 9 * idx[3] + 3 * idx[4] + idx[5]] = 1
    return M


DM = [dmat(m) for m in diagrams]
rank = sp.Matrix(np.array([M.reshape(-1) for M in DM])).rank()
assert rank == 15, 'diagrams not independent'
log('15 diagram operators built, linearly independent')


def find(M):
    hits = [b for b in range(15) if np.array_equal(M, DM[b])]
    assert len(hits) == 1
    return hits[0]


def ptrans(M, strand):
    T = M.reshape([d] * 6)
    T = np.swapaxes(T, strand, strand + 3)
    return T.reshape(DIM, DIM)


piA = [find(ptrans(M, 0)) for M in DM]      # Gamma_A(D_a) = D_{piA[a]}
piB = [find(ptrans(M, 1)) for M in DM]      # Gamma_B(D_a) = D_{piB[a]}
# swap of B and B'
S = np.zeros((DIM, DIM), dtype=np.int64)
for i, j, k in itertools.product(range(d), repeat=3):
    S[9 * i + 3 * k + j, 9 * i + 3 * j + k] = 1
muL = [find(S @ M) for M in DM]             # S D_a = D_{muL[a]}
muR = [find(M @ S) for M in DM]             # D_a S = D_{muR[a]}
iS = find(S)
iI = find(np.eye(DIM, dtype=np.int64))
log('partial transposes on A and on B, and left/right multiplication by SWAP(B,B\') act as permutations of diagrams')

# orbit-compressed certificates (lines QAO / QBO): one integer per orbit of diagram tuples under permutations of the
# copies and the transpose of every diagram; orbits are recomputed here from scratch (the transpose permutation of the
# diagrams is derived from the explicit matrices), and the tensor is the orbit value broadcast to all members.
if 'QAO' in cert:
    trp = np.array([find(M.T.copy()) for M in DM])
    allk_ = np.array(list(itertools.product(range(15), repeat=n)), dtype=np.int64)
    pw_ = 15 ** np.arange(n - 1, -1, -1)
    canon_ = np.minimum(np.sort(allk_, axis=1) @ pw_, np.sort(trp[allk_], axis=1) @ pw_)
    _, orb_ = np.unique(canon_, return_inverse=True)
    vA = np.array([int(v) for v in cert['QAO']], dtype=object)
    vB = np.array([int(v) for v in cert['QBO']], dtype=object)
    assert len(vA) == len(vB) == orb_.max() + 1
    log('orbit-compressed certificate:', len(vA), 'orbit values per operator')

# C_a (x) I_B' expanded in diagrams (exact solve)
F = np.zeros((9, 9), dtype=object)
Om = np.zeros((9, 9), dtype=object)
for i, j in itertools.product(range(d), repeat=2):
    F[3 * j + i, 3 * i + j] = 1
for i, j in itertools.product(range(d), repeat=2):
    Om[4 * i, 4 * j] = 1
Ca = np.empty((9, 9), dtype=object)
for r in range(9):
    for c in range(9):
        Ca[r, c] = (a if r == c else Fraction(0)) - Fraction(F[r, c] + Om[r, c], 4)
CaI = np.empty((DIM, DIM), dtype=object)
for r in range(DIM):
    for c in range(DIM):
        rA, rB, rB2 = r // 9, (r // 3) % 3, r % 3
        cA, cB, cB2 = c // 9, (c // 3) % 3, c % 3
        CaI[r, c] = Ca[3 * rA + rB, 3 * cA + cB] if rB2 == cB2 else Fraction(0)
Amat = sp.Matrix(np.array([M.reshape(-1) for M in DM]).T.tolist())
sol = Amat.solve_least_squares(sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in CaI.reshape(-1)]))
cvec = [Fraction(int(sp.fraction(v)[0]), int(sp.fraction(v)[1])) for v in sol]
chk = sum((np.array(DM[b], dtype=object) * cvec[b] for b in range(15)), np.zeros((DIM, DIM), dtype=object))
assert all(chk.reshape(-1)[k] == CaI.reshape(-1)[k] for k in range(DIM * DIM)), 'C_a (x) I not in diagram span'
log('C_a (x) I_B\' in diagram coordinates:', {b: str(cvec[b]) for b in range(15) if cvec[b] != 0})

# spin m = 0 bases, exact
eps = np.zeros((3, 3, 3), dtype=np.int64)
for p in itertools.permutations(range(3)):
    eps[p] = int(sp.combinatorics.Permutation(list(p)).signature())
Lg = [-eps[k] for k in range(3)]
I3 = np.eye(3, dtype=np.int64)


def tot(Lk):
    return (np.kron(np.kron(Lk, I3), I3) + np.kron(np.kron(I3, Lk), I3) + np.kron(np.kron(I3, I3), Lk))


Jt = [sp.Matrix(tot(Lk)) for Lk in Lg]
Jz = Jt[2]
Cas = Jt[0] ** 2 + Jt[1] ** 2 + Jt[2] ** 2
ker = Jz.nullspace()
assert len(ker) == 7
Bj = {}
tot_m = 0
for j in range(4):
    ns = (Cas + j * (j + 1) * sp.eye(DIM)).nullspace()
    # intersect with ker Jz: solve Jz * (N c) = 0
    Nm = sp.Matrix.hstack(*ns)
    cs = (Jz * Nm).nullspace()
    Bm = Nm * sp.Matrix.hstack(*cs)
    den = sp.ilcm(*[sp.fraction(x)[1] for x in Bm])
    Bm = Bm * den
    assert Jz * Bm == sp.zeros(DIM, Bm.shape[1])
    assert Cas * Bm == -j * (j + 1) * Bm
    assert Bm.rank() == Bm.shape[1]
    Bj[j] = np.array(Bm.tolist(), dtype=object)
    tot_m += Bm.shape[1]
    # dimension of the isotypic component must be (2j+1) * multiplicity
    assert len(ns) == (2 * j + 1) * Bm.shape[1]
assert tot_m == 7
mult = {j: Bj[j].shape[1] for j in range(4)}
log('spin m=0 bases: multiplicities', mult)
# every diagram commutes with the three total rotation generators (so with SO(3) acting as O(x)O(x)O)
for M in DM:
    Ms = sp.Matrix(M)
    for J in Jt:
        assert J * Ms == Ms * J
log('every diagram commutes with the total so(3) generators (exact)')
# per-copy compressions Psi_j(a) = B_j^T D_a B_j (integers)
Psi = {j: np.array([Bj[j].T.dot(np.array(DM[b], dtype=object)).dot(Bj[j]) for b in range(15)], dtype=object)
       for j in range(4)}
# the compression map (15 diagram coefficients -> the four blocks, 1+9+4+1 = 15 numbers) is injective, so an
# element of the span of diagram products is determined by its compressions (in particular it is symmetric
# as soon as all its compressions are symmetric)
PsiMat = sp.Matrix([[int(Psi[j][b][r, c]) for b in range(15)]
                    for j in range(4) for r in range(mult[j]) for c in range(mult[j])])
assert PsiMat.shape == (15, 15) and PsiMat.rank() == 15
log('per-copy compression map is a bijection (rank 15, exact)')

# float sanity check of the block criterion at one copy: spectrum of a random symmetric element of the algebra
rng = np.random.default_rng(1)
for _ in range(5):
    xr = rng.standard_normal(15)
    M = sum(xr[b] * DM[b] for b in range(15))
    M = (M + M.T) / 2
    ev = np.sort(np.linalg.eigvalsh(M))
    evb = []
    for j in range(4):
        Bf = np.array(Bj[j], dtype=float)
        G = Bf.T @ M @ Bf
        H = Bf.T @ Bf
        w = np.linalg.eigvals(np.linalg.solve(H, G)).real
        evb += list(w) * (2 * j + 1)
    assert np.allclose(np.sort(evb), ev, atol=1e-9)
log('float sanity check: spectrum of random symmetric algebra elements = union of block spectra (x(2j+1))')



# ------------------------------------------------------------------ S-adapted exact bases B_j' = [B_j^+ | B_j^-]
import flint
from sympy import prevprime
Sm = sp.Matrix(S)
Bp, sgn = {}, {}
for j in range(4):
    Bm = sp.Matrix(Bj[j].tolist())
    sig = (Bm.T * Bm).inv() * Bm.T * Sm * Bm
    assert Sm * Bm == Bm * sig                       # S preserves the m=0 spin-j space
    cols, sg = [], []
    for e in (1, -1):
        for v in (sig - e * sp.eye(sig.shape[0])).nullspace():
            w = Bm * v
            w = w * sp.ilcm(*[sp.fraction(x)[1] for x in w])
            assert Sm * w == e * w
            cols.append(w)
            sg.append(e)
    M_ = sp.Matrix.hstack(*cols)
    assert M_.shape[1] == Bm.shape[1] and M_.rank() == Bm.shape[1]
    Bp[j] = np.array(M_.tolist(), dtype=np.int64)
    sgn[j] = np.array(sg)
log('S-adapted spin bases: signs', {j: sgn[j].tolist() for j in range(4)})
PsiP = {j: np.array([Bp[j].T @ DM[b] @ Bp[j] for b in range(15)], dtype=np.int64) for j in range(4)}
PsiMat = sp.Matrix([[int(PsiP[j][b][r, c]) for b in range(15)]
                    for j in range(4) for r in range(mult[j]) for c in range(mult[j])])
assert PsiMat.rank() == 15
log('per-copy compression map in the S-adapted bases is a bijection (rank 15)')


# ======================================================================================================================
# MODULAR global algebra (memory-light): every global tensor is only ever formed modulo a prime p < 2^29.
# Exact integers enter only through (a) the certificate numerators (orbit values or full tensors), (b) explicit
# l1-bounds that fix how many primes are needed, and (c) the CRT reconstruction of the (small) compressed blocks.
# The operators checked are the positive multiples  L Q_B,  4L Q_A (Q_A = P Q_A_raw P),  4 D_Z Z  of the certificate.
# ======================================================================================================================
def perm_axes(X, pi):
    """Y[pi(a1),..,pi(an)] = X[a1,..,an]"""
    Y = np.empty_like(X)
    Y[np.ix_(*([np.array(pi)] * n))] = X
    return Y


def Pc4_mod(X, p_):
    """4 P X P = X + S X + X S + S X S  (mod p)"""
    SX = perm_axes(X, muL)
    return (X + SX + perm_axes(X, muR) + perm_axes(SX, muR)) % p_


Dc = 1
for c_ in cvec:
    Dc = math.lcm(Dc, c_.denominator)
cint = [int(c_ * Dc) for c_ in cvec]                       # W (x) I = Dc^-n * tpow(cint)
DZ = math.lcm(Dc ** n, t.denominator, 4 * L)
kW, kI, kA, kB = DZ // Dc ** n, t.numerator * (DZ // t.denominator), DZ // (4 * L), DZ // L
# inner = kW*tpow(cint) - kI*I - kA*Gamma_A(4L Q_A) - kB*Gamma_B(L Q_B);  4 D_Z Z = Pc4(inner)

if 'QAO' in cert:
    cnt_ = np.bincount(orb_)
    l1_A_raw = sum(abs(int(v)) * int(c) for v, c in zip(vA, cnt_))
    l1_B = sum(abs(int(v)) * int(c) for v, c in zip(vB, cnt_))
else:
    l1_A_raw = sum(abs(int(v)) for v in QA_rawI.reshape(-1))
    l1_B = sum(abs(int(v)) for v in QBI.reshape(-1))
l1_W = sum(abs(v) for v in cint) ** n
l1_A = 4 * l1_A_raw                                         # Pc4 at most quadruples the l1 norm
l1_Z = 4 * (kW * l1_W + abs(kI) + kA * l1_A + kB * l1_B)   # permutations preserve the l1 norm
maxPsi = max(int(np.abs(PsiP[j]).max()) for j in range(4))
glob_bound = max(l1_A, l1_B, l1_Z) * maxPsi ** n            # bound on every entry of every compressed block
primes = []
p = 2 ** 29
while math.prod(primes) <= 4 * glob_bound:
    p = prevprime(p)
    primes.append(p)
Mod = math.prod(primes)
assert 15 * max(primes) ** 2 < 2 ** 63
log('CRT with', len(primes), 'primes below 2^29 (entry bound', glob_bound.bit_length(), 'bits)')

beta_list = list(itertools.combinations_with_replacement(range(4), n))
log('blocks to check per operator (sorted spin labels; the certificate is copy-permutation invariant by'
    ' construction for orbit input, and is checked below otherwise):', len(beta_list))
res = {name: {beta: [] for beta in beta_list} for name in ('Q_B', 'Q_A', 'Z')}


def tpow_mod(v, p_):
    X = np.array([1], dtype=np.int64)
    for _ in range(n):
        X = (np.multiply.outer(X, np.array([x % p_ for x in v], dtype=np.int64)) % p_).reshape(-1)
    return X.reshape([15] * n)


def invariant(X):
    gens = [tuple([1, 0] + list(range(2, n)))] + ([tuple(list(range(1, n)) + [0])] if n > 2 else [])
    return all(np.array_equal(np.transpose(X, g), X) for g in gens)


if 'QAO' not in cert and not (invariant(QA_rawI) and invariant(QBI)):
    raise SystemExit('full-format certificate not invariant under permutations of the copies: use exact_verify_fast.py')
log('copy-permutation invariance: by construction (orbit format) or checked exactly on the integer tensors')
for ip, p_ in enumerate(primes):
    if 'QAO' in cert:
        A_raw = np.array([int(v) % p_ for v in vA], dtype=np.int64)[orb_].reshape([15] * n)
        Bm = np.array([int(v) % p_ for v in vB], dtype=np.int64)[orb_].reshape([15] * n)
    else:
        A_raw = np.array([int(v) % p_ for v in QA_rawI.reshape(-1)], dtype=np.int64).reshape([15] * n)
        Bm = np.array([int(v) % p_ for v in QBI.reshape(-1)], dtype=np.int64).reshape([15] * n)
    Am = Pc4_mod(A_raw, p_)                                  # 4L Q_A mod p
    inner = (tpow_mod(cint, p_) * (kW % p_)) % p_
    Iidx = tuple([iI] * n)
    inner[Iidx] = (int(inner[Iidx]) - kI) % p_
    inner = (inner - (perm_axes(Am, piA) * (kA % p_)) % p_ - (perm_axes(Bm, piB) * (kB % p_)) % p_) % p_
    Zm = Pc4_mod(inner, p_)                                  # 4 D_Z Z mod p
    if ip == 0 and not all(invariant(X) for X in (Bm, Am, Zm)):
        raise SystemExit('certificate not invariant under permutations of the copies (mod p): use the full checker')
    for name, X in (('Q_B', Bm), ('Q_A', Am), ('Z', Zm)):
        for beta in beta_list:
            cur = X
            for s in reversed(range(n)):
                cur = np.tensordot(cur, PsiP[beta[s]] % p_, axes=([s], [0])) % p_
            order = [2 * (n - 1 - s) for s in range(n)] + [2 * (n - 1 - s) + 1 for s in range(n)]
            b = int(np.prod([mult[j] for j in beta]))
            res[name][beta].append(np.transpose(cur, order).reshape(b, b).copy())
    del A_raw, Bm, Am, inner, Zm
    if ip % 10 == 9:
        log(f'  {ip + 1} primes done ({round(time.time() - T0, 1)} s)')


def crt(rs):
    out = np.zeros(rs[0].shape, dtype=object)
    for p_, r in zip(primes, rs):
        Mi = Mod // p_
        out = out + r.astype(object) * (Mi * pow(Mi, -1, p_))
    out = np.vectorize(lambda v: v % Mod, otypes=[object])(out)
    return np.vectorize(lambda v: v - Mod if v > Mod // 2 else v, otypes=[object])(out)


def pd_flint(G):
    k = G.shape[0]
    if k == 0:
        return True, None
    M = flint.fmpz_mat([[int(v) for v in row] for row in G])
    P, Lm, Dm, U = M.fflu()
    if P != flint.fmpz_mat([[1 if i == j else 0 for j in range(k)] for i in range(k)]):
        return None, 'row exchange'
    minors = [int(U[i, i]) for i in range(k)]
    return all(v > 0 for v in minors), min(minors)


def psd_flint(G):
    """G >= 0 iff (G definite) or (with N an exact kernel basis and coordinates `drop` on which N is invertible) the
    principal submatrix on the remaining coordinates is positive definite."""
    ok, info = pd_flint(G)
    if ok:
        return True, 0
    k = G.shape[0]
    M = flint.fmpz_mat([[int(v) for v in row] for row in G])
    Nm, nullity = M.nullspace()
    if nullity == 0:
        return (False, 0) if ok is False else (None, 0)
    N = flint.fmpz_mat([[Nm[i, j] for j in range(nullity)] for i in range(k)])
    assert M * N == flint.fmpz_mat(k, nullity)
    R, den_, rk = N.transpose().rref()
    assert rk == nullity
    drop = [next(j for j in range(k) if R[i, j] != 0) for i in range(rk)]
    keep = [i for i in range(k) if i not in drop]
    ok2, _ = pd_flint(G[np.ix_(keep, keep)])
    return ok2, nullity



def pd_congruence(G):
    k = G.shape[0]
    if k == 0:
        return True
    diag = [int(G[i, i]) for i in range(k)]
    if min(diag) <= 0:
        return False
    bl = [v.bit_length() for v in diag]
    Bmax = max(bl)
    c = [(Bmax - v) // 2 for v in bl]
    G1 = [[int(G[i, j]) << (c[i] + c[j]) for j in range(k)] for i in range(k)]
    e = max(abs(v).bit_length() for row in G1 for v in row) - 60
    e = max(e, 0)
    Fm = np.array([[float(v >> e) if v >= 0 else -float((-v) >> e) for v in row] for row in G1])
    w, V = np.linalg.eigh((Fm + Fm.T) / 2)
    if w.min() <= 0:
        return False
    Vi = np.rint(V * 2.0 ** 40).astype(np.int64)
    Vm = flint.fmpz_mat([[int(x) for x in row] for row in Vi])
    Hm = Vm.transpose() * flint.fmpz_mat(G1) * Vm
    for i in range(k):
        off = sum(abs(int(Hm[i, j])) for j in range(k) if j != i)
        if not int(Hm[i, i]) > off:
            return False
    return True


# forced kernel vectors in the S-adapted integer bases (from scratch): Bp_j c_j = Sym-projection of the m=0 spin-j space
Ssym_ = sp.zeros(27, 27)
for p3 in itertools.permutations(range(3)):
    for i3 in itertools.product(range(3), repeat=3):
        j3 = tuple(i3[p3[q]] for q in range(3))
        Ssym_[9 * j3[0] + 3 * j3[1] + j3[2], 9 * i3[0] + 3 * i3[1] + i3[2]] += 1
fco = {}
for j in (1, 3):
    Bj_ = sp.Matrix(Bp[j].tolist())
    img = Ssym_ * Bj_
    assert img.rank() == 1
    col = [q for q in range(img.shape[1]) if any(img[:, q])][0]
    cc = Bj_.solve_least_squares(img[:, col])
    assert Bj_ * cc == img[:, col]
    cc = cc * sp.ilcm(1, 1, *[sp.fraction(v)[1] for v in cc])
    fco[j] = [int(v) for v in cc]
log('forced vectors (S-adapted coordinates):', fco)


def pd_any(Gk):
    r = pd_flint(Gk)[0] if Gk.shape[0] <= 30 else None
    if r is None:
        r = pd_congruence(Gk)
    return r


def psd_check(G, beta, idx):
    """returns (ok, kernel dimension used); conservative: False means 'not certified'"""
    k = len(idx)
    if k == 0:
        return True, 0
    Gs = G[np.ix_(idx, idx)]
    if all(v == 0 for v in Gs.reshape(-1)):
        return True, k                                     # the zero matrix
    if all(j in (1, 3) for j in beta):
        f = np.array([1], dtype=object)
        for j in beta:
            f = np.multiply.outer(f, np.array(fco[j], dtype=object)).reshape(-1)
        fr = f[idx]
        if all(sum(Gs[i, q] * fr[q] for q in range(k)) == 0 for i in range(k)):
            d = max(range(k), key=lambda q: abs(fr[q]))
            assert fr[d] != 0
            keep = [q for q in range(k) if q != d]
            if pd_any(Gs[np.ix_(keep, keep)]):
                return True, 1
    elif pd_any(Gs):
        return True, 0
    if k <= 81:                                            # small blocks: general exact kernel method
        ok, nul = psd_flint(Gs)
        return bool(ok), nul
    return False, 0


allok = True
for name, symonly in (('Q_B', False), ('Q_A', True), ('Z', True)):
    ok_all, dims, kern = True, 0, 0
    for beta in beta_list:
        G = crt(res[name][beta])
        assert all(G[i, j] == G[j, i] for i in range(G.shape[0]) for j in range(i))
        rows = list(itertools.product(*[range(mult[j]) for j in beta]))
        rs = np.array([np.prod([sgn[beta[s]][r[s]] for s in range(n)]) for r in rows])
        if symonly:
            anti = np.where(rs == -1)[0]
            if len(anti) and any(G[i, j] != 0 for i in anti for j in range(G.shape[0])):
                log(f'  {name} block {beta}: antisymmetric part not zero')
                ok_all = False
                continue
            idx = np.where(rs == 1)[0]
        else:
            idx = np.arange(G.shape[0])
        ok, info = psd_check(G, beta, idx)
        if not ok:
            log(f'  {name} block {beta}: NOT positive semidefinite on the relevant part ({info})')
            ok_all = False
        else:
            kern += info
        dims += len(idx)
    log(f'{name}: positive semidefinite on all checked blocks (S=+1 part for Q_A, Z) = {ok_all}; '
        f'dimension checked {dims}, exact kernel dimension {kern}  ({round(time.time() - T0, 1)} s)')
    allok = allok and ok_all
if allok:
    log(f'PASS: for all x, y in C^{3 ** n}: <x(x)y| C_{a}^(x){n} |x(x)y> >= {t} |x|^2 |y|^2')
else:
    log('FAIL')
