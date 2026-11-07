"""Exact verification of a level-2 DPS certificate for block-positivity of W = C_a^{(x)n}.

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
    QA_raw = np.array([Fraction(int(v), L) for v in cert['QA']], dtype=object).reshape([15] * n)
    QB = np.array([Fraction(int(v), L) for v in cert['QB']], dtype=object).reshape([15] * n)
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
    QA_raw = np.array([Fraction(int(v), L) for v in vA], dtype=object)[orb_].reshape([15] * n)
    QB = np.array([Fraction(int(v), L) for v in vB], dtype=object)[orb_].reshape([15] * n)

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


# ------------------------------------------------------------------ global exact algebra (coefficient tensors)
def perm_axes(X, pi):
    """Y[pi(a1),..,pi(an)] = X[a1,..,an]"""
    Y = np.empty_like(X)
    for idx in itertools.product(range(15), repeat=n):
        Y[tuple(pi[i] for i in idx)] = X[idx]
    return Y


def Pc(X):
    """P X P = (X + S X + X S + S X S)/4 with S = SWAP^{(x)n}"""
    SX = perm_axes(X, muL)
    XS = perm_axes(X, muR)
    SXS = perm_axes(SX, muR)
    return (X + SX + XS + SXS) / 4


def tpow(v):
    X = np.array([Fraction(1)], dtype=object)
    for _ in range(n):
        X = np.multiply.outer(X, np.array(v, dtype=object)).reshape(-1)
    return X.reshape([15] * n)


Wc = tpow(cvec)
Ic = tpow([Fraction(1) if b == iI else Fraction(0) for b in range(15)])
QA = Pc(QA_raw)
Z = Pc(Wc - t * Ic - perm_axes(QA, piA) - perm_axes(QB, piB))
log('Z computed exactly', round(time.time() - T0, 1), 's')


# ------------------------------------------------------------------ exact PSD checks
def to_int(X):
    den = 1
    for v in X.reshape(-1):
        den = math.lcm(den, v.denominator)
    Y = np.array([int(v * den) for v in X.reshape(-1)], dtype=object).reshape(X.shape)
    return Y, den


def block(Xint, beta):
    """sum_a X[a] Psi_{j1}(a1) (x) ... (x) Psi_{jn}(an), integer"""
    T = Xint  # shape (15,)*n
    # contract last copy first: T[a1..a_{n-1}, a_n] -> T[a1..a_{n-1}, r_n, s_n] ...
    ms = [mult[j] for j in beta]
    cur = T
    for s in reversed(range(n)):
        P_ = Psi[beta[s]]                     # (15, m, m)
        cur = np.tensordot(cur, P_, axes=([s], [0]))   # removes axis s, appends (m, m)
    # now cur has axes (r_n, s_n, r_{n-1}, s_{n-1}, ..., r_1, s_1)
    order_r = [2 * (n - 1 - s) for s in range(n)]
    order_c = [2 * (n - 1 - s) + 1 for s in range(n)]
    cur = np.transpose(cur, order_r + order_c)
    b = int(np.prod(ms))
    return cur.reshape(b, b)


def psd_exact(M):
    """exact test of positive semidefiniteness (symmetric pivoting); returns (ok, rank, min positive pivot)"""
    A = [[Fraction(int(x)) for x in row] for row in M]
    k = len(A)
    for i in range(k):
        for j in range(k):
            if A[i][j] != A[j][i]:
                return False, None, 'not symmetric'
    active = list(range(k))
    r = 0
    minpiv = None
    while active:
        p = max(active, key=lambda i: A[i][i])
        piv = A[p][p]
        if piv < 0:
            return False, r, 'negative pivot'
        if piv == 0:
            for i in active:
                for j in active:
                    if A[i][j] != 0:
                        return False, r, 'zero pivot with nonzero row'
            return True, r, minpiv
        active.remove(p)
        for i in active:
            f = A[i][p] / piv
            if f != 0:
                for j in active:
                    A[i][j] -= f * A[p][j]
        r += 1
        minpiv = piv if minpiv is None or piv < minpiv else minpiv
    return True, r, minpiv


def is_copy_symmetric(X):
    for perm in itertools.permutations(range(n)):
        if not np.array_equal(np.transpose(X, perm), X):
            return False
    return True


# If Q_A, Q_B and Z are invariant under permutations of the copies, a block beta and any permutation of it are
# conjugate by a permutation matrix, so it suffices to check blocks with j_1 <= ... <= j_n.
symmetric = all(is_copy_symmetric(X) for X in (QB, QA, Z))
log('certificate invariant under permutations of the copies (exact check):', symmetric)
if symmetric:
    beta_list = list(itertools.combinations_with_replacement(range(4), n))
else:
    beta_list = list(itertools.product(range(4), repeat=n))
log('blocks to check per operator:', len(beta_list))
allok = True
for name, X in [('Q_B', QB), ('Q_A', QA), ('Z', Z)]:
    Xi, den = to_int(X)
    ok_all = True
    ranks = 0
    worst = None
    for beta in beta_list:
        G = block(Xi, beta)
        ok, rk, mp = psd_exact(G)
        if not ok:
            log(f'  {name} block {beta}: NOT PSD ({mp})')
            ok_all = False
            continue
        ranks += rk
        if mp is not None:
            v = Fraction(mp) / den
            worst = v if worst is None or v < worst else worst
    log(f'{name}: PSD on all blocks = {ok_all}; total rank of compressions {ranks}; '
        f'smallest pivot {float(worst) if worst is not None else float("nan"):.3e}  ({round(time.time() - T0, 1)} s)')
    allok = allok and ok_all
log('(total rank counts only the blocks checked; all blocks together have dimension sum_beta b =', 7 ** n, ')')
if allok:
    log(f'PASS: for all x, y in C^{3 ** n}: <x(x)y| C_{a}^(x){n} |x(x)y> >= {t} |x|^2 |y|^2')
else:
    log('FAIL')
