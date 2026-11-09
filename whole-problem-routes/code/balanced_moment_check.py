"""Checks for the balanced degree-four obstruction (OBSTRUCTIONS.md) and the Choi-norm conversion (Proposition M').
(a) M = [3(D0+D1+D2)^{(x)3} - 10 T_2]/135 on three sites (D_i the pairings delta_ab delta_cd, delta_ac delta_bd,
    delta_ad delta_bc; T_2 = sum over S_3 of D_pi0 (x) D_pi1 (x) D_pi2): fully symmetric in the four global indices,
    positive semidefinite as a 729 x 729 matrix, rank 168, trace 39; with L = M/39: L(|z|^4) = 1, L(p_s) = 23/39,
    L(T_tetra) = -1/39, L(s_3) = -3/832 at a = 3/4 (64 s_3 = 25|z|^4 - 14 sum_s p_s + 18 T_tetra, checked on random z).
(b) the four class sums of 135*64*L*39 over the code quadruples (twelve synchronous rotations, simplex-code flags), and the
    sign of R_l = -14580 23^{15l} + 9369540 23^{7l} 15^{8l} - 303480 23^{7l} + 3585384 23^{3l} 15^{4l}: R_1 > 0 > R_2.
(c) the Choi matrix C_a = aI - (F + Omega)/4 has eigenvalues a-1, a-1/4, a+1/4, the four-linear form equals a Choi matrix
    element, and E_3(7/10) (19/20)^7 < 217031/1280000000 (exact).
usage: python balanced_moment_check.py"""
import itertools
from fractions import Fraction as Fr
import numpy as np

ok = True
# (a) the moment
d = 3
I3 = np.eye(d)
D = [np.einsum('ab,cd->abcd', I3, I3), np.einsum('ac,bd->abcd', I3, I3), np.einsum('ad,bc->abcd', I3, I3)]
Dsum = D[0] + D[1] + D[2]


def site3(A, B, C):            # tensor product over three sites, then regroup to global indices (27^4)
    T = np.einsum('abcd,efgh,ijkl->aeibfjcgkdhl', A, B, C)
    return T.reshape(27, 27, 27, 27)


M135 = 3 * site3(Dsum, Dsum, Dsum) - 10 * sum(site3(D[p[0]], D[p[1]], D[p[2]]) for p in itertools.permutations(range(3)))
M135 = np.rint(M135).astype(np.int64)
sym = all((M135 == np.transpose(M135, perm)).all() for perm in itertools.permutations(range(4)))
Mat = M135.reshape(729, 729) / 135
ev = np.linalg.eigvalsh(Mat)
rank = int((ev > 1e-9).sum()); trace = np.trace(Mat)
ok &= sym and ev.min() > -1e-9 and rank == 168 and abs(trace - 39) < 1e-9
print(f"(a) M fully symmetric: {sym}; min eigenvalue {ev.min():.2e}; rank {rank}; trace {trace:.6f}")
L = lambda T4: float(np.einsum('abcd,abcd->', M135, T4)) / 135 / 39       # T4 = coefficient tensor of a quartic
norm4 = np.einsum('ab,cd->abcd', np.eye(27), np.eye(27))
zz = np.eye(27).reshape(3, 3, 3, 27)
# p_s = Tr(rho_s^2) for site s: rho_s[i,j] = sum_rest z_{i,rest} z_{j,rest}
def purity_tensor(s):
    T = np.zeros((27,) * 4)
    for alpha, beta, gamma, delta in itertools.product(range(27), repeat=4):
        pass
    return T
idx = np.array(list(itertools.product(range(3), repeat=3)))             # digits of the 27 indices
def same_except(s):
    keep = [i for i in range(3) if i != s]
    return (idx[:, None, keep] == idx[None, :, keep]).all(-1)
ps_vals = []
for s in range(3):
    E = same_except(s)                                                   # alpha, beta agree off site s
    # p_s = sum_{alpha,beta,gamma,delta} z_alpha z_beta z_gamma z_delta [alpha~beta off s][gamma~delta off s][alpha_s = gamma_s? ...]
    # rho_s[i,j] = sum_r z_{(i,r)} z_{(j,r)} ; Tr rho_s^2 = sum_{i,j,r,r'} z_{ir} z_{jr} z_{jr'} z_{ir'}
    T = np.zeros((27, 27, 27, 27))
    for i in range(3):
        for j in range(3):
            for r in itertools.product(range(3), repeat=2):
                for r2 in itertools.product(range(3), repeat=2):
                    def code(x, rr):
                        dig = list(rr); dig.insert(s, x); return dig[0] * 9 + dig[1] * 3 + dig[2]
                    T[code(i, r), code(j, r), code(j, r2), code(i, r2)] += 1
    ps_vals.append(L(T))
# tetra T = sum z_{abc} z_{ab'c'} z_{a'bc'} z_{a'b'c}
Tt = np.zeros((27,) * 4)
for a, b, c, a2, b2, c2 in itertools.product(range(3), repeat=6):
    Tt[a * 9 + b * 3 + c, a * 9 + b2 * 3 + c2, a2 * 9 + b * 3 + c2, a2 * 9 + b2 * 3 + c] += 1
vals = (L(norm4), ps_vals, L(Tt))
good = abs(vals[0] - 1) < 1e-12 and all(abs(p - 23 / 39) < 1e-12 for p in ps_vals) and abs(vals[2] + 1 / 39) < 1e-12
ok &= good
# 64 s_3 identity on random z, and L(s_3)
rng = np.random.default_rng(0)
def Pa(X, a=0.75):
    return a * np.trace(X) * np.eye(3) - (X + X.T) / 4
def s3(z):
    X = np.outer(z, z).reshape(3, 3, 3, 3, 3, 3)
    for site in range(3):
        X = np.moveaxis(X, [site, 3 + site], [0, 1])
        Xs = X.reshape(3, 3, -1); tr = np.trace(Xs, axis1=0, axis2=1)
        Xs = 0.75 * np.einsum('ij,k->ijk', np.eye(3), tr) - (Xs + Xs.transpose(1, 0, 2)) / 4
        X = np.moveaxis(Xs.reshape(X.shape), [0, 1], [site, 3 + site])
    return z @ X.reshape(27, 27) @ z
err = 0
for _ in range(5):
    z = rng.normal(size=27); Z = z.reshape(3, 3, 3)
    ps = [np.trace(np.moveaxis(Z, s, 0).reshape(3, 9) @ np.moveaxis(Z, s, 0).reshape(3, 9).T @ np.moveaxis(Z, s, 0).reshape(3, 9) @ np.moveaxis(Z, s, 0).reshape(3, 9).T) for s in range(3)]
    tet = np.einsum('abc,aBC,AbC,ABc->', Z, Z, Z, Z)
    err = max(err, abs(64 * s3(z) - (25 * (z @ z) ** 2 - 14 * sum(ps) + 18 * tet)))
Ls3 = (25 * vals[0] - 14 * sum(ps_vals) + 18 * vals[2]) / 64
good2 = err < 1e-9 and abs(Ls3 + 3 / 832) < 1e-12
ok &= good2
print(f"    L(|z|^4) = {vals[0]:.6f}, L(p_s) = {[round(p * 39, 6) for p in ps_vals]}/39, L(T_tetra) = {vals[2] * 39:.6f}/39; "
      f"identity error {err:.1e}; L(s_3) = {Ls3:.10f} (-3/832 = {-3 / 832:.10f}): {good and good2}")

# (b) code contraction: classes of quadruples
S = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]); signs = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
Os = [np.linalg.matrix_power(S, r) @ np.diag(e) for r in range(3) for e in signs]
Rs = [np.kron(np.kron(O, O), O) for O in Os]                            # synchronous rotations on three sites
vs = [v for v in itertools.product(range(2), repeat=4) if any(v)]
code = [[sum(x * y for x, y in zip([(j >> b) & 1 for b in range(4)], v)) % 2 for v in vs] for j in range(12)]
# kernel 64 P_{3/4}^{(x)3}: K1[mu nu, rho sigma] = 3 d_mu nu d_rho sigma - d_mu rho d_nu sigma - d_mu sigma d_nu rho
K1 = 3 * np.einsum('mn,rs->mnrs', I3, I3) - np.einsum('mr,ns->mnrs', I3, I3) - np.einsum('ms,nr->mnrs', I3, I3)
K3 = np.einsum('abcd,efgh,ijkl->aeibfjcgkdhl', K1, K1, K1).reshape(27, 27, 27, 27)   # [mu, nu, rho, sigma] global
def core(j, k, l, m):
    # f = sum K3[mu,nu,rho,sigma] (R_l z)_mu (R_m z)_nu (R_j z)_rho (R_k z)_sigma ; L-value times 135*39
    Tf = np.einsum('mnrs,ma,nb,rc,sd->abcd', K3, Rs[l], Rs[m], Rs[j], Rs[k], optimize=True)
    return int(np.rint(np.einsum('abcd,abcd->', M135, Tf)))
def cls(j, k, l, m):
    bits = [(code[j][p], code[k][p], code[l][p], code[m][p]) for p in range(15)]
    if any(sum(b) % 2 for b in bits):
        return None
    if j == k == l == m:
        return 'main'
    if j == k and l == m:
        return 'B'
    if (j == l and k == m) or (j == m and k == l):
        return 'K'
    if all(len(set(x)) > 1 for x in [(j, k, l, m)]) and len({j, k, l, m}) == 4:
        return 'four'
    return None
sums = {'main': 0, 'B': 0, 'K': 0, 'four': 0}; counts = {k: 0 for k in sums}
for j, k, l, m in itertools.product(range(12), repeat=4):
    c = cls(j, k, l, m)
    if c is None:
        continue
    counts[c] += 1; sums[c] += core(j, k, l, m)
want = {'main': -14580, 'B': 9369540, 'K': -303480, 'four': 3585384}
good = sums == want and counts == {'main': 12, 'B': 132, 'K': 264, 'four': 936}
ok &= good
R = lambda l: -14580 * 23 ** (15 * l) + 9369540 * 23 ** (7 * l) * 15 ** (8 * l) - 303480 * 23 ** (7 * l) + 3585384 * 23 ** (3 * l) * 15 ** (4 * l)
Dd = lambda l: 144 * 135 * 64 * 39 * 48 ** (15 * l)
val63 = Fr(R(2), Dd(2)) / Fr(3, 4) ** 63
signs_ok = R(1) > 0 > R(2)
ok &= signs_ok
print(f"(b) class counts {counts}; class sums {sums} (stated {want}): {good}; R_1 > 0 > R_2: {signs_ok}; "
      f"L_2(s_63)/(3/4)^63 = {float(val63):.6e}")

# (c) Choi norm and the four-linear identity; the 97-copy calibration
for a in (0.5, 0.75, 0.9):
    F = np.zeros((9, 9)); Om = np.zeros((9, 9))
    for i, j in itertools.product(range(3), repeat=2):
        F[3 * i + j, 3 * j + i] = 1
    om = np.eye(3).reshape(9); Om = np.outer(om, om)
    Ca = a * np.eye(9) - (F + Om) / 4
    evs = sorted(set(np.round(np.linalg.eigvalsh(Ca), 10)))
    x1, x2, y1, y2 = rng.normal(size=(4, 3))
    lhs = y1 @ Pa(np.outer(x1, x2), a) @ y2
    rhs = np.kron(x1, y1) @ Ca @ np.kron(x2, y2)
    ok &= abs(lhs - rhs) < 1e-12 and np.allclose(evs, sorted({a - 1, a - 0.25, a + 0.25}))
print(f"(c) Choi eigenvalues a-1, a-1/4, a+1/4 and the four-linear identity y1^T P_a(x1 x2^T) y2 = <x1 (x) y1, C_a (x2 (x) y2)>: checked")
a = Fr(7, 10); A = a * a - a / 3 + Fr(1, 6); B = a * a - a / 3; K = (5 - 8 * a) / 48
r, h = B / A, abs(K) / A
E3 = 11 * r ** 24 + 22 * h ** 24 + 78 * r ** 12 * h ** 24
cal = E3 * Fr(19, 20) ** 7 < Fr(217031, 1280000000)
ok &= cal
print(f"    E_3(7/10) (19/20)^7 = {float(E3 * Fr(19, 20) ** 7):.6e} < 1.69555e-4: {cal}  (97 = 7 + 30*3 copies)")
print("PASS" if ok else "FAIL")
