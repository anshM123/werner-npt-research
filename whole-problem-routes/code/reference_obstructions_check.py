"""Checks of the quadratic-reference obstructions for the trace-norm form of the map route (Section 5, "Norm bounds";
from an independent audit, re-checked here).  P = (3D - id - T)/7 is the normalised critical map.
 (a) Bell pair: Y = (P x P)(phi phi^*) = (15 I + 2F + 2 Omega)/147, eigenvalues 23, 17, 13 (/147), 9 Tr Y^2 = 2481/2401;
 (b) GHZ blocks: two-site diagonal marginal 19/147, 15/147, classical cost 2433/2401; GHZ_L outputs positive (L <= 6);
 (c) the extension C = (6D - id - 3T)/14: CP (Choi eigenvalues 0, 3/14, 9/14), equal to P on symmetric matrices,
     Q(Y|Z) = 116561/116424, ||Y - Z||^2 = 1/7203;
 (d) random CP maps Q equal to P on symmetric matrices: ||(Q x Q)(phi phi^*) - Y||^2 >= e^4/3 = 1/7203;
 (e) antisymmetric three-qutrit state h: Y_a = P_a^{x3}(hh^*) > 0, and for random CP extensions Q1, Q2, Q3,
     ||(Q1 x Q2 x Q3)(hh^*) - Y_a||^2 >= (3/8) r^2 e^4 (= 3/235298 at a = 3/4).
usage: python reference_obstructions_check.py"""
from fractions import Fraction as Fr
import itertools
import numpy as np
rng = np.random.default_rng(7)
I3 = np.eye(3)
def Pmap(a):
    t = 3 * a - 0.5
    return lambda X: (a * np.trace(X) * I3 - (X + X.T) / 4) / t
def apply_site(f, X, n, s):
    T = X.reshape((3,) * (2 * n))
    perm = [s, n + s] + [i for i in range(2 * n) if i not in (s, n + s)]
    Tp = np.transpose(T, perm).reshape(3, 3, -1)
    out = np.stack([f(Tp[:, :, k]) for k in range(Tp.shape[2])], axis=2)
    return np.transpose(out.reshape([3, 3] + [3] * (2 * n - 2)), np.argsort(perm)).reshape(3 ** n, 3 ** n)
def apply_all(fs, X):
    n = len(fs)
    for s, f in enumerate(fs): X = apply_site(f, X, n, s)
    return X
P = Pmap(0.75)
phi = np.zeros(9); phi[[0, 4, 8]] = 1 / np.sqrt(3)
Y = apply_all([P, P], np.outer(phi, phi).astype(complex))
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
Om = np.outer(phi, phi) * 3
print(f"(a) |Y - (15I+2F+2Om)/147| = {np.abs(Y - (15 * np.eye(9) + 2 * F + 2 * Om) / 147).max():.1e};"
      f" eigenvalues x147: {sorted(set(np.round(np.linalg.eigvalsh(Y).real * 147, 6)))}; 9TrY^2 = {9 * np.trace(Y @ Y).real:.10f}"
      f" vs 2481/2401 = {2481 / 2401:.10f}")
p = {True: Fr(1, 7), False: Fr(3, 7)}
pjj = sum(p[i == 0] * p[i == 0] for i in range(3)) / 3
pjk = sum(p[i == 0] * p[i == 1] for i in range(3)) / 3
cost = 9 * (3 * pjj ** 2 + 6 * pjk ** 2)
print(f"(b) two-site diagonal marginal {pjj}, {pjk}; classical cost {cost} (= 2433/2401: {cost == Fr(2433, 2401)})")
for L in range(2, 7):
    g = np.zeros(3 ** L)
    for i in range(3): g[sum(i * 3 ** k for k in range(L))] = 1 / np.sqrt(3)
    YL = apply_all([P] * L, np.outer(g, g).astype(complex))
    print(f"(b) GHZ_{L}: least output eigenvalue {np.linalg.eigvalsh(YL).min():.3e}")
def Cext(X):
    return (6 * np.trace(X) * I3 - X - 3 * X.T) / 14
def choi(f):
    return sum(np.kron(np.outer(I3[i], I3[j]), f(np.outer(I3[i], I3[j]).astype(complex))) for i in range(3) for j in range(3))
S = rng.normal(size=(3, 3)); S = S + S.T
print(f"(c) Choi eigenvalues of C: {sorted(set(np.round(np.linalg.eigvalsh(choi(Cext)).real, 6)))};"
      f" |C(S) - P(S)| on a random symmetric S: {np.abs(Cext(S) - P(S)).max():.1e}")
Z = apply_all([Cext, Cext], np.outer(phi, phi).astype(complex))
Qc = np.trace(Y @ np.linalg.inv(Z) @ Y).real
print(f"(c) Q(Y|Z) = {Qc:.12f} vs 116561/116424 = {116561 / 116424:.12f}; ||Y-Z||^2 = {np.linalg.norm(Y - Z) ** 2:.12f}"
      f" vs 1/7203 = {1 / 7203:.12f}")
# antisymmetric Hermitian basis A_k = i(E_ij - E_ji)/sqrt2 and random CP extensions
pairs = [(1, 2), (2, 0), (0, 1)]
A = [1j * (np.outer(I3[i], I3[j]) - np.outer(I3[j], I3[i])) / np.sqrt(2) for (i, j) in pairs]
def ext_from(a, B):
    Pa = Pmap(a)
    def Q(X):
        sym = (X + X.T) / 2
        return Pa(sym) + sum(np.trace(A[k] @ X) * B[k] for k in range(3))   # <A_k, X> = Tr(A_k X)
    return Q
def rand_herm_traceless():
    M = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); M = M + M.conj().T
    return M - np.trace(M) / 3 * I3
def random_extension(a):
    """bisect from the O(3)-covariant extension B_k = e A_k towards a random direction up to the CP boundary"""
    t = 3 * a - 0.5; e = (1 - a) / t
    B0 = [e * A[k] for k in range(3)]
    while True:
        Dl = [rand_herm_traceless() for _ in range(3)]
        if sum(np.trace(A[k] @ Dl[k]).real for k in range(3)) > 0: break
    lo, hi = 0.0, 1.0
    while np.linalg.eigvalsh(choi(ext_from(a, [B0[k] + hi * Dl[k] for k in range(3)]))).min() > -1e-12: hi *= 2
    for _ in range(40):
        mid = (lo + hi) / 2
        if np.linalg.eigvalsh(choi(ext_from(a, [B0[k] + mid * Dl[k] for k in range(3)]))).min() > -1e-12: lo = mid
        else: hi = mid
    s = lo * rng.uniform(0.3, 1.0)
    return ext_from(a, [B0[k] + s * Dl[k] for k in range(3)])
worst = np.inf
for _ in range(200):
    Q = random_extension(0.75)
    Zq = apply_all([Q, Q], np.outer(phi, phi).astype(complex))
    worst = min(worst, np.linalg.norm(Zq - Y) ** 2 * 7203)
print(f"(d) 200 random CP extensions: min ||(QxQ)(phi phi*) - Y||^2 / (1/7203) = {worst:.4f} (>= 1 expected)")
eps = np.zeros(27)
for (a_, b_, c_) in itertools.permutations(range(3)):
    sgn = np.linalg.det(np.eye(3)[[a_, b_, c_]])
    eps[9 * a_ + 3 * b_ + c_] = sgn
h = eps / np.sqrt(6)
for a in (0.75, 0.8, 0.9, 0.99):
    t = 3 * a - 0.5; r = 1 / (2 * t); e = (1 - a) / t
    Ya = apply_all([Pmap(a)] * 3, np.outer(h, h).astype(complex))
    gap = 3 / 8 * r * r * e ** 4
    worst = np.inf
    nt = 60 if a == 0.75 else 20
    for _ in range(nt):
        Qs = [random_extension(a) for _ in range(3)]
        Zq = apply_all(Qs, np.outer(h, h).astype(complex))
        worst = min(worst, np.linalg.norm(Zq - Ya) ** 2 / gap)
    Co = ext_from(a, [e * A[k] for k in range(3)])
    Zo = apply_all([Co] * 3, np.outer(h, h).astype(complex))
    print(f"(e) a={a}: least eigenvalue of Y_a {np.linalg.eigvalsh(Ya).min():.4e}; gap (3/8)r^2e^4 = {gap:.4e};"
          f" covariant extension ratio {np.linalg.norm(Zo - Ya) ** 2 / gap:.4f}; {nt} random triples min ratio {worst:.4f}")
print(f"(e) gap at a = 3/4 equals 3/235298: {abs(3 / 8 * (2 / 7) ** 2 * (1 / 7) ** 4 - 3 / 235298) < 1e-18}")
