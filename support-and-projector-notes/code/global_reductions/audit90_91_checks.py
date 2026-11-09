"""Checks of statements from an independent audit (OBSTRUCTIONS, map route and zero-free regions).
 (a) DDVV remainder: sum_I B_I(z)^2 = K_t^{(x)n}(zz^T), K_t(X) = (1+t^2)Tr(X) I + (5+2t) X + 2t X^T; Choi spectrum
     (t+4)^2, (t+1)^2, (t-1)^2; the fourth-power decomposition at the matching roots; Q_h <= c^{2n}||z||^4 beta^{n/2}.
 (b) P = G o W with W = (D - T)/2, Lambda = D - 2T, G = (D + 2 id + 2T)/7; the three-site seed
     343 G^{(x)3}(P_z) = B^{(x)3} + 32 (P_z - P_y), its spectrum, trace norm 353/343, the partial-inverse spectra of
     h = I - 2 P_y, and the integer inequalities used for every large N.
 (c) The local-gauge seed Y = P^{(x)3}(e e^T) (e the normalised Levi-Civita vector): spectrum and two-party marginal.
 (d) Zeros approaching an interior point x0 = -1/4 (n = 41) and the exact integer comparisons.
usage: python audit90_91_checks.py"""
import itertools
import numpy as np
rng = np.random.default_rng(9091)
I3 = np.eye(3)
F = np.zeros((9, 9)); Om = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
for i in range(3):
    for j in range(3): Om[4 * i, 4 * j] = 1
def apply_site(f, X, n, s):
    T = X.reshape([3] * (2 * n)); perm = [s, n + s] + [i for i in range(2 * n) if i not in (s, n + s)]
    Tp = np.transpose(T, perm).reshape(3, 3, -1)
    out = np.stack([f(Tp[:, :, k]) for k in range(Tp.shape[2])], axis=2)
    return np.transpose(out.reshape([3, 3] + [3] * (2 * n - 2)), np.argsort(perm)).reshape(3 ** n, 3 ** n)
def apply_all(f, X, n):
    for s in range(n): X = apply_site(f, X, n, s)
    return X
def choi(f):
    return sum(np.kron(np.outer(I3[i], I3[j]), f(np.outer(I3[i], I3[j]))) for i in range(3) for j in range(3))
# ---------------- (a)
def Bloc(i, k, t):
    return np.outer(I3[i], I3[k]) + np.outer(I3[k], I3[i]) + t * I3[i][k] * I3
err = 0
for t in (-1.3, -0.65, 0.4):
    Kt = lambda X, t=t: (1 + t * t) * np.trace(X) * I3 + (5 + 2 * t) * X + 2 * t * X.T
    for n in (1, 2):
        z = rng.normal(size=3 ** n)
        M = 0
        for I in itertools.product(range(3), repeat=n):
            Bm = 0
            for K in itertools.product(range(3), repeat=n):
                term = np.array([[1.0]])
                for s in range(n): term = np.kron(term, Bloc(I[s], K[s], t))
                Bm = Bm + z[np.ravel_multi_index(K, [3] * n)] * term
            M = M + Bm @ Bm
        err = max(err, np.abs(M - apply_all(Kt, np.outer(z, z), n)).max())
    ev = sorted(set(np.round(np.linalg.eigvalsh(choi(Kt)), 8)))
    print(f"(a) t = {t}: Choi eigenvalues {ev}, expected {sorted({round((t+4)**2, 8), round((t+1)**2, 8), round((t-1)**2, 8)})}")
print(f"(a) sum_I B_I^2 = K_t^(x)n (zz^T): max error {err:.1e}")
roots = np.roots([3, 8, 40, 68, 29]); roots = roots[np.abs(roots.imag) < 1e-9].real
for t in roots:
    c = 3 * t * t + 4 * t + 8; sig = (5 + 4 * t) / c; al = 5 / c
    Kn = lambda X, t=t, c=c: ((1 + t * t) * np.trace(X) * I3 + (5 + 2 * t) * X + 2 * t * X.T) / c
    K4 = lambda X: Kn(Kn(Kn(Kn(X))))
    E0 = lambda X: np.trace(X) * I3 / 3; Es = lambda X: (np.trace(X) * I3 + X + X.T) / 5; Ea = lambda X: (2 * np.trace(X) * I3 + X - X.T) / 6
    w0 = 1 - 2.5 * sig ** 4 - 3 * al ** 4
    X = rng.normal(size=(3, 3))
    dec = np.abs(K4(X) - (2.5 * sig ** 4 * Es(X) + 3 * al ** 4 * Ea(X) + w0 * E0(X))).max()
    beta = (2 + sig ** 4 + 3 * al ** 4) / 6
    worst = 0
    for n in (2, 3):
        for _ in range(20):
            z = rng.normal(size=3 ** n); z /= np.linalg.norm(z)
            Y = apply_all(Kn, np.outer(z, z), n)
            worst = max(worst, np.trace(Y @ Y) / beta ** (n / 2))
    print(f"(a) root t = {t:.6f}: |sigma| = {abs(sig):.4f} <= 11/30, alpha = {al:.4f} <= 3/4, fourth-power decomposition"
          f" error {dec:.1e}, weight of E0 {w0:.5f} > 0, beta = {beta:.5f} < 1/2;"
          f" max Q_h/(c^2n ||z||^4 beta^(n/2)) over random z, n = 2, 3: {worst:.4f} <= 1")
ok_sep = np.allclose(2 * np.eye(9) + Om - F,
                     2 * sum(np.outer(v, v.conj()) for i, j in itertools.combinations(range(3), 2) for sg in (1, -1)
                             for v in [np.kron((I3[i] + sg * 1j * I3[j]) / np.sqrt(2), (I3[i] - sg * 1j * I3[j]) / np.sqrt(2))])
                     + sum(np.outer(np.kron(I3[i], I3[j]), np.kron(I3[i], I3[j])) for i in range(3) for j in range(3) if i != j))
print(f"(a) product identity 2I + Omega - F = 2 sum |u x ubar><u x ubar| + sum_(i!=j) |ij><ij|: {ok_sep}")
# ---------------- (b)
D = lambda X: np.trace(X) * I3; Tm = lambda X: X.T
W = lambda X: (D(X) - X.T) / 2; Lam = lambda X: D(X) - 2 * X.T; G = lambda X: (D(X) + 2 * X + 2 * X.T) / 7
P = lambda X: (3 * D(X) - X - X.T) / 7
X = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
print(f"(b) Lambda W = id: {np.allclose(Lam(W(X)), X)}, W Lambda = id: {np.allclose(W(Lam(X)), X)}, G W = P: {np.allclose(G(W(X)), P(X))},"
      f" Choi of W = antisymmetric projector: {np.allclose(choi(W), (np.eye(9) - F) / 2)}")
z = np.zeros(27); y = np.zeros(27)
for word, s in {'000': 1, '011': -1, '101': -1, '110': -1}.items(): z[int(word, 3)] = s / 2
for word, s in {'001': 1, '010': 1, '100': 1, '111': -1}.items(): y[int(word, 3)] = s / 2
Pz = np.outer(z, z); Py = np.outer(y, y); Bd = np.diag([3.0, 3.0, 1.0])
GP = apply_all(G, Pz.astype(complex), 3).real
print(f"(b) 343 G^(x)3(P_z) = B^(x)3 + 32(P_z - P_y): error {np.abs(343 * GP - (np.kron(np.kron(Bd, Bd), Bd) + 32 * (Pz - Py))).max():.1e}")
ev = np.round(np.linalg.eigvalsh(343 * GP), 8)
print(f"(b) spectrum x343: {sorted(set(ev))} with multiplicities {[int((np.abs(ev - v) < 1e-6).sum()) for v in sorted(set(ev))]};"
      f" trace norm {np.abs(np.linalg.eigvalsh(GP)).sum():.10f} vs 353/343 = {353/343:.10f}")
h = np.eye(27) - 2 * Py
for k in range(4):
    sites = list(range(k)); Hk = h.astype(complex)
    for s in sites: Hk = apply_site(Lam, Hk, 3, s)
    evk = np.round(np.linalg.eigvalsh((Hk + Hk.conj().T).real / 2), 6)
    print(f"(b) Lambda on {k} site(s) of h: spectrum {sorted(set(evk))}, norm {np.abs(evk).max():.0f}")
print(f"(b) exact integer inequality 2*353^80 > 19*343^80: {2 * 353 ** 80 > 19 * 343 ** 80};"
      f" 2^3 < 9 and (19/2)^d > 9^d: {(19 / 2) > 9}")
# ---------------- (c)
e = np.zeros(27)
for (a_, b_, c_) in itertools.permutations(range(3)):
    e[9 * a_ + 3 * b_ + c_] = np.linalg.det(np.eye(3)[[a_, b_, c_]]) / np.sqrt(6)
Yc = apply_all(P, np.outer(e, e).astype(complex), 3).real
evc = np.round(np.linalg.eigvalsh(Yc) * 343, 6)
YAB = np.einsum('abcdec->abde', Yc.reshape(3, 3, 3, 3, 3, 3)).reshape(9, 9)
print(f"(c) Y = P^(x)3(ee^T): spectrum x343 {sorted(set(evc))} multiplicities {[int((np.abs(evc - v) < 1e-5).sum()) for v in sorted(set(evc))]};"
      f" Y_AB = (17I - F - Omega)/147: {np.allclose(YAB, (17 * np.eye(9) - F - Om) / 147)}; 17^12 > 2*16^12: {17 ** 12 > 2 * 16 ** 12}")
# ---------------- (d)
x0 = -0.25; R0 = 7 / 12; n = 41
def Gint(x):
    return (1 - R0 ** n) * ((x + 2) / 3) ** n + R0 ** n + (1 + R0 ** n) * ((x - 1) / 3) ** n
xz = -0.255134820851 + 0.133961133380j
for _ in range(50):
    d = (Gint(xz + 1e-8) - Gint(xz - 1e-8)) / 2e-8; xz = xz - Gint(xz) / d
print(f"(d) n = 41, x0 = -1/4: zero at {xz:.9f}; n (x - x0) = {n * (xz - x0):.4f} vs i pi (x0+2) = {1j * np.pi * (x0 + 2):.4f};"
      f" G(-1/2) = R0^n (1 - 2^(1-n)) > 0: {Gint(-0.5) > 0}")
print(f"(d) integer comparisons 12*1952^16 < 7*2023^16: {12 * 1952 ** 16 < 7 * 2023 ** 16}; 12*8^16 < 7*11^16: {12 * 8 ** 16 < 7 * 11 ** 16};"
      f" 3*1824^4 < 2*2047^4: {3 * 1824 ** 4 < 2 * 2047 ** 4}; 3*8^4 < 2*11^4: {3 * 8 ** 4 < 2 * 11 ** 4}")
