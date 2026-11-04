"""Checks for Proposition 9 (single-excitation-range/README.md): one copy of rank two, all others of rank one.
Copy 0: U_0 = diag(s1, s2); light copies t = 1..N: U_t = |1><conj(x_t)|, sum_t x_t x_t^* = I - Abar_0.
Reduced operators (y_t = Abar_0^{-1/2} x_t, light copies as qubits, levels k = number of excited light copies):
  X'_k = sum_t y_t y_t^* (x) n_t^hole,  b_k = sum_t y_t^* (x) sigma_t^- : V (x) H_k -> H_{k-1},
  M = 4 X'_{c-1} + b_{c-1}^* b_{c-1},  K = (1/4) b_c (I + X'_c)^{-1} b_c^*.
(a) reduction: 2Kt_c >= 0  <=>  (L) Tr_V (2 + M)^{-1} <= I/2 + (2 + 4K)^{-1}; we compare the sign of lambda_min(2Kt_c)
    (built independently with GeneralModel) with the sign of the (L) margin at random and nearly singular data;
(b) one copy carrying a full direction (s_max = 1, all x_t parallel): the commutator identity
    B B^* - B^* B = D - U, the bound m - 4K >= lambda I, and lambda_min(2Kt_c) >= 0 directly, with rotated frames,
    N <= 7, every c;
(c) the near-parallel condition nu (h + (c+3)^2 mu / 2) <= mu (mu = s_min^-2 - 1, nu = s_max^-2 - 1): at random data
    satisfying it, lambda_min(2Kt_c) >= 0;
(d) two routes that cannot prove (a) in general: a separable split M >= rho (x) 4K (rho a density matrix on V) fails at
    an instance with N = 3, c = 3, and (2 + tP)^{-1} + (2 + tQ)^{-1} <= 1/2 + (2 + t(P+Q))^{-1} fails for rank-one P, Q."""
import numpy as np, itertools, warnings
from threelevel_certificate import GeneralModel
warnings.filterwarnings('ignore')
rng = np.random.default_rng(12)

def levels(N, k):
    S = [frozenset(L) for L in itertools.combinations(range(N), k)]; return S, {L: i for i, L in enumerate(S)}
def ops(y, k):
    """X'_k (2m x 2m, index r*m + j) and b_k (m_{k-1} x 2m) for y of shape (N, 2)."""
    N = len(y); S, p = levels(N, k); m = len(S)
    X = np.zeros((2 * m, 2 * m), complex)
    for j, L in enumerate(S):
        Y = sum((np.outer(y[t], y[t].conj()) for t in range(N) if t not in L), np.zeros((2, 2)))
        for r in range(2):
            for rp in range(2): X[r * m + j, rp * m + j] = Y[r, rp]
    if k == 0: return X, np.zeros((0, 2 * m), complex)
    R, pr = levels(N, k - 1); b = np.zeros((len(R), 2 * m), complex)
    for j, L in enumerate(S):
        for t in L:
            for r in range(2): b[pr[L - {t}], r * m + j] += np.conj(y[t][r])
    return X, b
def L_margin(y, c):
    """lambda_max of Tr_V(2+M)^{-1} - I/2 - (2+4K)^{-1} on H_{c-1} (<= 0 iff (L))."""
    N = len(y); Xk, bk = ops(y, c - 1); m = Xk.shape[0] // 2
    M = 4 * Xk + bk.conj().T @ bk; Minv = np.linalg.inv(2 * np.eye(2 * m) + M); TrV = Minv[:m, :m] + Minv[m:, m:]
    if c <= N:
        Xc, bc = ops(y, c); K = 0.25 * bc @ np.linalg.inv(np.eye(Xc.shape[0]) + Xc) @ bc.conj().T
    else: K = np.zeros((m, m))
    Z = TrV - 0.5 * np.eye(m) - np.linalg.inv(2 * np.eye(m) + 4 * K)
    return np.linalg.eigvalsh((Z + Z.conj().T) / 2)[-1]
def model(s1, s2, x, rot=None):
    """GeneralModel with U_0 = diag(s1, s2) and light U_t = |1><conj(x_t)|, optionally in a rotated frame."""
    U = [np.diag([s1, s2]).astype(complex)] + [np.array([[xt[0], xt[1]]], complex) for xt in x]
    if rot is not None: U = [u @ rot for u in U]
    return GeneralModel(U)
def lights(s1, s2, N, parallel=False):
    W = np.diag([1 - s1 ** 2, 1 - s2 ** 2])
    v = rng.normal(size=(N, 2)) + 1j * rng.normal(size=(N, 2))
    if parallel: v[:, 1] = 0 if s2 == 1 else v[:, 1]; v[:, 0] = v[:, 0] if s2 == 1 else 0
    G = sum(np.outer(t, t.conj()) for t in v); w, V = np.linalg.eigh(G); Gm = V @ np.diag(np.where(w > 1e-14, w ** -0.5, 0)) @ V.conj().T
    return np.array([np.sqrt(W) @ Gm @ t for t in v])
def haar2():
    Z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)); Q, R = np.linalg.qr(Z); return Q * (np.diag(R) / abs(np.diag(R)))
ok = True
# (a)
agree = 0; total = 0; worst_gap = 0
for tr in range(80):
    N = int(rng.integers(2, 6)); c = int(rng.integers(1, N + 2))
    s = np.sort(rng.uniform(0.05, 0.999, 2))[::-1] if tr % 3 else np.sort([rng.uniform(0.9, 0.9999), rng.uniform(0.01, 0.3)])[::-1]
    x = lights(s[0], s[1], N); y = np.array([np.diag([1 / s[0], 1 / s[1]]) @ xt for xt in x])
    lm = 2 * np.linalg.eigvalsh(model(s[0], s[1], x, haar2()).ops(c)['K'])[0]; Lm = L_margin(y, c)
    total += 1; agree += (lm >= -1e-10) == (Lm <= 1e-10); worst_gap = max(worst_gap, Lm)
print(f"(a) sign of lambda_min(2Kt_c) matches (L) in {agree}/{total} random cases; largest (L) margin {worst_gap:+.2e} (<= 0 means (L) holds)")
ok &= agree == total
# (b)
werr = 0; wmin = np.inf; dmin = np.inf
for tr in range(50):
    N = int(rng.integers(2, 8)); c = int(rng.integers(1, N + 2)); d = rng.uniform(1e-3, 0.999)
    s1, s2 = 1.0, np.sqrt(1 - d)
    eta = (rng.uniform(0.05, 1, N) * np.exp(2j * np.pi * rng.uniform(size=N))); eta *= np.sqrt(d) / np.linalg.norm(eta)
    x = np.zeros((N, 2), complex); x[:, 1] = eta                      # all light copies along the second frame vector
    y = np.array([np.diag([1 / s1, 1 / s2]) @ xt for xt in x]); lam = np.sum(np.abs(y) ** 2)
    k = c - 1; S, p = levels(N, k); m = len(S)
    B = ops(y, k)[1][:, m:] if k >= 1 else np.zeros((0, m)); Bc = ops(y, c)[1][:, levels(N, c)[0].__len__():] if c <= N else np.zeros((m, 0))
    D = np.diag([sum(abs(y[t][1]) ** 2 for t in range(N) if t not in L) for L in S]); Uu = np.diag([sum(abs(y[t][1]) ** 2 for t in L) for L in S])
    werr = max(werr, np.abs(Bc @ Bc.conj().T - B.conj().T @ B - (D - Uu)).max())
    Dc = np.diag([sum(abs(y[t][1]) ** 2 for t in range(N) if t not in L) for L in levels(N, c)[0]]) if c <= N else np.zeros((0, 0))
    K4 = Bc @ np.linalg.inv(np.eye(Dc.shape[0]) + Dc) @ Bc.conj().T if c <= N else np.zeros((m, m))
    mm = 4 * D + B.conj().T @ B; wmin = min(wmin, np.linalg.eigvalsh(mm - K4 - lam * np.eye(m))[0])
    dmin = min(dmin, 2 * np.linalg.eigvalsh(model(s1, s2, x, haar2()).ops(c)['K'])[0])
print(f"(b) commutator identity error {werr:.1e}; min eig of m - 4K - lambda I {wmin:+.1e}; min lambda_min(2Kt_c) {dmin:+.3e} (50 POVMs, N <= 7)")
ok &= werr < 1e-12 and wmin > -1e-10 and dmin > -1e-10
# (c)
cnt = 0; cmin = np.inf
for tr in range(4000):
    N = int(rng.integers(2, 6)); c = int(rng.integers(2, N + 1)); h = N + 1 - c
    smin = rng.uniform(0.2, 0.9); smax = rng.uniform(0.95, 0.99999)
    mu, nu = smin ** -2 - 1, smax ** -2 - 1
    if nu * (h + (c + 3) ** 2 * mu / 2) > mu: continue
    x = lights(smax, smin, N); cnt += 1
    cmin = min(cmin, 2 * np.linalg.eigvalsh(model(smax, smin, x, haar2()).ops(c)['K'])[0])
    if cnt >= 120: break
print(f"(c) {cnt} random configurations satisfy the near-parallel condition; min lambda_min(2Kt_c) {cmin:+.3e}")
ok &= cnt >= 50 and cmin > -1e-10
# (d)
import cvxpy as cp
best = np.inf
for tr in range(40):
    y = rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2)); c = 3
    Xk, bk = ops(y, c - 1); m = Xk.shape[0] // 2; M = 4 * Xk + bk.conj().T @ bk
    K4 = np.zeros((m, m))                                   # c = N: no level-c creation, K = 0 ... use c = 3 with N = 4 instead
    y = rng.normal(size=(4, 2)) + 1j * rng.normal(size=(4, 2)) if tr % 2 else y
    N = len(y); c = 3
    Xk, bk = ops(y, c - 1); m = Xk.shape[0] // 2; M = 4 * Xk + bk.conj().T @ bk
    Xc, bc = ops(y, c); K4 = bc @ np.linalg.inv(np.eye(Xc.shape[0]) + Xc) @ bc.conj().T
    rho = cp.Variable((2, 2), hermitian=True); t = cp.Variable()
    cons = [rho >> 0, cp.real(cp.trace(rho)) == 1, M - cp.kron(rho, K4) - t * np.eye(2 * m) >> 0]
    cp.Problem(cp.Maximize(t), cons).solve(solver='CLARABEL')
    best = min(best, t.value / max(np.linalg.norm(K4, 2), 1e-12))
P = np.outer([1, 0], [1, 0]); q = np.array([1, 1]) / np.sqrt(2); Q = np.outer(q, q); tt = 3.0
Z = np.linalg.inv(2 * np.eye(2) + tt * P) + np.linalg.inv(2 * np.eye(2) + tt * Q) - 0.5 * np.eye(2) - np.linalg.inv(2 * np.eye(2) + tt * (P + Q))
sub = np.linalg.eigvalsh(Z)[-1]
print(f"(d) best separable-split margin / ||4K|| over 40 instances: {best:+.3e} (< 0: the split fails); subadditivity violated by {sub:+.3f}")
ok &= best < 0 and sub > 0
assert ok
print("PASS")
