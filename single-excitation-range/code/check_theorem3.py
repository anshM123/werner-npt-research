"""Checks for Theorem 3 (the co-range sector c = n - 1, every n) and Corollary 4 (n = 4).
(a) The identity  ||D_w G||^2 + ||Gamma||^2 = sum_t [ sum_{s != t} ||kappa_s G_t||^2 + ||w_t||^2 ||G_t||^2 ],
    (D_w G)_{st} = kappa_s G_t - kappa_t G_s, Gamma = sum_t w_t (x)_t G_t, kappa_s = <w_s| on slot s,
    for random complex w_s in C^m and G_t, n = 2..5, m = 2, 3.
(b) lambda_max(D^* D) <= 2 for the full map (two reference components) at random and adversarial U with
    sum_s U_s^* U_s = I, n = 3, 4, 5; and the sharp example (real tight frame) gives equality.
(d) For one reference index and every sector c, L_r = T + V - ||Gamma||^2 <= ||G||^2 (the remark after Corollary 4).
(e) Corollary 4 end to end at n = 4 from the dense operator, at random and minimised planes.
The compact form of the sector c = n - 1 (the hole form with h = 1) is checked against explicit partial traces by
check_wplanes.py, part (c)."""
import numpy as np, itertools
from scipy.optimize import minimize
rng = np.random.default_rng(3); ok = True
def cz(*s): return rng.normal(size=s) + 1j * rng.normal(size=s)
# ---------- (a)
def kappa(w, G, s, t, n):
    """contract slot s of G_t (hole at t) with conj(w): result has holes at s, t (slot order of the remaining copies)"""
    rem = [j for j in range(n) if j != t]; ax = rem.index(s)
    return np.tensordot(w.conj(), G, axes=([0], [ax]))
def insert(w, G, t, n):
    """place w at slot t of G (hole at t): result in E^{(x)n}"""
    out = np.tensordot(w, G, axes=0)                    # axis 0 = slot t, then the n-1 remaining slots in order
    rem = [j for j in range(n) if j != t]; order = [t] + rem; perm = np.argsort(order)
    return np.transpose(out, perm)
err_a = 0
for n in (2, 3, 4, 5):
    for m in (2, 3):
        w = [cz(m) for _ in range(n)]; G = [cz(*([m] * (n - 1))) for _ in range(n)]
        # compute ||D_w G||^2 carefully: both kappa_s G_t and kappa_t G_s live on the copies other than s, t (in increasing order)
        DD = 0
        for s, t in itertools.combinations(range(n), 2):
            a = kappa(w[s], G[t], s, t, n); b = kappa(w[t], G[s], t, s, n); DD += np.linalg.norm(a - b) ** 2
        Gam = sum(insert(w[t], G[t], t, n) for t in range(n)); right = 0
        for t in range(n):
            right += sum(np.linalg.norm(kappa(w[s], G[t], s, t, n)) ** 2 for s in range(n) if s != t) + np.linalg.norm(w[t]) ** 2 * np.linalg.norm(G[t]) ** 2
        err_a = max(err_a, abs(DD + np.linalg.norm(Gam) ** 2 - right) / right)
print(f"(a) identity, n = 2..5, m = 2, 3: largest relative error {err_a:.1e}"); ok &= err_a < 1e-12
# ---------- (b) the full operator D on (+)_u R (x) E^{(x)(n-1)}
def D_matrix(U, n):
    m = 2; dimF = 2 * m ** (n - 1); blocks = []
    pairs = list(itertools.combinations(range(n), 2)); rows = []
    for (s, t) in pairs:
        row = []
        for u in range(n):
            Mu = np.zeros((m ** (n - 2), dimF), complex)
            for k in range(dimF):
                e = np.zeros(dimF, complex); e[k] = 1; F = e.reshape([2] + [m] * (n - 1))     # F_u[r; slots]
                if u == t:      # + K_s F_t, K_s contracts slot s and r with U_s: sum_{r,z} conj(U_s[z,r]) F[r; z at s]
                    val = sum(kappa(U[s][:, r], F[r], s, t, n) for r in range(2)); Mu[:, k] = val.ravel()
                elif u == s:    # - K_t F_s
                    val = sum(kappa(U[t][:, r], F[r], t, s, n) for r in range(2)); Mu[:, k] = -val.ravel()
            row.append(Mu)
        rows.append(np.hstack(row))
    return np.vstack(rows)
def normalise(U):
    S = sum(u.conj().T @ u for u in U); w, V = np.linalg.eigh(S); T = V @ np.diag(w ** -0.5) @ V.conj().T; return [u @ T for u in U]
worst = 0
for n in (3, 4, 5):
    for k in range(4):
        U = normalise([cz(2, 2) for _ in range(n)]); Dm = D_matrix(U, n); worst = max(worst, np.linalg.eigvalsh(Dm.conj().T @ Dm)[-1])
    def negl(x):
        U = normalise([(x[4 * s:4 * s + 4] + 1j * x[4 * n + 4 * s:4 * n + 4 * s + 4]).reshape(2, 2) for s in range(n)]); Dm = D_matrix(U, n)
        return -np.linalg.eigvalsh(Dm.conj().T @ Dm)[-1]
    if n <= 4:
        res = minimize(negl, rng.normal(size=8 * n), method='L-BFGS-B', options={'maxiter': 150}); worst = max(worst, -res.fun)
print(f"(b) largest lambda_max(D^*D) over random and adversarial U at n = 3, 4, 5: {worst:.10f} (bound 2)"); ok &= worst <= 2 + 1e-9
# sharp example: real tight frame v_s (sum v v^T = I), unit z_s, U_s = |z_s><v_s|, F_t = (J v_t) (x) (x)_{j != t} z_j
for n in (3, 4, 5):
    ang = np.pi * np.arange(n) / n; v = [np.sqrt(2.0 / n) * np.array([np.cos(a), np.sin(a)]) for a in ang]
    z = [(lambda x: x / np.linalg.norm(x))(cz(2)) for _ in range(n)]; U = [np.outer(z[s], v[s]) for s in range(n)]
    J = np.array([[0, 1], [-1, 0]])
    F = []
    for t in range(n):
        vec = J @ v[t]
        for j in range(n):
            if j != t: vec = np.kron(vec, z[j])
        F.append(vec)
    Fv = np.concatenate(F); Dm = D_matrix(U, n); ratio = np.linalg.norm(Dm @ Fv) ** 2 / np.linalg.norm(Fv) ** 2
    SA = sum(u.conj().T @ u for u in U); print(f"    sharp example n={n}: ||D F||^2 / ||F||^2 = {ratio:.12f}; |sum A - I| = {np.abs(SA - np.eye(2)).max():.1e}")
    ok &= abs(ratio - 2) < 1e-9

# ---------- (d) one reference index, every sector c: with G_sigma (strings sigma with c letters, holes elsewhere),
# tau_t(rho) = kappa_t G_{rho + t} for holes t of rho (c - 1 letters), T = sum_rho sum_t |tau_t|^2,
# V = sum_sigma (sum_{t hole of sigma} ||w_t||^2) ||G_sigma||^2, Gamma(sigma') = sum_{t in sigma'} w_t (x)_t G_{sigma' - t}:
#   L_r = sum_rho [2 sum_t |tau_t|^2 - |sum_t tau_t|^2] = T + V - ||Gamma||^2 <= ||G||^2   (when sum_s ||w_s||^2 = 1)
def sector_check(n, c, trials=20):
    m = 2; worst_id, worst_ratio = 0, 0
    subsets = list(itertools.combinations(range(n), c))
    for _ in range(trials):
        w = [cz(m) for _ in range(n)]; nw = np.sqrt(sum(np.linalg.norm(x) ** 2 for x in w)); w = [x / nw for x in w]
        G = {S: cz(*([m] * c)) for S in subsets}                              # letters on the copies of S (in increasing order)
        normG = sum(np.linalg.norm(g) ** 2 for g in G.values())
        def contract(S, t):     # kappa_t G_S: contract the letter at copy t (t in S) with conj(w_t)
            ax = list(S).index(t); return np.tensordot(w[t].conj(), G[S], axes=([0], [ax]))
        L = 0; T = 0
        for R in itertools.combinations(range(n), c - 1):                       # rho: excited set R
            holes = [t for t in range(n) if t not in R]
            taus = [contract(tuple(sorted(R + (t,))), t) for t in holes]       # each lives on the letters of R
            L += 2 * sum(np.linalg.norm(x) ** 2 for x in taus) - np.linalg.norm(sum(taus)) ** 2; T += sum(np.linalg.norm(x) ** 2 for x in taus)
        V = sum(sum(np.linalg.norm(w[t]) ** 2 for t in range(n) if t not in S) * np.linalg.norm(G[S]) ** 2 for S in subsets)
        Gam2 = 0
        for Sp in itertools.combinations(range(n), c + 1):
            tot = 0
            for t in Sp:
                S = tuple(x for x in Sp if x != t); out = np.tensordot(w[t], G[S], axes=0)   # axis 0 = slot t
                order = [t] + list(S); perm = np.argsort(order); tot = tot + np.transpose(out, perm)
            Gam2 += np.linalg.norm(tot) ** 2
        worst_id = max(worst_id, abs(L - (T + V - Gam2)) / normG); worst_ratio = max(worst_ratio, L / normG)
    return worst_id, worst_ratio
for n in (3, 4, 5):
    for c in range(2, n):
        e, r = sector_check(n, c)
        print(f"(d) one reference index, n={n} c={c}: identity error {e:.1e}; largest L_r/||G||^2 over random w, G {r:.4f} (bound 1)")
        ok &= e < 1e-12 and r <= 1 + 1e-12

# ---------- (e) Corollary 4 end to end at n = 4, independent of the sector formula: m = min over ALL co-ranges of
# q_4/||C||^2 for single-excitation range planes, from the dense operator [Psi^(x)4(u_r u_r'^*)] (gauge sum A = I).
def psi_n(X, n):
    T = X.reshape([3] * (2 * n))
    for s_ in range(n):
        tr = np.trace(T, axis1=s_, axis2=n + s_)
        full = np.moveaxis(np.einsum('ab,...->...ab', np.eye(3), tr), [2 * n - 2, 2 * n - 1], [s_, n + s_]); T = full - 0.5 * T
    return T.reshape(3 ** n, 3 ** n)
def m_single(U, n):
    D = 3 ** n; st = list(itertools.product(range(3), repeat=n)); idx = {x: i for i, x in enumerate(st)}; u = np.zeros((2, D), complex)
    for r in range(2):
        for s_ in range(n):
            for z in (1, 2):
                t = [0] * n; t[s_] = z; u[r, idx[tuple(t)]] += U[s_][z - 1, r]
    K = np.zeros((2 * D, 2 * D), complex)
    for r in range(2):
        for rp in range(2):
            K[r * D:(r + 1) * D, rp * D:(rp + 1) * D] = psi_n(np.outer(u[r], u[rp].conj()), n)
    return np.linalg.eigvalsh((K + K.conj().T) / 2)[0]
n = 4; vals = [m_single(normalise([cz(2, 2) for _ in range(n)]), n) for _ in range(10)]
def mfun(x): return m_single(normalise([(x[4 * s_:4 * s_ + 4] + 1j * x[16 + 4 * s_:16 + 4 * s_ + 4]).reshape(2, 2) for s_ in range(n)]), n)
for k in range(3):
    res = minimize(mfun, rng.normal(size=32), method='L-BFGS-B', options={'maxiter': 400}); vals.append(res.fun)
print(f"(e) Corollary 4 at n = 4, dense evaluation: min over 10 random and 3 minimised single-excitation planes of the minimum "
      f"over all co-ranges = {min(vals):+.3e}")
ok &= min(vals) > -1e-10
assert ok
print("PASS")
