"""Checks for Proposition 8 (single-excitation-range/README.md): small rank-two parts.
(a) the exact identity 2Kt_c = 4 Vac + alpha^* alpha - 2 sum_t Delta_t^{exc}, Delta_t = P_t - Abar_t (x) I_letter,
    P_t = |omega_t><omega_t|, at random complex POVMs (n = 3, 4, 5, every c);
(b) Delta_t has exactly one positive eigenvalue, sqrt(det A_t), with a maximally entangled eigenvector;
(c) on every pattern, 2 sum_{t in S} Delta_t <= (sum_{t in S} sqrt(det A_t) + max_{t in S} sqrt(det A_t)) I
    (weighted central spin bound);
(d) the corollary: lambda_min(2Kt_c) >= min_S [4 lambda_min(A_{H(S)}) - sum_S sqrt(det) - max_S sqrt(det)] at random
    POVMs, and the sector is positive whenever the right side is >= 0;
(e) a POVM covered by Proposition 8 but not by the threshold of Proposition 5 (a heavy set of c copies)."""
import numpy as np, itertools, warnings
from threelevel_certificate import GeneralModel, gauge
warnings.filterwarnings('ignore')
rng = np.random.default_rng(8)
ok = True

def pieces(M, c):
    o = M.ops(c); S_list, m, n = o['S'], o['m'], M.n
    Vac = o['Vac']; AA = sum(o['B']).conj().T @ sum(o['B']); Dsum = np.zeros((2 * m, 2 * m), complex)
    for t in range(n):
        exc = np.diag([1.0 if S_list[j][t] else 0.0 for r in range(2) for j in range(m)])
        Ab = np.zeros((2 * m, 2 * m), complex)
        for j in range(m):
            for r in range(2):
                for rp in range(2): Ab[r * m + j, rp * m + j] = M.A[t].conj()[r, rp]
        Bt = o['B'][t]; Dsum += Bt.conj().T @ Bt - exc @ Ab @ exc
    return o, Vac, AA, Dsum

def delta(A): return np.sqrt(max(np.linalg.det(A).real, 0.0))

# (a)
worst = 0
for n in (3, 4, 5):
    for c in range(1, n + 1):
        M = GeneralModel(gauge([rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(n)]))
        o, Vac, AA, Dsum = pieces(M, c)
        worst = max(worst, np.abs(2 * o['K'] - (4 * Vac + AA - 2 * Dsum)).max())
print(f"(a) 2Kt = 4 Vac + alpha^*alpha - 2 sum Delta^exc: largest error {worst:.1e}"); ok &= worst < 1e-12
# (b)
worst = 0; ment = 0
for _ in range(200):
    U = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)); A = U.conj().T @ U
    om = np.array([U[z, r] for r in range(2) for z in range(2)])      # index (r, z): frame (x) letter
    Dl = np.outer(om, om.conj()) - np.kron(A.conj(), np.eye(2)); w, V = np.linalg.eigh(Dl)
    worst = max(worst, abs(w[-1] - delta(A)), max(w[-2], 0))
    v = V[:, -1].reshape(2, 2); ment = max(ment, np.abs(v @ v.conj().T - np.eye(2) / 2).max())
print(f"(b) top eigenvalue of Delta = sqrt(det A), the others <= 0: largest error {worst:.1e}; "
      f"top eigenvector maximally entangled to {ment:.1e}"); ok &= worst < 1e-10 and ment < 1e-10
# (c), (d)
worst_c = -np.inf; worst_d = np.inf; npos = 0
for trial in range(60):
    n = int(rng.integers(4, 6)); c = int(rng.integers(2, n))
    kind = trial % 3
    if kind == 0: X = [rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(n)]
    elif kind == 1: X = [np.outer(rng.normal(size=2), rng.normal(size=2)) + 0.03 * rng.normal(size=(2, 2)) for _ in range(n)]
    else: X = [np.outer(rng.normal(size=2) + 1j * rng.normal(size=2), rng.normal(size=2)) + 0.1 * rng.normal(size=(2, 2)) for _ in range(n)]
    M = GeneralModel(gauge(X)); o, Vac, AA, Dsum = pieces(M, c); S_list, m = o['S'], o['m']
    dl = [delta(a) for a in M.A]; bound = np.inf
    for Sx in itertools.combinations(range(n), c):
        idx = [r * m + j for r in range(2) for j, s in enumerate(S_list) if frozenset(t for t in range(n) if s[t]) == frozenset(Sx)]
        blk = 2 * Dsum[np.ix_(idx, idx)]
        lhs = np.linalg.eigvalsh((blk + blk.conj().T) / 2)[-1]; rhs = sum(dl[t] for t in Sx) + max(dl[t] for t in Sx)
        worst_c = max(worst_c, lhs - rhs)
        AH = sum((M.A[u] for u in range(n) if u not in Sx), np.zeros((2, 2)))
        bound = min(bound, 4 * np.linalg.eigvalsh(AH)[0] - rhs)
    lm = 2 * np.linalg.eigvalsh(o['K'])[0]; worst_d = min(worst_d, lm - bound); npos += bound >= 0
print(f"(c) largest violation of the central-spin bound: {worst_c:+.1e}")
print(f"(d) lambda_min(2Kt) - (Proposition 8 bound): smallest {worst_d:+.3e} over 60 POVMs ({npos} satisfy the condition)")
ok &= worst_c < 1e-10 and worst_d > -1e-10
# (e) every A_t of rank two, a heavy set of three copies, round hole sets: n = 5, c = 3 (found by a random search)
ang = [3.1005, 1.1748, 0.5977, 0.1559, 2.0937]; wts = [0.7368, 0.6727, 0.7926, 0.8653, 0.5961]; eps = 0.0308
X = [np.sqrt(w) * (np.outer([np.cos(a), np.sin(a)], [np.cos(a), np.sin(a)]) + eps * np.eye(2)) for a, w in zip(ang, wts)]
M = GeneralModel(gauge(X)); A = M.A; n, c = 5, 3; h = n - c
heavy = max(np.linalg.eigvalsh(sum(A[t] for t in Sx))[-1] for Sx in itertools.combinations(range(n), c))
dl = [delta(a) for a in A]; b8 = min(4 * np.linalg.eigvalsh(sum((A[u] for u in range(n) if u not in Sx), np.zeros((2, 2))))[0]
                                     - sum(dl[t] for t in Sx) - max(dl[t] for t in Sx) for Sx in itertools.combinations(range(n), c))
lm = 2 * np.linalg.eigvalsh(M.ops(c)['K'])[0]
print(f"(e) n=5 c=3, min sqrt(det A_t) = {min(dl):.3f}: max_S lambda_max(A_S) = {heavy:.3f} > Proposition 5 threshold {2 * h / (3 * h - 1):.3f}; "
      f"Proposition 8 bound {b8:+.3e} >= 0; lambda_min(2Kt) {lm:+.3e}")
ok &= heavy > 2 * h / (3 * h - 1) and b8 >= 0 and lm >= b8 - 1e-10
assert ok
print("PASS")
