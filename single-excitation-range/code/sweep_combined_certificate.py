"""Sweep for Conjecture 6: the combined certificate at structured hard POVMs (default n = 5, c = 3, 60 samples).
Kinds: rank-one frames with noise, nearly unused copies, nearly parallel directions, one heavy copy, two heavy copies,
near W-planes.  Prints the worst slack (re-checked eigenvalue) and lambda_min of the sector form at that POVM."""
import numpy as np, sys
from combined_certificate import gauge, sector_matrix, certificate
n, c = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (5, 3)
N = int(sys.argv[3]) if len(sys.argv) > 3 else 60
rng = np.random.default_rng(11)
def cz(*s): return rng.normal(size=s) + 1j * rng.normal(size=s)
def sample(kind):
    if kind == 'rank-one frame + noise':
        X = [np.outer(cz(2), cz(2)) + rng.choice([0.0, 0.01, 0.05, 0.2]) * cz(2, 2) for _ in range(n)]
    elif kind == 'nearly unused copies':
        X = [np.outer(cz(2), cz(2)) + 0.05 * cz(2, 2) for _ in range(n)]
        for i in rng.choice(n, rng.integers(1, n - 2), replace=False): X[i] = rng.choice([0.003, 0.02, 0.1]) * X[i]
    elif kind == 'nearly parallel':
        v = cz(2); X = [np.outer(cz(2), v + rng.choice([0.01, 0.1, 0.3]) * cz(2)) + 0.03 * cz(2, 2) for _ in range(n)]
    elif kind == 'one heavy copy':
        d = rng.choice([0.02, 0.1, 0.3])
        X = [np.diag([np.sqrt(1 - d), 1.0]) + 0.05 * cz(2, 2)] + [np.sqrt(d / (n - 1)) * (np.outer(cz(2), cz(2)) + 0.1 * cz(2, 2)) for _ in range(n - 1)]
    elif kind == 'two heavy copies':
        X = [np.eye(2) + 0.1 * cz(2, 2), np.eye(2) + 0.1 * cz(2, 2)] + [0.1 * cz(2, 2) for _ in range(n - 2)]
    else:
        X = [np.eye(2) + rng.choice([0.1, 0.5]) * cz(2, 2) for _ in range(n)]
    return gauge(X)
kinds = ['rank-one frame + noise', 'nearly unused copies', 'nearly parallel', 'one heavy copy', 'two heavy copies', 'near W-plane']
worst = (np.inf, None, None)
for i in range(N):
    kind = kinds[i % len(kinds)]; U = sample(kind); lam, chk, thmin = certificate(U, n, c)
    val = min(lam, chk) if lam is not None else -np.inf
    if val < worst[0]: worst = (val, kind, np.linalg.eigvalsh(sector_matrix(U, n, c))[0])
    if (i + 1) % 12 == 0: print(f"n={n} c={c}: {i + 1} samples, worst slack {worst[0]:+.3e} ({worst[1]}; lambda_min Q_c there {worst[2]:+.3e})", flush=True)
assert worst[0] > -1e-7
print("PASS")
