"""Numerical evidence for Conjecture 6 (single-excitation-range/README.md): the combined certificate exists.
(a) solves the certificate SDP at structured POVMs for n = 5, 6 and 3 <= c <= n - 1;
(b) re-checks every LMI of the returned solution with numpy eigenvalues (independent of the solver);
(c) checks lambda_min(Q_c) >= 0 from the explicit sector matrix at the same POVMs;
(d) shows that the purely local certificate (theta = 0) fails at c = n - 1, n = 6, one heavy copy.
"""
import numpy as np
from combined_certificate import gauge, sector_matrix, certificate
rng = np.random.default_rng(2026)

def heavy(N, d): return [np.array([[np.sqrt(1 - d), 0], [0, 1.0]], complex)] + [np.array([[np.sqrt(d / N), 0], [0, 0]], complex) for _ in range(N)]
def cz(*s): return rng.normal(size=s) + 1j * rng.normal(size=s)
def samples(n):
    out = [("uniform W-plane", gauge([np.eye(2, dtype=complex)] * n)), ("one heavy copy, d = 0.05", heavy(n - 1, 0.05)),
           ("random", gauge([cz(2, 2) for _ in range(n)]))]
    X = [np.outer(cz(2), cz(2)) + 0.05 * cz(2, 2) for _ in range(n)]; out.append(("near rank-one frame", gauge(X)))
    X = [np.outer(cz(2), cz(2)) + 0.05 * cz(2, 2) for _ in range(n)]; X[0] *= 0.02; X[1] *= 0.02; out.append(("two nearly unused copies", gauge(X)))
    X = [np.eye(2) + 0.1 * cz(2, 2), np.eye(2) + 0.1 * cz(2, 2)] + [0.1 * cz(2, 2) for _ in range(n - 2)]; out.append(("two heavy copies", gauge(X)))
    return out

ok = True
for n, c in [(5, 3), (5, 4), (6, 3), (6, 4), (6, 5)]:
    for name, U in samples(n):
        lam, worst, thmin = certificate(U, n, c); mu = np.linalg.eigvalsh(sector_matrix(U, n, c))[0]
        print(f"n={n} c={c} {name}: certificate slack {lam:+.3e} (re-checked min eigenvalue {worst:+.3e}, min theta {thmin:+.1e}); "
              f"lambda_min Q_c {mu:+.3e}", flush=True)
        ok &= worst > -1e-7 and thmin > -1e-9 and mu > -1e-9
lam, worst, _ = certificate(heavy(5, 0.2), 6, 5, local_only=True)
print(f"(d) n=6 c=5 one heavy copy (d = 0.2), local certificate only (theta = 0): slack {lam:+.3e} (< 0: infeasible)")
ok &= lam < -1e-4
assert ok
print("PASS")
