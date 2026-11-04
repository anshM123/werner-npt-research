"""Quasi-concavity of each sector function mu_c along excitation-squeeze geodesics U_s(tau) = exp(tau K_s) U_s
(re-gauged); counts curves with an interior value below both ends, per sector."""
import numpy as np, sys, warnings
warnings.filterwarnings("ignore")

from povm_form import sector_mu
from scipy.linalg import expm
n, trials, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); rng = np.random.default_rng(seed)
def povm_from_U(U):
    S = sum(u.conj().T @ u for u in U); w, V = np.linalg.eigh(S); T = V @ np.diag(w ** -0.5) @ V.conj().T
    return [T @ u.conj().T @ u @ T for u in U]
herm = lambda: (lambda g: (g + g.conj().T) / 2)(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
taus = np.linspace(-3, 3, 61); dips = {c: 0 for c in range(1, n + 1)}; worst = {c: -np.inf for c in range(1, n + 1)}
for k in range(trials):
    U = [rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(n)]; K = [herm() for _ in range(n)]
    As = [povm_from_U([expm(t * k_) @ u for k_, u in zip(K, U)]) for t in taus]
    for c in range(1, n + 1):
        prof = np.array([sector_mu(A, n, c) for A in As]); v = min(prof[0], prof[-1]) - prof[1:-1].min()
        worst[c] = max(worst[c], v); dips[c] += v > 1e-9
print(f"n={n}: curves with a dip, per sector: " + ", ".join(f"c={c}: {dips[c]}/{trials} (largest {worst[c]:+.1e})" for c in range(1, n + 1)))
