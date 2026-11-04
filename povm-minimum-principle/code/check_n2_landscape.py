"""n = 2: the POVM is (A, I - A), so m = phi(a, b) with (a, b) the eigenvalues of A (m is unitarily invariant).
Facts checked: phi = 0 at the four corners and 1/4 at the centre (the W-plane); phi is not concave (positive Hessian
eigenvalues near a corner); phi is not quasi-concave (an explicit segment dips below both ends); local minimisation
from 40 interior starts always ends at a corner (no interior local minimum)."""
import numpy as np, warnings
from scipy.optimize import minimize
from povm_form import m_povm
warnings.filterwarnings("ignore")
phi = lambda a, b: m_povm([np.diag([a, b]).astype(complex), np.diag([1 - a, 1 - b]).astype(complex)], 2)
corners = [phi(0, 0), phi(0, 1), phi(1, 0), phi(1, 1)]; centre = phi(0.5, 0.5)
N = 41; g = np.linspace(0, 1, N); P = np.array([[phi(a, b) for b in g] for a in g]); h = g[1] - g[0]
Haa = (P[2:, 1:-1] - 2 * P[1:-1, 1:-1] + P[:-2, 1:-1]) / h ** 2; Hbb = (P[1:-1, 2:] - 2 * P[1:-1, 1:-1] + P[1:-1, :-2]) / h ** 2
Hab = (P[2:, 2:] - P[2:, :-2] - P[:-2, 2:] + P[:-2, :-2]) / (4 * h * h)
lmax = 0.5 * (Haa + Hbb) + np.sqrt(0.25 * (Haa - Hbb) ** 2 + Hab ** 2)
p, r = np.array([0.1313, 0.1234]), np.array([0.2617, 0.0501]); ts = np.linspace(0, 1, 21)
prof = np.array([phi(*((1 - t) * p + t * r)) for t in ts]); dip = min(prof[0], prof[-1]) - prof.min()
rng = np.random.default_rng(7); ends = []
for k in range(40):
    res = minimize(lambda x: phi(*np.clip(x, 0, 1)) + 100 * (np.sum(np.minimum(x, 0) ** 2) + np.sum(np.maximum(x - 1, 0) ** 2)),
                   rng.uniform(0.02, 0.98, 2), method='Nelder-Mead', options={'xatol': 1e-7, 'fatol': 1e-10, 'maxiter': 2000})
    ends.append(np.clip(res.x, 0, 1))
at_corner = sum(1 for x in ends if np.all(np.minimum(x, 1 - x) < 1e-3))
print(f"corners {np.round(corners, 6)}, centre {centre:.6f}")
print(f"largest Hessian eigenvalue on the interior grid {lmax.max():+.3f} (positive: not concave); fraction of grid points with a positive one {(lmax > 1e-6).mean():.2f}")
print(f"segment {p} -> {r}: ends {prof[0]:.5f}, {prof[-1]:.5f}, interior minimum {prof.min():.5f}: dips {dip:+.2e} below both ends (not quasi-concave)")
print(f"local minimisations ending at a corner: {at_corner} of 40")
assert max(abs(c) for c in corners) < 1e-12 and abs(centre - 0.25) < 1e-12 and lmax.max() > 0 and dip > 1e-4 and at_corner == 40
print("PASS")
