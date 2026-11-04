"""Balanced-slice margins for Conjecture SC.  By scale invariance only the sign of nu_s is meaningful; balancing C at copy s
(C -> A C A^*, A in SL_3 at copy s, minimising ||A C A^*||, so that Tr_{-s}(CC^* + C^*C) is proportional to I) gives
a scale-free margin.  Samples random tangent operators in several regimes and reports the largest balanced nu at copy 0."""
import numpy as np, sys
from scipy.optimize import minimize
from scipy.linalg import expm
from common import tangent_range, squeeze_forms, pencil_margin, Bt, d
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
B8 = Bt[[0, 1, 3, 4, 5, 6, 7, 8]]


def act(A, X, n):
    D = d ** n; T = X.reshape(3, D // 3, 3, D // 3); return np.einsum('ab,bicj,cd->aidj', A, T, A.conj().T).reshape(D, D)


def balance(X, n):
    def hm(h):
        H = np.einsum('k,kab->ab', h, B8); return H - np.trace(H) / 3 * np.eye(3)
    r = minimize(lambda h: np.log(np.linalg.norm(act(expm(hm(h)), X, n)) ** 2), np.zeros(8), method='BFGS', options={'gtol': 1e-10, 'maxiter': 3000})
    Y = act(expm(hm(r.x)), X, n); return Y / np.linalg.norm(Y)


cz = lambda *s: rng.normal(size=s) + 1j * rng.normal(size=s)
overall = -np.inf
for n, samples in ((2, 12), (3, 10), (4, 6)):
    D = d ** n
    for regime in ("generic", "near-thin U_0", "a = 0", "range = co-range"):
        vals = []
        for k in range(samples):
            U = [cz(2, 2) for _ in range(n)]; a = cz(2) * (0 if regime == "a = 0" else 1)
            if regime == "near-thin U_0":
                Uu, S, Vh = np.linalg.svd(U[0]); U[0] = Uu @ np.diag([S[0], 1e-3 * S[0]]) @ Vh
            u = tangent_range(U, a, n)
            v = u.copy() if regime == "range = co-range" else cz(2, D)
            X = np.outer(u[0], v[0].conj()) + np.outer(u[1], v[1].conj()); X /= np.linalg.norm(X)
            Xb = balance(X, n); vals.append(pencil_margin(*squeeze_forms(Xb, n, 0))[0])
        overall = max(overall, max(vals))
        print(f"n={n} {regime:17s}: samples {samples}, largest balanced nu {max(vals):+.4f}, median {np.median(vals):+.4f}", flush=True)
print(f"RESULT: largest balanced nu over all samples {overall:+.4f}")
assert overall < 0
print("PASS")
