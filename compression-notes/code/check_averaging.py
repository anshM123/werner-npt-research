"""Section 2: averages of rank-one contractions that do not depend on C.

   For an operator C on n copies and a copy s put S = sum_{ij} q_{n-1}(X_ij) and T = Tr_s C, so that
   q_n(t, C) = S - t q_{n-1}(T). The script checks
       E_u          q_{n-1}(<u|_s C|u>_s) = (S + q(T)) / (d (d + 1))                 (u Haar),
       E_{u perp v} q_{n-1}(<u|_s C|v>_s) = (d S - q(T)) / ((d + 1) d (d - 1))      (orthonormal pairs).
   Both integrands have degree (2, 2) in the Haar vectors, so the averages are computed exactly over the twelve
   vectors of the four mutually unbiased bases of C^3 (a 2-design), and over ordered pairs inside each basis."""
import numpy as np
from werner import d, q, contract

rng = np.random.default_rng(2)
w = np.exp(2j * np.pi / 3)
bases = [np.eye(3, dtype=complex)]
for k in range(3):
    bases.append(np.array([[w ** (k * j * j + m * j) for j in range(3)] for m in range(3)]).T / np.sqrt(3))
for n in (2, 3):
    D = d ** n
    C = (rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2))) @ \
        (rng.standard_normal((2, D)) + 1j * rng.standard_normal((2, D)))
    for t in (0.3, 0.5):
        X = [[contract(C, n, 0, np.eye(d)[i], np.eye(d)[j]) for j in range(d)] for i in range(d)]
        S = sum(q(X[i][j], n - 1, t) for i in range(d) for j in range(d))
        qT = q(sum(X[a][a] for a in range(d)), n - 1, t)
        diag = np.mean([q(contract(C, n, 0, B[:, i], B[:, i]), n - 1, t) for B in bases for i in range(3)])
        off = np.mean([q(contract(C, n, 0, B[:, i], B[:, j]), n - 1, t) for B in bases for i in range(3) for j in range(3) if i != j])
        print(f"n={n} t={t}: E_u q(<u|C|u>) = {diag:.10f}  formula {(S + qT) / (d * (d + 1)):.10f};   "
              f"E_(u perp v) q(<u|C|v>) = {off:.10f}  formula {(d * S - qT) / ((d + 1) * d * (d - 1)):.10f}")
    print(f"  q_n(1/2, C) = {q(C, n, 0.5):.8f};  9 E_off - (3/2) E_diag at t = 1/2: {9 * off - 1.5 * diag:.8f}")
