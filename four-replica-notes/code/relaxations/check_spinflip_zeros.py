"""Spin-flip zeros on random qubit subspaces: for n = 2..5, Q_i = random 2-dim subspaces of C^3 (complex),
u a random vector of Q_1 (x) ... (x) Q_n, v = J^(x)n conj(u) (J = sigma_y in the Q_i bases).
n odd: q_n(P_u + P_v) = 0 for every u.  n even: q_n(P_u - P_v) = 0 once <u, J^(x)n conj u> = 0 (imposed by a root of
a quadratic along a second random direction).  usage: python check_spinflip_zeros.py"""
import numpy as np
from tm1_analyze import qform

SY = np.array([[0, -1j], [1j, 0]])
rng = np.random.default_rng(5)
for n in (2, 3, 4, 5):
    Y = SY
    for _ in range(n - 1): Y = np.kron(Y, SY)
    B = None
    for i in range(n):
        Q, _ = np.linalg.qr(rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))); R = Q[:, :2]
        B = R if B is None else np.kron(B, R)
    out = []
    for trial in range(3):
        u = rng.standard_normal(2 ** n) + 1j * rng.standard_normal(2 ** n); u /= np.linalg.norm(u)
        if n % 2 == 0:
            a = u; bb = rng.standard_normal(2 ** n) + 1j * rng.standard_normal(2 ** n)
            lam = np.roots([bb @ Y @ bb, 2 * (a @ Y @ bb), a @ Y @ a])[0]; u = a + lam * bb; u /= np.linalg.norm(u)
        v = Y @ u.conj()
        U, V = B @ u, B @ v
        t = 1.0 if n % 2 == 1 else -1.0
        C = np.outer(U, U.conj()) + t * np.outer(V, V.conj())
        out.append((abs(np.vdot(U, V)), qform(C, n)))
    print(f"n = {n}, C = P_u {'+' if n % 2 else '-'} P_v: " + '; '.join(f'|<u,v>| {a:.1e}, q_n {b:+.1e}' for a, b in out))
