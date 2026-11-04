"""Check of Lemma 2 (thin-copy removal inside the tangent class).  If C has range in T_e and its range is thin at copy s
(range block U_s of rank one, local support V = span{|0>, y}), then with the qubit 3-design (six Pauli eigenstates psi_k
of V) and y^perp the remaining direction,
   q_n(C) = 1/2 sum_k q_{n-1}(<psi_k|_s C |psi_k^perp>_s) + sum_{a in {0, y}} q_{n-1}(<a|_s C |y^perp>_s),
and every contraction has rank <= 2 and range in T_{e'} for e' = e with copy s removed."""
import numpy as np, itertools
from common import strings, q_direct, tangent_range, d
rng = np.random.default_rng(7)


def contract(X, n, s, x, y):
    """<x|_s X |y>_s (x on the range side, y on the co-range side)."""
    T = X.reshape([d] * (2 * n))
    T = np.tensordot(x.conj(), T, axes=([0], [s]))           # removes row index s
    T = np.tensordot(T, y, axes=([n - 1 + s], [0]))          # column index s (shifted by one)
    m = n - 1; return T.reshape(d ** m, d ** m)


def in_tangent(Y, m):
    """fraction of the range of Y outside T_{|0..0>} on m copies (should be 0)."""
    st = strings(m); allowed = np.array([sum(c != 0 for c in s_) <= 1 for s_ in st])
    Uq, S, _ = np.linalg.svd(Y); r = (S > 1e-10 * max(S[0], 1e-300)).sum()
    P = Uq[:, :r] @ Uq[:, :r].conj().T
    return np.real(np.trace(P[np.ix_(~allowed, ~allowed)])), r


worst = 0.0
for n in (2, 3, 4):
    for trial in range(3):
        U = [rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(n)]
        s = trial % n
        yv = rng.normal(size=2) + 1j * rng.normal(size=2); bv = rng.normal(size=2) + 1j * rng.normal(size=2)
        U[s] = np.outer(yv, bv.conj())                        # rank one: thin at copy s
        a = rng.normal(size=2) + 1j * rng.normal(size=2)
        u = tangent_range(U, a, n); D = d ** n
        v = rng.normal(size=(2, D)) + 1j * rng.normal(size=(2, D))
        X = np.outer(u[0], v[0].conj()) + np.outer(u[1], v[1].conj())
        y = np.zeros(3, complex); y[1:] = yv / np.linalg.norm(yv)
        yp = np.zeros(3, complex); yp[1:] = np.array([-np.conj(y[2]), np.conj(y[1])])
        e0 = np.array([1, 0, 0], complex)
        # six Pauli eigenstates of V = span{e0, y}
        basisV = np.column_stack([e0, y])
        pauli = [np.array([1, 0]), np.array([0, 1]), np.array([1, 1]) / np.sqrt(2), np.array([1, -1]) / np.sqrt(2),
                 np.array([1, 1j]) / np.sqrt(2), np.array([1, -1j]) / np.sqrt(2)]
        total = 0.0; bad = 0.0
        for c in pauli:
            psi = basisV @ c; perp = basisV @ np.array([-np.conj(c[1]), np.conj(c[0])])
            Ck = contract(X, n, s, psi, perp)
            total += 0.5 * q_direct(Ck, n - 1); frac, r = in_tangent(Ck, n - 1); bad = max(bad, frac, r - 2)
        for avec in (e0, y):
            Ck = contract(X, n, s, avec, yp)
            total += q_direct(Ck, n - 1); frac, r = in_tangent(Ck, n - 1); bad = max(bad, frac, r - 2)
        qn = q_direct(X, n); err = abs(total - qn) / max(1.0, abs(qn)); worst = max(worst, err)
        print(f"n={n} thin copy {s}: q_n = {qn:.10f}, sum of (n-1)-copy terms = {total:.10f}, rel. error {err:.1e}, "
              f"max(range outside T, rank-2) = {bad:.1e}")
print(f"RESULT: maximal relative error {worst:.1e}")
assert worst < 1e-10
print("PASS")
