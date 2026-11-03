"""Theorem 3: the block condition for a pair of qutrit copies with anchor Omega_3 = sum_a |aa>.
M = I - |Omega><Omega|/2 per copy (row (x) col), M_2 = M (x) M regrouped as (row1 row2) (x) (col1 col2), V = Omega_3 (x) I,
D = V^* M_2 V, R = M_2 V D^{-1} V^* M_2, Gamma = transpose on the row factor. Needs D > 0 and (M_2 +- R)^Gamma >= 0.
Exact: 352 (M_2 +- R)^Gamma are integer matrices; PSD is read off the signs of their characteristic polynomials
(python-flint)."""
import itertools
import sympy as sp

M1 = sp.eye(9) - sp.Rational(1, 2) * sp.Matrix(9, 9, lambda i, j: 1 if (i % 4 == 0 and j % 4 == 0) else 0)
full = sp.kronecker_product(M1, M1)                       # order (row1, col1, row2, col2)
perm = [None] * 81
for r1, c1, r2, c2 in itertools.product(range(3), repeat=4):
    perm[((r1 * 3 + r2) * 3 + c1) * 3 + c2] = ((r1 * 3 + c1) * 3 + r2) * 3 + c2
M2 = sp.Matrix(81, 81, lambda i, j: full[perm[i], perm[j]])
V = sp.Matrix(81, 9, lambda i, j: 1 if (i % 9 == j and i // 9 in (0, 4, 8)) else 0)
D = V.T * M2 * V
P_omega = sp.Matrix(9, 9, lambda i, j: 1 if (i % 4 == 0 and j % 4 == 0) else 0)
assert D == 2 * sp.eye(9) + P_omega / 4
print("D = 2 I + |Omega><Omega|/4 > 0")
R = M2 * V * D.inv() * V.T * M2


def gamma(X):
    return sp.Matrix(81, 81, lambda i, j: X[(j // 9) * 9 + i % 9, (i // 9) * 9 + j % 9])


A, B = gamma(M2), gamma(R)
try:
    import flint
except ImportError:
    flint = None
for sign, name in ((1, "352 (M_2 + R)^Gamma"), (-1, "352 (M_2 - R)^Gamma")):
    H = 352 * (A + sign * B)
    assert H == H.T and all(x.is_integer for x in H)
    if flint is None:
        import numpy as np
        ev = np.linalg.eigvalsh(np.array(H, dtype=float))
        print(f"{name}: smallest eigenvalue {ev.min():.3e} (floating point; install python-flint for the exact check)")
        continue
    cp = [int(c) for c in flint.fmpz_mat([[int(x) for x in H.row(i)] for i in range(81)]).charpoly().coeffs()]
    deg = len(cp) - 1                                    # det(tI - H) = sum_k cp[k] t^k
    psd = all(c == 0 or (c > 0) == ((deg - k) % 2 == 0) for k, c in enumerate(cp))
    kernel = next(k for k, c in enumerate(cp) if c != 0)
    print(f"{name}: integer, symmetric; characteristic polynomial has alternating signs (so PSD): {psd};"
          f" kernel dimension {kernel}")
