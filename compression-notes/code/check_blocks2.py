"""Section 10: antisymmetric and code blocks.

   (a) Exact (SymPy) partial-transpose spectrum of the one-sided compression (P_S (x) I) rho^{(x)b} (P_S (x) I) for
       S = Lambda (the antisymmetric subspace of two copies, b = 2) and S = Q (the qutrit [[3,1,2]] code, b = 3).
   (b) The same smallest eigenvalue in floating point for a few larger candidate blocks (Section 10, "Larger blocks").
   (c) Dense checks of the product bound for random operators with ranges in products of local-unitary images of the
       blocks (n = 2, 4, 5)."""
import itertools
import numpy as np
import sympy as sp
from werner import q

R1 = sp.zeros(9, 9)
for i in range(3):
    for j in range(3):
        R1[3 * i + j, 3 * i + j] += 1
        R1[3 * i + j, 3 * j + i] -= sp.Rational(1, 2)


def onesided_pt_exact(V, b):
    """V: list of k orthonormal vectors in (C^3)^{(x)b} (exact). Returns the partial transpose (on Bob) of the one-sided
       compression of rho^{(x)b} to span(V), as a (k 3^b) x (k 3^b) exact matrix."""
    N = 3 ** b; k = len(V); idx = list(itertools.product(range(3), repeat=b))

    def rb(a, bb, a2, b2):
        val = sp.Integer(1)
        for s in range(b):
            val *= R1[3 * a[s] + bb[s], 3 * a2[s] + b2[s]]
            if val == 0:
                return 0
        return val
    M = sp.zeros(k * N, k * N)
    for i, j in itertools.product(range(k), repeat=2):
        for (qi, qq), (pi, pp) in itertools.product(enumerate(idx), repeat=2):
            tot = 0
            for ai, a in enumerate(idx):
                if V[i][ai] == 0:
                    continue
                for a2i, a2 in enumerate(idx):
                    if V[j][a2i] != 0:
                        tot += sp.conjugate(V[i][ai]) * V[j][a2i] * rb(a, qq, a2, pp)
            M[i * N + pi, j * N + qi] = tot          # partial transpose: Bob's q and p exchanged
    return M


def vec(b, terms, norm):
    idx = list(itertools.product(range(3), repeat=b)); v = [sp.Integer(0)] * len(idx)
    for coef, ks in terms:
        v[idx.index(tuple(ks))] += sp.Integer(coef) / norm
    return v


Lam = [vec(2, [(1, (i, j)), (-1, (j, i))], sp.sqrt(2)) for i, j in ((0, 1), (0, 2), (1, 2))]
Q = [vec(3, [(1, (k, (k + j) % 3, (k + 2 * j) % 3)) for k in range(3)], sp.sqrt(3)) for j in range(3)]
print("(a) exact partial-transpose spectra of one-sided compressions:")
for name, V, b in (("Lambda (antisymmetric subspace of two copies)", Lam, 2), ("Q (qutrit [[3,1,2]] code, three copies)", Q, 3)):
    ev = onesided_pt_exact(V, b).eigenvals()
    print(f"  {name}: " + ", ".join(f"{sp.nsimplify(e)} (x{m})" for e, m in sorted(ev.items(), key=lambda t: float(t[0]))), flush=True)

# (b) floating point, larger candidates
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2


def tens(c):
    R = np.array([[1.0]])
    for _ in range(c):
        R = np.kron(R, rho)
    perm = list(range(0, 2 * c, 2)) + list(range(1, 2 * c, 2))
    return R.reshape([3] * (4 * c)).transpose(perm + [p + 2 * c for p in perm]).reshape(9 ** c, 9 ** c)


R2, R3 = tens(2), tens(3)


def ptmin(V, R, N):
    V = np.linalg.qr(V)[0]; k = V.shape[1]
    W = np.kron(V, np.eye(N)); s = W.conj().T @ R @ W
    s = s.reshape(k, N, k, N).transpose(0, 3, 2, 1).reshape(k * N, k * N)
    return np.linalg.eigvalsh((s + s.conj().T) / 2)[0]


e = np.eye(3)
k2 = lambda i, j: np.kron(e[i], e[j])
k3 = lambda i, j, k: np.kron(k2(i, j), e[k])
L2 = np.array([(k2(i, j) - k2(j, i)) / 2 ** .5 for i, j in ((0, 1), (0, 2), (1, 2))]).T
Qn = np.array([sum(k3(k, (k + j) % 3, (k + 2 * j) % 3) for k in range(3)) / 3 ** .5 for j in range(3)]).T
eps = sum(np.linalg.det(np.eye(3)[list(p)]) * k3(*p) for p in itertools.permutations(range(3))) / 6 ** .5
print("\n(b) smallest partial-transpose eigenvalue of other candidate blocks:")
for name, V, R, N in (("Lambda + |00>", np.column_stack([L2, k2(0, 0)]), R2, 9),
                      ("Lambda + sum_i |ii>/sqrt3", np.column_stack([L2, sum(k2(i, i) for i in range(3)) / 3 ** .5]), R2, 9),
                      ("Lambda + (|01> + |10>)/sqrt2", np.column_stack([L2, (k2(0, 1) + k2(1, 0)) / 2 ** .5]), R2, 9),
                      ("Sym^2(C^3)", np.array([k2(i, i) for i in range(3)] + [(k2(i, j) + k2(j, i)) / 2 ** .5 for i, j in ((0, 1), (0, 2), (1, 2))]).T, R2, 9),
                      ("Q + totally antisymmetric vector", np.column_stack([Qn, eps]), R3, 27),
                      ("C^3 (x) |00> (three copies)", np.kron(np.eye(3), k2(0, 0)[:, None]), R3, 27)):
    print(f"  {name}: {ptmin(V, R, N):+.4f}", flush=True)

# (c) dense checks
rng = np.random.default_rng(12)
cr = lambda *s: rng.standard_normal(s) + 1j * rng.standard_normal(s)


def haar(d):
    Qm, Rm = np.linalg.qr(cr(d, d)); return Qm * (np.diag(Rm) / abs(np.diag(Rm)))


print("\n(c) dense q_n / ||C||^2 for random C with range in a product of blocks (any rank, arbitrary co-range):")
for n, blocks, bound in ((2, ["L"], 1 / 4), (4, ["L", "L"], 1 / 16), (5, ["L", "Q"], 1 / 8)):
    worst = np.inf
    for t in range(40 if n < 5 else 15):
        S = np.eye(1)
        for bl in blocks:
            S = np.kron(S, np.kron(haar(3), haar(3)) @ L2 if bl == "L" else np.kron(np.kron(haar(3), haar(3)), haar(3)) @ Qn)
        r = rng.integers(1, S.shape[1] + 1)
        C = S @ cr(S.shape[1], r) @ cr(r, 3 ** n)
        worst = min(worst, q(C, n, 0.5) / np.vdot(C, C).real)
    print(f"  n = {n}, blocks {blocks}: smallest value {worst:.4f} (proved bound {bound:.4f})", flush=True)
