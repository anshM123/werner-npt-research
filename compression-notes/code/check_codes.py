"""Section 8: the endpoint form on the code space of a pure qutrit code.

   (a) Exact rationals A_n, B_n for n = 2d - 1 (every set of copies is either smaller than d or has a complement
       smaller than d):  A_n = sum_{l<d} C(n,l) (-1/6)^l,  B_n = sum_{l>=d} C(n,l) (-1/2)^l 3^{l-n},
       so that q_n(1/2, V x V^*) = A_n ||x||^2 + B_n |Tr x|^2 on a pure [[n, 1, d]]_3 code.
   (b) The cyclic five-qutrit code (stabilisers: cyclic shifts of X Z Z^dag X^dag I): dimension 3, purity of the
       stabiliser group (no nontrivial element of weight < 3), and the formula against the dense q_5."""
import itertools
from fractions import Fraction as Fr
from math import comb
import numpy as np
from werner import q

print("(a) exact coefficients:")
for n in range(3, 22, 2):
    d = (n + 1) // 2
    A = sum(comb(n, l) * Fr(-1, 2) ** l * Fr(1, 3 ** l) for l in range(d))
    B = sum(comb(n, l) * Fr(-1, 2) ** l * Fr(1, 3 ** (n - l)) for l in range(d, n + 1))
    print(f"  n = {n:2d}: A_n = {str(A):>22s} ({float(A):.6f}), B_n = {float(B):+.6f}, "
          f"A_n + 3 min(B_n, 0) = {float(A + 3 * min(B, 0)):+.6f}")

w = np.exp(2j * np.pi / 3)
X = np.roll(np.eye(3), 1, axis=0); Z = np.diag([1, w, w ** 2])
paulis = {(a, b): np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b) for a in range(3) for b in range(3)}
base = [(1, 0), (0, 1), (0, 2), (2, 0), (0, 0)]          # X, Z, Z^dag, X^dag, I  (as exponents of X^a Z^b, up to phase)
gens = [base[-k:] + base[:-k] for k in range(4)]


def op(lst):
    M = np.eye(1)
    for ab in lst:
        M = np.kron(M, paulis[ab])
    return M


G = [op(g) for g in gens]
print("\n(b) five-qutrit code: generators commute:", all(np.allclose(a @ b, b @ a) for a in G for b in G))
# purity: every element of the stabiliser group other than the identity has weight >= 3
weights = []
for e in itertools.product(range(3), repeat=4):
    if not any(e):
        continue
    exps = [((sum(e[k] * gens[k][s][0] for k in range(4))) % 3, (sum(e[k] * gens[k][s][1] for k in range(4))) % 3)
            for s in range(5)]
    weights.append(sum(1 for ab in exps if ab != (0, 0)))
print("  smallest weight of a nontrivial stabiliser element:", min(weights), "(purity needs >= 3)")
P = np.eye(243, dtype=complex)
for g in G:
    P = P @ sum(np.linalg.matrix_power(g, k) for k in range(3)) / 3
ev, U = np.linalg.eigh((P + P.conj().T) / 2)
V = U[:, ev > 0.5]
print("  code dimension:", V.shape[1])
# purity up to distance 3 on the whole code space: Tr_{S^c}(v_i v_j^*) = delta_ij I / 3^{|S|} for every |S| <= 2
dev = 0.0
for k in (1, 2):
    for S in itertools.combinations(range(5), k):
        rest = [s for s in range(5) if s not in S]
        for i in range(3):
            for j in range(3):
                T = np.multiply.outer(V[:, i].reshape([3] * 5), V[:, j].conj().reshape([3] * 5))
                T = np.transpose(T, list(S) + rest + [5 + s for s in S] + [5 + s for s in rest])
                M = T.reshape(3 ** k, 3 ** (5 - k), 3 ** k, 3 ** (5 - k))
                red = np.einsum('aibi->ab', M)
                dev = max(dev, np.abs(red - (i == j) * np.eye(3 ** k) / 3 ** k).max())
print(f"  max |Tr_(S^c)(v_i v_j^*) - delta_ij I/3^|S||| over |S| <= 2 and all i, j: {dev:.2e}")
A5, B5 = Fr(4, 9), Fr(-19, 288)
rng = np.random.default_rng(1)
for trial, r in itertools.product(range(2), (2, 3)):
    x = (rng.standard_normal((3, r)) + 1j * rng.standard_normal((3, r))) @ (rng.standard_normal((r, 3)) + 1j * rng.standard_normal((r, 3)))
    C = V @ x @ V.conj().T
    print(f"  rank-{r} x: dense q_5 = {q(C, 5, 0.5):+.10f}, formula {float(A5) * np.vdot(x, x).real + float(B5) * abs(np.trace(x)) ** 2:+.10f}")
print(f"  hence q_5 >= (A_5 + 3 B_5) ||C||^2 = {float(A5 + 3 * B5):.6f} ||C||^2 for every operator on the code space")

# (c) one-sided compression: range in the code space, co-range arbitrary.  Build (V^* (x) I)(rho^G)^{(x)5}(V (x) I) on
#     code (x) B^5 (729 x 729), rho^G = I - (3/2) P_Phi per copy, and print its spectrum range.
Phi = np.eye(3).reshape(-1) / np.sqrt(3)
rG4 = (np.eye(9) - 1.5 * np.outer(Phi, Phi)).reshape(3, 3, 3, 3)


def apply_rG5(T):
    for s in range(5):
        T = np.moveaxis(np.tensordot(rG4, T, axes=([2, 3], [s, 5 + s])), [0, 1], [s, 5 + s])
    return T


Vt = V.reshape([3] * 5 + [3])
M = np.zeros((729, 729), complex)
for c in range(3):
    for b in range(243):
        xi = np.zeros((3, 243), complex); xi[c, b] = 1
        T = apply_rG5(np.tensordot(Vt, xi.reshape([3] + [3] * 5), axes=([5], [0])))
        M[:, c * 243 + b] = np.tensordot(Vt.conj(), T, axes=(list(range(5)), list(range(5)))).reshape(-1)
ev = np.linalg.eigvalsh((M + M.conj().T) / 2)
print(f"\n(c) one-sided compression to the code space: partial-transpose spectrum in [{ev[0]:.10f}, {ev[-1]:.10f}];"
      f" 71/288 = {71 / 288:.10f}, 4/9 = {4 / 9:.10f}")
