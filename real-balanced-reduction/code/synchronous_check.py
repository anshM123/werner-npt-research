"""Checks for "Synchronous balancing" (README).  Exact unless stated.
(1) the twelve matrices O = S^r D_e form a group, (1/12) sum O^T H O = Tr(H) I/3 for a symbolic complex H, and with the
    same O on every copy, (1/12) sum (O^{(x)2})^T (H (x) I) O^{(x)2} = Tr(H)/3 I for a symbolic H on the first of two copies;
(2) the simplex code c(j)_v = j.v mod 2 (j = 0..11, v in F_2^4 \\ 0): 15 positions, all pairwise distances 8;
(3) the four-index flag kernel of the Werner form, product over the 15 positions of the one-position kernel
    (alpha if all four bits agree; beta if bits (j,k) agree, bits (l,m) agree and the pairs differ; chi if j = l != k = m;
    0 otherwise), summed over all 12^4 label quadruples by class: it is nonzero only for j = k = l = m (alpha^15),
    j = k != l = m (alpha^7 beta^8) and j = l != k = m (alpha^7 chi^8), 12, 132 and 132 times;
(4) the same classification for the kernel of the map P_a (A, B as above, K for both remaining two-plus-two patterns):
    classes (A^15) x 12, (A^7 B^8) x 132, (A^7 K^8) x 264, (A^3 B^4 K^8) x 936, and zero otherwise;
(5) the balance of the image at a flag site: for j != k the flags F_j, F_k stay orthogonal after one flag qutrit is
    traced (every pair of codewords differs in at least 7 further positions);
(6) floating point, for the vanishing-error Corollary: q_r(t, g_r g_r^*) = (2/3)(1 + (-t)^r) + (1 - t)^r/3 for
    g_r = (|0^r> + |1^r> + |2^r>)/sqrt3, r = 2..5; and for the balanced Latin-square plane on three copies, P (x) g_2 g_2^*
    is a balanced rank-two projector with q multiplicative.
usage: python synchronous_check.py"""
import itertools
import sympy as sp

ok = True
S = sp.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
signs = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
Os = [S ** r * sp.diag(*e) for r in range(3) for e in signs]
keys = {tuple(O): i for i, O in enumerate(Os)}
closed = all(tuple(A * B) in keys for A in Os for B in Os)
H = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'h{i}{j}'))
one = sp.simplify(sum((O.T * H * O for O in Os), sp.zeros(3, 3)) / 12 - H.trace() / 3 * sp.eye(3)) == sp.zeros(3, 3)
from sympy import kronecker_product as kp
two = sp.simplify(sum((kp(O, O).T * kp(H, sp.eye(3)) * kp(O, O) for O in Os), sp.zeros(9, 9)) / 12
                  - H.trace() / 3 * sp.eye(9)) == sp.zeros(9, 9)
ok &= closed and one and two
print(f"(1) group closure: {closed}; one-copy twirl exact: {one}; synchronous twirl on two copies exact: {two}")

vs = [v for v in itertools.product(range(2), repeat=4) if any(v)]
bits = lambda j: [(j >> b) & 1 for b in range(4)]
code = [[sum(x * y for x, y in zip(bits(j), v)) % 2 for v in vs] for j in range(12)]
dists = {sum(a != b for a, b in zip(code[j], code[k])) for j in range(12) for k in range(12) if j != k}
ok &= len(vs) == 15 and dists == {8}
print(f"(2) simplex code: {len(vs)} positions, pairwise distances {sorted(dists)}")

al, be, ch, A_, B_, K_ = sp.symbols('alpha beta chi A B K')


def werner_kernel(a, b, c, d):
    if a == b == c == d:
        return al
    if a == b and c == d:
        return be
    if a == c and b == d:
        return ch
    return 0


def map_kernel(a, b, c, d):
    if a == b == c == d:
        return A_
    if a == b and c == d:
        return B_
    if (a + b + c + d) % 2 == 0:
        return K_
    return 0


for name, ker, want in [('Werner', werner_kernel, 12 * al ** 15 + 132 * al ** 7 * be ** 8 + 132 * al ** 7 * ch ** 8),
                        ('map P_a', map_kernel, 12 * A_ ** 15 + 132 * A_ ** 7 * B_ ** 8 + 264 * A_ ** 7 * K_ ** 8
                         + 936 * A_ ** 3 * B_ ** 4 * K_ ** 8)]:
    tot = 0
    classes = {}
    for j, k, l, m in itertools.product(range(12), repeat=4):
        val = 1
        for p in range(15):
            val *= ker(code[j][p], code[k][p], code[l][p], code[m][p])
            if val == 0:
                break
        if val != 0:
            cls = ('all equal' if j == k == l == m else 'j=k!=l=m' if (j == k and l == m) else
                   'j=l!=k=m' if (j == l and k == m) else 'j=m!=k=l' if (j == m and k == l) else 'four distinct')
            classes[cls] = classes.get(cls, 0) + 1
        tot += val
    good = sp.expand(tot - want) == 0
    ok &= good
    print(f"(3/4) {name} flag kernel summed over 12^4 quadruples matches the stated class sum: {good}; nonzero classes {classes}")

# (5) flag orthogonality after tracing one flag qutrit: distinct codewords differ in >= 8 positions, so >= 7 untouched
ok &= min(dists) - 1 >= 1
print(f"(5) distinct labels differ at >= {min(dists)} flag pairs, so >= {min(dists) - 1} remain after one flag qutrit is traced")
import numpy as np


def qnum(C, n, t):
    tot = 0.0
    for r in range(n + 1):
        for K in itertools.combinations(range(n), r):
            X, rem = C.reshape([3] * 2 * n), list(range(n))
            for s_ in sorted(K, reverse=True):
                i = rem.index(s_); X = np.trace(X, axis1=i, axis2=len(rem) + i); rem.pop(i)
            tot += (-t) ** r * np.sum(np.abs(X) ** 2)
    return tot


err = 0.0
for r in (2, 3, 4, 5):
    g = np.zeros(3 ** r)
    for k in range(3):
        g[int(str(k) * r, 3)] = 1 / np.sqrt(3)
    for t in (0.35, 0.4, 0.5):
        err = max(err, abs(qnum(np.outer(g, g), r, t) - ((2 / 3) * (1 + (-t) ** r) + (1 - t) ** r / 3)))
v = np.zeros(27); w = np.zeros(27)
for (x, y, z) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    v[9 * x + 3 * y + z] = 1 / np.sqrt(3)
for (x, y, z) in [(0, 2, 1), (1, 0, 2), (2, 1, 0)]:
    w[9 * x + 3 * y + z] = 1 / np.sqrt(3)
g2 = np.zeros(9); g2[[0, 4, 8]] = 1 / np.sqrt(3)
P = np.outer(v, v) + np.outer(w, w)
for t in (0.4, 0.5):
    err = max(err, abs(qnum(np.kron(P, np.outer(g2, g2)), 5, t) - qnum(P, 3, t) * qnum(np.outer(g2, g2), 2, t)))
Vm = np.column_stack([np.kron(v, g2), np.kron(w, g2)])
rng = np.random.default_rng(1)
for s_ in range(5):
    Hc = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    ops = [np.eye(3)] * 5; ops[s_] = Hc
    Hs = ops[0]
    for o in ops[1:]:
        Hs = np.kron(Hs, o)
    err = max(err, np.abs(Vm.conj().T @ Hs @ Vm - np.trace(Hc) / 3 * np.eye(2)).max())
ok &= err < 1e-12
print(f"(6) GHZ padding formula, multiplicativity and balance of P (x) g_2: max error {err:.1e}; q_3 of the Latin-square plane at t = 1/2: {qnum(P, 3, 0.5):.6f}")
print("PASS" if ok else "FAIL")
