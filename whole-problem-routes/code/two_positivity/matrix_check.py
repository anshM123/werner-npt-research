"""Second check of the certificates of Proposition 7.11 (two_positivity/certificates.py) with explicit matrices.
(1) The four separable factors are recomputed as averages over the 60 rotations of the icosahedral group (which acts
    irreducibly on the 5-dim traceless symmetric and the 3-dim antisymmetric parts, so its twirl equals the O(3)
    twirl); this also exhibits each of them as an explicit separable mixture.
(2) For each certificate file (n = 2, 3) the operator C_a^{(x)n} - sum lam_i (M on copy p_i) (x) S ... is built as an
    explicit matrix and its least eigenvalue is computed."""
import numpy as np, itertools, sys, re
from fractions import Fraction as Fr
# icosahedral rotation group: closure of two generators
phi = (1 + 5 ** 0.5) / 2
def rot(axis, ang):
    axis = np.asarray(axis, float); axis /= np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K
g1 = rot([0, 1, phi], 2 * np.pi / 5); g2 = rot([1, 1, 1], 2 * np.pi / 3)
G = [np.eye(3)]
changed = True
while changed:
    changed = False
    for A in list(G):
        for g in (g1, g2):
            B = g @ A
            if not any(np.allclose(B, C, atol=1e-9) for C in G):
                G.append(B); changed = True
print('group order', len(G))
e = np.eye(3); om = sum(np.kron(e[k], e[k]) for k in range(3))
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        F[3 * i + j, 3 * j + i] = 1
PO = np.outer(om, om) / 3; Pan = (np.eye(9) - F) / 2; Ps0 = np.eye(9) - PO - Pan
M = np.eye(9) - 0.5 * np.outer(om, om)
w = np.array([1, 1j, 0]) / np.sqrt(2)
prods = {'e1e1': (e[0], e[0]), 'ww': (w, w), 'wwbar': (w, w.conj()), 'e1e2': (e[0], e[1])}
claimed = {'e1e1': (1/3, 2/15, 0), 'ww': (0, 1/5, 0), 'wwbar': (1/3, 1/30, 1/6), 'e1e2': (0, 1/10, 1/6)}
S = {}
for k, (u, v) in prods.items():
    T = np.zeros((9, 9), dtype=complex)
    for g in G:
        x = np.kron(g @ u, g @ v); T += np.outer(x, x.conj())
    T /= len(G)
    co = claimed[k]
    ok = np.allclose(T, co[0] * PO + co[1] * Ps0 + co[2] * Pan, atol=1e-12)
    S[k] = T.real if np.allclose(T.imag, 0) else T
    print(f'separable factor {k}: icosahedral twirl equals the claimed (P_Om, P_s0, P_an) coefficients {co}: {ok}')
def check(fname):
    lines = [l.split() for l in open(fname) if not l.startswith('#')]
    head = open(fname).readline()
    n = int(re.search(r'n = (\d+)', head).group(1)); a = float(Fr(re.search(r'a = ([0-9/]+)', head).group(1)))
    C = a * np.eye(9) - (F + np.outer(om, om)) / 4
    def kron_list(ops):
        out = np.ones((1, 1))
        for o in ops:
            out = np.kron(out, o)
        return out
    R = kron_list([C] * n)
    for l in lines:
        p = int(l[0]) - 1; others = l[1:-1]; lam = float(Fr(l[-1]))
        ops = [S[s] for s in others]; ops.insert(p, M)
        R = R - lam * kron_list(ops)
    ev = np.linalg.eigvalsh((R + R.conj().T) / 2)
    print(f'{fname}: n = {n}, a = {a}: least eigenvalue of C^n - sum lam_i blocks = {ev[0]:+.3e} (>= 0 up to rounding: {ev[0] > -1e-10})')
for f in sys.argv[1:]:
    check(f)
