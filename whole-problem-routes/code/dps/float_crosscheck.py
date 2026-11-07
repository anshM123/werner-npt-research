"""Floating-point end-to-end cross-check of a certificate on random complex product vectors (not part of the proof).

For random x, y in C^{27} it evaluates
   LHS = <x (x) y| C_a^{(x)3} |x (x) y> |y|^2            (C_a^{(x)3} built explicitly as a 729 x 729 matrix)
   RHS = t |x|^2 |y|^4 + <psi|Z|psi> + <conj(x) y y|Q_A|conj(x) y y> + <x conj(y) y|Q_B|x conj(y) y>
with psi = x (x) y (x) y written in copy-major order, Z computed (in floats) by the same formula as the exact
verifier, and Q_A = P Q_A_raw P.  The two must agree, and the three last terms must be >= 0.
usage: python float_crosscheck.py cert.txt
"""
import sys
import itertools
import numpy as np

cert = {}
with open(sys.argv[1]) as f:
    for line in f:
        p = line.split()
        cert[p[0]] = p[1:]
n = int(cert['n'][0])
assert n == 3
a = int(cert['a'][0]) / int(cert['a'][1])
t = int(cert['t'][0]) / int(cert['t'][1])
L = int(cert['L'][0])
diagrams = []
for k in range(15):
    pts = list(map(int, cert[f'D{k}']))
    diagrams.append([(pts[2 * i], pts[2 * i + 1]) for i in range(3)])
QA_raw = np.array([int(v) / L for v in cert['QA']]).reshape(15, 15, 15)
QB = np.array([int(v) / L for v in cert['QB']]).reshape(15, 15, 15)


def dmat(m):
    M = np.zeros((27, 27))
    for idx in itertools.product(range(3), repeat=6):
        if all(idx[p] == idx[q] for p, q in m):
            M[9 * idx[0] + 3 * idx[1] + idx[2], 9 * idx[3] + 3 * idx[4] + idx[5]] = 1
    return M


DM = np.array([dmat(m) for m in diagrams])


def find(M):
    h = [b for b in range(15) if np.array_equal(M, DM[b])]
    assert len(h) == 1
    return h[0]


def pt(M, s):
    return np.swapaxes(M.reshape([3] * 6), s, s + 3).reshape(27, 27)


piA = [find(pt(M, 0)) for M in DM]
piB = [find(pt(M, 1)) for M in DM]
S = np.zeros((27, 27))
for i, j, k in itertools.product(range(3), repeat=3):
    S[9 * i + 3 * k + j, 9 * i + 3 * j + k] = 1
muL = [find(S @ M) for M in DM]
muR = [find(M @ S) for M in DM]
# C_a (x) I in diagrams by least squares
F = np.zeros((9, 9))
Om = np.zeros((9, 9))
for i, j in itertools.product(range(3), repeat=2):
    F[3 * j + i, 3 * i + j] = 1
    Om[4 * i, 4 * j] = 1
C = a * np.eye(9) - (F + Om) / 4
CI = np.kron(C, np.eye(3))
cvec = np.linalg.lstsq(DM.reshape(15, -1).T, CI.reshape(-1), rcond=None)[0]
assert np.allclose(np.tensordot(cvec, DM, 1), CI)


def perm_axes(X, pi):
    Y = np.zeros_like(X)
    for idx in itertools.product(range(15), repeat=3):
        Y[pi[idx[0]], pi[idx[1]], pi[idx[2]]] = X[idx]
    return Y


def Pc(X):
    SX = perm_axes(X, muL)
    return (X + SX + perm_axes(X, muR) + perm_axes(SX, muR)) / 4


Ic = np.zeros(15)
Ic[find(np.eye(27))] = 1
Wc = np.einsum('a,b,c->abc', cvec, cvec, cvec)
I3 = np.einsum('a,b,c->abc', Ic, Ic, Ic)
QA = Pc(QA_raw)
Z = Pc(Wc - t * I3 - perm_axes(QA, piA) - perm_axes(QB, piB))


def expval(X, psi):
    """<psi| sum_a X[a] D_a1 (x) D_a2 (x) D_a3 |psi>, psi of shape (27,27,27)"""
    T1 = np.einsum('aij,jkl->aikl', DM, psi)                    # D_a1 on copy 1
    T2 = np.einsum('bkm,aiml->abikl', DM, T1)                   # D_a2 on copy 2
    R3 = np.einsum('cnl,ikn->cikl', DM, psi.conj())             # R3[c,i,k,l] = sum_n conj(psi[i,k,n]) D_c[n,l]
    # value[a,b,c] = <psi|(I(x)I(x)D_c)(D_a(x)D_b(x)I)|psi>
    V = np.einsum('cikl,abikl->abc', R3, T2)
    return np.sum(X * V)


rng = np.random.default_rng(7)
Wfull = np.kron(np.kron(C, C), C)          # on (A1 B1)(A2 B2)(A3 B3)
worst = 0
for trial in range(4):
    x = (rng.standard_normal(27) + 1j * rng.standard_normal(27)).reshape(3, 3, 3)
    y = (rng.standard_normal(27) + 1j * rng.standard_normal(27)).reshape(3, 3, 3)
    if trial == 3:      # a real product minimiser: x = y = e_1^{(x)3}
        x = np.zeros((3, 3, 3), complex)
        x[0, 0, 0] = 1
        y = x.copy()
    phi = np.einsum('ace,bdf->abcdef', x, y).reshape(729)
    lhs = np.real(np.vdot(phi, Wfull @ phi)) * np.vdot(y, y).real

    def psi_of(u, v, w):
        return np.einsum('adg,beh,cfi->abcdefghi', u, v, w).reshape(27, 27, 27)

    psi = psi_of(x, y, y)
    tz = expval(Z, psi).real
    tA = expval(QA, psi_of(x.conj(), y, y)).real
    tB = expval(QB, psi_of(x, y.conj(), y)).real
    tt = t * np.vdot(psi, psi).real
    rhs = tt + tz + tA + tB
    print(f'trial {trial}: LHS {lhs:.12e}  RHS {rhs:.12e}  diff {lhs - rhs:.2e}   terms Z {tz:.3e} QA {tA:.3e} QB {tB:.3e}')
    worst = max(worst, abs(lhs - rhs) / max(1, abs(lhs)))
    assert tz > -1e-9 and tA > -1e-9 and tB > -1e-9
print('max relative discrepancy', worst, 'PASS' if worst < 1e-9 else 'FAIL')
