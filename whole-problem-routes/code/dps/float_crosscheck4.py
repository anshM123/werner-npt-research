"""Floating-point cross-check of a 4-copy certificate (not part of the proof).

Test vectors: x = x12 (x) x34, y = y12 (x) y34 with x12, y12, x34, y34 random complex vectors on two copies each
(entangled between copies 1-2 and between copies 3-4), and the real product minimiser e_1^{(x)4}.
Checks LHS = <x(x)y|C_a^{(x)4}|x(x)y> |y|^2 against t|x|^2|y|^4 + <psi|Z|psi> + <conj(x)yy|Q_A|..> + <x conj(y) y|Q_B|..>.
usage: python float_crosscheck4.py cert.txt
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
assert n == 4
a = int(cert['a'][0]) / int(cert['a'][1])
t = int(cert['t'][0]) / int(cert['t'][1])
L = int(cert['L'][0])
diagrams = [[(int(cert[f'D{k}'][2 * i]), int(cert[f'D{k}'][2 * i + 1])) for i in range(3)] for k in range(15)]
QA_raw = np.array([int(v) / L for v in cert['QA']]).reshape([15] * 4)
QB = np.array([int(v) / L for v in cert['QB']]).reshape([15] * 4)


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


piA = np.array([find(pt(M, 0)) for M in DM])
piB = np.array([find(pt(M, 1)) for M in DM])
S = np.zeros((27, 27))
for i, j, k in itertools.product(range(3), repeat=3):
    S[9 * i + 3 * k + j, 9 * i + 3 * j + k] = 1
muL = np.array([find(S @ M) for M in DM])
muR = np.array([find(M @ S) for M in DM])
F = np.zeros((9, 9))
Om = np.zeros((9, 9))
for i, j in itertools.product(range(3), repeat=2):
    F[3 * j + i, 3 * i + j] = 1
    Om[4 * i, 4 * j] = 1
C = a * np.eye(9) - (F + Om) / 4
cvec = np.linalg.lstsq(DM.reshape(15, -1).T, np.kron(C, np.eye(3)).reshape(-1), rcond=None)[0]


def perm_axes(X, pi):
    Y = np.zeros_like(X)
    Y[np.ix_(pi, pi, pi, pi)] = X
    return Y


def Pc(X):
    SX = perm_axes(X, muL)
    return (X + SX + perm_axes(X, muR) + perm_axes(SX, muR)) / 4


Ic = np.zeros(15)
Ic[find(np.eye(27))] = 1
Wc = np.einsum('a,b,c,d->abcd', cvec, cvec, cvec, cvec)
I4 = np.einsum('a,b,c,d->abcd', Ic, Ic, Ic, Ic)
QA = Pc(QA_raw)
Z = Pc(Wc - t * I4 - perm_axes(QA, piA) - perm_axes(QB, piB))


def pair_vals(psi2):
    """V[a,b] = <psi2| D_a (x) D_b |psi2>, psi2 of shape (27, 27)"""
    T = np.einsum('aij,jk->aik', DM, psi2)            # D_a on first copy
    T = np.einsum('bkl,ail->abik', DM, T)             # D_b on second copy
    return np.einsum('ik,abik->ab', psi2.conj(), T)


def psi2_of(u, v, w):
    # u, v, w of shape (3,3) = (copy1, copy2) ; psi2[(i1 j1 k1),(i2 j2 k2)] = u[i1,i2] v[j1,j2] w[k1,k2]
    return np.einsum('ad,be,cf->abcdef', u, v, w).reshape(27, 27)


def expval(X, pairs):
    V12 = pair_vals(pairs[0])
    V34 = pair_vals(pairs[1])
    return np.einsum('abcd,ab,cd->', X, V12, V34)


rng = np.random.default_rng(11)
worst = 0
for trial in range(4):
    if trial < 3:
        xs = [(rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))) for _ in range(2)]
        ys = [(rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))) for _ in range(2)]
    else:
        e = np.zeros((3, 3), complex)
        e[0, 0] = 1
        xs, ys = [e, e], [e, e]
    lhs = 1.0
    for xp, yp in zip(xs, ys):
        phi = np.einsum('ac,bd->abcd', xp, yp).reshape(81)       # (A1 B1)(A2 B2)
        lhs *= np.vdot(phi, np.kron(C, C) @ phi).real
    ny2 = np.prod([np.vdot(yp, yp).real for yp in ys])
    nx2 = np.prod([np.vdot(xp, xp).real for xp in xs])
    lhs *= ny2
    tz = expval(Z, [psi2_of(xs[k], ys[k], ys[k]) for k in range(2)]).real
    tA = expval(QA, [psi2_of(xs[k].conj(), ys[k], ys[k]) for k in range(2)]).real
    tB = expval(QB, [psi2_of(xs[k], ys[k].conj(), ys[k]) for k in range(2)]).real
    rhs = t * nx2 * ny2 * ny2 + tz + tA + tB
    print(f'trial {trial}: LHS {lhs:.12e}  RHS {rhs:.12e}  diff {lhs - rhs:.2e}   terms Z {tz:.3e} QA {tA:.3e} QB {tB:.3e}')
    worst = max(worst, abs(lhs - rhs) / max(1, abs(lhs)))
    assert tz > -1e-9 and tA > -1e-9 and tB > -1e-9
print('max relative discrepancy', worst, 'PASS' if worst < 1e-9 else 'FAIL')
