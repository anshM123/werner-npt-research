"""Completeness of the reduction at small k.
For each (k, j, i) the class sum C2 = sum of the transpositions of S_u (u = k - 1) acts on the explicit twisted space
V_i; on the isotypic part of an irrep lambda it is the content sum c(lambda).  Exact part: C2 preserves V_i, its only
eigenvalues there are c((u-i, i)) and c((u-i-1, i, 1)), and with m_lambda(i) the dimension of the eigenspace,
sum_i sum_lambda dim(lambda) m_lambda(i) = number of states of the block (exact ranks over Q).  Floating-point part
(not used in any proof): the spectrum of the rescaled block Y~ built by j1lib.py equals the union over (i, lambda) of
the generalized eigenvalues of (Q_lambda, G_lambda), the reduced form of mlib restricted to the eigenspace, each with
multiplicity dim(lambda).
Usage: python check_spectrum.py KMAX"""
import sys, time
from fractions import Fraction as Fr
import numpy as np
import sympy
import mlib, explicit

KMAX = int(sys.argv[1])
T0 = time.time()
ok = True


def content(lam):
    return sum(c - r for r, row in enumerate(lam) for c in range(row))


def hookdim(lam):
    lam = [x for x in lam if x > 0]
    n = sum(lam)
    conj = [sum(1 for x in lam if x > c) for c in range(lam[0])] if lam else []
    prod = 1
    for r, row in enumerate(lam):
        for c in range(row):
            prod *= (row - c - 1) + (conj[c] - r - 1) + 1
    from math import factorial
    return factorial(n) // prod


def act(state, a, b):
    row, y = state
    y = list(y)
    y[a], y[b] = y[b], y[a]
    if row is not None:
        r = row[1]
        r = b if r == a else a if r == b else r
        row = (row[0], r)
    return (row, tuple(y))


nchk = 0
for k in range(2, KMAX + 1):
    u = k - 1
    for w0 in (Fr(1, 2), Fr(3)):
        sig = w0 * w0
        w = [w0] + [Fr(1)] * (k - 1)
        for j in range(-1, k + 1):
            st, M = explicit.block(k, j, w)
            idx = {s: n for n, s in enumerate(st)}
            total = 0
            parts = []
            for i in range(0, u // 2 + 1):
                st2, vecs = explicit.twisted(k, j, i)
                names, norms, Q = mlib.reduced(mlib.Calc(k, j, i, sigma=sig))
                if not names:
                    continue
                n = len(names)
                # C2 on V_i in the basis of the twisted vectors
                C2 = [[Fr(0)] * n for _ in range(n)]
                for bcol, nb in enumerate(names):
                    img = {}
                    for x, c in vecs[nb].items():
                        for a in range(1, k):
                            for b in range(a + 1, k):
                                y = idx[act(st[x], a, b)]
                                img[y] = img.get(y, 0) + c
                    # C2 v_b must lie in V_i: expand and check
                    rec = {}
                    for arow, na in enumerate(names):
                        co = sum(vecs[na].get(x, 0) * c for x, c in img.items()) / norms[arow]
                        C2[arow][bcol] = co
                        for x, cc in vecs[na].items():
                            rec[x] = rec.get(x, 0) + co * cc
                    if any(rec.get(x, 0) != c for x, c in img.items()) or any(img.get(x, 0) != c for x, c in rec.items()):
                        ok = False
                        print('C2 DOES NOT PRESERVE V_i', k, j, i)
                lams = [(u - i, i), (u - i - 1, i, 1)]
                Cm = sympy.Matrix(C2)
                Qm = sympy.Matrix([[Fr(int(x.p), int(x.q)) for x in row] for row in Q])
                Gm = sympy.diag(*norms)
                mult_sum = 0
                for lam in lams:
                    if min(lam) < 0 or any(lam[t] < lam[t + 1] for t in range(len(lam) - 1)):
                        continue  # not a partition
                    cval = content(lam)
                    ns = (Cm - cval * sympy.eye(n)).nullspace()
                    mlam = len(ns)
                    mult_sum += mlam
                    if mlam == 0:
                        continue
                    B = sympy.Matrix.hstack(*ns)
                    Ql = B.T * Qm * B
                    Gl = B.T * Gm * B
                    ev = np.linalg.eigvals(np.array((Gl.inv() * Ql).evalf(30).tolist(), dtype=float)).real
                    d = hookdim(list(lam))
                    total += d * mlam
                    parts += list(ev) * d
                if mult_sum != n:
                    ok = False
                    print('C2 EIGENSPACES DO NOT FILL V_i', k, j, i, mult_sum, n)
            if total != len(st):
                ok = False
                print('DIMENSION COUNT FAILS', k, j, total, len(st))
            full = np.linalg.eigvalsh(np.array([[float(x) for x in row] for row in M]))
            parts = np.sort(np.array(parts))
            if len(parts) != len(full) or np.max(np.abs(parts - full)) > 1e-8 * max(1.0, np.max(np.abs(full))):
                ok = False
                print('SPECTRA DIFFER', k, j, sig)
            nchk += 1
    print(f'k={k}: done ({time.time() - T0:.0f}s)', flush=True)
print(f'{nchk} blocks: dimension counts exact, spectra equal (floating point)')
print('PASS' if ok else 'FAIL', f'{time.time() - T0:.0f}s')
