"""States of the torus blocks of (R_inf)_2 (rows in H_1, columns in H_2, and the vacuum row with columns in H_3),
for the Lemma B builder j1lib.J1 (whose entry is Lemma B for every j; units 2^{k-2} here).

  family 1 (type I, T = {t, t2}):  R_c = (2@t, 2@t2 1@c), c != t2;  R'_c = (2@t2, 2@t 1@c), c != t;
                                    P_r = (1@r, 2@t 2@t2);  Q_c = (vac, 2@t 2@t2 1@c), c not in T.   4k - 4 states.
  family 2 (type I, T = {t}):       R_P = (2@t, 1@P), P a pair of [k];  P_{r,c} = (1@r, 2@t 1@c), c != t;
                                    Q_P = (vac, 2@t 1@P), P a pair of [k] minus t.
  family 3 (letter-1 block):        P_{r,P} = (1@r, 1@P), P a pair;  Q_P = (vac, 1@P), P a triple.
The kernel vectors p (x) y: for every column y that carries the rows (1@r, y), the vector with entry w_r on (1@r, y).
"""
import itertools
from fractions import Fraction as Fr
from j1lib import J1, onehot


def states_T2(k, t, t2):
    st = [((2, t), onehot(k, {t2: 2, c: 1})) for c in range(k) if c != t2]
    st += [((2, t2), onehot(k, {t: 2, c: 1})) for c in range(k) if c != t]
    st += [((1, r), onehot(k, {t: 2, t2: 2})) for r in range(k)]
    st += [(None, onehot(k, {t: 2, t2: 2, c: 1})) for c in range(k) if c not in (t, t2)]
    return st


def states_T1(k, t):
    st = [((2, t), onehot(k, {a: 1, b: 1})) for a, b in itertools.combinations(range(k), 2)]
    st += [((1, r), onehot(k, {t: 2, c: 1})) for r in range(k) for c in range(k) if c != t]
    U = [c for c in range(k) if c != t]
    st += [(None, onehot(k, {t: 2, a: 1, b: 1})) for a, b in itertools.combinations(U, 2)]
    return st


def states_L1(k):
    st = [((1, r), onehot(k, {a: 1, b: 1})) for r in range(k) for a, b in itertools.combinations(range(k), 2)]
    st += [(None, onehot(k, {a: 1, b: 1, c: 1})) for a, b, c in itertools.combinations(range(k), 3)]
    return st


def kernel_vectors(states, w):
    cols = []
    for (a, y) in states:
        if a is not None and a[0] == 1 and y not in cols:
            cols.append(y)
    out = []
    for y in cols:
        out.append([w[a[1]] if (a is not None and a[0] == 1 and yy == y) else 0 for (a, yy) in states])
    return cols, out


def rank_exact(M):
    A = [list(r) for r in M]
    n, m = len(A), len(A[0])
    rk = 0
    for c in range(m):
        piv = next((i for i in range(rk, n) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        p = A[rk][c]
        for i in range(n):
            if i != rk and A[i][c] != 0:
                f = A[i][c] / p
                A[i] = [x - f * y for x, y in zip(A[i], A[rk])]
        rk += 1
    return rk


def ldl_pd(A):
    """exact symmetric elimination; True iff every pivot is > 0 (A positive definite)."""
    A = [list(r) for r in A]
    n = len(A)
    for c in range(n):
        p = A[c][c]
        if p <= 0:
            return False
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / p
                for j in range(c + 1, n):
                    A[i][j] -= f * A[c][j]
    return True
