"""Fast floating-point direct construction of Y = A - B (unrescaled) for the group of k copies, p_G = sum_t w_t |1_t>.

Same definitions as check_bath_reduction.group_ops (single-excitation-range/code):
  Mel(b, y, a, z) = prod_c m(b_c, y_c, a_c, z_c),  m(b, y, a, z) = [a = b][y = z] - 1/2 [b = y][a = z];
  D[x, x'] = sum_{t,t'} w_t w_t' Mel(e_t, x, e_t', x');  MV_b[y, z] = sum_t w_t Mel(b, y, e_t, z);
  A(a, y', b, y) = Mel(b, y', a, y);  B(a, y', b, y) = MV_b[y', :] D^{-1} MV_a[y, :]^T.
The tensor-product structure is used (Kronecker products of 3 x 3 factors), so k <= 7 is cheap."""
import itertools
import numpy as np


def m1(b, y, a, z):
    return (1.0 if (a == b and y == z) else 0.0) - (0.5 if (b == y and a == z) else 0.0)


def fac(b, a):
    """3x3 matrix [y, z] -> m(b, y, a, z) for fixed letters b, a."""
    return np.array([[m1(b, y, a, z) for z in range(3)] for y in range(3)])


def kron_all(mats):
    out = np.ones((1, 1))
    for M in mats:
        out = np.kron(out, M)
    return out


class Direct:
    def __init__(self, w):
        self.w = np.asarray(w, float)
        k = self.k = len(w)
        self.e = [tuple(1 if c == t else 0 for c in range(k)) for t in range(k)]
        self.cols = list(itertools.product(range(3), repeat=k))
        self.pos = {y: i for i, y in enumerate(self.cols)}
        D = np.zeros((3 ** k, 3 ** k))
        for t in range(k):
            for tp in range(k):
                D += self.w[t] * self.w[tp] * kron_all([fac(self.e[t][c], self.e[tp][c]) for c in range(k)])
        self.Di = np.linalg.inv(D)
        self._mv = {}
        self._mvd = {}

    def MVrow(self, b, y):
        key = (b, y)
        if key not in self._mv:
            k = self.k
            v = np.zeros(3 ** k)
            for t in range(k):
                vec = np.ones(1)
                for c in range(k):
                    vec = np.kron(vec, fac(b[c], self.e[t][c])[y[c]])
                v += self.w[t] * vec
            self._mv[key] = v
        return self._mv[key]

    def A(self, a, yp, b, y):
        v = 1.0
        for c in range(self.k):
            v *= m1(b[c], yp[c], a[c], y[c])
            if v == 0.0:
                return 0.0
        return v

    def B(self, a, yp, b, y):
        key = (b, yp)
        if key not in self._mvd:
            self._mvd[key] = self.MVrow(b, yp) @ self.Di
        return self._mvd[key] @ self.MVrow(a, y)

    def Y(self, X1, X2):
        (a, yp), (b, y) = X1, X2
        return self.A(a, yp, b, y) - self.B(a, yp, b, y)

    def matrix(self, states):
        return np.array([[self.Y(X1, X2) for X2 in states] for X1 in states])


def typeI_states(k, T):
    """States of the type I block of (R_inf)_{k-1} for the set T, in the order
    P_{r,h} (r in [k], h in U), R_t (t in T), Z_{t,h} (t in T, h in U), Q."""
    T = sorted(T); U = [c for c in range(k) if c not in T]
    one = lambda r, l: tuple(l if c == r else 0 for c in range(k))
    vac = (0,) * k
    def col(two, zero):   # 2 on 'two', 0 at 'zero', 1 elsewhere
        return tuple(2 if c in two else (0 if c == zero else 1) for c in range(k))
    st, lab = [], []
    for r in range(k):
        for h in U:
            st.append((one(r, 1), col(T, h))); lab.append(('P', r, h))
    for t in T:
        st.append((one(t, 2), col([x for x in T if x != t], t))); lab.append(('R', t))
    for t in T:
        for h in U:
            st.append((one(t, 2), tuple(2 if (c in T and c != t) else (1 if c == t else (0 if c == h else 1)) for c in range(k))))
            lab.append(('Z', t, h))
    st.append((vac, tuple(2 if c in T else 1 for c in range(k)))); lab.append(('Q',))
    return st, lab
