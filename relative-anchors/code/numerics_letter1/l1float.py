"""Floating-point builder of the letter-1 block of (R_inf)_j (Lemma B, qubit letters only), for exploration.
States: G1 = (a, y) with row letter 1 at copy a and column y a j-subset; G0 = (vac, y) with y a (j+1)-subset.
Returns 2^Z Y (Z = k - j) as a dense numpy array, the state list, and the known kernel vectors p (x) y."""
import itertools
import numpy as np


def subsets(k, c):
    return [frozenset(P) for P in itertools.combinations(range(k), c)]


class L1Block:
    def __init__(self, s, j, bscale=1.0):
        self.s = np.array(s, dtype=float)
        self.w = np.sqrt(self.s)
        self.k = k = len(s)
        self.j = j
        self.bscale = bscale
        self.H = {}
        self.idx = {}
        for c in range(0, k + 1):
            self.H[c] = subsets(k, c)
            self.idx[c] = {y: i for i, y in enumerate(self.H[c])}
        self._Ei = {}

    def dim(self, c):
        return len(self.H[c]) if 0 <= c <= self.k else 0

    def Jp(self, c):
        """J^+ : H_c -> H_{c+1}"""
        M = np.zeros((self.dim(c + 1), self.dim(c)))
        for i, y in enumerate(self.H[c]):
            for t in range(self.k):
                if t not in y:
                    M[self.idx[c + 1][y | {t}], i] += self.w[t]
        return M

    def sig(self, r, c):
        """sigma_r^- : H_c -> H_{c-1} (removes the letter 1 at r)"""
        M = np.zeros((self.dim(c - 1), self.dim(c)))
        for i, y in enumerate(self.H[c]):
            if r in y:
                M[self.idx[c - 1][y - {r}], i] = 1.0
        return M

    def Cdiag(self, c):
        return np.array([2.0 * sum(self.s[t] for t in range(self.k) if t not in y) for y in self.H[c]])

    def E(self, c):
        Jm = self.Jp(c - 1).T if c >= 1 else None
        M = np.diag(self.Cdiag(c))
        if c >= 1:
            Jp = self.Jp(c - 1)
            M = M + 0.5 * Jp @ Jp.T
        return M

    def Einv(self, c):
        if c not in self._Ei:
            self._Ei[c] = np.linalg.inv(self.E(c))
        return self._Ei[c]

    def T(self, r, c):
        """T_r = 1/2 J^+ sigma_r^- + 2 w_r |0><0|_r on H_c"""
        M = 0.5 * self.Jp(c - 1) @ self.sig(r, c) if c >= 1 else np.zeros((self.dim(c), self.dim(c)))
        d = np.array([2.0 * self.w[r] if r not in y else 0.0 for y in self.H[c]])
        return M + np.diag(d)

    def build(self):
        k, j = self.k, self.j
        nH, nH1 = self.dim(j), self.dim(j + 1)
        n1 = k * nH
        N = n1 + nH1
        Y = np.zeros((N, N))
        Ts = [self.T(a, j) for a in range(k)]
        Ss = [self.sig(a, j) for a in range(k)]
        Ei = self.Einv(j)
        ETs = [Ei @ Ts[a] for a in range(k)]
        for a in range(k):
            for b in range(k):
                blk = (0.5 - (1.0 if a == b else 0.0)) * (Ss[b].T @ Ss[a])
                if a == b:
                    blk = blk + np.diag([2.0 - (1.0 if a in y else 0.0) for y in self.H[j]])
                blk = blk - self.bscale * (Ts[b].T @ ETs[a])
                # rows: (a, y') ; columns: (b, y)  -> entry [(a,y'),(b,y)] = <y'| ... |y> with the roles of Lemma B:
                # Y[(a,y'),(b,y)] = <y'|(1/2-[a=b]) s_b^+ s_a^- + [a=b](2-n_a) - T_b^T E^-1 T_a |y>
                Y[a * nH:(a + 1) * nH, b * nH:(b + 1) * nH] = blk
        if nH1 > 0:
            T1 = [self.T(a, j + 1) for a in range(k)]
            S1 = [self.sig(a, j + 1) for a in range(k)]
            Jm = self.Jp(j).T  # H_{j+1} -> H_j
            Ei1 = self.Einv(j + 1)
            for a in range(k):
                blk = -(S1[a] - self.bscale * (Jm @ Ei1 @ T1[a]))   # [y' in H_j, y in H_{j+1}]
                Y[a * nH:(a + 1) * nH, n1:] = blk
                Y[n1:, a * nH:(a + 1) * nH] = blk.T
            G = 2.0 * np.eye(nH1)
            if j + 2 <= k:
                Jp1 = self.Jp(j + 1)
                G = G - self.bscale * (Jp1.T @ self.Einv(j + 2) @ Jp1)
            Y[n1:, n1:] = G
        self.states = [(a, y) for a in range(k) for y in self.H[j]] + [(None, y) for y in self.H[j + 1]]
        return Y

    def kernel(self):
        k, j = self.k, self.j
        nH = self.dim(j)
        N = k * nH + self.dim(j + 1)
        K = []
        for i in range(nH):
            v = np.zeros(N)
            for a in range(k):
                v[a * nH + i] = self.w[a]
            K.append(v)
        return np.array(K).T
