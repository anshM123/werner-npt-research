"""Closed form of the letter-1 torus block of (R_inf)_{k-2} (Lemma B with j = k-2, Z = 2), any number type.

Hole notation: a letter-1 column string is identified with its set of holes (zeros).
  level k-2: pairs P = (h, h'), h < h';   level k-1: single holes h;   level k: no hole.
States: P_{r,P} = (1@r, holes P), r in [k], P a pair;  Q_h = (vac, hole h).
Units: 2^Z <Psi, Y Psi> with Z = 2 (multiply by 1/4 to get Y).

Operators (Lemma N, README Section 1):
  T_r |P> = 1/2 w_r |P> + 1/2 sum_{t in P} w_t |P - t + r>   (r not in P),   T_r |P> = 2 w_r |P>   (r in P)
  T_r |h> = 1/2 w_h |r> + 1/2 w_r |h> (h != r),   T_r |r> = 2 w_r |r>
  E2 = E_{k-2} = diag(|s|/2 + s_P) + 1/2 J^- J^+  on pairs,  E1 = E_{k-1} = diag(|s|/2 + s_h) + 1/2 w w^T on holes
  J^- |h> = sum_{t != h} w_t |{h,t}>.
"""
import itertools
from fractions import Fraction as Fr


def pairs(k):
    return list(itertools.combinations(range(k), 2))


def solve_sym(M, B, F):
    """Solve M X = B exactly (Gauss-Jordan), M n x n list, B n x m list."""
    n = len(M)
    m = len(B[0]) if B else 0
    A = [list(M[i]) + list(B[i]) for i in range(n)]
    for c in range(n):
        r = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[r] = A[r], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]


class Mid:
    """Letter-1 block of (R_inf)_{k-2}; w a list of weights (Fractions or floats); F the number type."""

    def __init__(self, w, F=Fr):
        self.F = F
        self.w = [F(x) for x in w]
        k = self.k = len(w)
        self.s = [x * x for x in self.w]
        self.S = sum(self.s)
        self.P = pairs(k)
        self.pix = {P: i for i, P in enumerate(self.P)}
        self.np = len(self.P)

    # ---- operators as dense matrices (lists) ----
    def T2(self, r):
        """T_r on pairs: matrix M[P', P] = <P'|T_r|P>."""
        F, w, n = self.F, self.w, self.np
        M = [[F(0)] * n for _ in range(n)]
        for i, P in enumerate(self.P):
            if r in P:
                M[i][i] += 2 * w[r]
            else:
                M[i][i] += w[r] / 2
                for t in P:
                    Q = tuple(sorted([x for x in P if x != t] + [r]))
                    M[self.pix[Q]][i] += w[t] / 2
        return M

    def T1(self, r):
        F, w, k = self.F, self.w, self.k
        M = [[F(0)] * k for _ in range(k)]
        for h in range(k):
            if h == r:
                M[r][r] += 2 * w[r]
            else:
                M[r][h] += w[h] / 2
                M[h][h] += w[r] / 2
        return M

    def Jm(self):
        """J^- : holes -> pairs, matrix [P][h]."""
        F, w = self.F, self.w
        M = [[F(0)] * self.k for _ in range(self.np)]
        for i, P in enumerate(self.P):
            h, hp = P
            M[i][h] += w[hp]
            M[i][hp] += w[h]
        return M

    def sig(self, r):
        """sigma_r^- : holes -> pairs."""
        F = self.F
        M = [[F(0)] * self.k for _ in range(self.np)]
        for h in range(self.k):
            if h != r:
                M[self.pix[tuple(sorted((h, r)))]][h] = F(1)
        return M

    def E2(self):
        F, s, S = self.F, self.s, self.S
        Jm = self.Jm()
        n = self.np
        M = [[sum(Jm[i][h] * Jm[j][h] for h in range(self.k)) / 2 for j in range(n)] for i in range(n)]
        for i, P in enumerate(self.P):
            M[i][i] += S / 2 + s[P[0]] + s[P[1]]
        return M

    def E1(self):
        F, s, S, w, k = self.F, self.s, self.S, self.w, self.k
        return [[(S / 2 + s[h] if h == g else F(0)) + w[h] * w[g] / 2 for g in range(k)] for h in range(k)]

    # ---- the quadratic form ----
    def states(self):
        return [('P', r, P) for r in range(self.k) for P in self.P] + [('Q', h) for h in range(self.k)]

    def g1index(self, r, i):
        return r * self.np + i

    def full(self):
        """Full block (units 2^Z = 4): returns matrix over states() order."""
        F, k, n, w, s, S = self.F, self.k, self.np, self.w, self.s, self.S
        N1 = k * n
        N = N1 + k
        Y = [[F(0)] * N for _ in range(N)]
        # A part
        for r in range(k):
            for i, P in enumerate(self.P):
                a = self.g1index(r, i)
                Y[a][a] += F(1) / 2 if r not in P else F(2)
        # crossed swap terms: + 1/2 sum_{a != b} sum_{x not in {a,b}} X[a,{a,x}] X[b,{b,x}]   (as symmetric matrix)
        # and the -1/2 [r_a = r_b] |sigma_r psi_r|^2 term was absorbed: diag for P not containing r is 1 - 1/2 = 1/2,
        # for P containing r it is 2 (sigma_r kills).  Swap term for a != b:
        for x in range(k):
            for a_ in range(k):
                if a_ == x:
                    continue
                for b_ in range(k):
                    if b_ == x or b_ == a_:
                        continue
                    ia = self.g1index(a_, self.pix[tuple(sorted((a_, x)))])
                    ib = self.g1index(b_, self.pix[tuple(sorted((b_, x)))])
                    Y[ia][ib] += F(1) / 2
        # B part: - sum_{a,b} <T_b psi_a, E^{-1} T_a psi_b>
        T = [self.T2(r) for r in range(k)]
        E = self.E2()
        # E^{-1} T_a for all a: solve
        Einv = solve_sym(E, [[F(1) if i == j else F(0) for j in range(n)] for i in range(n)], F)
        # B[(a,P'),(b,P)] = (T_b e_P')^T Einv (T_a e_P)
        TbT = [[[T[b][q][p] for q in range(n)] for p in range(n)] for b in range(k)]  # TbT[b][p] = column p of T_b
        EinvT = []
        for a_ in range(k):
            # Einv @ T_a  (n x n)
            EinvT.append([[sum(Einv[i][q] * T[a_][q][p] for q in range(n)) for p in range(n)] for i in range(n)])
        for a_ in range(k):
            for b_ in range(k):
                for pp in range(n):          # psi_a index P'
                    col_b = TbT[b_][pp]      # T_b e_{P'}
                    for p in range(n):       # psi_b index P
                        v = sum(col_b[i] * EinvT[a_][i][p] for i in range(n))
                        if v != 0:
                            Y[self.g1index(a_, pp)][self.g1index(b_, p)] -= v
        # coupling: -2 sum_a <psi_a, (sigma_a - J^- E1^{-1} T_a) phi>  -> off-diagonal entries -1 times (each side)
        Jm = self.Jm()
        E1 = self.E1()
        for a_ in range(k):
            T1a = self.T1(a_)
            X = solve_sym(E1, T1a, F)            # E1^{-1} T_a  (k x k)
            sg = self.sig(a_)
            Cm = [[sg[i][h] - sum(Jm[i][g] * X[g][h] for g in range(k)) for h in range(k)] for i in range(n)]
            for i in range(n):
                for h in range(k):
                    if Cm[i][h] != 0:
                        Y[self.g1index(a_, i)][N1 + h] -= Cm[i][h]
                        Y[N1 + h][self.g1index(a_, i)] -= Cm[i][h]
        # G0 block: 2 (I - w w^T/|s|)
        for h in range(k):
            for g in range(k):
                Y[N1 + h][N1 + g] += 2 * ((F(1) if h == g else F(0)) - w[h] * w[g] / S)
        return Y
