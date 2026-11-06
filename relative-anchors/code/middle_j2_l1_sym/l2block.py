"""The letter-1 block of (R_inf)_2 in the rescaled frame, from Lemma B (units 2^{k-2}), for any field of numbers.

States (same order as r2lib.states_L1): P-states (r, X) = (1@r, 1@X), r in [k], X a pair (r outer, pairs in
lexicographic order inner); Q-states T = (vac, 1@T), T a triple.  Squared weights s (field elements).  For a set Z:
pi_Z = prod_{z in Z} s_z, sig_Z = sum_{z in Z} s_z, L = sum s.

Rescaled frame: each state is multiplied by the product of the weights over its letters 1 (mu of j1lib); then, with
the conjugated operators (D = diag(w_Z) on letter-1 strings):
  Ehat_m = D^{-1} E_m D^{-1} on the m-sets:  Ehat[Z, Z] = 2(L - sig_Z)/pi_Z + 1/2 sum_{t in Z} 1/pi_{Z - t},
           Ehat[Z, Z'] = 1/(2 pi_{Z cap Z'}) if |Z cap Z'| = m - 1;
  t_r e_X = 1/2 sum_{t not in X - r} e_{X - r + t}  if r in X,   2 e_X  if r not in X   (t_r = w_r^{-1} D^{-1} T_(1,r) D);
  jbar^- e_Z = sum_{t in Z} e_{Z - t},  jbar^+ = (jbar^-)^T  (unweighted).
Entries (b = bscale multiplies every term of Lemma B that contains an E^{-1}; b = 1 is the true block):
  PP: Y[(r, X'), (r', X)] = (1/2 - [r = r'])[r in X][X' = X - r + r'] s_r' pi_X + [r = r'][X = X'](2 - [r in X]) s_r pi_X
                            - b s_r s_r' <t_r' e_X', Ehat_2^{-1} t_r e_X>
  PQ: Y[(r, X'), T]       = -[r in T][X' = T - r] pi_T + b s_r (jbar^- Ehat_3^{-1} t_r e_T)[X']
  QQ: Y[T', T]            = 2 [T = T'] pi_T - b (jbar^- Ehat_4^{-1} jbar^+)[T', T]
and, for b = 1 and k >= 5 (Lemma G with Woodbury), Y_QQ = 2 What_3^{-1} with the explicit matrix
  What_3[T, T] = (1 + 1/2 sum_{t not in T} s_t/(2(L - sig_T - s_t)))/pi_T,
  What_3[T, T'] = 1/(4 pi_{T cap T'} (L - sig_{T cup T'})) if |T cap T'| = 2, and 0 otherwise."""
import itertools


def msets(k, m):
    return list(itertools.combinations(range(k), m))


def ginv(M, one):
    """inverse by Gauss-Jordan over a field (elements support + - * /, == 0)"""
    n = len(M)
    A = [list(r) + [one if i == j else one * 0 for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        piv = next(i for i in range(c, n) if not (A[i][c] == 0))
        A[c], A[piv] = A[piv], A[c]
        p = A[c][c]
        pinv = one / p
        A[c] = [x * pinv for x in A[c]]
        for i in range(n):
            if i != c and not (A[i][c] == 0):
                f = A[i][c]
                rc = A[c]
                A[i] = [x - f * y for x, y in zip(A[i], rc)]
    return [r[n:] for r in A]


class Block:
    def __init__(self, s, one, bscale=1, need_qq=True):
        self.s = list(s)
        self.one = one
        self.zero = one * 0
        self.b = bscale
        k = self.k = len(s)
        self.L = sum(self.s[1:], self.s[0])
        self.sets = {m: msets(k, m) for m in range(0, 5)}
        self.idx = {m: {Z: i for i, Z in enumerate(self.sets[m])} for m in self.sets}
        self.P = [(r, X) for r in range(k) for X in self.sets[2]]
        self.Qs = self.sets[3]
        self.Einv = {2: ginv(self.Ehat(2), one), 3: ginv(self.Ehat(3), one)}
        if need_qq and k >= 4:
            self.Einv[4] = ginv(self.Ehat(4), one)

    def pi(self, Z):
        v = self.one
        for z in Z:
            v = v * self.s[z]
        return v

    def sig(self, Z):
        v = self.zero
        for z in Z:
            v = v + self.s[z]
        return v

    def Ehat(self, m):
        S = self.sets[m]
        n = len(S)
        M = [[self.zero] * n for _ in range(n)]
        for i, Z in enumerate(S):
            d = 2 * (self.L - self.sig(Z)) / self.pi(Z)
            for t in Z:
                d = d + self.one / (2 * self.pi(tuple(x for x in Z if x != t)))
            M[i][i] = d
            for j, Z2 in enumerate(S):
                if j != i:
                    I = tuple(x for x in Z if x in Z2)
                    if len(I) == m - 1:
                        M[i][j] = self.one / (2 * self.pi(I))
        return M

    def t_vec(self, r, X):
        """t_r e_X as a dict over sets of the size of X"""
        if r in X:
            Y = tuple(x for x in X if x != r)
            out = {}
            for t in range(self.k):
                if t not in Y:
                    Z = tuple(sorted(Y + (t,)))
                    out[Z] = out.get(Z, 0) + self.one / 2
            return out
        return {tuple(X): 2 * self.one}

    def apply_inv(self, m, v):
        """Ehat_m^{-1} v for a dict v over m-sets; returns a dense list"""
        Ei = self.Einv[m]
        idx = self.idx[m]
        n = len(self.sets[m])
        out = []
        cols = [(idx[Z], c) for Z, c in v.items()]
        for i in range(n):
            row = Ei[i]
            acc = self.zero
            for j, c in cols:
                acc = acc + row[j] * c
            out.append(acc)
        return out

    # ---- blocks ----
    def PP(self):
        k, s, b = self.k, self.s, self.b
        P = self.P
        n = len(P)
        # G[(a, Z)][(c, W)] = <t_a e_Z, Ehat^{-1} t_c e_W>
        tv = {(a, Z): self.t_vec(a, Z) for a in range(k) for Z in self.sets[2]}
        idx2 = self.idx[2]
        Ev = {key: self.apply_inv(2, v) for key, v in tv.items()}
        M = [[self.zero] * n for _ in range(n)]
        pos = {st: i for i, st in enumerate(P)}
        for i, (r, Xp) in enumerate(P):
            for j, (rp, X) in enumerate(P):
                if j < i:
                    continue
                # B-term: s_r s_r' <t_r' e_X', Ehat^{-1} t_r e_X>
                u = tv[(rp, Xp)]
                ev = Ev[(r, X)]
                g = self.zero
                for Z, c in u.items():
                    g = g + c * ev[idx2[Z]]
                val = -(b * s[r] * s[rp] * g) if not (g == 0) else self.zero
                # A-part
                if r in X:
                    rest = tuple(x for x in X if x != r)
                    if rp == r or rp not in rest:
                        Xn = tuple(sorted(rest + (rp,)))
                        if Xn == Xp:
                            coef = (self.one / 2 - (1 if r == rp else 0))
                            val = val + coef * s[rp] * self.pi(X)
                if r == rp and X == Xp:
                    val = val + (2 - (1 if r in X else 0)) * s[r] * self.pi(X)
                M[i][j] = val
                M[j][i] = val
        return M

    def PQ(self):
        k, s, b = self.k, self.s, self.b
        P, Qs = self.P, self.Qs
        idx2 = self.idx[2]
        M = [[self.zero] * len(Qs) for _ in P]
        pos = {st: i for i, st in enumerate(P)}
        for jT, T in enumerate(Qs):
            for r in range(k):
                v = self.apply_inv(3, self.t_vec(r, T))      # dense over triples
                # jbar^- v on pairs: (jbar^- v)[X'] = sum_{t not in X'} v[X' + t]
                for X in self.sets[2]:
                    acc = self.zero
                    for t in range(k):
                        if t not in X:
                            acc = acc + v[self.idx[3][tuple(sorted(X + (t,)))]]
                    val = b * s[r] * acc if not (acc == 0) else self.zero
                    if r in T and X == tuple(x for x in T if x != r):
                        val = val - self.pi(T)
                    M[pos[(r, X)]][jT] = val
        return M

    def QQ(self):
        """direct form 2 pi_T delta - b jbar^- Ehat_4^{-1} jbar^+ (needs Einv[4])"""
        k, b = self.k, self.b
        Qs = self.Qs
        n = len(Qs)
        M = [[self.zero] * n for _ in range(n)]
        cols = []
        for T in Qs:
            v = {}
            for t in range(k):
                if t not in T:
                    v[tuple(sorted(T + (t,)))] = self.one
            cols.append(self.apply_inv(4, v))
        for i, Tp in enumerate(Qs):
            for j, T in enumerate(Qs):
                acc = self.zero
                for t in range(k):
                    if t not in Tp:
                        acc = acc + cols[j][self.idx[4][tuple(sorted(Tp + (t,)))]]
                M[i][j] = (2 * self.pi(T) if i == j else self.zero) - b * acc
        return M

    def What3(self):
        """the explicit matrix of Lemma G (k >= 5): Y_QQ = 2 What3^{-1} at b = 1"""
        k, s, L = self.k, self.s, self.L
        Qs = self.Qs
        n = len(Qs)
        M = [[self.zero] * n for _ in range(n)]
        for i, T in enumerate(Qs):
            N = self.zero
            for t in range(k):
                if t not in T:
                    N = N + s[t] / (2 * (L - self.sig(T) - s[t]))
            M[i][i] = (self.one + N / 2) / self.pi(T)
            for j, T2 in enumerate(Qs):
                if j != i:
                    I = tuple(x for x in T if x in T2)
                    if len(I) == 2:
                        U = tuple(sorted(set(T) | set(T2)))
                        M[i][j] = self.one / (4 * self.pi(I) * (L - self.sig(U)))
        return M

    def full(self):
        PP, PQ, QQ = self.PP(), self.PQ(), self.QQ()
        nP, nQ = len(PP), len(QQ)
        M = [PP[i] + PQ[i] for i in range(nP)]
        M += [[PQ[i][j] for i in range(nP)] + QQ[j] for j in range(nQ)]
        return M


def kernel_indicators(k):
    """the rescaled kernel vectors p (x) |1@P>: indicators of the column groups {(r, P) : r in [k]} (P-coordinates)"""
    P = [(r, X) for r in range(k) for X in msets(k, 2)]
    out = []
    for X in msets(k, 2):
        out.append([1 if XX == X else 0 for (r, XX) in P])
    return out
