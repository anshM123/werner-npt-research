"""Closed form of the type I block of (R_inf)_{k-1} (units of 2Y), derived from Lemma B and Lemma N.

Variables (order of fastdirect.typeI_states): X[r,h] on P_{r,h} (r in [k], h in U), rho_t on R_t, zeta[t,h] on
Z_{t,h}, q on Q.  Column vectors: pattern T lives on R^U (one-hole strings f_h, h in U); pattern T\\t on
R^{U+t}; pattern T\\{t,t'} on R^{U+t+t'}.  E on pattern P with free set F: diag(s_x + s_P + s_F/2) + w_F w_F^T/2.

  2<Psi,Y Psi> = A-part + coupling + 2 q^2 - Bterm, with
  A-part   = 1/2 sum_{r!=h in U} X_rh^2 + 3/2 sum_h X_hh^2 + 1/2 (tr_U X + R)^2 + 3/2 |rho|^2
             + sum_{t,h} (X_th^2 - X_th zeta_th + zeta_th^2)
  coupling = 2 q [ -tr_U X + (w_T^T X_T w_U + 1/2 w_U^T X_U w_U)/K - sum_t (2K rho_t + w_t w_U.zeta_t)/(2K - s_t) ],
             K = s_T + s_U/2
  Bterm    = sum_{r,r'} <T_{r'} psi_r, E_T^{-1} T_r psi_{r'}> + 2 sum_{r,t} <T_{(t,2)} psi_r, E_{T\\t}^{-1} T_r y_t>
             + sum_{t != t'} <T_{(t',2)} y_t, E_{T\\tt'}^{-1} T_{(t,2)} y_t'>.
Generic in the number type (Fraction, or flint fmpq, or float)."""
from fractions import Fraction as Fr


def inv_mat(M, F):
    n = len(M)
    A = [list(M[i]) + [F(1) if i == j else F(0) for j in range(n)] for i in range(n)]
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


class TypeI:
    def __init__(self, w, T, F=Fr):
        self.F = F
        self.w = [F(x) for x in w]
        self.s = [x * x for x in self.w]
        k = self.k = len(w)
        self.T = sorted(T)
        self.U = [c for c in range(k) if c not in self.T]
        U, Tl = self.U, self.T
        idx = {}
        for r in range(k):
            for h in U:
                idx[('P', r, h)] = len(idx)
        for t in Tl:
            idx[('R', t)] = len(idx)
        for t in Tl:
            for h in U:
                idx[('Z', t, h)] = len(idx)
        idx[('Q',)] = len(idx)
        self.idx = idx
        self.n = len(idx)
        self.labels = sorted(idx, key=lambda l: idx[l])

    def Emat(self, P):
        F = self.F
        Fr_ = [x for x in range(self.k) if x not in P]
        sP = sum((self.s[x] for x in P), F(0))
        sF = sum((self.s[x] for x in Fr_), F(0))
        E = [[(self.s[x] + sP + sF / 2 if x == y else F(0)) + self.w[x] * self.w[y] / 2 for y in Fr_] for x in Fr_]
        return Fr_, E

    def build(self):
        F = self.F
        k, U, T, w, s, idx = self.k, self.U, self.T, self.w, self.s, self.idx
        n = self.n
        M = [[F(0)] * n for _ in range(n)]

        def addq(L1, L2, c):
            """M += c * (L1 L2^T + L2 L1^T)/2  (so that v^T M v gains c * L1(v) L2(v))"""
            for i, x in L1.items():
                for j, y in L2.items():
                    v = c * x * y / 2
                    M[i][j] += v
                    M[j][i] += v

        def lin(*pairs):
            d = {}
            for c, i in pairs:
                d[i] = d.get(i, F(0)) + c
            return d

        def ladd(a, b, cb=1):
            d = dict(a)
            for i, x in b.items():
                d[i] = d.get(i, F(0)) + cb * x
            return d

        def lscale(a, c):
            return {i: c * x for i, x in a.items()}

        X = lambda r, h: {idx[('P', r, h)]: F(1)}
        sT = sum((s[t] for t in T), F(0))
        sU = sum((s[h] for h in U), F(0))
        K = sT + sU / 2
        # ---- A-part
        for r in U:
            for h in U:
                if r != h:
                    addq(X(r, h), X(r, h), F(1) / 2)
                else:
                    addq(X(h, h), X(h, h), F(3) / 2)
        trR = {}
        for h in U:
            trR = ladd(trR, X(h, h))
        for t in T:
            trR = ladd(trR, {idx[('R', t)]: F(1)})
        addq(trR, trR, F(1) / 2)
        for t in T:
            addq({idx[('R', t)]: F(1)}, {idx[('R', t)]: F(1)}, F(3) / 2)
            for h in U:
                z = {idx[('Z', t, h)]: F(1)}
                addq(X(t, h), X(t, h), F(1))
                addq(X(t, h), z, F(-1))
                addq(z, z, F(1))
        # ---- coupling with q and the q diagonal
        L = {}
        for h in U:
            L = ladd(L, X(h, h), -1)
        for t in T:
            for h in U:
                L = ladd(L, X(t, h), w[t] * w[h] / K)
        for r in U:
            for h in U:
                L = ladd(L, X(r, h), w[r] * w[h] / (2 * K))
        for t in T:
            L = ladd(L, {idx[('R', t)]: -2 * K / (2 * K - s[t])})
            for h in U:
                L = ladd(L, {idx[('Z', t, h)]: -w[t] * w[h] / (2 * K - s[t])})
        q = {idx[('Q',)]: F(1)}
        addq(q, L, F(2))
        addq(q, q, F(2))
        # ---- B-term
        # psi_r on pattern T: vector over U, entries X[r,h]
        # T_{(r',1)} psi_r (pattern T, free U):
        def T1_patT(rp, r):
            out = {h: {} for h in U}
            if rp in U:
                xw = {}
                for h in U:
                    xw = ladd(xw, X(r, h), w[h])
                out[rp] = ladd(out[rp], xw, F(1) / 2)
                for h in U:
                    out[h] = ladd(out[h], X(r, h), w[rp] / 2)
                out[rp] = ladd(out[rp], X(r, rp), w[rp])
            else:
                for h in U:
                    out[h] = ladd(out[h], X(r, h), w[rp])
            return out
        FT, ET = self.Emat(T)
        ETi = inv_mat(ET, F)
        pos = {x: i for i, x in enumerate(FT)}

        def bil(v1, v2, Ei, posd, c):
            for x, L1 in v1.items():
                for y, L2 in v2.items():
                    g = Ei[posd[x]][posd[y]]
                    if g != 0 and L1 and L2:
                        addq(L1, L2, c * g)
        for r in range(k):
            for rp in range(k):
                bil(T1_patT(rp, r), T1_patT(r, rp), ETi, pos, F(-1))
        # cross: T_{(t,2)} psi_r (pattern T\t, free U+t) and T_{(r,1)} y_t
        def y_vec(t):
            v = {t: {idx[('R', t)]: F(1)}}
            for h in U:
                v[h] = {idx[('Z', t, h)]: F(1)}
            return v

        def T2_on_patT(t, r):
            out = {x: {} for x in U + [t]}
            xw = {}
            for h in U:
                xw = ladd(xw, X(r, h), w[h])
            out[t] = lscale(xw, F(1) / 2)
            for h in U:
                out[h] = ladd(out[h], X(r, h), -w[t] / 2)
            return out

        def T1_on_vec(r, v, free):
            """T_{(r,1)} on a vector v (dict x -> linear form) of one-hole strings with free set `free`."""
            out = {x: {} for x in free}
            if r in free:
                wv = {}
                for x in free:
                    wv = ladd(wv, v[x], w[x])
                out[r] = ladd(out[r], wv, F(1) / 2)
                for x in free:
                    out[x] = ladd(out[x], v[x], w[r] / 2)
                out[r] = ladd(out[r], v[r], w[r])
            else:
                for x in free:
                    out[x] = ladd(out[x], v[x], w[r])
            return out

        for t in T:
            Pt = [x for x in T if x != t]
            Ft, Et = self.Emat(Pt)
            Eti = inv_mat(Et, F)
            post = {x: i for i, x in enumerate(Ft)}
            yt = y_vec(t)
            for r in range(k):
                bil(T2_on_patT(t, r), T1_on_vec(r, yt, Ft), Eti, post, F(-2))
        # letter-2 / letter-2
        def T2_on_vec(tp, v, free):
            """T_{(t',2)} on v (pattern containing t'), result on free + t'."""
            out = {x: {} for x in free + [tp]}
            wv = {}
            for x in free:
                wv = ladd(wv, v[x], w[x])
            out[tp] = lscale(wv, F(1) / 2)
            for x in free:
                out[x] = ladd(out[x], v[x], -w[tp] / 2)
            return out
        for t in T:
            for tp in T:
                if t == tp:
                    continue
                Ptt = [x for x in T if x not in (t, tp)]
                Ftt, Ett = self.Emat(Ptt)
                Etti = inv_mat(Ett, F)
                postt = {x: i for i, x in enumerate(Ftt)}
                v1 = T2_on_vec(tp, y_vec(t), U + [t])
                v2 = T2_on_vec(t, y_vec(tp), U + [tp])
                bil(v1, v2, Etti, postt, F(-1))
        return M
