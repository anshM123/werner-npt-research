"""Explicit reduced form for u = 1 (U = {0}, T = {1..m}): the kernel coordinate x_0 = 0 and the coordinates x_t
(P_{t,0}, t in T) eliminated exactly.  Units of 2Y.  Variables v = (rho_1..rho_m, zeta_1..zeta_m, q) for the states
R_t, Z_{t,0}, Q.  Generic in the number type F (float or Fraction); weights w given (s = w^2)."""


def phi_matrix(w, F=float):
    w = [F(x) for x in w]
    s = [x * x for x in w]
    om, sig = w[0], s[0]
    wt, st = w[1:], s[1:]
    m = len(wt)
    ST = sum(st, F(0))
    e = ST + 2 * sig
    K = ST + sig / 2
    n = 2 * m + 1
    iR = lambda t: t
    iZ = lambda t: m + t
    iq = 2 * m
    S = [[F(0)] * n for _ in range(n)]

    def add_outer(a, b, c):   # S += c * (a b^T + b a^T)/2 ... used as S += c * a a^T when a is b
        for i, x in a.items():
            for j, y in b.items():
                S[i][j] += c * x * y

    # 1/2 R^2 + 3/2 |rho|^2 + sum zeta^2 + 2 q^2
    R = {iR(t): F(1) for t in range(m)}
    add_outer(R, R, F(1, 2) if F is not float else 0.5)
    for t in range(m):
        S[iR(t)][iR(t)] += F(3) / 2
        S[iZ(t)][iZ(t)] += F(1)
    S[iq][iq] += F(2)
    # -2 q sum_t (2K rho_t + w_t om zeta_t)/(2K - s_t)
    for t in range(m):
        c1 = -(2 * K) / (2 * K - st[t]); c2 = -(wt[t] * om) / (2 * K - st[t])
        S[iq][iR(t)] += c1; S[iR(t)][iq] += c1
        S[iq][iZ(t)] += c2; S[iZ(t)][iq] += c2
    # b_t = 1/2 (alpha_t + w_t gamma)
    det = [(ST - st[t] / 2 + 2 * sig) * (ST + st[t] + sig / 2) - sig * st[t] / 4 for t in range(m)]
    alpha = []
    for t in range(m):
        a = {iZ(t): -F(1) - ((st[t] * (ST + st[t] + sig) + sig * e) / 2) / det[t], iR(t): -(om * wt[t] * e) / det[t]}
        alpha.append(a)
    gamma = {iq: 2 * om / K}
    for t in range(m):
        gamma[iZ(t)] = gamma.get(iZ(t), F(0)) + wt[t] * (ST + st[t] + sig) / det[t]
        gamma[iR(t)] = gamma.get(iR(t), F(0)) - om * e / det[t]
    b = []
    for t in range(m):
        bt = {}
        for i, x in alpha[t].items():
            bt[i] = bt.get(i, F(0)) + x / 2
        for i, x in gamma.items():
            bt[i] = bt.get(i, F(0)) + wt[t] * x / 2
        b.append(bt)
    for t in range(m):
        add_outer(b[t], b[t], -F(1))
    wb = {}
    for t in range(m):
        for i, x in b[t].items():
            wb[i] = wb.get(i, F(0)) + wt[t] * x
    add_outer(wb, wb, -F(1) / (2 * sig))
    # pair terms: for t != t' (ordered), - < T_(t',2) y_t , E_(tt')^{-1} T_(t,2) y_t' >
    for t in range(m):
        for tp in range(m):
            if t == tp:
                continue
            wF = [om, wt[t], wt[tp]]; sF = [sig, st[t], st[tp]]
            base = K - st[t] / 2 - st[tp] / 2
            E = [[(sF[i] + base if i == j else F(0)) + wF[i] * wF[j] / 2 for j in range(3)] for i in range(3)]
            Ei = inv3(E)
            # T_(t',2) y_t in coords (0, t, t'):  y_t = zeta_t e_0 + rho_t e_t
            # = 1/2 ( <w,y_t> e_t' - w_t' y_t )
            M1 = {0: {iZ(t): -wt[tp] / 2}, 1: {iR(t): -wt[tp] / 2}, 2: {iZ(t): om / 2, iR(t): wt[t] / 2}}
            # T_(t,2) y_t' in coords (0, t, t'):  y_t' = zeta_t' e_0 + rho_t' e_t'
            M2 = {0: {iZ(tp): -wt[t] / 2}, 1: {iZ(tp): om / 2, iR(tp): wt[tp] / 2}, 2: {iR(tp): -wt[t] / 2}}
            for a in range(3):
                for c in range(3):
                    if Ei[a][c] == 0:
                        continue
                    for i, x in M1[a].items():
                        for j, y in M2[c].items():
                            S[i][j] -= Ei[a][c] * x * y
    # symmetrize (the pair loop adds S[i][j] for ordered pairs; together they are symmetric)
    return S


def inv3(E):
    a, b, c = E[0]; d, e, f = E[1]; g, h, i = E[2]
    A = e * i - f * h; B = -(d * i - f * g); C = d * h - e * g
    D = -(b * i - c * h); E_ = a * i - c * g; F_ = -(a * h - b * g)
    G = b * f - c * e; H = -(a * f - c * d); I = a * e - b * d
    det = a * A + b * B + c * C
    return [[A / det, D / det, G / det], [B / det, E_ / det, H / det], [C / det, F_ / det, I / det]]
