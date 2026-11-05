"""Exact checks (Fractions, rational weights w, so s = w^2 is rational) of the identities in NOTES.md, k = 4, 5, 6.

(R1) G0 elimination: Y00 = 2(I - w w^T/|s|), Y10 w = 0 (Proposition K2), hence
     Y(X, phi) = F(X) + 2|phi_perp + 1/2 Y01 X|^2 with F = Y11 - 1/2 Y10 Y01  (units 2^Z = 4),
     and c(X) := -Y01 X equals sum_r (sigma_r^+ - T_r^T E1^{-1} J^+) psi_r, with c(X) perp w.
(J)  Commutation: T_r J^- = E2 sigma_r^- - sigma_r^- C1 + 2 w_r J^- Pi_r  (maps H_{k-1} -> H_{k-2}, letter-1 strings).
(W)  Woodbury: E2^{-1} = N^{-1} - 1/2 N^{-1} J^- M^{-1} J^+ N^{-1}, M = I + 1/2 J^+ N^{-1} J^-, with
     M_hg = 1/2 w_h w_g/(|s|/2 + s_h + s_g) (h != g), M_hh = 1 + 1/2 sum_{t!=h} s_t/(|s|/2 + s_h + s_t);
     and E2^{-1} J^- = N^{-1} J^- M^{-1}.
(BJ) On the J^- part (psi_r = sum_h X_rh J^-|h>): the B-term equals sum_{h,g} <G_h X_.g, E2^{-1} G_g X_.h> with
     G_h b = (E2 - 2 s_h) sigma_b^-|h> + 2 w_h b_h J^-|h>.
(S)  Soft modes: for a perp w, a_h = 0 and Psi(a,h) = sum_r a_r e_r (x) J^-|h>:
     F(Psi) = 2 s_h |a|^2 - 4 s_h^2 <u, E2^{-1} u> - 1/2 |c|^2, u = sigma_a^-|h>,
     c = 4 s_h [ 1/2 <E1^{-1} e_h, a> w + (E1^{-1} e_h) o a o w ];  |Psi|^2 = (|s| - s_h)|a|^2;  Psi perp p (x) H_{k-2}.
"""
import random
from fractions import Fraction as Fr
from mid import Mid, solve_sym

random.seed(5)


def mat(A, B):
    return [[sum(A[i][l] * B[l][j] for l in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def matvec(A, v):
    return [sum(A[i][j] * v[j] for j in range(len(v))) for i in range(len(A))]


def dot(u, v):
    return sum(x * y for x, y in zip(u, v))


def T(A):
    return [list(r) for r in zip(*A)]


def eye(n):
    return [[Fr(int(i == j)) for j in range(n)] for i in range(n)]


def inv(A):
    return solve_sym(A, eye(len(A)), Fr)


ok = True
for k in (4, 5, 6):
    for trial in range(2):
        w = [Fr(random.randint(1, 30), random.randint(1, 12)) for _ in range(k)]
        if trial == 1:
            w[0] = Fr(1, 40)
        m = Mid(w)
        s, S, n = m.s, m.S, m.np
        N1 = k * n
        Y = m.full()
        Y11 = [row[:N1] for row in Y[:N1]]
        Y10 = [row[N1:] for row in Y[:N1]]
        Y01 = T(Y10)
        Y00 = [row[N1:] for row in Y[N1:]]
        # (R1)
        r1 = all(Y00[h][g] == 2 * (Fr(int(h == g)) - w[h] * w[g] / S) for h in range(k) for g in range(k))
        r1 &= all(v == 0 for v in matvec(Y10, w))
        Jm = m.Jm()
        E1 = m.E1()
        E1i = inv(E1)
        E2 = m.E2()
        E2i = inv(E2)
        Jp = T(Jm)
        # explicit c operator: c = sum_r (sigma_r^+ - T_r^T E1^{-1} J^+) psi_r  -> matrix k x N1
        Cop = [[Fr(0)] * N1 for _ in range(k)]
        for r in range(k):
            sg = m.sig(r)
            T1r = m.T1(r)
            Mr = [[sg[i][h] - sum(Jm[i][g] * sum(E1i[g][l] * T1r[l][h] for l in range(k)) for g in range(k))
                   for h in range(k)] for i in range(n)]  # pairs x holes: sigma_r - J^- E1^{-1} T_r
            for h in range(k):
                for i in range(n):
                    Cop[h][r * n + i] = Mr[i][h]
        r1 &= all(Cop[h][j] == -Y01[h][j] for h in range(k) for j in range(N1))
        r1 &= all(dot(w, [Cop[h][j] for h in range(k)]) == 0 for j in range(N1))
        # (J)
        C1 = [[Fr(2) * s[h] if h == g else Fr(0) for g in range(k)] for h in range(k)]
        rj = True
        for r in range(k):
            L = mat(m.T2(r), Jm)
            Pi = [[Fr(int(h == g == r)) for g in range(k)] for h in range(k)]
            R = mat(E2, m.sig(r))
            R2 = mat(m.sig(r), C1)
            R3 = mat(Jm, Pi)
            rj &= all(L[i][h] == R[i][h] - R2[i][h] + 2 * w[r] * R3[i][h] for i in range(n) for h in range(k))
        # (W)
        Nd = [S / 2 + s[P[0]] + s[P[1]] for P in m.P]
        Ni = [[(1 / Nd[i] if i == j else Fr(0)) for j in range(n)] for i in range(n)]
        Mc = [[Fr(int(h == g)) + sum(Jm[i][h] * Jm[i][g] / Nd[i] for i in range(n)) / 2 for g in range(k)] for h in range(k)]
        rw = all(Mc[h][g] == (w[h] * w[g] / (2 * (S / 2 + s[h] + s[g])) if h != g else
                              1 + sum(s[t] / (S / 2 + s[h] + s[t]) for t in range(k) if t != h) / 2)
                 for h in range(k) for g in range(k))
        Mci = inv(Mc)
        NJ = mat(Ni, Jm)
        Wd = mat(mat(NJ, Mci), T(NJ))
        rw &= all(E2i[i][j] == Ni[i][j] - Wd[i][j] / 2 for i in range(n) for j in range(n))
        L = mat(E2i, Jm)
        R = mat(NJ, Mci)
        rw &= all(L[i][h] == R[i][h] for i in range(n) for h in range(k))
        # reduced form F
        F = [[Y11[i][j] - sum(Y10[i][h] * Y01[h][j] for h in range(k)) / 2 for j in range(N1)] for i in range(N1)]

        def Fq(x):
            return dot(x, matvec(F, x))

        # (BJ)
        X = [[Fr(random.randint(-5, 5), random.randint(1, 4)) for _ in range(k)] for _ in range(k)]
        psi = [[sum(X[r][h] * Jm[i][h] for h in range(k)) for i in range(n)] for r in range(k)]
        Tm = [m.T2(r) for r in range(k)]
        V = [[matvec(Tm[b], psi[a]) for b in range(k)] for a in range(k)]   # V[a][b] = T_b psi_a
        Bterm = sum(dot(V[a][b], matvec(E2i, V[b][a])) for a in range(k) for b in range(k))

        def G(h, b):
            u = [Fr(0)] * n
            for t in range(k):
                if t != h:
                    u[m.pix[tuple(sorted((h, t)))]] += b[t]
            Eu = matvec(E2, u)
            return [Eu[i] - 2 * s[h] * u[i] + 2 * w[h] * b[h] * Jm[i][h] for i in range(n)]

        col = lambda g: [X[r][g] for r in range(k)]
        Bform = sum(dot(G(h, col(g)), matvec(E2i, G(g, col(h)))) for h in range(k) for g in range(k))
        rbj = Bterm == Bform
        # (S)
        rs = True
        for h in range(k):
            # random a with a_h = 0, a.w = 0
            a = [Fr(random.randint(-6, 6), random.randint(1, 5)) if t != h else Fr(0) for t in range(k)]
            j0 = next(t for t in range(k) if t != h)
            a[j0] = 0
            a[j0] = -sum(a[t] * w[t] for t in range(k)) / w[j0]
            x = [Fr(0)] * N1
            for r in range(k):
                for i in range(n):
                    x[r * n + i] = a[r] * Jm[i][h]
            u = [Fr(0)] * n
            for t in range(k):
                if t != h:
                    u[m.pix[tuple(sorted((h, t)))]] += a[t]
            eh = matvec(E1i, [Fr(int(g == h)) for g in range(k)])
            c = [4 * s[h] * (dot(eh, a) * w[g] / 2 + eh[g] * a[g] * w[g]) for g in range(k)]
            cdirect = matvec(Cop, x)
            rs &= c == cdirect
            aa = dot(a, a)
            rs &= Fq(x) == 2 * s[h] * aa - 4 * s[h] ** 2 * dot(u, matvec(E2i, u)) - dot(c, c) / 2
            rs &= dot(x, x) == (S - s[h]) * aa
            # orthogonal to the kernel p (x) H_{k-2}
            rs &= all(sum(x[r * n + i] * w[r] for r in range(k)) == 0 for i in range(n))
        print(f"k={k} w={[str(v) for v in w]}: R1 {r1}  J {rj}  W {rw}  BJ {rbj}  S {rs}")
        ok &= r1 and rj and rw and rbj and rs
print("PASS" if ok else "FAIL")
