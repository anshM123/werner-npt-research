"""Checks for Section 7 of whole-problem-routes/README.md: the n-copy positivity thresholds of
P_a(X) = a Tr(X) I_3 - (X + X^T)/4.

usage: python map_thresholds.py witness   exact: Theorem 7.1 (the Y-basis GHZ witness) for n <= 8 at rational a
       python map_thresholds.py two       exact: Theorem 7.2 (two copies: blocks, the A-block identity, the PPT
                                          witness of non-decomposability)
       python map_thresholds.py lp        floating point: decomposability thresholds by the LP of Proposition 7.3
       python map_thresholds.py search    floating point: Observation 7.4 (searches around a_n for n = 2..5)
       python map_thresholds.py planes    floating point: the plane identity of Remark 7.5 (n = 1, 2)
       python map_thresholds.py thin      floating point: Lemma 7.7 and Proposition 7.8 (thin copies)
       python map_thresholds.py collision floating point: the icosahedral form of Remark 7.6 (n = 1, 2)
Each mode ends with PASS or FAIL."""
import sys
from fractions import Fraction as Fr
import numpy as np

mode = sys.argv[1]


# ---------------------------------------------------------------------------------------------------------------
# the map on (C^d)^{(x)n}, d = 2 or 3, floating point
def phi_float(a, n, d):
    e = np.eye(d); M = np.zeros((d * d, d * d))
    for i in range(d):
        for j in range(d):
            E = np.outer(e[i], e[j])
            M[:, d * i + j] = (a * np.trace(E) * np.eye(d) - (E + E.T) / 4).reshape(d * d)
    M4 = M.reshape(d, d, d, d); D = d ** n

    def Phi(X):
        T = X.reshape([d] * (2 * n))
        for c in range(n):
            T = np.moveaxis(np.tensordot(M4, T, axes=([2, 3], [c, n + c])), [0, 1], [c, n + c])
        return T.reshape(D, D)
    return Phi


ok = True

if mode == 'witness':
    # x = (p^n + conj(p)^n)/sqrt2, p = (e1 + i e2)/sqrt2: x(b) = 2^{(1-n)/2} cos(pi |b| / 2) on b in {0,1}^n, inside
    # span{e1, e2}^n.  y = x (n odd) or (Y (x) I) x (n even), Y = e1 e2^T - e2 e1^T.
    # Claim: <y|P_a^n(x x^T)|y> = 2^{-n}((2a - 1/2)^n - 1/2).
    # (1) the compression of P_a^n to span{e1, e2}^n is Q_a^n, Q_a(X) = a Tr(X) I_2 - (X + X^T)/4 (floating, n <= 4);
    goodc = True
    for n in (1, 2, 3, 4):
        for a in (0.6, 0.8):
            P3 = phi_float(a, n, 3); P2 = phi_float(a, n, 2)
            emb = [sum(((b >> (n - 1 - c)) & 1) * 3 ** (n - 1 - c) for c in range(n)) for b in range(2 ** n)]
            for _ in range(3):
                X2 = np.random.default_rng(n).normal(size=(2 ** n, 2 ** n))
                X3 = np.zeros((3 ** n, 3 ** n)); X3[np.ix_(emb, emb)] = X2
                goodc &= np.allclose(P3(X3)[np.ix_(emb, emb)], P2(X2))
    print('compression of P_a^n to span{e1, e2}^n equals Q_a^n (n <= 4):', goodc)
    ok &= goodc
    # (2) exact value on the qubit planes for n = 1..8 at rational a: integer arithmetic with 4q Q_{p/q}
    def q_int(p, q, n, Z):
        # applies (4q Q_{p/q})^{(x)n} = (4p Tr(.) I - q (. + .^T))^{(x)n} to an integer 2^n x 2^n matrix Z (object array)
        T = Z.reshape([2] * (2 * n))
        for c in range(n):
            T = np.moveaxis(T, [c, n + c], [0, 1])
            tr = T[0, 0] + T[1, 1]
            new = np.empty_like(T)
            for i in range(2):
                for j in range(2):
                    new[i, j] = (4 * p * tr if i == j else 0) - q * (T[i, j] + T[j, i])
            T = np.moveaxis(new, [0, 1], [c, n + c])
        return T.reshape(2 ** n, 2 ** n)
    for n in range(1, 9):
        zx = np.array([[1, 0, -1, 0][bin(b).count('1') % 4] for b in range(2 ** n)], dtype=object)
        if n % 2:
            zy = zx.copy()
        else:
            zy = np.zeros(2 ** n, dtype=object)
            for b in range(2 ** n):
                first = (b >> (n - 1)) & 1
                nb = b ^ (1 << (n - 1))
                zy[nb] += (-1 if first == 0 else 1) * zx[b]       # Y e1 = -e2, Y e2 = e1 on the first copy
        good = True
        for (p, q) in ((1, 2), (3, 5), (13, 20), (7, 10), (3, 4), (9, 10)):
            out = q_int(p, q, n, np.outer(zx, zx))
            num = int(zy @ out @ zy)
            val = Fr(num, (4 * q) ** n * 2 ** (n - 1) * 2 ** (n - 1))
            a = Fr(p, q)
            good &= (val == Fr(1, 2 ** n) * ((2 * a - Fr(1, 2)) ** n - Fr(1, 2)))
        ok &= good
        print(f'n = {n}: exact value 2^-n((2a - 1/2)^n - 1/2) at a = 1/2, 3/5, 13/20, 7/10, 3/4, 9/10: {good}; '
              f'a_n = (1 + 2^(1-1/n))/4 = {(1 + 2 ** (1 - 1 / n)) / 4:.8f}')
    # (3) floating point: on the qutrit space this y is a bottom eigenvector of P_a^n(xx^T) for a < 3/4, n <= 6
    for n in range(2, 7):
        D = 3 ** n; x = np.zeros(D)
        for b in range(2 ** n):
            x[sum(((b >> (n - 1 - c)) & 1) * 3 ** (n - 1 - c) for c in range(n))] = [1, 0, -1, 0][bin(b).count('1') % 4]
        x /= np.linalg.norm(x)
        good = all(abs(np.linalg.eigvalsh(phi_float(a, n, 3)(np.outer(x, x)))[0] - 2.0 ** -n * ((2 * a - 0.5) ** n - 0.5))
                   < 1e-12 for a in (0.6, 0.65, 0.7))
        ok &= good
        print(f'n = {n}: least eigenvalue of the qutrit output P_a^n(xx^T) equals the formula at a = 0.6, 0.65, 0.7: {good}')

elif mode == 'two':
    import sympy as sp
    a, s1, s2, s3 = sp.symbols('a s1 s2 s3', real=True)
    s = [s1, s2, s3]
    # P_a^{(x)2}(x x^T) for x = sum_i s_i e_i (x) e_i, as a 9 x 9 symbolic matrix
    e = sp.eye(3)
    def Pm(i, j):
        Eij = e[:, i] * e[:, j].T
        return (a * (1 if i == j else 0)) * sp.eye(3) - (Eij + Eij.T) / 4
    out = sp.zeros(9, 9)
    for i in range(3):
        for j in range(3):
            out += s[i] * s[j] * sp.kronecker_product(Pm(i, j), Pm(i, j))
    unit = {s3: sp.sqrt(1 - s1 ** 2 - s2 ** 2)}
    # block A on |kk>, blocks B on {|ij>, |ji>}
    Aidx = [0, 4, 8]
    MA = sp.Matrix(3, 3, lambda p, q: out[Aidx[p], Aidx[q]])
    MA_claim = sp.Matrix(3, 3, lambda p, q: (a ** 2 - (a - sp.Rational(1, 8)) * s[p] ** 2 if p == q else 0) + s[p] * s[q] / 8)
    # the entries use sum s_i^2 = 1, so the comparison is made with s3 = sqrt(1 - s1^2 - s2^2)
    diffA = (MA - MA_claim).applyfunc(lambda z: sp.simplify(sp.expand(z).subs(s3, unit[s3])))
    goodA = all(z == 0 for z in diffA)
    goodB = True; goodZero = True
    for (i, j) in [(0, 1), (0, 2), (1, 2)]:
        u, v = 3 * i + j, 3 * j + i
        blk = sp.Matrix([[out[u, u], out[u, v]], [out[v, u], out[v, v]]])
        claim = sp.Matrix([[a ** 2 - a / 2 * (s[i] ** 2 + s[j] ** 2), s[i] * s[j] / 8],
                           [s[i] * s[j] / 8, a ** 2 - a / 2 * (s[i] ** 2 + s[j] ** 2)]])
        goodB &= all(sp.simplify(sp.expand(z).subs(s3, unit[s3])) == 0 for z in (blk - claim))
    # every entry outside the blocks vanishes
    inblock = set()
    for p in Aidx:
        for q in Aidx:
            inblock.add((p, q))
    for (i, j) in [(0, 1), (0, 2), (1, 2)]:
        u, v = 3 * i + j, 3 * j + i
        for p in (u, v):
            for q in (u, v):
                inblock.add((p, q))
    for p in range(9):
        for q in range(9):
            if (p, q) not in inblock and sp.simplify(out[p, q]) != 0:
                goodZero = False
    print('block decomposition of P_a^2(xx^T) at x = sum s_i e_i e_i (sum s_i^2 = 1): block A as claimed', goodA,
          '; blocks B as claimed', goodB, '; no entries outside the blocks', goodZero)
    ok &= goodA and goodB and goodZero
    # the A-block identity: z^T (M_A - (a - 1/2)^2 I) z = sum_{k<l} [b (s_l^2 z_k^2 + s_k^2 z_l^2) + s_k s_l z_k z_l / 4]
    z = sp.symbols('z1:4', real=True); b = a - sp.Rational(1, 4)
    zv = sp.Matrix(z)
    lhs = (zv.T * (MA_claim - (a - sp.Rational(1, 2)) ** 2 * sp.eye(3)) * zv)[0, 0]
    rhs = sum(b * (s[l] ** 2 * z[k] ** 2 + s[k] ** 2 * z[l] ** 2) + s[k] * s[l] * z[k] * z[l] / 4
              for k in range(3) for l in range(k + 1, 3))
    idA = sp.simplify(sp.expand(lhs - rhs).subs(s3, unit[s3])) == 0
    print('A-block identity (with sum s_i^2 = 1):', idA)
    ok &= idA
    # the B-block bound: a^2 - a r/2 - r/16 >= a^2 - a/2 - 1/16 = (a - 1/4)^2 - 1/8 for r <= 1, root (1 + sqrt2)/4
    r = sp.symbols('r', nonnegative=True)
    idB = sp.expand(a ** 2 - a / 2 - sp.Rational(1, 16) - ((a - sp.Rational(1, 4)) ** 2 - sp.Rational(1, 8))) == 0
    root = sp.nsimplify(sp.solve(sp.Eq((a - sp.Rational(1, 4)) ** 2, sp.Rational(1, 8)), a)[1])
    print('B-block: a^2 - a/2 - 1/16 = (a - 1/4)^2 - 1/8:', idB, '; larger root', root, '=', sp.nsimplify(root - (1 + sp.sqrt(2)) / 4) == 0)
    ok &= idB and sp.simplify(root - (1 + sp.sqrt(2)) / 4) == 0
    # sharpness: at the Bell point s = (1, 1, 0)/sqrt2 the block B_12 has eigenvalue a^2 - a/2 - 1/16
    ev = sp.Matrix([[a ** 2 - a / 2, sp.Rational(1, 16)], [sp.Rational(1, 16), a ** 2 - a / 2]]).eigenvals()
    print('Bell point, block B_12 eigenvalues:', ev)
    # the PPT witness of non-decomposability: rho = (5 P_O(x)P_O + 10 P_O(x)P_an + 4 P_s0(x)P_s0 + 10 P_an(x)P_O)/29
    G = [[Fr(1, 3), Fr(5, 3), Fr(-1)], [Fr(1, 3), Fr(1, 6), Fr(1, 2)], [Fr(-1, 3), Fr(5, 6), Fr(1, 2)]]
    # check G = partial transpose in the projector basis, numerically from the definitions
    om = np.zeros(9); om[[0, 4, 8]] = 1
    F = np.zeros((9, 9))
    for i in range(3):
        for j in range(3):
            F[3 * i + j, 3 * j + i] = 1
    PO = np.outer(om, om) / 3; Pan = (np.eye(9) - F) / 2; Ps0 = np.eye(9) - PO - Pan
    def gam(M):
        return M.reshape(3, 3, 3, 3).transpose(0, 3, 2, 1).reshape(9, 9)
    basis = [PO, Ps0, Pan]; goodG = True
    for c in range(3):
        img = sum(float(G[rr][c]) * basis[rr] for rr in range(3))
        goodG &= np.allclose(gam(basis[c]), img)
    print('G is the partial transpose in the basis (P_Om, P_s0, P_an):', goodG)
    ok &= goodG
    rho = {(0, 0): Fr(5, 29), (0, 2): Fr(10, 29), (1, 1): Fr(4, 29), (2, 0): Fr(10, 29)}
    rG = {(p, q): sum(G[p][i] * G[q][j] * rho.get((i, j), 0) for i in range(3) for j in range(3))
          for p in range(3) for q in range(3)}
    ppt = all(v >= 0 for v in rG.values())
    mult = [1, 5, 3]
    A = sp.Symbol('a')
    cC = [A - 1, A - sp.Rational(1, 4), A + sp.Rational(1, 4)]
    tr = sp.expand(sum(mult[i] * mult[j] * cC[i] * cC[j] * sp.Rational(rho[(i, j)].numerator, rho[(i, j)].denominator)
                       for (i, j) in rho))
    rt = sp.solve(sp.Eq(tr, 0), A)
    print('PPT witness: rho >= 0 True; rho^Gamma >= 0', ppt, '; Tr(C_a^2 rho) =', sp.factor(tr), '; roots', rt)
    ok &= ppt and sp.simplify(max(rt, key=lambda z: float(z)) - (7 + 2 * sp.sqrt(15)) / 22) == 0
    print('so P_a^2 is positive but not decomposable for (1 + sqrt2)/4 <= a < (7 + 2 sqrt15)/22 =',
          float((7 + 2 * sp.sqrt(15)) / 22))

elif mode == 'lp':
    from scipy.optimize import linprog
    def run(d, nmax):
        if d == 3:
            G = np.array([[1 / 3, 5 / 3, -1.0], [1 / 3, 1 / 6, 1 / 2], [-1 / 3, 5 / 6, 1 / 2]]); CC = (-1.0, -0.25, 0.25)
        else:
            G = np.array([[0.5, 1.0, -0.5], [0.5, 0.0, 0.5], [-0.5, 1.0, 0.5]]); CC = (-0.75, -0.25, 0.25)
        assert np.allclose(G @ G, np.eye(3))
        res_all = []
        for n in range(1, nmax + 1):
            T = [(i, j, n - i - j) for i in range(n + 1) for j in range(n + 1 - i)]; ix = {t: k for k, t in enumerate(T)}
            H = np.zeros((len(T), len(T)))
            for t in T:
                poly = {(0, 0, 0): 1.0}
                for rrow, k in zip(G, t):
                    for _ in range(k):
                        new = {}
                        for ex, c in poly.items():
                            for j in range(3):
                                f = list(ex); f[j] += 1; f = tuple(f); new[f] = new.get(f, 0.0) + c * rrow[j]
                        poly = new
                for ex, c in poly.items():
                    H[ix[t], ix[ex]] = c
            def slack(av):
                v = np.array([av + CC[0], av + CC[1], av + CC[2]]) / (av + 0.25)
                cC = np.array([v[0] ** t[0] * v[1] ** t[1] * v[2] ** t[2] for t in T]); m = len(T)
                cobj = np.zeros(m + 1); cobj[-1] = -1
                A = np.vstack([np.hstack([-np.eye(m), np.ones((m, 1))]), np.hstack([H, np.ones((m, 1))])])
                bb = np.concatenate([np.zeros(m), cC])
                r = linprog(cobj, A_ub=A, b_ub=bb, bounds=[(None, None)] * m + [(None, 1.0)], method='highs')
                return -r.fun
            lo, hi = 0.5, 1.0
            for _ in range(40):
                mid = (lo + hi) / 2
                if slack(mid) >= 0:
                    hi = mid
                else:
                    lo = mid
            an = (1 + 2 ** (1 - 1 / n)) / 4
            res_all.append((n, hi, an, slack(0.75)))
            print(f'd = {d}, n = {n:2d}: decomposability threshold {hi:.8f}; a_n = {an:.8f}; LP slack at a = 3/4: '
                  f'{slack(0.75):+.3e}', flush=True)
        return res_all
    # with a qubit ancilla: P_a^k (x) id_2 decomposable would make P_a^k 2-positive, hence P_a^{k+1} positive
    def Gd(d):
        G = np.zeros((3, 3))
        G[:, 0] = [1 / d, 1 / d, -1 / d]
        G[:, 1] = [(1 + d) / 2 - 1 / d, 0.5 - 1 / d, 0.5 + 1 / d]
        G[:, 2] = [(1 - d) / 2, 0.5, 0.5]
        return G
    def anc_slack(k, av):
        G3 = Gd(3); T = [(i, j, k - i - j) for i in range(k + 1) for j in range(k + 1 - i)]; ix = {t: q for q, t in enumerate(T)}
        H = np.zeros((len(T), len(T)))
        for t in T:
            poly = {(0, 0, 0): 1.0}
            for rrow, cnt in zip(G3, t):
                for _ in range(cnt):
                    new = {}
                    for ex, c in poly.items():
                        for j in range(3):
                            f = list(ex); f[j] += 1; f = tuple(f); new[f] = new.get(f, 0.0) + c * rrow[j]
                    poly = new
            for ex, c in poly.items():
                H[ix[t], ix[ex]] = c
        v = np.array([av - 1, av - 0.25, av + 0.25]) / (av + 0.25)
        cC = np.kron(np.array([v[0] ** t[0] * v[1] ** t[1] * v[2] ** t[2] for t in T]), np.array([1.0, 0.0, 0.0]))
        HH = np.kron(H, Gd(2)); m = len(cC)
        cobj = np.zeros(m + 1); cobj[-1] = -1
        A = np.vstack([np.hstack([-np.eye(m), np.ones((m, 1))]), np.hstack([HH, np.ones((m, 1))])])
        bb = np.concatenate([np.zeros(m), cC])
        r = linprog(cobj, A_ub=A, b_ub=bb, bounds=[(None, None)] * m + [(None, 1.0)], method='highs')
        return -r.fun
    goodanc = True
    for k in range(1, 7):
        s_vals = [anc_slack(k, av) for av in (0.75, 0.9, 0.99, 1.0)]
        goodanc &= all(s < 0 for s in s_vals[:3]) and s_vals[3] >= -1e-12
        print(f'P_a^{k} (x) id_2: LP slack at a = 0.75, 0.9, 0.99, 1: ' + ', '.join(f'{s:+.2e}' for s in s_vals))
    print('P_a^k (x) id_2 is decomposable only at a = 1 (k <= 6):', goodanc)
    ok &= goodanc
    r3 = run(3, 12)
    r2 = run(2, 12)
    good = all(abs(h - an) < 1e-7 for (n, h, an, s) in r2)
    good &= abs(r3[1][1] - (7 + 2 * 15 ** 0.5) / 22) < 1e-7
    good &= all(r3[k][1] < r3[k + 1][1] for k in range(len(r3) - 1)) and r3[2][3] < 0
    print('qubit thresholds equal a_n for n <= 12:', all(abs(h - an) < 1e-7 for (n, h, an, s) in r2),
          '; qutrit n = 2 threshold equals (7 + 2 sqrt15)/22:', abs(r3[1][1] - (7 + 2 * 15 ** 0.5) / 22) < 1e-7,
          '; qutrit thresholds increasing, and not decomposable at a = 3/4 for n = 3:', r3[2][3] < 0)
    ok &= good

elif mode == 'search':
    rng = np.random.default_rng(2026)
    for n in (2, 3, 4, 5):
        d = 3; D = d ** n; an = (1 + 2 ** (1 - 1 / n)) / 4
        p = np.array([1, 1j, 0]) / np.sqrt(2); g = np.ones(1, dtype=complex)
        for _ in range(n):
            g = np.kron(g, p)
        ghz = np.real(g + np.conj(g)); ghz /= np.linalg.norm(ghz)
        def rot(v):
            T = v.reshape([d] * n)
            for c in range(n):
                Q, _ = np.linalg.qr(rng.normal(size=(d, d)))
                T = np.moveaxis(np.tensordot(Q, T, axes=([1], [c])), 0, c)
            return T.reshape(D)
        carry = None; res = []
        for a in (an - 0.004, an - 0.001, an + 0.001, an + 0.004):
            Phi = phi_float(a, n, d)
            seeds = [ghz, rot(ghz)] + [rot(ghz) + eps * rng.normal(size=D) for eps in (0.05, 0.3) for _ in range(4)]
            seeds += [rng.normal(size=D) for _ in range(12)]
            if carry is not None:
                seeds = [carry, carry + 0.05 * rng.normal(size=D)] + seeds
            best = (np.inf, None)
            for x in seeds:
                x = x / np.linalg.norm(x); val = np.inf
                for it in range(800):
                    try:
                        w, V = np.linalg.eigh(Phi(np.outer(x, x))); y = V[:, 0]
                        w2, V2 = np.linalg.eigh(Phi(np.outer(y, y)))
                    except np.linalg.LinAlgError:
                        break
                    nv = w2[0]; x = V2[:, 0]
                    if abs(nv - val) < 1e-15:
                        val = nv; break
                    val = nv
                if val < best[0]:
                    best = (val, x)
            carry = best[1]; res.append(best[0])
            T = best[1].reshape([d] * n)
            rk = [int(np.sum(np.linalg.eigvalsh(np.moveaxis(T, c, 0).reshape(d, -1) @ np.moveaxis(T, c, 0).reshape(d, -1).T) > 1e-7)) for c in range(n)]
            print(f'n = {n}, a = a_n {a - an:+.3f}: least value found {best[0]:+.4e}; one-copy ranks of the best x {rk}', flush=True)
        good = res[0] < 0 and res[1] < 0 and res[2] > 0 and res[3] > 0
        print(f'n = {n}: negative just below a_n, none found just above: {good}')
        ok &= good

elif mode == 'planes':
    # Remark 7.5: P_a = sum over the three coordinate planes of iota o Q'_a o c, as bilinear forms on real vectors
    rng = np.random.default_rng(3)
    good = True
    for a in (0.6, 0.7, 0.75, 0.9):
        def Qp(X):
            return a * np.trace(X) * np.eye(2) - (X + X.T) / 4 - 0.5 * (a - 0.5) * np.diag(np.diag(X))
        planes = [(0, 1), (0, 2), (1, 2)]
        def Cp(pi):
            M = np.zeros((2, 3)); M[0, pi[0]] = 1; M[1, pi[1]] = 1
            return M
        for n in (1, 2):
            Pn = phi_float(a, n, 3)
            e = np.eye(2); M = np.zeros((4, 4))
            for i in range(2):
                for j in range(2):
                    M[:, 2 * i + j] = Qp(np.outer(e[i], e[j])).reshape(4)
            M4 = M.reshape(2, 2, 2, 2)
            def Qn(X):
                T = X.reshape([2] * (2 * n))
                for c in range(n):
                    T = np.moveaxis(np.tensordot(M4, T, axes=([2, 3], [c, n + c])), [0, 1], [c, n + c])
                return T.reshape(2 ** n, 2 ** n)
            for _ in range(4):
                x = rng.normal(size=3 ** n); y = rng.normal(size=3 ** n)
                lhs = y @ Pn(np.outer(x, x)) @ y
                rhs = 0.0
                for assign in np.ndindex(*([3] * n)):
                    Cm = np.ones((1, 1))
                    for k in assign:
                        Cm = np.kron(Cm, Cp(planes[k]))
                    xc, yc = Cm @ x, Cm @ y
                    rhs += yc @ Qn(np.outer(xc, xc)) @ yc
                good &= abs(lhs - rhs) < 1e-10
        Ch = sum(np.kron(np.outer(np.eye(2)[i], np.eye(2)[j]), Qp(np.outer(np.eye(2)[i], np.eye(2)[j])))
                 for i in range(2) for j in range(2))
        lmin = np.linalg.eigvalsh(Ch)[0]
        good &= abs(lmin + (1 - a) / 2) < 1e-12
        print(f"a = {a}: plane identity at n = 1, 2 holds; least eigenvalue of the Choi matrix of Q'_a = {lmin:+.4f} "
              f"(= -(1 - a)/2)")
    print('plane decomposition identity and the Choi eigenvalue:', good)
    ok &= good

elif mode == 'thin':
    # Lemma 7.7: for every isometry V: C^2 -> C^3 the Choi matrix of Y -> P_a(V Y V^*) equals
    # a I - (|omega_V><omega_V| + S_V)/4, is invariant under the partial transpose of the output, and is PSD for a >= 3/4.
    rng = np.random.default_rng(77)
    def P(a, X):
        return a * np.trace(X) * np.eye(3) - (X + X.T) / 4
    def choi(a, V):
        C = np.zeros((6, 6), dtype=complex)
        for i in range(2):
            for j in range(2):
                E = np.zeros((2, 2)); E[i, j] = 1
                C[3 * i:3 * i + 3, 3 * j:3 * j + 3] = P(a, V @ E @ V.conj().T)
        return C
    def ptr(C):
        return C.reshape(2, 3, 2, 3).transpose(0, 3, 2, 1).reshape(6, 6)
    good = True; worst = np.inf
    for _ in range(5000):
        V, _r = np.linalg.qr(rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2)))
        C = choi(0.75, V)
        om = sum(np.kron(np.eye(2)[i], V[:, i]) for i in range(2))
        # S_V = sum_ij |e_i><e_j| (x) |conj(V e_j)><conj(V e_i)|
        S = sum(np.kron(np.outer(np.eye(2)[i], np.eye(2)[j]), np.outer(V[:, j].conj(), V[:, i])) for i in range(2) for j in range(2))
        good &= np.allclose(C, 0.75 * np.eye(6) - (np.outer(om, om.conj()) + S) / 4)
        good &= np.allclose(ptr(C), C)
        worst = min(worst, np.linalg.eigvalsh(C).min())
    good &= worst > -1e-12
    Vr = np.eye(3)[:, :2]
    ctrl = np.linalg.eigvalsh(choi(0.74, Vr)).min()
    good &= ctrl < 0
    print(f'Lemma 7.7 at a = 3/4 over 5000 random complex planes: Choi formula, invariance under the partial transpose: '
          f'True; least eigenvalue {worst:+.3e}; control a = 0.74, real plane: {ctrl:+.3e} (negative)')
    # Proposition 7.8 at n = 3, 4: inputs thin (one-copy support of dimension 2) at all but two copies
    for n in (3, 4):
        Phi = phi_float(0.75, n, 3)
        lmin = np.inf
        for _ in range(60):
            thick = rng.choice(n, size=2, replace=False)
            factors = []
            for s in range(n):
                if s in thick:
                    factors.append(np.eye(3))
                else:
                    V, _r = np.linalg.qr(rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2)))
                    factors.append(V)
            W = np.ones((1, 1))
            for F in factors:
                W = np.kron(W, F)
            z = rng.normal(size=W.shape[1]) + 1j * rng.normal(size=W.shape[1])
            x = W @ z; x /= np.linalg.norm(x)
            X = np.outer(x, x.conj())
            out = Phi(np.real(X)) + 1j * Phi(np.imag(X))
            out = (out + out.conj().T) / 2
            lmin = min(lmin, np.linalg.eigvalsh(out).min())
        good &= lmin >= 0.25 ** n - 1e-12
        print(f'Proposition 7.8, n = {n}, a = 3/4: inputs thin at all but two copies (60 random ones): least output eigenvalue '
              f'{lmin:+.3e} >= 4^-n = {0.25 ** n:.3e}: {lmin >= 0.25 ** n - 1e-12}')
        # adversarial: alternating descent inside the class (real planes at the thin copies), from random starts
        dmin = np.inf
        for _ in range(8):
            thick = rng.choice(n, size=2, replace=False)
            W = np.ones((1, 1))
            for s_ in range(n):
                if s_ in thick:
                    W = np.kron(W, np.eye(3))
                else:
                    Q, _r = np.linalg.qr(rng.normal(size=(3, 3)))
                    W = np.kron(W, Q[:, :2])
            z = rng.normal(size=W.shape[1]); val = np.inf
            for it in range(300):
                x = W @ z; x /= np.linalg.norm(x)
                w_, V_ = np.linalg.eigh(Phi(np.outer(x, x))); y = V_[:, 0]
                Ly = W.T @ Phi(np.outer(y, y)) @ W     # P is self-adjoint on real symmetric matrices; W is an isometry
                ev, Z = np.linalg.eigh(Ly)
                z = Z[:, 0]; nv = ev[0]
                if abs(nv - val) < 1e-14:
                    val = nv; break
                val = nv
            dmin = min(dmin, val)
        good &= dmin >= 0.25 ** n - 1e-10
        print(f'Proposition 7.8, n = {n}, a = 3/4: descent inside the class from 8 random starts: least value {dmin:+.6e} '
              f'(4^-n = {0.25 ** n:.6e})')
    ok &= good

elif mode == 'collision':
    gr = (1 + 5 ** 0.5) / 2
    V = np.array([[0, 1, gr], [0, 1, -gr], [1, gr, 0], [1, -gr, 0], [gr, 0, 1], [-gr, 0, 1]], float)
    V /= np.linalg.norm(V, axis=1)[:, None]
    Pl = [np.outer(v, v) for v in V]
    good = np.allclose(sum(Pl), 2 * np.eye(3))
    rng = np.random.default_rng(5)
    for _ in range(5):
        X = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        for a in (0.6, 0.75, 1.0):
            good &= np.allclose(a * np.trace(X) * np.eye(3) - (X + X.T) / 4,
                                sum(np.trace(Pj @ X) / 2 * ((a + 0.25) * np.eye(3) - 1.25 * Pj) for Pj in Pl))
    print('six icosahedral lines: sum of projectors = 2I and P_a(X) = sum_j Tr(E_j X) Q_j:', good)
    for n in (1, 2):
        for _ in range(4):
            x = rng.normal(size=3 ** n) + 1j * rng.normal(size=3 ** n); x /= np.linalg.norm(x)
            y = rng.normal(size=3 ** n) + 1j * rng.normal(size=3 ** n); y /= np.linalg.norm(y)
            rho = np.outer(x, x.conj()); sig = np.outer(y, y.conj())
            E = [Pj / 2 for Pj in Pl]
            strings = list(np.ndindex(*([6] * n)))
            def op(s):
                M = np.ones((1, 1))
                for j in s:
                    M = np.kron(M, E[j])
                return M
            pr = {s: np.real(np.trace(op(s) @ rho)) for s in strings}
            ps = {s: np.real(np.trace(op(s) @ sig)) for s in strings}
            for a in (0.6, 0.75, 0.9):
                lhs = np.real(np.trace(sig @ phi_float(a, n, 3)(rho)))
                r = (9 / 4 - a) / (a + 0.25)
                val = sum(pr[s] * ps[t] * (-r) ** sum(1 for u, v in zip(s, t) if u == v) for s in strings for t in strings)
                good &= abs(lhs - (a + 0.25) ** n * val) < 1e-12
    print('<sigma, P_a^n(rho)> = (a + 1/4)^n E[(-r)^N], r = (9/4 - a)/(a + 1/4), random entangled states, n = 1, 2:', good)
    ok &= good

print('PASS' if ok else 'FAIL')
