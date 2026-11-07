"""Exact check (python-flint, rational amplitudes) of Lemma Q (common-mode shift) for type I blocks of (R_inf)_{k-1}.

With delta_t = C_t + (1/2) C e_q e_rho^T (C the (C_u) block on (X_U, q), C_t the coupling of y_t = (rho_t, zeta_t)),
the claims checked are
  (a) for t != t':  G_tt' = -Ptil_tt' + (w_t w_t'/(4 sigma)) L_t^T L_t' - delta_t^T C^{-1} delta_t',
      with -Ptil given by the closed form (all a = s_t, b = s_t', S = sum_U s_h, h_t = w_U . zeta_t,
      k_0 = K - (a+b)/2, D_a = k_0 + a, D_b = k_0 + b, k_a = K - a/2, k_b = K - b/2, g_0(x) = 1/(x + k_0)):
        y^T(-Ptil)y' = -(a w_t'/(8 D_a k_b)) rho (w_t' rho' + h') - (b w_t/(8 D_b k_a)) rho' (w_t rho + h)
                       - (w_t w_t'/4) sum_h zeta_h zeta'_h/(s_h + k_0) + (w_t w_t'/4) gamma_0 T T',
        T  = w_t rho (1/D_b - 1/D_a) + h/D_b - sum_h w_h zeta_h/(s_h + k_0),
        T' = w_t' rho' (1/D_a - 1/D_b) + h'/D_a - sum_h w_h zeta'_h/(s_h + k_0),
        gamma_0 = 1/(2 + a/D_a + b/D_b + sum_h s_h/(s_h + k_0));
  (b) G itself (all blocks) equals the exact Schur complement of the whole block computed by direct elimination in
      another gauge (coordinates P_{r0,h} deleted), so the formula is the true G;
  (c) the vanishing at w_t = 0: building the block with w_t = 0 (other weights positive), C_t = -(1/2) C e_q e_rho^T
      exactly, i.e. delta_t = 0; and P_tt' + (1/2) e_rho (delta_t')_q + (1/2)(delta_t)_q^T e_rho^T = 0 for every t'.
Usage: python lemmaq_exact.py [seed] [npoints]"""
import sys, os, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'typeI_G_u3'))
import flint
from cf import TypeI

Q = flint.fmpq


def mat(rows):
    return flint.fmpq_mat(rows)


def zeros(r, c):
    return flint.fmpq_mat(r, c)


def eye(n):
    M = flint.fmpq_mat(n, n)
    for i in range(n):
        M[i, i] = 1
    return M


def Emat(s, w, free, sP):
    sF = sum((s[x] for x in free), Q(0))
    n = len(free)
    E = zeros(n, n)
    for i, x in enumerate(free):
        E[i, i] = s[x] + sP + sF / 2
        for j, y in enumerate(free):
            E[i, j] += w[x] * w[y] / 2
    return E


def T1(r, w, free):
    n = len(free)
    A = zeros(n, n)
    if r in free:
        ir = free.index(r)
        for j, x in enumerate(free):
            A[ir, j] += w[x] / 2
            A[j, j] += w[r] / 2
        A[ir, ir] += w[r]
    else:
        for j in range(n):
            A[j, j] = w[r]
    return A


def T2(tp, w, free_in):
    n = len(free_in)
    A = zeros(n + 1, n)
    for j, x in enumerate(free_in):
        A[n, j] = w[x] / 2
        A[j, j] = -w[tp] / 2
    return A


def rows_perm(A, order):
    return mat([[A[i, j] for j in range(A.ncols())] for i in order])


def sub(M, ri, ci):
    return mat([[M[i][j] for j in ci] for i in ri])


def build(wT, wU):
    w = [Q(x) for x in list(wT) + list(wU)]
    s = [x * x for x in w]
    m, u = len(wT), len(wU)
    T = list(range(m)); U = list(range(m, m + u))
    ti = TypeI(w, T, F=Q)
    M = ti.build()
    return ti, M, w, s, m, u, T, U


def formula(wT, wU):
    ti, M, w, s, m, u, T, U = build(wT, wU)
    idx = ti.idx
    sig = sum((s[t] for t in T), Q(0))
    S = sum((s[h] for h in U), Q(0))
    K = sig + S / 2
    zi = [idx[('P', r, h)] for r in U for h in U] + [idx[('Q',)]]
    yi = {t: [idx[('R', t)]] + [idx[('Z', t, h)] for h in U] for t in T}
    C = sub(M, zi, zi)
    Ct = {t: sub(M, zi, yi[t]) for t in T}
    d = u + 1
    L = {}
    for t in T:
        free = [t] + U
        Et = Emat(s, w, free, sig - s[t])
        A2 = T2(t, w, U)
        A2 = rows_perm(A2, [u] + list(range(u)))
        sel = zeros(u, d)
        for i in range(u):
            sel[i, i + 1] = 1
        L[t] = -sel - 2 * A2.transpose() * Et.solve(T1(t, w, free) - w[t] * eye(d))
    nz = C.nrows()
    eq = zeros(nz, 1); eq[nz - 1, 0] = 1
    er = zeros(1, d); er[0, 0] = 1
    Ceq = C * eq
    assert C[nz - 1, nz - 1] == 2
    dl = {t: Ct[t] + Ceq * er / 2 for t in T}
    G = {}
    Ptil = {}
    for t in T:
        for tp in T:
            if t == tp:
                base = zeros(d, d); base[0, 0] = 2
                for i in range(1, d):
                    base[i, i] = 1
                G[(t, t)] = base - (1 - s[t] / sig) / 4 * L[t].transpose() * L[t] - Ct[t].transpose() * C.solve(Ct[t])
            else:
                free2 = [t, tp] + U
                Ett = Emat(s, w, free2, sig - s[t] - s[tp])
                B1 = rows_perm(T2(tp, w, [t] + U), [0, d] + list(range(1, d)))
                B2 = rows_perm(T2(t, w, [tp] + U), [d, 0] + list(range(1, d)))
                P = B1.transpose() * Ett.solve(B2)
                dqt = mat([[dl[t][nz - 1, j] for j in range(d)]])
                dqtp = mat([[dl[tp][nz - 1, j] for j in range(d)]])
                Pt = P - er.transpose() * dqtp / 2 - dqt.transpose() * er / 2
                Ptil[(t, tp)] = Pt
                gau = w[t] * w[tp] / (4 * sig) * L[t].transpose() * L[tp]
                G[(t, tp)] = -Pt + gau - dl[t].transpose() * C.solve(dl[tp])
                # sanity: same as the four-term formula
                G4 = -P - Ct[t].transpose() * C.solve(Ct[tp]) + gau
                G4[0, 0] += Q(1, 2)
                assert G4 == G[(t, tp)]
    return G, Ptil, dl, dict(w=w, s=s, m=m, u=u, T=T, U=U, sig=sig, S=S, K=K, d=d, C=C, Ct=Ct, eq=eq, er=er)


def Ptil_closed(info, t, tp):
    """the closed form of -Ptil_tt' as a d x d matrix (rows: y_t coords, cols: y_t' coords)."""
    w, s, U, K, d = info['w'], info['s'], info['U'], info['K'], info['d']
    a, b = s[t], s[tp]
    k0 = K - (a + b) / 2
    Da, Db = k0 + a, k0 + b
    ka, kb = K - a / 2, K - b / 2
    g0 = [1 / (s[h] + k0) for h in U]
    gam0 = 1 / (2 + a / Da + b / Db + sum((s[h] * g0[i] for i, h in enumerate(U)), Q(0)))
    # linear forms: rho = e0, h = sum_h w_h zeta_h, T, T'
    def vec(f):
        v = [Q(0)] * d
        for i, c in f.items():
            v[i] += c
        return v
    rho = vec({0: Q(1)})
    hvec = vec({i + 1: w[h] for i, h in enumerate(U)})
    Tv = [w[t] * (1 / Db - 1 / Da) * rho[i] + hvec[i] / Db - (w[U[i - 1]] * g0[i - 1] if i > 0 else 0) for i in range(d)]
    Tpv = [w[tp] * (1 / Da - 1 / Db) * rho[i] + hvec[i] / Da - (w[U[i - 1]] * g0[i - 1] if i > 0 else 0) for i in range(d)]
    A = zeros(d, d)
    for i in range(d):
        for j in range(d):
            v = -(a * w[tp] / (8 * Da * kb)) * rho[i] * (w[tp] * rho[j] + hvec[j])
            v += -(b * w[t] / (8 * Db * ka)) * rho[j] * (w[t] * rho[i] + hvec[i])
            if i > 0 and i == j:
                v += -(w[t] * w[tp] / 4) * g0[i - 1]
            v += (w[t] * w[tp] / 4) * gam0 * Tv[i] * Tpv[j]
            A[i, j] = v
    return A


def direct_G(wT, wU):
    ti, M, w, s, m, u, T, U = build(wT, wU)
    r0 = U[0]
    tb = []
    for t in T:
        tb.append(ti.idx[('R', t)])
        for h in U:
            tb.append(ti.idx[('Z', t, h)])
    el = [i for i in range(ti.n) if i not in tb and not (ti.labels[i][0] == 'P' and ti.labels[i][1] == r0)]
    Cm = sub(M, el, el)
    Bm = sub(M, el, tb)
    Am = sub(M, tb, tb)
    return Am - Bm.transpose() * Cm.solve(Bm), m, u


def randq(rng):
    return Q(rng.randint(1, 40), rng.randint(1, 40))


if __name__ == '__main__':
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    npts = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    rng = random.Random(seed)
    nblk = 0
    for it in range(npts):
        m = rng.randint(2, 4); u = rng.randint(1, 4)
        wT = [randq(rng) for _ in range(m)]; wU = [randq(rng) for _ in range(u)]
        G, Ptil, dl, info = formula(wT, wU)
        Gd, _, _ = direct_G(wT, wU)
        d = info['d']
        for t in range(m):
            for tp in range(m):
                blk = mat([[Gd[t * d + i, tp * d + j] for j in range(d)] for i in range(d)])
                assert blk == G[(t, tp)], ('G mismatch', it, t, tp)
                nblk += 1
                if t != tp:
                    assert Ptil_closed(info, t, tp) == -Ptil[(t, tp)], ('Ptil closed form mismatch', it, t, tp)
        print(f'point {it}: m={m} u={u} wT={[str(x) for x in wT]} wU={[str(x) for x in wU]}: G formula = direct '
              f'elimination, -Ptil closed form exact', flush=True)
    print(f'(a),(b): {nblk} blocks at {npts} rational points: PASS')
    # (c) vanishing at w_t = 0
    for it in range(4):
        m = rng.randint(2, 4); u = rng.randint(1, 4)
        wT = [Q(0)] + [randq(rng) for _ in range(m - 1)]; wU = [randq(rng) for _ in range(u)]
        G, Ptil, dl, info = formula(wT, wU)
        assert dl[0] == zeros(dl[0].nrows(), dl[0].ncols()), 'delta_0 != 0 at w_0 = 0'
        for tp in range(1, m):
            assert Ptil[(0, tp)] == zeros(info['d'], info['d']) and Ptil[(tp, 0)] == zeros(info['d'], info['d'])
            assert G[(0, tp)] == zeros(info['d'], info['d'])
        print(f'w_0 = 0 point {it} (m={m}, u={u}): delta_0 = 0, Ptil_0t\' = 0, G_0t\' = 0 exactly', flush=True)
    print('(c): PASS')
