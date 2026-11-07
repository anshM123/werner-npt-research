"""Closed form of the raw pair R_raw,tt' (psi, q and the common mean p eliminated) in measure coordinates.

Normalisation S = sum_U s_h = 1.  Free data: a finite measure given by locations x_i in [0, 1] and masses m_i >= 0
with sum m_i = 1 (physical case: one atom per free copy, m_h = x_h = s_h).  Letter-2 pair: a = s_t, b = s_t',
sigma = sum_T s >= a + b, wa = sqrt(a), wb = sqrt(b).  Coordinates of one letter-2 copy: (rho, Z_1, ..., Z_n) with
zeta_h = w_h Z_h in the physical case.  E[phi] = sum_i m_i phi(x_i).

The bilinear form (convention: the cross term of the quadratic form is 2 R(y, y')) is

  R(y, y') = (wa wb/(4 sigma)) E[l_t l_t'] - (wa wb/4) E[Z Z' g0] + rho wb (wb rho' + h')/(4 Da)
             + rho' wa (wa rho + h)/(4 Db) + (wa wb/4) gam0 T T' + rho rho'/2 - Lt Lt'/2 - bt bt'/cp,

with every symbol defined in `pieces` below.  Generic in the number type: pass `F` (float, Fraction, flint.fmpq) and,
for exact work, wa, wb as exact numbers whose squares are a, b."""


def stats(x, m, a, b, sig, F=float):
    K = sig + F(1) / 2
    ka, kb, k0 = K - a / 2, K - b / 2, K - (a + b) / 2
    st = dict(K=K, ka=ka, kb=kb, k0=k0, da=K + a / 2, db=K + b / 2, Da=k0 + a, Db=k0 + b)
    st['Ff'] = sum((mi / (K + xi) for xi, mi in zip(x, m)), F(0))
    st['Ba'] = sum((mi / (ka + xi) for xi, mi in zip(x, m)), F(0))
    st['Bb'] = sum((mi / (kb + xi) for xi, mi in zip(x, m)), F(0))
    st['B0'] = sum((mi / (k0 + xi) for xi, mi in zip(x, m)), F(0))
    st['Cab'] = sum((mi / ((ka + xi) * (kb + xi)) for xi, mi in zip(x, m)), F(0))
    st['alpha'] = 1 / (2 + st['Ff'])
    st['gama'] = 1 / (2 + a / st['da'] + st['Ba'])
    st['gamb'] = 1 / (2 + b / st['db'] + st['Bb'])
    st['gam0'] = 1 / (2 + a / st['Da'] + b / st['Db'] + st['B0'])
    st['cp'] = sig - 2 * sig ** 2 * st['alpha'] * st['Ff'] - sig ** 2 / (2 * K ** 2)
    return st


def lin_features(x, m, st, F=float):
    """linear functionals of a field Z as coefficient vectors over the atoms (length n): h, z^a, z^b, z^ab, z^0."""
    ka, kb, k0 = st['ka'], st['kb'], st['k0']
    h = [mi for mi in m]
    za = [mi / (ka + xi) for xi, mi in zip(x, m)]
    zb = [mi / (kb + xi) for xi, mi in zip(x, m)]
    zab = [mi / ((ka + xi) * (kb + xi)) for xi, mi in zip(x, m)]
    z0 = [mi / (k0 + xi) for xi, mi in zip(x, m)]
    return dict(h=h, za=za, zb=zb, zab=zab, z0=z0)


def copy_forms(st, feats, a, wa, sig, kself, dself, gamself, Bself, zself_key, Dother, Dself, F=float):
    """For one copy (weight a, root wa; kself = k_a, dself = delta_a, gamself = gamma_a, Bself = B_a), return the
    linear forms (as vectors over (rho, Z_1..Z_n)) of: u, nu, A, B (of l), L (q source), beta, b (p half-source), T,
    and the vectors V = E[e_a Z], W = E[e_a g_other Z] are built by the caller."""
    n = len(feats['h'])
    K = st['K']
    def vec(r, zcoef):
        return [r] + list(zcoef)
    def comb(*terms):
        out = [F(0)] * (n + 1)
        for c, v in terms:
            for i in range(n + 1):
                out[i] += c * v[i]
        return out
    e_rho = vec(F(1), [F(0)] * n)
    h = vec(F(0), feats['h'])
    zs = vec(F(0), feats[zself_key])
    z0 = vec(F(0), feats['z0'])
    u = comb((wa, e_rho), (F(1) / 2, h))
    nu = comb((gamself / dself, u), (-gamself / 2, zs))
    A = comb((-1 / dself, u), (a / dself, nu))
    B = comb((-a, nu))
    L = comb((-K / kself, e_rho), (-wa / (2 * kself), h))
    # beta = rho/2 - <lv, E^{-1} v>
    # P1 = rho/(4 d) - (wa/2)[ (wa rho + h) Bself/2 + h + (1/2 - k) z^self ]
    inner = comb((wa * Bself / 2, e_rho), (Bself / 2, h), (F(1), h), (F(1) / 2 - kself, zs))
    P1 = comb((1 / (4 * dself), e_rho), (-wa / 2, inner))
    Wv = comb((wa / (2 * dself), e_rho), (F(1), inner))
    Wl = (wa / 2) * (1 / dself - Bself)
    beta = comb((F(1) / 2, e_rho), (-1, P1), (gamself * Wl, Wv))
    bvec = comb((F(1), beta), (sig / (2 * K), L))
    T = comb((wa * (1 / Dother - 1 / Dself), e_rho), (1 / Dother, h), (-1, z0))
    return dict(u=u, nu=nu, A=A, B=B, L=L, beta=beta, b=bvec, T=T, h=h, rho=e_rho)


def raw_pair(x, m, a, b, sig, wa, wb, F=float):
    """matrix R (n+1) x (n+1) with R(y, y') = y^T R y' in coordinates (rho, Z_i)."""
    n = len(x)
    st = stats(x, m, a, b, sig, F)
    fe = lin_features(x, m, st, F)
    cj = copy_forms(st, fe, a, wa, sig, st['ka'], st['da'], st['gama'], st['Ba'], 'za', st['Db'], st['Da'], F)
    ck = copy_forms(st, fe, b, wb, sig, st['kb'], st['db'], st['gamb'], st['Bb'], 'zb', st['Da'], st['Db'], F)
    M = [[F(0)] * (n + 1) for _ in range(n + 1)]

    def add_outer(c, u, v):
        for i in range(n + 1):
            if u[i] != 0:
                for j in range(n + 1):
                    if v[j] != 0:
                        M[i][j] += c * u[i] * v[j]

    def lv(rc, zc):
        return [rc] + list(zc)
    zero = [F(0)] * n
    ka, kb, k0, K = st['ka'], st['kb'], st['k0'], st['K']
    # local terms: (wa wb/(4 sig)) E[e_a e_b Z Z'] - (wa wb/4) E[g0 Z Z']
    for i, (xi, mi) in enumerate(zip(x, m)):
        ea = (K + xi) / (ka + xi); eb = (K + xi) / (kb + xi)
        M[i + 1][i + 1] += wa * wb / (4 * sig) * mi * ea * eb - wa * wb / 4 * mi / (k0 + xi)
    # finite part of E[l_j l_k]:  -A_k V_j - B_k W_j - A_j V_k - B_j W_k + A_j A_k + A_j B_k Bb + B_j A_k Ba + B_j B_k Cab
    Vj = lv(F(0), [hh + a / 2 * z for hh, z in zip(fe['h'], fe['za'])])
    Wj = lv(F(0), [z + a / 2 * zz for z, zz in zip(fe['zb'], fe['zab'])])
    Vk = lv(F(0), [hh + b / 2 * z for hh, z in zip(fe['h'], fe['zb'])])
    Wk = lv(F(0), [z + b / 2 * zz for z, zz in zip(fe['za'], fe['zab'])])
    c = wa * wb / (4 * sig)
    add_outer(-c, Vj, ck['A']); add_outer(-c, Wj, ck['B'])
    add_outer(-c, cj['A'], Vk); add_outer(-c, cj['B'], Wk)
    add_outer(c, cj['A'], ck['A']); add_outer(c * st['Bb'], cj['A'], ck['B'])
    add_outer(c * st['Ba'], cj['B'], ck['A']); add_outer(c * st['Cab'], cj['B'], ck['B'])
    # pair term pieces
    rho = cj['rho']; rhop = ck['rho']
    add_outer(wb / (4 * st['Da']), rho, [wb * rhop[i] + ck['h'][i] for i in range(n + 1)])
    add_outer(wa / (4 * st['Db']), [wa * rho[i] + cj['h'][i] for i in range(n + 1)], rhop)
    add_outer(wa * wb / 4 * st['gam0'], cj['T'], ck['T'])
    # rho rho'/2 - L L'/2 - b b'/cp
    add_outer(F(1) / 2, rho, rhop)
    add_outer(-F(1) / 2, cj['L'], ck['L'])
    add_outer(-1 / st['cp'], cj['b'], ck['b'])
    return M, st


def metric(x, m, sig, F=float):
    """Gram matrix of N(y) = rho^2 + E Z^2 - (E Z)^2/(sigma + 1) in coordinates (rho, Z_i)."""
    n = len(x)
    N = [[F(0)] * (n + 1) for _ in range(n + 1)]
    N[0][0] = F(1)
    for i in range(n):
        N[i + 1][i + 1] += m[i]
        for j in range(n):
            N[i + 1][j + 1] -= m[i] * m[j] / (sig + 1)
    return N
