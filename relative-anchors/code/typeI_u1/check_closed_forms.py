"""EXACT check of the closed forms in gstruct.pkl against u1form.phi_matrix.

For rational weight vectors w = (omega, w_1, ..., w_m) we build Phi exactly (Fractions), symmetrize, eliminate q
exactly (G = Phi_xx - Phi_xq Phi_qx / Phi_qq), and compare the scaled blocks
    Gt_tt' = D_t G_tt' D_t',   D_t = diag(1, w_t * omega)   (coordinates (rho_t, zeta_t))
entry by entry with the rational functions Gt_off(s_t, s_t', S_T, sigma) and Gt_diag(s_t, S_T, sigma).
Equality is tested exactly (Fraction == Fraction).  Also checks Phi_qq = 2 - (S/4) gq^2 - S^2 gq^2/(8 sigma) > 0.
"""
import pickle, sys, time
from fractions import Fraction as Fr
import sympy as sp
from u1form import phi_matrix

t0 = time.time()
Go, Gd = pickle.load(open('gstruct.pkl', 'rb'))
s, s2, S, sg = sp.symbols('s s2 S sigma', positive=True)


def to_frac_fun(expr, vars_):
    num, den = sp.fraction(sp.together(expr))
    pn = sp.Poly(sp.expand(num), *vars_)
    pd = sp.Poly(sp.expand(den), *vars_)
    tn = [(mon, Fr(int(sp.numer(c)), int(sp.denom(c)))) for mon, c in pn.terms()]
    td = [(mon, Fr(int(sp.numer(c)), int(sp.denom(c)))) for mon, c in pd.terms()]

    def ev(terms, x):
        tot = Fr(0)
        for mon, c in terms:
            v = c
            for xi, e in zip(x, mon):
                if e:
                    v *= xi ** e
            tot += v
        return tot

    return lambda *x: ev(tn, x) / ev(td, x)


V4 = (s, s2, S, sg)
V3 = (s, S, sg)
foff = [[to_frac_fun(Go[i, j], V4) for j in range(2)] for i in range(2)]
fdiag = [[to_frac_fun(Gd[i, j], V3) for j in range(2)] for i in range(2)]
print(f"closed forms loaded and converted ({time.time() - t0:.1f} s)", flush=True)


def schur_G(w):
    P = phi_matrix(w, Fr)
    n = len(P)
    P = [[(P[i][j] + P[j][i]) / 2 for j in range(n)] for i in range(n)]
    iq = n - 1
    G = [[P[i][j] - P[i][iq] * P[iq][j] / P[iq][iq] for j in range(iq)] for i in range(iq)]
    return G, P[iq][iq]


import random
random.seed(20261104)
ok = True
count_blocks = 0
vecs = []
for m in range(1, 8):
    for trial in range(4 if m <= 4 else 2):
        w = [Fr(random.randint(1, 15), random.randint(1, 15)) for _ in range(m + 1)]
        vecs.append(w)
# a few extreme ones
vecs.append([Fr(1, 100), Fr(1), Fr(1, 50), Fr(3)])
vecs.append([Fr(40), Fr(1, 30), Fr(2), Fr(1, 7), Fr(5)])
vecs.append([Fr(1), Fr(1000), Fr(1, 1000)])
for w in vecs:
    m = len(w) - 1
    om = w[0]; sig = om * om
    st = [x * x for x in w[1:]]
    ST = sum(st, Fr(0))
    G, Gqq = schur_G(w)
    K = ST + sig / 2
    gq = 2 * om / K
    kappa = 2 - (ST / 4) * gq ** 2 - ST ** 2 * gq ** 2 / (8 * sig)
    ok &= (Gqq == kappa) and kappa > 0 and kappa == (3 * ST ** 2 + 2 * ST * sig + sig ** 2) / (2 * K ** 2)
    D = [[Fr(1), w[1 + t] * om] for t in range(m)]
    for t in range(m):
        for u in range(m):
            idx_t = [t, m + t]; idx_u = [u, m + u]
            for a in range(2):
                for b in range(2):
                    val = D[t][a] * G[idx_t[a]][idx_u[b]] * D[u][b]
                    if t == u:
                        ref = fdiag[a][b](st[t], ST, sig)
                    else:
                        ref = foff[a][b](st[t], st[u], ST, sig)
                    if val != ref:
                        ok = False
                        print('MISMATCH', w, t, u, a, b, float(val), float(ref))
            count_blocks += 1
print(f"exact comparison: {len(vecs)} weight vectors, m = 1..7 (+3 extreme), {count_blocks} 2x2 blocks, "
      f"Phi_qq = kappa > 0 checked ({time.time() - t0:.1f} s)")
print("PASS" if ok else "FAIL")
