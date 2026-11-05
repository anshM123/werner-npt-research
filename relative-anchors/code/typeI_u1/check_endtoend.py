"""EXACT end-to-end check: u1form.phi_matrix  ->  the polynomials FA, FB1, FB2 of polys.py.

For rational weight vectors (omega, w_1, ..., w_m) (m = 1..6, also spread-out weights) we compute G (Schur complement
of the symmetrised Phi onto (rho, zeta)) exactly, and with c^2 = 1 + sigma/S_T, H = 2/5:
  (A) for t != t':  H^2 - ||Ghat_tt'||_F^2 S_T^2/(s_t s_t')            ==  FA/(4 Den^2)
  (B) for every t:  (Ghat_tt)_{zeta zeta} - h                          ==  FB1/(2 Dd)
                    det(Ghat_tt - h I)                                 ==  FB2/(2 Dd)^2,   h = H (1 - s_t/S_T)
where Ghat = diag(1, c) G diag(1, c) blockwise, the right-hand sides are evaluated at the normalised point
(s_t/S_T, s_t'/S_T, sigma/S_T).  Only squares of c and of the off-diagonal entries occur, so everything is rational.
"""
import pickle, random, time
from fractions import Fraction as Fr
import sympy as sp
from u1form import phi_matrix
from polys import s, s2, sg

t0 = time.time()
d = pickle.load(open('polys.pkl', 'rb'))
H = Fr(int(sp.numer(d['H'])), int(sp.denom(d['H'])))


def frac_fun(expr, gens):
    P = sp.Poly(expr, *gens)
    terms = [(m, Fr(int(sp.numer(c)), int(sp.denom(c)))) for m, c in P.terms()]

    def f(*x):
        tot = Fr(0)
        for m, c in terms:
            v = c
            for xi, e in zip(x, m):
                if e:
                    v *= xi ** e
            tot += v
        return tot
    return f


fFA = frac_fun(d['FA'], (s, s2, sg)); fDen = frac_fun(d['Den'], (s, s2, sg))
fFB1 = frac_fun(d['FB1'], (s, sg)); fFB2 = frac_fun(d['FB2'], (s, sg)); fDd = frac_fun(d['Dd'], (s, sg))


def schur_G(w):
    P = phi_matrix(w, Fr)
    n = len(P)
    P = [[(P[i][j] + P[j][i]) / 2 for j in range(n)] for i in range(n)]
    iq = n - 1
    return [[P[i][j] - P[i][iq] * P[iq][j] / P[iq][iq] for j in range(iq)] for i in range(iq)]


random.seed(7)
vecs = []
for m in range(1, 7):
    for _ in range(3):
        vecs.append([Fr(random.randint(1, 20), random.randint(1, 20)) for _ in range(m + 1)])
vecs += [[Fr(1, 200), Fr(3), Fr(1, 90), Fr(7, 3)], [Fr(50), Fr(1, 40), Fr(1), Fr(2, 9), Fr(9)], [Fr(1), Fr(1, 3)]]
ok = True
nA = nB = 0
for w in vecs:
    m = len(w) - 1
    sig = w[0] ** 2
    st = [x * x for x in w[1:]]
    ST = sum(st, Fr(0))
    c2 = 1 + sig / ST
    G = schur_G(w)
    blk = lambda t, u, a, b: G[[t, m + t][a]][[u, m + u][b]]
    for t in range(m):
        # (B)
        h = H * (1 - st[t] / ST)
        g00, g01, g11 = blk(t, t, 0, 0), blk(t, t, 0, 1), blk(t, t, 1, 1)
        lhs1 = c2 * g11 - h
        lhs2 = (g00 - h) * (c2 * g11 - h) - c2 * g01 ** 2
        x, y = st[t] / ST, sig / ST
        ok &= lhs1 == fFB1(x, y) / (2 * fDd(x, y))
        ok &= lhs2 == fFB2(x, y) / (2 * fDd(x, y)) ** 2
        nB += 1
        for u in range(m):
            if u == t:
                continue
            fro2 = blk(t, u, 0, 0) ** 2 + c2 * (blk(t, u, 0, 1) ** 2 + blk(t, u, 1, 0) ** 2) + c2 ** 2 * blk(t, u, 1, 1) ** 2
            lhs = H ** 2 - fro2 * ST ** 2 / (st[t] * st[u])
            x1, x2 = st[t] / ST, st[u] / ST
            ok &= lhs == fFA(x1, x2, y) / (4 * fDen(x1, x2, y) ** 2)
            nA += 1
print(f"exact end-to-end check: {len(vecs)} weight vectors (m = 1..6), {nA} ordered pairs (A), {nB} diagonal blocks (B) "
      f"({time.time() - t0:.1f} s)")
print("PASS" if ok else "FAIL")
