"""SYMBOLIC identity check of the closed forms (gstruct.pkl) against u1form.phi_matrix, in the rational function
field Q(w_0, ..., w_m) (sparse SymPy FracField arithmetic; u1form.phi_matrix is run unchanged on field elements).

Checked, for each m given on the command line (default 1 2 3 4): Phi_qq = (3 S^2 + 2 S sigma + sigma^2)/(2 K^2), and
the scaled blocks Gt_tu = D_t G_tu D_u (D_t = diag(1, w_t w_0)) for (t, u) in {1, 2, 3}^2 (those that exist) equal
the closed forms Gt_diag(s_t, S_T, sigma), Gt_off(s_t, s_u, S_T, sigma) AS RATIONAL FUNCTIONS of the weights.
For m >= 3 the remaining weights are free, so S_T is independent of (s_t, s_u, sigma) and the full four-parameter
family is covered.  (Why this gives every m: see NOTES.md, Step 1.)
"""
import sys, time, pickle
import sympy as sp
from sympy import QQ
from sympy.polys.fields import field
from u1form import phi_matrix

Go, Gd = pickle.load(open('gstruct.pkl', 'rb'))
s, s2, S, sg = sp.symbols('s s2 S sigma', positive=True)
ms = [int(x) for x in sys.argv[1:]] or [1, 2, 3, 4]
ok_all = True
t0 = time.time()
for m in ms:
    names = ','.join(['w%d' % i for i in range(m + 1)])
    Fld, *W = field(names, QQ)
    Wsym = sp.symbols(names)
    if m == 0:
        continue

    def F(x, y=None):
        if y is not None:
            return Fld(QQ(x, y))
        return x if hasattr(x, 'field') else Fld(x)

    P = phi_matrix(W, F=F)
    iq = 2 * m
    Ps = lambda i, j: (P[i][j] + P[j][i]) / 2
    Pqq = Ps(iq, iq)
    s_ = [w * w for w in W]
    ST = sum(s_[1:], Fld(0)); sig = s_[0]; K = ST + sig / 2
    ok = (Pqq - (3 * ST ** 2 + 2 * ST * sig + sig ** 2) / (2 * K ** 2)) == 0
    STsym = sum(Wsym[i] ** 2 for i in range(1, m + 1))

    def to_field(expr, subs):
        num, den = sp.fraction(sp.together(expr))
        return Fld.from_expr(sp.expand(num.subs(subs))) / Fld.from_expr(sp.expand(den.subs(subs)))

    idx = range(1, min(m, 3) + 1)
    nblk = 0
    for t in idx:
        for u in idx:
            for a in range(2):
                for b in range(2):
                    i = (t - 1) if a == 0 else (m + t - 1)
                    j = (u - 1) if b == 0 else (m + u - 1)
                    val = Ps(i, j) - Ps(i, iq) * Ps(iq, j) / Pqq
                    if a == 1: val = val * W[t] * W[0]
                    if b == 1: val = val * W[u] * W[0]
                    if t == u:
                        ref = to_field(Gd[a, b], {s: Wsym[t] ** 2, S: STsym, sg: Wsym[0] ** 2})
                    else:
                        ref = to_field(Go[a, b], {s: Wsym[t] ** 2, s2: Wsym[u] ** 2, S: STsym, sg: Wsym[0] ** 2})
                    ok &= (val - ref) == 0
            nblk += 1
    print(f"m = {m}: Phi_qq formula and {nblk} blocks (t, u in {list(idx)}) identical as rational functions: {ok} "
          f"({time.time() - t0:.1f} s)", flush=True)
    ok_all &= ok
print('PASS' if ok_all else 'FAIL')
