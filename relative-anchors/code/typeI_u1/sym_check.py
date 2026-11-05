"""SYMBOLIC identity check of the closed forms (gstruct.pkl) against u1form.phi_matrix for m = 3 (k = 4), with
symbolic weights omega = sqrt(sigma), w_t = sqrt(s_t) (s_1, s_2, s_3, sigma independent positive symbols).

Blocks checked: Gt_11 (diagonal, copy 1) and Gt_12 (off-diagonal, copies 1, 2), both in scaled coordinates
Gt = D G D, D_t = diag(1, w_t omega).  With m = 3 the third weight is free, so S_T = s_1 + s_2 + s_3 is an
independent variable and the check covers the full four-parameter family (s_t, s_t', sigma, S_T).

Why m = 3 suffices for every m: in u1form.phi_matrix every entry of Phi in the rows/columns (rho_t, zeta_t, rho_t',
zeta_t', q) is a fixed function of (omega, w_t, w_t', S_T): the other copies tau enter only through
sum_tau b_tau[t] b_tau[t'] = (1/4) sum_{tau != t, t'} s_tau gamma[t] gamma[t'] + (terms of t, t'), through
wb = sum_tau w_tau b_tau (coefficients (1/2) w_t alpha_t + (1/2) S_T gamma), and through e, K (functions of S_T, sigma).
"""
import sympy as sp, pickle, time, sys
from u1form import phi_matrix

t0 = time.time()
m = int(sys.argv[1]) if len(sys.argv) > 1 else 3
sig = sp.Symbol('sigma', positive=True)
ss = sp.symbols('s1:%d' % (m + 1), positive=True)
w = [sp.sqrt(sig)] + [sp.sqrt(x) for x in ss]


def F(x, y=None):
    if y is not None:
        return sp.Rational(x, y)
    return sp.Integer(x) if isinstance(x, int) else x


P = phi_matrix(w, F=F)
n = 2 * m + 1
iq = 2 * m
Ps = lambda i, j: (P[i][j] + P[j][i]) / 2
idx = [0, m, 1, m + 1]          # rho_1, zeta_1, rho_2, zeta_2
Pqq = sp.cancel(sp.together(Ps(iq, iq)))
print('Phi_qq =', sp.factor(Pqq), f'({time.time() - t0:.1f} s)', flush=True)
Go, Gd = pickle.load(open('gstruct.pkl', 'rb'))
s, s2, S, sg = sp.symbols('s s2 S sigma', positive=True)
ST = sum(ss)
sub_off = {s: ss[0], s2: ss[1], S: ST, sg: sig}
sub_dg = {s: ss[0], S: ST, sg: sig}
D1 = [1, w[1] * w[0]]
D2 = [1, w[2] * w[0]]
ok = True
for a in range(2):
    for b in range(2):
        i, j = idx[a], idx[b]
        Gij = Ps(i, j) - Ps(i, iq) * Ps(iq, j) / Pqq
        val = D1[a] * Gij * D1[b]
        diff = sp.cancel(sp.together(val - Gd[a, b].subs(sub_dg)))
        print(f'diag ({a},{b}): difference = {diff}  ({time.time() - t0:.1f} s)', flush=True)
        ok &= diff == 0
        i2 = idx[2 + b]
        Gij = Ps(i, i2) - Ps(i, iq) * Ps(iq, i2) / Pqq
        val = D1[a] * Gij * D2[b]
        diff = sp.cancel(sp.together(val - Go[a, b].subs(sub_off)))
        print(f'off  ({a},{b}): difference = {diff}  ({time.time() - t0:.1f} s)', flush=True)
        ok &= diff == 0
print(f'm = {m}: symbolic identity of all 8 entries: {ok}  ({time.time() - t0:.1f} s)')
print('PASS' if ok else 'FAIL')
