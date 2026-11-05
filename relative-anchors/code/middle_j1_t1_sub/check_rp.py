"""Exact checks of Proposition RP (NOTES.md, Section 3): the R- and P-states of the type I block T = {t} of (R_inf)_1.

(a) Closed form.  At rational weight vectors with |w| = 1 (j1lib.py), with Pi_R = Pi - Y_RP^T D_R^{-1} Y_RP,
        pi^T Pi_R pi = sum_{c in U} a_c x_c e_c^2 + b_t e_t^2,   e_r = pi_r/w_r - sum_s w_s pi_s,
    a_c = 2 - tau^2 lam_c^2/8, b_t = 2 sigma (1 - sigma)/(2 - sigma) - (2 - tau)^2 sigma/(4 (1 - sigma)^2),
    as quadratic forms (polarisation), and sum_U x_c e_c + sigma e_t = 0 identically.
(b) The two scalar bounds of the proof, a_c >= A(sigma) and (2 - tau) <= 2 (1 - sigma^2)/(1 + 2 sigma - 2 sigma^2), and
    the final inequality, at the same points; the (R, P) block from j1lib with P_t deleted is positive definite and
    (0, w) is in the kernel of the undeleted block.
(c) The one-variable inequality: F(s) = s A(s) - (1 - s) br(s) equals 3 s p(s)/(2 (2 - s)(2 s^2 - 2 s - 1)^2) with
    p(s) = 5 s^3 - 12 s^2 + 6 s + 2, whose Bernstein coefficients on [0, 1] are 2, 4, 2, 1 (sympy, exact).
(d) Negative control: with the B-part scaled by 1001/1000, the (R, P) block with P_t deleted is not positive definite
    at explicit rational points next to the corner sigma -> 0 with one dominant copy of U (where a_c -> 0), while the
    true one is.
Usage: python check_rp.py"""
import random
from fractions import Fraction as Fr
import sympy as sp
from j1lib import J1, block_states

random.seed(11)


def unit_weights(k):
    v = [Fr(random.randint(1, 9), random.randint(1, 9)) for _ in range(k - 1)]
    n2 = sum(a * a for a in v)
    w = [abs(2 * a / (1 + n2)) for a in v] + [abs((1 - n2) / (1 + n2))]
    if any(a == 0 for a in w):
        return unit_weights(k)
    return w


def ldl_pd(A):
    A = [row[:] for row in A]
    n = len(A)
    for c in range(n):
        if A[c][c] <= 0:
            return False
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / A[c][c]
                for j in range(c + 1, n):
                    A[i][j] -= f * A[c][j]
    return True


def rp_block(w, t, bscale=1):
    k = len(w)
    st = block_states(k, (t,))
    M = J1(w, bscale=bscale).matrix(st)
    R = sorted([i for i, (a, y) in enumerate(st) if a is not None and a[0] == 2],
               key=lambda i: [c for c in range(k) if st[i][1][c] == 1][0])
    P = sorted([i for i, (a, y) in enumerate(st) if a is not None and a[0] == 1], key=lambda i: st[i][0][1])
    idx = R + P
    return [[M[i][j] for j in idx] for i in idx], k


ok_all = True


def report(name, ok):
    global ok_all
    ok_all &= ok
    print(f'{name}: {"PASS" if ok else "FAIL"}', flush=True)


oka = okb = True
cnt = 0
for k in range(2, 9):
    for trial in range(3):
        w = unit_weights(k)
        t = (trial + 2 * k) % k
        x = [a * a for a in w]
        sig = x[t]
        U = [c for c in range(k) if c != t]
        Mb, _ = rp_block(w, t)
        DR = [Mb[i][i] for i in range(k)]
        assert all(Mb[i][j] == 0 for i in range(k) for j in range(k) if i != j)
        YRP = [[Mb[i][k + j] for j in range(k)] for i in range(k)]
        Pi = [[Mb[k + i][k + j] for j in range(k)] for i in range(k)]
        PiR = [[Pi[i][j] - sum(YRP[c][i] * YRP[c][j] / DR[c] for c in range(k)) for j in range(k)] for i in range(k)]
        g = [xc / (2 * (1 - xc)) for xc in x]
        tau = (1 + g[t]) / (1 + sum(g) / 2)
        a = {c: 2 - tau ** 2 / (1 - x[c]) ** 2 / 8 for c in U}
        bt = 2 * sig * (1 - sig) / (2 - sig) - (2 - tau) ** 2 * sig / (4 * (1 - sig) ** 2)

        def e_of(pi):
            m1 = sum(w[s] * pi[s] for s in range(k))
            return [pi[r] / w[r] - m1 for r in range(k)]

        def form(pi):
            e = e_of(pi)
            return sum(a[c] * x[c] * e[c] ** 2 for c in U) + bt * e[t] ** 2
        good = True
        for i in range(k):
            for j in range(k):
                ei = [Fr(1) if r == i else Fr(0) for r in range(k)]
                ej = [Fr(1) if r == j else Fr(0) for r in range(k)]
                eij = [p + q for p, q in zip(ei, ej)]
                Bij = (form(eij) - form(ei) - form(ej)) / 2 if i != j else form(ei)
                good &= Bij == PiR[i][j]
                e = e_of(ei)
                good &= sum(x[c] * e[c] for c in U) + sig * e[t] == 0
        oka &= good
        # (b)
        Asig = 2 - (2 - sig) ** 2 / (2 * (1 + 2 * sig - 2 * sig ** 2) ** 2)
        br = (1 + sig) ** 2 / (1 + 2 * sig - 2 * sig ** 2) ** 2 - 2 * (1 - sig) / (2 - sig)
        okb &= all(a[c] >= Asig > 0 for c in U) if U else True
        okb &= (2 - tau) <= 2 * (1 - sig ** 2) / (1 + 2 * sig - 2 * sig ** 2)
        okb &= bt >= -sig * br
        okb &= sig * Asig - (1 - sig) * br > 0
        keep = [i for i in range(2 * k) if i != k + t]
        okb &= ldl_pd([[Mb[i][j] for j in keep] for i in keep])
        vec = [Fr(0)] * k + list(w)
        okb &= all(sum(Mb[i][j] * vec[j] for j in range(2 * k)) == 0 for i in range(2 * k))
        cnt += 1
print(f'{cnt} weight vectors with |w| = 1, k = 2..8')
report('(a) closed form of Pi_R and the constraint, exactly', oka)
report('(b) bounds a_c >= A(sigma), 2 - tau bound, b_t >= -sigma br, final inequality; (R,P) block PD with P_t deleted', okb)

S = sp.Symbol('s')
br = (1 + S) ** 2 / (1 + 2 * S - 2 * S ** 2) ** 2 - 2 * (1 - S) / (2 - S)
A = 2 - (2 - S) ** 2 / (2 * (1 + 2 * S - 2 * S ** 2) ** 2)
F = sp.together(S * A - (1 - S) * br)
p = 5 * S ** 3 - 12 * S ** 2 + 6 * S + 2
target = 3 * S * p / (2 * (2 - S) * (2 * S ** 2 - 2 * S - 1) ** 2)
okc = sp.simplify(F - target) == 0
coeffs = sp.Poly(p, S).all_coeffs()[::-1]          # a0, a1, a2, a3
a0, a1, a2, a3 = coeffs
bern = [a0, a0 + sp.Rational(1, 3) * a1, a0 + sp.Rational(2, 3) * a1 + sp.Rational(1, 3) * a2, a0 + a1 + a2 + a3]
print(f'  F(s) = 3 s p(s)/(2 (2 - s)(2 s^2 - 2 s - 1)^2): {okc}; Bernstein coefficients of p on [0, 1]: {bern}')
okc &= all(b > 0 for b in bern) and bern == [2, 4, 2, 1]
report('(c) one-variable inequality F > 0 on (0, 1), exact', okc)

okd = True
for w, t in [([Fr(1, 100), Fr(1), Fr(1, 100)], 0), ([Fr(1, 60), Fr(1), Fr(1, 50), Fr(1, 70)], 0),
             ([Fr(1), Fr(1, 80), Fr(1, 90), Fr(1, 100), Fr(30)], 1)]:
    k = len(w)
    res = []
    for bs in (1, Fr(1001, 1000)):
        Mb, _ = rp_block(w, t, bscale=bs)
        keep = [i for i in range(2 * k) if i != k + t]
        res.append(ldl_pd([[Mb[i][j] for j in keep] for i in keep]))
    print(f'  w = {[str(v) for v in w]}, t = {t}: true (R,P) block with P_t deleted PD {res[0]}; scaled PD {res[1]}')
    okd &= res[0] and not res[1]
report('(d) negative controls (B scaled by 1001/1000) at 3 explicit points', okd)
print('ALL PASS' if ok_all else 'SOME FAIL')
