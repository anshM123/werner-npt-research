"""Exact checks and negative controls for Theorem D (condition (i), NOTES.md Section 12).

(a) At rational weight vectors (u = 3..8): every D(s_h) of the reduced system (reduced_c.py) is > 0; the parameters
    (t, xi, th, th0) of prove_D.py computed from the weights lie in the certified domain; and D(s_h) - Dt(s_h) is a
    diagonal positive semidefinite matrix, where Dt is D_natural.pkl evaluated at N = g/c, N0 = g0/c0.
(b) Negative controls at explicit rational points of the domain: Dt - (1/10) diag(Dt) is not positive definite (the
    margin of the certificate is below 1/10 there), and D evaluated at N = 2 g/c, N0 = 2 g0/c0 (a false bound) is not
    positive definite."""
import random, pickle
import flint
from reduced_c import ReducedC
from certify_points import family

Q = flint.fmpq
P = pickle.load(open('D_natural.pkl', 'rb'))
names = P['names']


def evald(d, vals):
    tot = Q(0)
    for m, c in d.items():
        t = Q(c)
        for v, e in zip(vals, m):
            if e:
                t *= v ** int(e)
        tot += t
    return tot


def Dnat(sv, sg, SS, g, N, g0, N0):
    dt0 = (3 * sg + SS) / 2
    vals = dict(s=sv, sigma=sg, S=SS, IG=1 / (2 + g), IGT=1 / (2 + sg / dt0 + g0), g=g, N=N, g0=g0, N0=N0, u=Q(1))
    vv = [vals[n] for n in names]
    M = [[Q(0)] * 3 for _ in range(3)]
    idx = {'ah': 0, 'ch': 1, 'zh': 2}
    for (a, b), (nd, dd) in P['D'].items():
        val = evald(nd, vv) / evald(dd, vv)
        i, j = idx[a], idx[b]
        if i == j:
            M[i][i] += val
        else:
            M[i][j] += val / 2
            M[j][i] += val / 2
    return M


def posdef3(M):
    m1 = M[0][0]
    m2 = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    m3 = flint.fmpq_mat(3, 3, [M[i][j] for i in range(3) for j in range(3)]).det()
    return m1 > 0 and m2 > 0 and m3 > 0


R = ReducedC()
rng = random.Random(17)
ok = True
cnt = 0
for u in range(3, 9):
    for k in range(4):
        wt, wU = family(u, k, rng)
        num = R.numeric(wt, wU)
        sv = num['sv']
        sg = num['sg']
        SS = sum(sv, Q(0))
        c = SS / 2 + sg
        c0 = (sg + SS) / 2
        g = sum((x / (x + c) for x in sv), Q(0))
        N = sum((x / (x + c) ** 2 for x in sv), Q(0))
        g0 = sum((x / (x + c0) for x in sv), Q(0))
        N0 = sum((x / (x + c0) ** 2 for x in sv), Q(0))
        for h in range(u):
            Dh = num['D'][h]
            a1 = posdef3(Dh)
            # D_natural at the true sums equals D
            Dn = Dnat(sv[h], sg, SS, g, N, g0, N0)
            a0 = all(Dn[i][j] == Dh[i][j] for i in range(3) for j in range(3))
            # parameters (K = 1 normalisation: divide all squared weights by c)
            t = sg / c
            x = sv[h] / c
            Sn = SS / c
            b = 1 - t / 2
            xi = x / Sn
            gp = g - x / (1 + x)
            g0p = g0 - x / (x + b)
            th = gp / (Sn - x) if Sn != x else Q(0)
            th0 = g0p * b / (Sn - x) if Sn != x else Q(0)
            a2 = (0 < t < 1) and (0 < xi <= 1) and (0 <= th <= 1) and (0 <= th0 <= 1)
            Dt = Dnat(sv[h], sg, SS, g, g / c, g0, g0 / c0)
            diff = [[Dh[i][j] - Dt[i][j] for j in range(3)] for i in range(3)]
            a3 = all(diff[i][j] == 0 for i in range(3) for j in range(3) if i != j) and all(diff[i][i] >= 0
                                                                                        for i in range(3))
            ok &= a0 and a1 and a2 and a3
            cnt += 1
print(f'(a) {cnt} blocks D(s_h) at 24 weight vectors (u = 3..8): D_natural = D, D > 0, parameters in the domain, '
      f'D - Dt diagonal >= 0: {ok}', flush=True)


def Dt_at(t, xi, th, th0, factorN=1):
    sg = t
    SS = 2 * (1 - t)
    x = SS * xi
    b = 1 - t / 2
    g = x / (1 + x) + th * (SS - x)
    g0 = x / (x + b) + th0 * (SS - x) / b
    return Dnat(x, sg, SS, g, factorN * g, g0, factorN * g0 / b)


ctl1 = False
ctl2 = False
for t in [Q(9999, 10000), Q(99, 100), Q(9, 10), Q(1, 2)]:
    for xi in [Q(1), Q(9, 10), Q(1, 2), Q(1, 10)]:
        for th in [Q(0), Q(1, 2), Q(7, 10), Q(1)]:
            for th0 in [Q(0), Q(4, 5), Q(1)]:
                M = Dt_at(t, xi, th, th0)
                Mm = [[M[i][j] - (M[i][j] / 10 if i == j else 0) for j in range(3)] for i in range(3)]
                if not ctl1 and not posdef3(Mm):
                    ctl1 = True
                    print(f'(b1) Dt - diag(Dt)/10 not positive definite at t = {t}, xi = {xi}, th = {th}, th0 = {th0}',
                          flush=True)
                M2 = Dt_at(t, xi, th, th0, factorN=2)
                if not ctl2 and not posdef3(M2):
                    ctl2 = True
                    print(f'(b2) D at N = 2g/c, N0 = 2g0/c0 not positive definite at t = {t}, xi = {xi}, th = {th}, '
                          f'th0 = {th0}', flush=True)
print(f'controls found: (b1) {ctl1}, (b2) {ctl2}')
print('PASS' if ok else 'FAIL')
