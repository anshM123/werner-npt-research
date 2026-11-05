import numpy as np
from explore10 import entries
from rowbounds import bounds
beta = 12 / 5
def exact_rows(x, t):
    Mc, Mca, q = entries(x, t)
    w, d, rho, g, eta2, N, Lam, e, m, b, B, Pm = [q[k] for k in ('w', 'd', 'rho', 'g', 'eta2', 'N', 'Lam', 'e', 'm', 'b', 'B', 'Pm')]
    u = len(x); out = []
    for h in range(u):
        R = 0; R1 = R3 = R4 = 0
        for k in range(u):
            if k == h: continue
            P1 = 0.5 - 1 / ((2 + g) * d[h] * d[k])
            P3 = 0.25 * (B + 2 * b[h] + 2 * b[k] + Pm - m[h] ** 2 / Lam[h] - m[k] ** 2 / Lam[k])
            P4 = 0.25 * (m[h] * (2 * e[h] - m[h]) / Lam[h] + m[k] * (2 * e[k] - m[k]) / Lam[k])
            R1 += x[k] * P1; R3 += x[k] * P3; R4 += x[k] * abs(P4)
            pi = (1 / ((2 + g) * d[h] * d[k]) - 0.5 - 0.25 * (B + 2 * b[h] + 2 * b[k] + Pm - 2 * e[h] * m[h] / Lam[h] - 2 * e[k] * m[k] / Lam[k]))
            R += x[k] * abs(pi)
        out.append(dict(Mchh=Mc[h, h], R=t * R, R1=t * R1, R3=t * R3, R4=t * R4, B=B, bh=b[h], Pm=Pm, N=N, g=g))
    return out
S = 1.025; t = 1 - S / 2
for xs in ([1e-4, S - 1e-4], [1e-4] + [(S - 1e-4) / 20] * 20):
    x = np.array(xs)
    r = exact_rows(x, t)[0]
    print('config', 'one big' if len(xs) == 2 else 'many', {k: round(v, 4) for k, v in r.items()})
    print('   exact LB/(1-t):', (r['Mchh'] - r['R']) / (1 - t), ' triangle LB/(1-t):', (r['Mchh'] - r['R1'] - r['R3'] - r['R4']) / (1 - t))
    xi = x[0] / S; Sp = S - x[0]; rSp = Sp / (1 + Sp); gp = r['g'] - x[0] / (1 + x[0])
    theta = (gp - rSp) / (Sp - rSp)
    for nu in (0.0, 1.0):
        bb = bounds(S, xi, theta, nu, one=1.0)
        print('   bound nu=%g: LBc %.4f' % (nu, bb['LBc_a'] if bb['mh'] >= 0 else bb['LBc_b']), ' theta', round(theta, 3))
