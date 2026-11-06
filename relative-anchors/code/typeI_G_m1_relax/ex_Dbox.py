"""Exploration (floating point from exact formulas; no claim): the local block D in natural sums (D_natural.pkl),
normalised K = sigma + S/2 = 1 (t = sigma, S = 2(1 - t), x = s_h, b = 1 - t/2), with the other copies entering
through g' in [rho(S'), S'], N' in [g'/(1+S'), g'], g0' in [rho0(S'), S'/b], N0' in [(g0'/b)(1 - rho0(S')), g0'/b]
(S' = S - x).  Reports the coefficients of N and N0 (sign), and the smallest Jacobi-normalised eigenvalue of D over
random points of this box and its vertices."""
import pickle, random, itertools
import numpy as np
P = pickle.load(open('D_natural.pkl', 'rb'))
names = P['names']
def evald(d, vals):
    tot = 0.0
    for m, c in d.items():
        t = float(int(c.p)) / float(int(c.q)) if hasattr(c, "p") else float(c)
        for v, e in zip(vals, m):
            if e:
                t *= v ** int(e)
        tot += t
    return tot
def Dm(x, t, gp, Np, g0p, N0p):
    S = 2 * (1 - t); sg = t; b = 1 - t / 2
    rho = lambda y: y / (1 + y)
    rho0 = lambda y: y / (y + b)
    g = rho(x) + gp; N = rho(x) / (1 + x) + Np
    g0 = rho0(x) + g0p; N0 = rho0(x) / (x + b) + N0p
    IG = 1 / (2 + g); IGT = 1 / (2 + sg / (1 + t / 2) + g0)
    vals = dict(s=x, sigma=sg, S=S, IG=IG, IGT=IGT, g=g, N=N, g0=g0, N0=N0, u=1.0)
    vv = [vals[n] for n in names]
    M = np.zeros((3, 3)); idx = {'ah': 0, 'ch': 1, 'zh': 2}
    for (a, b_), (nd, dd) in P['D'].items():
        val = evald(nd, vv) / evald(dd, vv)
        i, j = idx[a], idx[b_]
        if i == j:
            M[i, i] += val
        else:
            M[i, j] += val / 2; M[j, i] += val / 2
    return M
rng = random.Random(1)
worst = (9, None)
for it in range(20000):
    t = rng.uniform(0.001, 0.999)
    S = 2 * (1 - t); x = S * rng.uniform(0.0005, 1.0); Sp = S - x; b = 1 - t / 2
    rho = lambda y: y / (1 + y); rho0 = lambda y: y / (y + b)
    th = [rng.random() for _ in range(4)]
    if it % 4 == 0:
        th = [float(rng.random() < 0.5) for _ in range(4)]
    gp = rho(Sp) + th[0] * (Sp - rho(Sp))
    Np = gp / (1 + Sp) + th[1] * (gp - gp / (1 + Sp))
    g0p = rho0(Sp) + th[2] * (Sp / b - rho0(Sp))
    N0p = (g0p / b) * (1 - rho0(Sp)) + th[3] * (g0p / b - (g0p / b) * (1 - rho0(Sp)))
    M = Dm(x, t, gp, Np, g0p, N0p)
    dsc = 1 / np.sqrt(np.abs(np.diag(M)))
    lam = np.linalg.eigvalsh(M * np.outer(dsc, dsc)).min()
    if np.any(np.diag(M) <= 0):
        lam = -1
    if lam < worst[0]:
        worst = (lam, (round(t, 4), round(x / S, 4), [round(v, 3) for v in th]))
print('smallest Jacobi-normalised eigenvalue of D over the box:', worst)
# coefficient of N in D_aa, D_cc and of N0 in D_zz: finite differences (D is affine in N, N0)
for it in range(5):
    t = rng.uniform(0.01, 0.99); S = 2 * (1 - t); x = S * rng.uniform(0.01, 1); Sp = S - x; b = 1 - t / 2
    gp = Sp * 0.5; g0p = Sp * 0.6
    M1 = Dm(x, t, gp, 0.0, g0p, 0.0); M2 = Dm(x, t, gp, 1.0, g0p, 1.0)
    print('d/dN of (D_aa, D_cc, D_zz) at a random point:', np.round(np.diag(M2 - M1), 6),
          ' off-diagonal change:', np.round(M2 - M1 - np.diag(np.diag(M2 - M1)), 9).max())
