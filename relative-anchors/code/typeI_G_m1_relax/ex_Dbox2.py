"""Exploration (floating point): D with N, N0 at their upper bounds N = g/c, N0 = g0/c0 (c = S/2 + sigma = 1,
c0 = (sigma + S)/2 = b), and g' in [0, S'], g0' in [0, S'/b] (the relaxed box), normalisation K = 1.
Smallest Jacobi-normalised eigenvalue over random points and vertices."""
import pickle, random
import numpy as np
P = pickle.load(open('D_natural.pkl', 'rb'))
names = P['names']
def evald(d, vals):
    tot = 0.0
    for m, c in d.items():
        t = float(int(c.p)) / float(int(c.q)) if hasattr(c, 'p') else float(c)
        for v, e in zip(vals, m):
            if e:
                t *= v ** int(e)
        tot += t
    return tot
def Dt(x, t, th, th0):
    S = 2 * (1 - t); sg = t; b = 1 - t / 2
    g = x / (1 + x) + th * (S - x)
    g0 = x / (x + b) + th0 * (S - x) / b
    N = g; N0 = g0 / b
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
rng = random.Random(2)
worst = (9, None)
for it in range(40000):
    t = rng.uniform(1e-4, 1 - 1e-4) if it % 3 else rng.choice([1e-4, 1 - 1e-4, 0.5])
    S = 2 * (1 - t)
    xi = rng.uniform(1e-4, 1.0) if it % 5 else rng.choice([1e-4, 1.0])
    th = rng.random() if it % 2 else float(rng.random() < 0.5)
    th0 = rng.random() if it % 2 else float(rng.random() < 0.5)
    M = Dt(S * xi, t, th, th0)
    if np.any(np.diag(M) <= 0):
        worst = (-1, (t, xi, th, th0)); break
    dsc = 1 / np.sqrt(np.diag(M))
    lam = np.linalg.eigvalsh(M * np.outer(dsc, dsc)).min()
    if lam < worst[0]:
        worst = (lam, (round(t, 5), round(xi, 5), round(th, 3), round(th0, 3)))
print('relaxed box, N and N0 at their maxima: smallest Jacobi-normalised eigenvalue', worst)
