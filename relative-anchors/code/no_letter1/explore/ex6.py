"""Pointwise sign of the symmetrised 4-kernel of the homogenised numerator N_H (degree 4 in the measure)."""
import numpy as np, itertools
rng = np.random.default_rng(6)
perms = list(itertools.permutations(range(4)))
worst = (1e9, None); negcount = 0; tot = 0
for trial in range(4000):
    v = rng.exponential(1, 4) * np.exp(rng.normal(0, 2, 4))
    if rng.random() < 0.3: v[2] = 0
    v = v / v.sum(); p, q, r, e = v
    if min(p, q, e) <= 0: continue
    xT = p + q + r; c = xT + 2 * e
    gi = p / (c - p); gl = q / (c - q); Gi = p / (c - p + q); Gl = q / (c - q + p)
    dl = (1 + gi) * (1 + gl); ep = (1 + Gi) * (1 + Gl)
    for z in range(10):
        ts = xT + np.exp(rng.normal(np.log(e + 1e-300), 2, 4)) * (rng.random(4) < 0.8)
        if rng.random() < 0.3: ts = xT + 2 * e * rng.random(4)
        Xv = 2 / e + 1 / ts; Yv = (2 + Gi + Gl) / e + 1 / (ts + p + q); Zv = (2 + gi) / e + 1 / (ts + p); Wv = (2 + gl) / e + 1 / (ts + q)
        acc = 0.0
        for s in perms:
            i, j, k, l = s
            acc += 2 * q * Zv[i] * Zv[j] * Wv[k] * Yv[l] - (q * (1 + gi) ** 2 * Xv[i] * Wv[j] * Yv[k] + (q + r) * dl * Xv[i] * Yv[j] * Zv[k] - (q + r) * ep * Zv[i] * Zv[j] * Wv[k]) / e
        acc /= 24
        # natural scale: N_H(e delta_tau) / e^4 at tau = mean
        tot += 1
        if acc < 0: negcount += 1
        sc = q * (c - p) * np.prod(np.minimum(Xv, 1e300)) ** 0  # placeholder scale
        val = acc / (q * (c - p) + 1e-300) * e ** 4 / max(1e-300, (xT + e) ** 0)
        if val < worst[0]: worst = (val, (p, q, r, e), ts)
print(tot, 'kernel evaluations, negative', negcount, 'worst', worst)
