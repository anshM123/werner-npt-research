"""K0(t) >= 0 and K1(t, s) >= 0 for t, s in (0, inf)?  and for t, s >= xT?"""
import numpy as np
rng = np.random.default_rng(5)
res = {}
for trial in range(20000):
    v = rng.exponential(1, 4) * np.exp(rng.normal(0, 2, 4))
    if rng.random() < 0.3: v[2] = 0
    v = v / v.sum(); a, b, xS, e = v
    if min(a, b, e) <= 0: continue
    xT = a + b + xS; c = 2 - xT
    Gi = a / (c - a + b); Gj = b / (c - b + a); gi = a / (c - a); gj = b / (c - b)
    dl = (1 + gi) * (1 + gj); ep = (1 + Gi) * (1 + Gj)
    for name, ts in (('all_t', np.exp(np.linspace(-14, 8, 120))), ('t>=xT', xT + np.exp(np.linspace(-14, 8, 120)))):
        T, S = np.meshgrid(ts, ts, indexing='ij')
        K0 = dl * (2 + Gi + Gj) / ts ** 2 + 2 * dl / (ts + a + b) ** 2 - ep * (2 + gj) / (ts + a) ** 2 - ep * (2 + gi) / (ts + b) ** 2
        K1 = dl / ((S + a + b) * T ** 2) + dl / (S * (T + a + b) ** 2) - ep / ((S + b) * (T + a) ** 2) - ep / ((S + a) * (T + b) ** 2)
        for nm, val in ((name + ' K0', (K0 * ts ** 2).min() / (a * b)), (name + ' K1', (K1 * T ** 2 * S).min() / (a * b))):
            if nm not in res or val < res[nm][0]:
                res[nm] = (val, (a, b, xS, e))
for k, v in res.items(): print(k, v)
