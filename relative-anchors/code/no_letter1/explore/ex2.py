"""Test K(t, s) >= 0: the derivative of D in the position t of an atom (sign reversed), with the measure nu = e delta_s
(extreme points).  If K >= 0 on [xT, c]^2 then D(mu) >= D(e delta_c) for all measures of mass e on [xT, c]."""
import numpy as np
rng = np.random.default_rng(2)
worst = (1e9, None); cnt = 0
for trial in range(20000):
    v = rng.exponential(1, 4) * np.exp(rng.normal(0, 2, 4))
    if rng.random() < 0.3: v[2] = 0
    v = v / v.sum(); a, b, xS, e = v
    if min(a, b, e) <= 0: continue
    xT = a + b + xS; c = 2 - xT
    Gi = a / (c - a + b); Gj = b / (c - b + a); gi = a / (c - a); gj = b / (c - b)
    dl = (1 + gi) * (1 + gj); ep = (1 + Gi) * (1 + Gj)
    ts = xT + (c - xT) * np.linspace(0, 1, 60) ** 2
    T, S = np.meshgrid(ts, ts, indexing='ij')
    X = 2 + e / S; Y = 2 + Gi + Gj + e / (S + a + b); Z = 2 + gi + e / (S + a); W = 2 + gj + e / (S + b)
    K = dl * Y / T ** 2 + dl * X / (T + a + b) ** 2 - ep * W / (T + a) ** 2 - ep * Z / (T + b) ** 2
    # natural scale: everything ~ ab?  normalise by a*b
    m = K.min() / (a * b)
    cnt += 1
    if m < worst[0]:
        i, j = np.unravel_index(np.argmin(K), K.shape)
        worst = (m, (a, b, xS, e), ts[i], ts[j])
print(cnt, 'configs; min K/(ab):', worst)
