"""K(t, s) = K0(t) + e K1(t, s): signs of the two parts separately; H at mixtures theta e delta_xT + (1-theta) e delta_c."""
import numpy as np
rng = np.random.default_rng(4)
w0 = (1e9, None); w1 = (1e9, None); mixbad = 0; w_mix = 1e9
for trial in range(20000):
    v = rng.exponential(1, 4) * np.exp(rng.normal(0, 2, 4))
    if rng.random() < 0.3: v[2] = 0
    v = v / v.sum(); a, b, xS, e = v
    if min(a, b, e) <= 0: continue
    xT = a + b + xS; c = 2 - xT
    Gi = a / (c - a + b); Gj = b / (c - b + a); gi = a / (c - a); gj = b / (c - b)
    dl = (1 + gi) * (1 + gj); ep = (1 + Gi) * (1 + Gj)
    ts = xT + (c - xT) * np.linspace(0, 1, 50) ** 2
    T, S = np.meshgrid(ts, ts, indexing='ij')
    K0 = dl * (2 + Gi + Gj) / ts ** 2 + 2 * dl / (ts + a + b) ** 2 - ep * (2 + gj) / (ts + a) ** 2 - ep * (2 + gi) / (ts + b) ** 2
    K1 = dl / ((S + a + b) * T ** 2) + dl / (S * (T + a + b) ** 2) - ep / ((S + b) * (T + a) ** 2) - ep / ((S + a) * (T + b) ** 2)
    m0 = K0.min() / (a * b); m1 = K1.min() / (a * b)
    if m0 < w0[0]: w0 = (m0, (a, b, xS, e))
    if m1 < w1[0]: w1 = (m1, (a, b, xS, e))
print('min K0/(ab)', w0); print('min K1/(ab)', w1)
