"""Float exploration: dependence of the pair quantity D (inequality (i)) and of the reduced row quantity H on the
U-configuration, at fixed (p, q, r, e).  Measures with the mass-position link: atoms m_u (sum e) at t_u = c - 2 m_u."""
import numpy as np
rng = np.random.default_rng(1)

def Phi(ms, c, d):
    ms = np.asarray(ms, float)
    return np.sum(ms / (c - 2 * ms + d))

def pieces(p, q, r, ms):
    xT = p + q + r; c = 2 - xT
    Gi = p / (c - p + q); Gj = q / (c - q + p)
    gi = p / (c - p); gj = q / (c - q)
    T0 = Phi(ms, c, 0); P = Phi(ms, c, p); Q = Phi(ms, c, q); R = Phi(ms, c, p + q)
    X = 2 + T0; Y = 2 + Gi + Gj + R; Z = 2 + gi + P; W = 2 + gj + Q
    dl = (1 + gi) * (1 + gj); ep = (1 + Gi) * (1 + Gj)
    D = dl * X * Y - ep * Z * W
    Fii = 2 - X * (1 + gi) ** 2 / Z ** 2
    phi = D / (Z * W * Y)
    H = q * Fii - (q + r) * phi
    return D, H, Fii, phi

def configs(e, N=8):
    out = {'one': [e]}
    for n in (2, 3, 5, 20, 200):
        out[f'eq{n}'] = [e / n] * n
    for s in (0.1, 0.3, 0.5):
        out[f'two{s}'] = [e * s, e * (1 - s)]
    for z in range(6):
        w = rng.exponential(1, rng.integers(2, N)) ** 2; out[f'rnd{z}'] = list(e * w / w.sum())
    out['tiny'] = [e / 5000] * 5000
    return out

stats = {}
for trial in range(3000):
    v = rng.exponential(1, 4) * np.exp(rng.normal(0, 2, 4))
    if rng.random() < 0.3: v[2] = 0
    v = v / v.sum(); p, q, r, e = v
    if min(p, q, e) <= 0: continue
    res = {name: pieces(p, q, r, ms) for name, ms in configs(e).items()}
    Dn = {k: val[0] / (p * q) for k, val in res.items()}
    Hn = {k: val[1] for k, val in res.items()}
    kD = min(Dn, key=Dn.get); kH = min(Hn, key=Hn.get)
    stats.setdefault('D', {}).setdefault(kD, 0); stats['D'][kD] += 1
    stats.setdefault('H', {}).setdefault(kH, 0); stats['H'][kH] += 1
    for nm, d in (('Dmin', min(Dn.values())), ('Hmin', min(Hn.values()))):
        if nm not in stats or d < stats[nm][0]:
            stats[nm] = (d, (p, q, r, e), kD if nm == 'Dmin' else kH)
print(stats)
