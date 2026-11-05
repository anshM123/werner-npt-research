"""H (reduced row quantity): compare admissible configurations with the two single-atom extremes and with single
atoms e delta_tau at intermediate tau, and with two-atom measures of arbitrary masses (relaxation)."""
import numpy as np
rng = np.random.default_rng(3)

def Hval(p, q, r, e, atoms):  # atoms: list of (mass, position)
    xT = p + q + r; c = 2 - xT
    m = np.array([a[0] for a in atoms]); t = np.array([a[1] for a in atoms])
    Ph = lambda d: np.sum(m / (t + d))
    Gi = p / (c - p + q); Gj = q / (c - q + p); gi = p / (c - p); gj = q / (c - q)
    X = 2 + Ph(0); Y = 2 + Gi + Gj + Ph(p + q); Z = 2 + gi + Ph(p); W = 2 + gj + Ph(q)
    dl = (1 + gi) * (1 + gj); ep = (1 + Gi) * (1 + Gj)
    D = dl * X * Y - ep * Z * W
    Fii = 2 - X * (1 + gi) ** 2 / Z ** 2
    return q * Fii - (q + r) * D / (Z * W * Y)

res = {'adm_below_both_ext': 0, 'single_tau_below_ext': 0, 'two_atom_neg': 0, 'single_neg': 0}
worst_rel = 1e9; worst_two = (1e9, None); worst_single = (1e9, None)
for trial in range(4000):
    v = rng.exponential(1, 4) * np.exp(rng.normal(0, 2, 4))
    if rng.random() < 0.3: v[2] = 0
    v = v / v.sum(); p, q, r, e = v
    if min(p, q, e) <= 0: continue
    xT = p + q + r; c = 2 - xT
    Hone = Hval(p, q, r, e, [(e, xT)]); Htiny = Hval(p, q, r, e, [(e, c)])
    ext = min(Hone, Htiny)
    scale = q * (c - p)
    # admissible random
    for z in range(5):
        w = rng.exponential(1, rng.integers(2, 8)) ** 2; ms = e * w / w.sum()
        Ha = Hval(p, q, r, e, [(mm, c - 2 * mm) for mm in ms])
        if Ha < ext * (1 - 1e-9): res['adm_below_both_ext'] += 1
    taus = xT + (c - xT) * np.linspace(0, 1, 41)
    Hs = np.array([Hval(p, q, r, e, [(e, tau)]) for tau in taus])
    if Hs.min() < ext * (1 - 1e-9): res['single_tau_below_ext'] += 1
    if Hs.min() <= 0: res['single_neg'] += 1
    if Hs.min() / scale < worst_single[0]: worst_single = (Hs.min() / scale, (p, q, r, e), taus[np.argmin(Hs)])
    for z in range(20):
        t1, t2 = xT + (c - xT) * rng.random(2); th = rng.random()
        H2 = Hval(p, q, r, e, [(e * th, t1), (e * (1 - th), t2)])
        if H2 <= 0: res['two_atom_neg'] += 1
        if H2 / scale < worst_two[0]: worst_two = (H2 / scale, (p, q, r, e), t1, t2, th)
print(res); print('worst single/scale', worst_single); print('worst two-atom/scale', worst_two)
