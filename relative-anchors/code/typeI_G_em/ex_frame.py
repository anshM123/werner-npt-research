"""Exploration (no claim rests on it): the weighted block Gershgorin quantities at equal U-weights, from the closed
forms of equal_m.py, evaluated exactly at rational points and then rounded.

Frame (Proposition U1, Theorem U2): invariant part, coordinates (rho_t, zhat_t) with zhat_t = kappa_t z_t,
kappa_t^2 = s_t S sigma/(sigma + S); standard part, zhat = w_t om |f| z.  Normalisation sigma + S = 1, S = mu,
sigma = 1 - mu, s = mu/u, s_t = theta sigma, s_t' = phi sigma.
  A_inv = ||Ghat_tt'||_F sigma/sqrt(s_t s_t'),  B_inv = lambda_min(Ghat_tt)/(1 - theta),
  A_std = |Ghat_tt'| sigma/sqrt(s_t s_t'),     B_std = Ghat_tt/(1 - theta)."""
import sys, os, math, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from emlib import load, ev, Q

R = load()


def blocks(u, mu, th, ph):
    sig = 1 - mu
    S = mu
    s = mu / u
    st, stp = th * sig, ph * sig
    vd = [u, sig, s, st, 0]
    vo = [u, sig, s, st, stp]
    Gd = np.array([[float(ev(R['inv'][i][j], vd)) for j in range(2)] for i in range(2)])
    Go = np.array([[float(ev(R['inv'][i][2 + j], vo)) for j in range(2)] for i in range(2)])
    Gd2 = np.array([[float(ev(R['inv'][2 + i][2 + j], vo)) for j in range(2)] for i in range(2)])
    sd = float(ev(R['std'][0][0], vd))
    so = float(ev(R['std'][0][1], vo))
    kap = lambda x: math.sqrt(float(x * S * sig / (sig + S)))
    D1 = np.diag([1, 1 / kap(st)])
    D2 = np.diag([1, 1 / kap(stp)]) if stp != 0 else np.eye(2)
    Gh = D1 @ Gd @ D1
    Gho = D1 @ Go @ D2
    sF = float(s)
    stdd = sd / (float(st) * sF)
    stdo = so / (math.sqrt(float(st * stp)) * sF) if stp != 0 else 0.0
    return Gh, Gho, stdd, stdo


fr = lambda x: Q(x).limit_denominator(10 ** 6) if not isinstance(x, int) else Q(x)
worst = dict(Ainv=0, Binv=9, Astd=0, Bstd=9, B1inv=9, B1std=9)
loc = {}
mus = [Q(1, 10 ** 4), Q(1, 100), Q(1, 10), Q(3, 10), Q(1, 2), Q(7, 10), Q(9, 10), Q(99, 100), Q(9999, 10 ** 4)]
ths = [Q(1, 1000), Q(1, 20), Q(1, 4), Q(1, 2), Q(3, 4), Q(19, 20)]
for u in (3, 4, 6, 10, 30, 100, 1000):
    for mu in mus:
        for th in ths:
            for frac in (Q(1, 100), Q(1, 2), Q(1)):
                ph = (1 - th) * frac
                Gh, Gho, sd, so = blocks(u, mu, th, ph)
                A = np.linalg.norm(Gho) / math.sqrt(float(th * ph))
                B = np.linalg.eigvalsh(Gh).min() / float(1 - th)
                As = abs(so) / math.sqrt(float(th * ph))
                Bs = sd / float(1 - th)
                for key, val, better in (('Ainv', A, max), ('Binv', B, min), ('Astd', As, max), ('Bstd', Bs, min)):
                    nv = better(worst[key], val)
                    if nv != worst[key]:
                        worst[key] = nv
                        loc[key] = (u, str(mu), str(th), str(ph))
        # m = 1: theta = 1
        Gh, _, sd, _ = blocks(u, mu, Q(1), Q(0))
        lm = np.linalg.eigvalsh(Gh).min()
        if lm < worst['B1inv']:
            worst['B1inv'] = lm
            loc['B1inv'] = (u, str(mu))
        if sd < worst['B1std']:
            worst['B1std'] = sd
            loc['B1std'] = (u, str(mu))
    print(u, {k: round(v, 5) for k, v in worst.items()}, flush=True)
print(worst)
print(loc)
print('PASS' if worst['Ainv'] < 0.5 < worst['Binv'] and worst['Astd'] < 0.5 < worst['Bstd'] else 'FAIL')
