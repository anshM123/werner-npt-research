"""Exploration (floating point; no claim): the per-copy lower bound Q_h of Psi on the perpendicular subspace p
(NOTES.md, Section 13), for given delta, theta:
  Q_h = sig (s - sig rho)(a^2 + c^2) + 4 sig IG s rho a c + sig s z^2
        + sig s^2 z (a + c)/d0 + sig IGT (S/dt0) s^2 z a/d0 - sig IGT g0 s^2 z a/d0 + (S sig/2) s z c/d0
        - (3/5) sig^2 (1 + delta) s^2 mu^2,  mu = 2(a + c)/d + z/d0
        - (9/5) sig^2 (1 + 1/delta) [ (N - g^2/S) 4 IG^2 s^3 (a^2 + c^2)/d^2 + (N0 - g0^2/S) IGT^2 s^3 z^2/d0^2 ]
        - (sig^2 IGT/4)(N0 - g0^2/S)(theta s c^2 + s z^2/theta).
Reports the smallest Jacobi-normalised eigenvalue of Q_h over random weight vectors and copies."""
import random, math, sys
import numpy as np


def Qh(sh, sg, S, g, N, g0, N0, delta, theta):
    d = sh + S / 2 + sg; d0 = sh + (sg + S) / 2; dt0 = (3 * sg + S) / 2; rho = sh / d
    IG = 1 / (2 + g); IGT = 1 / (2 + sg / dt0 + g0)
    M = np.zeros((3, 3))
    def add(i, j, v):
        if i == j:
            M[i, i] += v
        else:
            M[i, j] += v / 2; M[j, i] += v / 2
    add(0, 0, sg * (sh - sg * rho)); add(1, 1, sg * (sh - sg * rho)); add(0, 1, 4 * sg * IG * sh * rho)
    add(2, 2, sg * sh)
    add(0, 2, sg * sh * sh / d0 + sg * IGT * (S / dt0) * sh * sh / d0 - sg * IGT * g0 * sh * sh / d0)
    add(1, 2, sg * sh * sh / d0 + (S * sg / 2) * sh / d0)
    mu = np.array([2 / d, 2 / d, 1 / d0])
    M -= 0.6 * sg * sg * (1 + delta) * sh * sh * np.outer(mu, mu)
    VN = N - g * g / S; VN0 = N0 - g0 * g0 / S
    k = 1.8 * sg * sg * (1 + 1 / delta)
    add(0, 0, -k * VN * 4 * IG * IG * sh ** 3 / d ** 2); add(1, 1, -k * VN * 4 * IG * IG * sh ** 3 / d ** 2)
    add(2, 2, -k * VN0 * IGT * IGT * sh ** 3 / d0 ** 2)
    add(1, 1, -(sg * sg * IGT / 4) * VN0 * theta * sh); add(2, 2, -(sg * sg * IGT / 4) * VN0 * sh / theta)
    return M


rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
for delta in (0.05, 0.1, 0.2, 0.4):
    for theta in (0.5, 1.0, 2.0):
        worst = 9
        for it in range(3000):
            u = rng.choice([3, 4, 5, 6, 8, 12, 20])
            sg = math.exp(rng.uniform(-9, 9))
            fam = it % 4
            if fam == 0:
                sv = [math.exp(rng.uniform(-3, 3)) for _ in range(u)]
            elif fam == 1:
                sv = [math.exp(rng.uniform(-9, 9)) for _ in range(u)]
            elif fam == 2:
                sv = [1.0] * u; sv[0] = math.exp(rng.uniform(-10, 10))
            else:
                r = math.exp(rng.uniform(-3, 3)); sv = [r ** j for j in range(u)]
            S = sum(sv); c = S / 2 + sg; c0 = (sg + S) / 2
            g = sum(x / (x + c) for x in sv); N = sum(x / (x + c) ** 2 for x in sv)
            g0 = sum(x / (x + c0) for x in sv); N0 = sum(x / (x + c0) ** 2 for x in sv)
            for sh in sv:
                M = Qh(sh, sg, S, g, N, g0, N0, delta, theta)
                dg = np.diag(M)
                if np.any(dg <= 0):
                    worst = -1; break
                sc = 1 / np.sqrt(dg)
                worst = min(worst, np.linalg.eigvalsh(M * np.outer(sc, sc)).min())
            if worst < 0:
                break
        print(f'delta = {delta}, theta = {theta}: smallest Jacobi-normalised eigenvalue of Q_h: {worst:.4f}', flush=True)
