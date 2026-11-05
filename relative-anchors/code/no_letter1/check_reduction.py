"""Exact checks that tie the certified objects to the block itself, at random rational configurations.
For the closed-form block M (djlib.block, itself checked against ylib in check_closed.py) and F = M with the vacuum
eliminated, for ordered pairs (i, l), with p = a_i, q = a_l, r = x_T - p - q, e = x_U:
 (1) F_ii = 2 - X (1 + g_i)^2/Z^2,  F_il = -D/(Z W Y),  a_l F_ii + (x_T - a_i) F_il = N_H/(Z^2 W Y);
 (2) D = (1/2) sum_{u,v} x_u x_v Dk(t_u, t_v), with Dk the very rational function certified in cert_pair.py,
     evaluated at s = (t - x_T)/(c - t);
 (3) Lambda N_H = sum over ordered 4-tuples of atoms x_u1 .. x_u4 P(s_u1, .., s_u4)/(24 e^4 prod den(t_ui)), with P the
     very polynomial certified in cert_row.py (s = t - x_T); P is symmetric in (s1, .., s4) (checked at points).
Also the U-empty formulas (e = 0)."""
import random, itertools, time
from math import factorial
from collections import Counter
import flint
from djlib import block, eliminate_vacuum
from cert_pair import build_kernel
from cert_row import build_P
Q = flint.fmpq
random.seed(5112026)
t0 = time.time()
_, Dk, _ = build_kernel()
ctxP, P = build_P()
print(f'objects built [{time.time() - t0:.1f}s]', flush=True)


def quantities(aT, xU, i, l):
    xT = sum(aT, Q(0)); c = 2 - xT; e = sum(xU, Q(0))
    p, q = aT[i], aT[l]; r = xT - p - q
    Ph = lambda d: sum((x / (c - 2 * x + d) for x in xU), Q(0))
    gi = p / (c - p); gl = q / (c - q); Gi = p / (c - p + q); Gl = q / (c - q + p)
    dl = (1 + gi) * (1 + gl); ep = (1 + Gi) * (1 + Gl)
    X = 2 + Ph(Q(0)); Y = 2 + Gi + Gl + Ph(p + q); Z = 2 + gi + Ph(p); W = 2 + gl + Ph(q)
    D = dl * X * Y - ep * Z * W
    NH = 2 * q * Z * Z * W * Y - q * (1 + gi) ** 2 * X * W * Y - (q + r) * dl * X * Y * Z + (q + r) * ep * Z * Z * W
    return dict(xT=xT, c=c, e=e, p=p, q=q, r=r, gi=gi, X=X, Y=Y, Z=Z, W=W, D=D, NH=NH)


ok = True
nconf = 0
configs = [(2, 0), (3, 0), (5, 0), (2, 1), (3, 1), (2, 2), (4, 2), (3, 3), (2, 3), (5, 3), (2, 4), (3, 4)]
for n, N in configs:
    w = [Q(random.randint(1, 40), random.randint(1, 40)) for _ in range(n + N)]
    if random.random() < 0.5:
        w[0] = w[0] * 50
    S = sum(w, Q(0)); w = [v / S for v in w]
    aT, xU = w[:n], w[n:]
    M = block(aT, xU); F = eliminate_vacuum(M)
    pairs = list(itertools.permutations(range(n), 2))
    random.shuffle(pairs)
    good = True
    for idx, (i, l) in enumerate(pairs):
        V = quantities(aT, xU, i, l)
        good &= F[i][i] == 2 - V['X'] * (1 + V['gi']) ** 2 / V['Z'] ** 2
        good &= F[i][l] == -V['D'] / (V['Z'] * V['W'] * V['Y'])
        good &= aT[l] * F[i][i] + (V['xT'] - aT[i]) * F[i][l] == V['NH'] / (V['Z'] ** 2 * V['W'] * V['Y'])
        if N == 0 or idx >= 2:
            continue
        p, q, r, e, c, xT = V['p'], V['q'], V['r'], V['e'], V['c'], V['xT']
        ts = [c - 2 * x for x in xU]
        # (2) pair kernel
        sc = [(t - xT) / (c - t) for t in ts]
        tot = Q(0)
        for u in range(N):
            for v in range(N):
                vals = [p, q, r, e, sc[u], sc[v]]
                tot += xU[u] * xU[v] * Dk.n(*vals) / Dk.d(*vals)
        good &= V['D'] == tot / 2
        # (3) row kernel
        al1 = c - p; al2 = c - q; be1 = c - p + q; be2 = c - q + p
        Lam = al1 ** 2 * al2 * be1 * be2
        den = lambda t: t * (t + p) * (t + q) * (t + p + q)
        tot = Q(0)
        for ms in itertools.combinations_with_replacement(range(N), 4):
            cnt = factorial(4)
            for m in Counter(ms).values():
                cnt //= factorial(m)
            vals = [p, q, r, e] + [ts[u] - xT for u in ms]
            pr = Q(1)
            for u in ms:
                pr *= xU[u] / den(ts[u])
            tot += cnt * pr * P(*vals)
        good &= Lam * V['NH'] == tot / (24 * e ** 4)
        # symmetry of P at one permuted point
        vals = [p, q, r, e] + [Q(random.randint(0, 9), random.randint(1, 9)) for _ in range(4)]
        perm = vals[:4] + [vals[5], vals[7], vals[4], vals[6]]
        good &= P(*vals) == P(*perm)
    nconf += 1
    ok &= good
    print(f'|T| = {n}, |U| = {N}: {len(pairs)} ordered pairs, formulas (1){" and kernel identities (2), (3)" if N else ""}: '
          f'{good} [{time.time() - t0:.1f}s]', flush=True)
print(nconf, 'configurations'); print('PASS' if ok else 'FAIL')
