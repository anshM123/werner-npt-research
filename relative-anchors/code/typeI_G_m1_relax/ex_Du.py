"""Exact test: does the local block D(s) depend on u once g = sum s/d, N = sum s/d^2, g0 = sum s/d0, N0 = sum s/d0^2
(and IG, IGT) are fixed?  Evaluate D at random rational values of (s, sigma, S, g, N, g0, N0) for two values of u,
with M_1 = (u - g)/c, M_3 = (M_1 - N)/c, M_0 = (u - g0)/c0, M_2 = (M_0 - N0)/c0, c = S/2 + sigma, c0 = (sigma + S)/2."""
import random, pickle
import flint
Q = flint.fmpq
D = pickle.load(open('reduced_compact.pkl', 'rb'))
consts = D['consts']
def ev(terms, vals):
    tot = Q(0)
    for m, p, q in terms:
        t = Q(p, q)
        for v, e in zip(vals, m):
            if e:
                t *= v ** e
        tot += t
    return tot
def Dmat(sv, sg, SS, IGv, IGTv, u, g, N, g0, N0):
    c = SS / 2 + sg; c0 = (sg + SS) / 2
    M1 = (u - g) / c; M3 = (M1 - N) / c; M0 = (u - g0) / c0; M2 = (M0 - N0) / c0
    val = {'IG': IGv, 'IGT': IGTv, 'u': u, 'M_0': M0, 'M_1': M1, 'M_2': M2, 'M_3': M3}
    out = {}
    for a, b, ex, nt, dt in D['loc']:
        t = Q(1)
        for cn, e in zip(consts, ex):
            if e:
                t *= val[cn] ** e
        t *= ev(nt, (sv, sg, SS)) / ev(dt, (sv, sg, SS))
        out[(a, b)] = out.get((a, b), Q(0)) + t
    return out
random.seed(5)
for trial in range(4):
    r = lambda: Q(random.randint(1, 50), random.randint(1, 50))
    sv, sg, SS, IGv, IGTv, g, N, g0, N0 = r(), r(), r(), r(), r(), r(), r(), r(), r()
    d1 = Dmat(sv, sg, SS, IGv, IGTv, Q(7), g, N, g0, N0)
    d2 = Dmat(sv, sg, SS, IGv, IGTv, Q(19), g, N, g0, N0)
    print('trial', trial, 'D independent of u at fixed natural sums:', all(d1[k] == d2[k] for k in d1),
          {k: float(d1[k] - d2[k]) for k in d1})
