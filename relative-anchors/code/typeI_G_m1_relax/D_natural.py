"""The local block D(s) of the reduced system in natural sums (exact, python-flint rational functions).

With c = S/2 + sigma and c0 = (sigma + S)/2, the moment symbols of D are rewritten through
  g = sum s/d = u - c M_1,  N = sum s/d^2 = M_1 - c M_3,  g0 = sum s/d0 = u - c0 M_0,  N0 = sum s/d0^2 = M_0 - c0 M_2.
The script checks symbolically that u cancels, writes D_natural.pkl (entries as rational functions of
(s, sigma, S, IG, IGT, g, N, g0, N0)), and prints the degree of each entry in each global quantity."""
import pickle, time
import flint
from flintfield import Field

t0 = time.time()
D = pickle.load(open('reduced_compact.pkl', 'rb'))
names = ['s', 'sigma', 'S', 'IG', 'IGT', 'g', 'N', 'g0', 'N0', 'u']
Kf = Field(names)
s, sg, S, IG, IGT, g, N, g0, N0, u = Kf.gens


def pw(v, e):
    r = Kf(1)
    for _ in range(e):
        r = r * v
    return r


def ev(terms, vals):
    tot = Kf(0)
    for m, p, q in terms:
        t = Kf(flint.fmpq(p, q))
        for v, e in zip(vals, m):
            if e:
                t = t * pw(v, e)
        tot = tot + t
    return tot


c = S / 2 + sg
c0 = (sg + S) / 2
M1 = (u - g) / c
M3 = (M1 - N) / c
M0 = (u - g0) / c0
M2 = (M0 - N0) / c0
val = {'IG': IG, 'IGT': IGT, 'u': u, 'M_0': M0, 'M_1': M1, 'M_2': M2, 'M_3': M3}
ent = {}
for a, b, ex, nt, dt in D['loc']:
    t = Kf(1)
    for cn, e in zip(D['consts'], ex):
        if e:
            t = t * pw(val[cn], e)
    t = t * ev(nt, (s, sg, S)) / ev(dt, (s, sg, S))
    ent[(a, b)] = ent.get((a, b), Kf(0)) + t
out = {}
for k, e in ent.items():
    num, den = e.n, e.d
    degs = {}
    for i, nm in enumerate(names):
        dn = max((m[i] for m in num.to_dict()), default=0)
        dd = max((m[i] for m in den.to_dict()), default=0)
        degs[nm] = (dn, dd)
    print(k, 'numerator terms', len(num.to_dict()), 'degrees (num, den):',
          {nm: degs[nm] for nm in ['u', 'IG', 'IGT', 'g', 'N', 'g0', 'N0', 's']}, flush=True)
    out[k] = (num.to_dict(), den.to_dict())
pickle.dump(dict(names=names, D=out), open('D_natural.pkl', 'wb'))
print(f'saved D_natural.pkl ({time.time() - t0:.1f} s)')
