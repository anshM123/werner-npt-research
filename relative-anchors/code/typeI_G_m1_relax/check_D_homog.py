"""Exact check of the homogeneity used in the normalisation K = 1 of prove_D.py: every entry of D_natural satisfies
D(l s, l sigma, l S, g, N/l, g0, N0/l) = l^2 D(s, sigma, S, g, N, g0, N0) (IG = 1/(2+g), IGT = 1/(2 + sigma/d0_t + g0)
are invariant).  Tested at random rational points and scalings; it also follows from the degrees of the monomials."""
import random, pickle
import flint
Q = flint.fmpq
P = pickle.load(open('D_natural.pkl', 'rb'))
names = P['names']
def evald(d, vals):
    tot = Q(0)
    for m, c in d.items():
        t = Q(c)
        for v, e in zip(vals, m):
            if e:
                t *= v ** int(e)
        tot += t
    return tot
def D(sv, sg, SS, g, N, g0, N0):
    dt0 = (3 * sg + SS) / 2
    vals = dict(s=sv, sigma=sg, S=SS, IG=1 / (2 + g), IGT=1 / (2 + sg / dt0 + g0), g=g, N=N, g0=g0, N0=N0, u=Q(1))
    vv = [vals[n] for n in names]
    return {k: evald(nd, vv) / evald(dd, vv) for k, (nd, dd) in P['D'].items()}
random.seed(3)
ok = True
for it in range(10):
    r = lambda: Q(random.randint(1, 40), random.randint(1, 40))
    sv, sg, SS, g, N, g0, N0, lam = r(), r(), r(), r(), r(), r(), r(), r()
    A = D(sv, sg, SS, g, N, g0, N0)
    B = D(lam * sv, lam * sg, lam * SS, g, N / lam, g0, N0 / lam)
    ok &= all(B[k] == lam ** 2 * A[k] for k in A)
print('homogeneity of degree 2 (N, N0 of degree -1):', ok)
print('PASS' if ok else 'FAIL')
