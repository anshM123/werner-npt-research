"""Three levels (two single U-copies with weights s1, s2, and u3 >= 2 copies of weight s3): exact validation of
three_level_fast.py against direct elimination, and sign certificates with u3 = 2 + x3 for
  (Z1) G_z[rho,h] < 0 for every h;  (Z3) G_z[h,h'] < 0 for h, h' in different levels;  (D) G_z v > 0.
Not covered: G_z[h,h'] with h, h' both among the u3 equal copies (that needs the standard part of level 3).
sigma = s_t = 1 throughout (homogeneity: G_z(l s) = D_l G_z(s) D_l with D_l = diag(1, l, ..., l), and v scales
as l^2 D_l^{-1} v, so (Z) and (D) at sigma = 1 give them for every sigma)."""
import sys, os, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from gzexact import Gz_exact

Q = flint.fmpq
D = pickle.load(open('three_level.pkl', 'rb'))
names = D["names"]                       # s1, s2, s3, u3 (sigma = 1)
ctx = flint.fmpq_mpoly_ctx.get(tuple(names), 'lex')
Q4 = [[(ctx.from_dict(D['Q4'][i][j][0]), ctx.from_dict(D['Q4'][i][j][1])) for j in range(4)] for i in range(4)]
ok = True

t0 = time.time()
nv = 0
for a3 in (2, 3, 5):
    for wt, w1, w2, w3 in [(Q(1), Q(1, 2), Q(5, 3), Q(2, 7)), (Q(1), Q(1, 30), Q(1, 3), Q(4)),
                           (Q(1), Q(5, 2), Q(1, 100), Q(3, 5)), (Q(1), Q(7, 4), Q(2, 9), Q(1, 50))]:
        G = Gz_exact(wt, [w1, w2] + [w3] * a3)
        idx = [[0], [1], [2], list(range(3, 3 + a3))]
        direct = [[sum((G[i, j] for i in idx[p] for j in idx[q]), Q(0)) for q in range(4)] for p in range(4)]
        vals = [w1 * w1, w2 * w2, w3 * w3, Q(a3)]
        sym = [[Q4[p][q][0](*vals) / Q4[p][q][1](*vals) for q in range(4)] for p in range(4)]
        good = sym == direct
        ok &= good
        nv += 1
        if not good:
            print('   MISMATCH', a3, wt, w1, w2, w3)
print(f'validation: symbolic Q4 = exact elimination at {nv} points: {"all equal" if ok else "MISMATCH"} '
      f'({time.time() - t0:.1f} s)', flush=True)

TX = flint.fmpq_mpoly_ctx.get(("s1", "s2", "s3", "x3"), "lex")
g = TX.gens()
mp = {"s1": g[0], "s2": g[1], "s3": g[2], "u3": g[3] + 2}


def transport(p):
    out = TX.from_dict({})
    cache = {}
    for m, c in p.to_dict().items():
        term = TX.from_dict({(0,) * 4: c})
        for nm, e in zip(names, m):
            if e:
                if (nm, e) not in cache:
                    cache[(nm, e)] = mp[nm] ** e
                term = term * cache[(nm, e)]
        out = out + term
    return out


def sign(p):
    c = list(p.coeffs())
    return 1 if all(x > 0 for x in c) else (-1 if all(x < 0 for x in c) else 0)


T = [[(transport(Q4[i][j][0]), transport(Q4[i][j][1])) for j in range(4)] for i in range(4)]
s1, s2, s3, x3 = g
S = s1 + s2 + (x3 + 2) * s3
v = [4 * (1 + S)] + [3 * (2 + S)] * 3          # sigma = 1
res = []
for (i, j, nm) in [(0, 1, 'rho,z1'), (0, 2, 'rho,z2'), (0, 3, 'rho,z3'), (1, 2, 'z1,z2'), (1, 3, 'z1,z3'), (2, 3, 'z2,z3')]:
    n, d = T[i][j]
    sgn = sign(n) * sign(d)
    res.append(sgn == -1)
    print(f'Q4[{nm}]: numerator {len(n.coeffs())} terms, sign {sign(n):+d}; denominator sign {sign(d):+d}', flush=True)
for i, nm in enumerate(['rho', 'z1', 'z2', 'z3']):
    # row i of Q4 v, over a common denominator
    dens = [T[i][j][1] for j in range(4)]
    if all(dd == dens[0] for dd in dens):
        den = dens[0]
        num = sum((T[i][j][0] * v[j] for j in range(4)), TX.from_dict({}))
    else:
        num = TX.from_dict({})
        den = TX.from_dict({(0,) * 4: 1})
        for j in range(4):
            n, d = T[i][j]
            num = num * d + n * v[j] * den
            den = den * d
    sgn = sign(num) * sign(den)
    res.append(sgn == 1)
    print(f'(D) row {nm}: numerator {len(num.coeffs())} terms, sign {sign(num):+d}; denominator sign {sign(den):+d}',
          flush=True)
ok &= all(res)
print(f'certificates: {res} ({time.time() - t0:.1f} s)')
print('PASS' if ok else 'FAIL')
