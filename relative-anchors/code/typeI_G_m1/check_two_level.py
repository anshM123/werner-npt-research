"""Two-level U-weights (u1 copies of weight s1, u2 copies of weight s2): exact validation of the symbolic reduction
(two_level.py) against the exact engine, and sign certificates.

Q3 is the matrix of G_z on the invariant letter-2 data (rho, z1, z2):
  Q3[0,0] = G_z[rho,rho], Q3[0,i] = sum_{h in U_i} G_z[rho,h], Q3[i,j] = sum_{h in U_i, h' in U_j} G_z[h,h'].
Certificates (u_i = u_i0 + x_i, x_i >= 0; every polynomial below has coefficients of one sign):
  (Z1) Q3[0,1] < 0 and Q3[0,2] < 0: G_z[rho,h] < 0 for every h;
  (Z3) Q3[1,2] < 0: G_z[h,h'] < 0 for h, h' in different levels;
  (D)  Q3 (v0, v1, v1) > 0 entrywise: (G_z v)_rho > 0 and (G_z v)_h > 0 for every h (u_i (G_z v)_h = row i).
The off-diagonal entries G_z[h,h'] with h, h' in the same level are not covered here.
Usage: python check_two_level.py case   (case = one: u1 = 1, u2 >= 2;  gen: u1, u2 >= 2)"""
import sys, os, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from flintfield import Field, F as FF
from gzexact import Gz_exact

Q = flint.fmpq
case = sys.argv[1]
D = pickle.load(open(f'two_level_{case}.pkl', 'rb'))
names = D['names']            # [sigma, s1, s2, u2] (+ u1 for gen)
K = Field(names)
Q3 = [[FF(K, K.ctx.from_dict(D['Q3'][i][j][0]), K.ctx.from_dict(D['Q3'][i][j][1])) for j in range(3)] for i in range(3)]
ok = True


def ev(f, vals):
    return f.n(*vals) / f.d(*vals)


# ---- validation
t0 = time.time()
pts = [(Q(3, 4), Q(1, 2), Q(5, 3)), (Q(2), Q(1, 30), Q(1, 3)), (Q(1, 7), Q(5, 2), Q(1, 100))]
uus = [(1, 2), (1, 3), (1, 5), (1, 7)] if case == 'one' else [(2, 2), (2, 3), (3, 2), (3, 4), (4, 3)]
nv = 0
for (a1, a2) in uus:
    for wt, om1, om2 in pts:
        wU = [om1] * a1 + [om2] * a2
        G = Gz_exact(wt, wU)
        I1 = list(range(1, a1 + 1))
        I2 = list(range(a1 + 1, a1 + a2 + 1))
        idx = [[0], I1, I2]
        direct = [[sum((G[i, j] for i in idx[p] for j in idx[q]), Q(0)) for q in range(3)] for p in range(3)]
        vals = [wt * wt, om1 * om1, om2 * om2, Q(a2)] + ([] if case == 'one' else [Q(a1)])
        sym = [[ev(Q3[p][q], vals) for q in range(3)] for p in range(3)]
        good = sym == direct
        ok &= good
        nv += 1
        if not good:
            print('   MISMATCH', a1, a2, wt, om1, om2)
print(f'validation: symbolic Q3 = exact elimination at {nv} points: {"all equal" if ok else "MISMATCH"} '
      f'({time.time() - t0:.1f} s)', flush=True)

# ---- certificates
sig, s1, s2 = K.gens[0], K.gens[1], K.gens[2]
U2 = K.gens[3]
U1 = K.gens[4] if case == 'gen' else K(1)
S = U1 * s1 + U2 * s2
v0 = 4 * sig * (sig + S)
v1 = 3 * (2 * sig + S)
Kx = Field(['sigma', 's1', 's2', 'x2'] + ([] if case == 'one' else ['x1']))
shift = {3: 2, 4: 2}            # u2 = 2 + x2, u1 = 2 + x1


def subst(p):
    """polynomial in K.ctx -> polynomial in Kx.ctx with u_i = 2 + x_i."""
    gx = Kx.ctx.gens()
    out = Kx.ctx.from_dict({})
    pw = {}
    for m, c in p.to_dict().items():
        term = Kx.ctx.from_dict({tuple(m[:3]) + (0,) * (len(m) - 3): c})
        for vi in (3, 4):
            if vi < len(m) and m[vi]:
                key = (vi, m[vi])
                if key not in pw:
                    pw[key] = (gx[vi] + shift[vi]) ** m[vi]
                term = term * pw[key]
        out = out + term
    return out


def cert(f, name):
    num, den = subst(f.n), subst(f.d)
    cn = [c for c in num.coeffs()]
    cd = [c for c in den.coeffs()]
    sn = 1 if all(c > 0 for c in cn) else (-1 if all(c < 0 for c in cn) else 0)
    sd = 1 if all(c > 0 for c in cd) else (-1 if all(c < 0 for c in cd) else 0)
    print(f'{name}: numerator {len(cn)} terms sign {sn:+d}, denominator {len(cd)} terms sign {sd:+d}', flush=True)
    return sn * sd


t0 = time.time()
r = []
r.append(cert(Q3[0][1], '(Z1) Q3[rho,z1]') == -1)
r.append(cert(Q3[0][2], '(Z1) Q3[rho,z2]') == -1)
r.append(cert(Q3[1][2], '(Z3) Q3[z1,z2]') == -1)
for i, nm in enumerate(['rho', 'z1', 'z2']):
    e = Q3[i][0] * v0 + (Q3[i][1] + Q3[i][2]) * v1
    r.append(cert(e, f'(D) row {nm}') == 1)
ok &= all(r)
print(f'certificates: {r} ({time.time() - t0:.1f} s)')
print('PASS' if ok else 'FAIL')
