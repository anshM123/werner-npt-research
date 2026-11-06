"""The fixed-u1 forms tl_inv_fN_T11/T22 (u1 = N, u2 symbolic) against the symbolic forms tl_inv_gen_T11/T22 (u1, u2
symbolic; validated by check_tl.py against direct elimination): every entry of the blocks (A, A) and (A, B) agrees
exactly at 4 random rational points (u2, sigma, s1) for every N = 2, ..., NMAX, after putting u1 = N in the symbolic
form.  (An exact check, pointwise; the certificates use the fixed-u1 forms.)
Part 2 (exact evaluation, not a proof): the three LDL^T pivots of (B_inv) with H = 1/2, computed from the symbolic forms
tl_inv_gen_T11/T22 at rational points with u1, u2 >= 17 (outside the range of Theorem TV), s_t = s and s_t = 1.
python check_fN.py NMAX"""
import sys, os, random, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint

Q = flint.fmpq
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 16


def load(fn):
    D = pickle.load(open(fn, 'rb'))
    ctx = flint.fmpq_mpoly_ctx.get(tuple(D['names']), 'lex')
    P = lambda d: ctx.from_dict({m: Q(p, q) for m, (p, q) in d.items()})
    return D['names'], {k: (P(n), P(d)) for k, (n, d) in D['G'].items()}


rng = random.Random(61106)
ok = True
cnt = 0
gen = {md: load('tl_inv_gen_%s.pkl' % md) for md in ('T11', 'T22')}
for N in range(2, NMAX + 1):
    for md in ('T11', 'T22'):
        gnames, G = gen[md]
        fnames, Fm = load('tl_inv_f%d_%s.pkl' % (N, md))
        good = set(G) == set(Fm)
        for _ in range(4):
            vals = {'u1': Q(N), 'u2': Q(rng.randint(2, 9)), 'sig': Q(rng.randint(1, 60), rng.randint(1, 13)),
                    's1': Q(rng.randint(1, 40), rng.randint(1, 17))}
            for key in G:
                gn, gd = G[key]
                fn_, fd = Fm[key]
                ga = [vals[x] for x in gnames]
                fa = [vals[x] for x in fnames]
                good &= gn(*ga) / gd(*ga) == fn_(*fa) / fd(*fa)
                cnt += 1
        ok &= good
        print('N = %d, %s: %s' % (N, md, 'equal' if good else 'MISMATCH'), flush=True)
print('%d entry values compared' % cnt)
# ---- part 2: (B_inv) pivots at large sizes, from the symbolic forms
H = Q(1, 2)
npt = 0
minp = None
for (u1, u2) in ((17, 17), (17, 40), (25, 19), (50, 50), (100, 30), (200, 200)):
    for md, st_is_s in (('T11', True), ('T22', False)):
        gnames, G = gen[md]
        for s1 in (Q(1, 1000), Q(1, 10), Q(1), Q(7, 3), Q(10), Q(1000)):
            for rr in (Q(1, 1000), Q(1, 3), Q(2), Q(50)):
                st = s1 if st_is_s else Q(1)
                sig = st + rr
                vals = {'u1': Q(u1), 'u2': Q(u2), 'sig': sig, 's1': s1}
                args = [vals[x] for x in gnames]
                Gz = [[G[(min(i, j), max(i, j))][0](*args) / G[(min(i, j), max(i, j))][1](*args) for j in range(3)]
                      for i in range(3)]
                S1, S2 = u1 * s1, Q(u2)
                S = S1 + S2
                Nm = [[Q(1), Q(0), Q(0)], [Q(0), 1 / S, S2 / S], [Q(0), 1 / S, -S1 / S]]
                Gy = [[sum((Nm[k][i] * Gz[k][l] * Nm[l][j] for k in range(3) for l in range(3)), Q(0)) for j in range(3)]
                      for i in range(3)]
                kap = [Q(1), st * sig / (S * (sig + S)), st * S1 * S2 / S]
                M = [[Gy[i][j] - (H * (rr / sig) * kap[i] if i == j else 0) for j in range(3)] for i in range(3)]
                A = [row[:] for row in M]
                piv = []
                for p in range(3):
                    piv.append(A[p][p])
                    for i in range(p + 1, 3):
                        f = A[i][p] / A[p][p]
                        for j in range(p + 1, 3):
                            A[i][j] -= f * A[p][j]
                good = all(p > 0 for p in piv)
                ok &= good
                npt += 1
                rel = piv[2] / (H * rr / sig * kap[2])
                minp = rel if minp is None or rel < minp else minp
print('part 2: (B_inv) pivots positive at %d rational points with u1, u2 >= 17: %s; smallest third pivot relative to '
      'H (1 - theta) kappa_2^2: %.4g' % (npt, ok, float(minp)))
print('PASS' if ok else 'FAIL')
