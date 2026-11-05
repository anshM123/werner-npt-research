"""Leading-principal-minor certificate by fraction-free (Bareiss) elimination with the denominator factors pulled out.

Input: blocks_sym.pkl (build_sym.py).  Let Mk be the gauge-fixed block (entries in Q(s1, s2, s3), s0 = 1), reordered
by the chosen order, den the lcm of the denominators of its entries (times a positive integer so that Nt = den * Mk has
integer coefficients), den = cst * prod_i f_i^{E_i} with f_i irreducible.  Bareiss elimination of Nt: after pivot p
the entry (i, l), i, l > p, equals the (p + 2)-minor of Nt on rows {0..p, i} and columns {0..p, l} (Sylvester's
identity), and the division by the previous pivot is exact.  Every entry is stored as prod_i f_i^{x_i} * core, the
factors f_i being divided out exactly after every step.  The leading minors are
        D_r(Mk) = D_r(Nt) / den^r = core_r * prod_i f_i^{x_ir - r E_i} / cst^r.
If cst > 0, every f_i and every core_r has only positive coefficients, every leading principal minor of Mk is positive on
the open orthant and Mk > 0 there (Sylvester's criterion).  Exact checks: every D_r(Nt) is compared with a direct
determinant of the leading r x r part of Nt at a random integer point.  The leading cores are saved to
cert_j{j}_T{T}_{order}.pkl.  Usage: python bareiss_cert.py j T order   (T: '-' or e.g. 0, 01, 012; order as in
elim_minors.py, or 'file:NAME' to read a comma-separated order from a file)."""
import sys, time, pickle, random
from math import lcm
from fractions import Fraction as Fr
import flint

j = int(sys.argv[1])
Tstr = sys.argv[2]
T = () if Tstr == '-' else tuple(int(c) for c in Tstr)
ordspec = sys.argv[3] if len(sys.argv) > 3 else 'nat'
D = pickle.load(open('blocks_sym.pkl', 'rb'))
names = tuple(D['names'])
blk = D['blocks'][(j, T)]
S = blk['states']
keep = blk['keep']
n = len(keep)
labs = [blk['labels'][a] for a in keep]
vac = (0, 0, 0, 0)
isvac = [S[a][0] == vac for a in keep]
if ordspec.startswith('file:'):
    ordspec_s = open(ordspec[5:]).read().strip()
else:
    ordspec_s = ordspec
if ordspec_s == 'nat':
    order = list(range(n))
elif ordspec_s == 'rev':
    order = list(range(n))[::-1]
elif ordspec_s == 'vacfirst':
    order = [i for i in range(n) if isvac[i]] + [i for i in range(n) if not isvac[i]]
elif ordspec_s == 'vaclast':
    order = [i for i in range(n) if not isvac[i]] + [i for i in range(n) if isvac[i]]
elif ordspec_s == 'col':
    order = sorted(range(n), key=lambda i: (isvac[i], S[keep[i]][1], S[keep[i]][0]))
else:
    order = [int(x) for x in ordspec_s.split(',')]
assert sorted(order) == list(range(n))
ctxQ = flint.fmpq_mpoly_ctx.get(names, 'lex')
ctxZ = flint.fmpz_mpoly_ctx.get(names, 'lex')
t0 = time.time()
E = [[blk['M'][keep[order[a]]][keep[order[b]]] for b in range(n)] for a in range(n)]
L = 1
for row in E:
    for (nn, dd) in row:
        for c in list(nn.values()) + list(dd.values()):
            L = lcm(L, int(c.q))


def toZ(dct):
    return ctxZ.from_dict({m: int(c * L) for m, c in dct.items()})


NUMS = [[toZ(E[a][b][0]) for b in range(n)] for a in range(n)]
DENS = [[toZ(E[a][b][1]) for b in range(n)] for a in range(n)]
den = DENS[0][0]
for a in range(n):
    for b in range(n):
        g = den.gcd(DENS[a][b])
        den = den * divmod(DENS[a][b], g)[0]
Nt = [[NUMS[a][b] * divmod(den, DENS[a][b])[0] for b in range(n)] for a in range(n)]
for a in range(n):
    for b in range(n):
        assert divmod(den, DENS[a][b])[1] == 0
if den.leading_coefficient() < 0:
    den = -den
    Nt = [[-x for x in row] for row in Nt]
assert all(Nt[a][b] == Nt[b][a] for a in range(n) for b in range(n))
cst, fac = den.factor()
FACS = [f for f, e in fac]
EXPS = [e for f, e in fac]
nf = len(FACS)
print(f'j = {j}, T = {T}, size {n}, order {ordspec}: {[labs[i] for i in order]}', flush=True)
print(f'  den = {cst} * ' + ' * '.join(f'[deg {f.total_degree()}, {len(f)} t]^{e}' for f, e in fac), flush=True)
denpos = cst > 0 and all(all(c > 0 for c in f.coeffs()) for f in FACS)
print('  den constant > 0 and den factors with positive coefficients:', denpos, flush=True)


def exact_div(p, q):
    """p / q if q divides p exactly (flint's exact division, which stops early when it is not exact), else None."""
    try:
        return p / q
    except Exception:
        return None


def strip(p):
    ex = [0] * nf
    if p == 0:
        return p, ex
    for i, f in enumerate(FACS):
        while True:
            q = exact_div(p, f)
            if q is None:
                break
            p = q
            ex[i] += 1
    return p, ex


ZERO = ctxZ.from_dict({})
b = [[None] * n for _ in range(n)]
for a in range(n):
    for c in range(a, n):
        b[a][c] = strip(Nt[a][c])
leading = [b[0][0]]
prev = (ctxZ.from_dict({(0,) * len(names): 1}), [0] * nf)
for p in range(n - 1):
    cp, ep = b[p][p]
    if cp == 0:
        print('  zero pivot', p, flush=True)
        break
    cprev, eprev = prev
    for i in range(p + 1, n):
        ci, ei = b[p][i] if p <= i else b[i][p]
        for l in range(i, n):
            cl, el = b[p][l]
            cil, eil = b[i][l]
            # t1 = b[p][p] * b[i][l],  t2 = b[i][p] * b[p][l]  (b[i][p] = b[p][i])
            e1 = [x + y for x, y in zip(ep, eil)]
            e2 = [x + y for x, y in zip(ei, el)]
            if cil == 0:
                g = e2
            elif ci == 0 or cl == 0:
                g = e1
            else:
                g = [min(x, y) for x, y in zip(e1, e2)]
            t = ZERO
            if cil != 0:
                m1 = cp * cil
                for f, x, gg in zip(FACS, e1, g):
                    if x > gg:
                        m1 = m1 * f ** (x - gg)
                t = t + m1
            if ci != 0 and cl != 0:
                m2 = ci * cl
                for f, x, gg in zip(FACS, e2, g):
                    if x > gg:
                        m2 = m2 * f ** (x - gg)
                t = t - m2
            if t == 0:
                b[i][l] = (ZERO, [0] * nf)
                continue
            q = exact_div(t, cprev)
            assert q is not None, 'inexact division by previous core'
            ex = [gg - x for gg, x in zip(g, eprev)]
            for k_, f in enumerate(FACS):
                while ex[k_] < 0:
                    q2 = exact_div(q, f)
                    assert q2 is not None, 'inexact division by a factor'
                    q = q2
                    ex[k_] += 1
            c2, e3 = strip(q)
            b[i][l] = (c2, [x + y for x, y in zip(ex, e3)])
    prev = b[p][p]
    leading.append(b[p + 1][p + 1])
    lc = b[p + 1][p + 1][0]
    print(f'  step {p + 1}: D_{p + 2} core deg {lc.total_degree()}, {len(lc)} terms, exps {b[p + 1][p + 1][1]}, '
          f'negative coeffs {sum(1 for c in lc.coeffs() if c < 0)}; max core terms in the rest '
          f'{max(len(b[i][l][0]) for i in range(p + 1, n) for l in range(i, n))}  [{time.time() - t0:.1f}s]', flush=True)
    # free the row/col p
    for l in range(p, n):
        b[p][l] = None if l != p else b[p][l]
# exact check at a random integer point
random.seed(11)
pt = [random.randint(2, 30) for _ in names]
Nv = [[Fr(int(Nt[a][c](*pt))) for c in range(n)] for a in range(n)]


def detF(A):
    A = [row[:] for row in A]
    k = len(A)
    det_ = Fr(1)
    for p in range(k):
        piv = next((i for i in range(p, k) if A[i][p] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != p:
            A[p], A[piv] = A[piv], A[p]
            det_ = -det_
        det_ *= A[p][p]
        for i in range(p + 1, k):
            f = A[i][p] / A[p][p]
            for l in range(p, k):
                A[i][l] -= f * A[p][l]
    return det_


ok = denpos
allpos = True
for r, (c, ex) in enumerate(leading, 1):
    val = Fr(int(c(*pt)))
    for f, e in zip(FACS, ex):
        val *= Fr(int(f(*pt))) ** e
    good = val == detF([row[:r] for row in Nv[:r]])
    ok &= good
    coeffs = c.coeffs()
    pos = all(x > 0 for x in coeffs)
    allpos &= pos
    print(f'  leading minor {r}: point check {good}; core deg {c.total_degree()}, {len(c)} terms, exps {ex}, '
          f'{"all coefficients positive" if pos else "MIXED (" + str(sum(1 for x in coeffs if x < 0)) + " negative)"}',
          flush=True)
ok &= allpos and len(leading) == n
tag = f'cert_j{j}_T{Tstr}_{ordspec.replace(":", "_").replace(",", "")[:40]}'
pickle.dump(dict(names=names, order=order, labels=[labs[i] for i in order], den_const=int(cst),
                 den_factors=[(f.to_dict(), e) for f, e in fac],
                 cores=[(c.to_dict(), ex) for c, ex in leading]), open(tag + '.pkl', 'wb'))
if ok:
    print(f'CONCLUSION: all {n} leading principal minors of the gauge-fixed block are positive on the open orthant; '
          f'the block j = {j}, T = {T} is >= 0 with kernel exactly the known one, for all positive weights')
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
