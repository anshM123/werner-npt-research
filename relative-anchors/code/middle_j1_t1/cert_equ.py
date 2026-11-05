"""Certificate for equal U-weights (every k >= 3, every sigma): the invariant form Finv (3 x 3) and the standard form
Fstd (2 x 2) of equ_sym.py are positive definite for all p, q > 0 and every real u >= 2 (u = 2 + x, x >= 0).

Sylvester: every leading principal minor is a rational function of (p, q, x); it is factored over Q and every
irreducible factor of numerator and denominator is checked to have coefficients of one sign; the overall sign is then
read off (positive).  Negative control: the same check on Finv - c E00 and Fstd - c E00 with c = ..., and on the
minors of the forms with the B-part scaled is in check_equ.py (exact points)."""
import pickle, time
import flint
from flintfield import Field, F

t0 = time.time()
D = pickle.load(open('logs/equ_sym.pkl', 'rb'))
K = Field(D['names'])
ctx = K.ctx


def rf(nd):
    n, d = nd
    return F(K, ctx.from_dict(n), ctx.from_dict(d))


Finv = [[rf(v) for v in r] for r in D['Finv']]
Fstd = [[rf(v) for v in r] for r in D['Fstd']]


def det(A):
    n = len(A)
    if n == 1:
        return A[0][0]
    if n == 2:
        return A[0][0] * A[1][1] - A[0][1] * A[1][0]
    return (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
            + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))


def signcheck(f):
    """returns +1 / -1 if f = c * prod g_i^e_i with every g_i one-signed in coefficients (sign of f on the open
    orthant then known), else None; also a short description."""
    out = 1
    desc = []
    for poly in (f.n, f.d):
        c, fac = poly.factor()
        sg = 1 if c > 0 else -1
        for g, e in fac:
            cs = list(g.coeffs())
            if all(v > 0 for v in cs):
                pass
            elif all(v < 0 for v in cs):
                sg *= (-1) ** e
            else:
                return None, f'factor with mixed signs: deg {g.total_degree()}, {len(g)} terms'
            desc.append(f'deg {g.total_degree()}/{len(g)}t^{e}')
        out *= sg
    return out, ', '.join(desc)


ok = True
for name, A in (('Finv', Finv), ('Fstd', Fstd)):
    for m in range(1, len(A) + 1):
        dm = det([row[:m] for row in A[:m]])
        sg, desc = signcheck(dm)
        good = sg == 1
        ok &= good
        print(f'{name} leading minor {m}: sign {sg}; factors: {desc}', flush=True)
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
