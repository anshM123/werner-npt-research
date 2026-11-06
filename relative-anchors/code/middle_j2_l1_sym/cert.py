"""Certificate for one family and one size regime (NOTES Section 4).

python cert.py FAMILY U0 U1 [save=FILE]
  FAMILY: a name of r2v.families() (inv, mix, std0, two0, anti0, three0, tto0, twostd0, antistd0, and the same with 1);
  U0, U1: a fixed level size (an integer m >= 2) or 'o+' for every real size u = o + x, x >= 0 (o >= 4).
Squared weights 1 (level 0) and q (level 1), q > 0 free.
Steps: bases of e'_lambda K, e'_lambda W and e'_lambda M_Q chosen at a sample size of the regime; the reduced forms over
Q(q, x...) (python-flint); the Gram determinants det G_P, det G_Q are positive on the regime (sign test), so the
vectors are a basis on the whole regime; S' (the congruent Schur complement) vanishes identically on the K-rows;
then every pivot of the symmetric elimination of the W-block (= ratio of consecutive leading principal minors) is a
ratio of polynomials that are positive for q > 0, x >= 0 (sign test: every irreducible factor one-signed, with a
monomial free of the x-variables).  PASS iff all of this holds."""
import sys, time, pickle
from fractions import Fraction as Fr
from flintfield import Field
import r2v

t0 = time.time()


def signcheck(f):
    """(sign, description) if f = c prod g_i^e_i with every irreducible g_i one-signed in coefficients and with a
    monomial free of the x-variables (so g_i != 0 for q > 0, x >= 0); sign is then the sign of f there.
    (None, reason) otherwise.  (As in middle_j1_l1/cert.py.)"""
    out = 1
    desc = []
    for poly in (f.n, f.d):
        names = poly.context().names()
        xi = [i for i, nm in enumerate(names) if nm.startswith('x')]
        cs = list(poly.coeffs())
        if (all(v > 0 for v in cs) or all(v < 0 for v in cs)) and any(all(m[i] == 0 for i in xi) for m in poly.monoms()):
            out *= 1 if cs[0] > 0 else -1
            desc.append(f'[one-signed {poly.total_degree()}/{len(poly)}]')
            continue
        c, fac = poly.factor()
        sg = 1 if c > 0 else -1
        for gg, e in fac:
            cs = list(gg.coeffs())
            if all(v > 0 for v in cs):
                pass
            elif all(v < 0 for v in cs):
                sg *= (-1) ** e
            else:
                return None, f'factor with mixed signs: deg {gg.total_degree()}, {len(gg)} terms'
            if not any(all(m[i] == 0 for i in xi) for m in gg.monoms()):
                return None, f'factor {gg} vanishes at x = 0'
            desc.append(f'{gg.total_degree()}/{len(gg)}^{e}')
        out *= sg
    return out, ' '.join(desc)


def det(M, one):
    n = len(M)
    A = [list(r) for r in M]
    d = one
    for c in range(n):
        piv = next((i for i in range(c, n) if not (A[i][c] == 0)), None)
        if piv is None:
            return one * 0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            d = -d
        p = A[c][c]
        d = d * p
        for i in range(c + 1, n):
            if not (A[i][c] == 0):
                f = A[i][c] / p
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return d


def parse(spec):
    if spec.endswith('+'):
        return None, int(spec[:-1])
    return int(spec), None


def build(fam, s0, s1, log=print):
    (f0, o0), (f1, o1) = parse(s0), parse(s1)
    names = ['q'] + (['x0'] if f0 is None else []) + (['x1'] if f1 is None else [])
    K = Field(names)
    g = dict(zip(names, K.gens))
    one = K(1)
    u = [g['x0'] + o0 if f0 is None else K(f0), g['x1'] + o1 if f1 is None else K(f1)]
    ctx = r2v.Ctx(one, g['q'], u[0], u[1])
    F = r2v.fam_by_name(fam)
    sample = (o0 + 3 if f0 is None else f0, o1 + 3 if f1 is None else f1)
    bases = r2v.choose(F, sample)
    red = r2v.reduce_family(F, ctx, bases, full=True, log=log)
    return K, ctx, F, bases, red


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if '=' not in a]
    opts = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
    fam, s0, s1 = args
    K, ctx, F, bases, red = build(fam, s0, s1, log=lambda m: print(m, flush=True))
    one = K(1)
    nK, nW, nQ = (len(b) for b in bases)
    print(f'{fam}, u0 = {s0}, u1 = {s1}: |K| = {nK}, |W| = {nW}, |Q| = {nQ} [{time.time() - t0:.1f}s]', flush=True)
    ok = True
    if nW == 0:
        print('  the irrep does not occur in the W-part here (nothing to certify)')
        print('PASS')
        sys.exit(0)
    for gname in ('GP', 'GQ'):
        if gname in red:
            sg, desc = signcheck(det(red[gname], one))
            print(f'  det {gname}: sign {sg} ({desc})', flush=True)
            ok &= sg == 1
    Sf = red['Sfull']
    n = len(Sf)
    sym = all(Sf[i][j] == Sf[j][i] for i in range(n) for j in range(i + 1, n))
    kz = all(Sf[i][j] == 0 for i in range(nK) for j in range(n))
    print(f'  reduced S\' symmetric: {sym}; K-rows vanish identically: {kz} [{time.time() - t0:.1f}s]', flush=True)
    ok &= sym and kz
    S = red['S']
    A = [list(r) for r in S]
    m = len(A)
    for c in range(m):
        p = A[c][c]
        sg, desc = signcheck(p)
        print(f'  pivot {c + 1}/{m} ({red["names"][1][c]}): sign {sg}; numerator deg {p.n.total_degree()} '
              f'({len(p.n)} terms), denominator deg {p.d.total_degree()} ({len(p.d)} terms); {desc} '
              f'[{time.time() - t0:.1f}s]', flush=True)
        if sg != 1:
            ok = False
            break
        for i in range(c + 1, m):
            if not (A[i][c] == 0):
                f = A[i][c] / p
                for j in range(c + 1, m):
                    A[i][j] = A[i][j] - f * A[c][j]
    if 'save' in opts:
        with open(opts['save'], 'wb') as fh:
            pickle.dump(dict(fam=fam, s0=s0, s1=s1, names=K.ctx.names(), basis=red['names'],
                             S=[[(str(x.n), str(x.d)) for x in row] for row in S]), fh)
    print(f'time {time.time() - t0:.1f}s')
    print('PASS' if ok else 'FAIL')
