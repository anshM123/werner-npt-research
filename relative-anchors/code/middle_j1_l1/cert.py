"""Leading-minor sign certificates for the reduced forms of l1sym.py (letter-1 block of (R_inf)_1, two weight levels).

python cert.py PART U0 U1 [bscale] [order=i,j,...] [save=file]
  PART in inv, std, two, anti, mix (irrep on level 0; 'mix' on both levels);
  U0, U1: a fixed level size (an integer, e.g. 1 or 2) or 'o+' for every size u = o + x, x >= 0 real (o >= 3).
Squared weights Q0 = 1, Q1 = q (the rescaled forms are homogeneous in the weights, so this is no restriction).
Steps: build the reduced forms over Q(q, x0, x1) (flint), check that the Gram determinants of the test vectors are
positive on the regime (so the vectors are a basis of the multiplicity space), eliminate the vacuum states (Schur complement S onto the
P-vectors; the vacuum block is positive definite by Lemma G), check S kappa = 0 for the kernel vectors kappa of the part,
delete the gauge coordinates (one per kernel vector, with coefficient 1 in it and 0 in the others), and factor every
leading principal minor of the gauge-fixed S over Q.  A minor passes if every irreducible factor of its numerator and
denominator has coefficients of one sign and a monomial free of the x_i, and the resulting sign is +1; then the minor
is positive for all q > 0 and all real x_i >= 0, and by Sylvester's criterion the gauge-fixed S is positive definite.
Prints PASS if every minor passes (and the kernel identities hold), FAIL otherwise."""
import sys, time, pickle
from fractions import Fraction as Fr
import flint
from flintfield import Field
from l1sym import Model, reduce_part, keep_for, det, num

t0 = time.time()


def signcheck(f):
    """(sign, description) if f = c prod g_i^e_i with every irreducible g_i one-signed in coefficients and with a
    monomial free of the x-variables (so g_i != 0 for q > 0, x >= 0); sign is then the sign of f there.
    (None, reason) otherwise.  (Copied from middle_j1_t1/lev.py.)"""
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


def parse_u(spec):
    if spec.endswith('+'):
        return None, int(spec[:-1])
    return int(spec), None


def build(part, s0, s1, bsc=Fr(1)):
    (f0, o0), (f1, o1) = parse_u(s0), parse_u(s1)
    names = ['q'] + (['x0'] if f0 is None else []) + (['x1'] if f1 is None else [])
    K = Field(names)
    g = dict(zip(names, K.gens))
    one = K(1)
    U = [g['x0'] + o0 if f0 is None else K(f0), g['x1'] + o1 if f1 is None else K(f1)]
    M = Model(one, [one, g['q']], U, bsc)
    kp, kq, gauge = keep_for(part, f0 if f0 is not None else o0, f1 if f1 is not None else o1)
    red = reduce_part(M, part, kp, kq)
    return K, M, red, gauge


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if '=' not in a]
    opts = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
    part, s0, s1 = args[0], args[1], args[2]
    bsc = Fr(args[3]) if len(args) > 3 else Fr(1)
    K, M, red, gauge = build(part, s0, s1, bsc)
    S = red['S']
    n = len(S)
    one = K(1)
    print(f'{part}, u0 = {s0}, u1 = {s1}, bscale {bsc}: S of size {n} ({red["names"][0]}), Q-vectors '
          f'{red["names"][1]} eliminated [{time.time() - t0:.1f}s]', flush=True)
    ok = all(S[i][j] == S[j][i] for i in range(n) for j in range(n))
    print(f'  S symmetric: {ok}', flush=True)
    # the test vectors are a basis of their multiplicity space on the whole regime: Gram determinants > 0
    for gname in ('GP', 'GQ'):
        if gname in red:
            sg, desc = signcheck(det(red[gname], one))
            print(f'  det of the Gram matrix {gname}: sign {sg} ({desc})', flush=True)
            ok &= sg == 1
    if bsc == 1:
        for kv in red['ker']:
            z = all(sum((S[i][j] * kv[j] for j in range(n)), K(0)) == 0 for i in range(n))
            print(f'  kernel vector {kv}: S kappa = 0 identically: {z}', flush=True)
            ok &= z
    keep = [i for i in range(n) if i not in gauge]
    if 'order' in opts:
        o = [int(v) for v in opts['order'].split(',')]
        assert sorted(o) == sorted(keep), (o, keep)
        keep = o
    Sg = [[S[i][j] for j in keep] for i in keep]
    print(f'  gauge-fixed coordinates (order): {[red["names"][0][i] for i in keep]}', flush=True)
    minors = []
    for m in range(1, len(keep) + 1):
        dm = det([row[:m] for row in Sg[:m]], one)
        minors.append(dm)
        sg, desc = signcheck(dm)
        good = sg == 1
        ok &= good
        print(f'  minor {m}: sign {sg}; numerator deg {dm.n.total_degree()} ({len(dm.n)} terms), denominator deg '
              f'{dm.d.total_degree()} ({len(dm.d)} terms); factors (deg/terms^mult): {desc} [{time.time() - t0:.1f}s]',
              flush=True)
        if not good:
            break
    if 'save' in opts:
        with open(opts['save'], 'wb') as fh:
            pickle.dump(dict(names=K.ctx.names() if hasattr(K.ctx, 'names') else None,
                             minors=[(str(d.n), str(d.d)) for d in minors]), fh)
    print(f'time {time.time() - t0:.1f}s')
    print('PASS' if ok else 'FAIL')
