"""Worker for the determinant route of the generic two-row family [q0 >= 4, t >= 4, i >= 7].

At a slice (t, i) = (4 + y, 7 + z), y, z integers >= 0, compute det Q over Q(q0, s) (Q the gauge-fixed two-row form of
g2lib.build, specialised at t and i), multiply by Den_c (the explicit candidate denominator, DEN_EXP below), check
that the product is a polynomial, substitute q0 = 4 + x, and either
  - mode 'save': write its coefficients {(a, b): (p, q)} (x^a s^b, coefficient p/q) to OUT/s_y_z.pkl, or
  - mode 'check': compare it with the interpolated polynomial (OUT/interp.pkl) at (y, z) and print PASS / FAIL.
Slices are taken in column order y, then z; worker W of NW takes the columns y with y % NW == W.
python g2slices.py MODE OUT Y0 Y1 Z0 Z1 W NW"""
import sys, os, time, pickle
import flint
import g2core
from flintfield import F

# candidate exponents of the 12 denominator factors of det Q, in the order of g2core.Slicer.factors (observed on
# slices; the verification in mode 'check' proves that det Q * Den_c is a polynomial)
DEN_EXP = [6, 4, 4, 4, 4, 2, 1, 2, 2, 1, 1, 1]
ZLOW = int(os.environ.get('G2_ZLOW', '60'))


def den_c(S):
    D = S.C.K.ctx.from_dict({(0,) * 4: 1})
    for f, e in zip(S.factors, DEN_EXP):
        D = D * f ** e
    return D


def slice_poly(S, Dc, y, z):
    """the polynomial det Q * Den_c at t = 4 + y, i = 7 + z, after q0 = 4 + x, as {(a, b): fmpq}"""
    free = ('q0', 's')
    vals = {'t': 4 + y, 'i': 7 + z}
    d = S.det(free, vals)
    Dc2 = S.spec(Dc, free, vals)
    q, r = divmod(Dc2, d.d)
    if r != 0:
        return None
    P = d.n * q
    K = S.field(free)
    x, s = K.ctx.gens()
    P = P.compose(x + 4, s, ctx=K.ctx)
    return {tuple(e): c for e, c in P.to_dict().items()}


def main():
    mode, out = sys.argv[1], sys.argv[2]
    y0, y1, z0, z1, w, nw = (int(v) for v in sys.argv[3:9])
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    S = g2core.Slicer(check=False, log=lambda m: None)
    assert len(S.factors) == len(DEN_EXP)
    Dc = den_c(S)
    interp = None
    if mode == 'check':
        import g2interp
        interp = g2interp.Evaluator(os.path.join(out, 'interp.pkl'))
    ok = True
    nd = 0
    ys = range(y0, y1 + 1)
    if os.environ.get('G2_REVERSE') == '1':
        ys = range(y1, y0 - 1, -1)
    for y in ys:
        if y % nw != w:
            continue
        for z in range(z0, z1 + 1):
            fn = os.path.join(out, f's_{y}_{z}.pkl')
            if mode == 'save' and os.path.exists(fn):
                continue
            if mode == 'check' and interp.in_grid(y, z):
                continue
            t1 = time.time()
            P = slice_poly(S, Dc, y, z)
            if P is None:
                print(f'y={y} z={z}: det * Den_c NOT A POLYNOMIAL', flush=True)
                ok = False
                continue
            if mode == 'save':
                data = {e: (int(c.p), int(c.q)) for e, c in P.items()}
                with open(fn + '.tmp', 'wb') as fh:
                    pickle.dump(data, fh, protocol=4)
                os.replace(fn + '.tmp', fn)
                res = 'saved'
            else:
                good = interp.equals(y, z, P)
                ok &= good
                res = 'PASS' if good else 'FAIL'
            nd += 1
            print(f'y={y} z={z}: {len(P)} terms, {res} [{time.time()-t1:.1f}s, total {time.time()-t0:.0f}s]',
                  flush=True)
        if mode == 'save' and z0 == 0 and not os.path.exists(os.path.join(out, f'col_{y}.pkl')):
            import g2interp
            ok &= g2interp.stage_a(out, y, z1, ZLOW)
    print(f'worker {w}/{nw}: {nd} slices, {time.time()-t0:.0f}s')
    print('PASS' if ok else 'FAIL')


if __name__ == '__main__':
    main()
