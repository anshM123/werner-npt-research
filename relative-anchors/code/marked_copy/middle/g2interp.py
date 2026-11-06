"""Exact interpolation for the determinant route (see g2slices.py).

Stage A (one column y): the slices s_y_z.pkl, z = 0..ZMAX, give for every monomial x^a s^b its value at z; the
interpolant in z of degree <= ZMAX is computed exactly (inverse Vandermonde over Q); its coefficients of z^l with
l > ZLOW must vanish; the coefficients l = 0..ZLOW are written to col_y.pkl.
Stage B (all columns): for every monomial and every l <= ZLOW the values at y = 0..YMAX are interpolated in y; the
coefficients of y^k with k > YLOW must vanish; the coefficients k <= YLOW are the coefficients of x^a y^k z^l s^b of
the polynomial N~(x, y, z, s); every one must be >= 0, and some coefficient with a = k = l = 0 must be > 0.
python g2interp.py A OUT Y ZMAX ZLOW      python g2interp.py B OUT YMAX YLOW ZLOW [CHUNK]"""
import sys, os, time, pickle
import flint


def inv_vandermonde(D):
    V = flint.fmpq_mat(D + 1, D + 1, [flint.fmpq(z) ** l for z in range(D + 1) for l in range(D + 1)])
    return V.inv()


def stage_a(out, y, zmax, zlow):
    t0 = time.time()
    cols = []
    mons = set()
    for z in range(zmax + 1):
        with open(os.path.join(out, f's_{y}_{z}.pkl'), 'rb') as fh:
            d = pickle.load(fh)
        cols.append(d)
        mons.update(d.keys())
    mons = sorted(mons)
    nm = len(mons)
    Y = flint.fmpq_mat(zmax + 1, nm)
    for z, d in enumerate(cols):
        for j, mo in enumerate(mons):
            if mo in d:
                p, q = d[mo]
                Y[z, j] = flint.fmpq(p, q)
    Cz = inv_vandermonde(zmax) * Y
    bad = 0
    for l in range(zlow + 1, zmax + 1):
        for j in range(nm):
            if Cz[l, j] != 0:
                bad += 1
    data = {mo: [(int(Cz[l, j].p), int(Cz[l, j].q)) for l in range(zlow + 1)] for j, mo in enumerate(mons)}
    with open(os.path.join(out, f'col_{y}.pkl.tmp'), 'wb') as fh:
        pickle.dump({'ok': bad == 0, 'bad': bad, 'data': data, 'zmax': zmax, 'zlow': zlow}, fh, protocol=4)
    os.replace(os.path.join(out, f'col_{y}.pkl.tmp'), os.path.join(out, f'col_{y}.pkl'))
    print(f'column y={y}: {nm} monomials, interpolated in z (degree <= {zmax}); nonzero coefficients above '
          f'z^{zlow}: {bad} [{time.time()-t0:.1f}s]', flush=True)
    return bad == 0


def stage_b(out, ymax, ylow, zlow, chunk=200):
    t0 = time.time()
    Vy = inv_vandermonde(ymax)
    # the monomial set: union over the columns
    mons = set()
    oks = True
    for y in range(ymax + 1):
        with open(os.path.join(out, f'col_{y}.pkl'), 'rb') as fh:
            c = pickle.load(fh)
        assert c['zlow'] == zlow
        oks &= c['ok']
        mons.update(c['data'].keys())
        del c
    mons = sorted(mons)
    print(f'{len(mons)} monomials x^a s^b; every column vanishes above z^{zlow}: {oks} [{time.time()-t0:.0f}s]',
          flush=True)
    nbad = nneg = npos = nzero = 0
    base_pos = False
    maxdeg = [0, 0, 0, 0]
    # values of N~ at test points (x, y, z, s) outside the grid, for an independent comparison (g2points.py)
    tpts = [(flint.fmpq(3), flint.fmpq(131), flint.fmpq(157), flint.fmpq(2, 7)),
            (flint.fmpq(40), flint.fmpq(250), flint.fmpq(9), flint.fmpq(5)),
            (flint.fmpq(0), flint.fmpq(0), flint.fmpq(400), flint.fmpq(1, 1000)),
            (flint.fmpq(17), flint.fmpq(1000), flint.fmpq(1000), flint.fmpq(1000))]
    tval = [flint.fmpq(0) for _ in tpts]
    for c0 in range(0, len(mons), chunk):
        part = mons[c0:c0 + chunk]
        # values[l] : (ymax+1) x len(part)
        vals = [flint.fmpq_mat(ymax + 1, len(part)) for _ in range(zlow + 1)]
        for y in range(ymax + 1):
            with open(os.path.join(out, f'col_{y}.pkl'), 'rb') as fh:
                data = pickle.load(fh)['data']
            for j, mo in enumerate(part):
                lst = data.get(mo)
                if lst is None:
                    continue
                for l, (p, q) in enumerate(lst):
                    if p:
                        vals[l][y, j] = flint.fmpq(p, q)
            del data
        for l in range(zlow + 1):
            Cy = Vy * vals[l]
            for j, mo in enumerate(part):
                for k in range(ymax + 1):
                    c = Cy[k, j]
                    if c == 0:
                        if k <= ylow:
                            nzero += 1
                        continue
                    if k > ylow:
                        nbad += 1
                        continue
                    a, b = mo
                    for ip, (px, py, pz, ps) in enumerate(tpts):
                        tval[ip] += c * px ** a * py ** k * pz ** l * ps ** b
                    if c < 0:
                        nneg += 1
                        if nneg <= 20:
                            print(f'  NEGATIVE coefficient of x^{a} y^{k} z^{l} s^{b}: {float(c):.3e}', flush=True)
                    else:
                        npos += 1
                        maxdeg = [max(maxdeg[0], a), max(maxdeg[1], k), max(maxdeg[2], l), max(maxdeg[3], b)]
                        if a == 0 and k == 0 and l == 0:
                            base_pos = True
        print(f'  monomials {c0}..{c0 + len(part) - 1}: so far {npos} positive, {nneg} negative, {nbad} nonzero above '
              f'y^{ylow} [{time.time()-t0:.0f}s]', flush=True)
    ok = oks and nbad == 0 and nneg == 0 and base_pos
    with open(os.path.join(out, 'testpoints.pkl'), 'wb') as fh:
        pickle.dump([(tuple((int(v.p), int(v.q)) for v in pt), (int(v.p), int(v.q))) for pt, v in zip(tpts, tval)], fh)
    for pt, v in zip(tpts, tval):
        print(f'  N~ at (x, y, z, s) = {tuple(str(c) for c in pt)}: {float(v):.6e}')
    print(f'N~(x, y, z, s): {npos} positive coefficients, {nneg} negative, {nzero} zero (inside the box); degrees '
          f'(x, y, z, s) <= {tuple(maxdeg)}; nonzero above y^{ylow}: {nbad}; a positive coefficient free of x, y, z: '
          f'{base_pos} [{time.time()-t0:.0f}s]')
    return ok


class Evaluator:
    """unused placeholder for g2slices mode 'check'"""
    def __init__(self, fn):
        raise NotImplementedError


if __name__ == '__main__':
    if sys.argv[1] == 'A':
        out, y, zmax, zlow = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        print('PASS' if stage_a(out, y, zmax, zlow) else 'FAIL')
    else:
        out, ymax, ylow, zlow = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        chunk = int(sys.argv[6]) if len(sys.argv) > 6 else 200
        print('PASS' if stage_b(out, ymax, ylow, zlow, chunk) else 'FAIL')
