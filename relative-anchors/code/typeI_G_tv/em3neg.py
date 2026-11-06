"""Negative controls for Section 10 (exact).
 1. A larger Gershgorin constant: the third pivot fails.  For the first H in H_LIST with such a point, an explicit
    rational point with u1, u2 >= 17 where the first two pivots of A are positive and det A < 0 (closed form evaluated
    over Q, as in em3points.py; at H = 1/2 all three pivots are positive there).  On the slice through that point,
    det A at that H has a numerator with coefficients of both signs, so the sign certificate cannot hold.
 2. The vanishing test of the interpolation is not vacuous: if one slice (x1 = 38) is multiplied by 1 + 1/1000, the
    interpolant of the first 1500 monomials has nonzero coefficients above x1^L; unperturbed, it has none.
 3. Every exponent of Den_c is needed: lowering any one of them by 1, Den_c det A is no longer a polynomial on the
    slice x1 = 0.
python em3neg.py"""
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
import em3lib as L
from em3det import interpolate, load_setup, OUT
from em3points import load_closed

Q = flint.fmpq
H_LIST = (Q(1), Q(3, 2), Q(2), Q(3, 4))


def pivots(names, G, u1, u2, s1, st, rr, h):
    sig = st + rr
    vals = {'u1': u1, 'u2': u2, 'sig': sig, 's1': s1}
    args = [vals[x] for x in names]
    Gz = [[G[(min(i, j), max(i, j))][0](*args) / G[(min(i, j), max(i, j))][1](*args) for j in range(3)]
          for i in range(3)]
    S1, S2 = u1 * s1, u2
    S = S1 + S2
    Nm = [[Q(1), Q(0), Q(0)], [Q(0), 1 / S, S2 / S], [Q(0), 1 / S, -S1 / S]]
    Gy = [[sum((Nm[k][i] * Gz[k][l] * Nm[l][j] for k in range(3) for l in range(3)), Q(0)) for j in range(3)]
          for i in range(3)]
    kap = [Q(1), st * sig / (S * (sig + S)), st * S1 * S2 / S]
    M = [[Gy[i][j] - (h * (rr / sig) * kap[i] if i == j else 0) for j in range(3)] for i in range(3)]
    p1 = M[0][0]
    b11 = M[1][1] - M[0][1] * M[0][1] / p1
    b12 = M[1][2] - M[0][1] * M[0][2] / p1
    b22 = M[2][2] - M[0][2] * M[0][2] / p1
    det = p1 * (b11 * b22 - b12 * b12)
    return p1, b11, det


def control_H():
    """for each type, the first H in H_LIST (all > 1/2) with a rational point (u1, u2 >= 17) where pivots 1 and 2 are
    positive and det A < 0"""
    ok = True
    for ka, md in (('s', 'T11'), ('1', 'T22')):
        names, G = load_closed(md)
        wit = None
        for h in H_LIST:
            for (u1, u2) in ((17, 17), (17, 50), (50, 17), (50, 50), (200, 200)):
                for s1 in (Q(1, 1000), Q(1, 10), Q(1), Q(10), Q(1000)):
                    for rr in (Q(1, 1000), Q(2), Q(50), Q(1000)):
                        st = s1 if ka == 's' else Q(1)
                        p1, p2, det = pivots(names, G, Q(u1), Q(u2), s1, st, rr, h)
                        if p1 > 0 and p2 > 0 and det < 0:
                            wit = (h, u1, u2, s1, rr, p1, p2, det)
                            break
                    if wit:
                        break
                if wit:
                    break
            if wit:
                break
        if not wit:
            print('s_t = %s: no witness found for H in %s' % (ka, [str(h) for h in H_LIST]), flush=True)
            ok = False
            continue
        h, u1, u2, s1, rr, p1, p2, det = wit
        st = s1 if ka == 's' else Q(1)
        q1, q2, qd = pivots(names, G, Q(u1), Q(u2), s1, st, rr, Q(1, 2))
        good_half = q1 > 0 and q2 > 0 and qd > 0
        print('s_t = %s: at (u1, u2, s, r) = (%d, %d, %s, %s) with H = %s: pivots 1, 2 = %.4e, %.4e > 0 and det A = '
              '%.4e < 0 (third pivot negative); with H = 1/2 all pivots positive there: %s'
              % (ka, u1, u2, s1, rr, h, float(p1), float(p2), float(det), good_half), flush=True)
        ok &= good_half
        # the slice through the witness at this H: numerator coefficients of both signs
        t0 = time.time()
        A, (nm, Ao) = L.build_A(ka, H=h, orig=True, log=lambda *a, **k: None)
        OS = L.OrigSlicer(ka, nm, Ao)
        c = u1 - 2
        dt = OS.det(c)
        cf = list(dt.n.coeffs())
        neg = sum(1 for x in cf if x < 0)
        pos = len(cf) - neg
        val = dt.n(Q(u2 - 2), rr, s1) / dt.d(Q(u2 - 2), rr, s1)
        print('  slice x1 = %d at H = %s: numerator %d terms, %d positive and %d negative coefficients; its value at '
              'the witness equals the direct det: %s [%.0fs]' % (c, h, len(cf), pos, neg, val == det,
                                                                 time.time() - t0), flush=True)
        ok &= pos > 0 and neg > 0 and val == det
        del A, Ao, OS, dt
    return ok


def control_interp(ka='s', cpert=38, nmon=1500):
    st = load_setup(ka)
    B, low = st['B'], st['L']
    slices = []
    for c in range(B + 1):
        with open(os.path.join(OUT, ka, 's_%d.pkl' % c), 'rb') as fh:
            slices.append(pickle.load(fh))
    mons = set()
    for d in slices:
        mons.update(d.keys())
    part = set(sorted(mons)[:nmon])
    sub = [{m: v for m, v in d.items() if m in part} for d in slices]
    _, nbad0, _ = interpolate(sub, B, low, log=lambda *a, **k: None)
    pert = [dict(d) for d in sub]
    f = Q(1001, 1000)
    for m, (p, q) in list(pert[cpert].items()):
        v = Q(p, q) * f
        pert[cpert][m] = (int(v.p), int(v.q))
    _, nbad1, _ = interpolate(pert, B, low, log=lambda *a, **k: None)
    print('s_t = %s, first %d monomials: nonzero coefficients above x1^%d: %d as computed, %d with the slice x1 = %d '
          'multiplied by 1 + 1/1000' % (ka, len(part), low, nbad0, nbad1, cpert), flush=True)
    return nbad0 == 0 and nbad1 > 0


def control_den(ka='s'):
    st = load_setup(ka)
    with open(os.path.join(OUT, ka, 's_0.pkl'), 'rb') as fh:
        d = pickle.load(fh)
    P0 = L.K3.ctx.from_dict({m: Q(p, q) for m, (p, q) in d.items()})
    facs = [L.KX.ctx.from_dict({m: Q(p, q) for m, (p, q) in fd.items()}) for fd in st['fac_dicts']]
    res = []
    for k, (f, e) in enumerate(zip(facs, st['dexp'])):
        if e == 0:
            continue
        q, rem = divmod(P0, L.spec(f, 0))
        res.append(rem != 0)
    print('s_t = %s: lowering the exponent of each of the %d factors of Den_c by 1, Den_c det A is not a polynomial on '
          'the slice x1 = 0: %s' % (ka, len(res), all(res)), flush=True)
    return all(res)


def main():
    t0 = time.time()
    ok = control_H()
    for ka in ('s', '1'):
        ok &= control_interp(ka)
        ok &= control_den(ka)
    print('every control behaves as expected: %s [%.0fs]' % (ok, time.time() - t0))
    return ok


if __name__ == '__main__':
    print('PASS' if main() else 'FAIL')
