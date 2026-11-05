"""Three levels at fixed rational weights, every u3 at once: two single U-copies (weights s1 = w1^2, s2 = w2^2) and u3
copies of weight s3 = w3^2, sigma = s_t = 1.  For each sampled rational (w1, w2, w3) the invariant reduction of
three_level_fast.py is carried out exactly over Q(u3) (one variable, fast), and with u3 = 2 + x every numerator and
denominator below is checked to have coefficients of one sign:
  (Z1) G_z[rho,h] < 0;  (Z3) G_z[h,h'] < 0 for h, h' in different levels;  (D) G_z v > 0.
So at each sampled weight triple these hold for every u3 >= 2 (pairs inside the u3 equal copies are not covered).
This is a family of exact checks (pointwise in the weights, uniform in u3), not a proof.
The reduction is also compared with direct elimination (gzexact.py) at u3 = 2, 3, 5 for every triple."""
import sys, os, time, random
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sympy.polys.domains import QQ
from sympy.polys.fields import field
import flint
import three_level_fast as T3
from gzexact import Gz_exact

Q = flint.fmpq


def reduce_at(w1, w2, w3):
    """Q4 as 4 x 4 list of flint fmpq_poly pairs (num, den) in u3, at sigma = 1 and the given rational weights."""
    Kf, U3s = field('u3', QQ)
    T3.Kf, T3.U3s = Kf, U3s
    T3.HALF = Kf.one / 2
    T3.WT = Kf.one
    T3.O1, T3.O2, T3.O3 = Kf(QQ(w1.p, w1.q)), Kf(QQ(w2.p, w2.q)), Kf(QQ(w3.p, w3.q))
    T3.OM = {1: T3.O1, 2: T3.O2, 3: T3.O3}
    T3.SS = {i: T3.OM[i] ** 2 for i in T3.LV}
    T3.UU = {1: Kf.one, 2: Kf.one, 3: U3s}
    T3.SIG = T3.WT ** 2
    R, Phi, elim, keep = T3.make()
    allv = elim + keep
    idx = {R.gens.index(v): k for k, v in enumerate(allv)}
    n = len(allv)
    H = [[Kf.zero] * n for _ in range(n)]
    for mon, co in Phi.terms():
        vs = []
        for pos, e in enumerate(mon):
            vs += [idx[pos]] * e
        a, b = vs
        if a == b:
            H[a][a] += co
        else:
            H[a][b] += co / 2
            H[b][a] += co / 2
    # elimination over Q(u3) in SymPy's field (univariate: cheap)
    ne = len(elim)
    for p in range(ne):
        pv = H[p][p]
        assert pv != 0
        for i in range(p + 1, n):
            if H[i][p] == 0:
                continue
            f = H[i][p] / pv
            for j in range(p + 1, n):
                if H[p][j] != 0:
                    H[i][j] = H[i][j] - f * H[p][j]
    return [[H[ne + i][ne + j] for j in range(4)] for i in range(4)], Kf, U3s


def to_poly(p):
    """univariate SymPy PolyElement in u3 -> list of Fractions (coefficient of u3^k at index k)."""
    d = {}
    for mon, co in p.terms():
        d[mon[0]] = Fraction(str(co))
    deg = max(d) if d else 0
    return [d.get(k, Fraction(0)) for k in range(deg + 1)]


def shift2(c):
    """coefficients of p(2 + x) from those of p(u)."""
    out = [Fraction(0)] * len(c)
    # Horner: p(2 + x)
    for a in reversed(c):
        # out = out * (2 + x) + a
        new = [Fraction(0)] * len(c)
        for k, b in enumerate(out):
            if b:
                new[k] += 2 * b
                if k + 1 < len(new):
                    new[k + 1] += b
        new[0] += a
        out = new
    return out


def sign(c):
    nz = [x for x in c if x != 0]
    return 1 if nz and all(x > 0 for x in nz) else (-1 if nz and all(x < 0 for x in nz) else 0)


if __name__ == '__main__':
    npts = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    rng = random.Random(7)
    ok = True
    t0 = time.time()
    for it in range(npts):
        def rw():
            e = rng.uniform(-2.5, 2.5)
            if it % 4 == 1:
                e = rng.uniform(-6, -3)     # first single copy tiny
            return Q(round(10 ** e * 1000) or 1, 1000)
        w1 = rw()
        w2 = Q(round(10 ** rng.uniform(-2, 2) * 1000) or 1, 1000)
        w3 = Q(round(10 ** rng.uniform(-2, 2) * 1000) or 1, 1000)
        if it % 4 == 2:
            w3 = w3 / 1000                  # the equal copies tiny
        Q4, Kf, U3s = reduce_at(w1, w2, w3)
        # validation at u3 = 2, 3, 5
        good = True
        for a3 in (2, 3, 5):
            G = Gz_exact(Q(1), [w1, w2] + [w3] * a3)
            idx = [[0], [1], [2], list(range(3, 3 + a3))]
            for p in range(4):
                for q in range(4):
                    direct = sum((G[i, j] for i in idx[p] for j in idx[q]), Q(0))
                    num = sum((Fraction(str(co)) * a3 ** mon[0] for mon, co in Q4[p][q].numer.terms()), Fraction(0))
                    den = sum((Fraction(str(co)) * a3 ** mon[0] for mon, co in Q4[p][q].denom.terms()), Fraction(0))
                    good &= Fraction(int(direct.p), int(direct.q)) == num / den
        # certificates in u3 = 2 + x
        S = w1 * w1 + w2 * w2
        s3 = w3 * w3
        signs = []
        for (p, q) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]:
            e = Q4[p][q]
            signs.append(sign(shift2(to_poly(e.numer))) * sign(shift2(to_poly(e.denom))) == -1)
        # v = (4 (1 + S_U), 3 (2 + S_U), ...), S_U = s1 + s2 + u3 s3  (linear in u3)
        SU = Kf(QQ(int(S.p), int(S.q))) + U3s * Kf(QQ(int(s3.p), int(s3.q)))
        v = [4 * (1 + SU)] + [3 * (2 + SU)] * 3
        for p in range(4):
            row = sum((Q4[p][q] * v[q] for q in range(4)), Kf.zero)
            signs.append(sign(shift2(to_poly(row.numer))) * sign(shift2(to_poly(row.denom))) == 1)
        allg = good and all(signs)
        ok &= allg
        print(f'point {it}: w = ({w1}, {w2}, {w3}): reduction = direct elimination at u3 = 2, 3, 5: {good}; '
              f'(Z1), (Z3), (D) for every u3 >= 2: {all(signs)}  ({time.time() - t0:.1f} s)', flush=True)
    print('PASS' if ok else 'FAIL')
