"""Route (i) of NOTES Section 6: an explicit lift of a whole block, and the test whether it is a fixed-vector matrix.

In the rescaled basis (ylib) the B-part of a block is crossed:
    x^T B' x = sum_{a,b in rows} <V_ab(x), D'^{-1} V_ba(x)>,   V_ab(x) = sum_{y'} x_(a,y') MV'(b, y'),
with MV' and D' = mu D mu as in ylib.  Split  sum_{a,b} = sum_a (a = b) + sum_{a<b} [ (1/2)<S, D'^{-1} S> - (1/2)<Dl, D'^{-1} Dl> ],
S = V_ab + V_ba, Dl = V_ab - V_ba.  The terms -<V_aa, D'^{-1} V_aa> and -(1/2)<S, D'^{-1} S> are Schur complements
(variables z_a, zp_ab with matrix D').  The wrong-sign term +(1/2)<Dl, D'^{-1} Dl> is rewritten with Woodbury,
D' = C' + H K H^T (C' = 2^{-(k-c)} s^{N(z)} C_z diagonal, H[u + 1_t, u] = s_t, K = diag(2^{-(k-c)-1} s^{N(u)})):
    <Dl, D'^{-1} Dl> = min_zeta [ <Dl + H zeta, C'^{-1} (Dl + H zeta)> + <zeta, K^{-1} zeta> ].
So Y' = Schur complement onto x of the lift
    Q = A'(x) + sum_a [<z_a, D' z_a> + 2<z_a, V_aa>] + (1/2) sum_{a<b} [<zp, D' zp> + 2<zp, S>]
        + (1/2) sum_{a<b} [ <Dl + H zeta, C'^{-1}(Dl + H zeta)> + <zeta, K^{-1} zeta> ]
(the minus sign of -B' is absorbed: min over z_a of <z,D'z> + 2<z,V> is -<V, D'^{-1} V>).
Checks: (1) exact, at a rational point (w rational): the Schur complement of the lift onto the kept coordinates x equals
the gauge-fixed block; D' = C' + H K H^T.  (2) symbolic: the lift matrix L(s); exact obstruction (newton test) and, if
unobstructed or rescalable, the monomial-PSD test.
Usage: python lift_test.py k blockname"""
import sys, itertools, time
import flint
from ylib import Builder, Ring, point_ring, Mel, Nv, halfex, wcol
from flintfield import Field, F
from psdtools import load_blocks
from analyze import obstruction, mono_psd, farkas, rescale_lp, apply_rescale

k = int(sys.argv[1])
bname = sys.argv[2]
t0 = time.time()
blk = [b for b in load_blocks(k) if b['name'] == bname][0]
states = blk['states']           # kept states
rows = sorted({st[0] for st in states})
vac = (0,) * k


def Dprime_blocks(Bd):
    """D' blocks (as in ylib.Builder.Dinv, without inversion): dict wcol -> (list of strings, matrix)."""
    R = Bd.R
    out = {}
    blocks = {}
    for y in Bd.cols:
        blocks.setdefault(wcol(y), []).append(y)
    for key, ys in blocks.items():
        Bm = []
        for y in ys:
            row = []
            for yp in ys:
                v = R.zero
                for t in range(k):
                    for tp in range(k):
                        mv = Mel(Bd.e[t], y, Bd.e[tp], yp)
                        if mv:
                            v += R.one * mv * R.mono(halfex(Nv(y), Nv(yp), Bd.e[t], Bd.e[tp]))
                row.append(v)
            Bm.append(row)
        out[key] = (ys, Bm)
    return out


def build_lift(Bd, s, zero, one):
    """returns (index list, matrix Lmat as dict of dicts, x-indices) of the lift in the ring of Bd."""
    Dp = Dprime_blocks(Bd)
    idx = {}
    def var(name):
        if name not in idx:
            idx[name] = len(idx)
        return idx[name]
    for st in states:
        var(('x', st))
    L = {}
    def add(i, j, v):
        if v == 0:
            return
        L.setdefault(i, {})
        L[i][j] = L[i].get(j, zero) + v
    # A' part on x
    for X1 in states:
        for X2 in states:
            v = Bd.A(X1[0], X1[1], X2[0], X2[1])
            add(idx[('x', X1)], idx[('x', X2)], v)
    # V_ab as dict: z -> dict(x index -> coefficient)
    def V(a, b):
        out = {}
        for st in states:
            if st[0] != a:
                continue
            for z, c in Bd.MV(b, st[1]).items():
                out.setdefault(z, {})
                out[z][idx[('x', st)]] = out[z].get(idx[('x', st)], zero) + c
        return out

    def sector_strings(vec):
        keys = {wcol(z) for z in vec}
        res = []
        for key in keys:
            res.extend(Dp[key][0])
        return keys

    def add_schur_pair(tag, vec, scale):
        """adds scale * ( <z, D' z> + 2 <z, vec(x)> ) with new variables z over the D'-blocks touched by vec."""
        for key in sector_strings(vec):
            ys, Bm = Dp[key]
            zi = [var((tag, key, y)) for y in ys]
            for p in range(len(ys)):
                for q in range(len(ys)):
                    add(zi[p], zi[q], scale * Bm[p][q])
            for p, y in enumerate(ys):
                for xi, c in vec.get(y, {}).items():
                    add(zi[p], xi, scale * c)
                    add(xi, zi[p], scale * c)

    def add_woodbury(tag, vec, scale):
        """adds scale * ( <vec + H zeta, C'^{-1} (vec + H zeta)> + <zeta, K^{-1} zeta> ), zeta over H_{c-1}."""
        for key in sector_strings(vec):
            ys, Bm = Dp[key]
            c = key[0] + sum(key[1:])            # number of nonzero letters
            two = one * 2 ** (k - c)
            for y in ys:
                Cy = sum((s[t] * (2 * (y[t] == 0) + (y[t] == 2)) for t in range(k)), zero)
                mon = one
                for t in range(k):
                    if y[t] == 1:
                        mon = mon * s[t]
                Cp = Cy * mon / two
                # the form f = vec_y(x) + sum_{t: y_t = 1} s_t zeta_{y - 1_t}
                f = {}
                for xi, cc in vec.get(y, {}).items():
                    f[xi] = f.get(xi, zero) + cc
                for t in range(k):
                    if y[t] == 1:
                        u = y[:t] + (0,) + y[t + 1:]
                        zi = var((tag, 'zeta', u))
                        f[zi] = f.get(zi, zero) + s[t]
                w = scale / Cp
                for i, ci in f.items():
                    for j, cj in f.items():
                        add(i, j, w * ci * cj)
            # K^{-1} on zeta
            for name, i in list(idx.items()):
                if name[0] == tag and name[1] == 'zeta':
                    u = name[2]
                    if wcol(u)[1:] != key[1:] or wcol(u)[0] != key[0] - 1:
                        continue
                    mon = one
                    for t in range(k):
                        if u[t] == 1:
                            mon = mon * s[t]
                    Kinv = 2 * two / mon
                    add(i, i, scale * Kinv)

    for a in rows:
        add_schur_pair(('z', a), V(a, a), one)
    for ia, a in enumerate(rows):
        for b in rows[ia + 1:]:
            Vab, Vba = V(a, b), V(b, a)
            Sv, Dv = {}, {}
            for z in set(Vab) | set(Vba):
                Sv[z], Dv[z] = {}, {}
                for xi in set(Vab.get(z, {})) | set(Vba.get(z, {})):
                    p, q = Vab.get(z, {}).get(xi, zero), Vba.get(z, {}).get(xi, zero)
                    Sv[z][xi] = p + q
                    Dv[z][xi] = p - q
            half = one / 2
            add_schur_pair(('zp', a, b), Sv, half)
            add_woodbury(('zm', a, b), Dv, half)
    return idx, L, Dp


def dense(idx, L, zero):
    n = len(idx)
    return [[L.get(i, {}).get(j, zero) for j in range(n)] for i in range(n)]


# (1) exact check at a rational point
w = [flint.fmpq(p, q) for p, q in [(2, 3), (5, 4), (3, 7), (7, 5)][:k]]
sp = [x * x for x in w]
Bp = Builder(point_ring(sp))
Qz, Qo = flint.fmpq(0), flint.fmpq(1)
idx, L, Dp = build_lift(Bp, sp, Qz, Qo)
n = len(idx)
xs = [idx[('x', st)] for st in states]
zs = [i for i in range(n) if i not in xs]
Mfull = dense(idx, L, Qz)
Mq = flint.fmpq_mat(n, n, [v for r in Mfull for v in r])
A = flint.fmpq_mat(len(xs), len(xs), [Mfull[i][j] for i in xs for j in xs])
W = flint.fmpq_mat(len(zs), len(xs), [Mfull[i][j] for i in zs for j in xs])
Zm = flint.fmpq_mat(len(zs), len(zs), [Mfull[i][j] for i in zs for j in zs])
S = A - W.transpose() * Zm.solve(W)
Y = [[Bp.Y(X1, X2) for X2 in states] for X1 in states]
ok1 = all(S[i, j] == Y[i][j] for i in range(len(xs)) for j in range(len(xs)))
# D' split check
okD = True
for key, (ys, Bm) in Dp.items():
    c = key[0] + sum(key[1:])
    for p, y in enumerate(ys):
        for q, yp in enumerate(ys):
            pred = flint.fmpq(0)
            if y == yp:
                Cy = sum((sp[t] * (2 * (y[t] == 0) + (y[t] == 2)) for t in range(k)), flint.fmpq(0))
                mon = flint.fmpq(1)
                for t in range(k):
                    if y[t] == 1:
                        mon *= sp[t]
                pred += Cy * mon / 2 ** (k - c)
            for t in range(k):
                for tp in range(k):
                    if y[t] == 1 and yp[tp] == 1:
                        u = y[:t] + (0,) + y[t + 1:]
                        if yp[:tp] + (0,) + yp[tp + 1:] == u:
                            mon = flint.fmpq(1)
                            for r in range(k):
                                if u[r] == 1:
                                    mon *= sp[r]
                            pred += mon * sp[t] * sp[tp] / 2 ** (k - c + 1)
            okD &= pred == Bm[p][q]
print(f'{bname}: lift size {n} ({len(xs)} block coordinates, {len(zs)} lifted); D\' = C\' + H K H^T exact: {okD}; '
      f'Schur complement of the lift = block (exact, rational point): {ok1}  [{time.time() - t0:.1f}s]', flush=True)
# eigen-check: the lift must be PSD at the point
pd = True
Mt = [r[:] for r in Mfull]
act = list(range(n))
while act:
    cand = [i for i in act if Mt[i][i] > 0]
    if not cand:
        break
    i = cand[0]
    rest = [r for r in act if r != i]
    for r in rest:
        if Mt[r][i] != 0:
            f = Mt[r][i] / Mt[i][i]
            for l in rest:
                Mt[r][l] -= f * Mt[i][l]
    act = rest
pd = all(Mt[i][i] == 0 for i in act) and all(Mt[i][j] == 0 for i in act for j in act)
print(f'   lift PSD at the point: {pd}', flush=True)
# (2) symbolic
K = Field(list(blk['names']))
if k == 4:
    ssym = [K(1)] + list(K.gens)
    Rr = Ring(ssym, K(0), K(1))
else:
    ssym = list(K.gens)
    Rr = Ring(ssym, K(0), K(1))
Bs = Builder(Rr)
idx2, L2, _ = build_lift(Bs, ssym, K(0), K(1))
assert idx2 == idx
Ms = dense(idx2, L2, K(0))
Mg = [[(x.n, x.d) for x in row] for row in Ms]
ob = obstruction(Mg)
if ob is None:
    print('   symbolic lift: a diagonal entry is not subtraction-free (test not applicable)')
else:
    print(f'   symbolic lift: {len(ob[0])} obstructed ordered pairs ({len(ob[1])} with a width certificate)', flush=True)
    nm, nf = mono_psd(Mg)
    print(f'   monomial-PSD test of the lift (q = 1): {nf} non-PSD of {nm} coefficient matrices  [{time.time() - t0:.1f}s]')
    homog = False  # lift entries are not homogeneous of one degree: keep every exponent
    if ob[0] and not ob[1]:
        fk = farkas(Mg, homog)
        print(f'   monomial rescalings of the lift: {fk[0]} {fk[1] if fk[1] is not None else ""}', flush=True)
        if fk[0] != 'certificate':
            d2 = rescale_lp(Mg, homog)
            if isinstance(d2, str) or d2 is None:
                print(f'   LP rescaling: {d2}')
            else:
                Mr = apply_rescale(Mg, d2, homog)
                ob2 = obstruction(Mr)
                nm2, nf2 = mono_psd(Mr)
                print(f'   after the LP rescaling: {len(ob2[0]) if ob2 else "n/a"} obstructed pairs; monomial-PSD: '
                      f'{nf2} non-PSD of {nm2}  [{time.time() - t0:.1f}s]', flush=True)
    elif ob[1]:
        print('   width certificate: no monomial rescaling of the lift is a fixed-vector matrix (exact)')
print('PASS' if (ok1 and okD and pd) else 'FAIL')
