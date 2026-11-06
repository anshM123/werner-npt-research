"""Family 2 (the |T| = 1 blocks of (R_inf)_2), part inv, regime (3+, 3+) (u0 = 3 + x0, u1 = 3 + x1), t-weight 0:
the Schur complement S' of qfirst_f2.py is positive definite on {q > 0, x0 >= 0, x1 >= 0}, by the determinant route.

The steps of qfirst_f2.py before its elimination are kept (Gram determinants, kernel vectors, R_QQ = 2 G_Q R_W^-1 G_Q,
R_W > 0); then:
 1. every entry denominator of S' is a product of factors one-signed in coefficients with a monomial free of x0, x1
    (cert.signcheck), so S' is continuous on the region;
 2. D*: for every irreducible factor f of an entry denominator, exponent max_p sum_a mult_f(den S'_{a,p(a)}); then
    T = D* det S' is a polynomial with deg_x1 T <= B = deg_x1 D* + max_p sum_a (deg_x1 num - deg_x1 den);
 3. Den_c, the candidate denominator of det S' (exponents observed on slices; Den_c | D* is checked), R = D*/Den_c;
 4. slices x1 = c, c = 0..B: det S' over Q(q, x0) exactly, times Den_c, must be a polynomial (saved);
 5. interpolation in x1 over c = 0..B: the coefficients above x1^(B - deg_x1 R) vanish, so T and N~ R (N~ the
    interpolant) have degree <= B in x1 and agree at B + 1 values: T = N~ R and det S' = N~/Den_c; every coefficient
    of N~ is >= 0 and one free of x0, x1 is > 0, so det S' > 0 on the region (Den_c > 0 by 1);
 6. S' is positive definite at (q, x0, x1) = (1, 0, 0) (exact elimination).
Since the region is connected and S' is continuous there with det S' > 0, S' > 0 on the whole region.
python f2det.py bounds | slices W NW | interp | final       (data in f2det_out/)"""
import sys, os, time, pickle
from fractions import Fraction as Fr
import flint
from flintfield import Field, F as FF
from f2sym import Model2, reduce_part2, assemble, part_spec2, keep_for2, pair_sums2
from l1sym import det, minv, mm
from cert import signcheck
import t2lev

OUT = 'f2det_out'


def maxperm(W):
    n = len(W)
    best = {0: 0}
    for a in range(n):
        nb = {}
        for mask, val in best.items():
            for b in range(n):
                if mask >> b & 1 or W[a][b] is None:
                    continue
                m2 = mask | (1 << b)
                v = val + W[a][b]
                if m2 not in nb or v > nb[m2]:
                    nb[m2] = v
        best = nb
    return best.get((1 << n) - 1)


def build(log=print):
    """S' exactly as in qfirst_f2.py (part inv, (3+, 3+), tw = 0, order PR) with its checks; returns (K, S, ok)"""
    t0 = time.time()
    part, o0, o1 = 'inv', 3, 3
    K = Field(['q', 'x0', 'x1'])
    q, x0, x1 = K.gens
    one, zero = K(1), K(0)
    M = Model2(one, one, [one, q], [x0 + o0, x1 + o1], Fr(1))
    red = reduce_part2(M, part, o0, o1)
    ok = True
    for gname in ('GR', 'GP', 'GQ'):
        if gname in red and len(red[gname]):
            sg, desc = signcheck(det(red[gname], one))
            ok &= sg == 1
    F, iR, iP, iQ = assemble(red)
    n = len(F)
    ok &= all(F[i][j] == F[j][i] for i in range(n) for j in range(n))
    for kv in red['ker']:
        vec = [zero] * n
        for j, c in enumerate(kv):
            vec[iP[j]] = K(c)
        ok &= all(sum((F[i][j] * vec[j] for j in range(n)), zero) == 0 for i in range(n))
    named, R, P, Qv, ker0, gauge0, names = part_spec2(part)
    kR, kP, kQ, gauge, ker = keep_for2(part, o0, o1)
    Qk = [Qv[i] for i in kQ]
    RW = pair_sums2(M, named, Qk, Qk, 'pair', 'pair', M.W4)
    GQ = red['GQ']
    nQ = len(Qk)
    GQi = minv(GQ, one)
    QQ2 = mm(mm(GQ, minv(RW, one), zero), GQ, zero)
    ok &= all(red['QQ'][i][j] == 2 * QQ2[i][j] for i in range(nQ) for j in range(nQ))
    res = []
    t2lev.leading_minors([row[:] for row in RW], one, lambda m_, d: res.append(signcheck(d)[0] == 1) or res[-1])
    ok &= all(res) and len(res) == nQ
    iPg = [iP[j] for j in range(len(iP)) if j not in red['gauge']]
    X = iPg + iR
    Wm = mm(mm(GQi, RW, zero), GQi, zero)
    FXQ = [[F[i][j] for j in iQ] for i in X]
    T = mm(FXQ, Wm, zero)
    S = [[F[X[a]][X[b]] - sum((T[a][l] * F[iQ[l]][X[b]] for l in range(nQ)), zero) / 2 for b in range(len(X))]
         for a in range(len(X))]
    ok &= all(S[a][b] == S[b][a] for a in range(len(X)) for b in range(len(X)))
    log(f"  S' built ({len(X)} x {len(X)}), checks of qfirst_f2.py: {ok} [{time.time()-t0:.1f}s]")
    return K, S, ok


def factors_of(K, S):
    facs, mult = {}, {}
    n = len(S)
    for a in range(n):
        for b in range(n):
            if S[a][b].n.is_zero():
                continue
            c, fl = S[a][b].d.factor()
            for f, m in fl:
                facs[str(f)] = f
                mult[(str(f), a, b)] = m
    return list(facs.values()), mult


def spec(K2, P, var, c):
    """substitute variable var (1 = x0, 2 = x1) = c; the remaining two become the generators of K2"""
    g = K2.ctx.gens()
    cst = K2.ctx.from_dict({(0, 0): c})
    args = [g[0], g[1], cst] if var == 2 else [g[0], cst, g[1]]
    return P.compose(*args, ctx=K2.ctx)


def det_slice(K2, S, var, c):
    n = len(S)
    A = [[(FF(K2, spec(K2, S[a][b].n, var, c), spec(K2, S[a][b].d, var, c)) if not S[a][b].n.is_zero() else K2(0))
          for b in range(n)] for a in range(n)]
    dt = K2(1)
    for m_ in range(n):
        p = A[m_][m_]
        dt = dt * p
        for i in range(m_ + 1, n):
            if A[i][m_] != 0:
                f = A[i][m_] / p
                for j in range(m_ + 1, n):
                    A[i][j] = A[i][j] - f * A[m_][j]
    return dt


def multiplicity(f, P):
    e = 0
    while True:
        qq, r = divmod(P, f)
        if r != 0:
            return e
        P = qq
        e += 1


def setup(log=print):
    K, S, ok = build(log)
    facs, mult = factors_of(K, S)
    n = len(S)
    estar = []
    for f in facs:
        W = [[(mult.get((str(f), a, b), 0) if not S[a][b].n.is_zero() else None) for b in range(n)] for a in range(n)]
        estar.append(maxperm(W))
    return K, S, ok, facs, estar


def bounds(K, S, facs, estar, dexp):
    n = len(S)
    out = {}
    for iv, v in enumerate(('q', 'x0', 'x1')):
        W = [[(S[a][b].n.degrees()[iv] - S[a][b].d.degrees()[iv]) if not S[a][b].n.is_zero() else None
              for b in range(n)] for a in range(n)]
        dD = sum(e * f.degrees()[iv] for f, e in zip(facs, estar))
        dR = dD - sum(e * f.degrees()[iv] for f, e in zip(facs, dexp))
        out[v] = (dD + maxperm(W), dR)
    return out


def main():
    mode = sys.argv[1]
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    if mode == 'bounds':
        K, S, ok, facs, estar = setup()
        n = len(S)
        dens_ok = all(signcheck(FF(K, S[a][b].d, K.ctx.from_dict({(0, 0, 0): 1})))[0] == 1
                      for a in range(n) for b in range(n) if not S[a][b].n.is_zero())
        print(f'  every entry denominator of S\' one-signed positive on the region: {dens_ok}', flush=True)
        # observed exponents on slices: (q, x0) at x1 = c and (q, x1) at x0 = c
        obs = [0] * len(facs)
        for var, cs in ((2, (0, 3, 7)), (1, (0, 5))):
            K2 = Field(['q', 'x1' if var == 1 else 'x0'])
            for c in cs:
                dt = det_slice(K2, S, var, c)
                row = []
                for k_, f in enumerate(facs):
                    fs = spec(K2, f, var, c)
                    if fs.is_constant():
                        row.append('-')
                        continue
                    m = multiplicity(fs, dt.d)
                    assert multiplicity(fs, dt.n) == 0
                    obs[k_] = max(obs[k_], m)
                    row.append(str(m))
                print(f'  slice {"x1" if var == 2 else "x0"} = {c}: observed exponents {" ".join(row)} '
                      f'[{time.time()-t0:.0f}s]', flush=True)
        print(f'  D* exponents {estar}; observed {obs}', flush=True)
        bd = bounds(K, S, facs, estar, obs)
        print(f'  degree bounds of T and degrees of R (q, x0, x1): {bd}', flush=True)
        with open(os.path.join(OUT, 'setup.pkl'), 'wb') as fh:
            pickle.dump({'facs': [str(f) for f in facs], 'estar': estar, 'dexp': obs, 'bounds': bd,
                         'dens_ok': dens_ok, 'ok': ok}, fh)
        print('PASS' if ok and dens_ok and all(o <= e for o, e in zip(obs, estar)) else 'FAIL')
    elif mode == 'slices':
        w, nw = int(sys.argv[2]), int(sys.argv[3])
        with open(os.path.join(OUT, 'setup.pkl'), 'rb') as fh:
            st = pickle.load(fh)
        K, S, ok, facs, estar = setup(log=lambda m: None)
        assert [str(f) for f in facs] == st['facs'] and estar == st['estar']
        Dc = K.ctx.from_dict({(0, 0, 0): 1})
        for f, e in zip(facs, st['dexp']):
            Dc = Dc * f ** e
        B = st['bounds']['x1'][0]
        K2 = Field(['q', 'x0'])
        good = True
        for c in range(B + 1):
            if c % nw != w:
                continue
            fn = os.path.join(OUT, f's_{c}.pkl')
            if os.path.exists(fn):
                continue
            t1 = time.time()
            dt = det_slice(K2, S, 2, c)
            Dc2 = spec(K2, Dc, 2, c)
            qq, r = divmod(Dc2, dt.d)
            if r != 0:
                print(f'x1={c}: det * Den_c NOT A POLYNOMIAL', flush=True)
                good = False
                continue
            P = dt.n * qq
            data = {tuple(e): (int(v.p), int(v.q)) for e, v in P.to_dict().items()}
            with open(fn + '.tmp', 'wb') as fh:
                pickle.dump(data, fh, protocol=4)
            os.replace(fn + '.tmp', fn)
            print(f'x1={c}: {len(data)} terms saved [{time.time()-t1:.0f}s, total {time.time()-t0:.0f}s]', flush=True)
        print('PASS' if good else 'FAIL')
    elif mode == 'interp':
        with open(os.path.join(OUT, 'setup.pkl'), 'rb') as fh:
            st = pickle.load(fh)
        B, dR = st['bounds']['x1']
        L = B - dR
        vals = []
        mons = set()
        for c in range(B + 1):
            with open(os.path.join(OUT, f's_{c}.pkl'), 'rb') as fh:
                d = pickle.load(fh)
            vals.append(d)
            mons.update(d.keys())
        mons = sorted(mons)
        V = flint.fmpq_mat(B + 1, B + 1, [flint.fmpq(c) ** l for c in range(B + 1) for l in range(B + 1)]).inv()
        Y = flint.fmpq_mat(B + 1, len(mons))
        for c, d in enumerate(vals):
            for j, mo in enumerate(mons):
                if mo in d:
                    Y[c, j] = flint.fmpq(*d[mo])
        Cf = V * Y
        nbad = nneg = npos = 0
        base = False
        maxl = 0
        for j, mo in enumerate(mons):
            for l in range(B + 1):
                v = Cf[l, j]
                if v == 0:
                    continue
                if l > L:
                    nbad += 1
                elif v < 0:
                    nneg += 1
                else:
                    npos += 1
                    maxl = max(maxl, l)
                    if mo[1] == 0 and l == 0:
                        base = True
        print(f'N~(q, x0, x1): {len(mons)} monomials q^a x0^b, interpolated over x1 = 0..{B}; nonzero coefficients '
              f'above x1^{L}: {nbad}; {npos} positive, {nneg} negative coefficients; degree in x1 {maxl}; a positive '
              f'coefficient free of x0, x1: {base} [{time.time()-t0:.0f}s]')
        print('PASS' if nbad == 0 and nneg == 0 and base else 'FAIL')
    elif mode == 'final':
        K, S, ok, facs, estar = setup(log=lambda m: None)
        pt = [flint.fmpq(1), flint.fmpq(0), flint.fmpq(0)]
        A = [[(S[a][b].n(*pt) / S[a][b].d(*pt)) if not S[a][b].n.is_zero() else flint.fmpq(0) for b in range(len(S))]
             for a in range(len(S))]
        n = len(A)
        pd = True
        for c in range(n):
            p = A[c][c]
            if p <= 0:
                pd = False
                break
            for r in range(c + 1, n):
                f = A[r][c] / p
                for cc in range(c + 1, n):
                    A[r][cc] -= f * A[c][cc]
        print(f"  S' positive definite at (q, x0, x1) = (1, 0, 0): {pd}")
        print('PASS' if pd and ok else 'FAIL')


if __name__ == '__main__':
    main()
