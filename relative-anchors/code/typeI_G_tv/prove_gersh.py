"""Weighted block Gershgorin certificates (Lemma W of typeI_G_em/NOTES.md) at two U-weight levels, for T-weights in the
same two values.

Squared weights: s (level 1, u1 copies) and 1 (level 2, u2 copies; homogeneity).  The copies of T have squared weight
s or 1, sigma = s_T.  Frame per copy t (as in Theorem EM, two directions in the invariant part):
  rho_t;  y^1 = S1 z^1 + S2 z^2,  y^2 = z^1 - z^2  (z^i the invariant z-coordinate of level i; S_i = u_i s_i, S = S1 + S2),
  zhat^1 = kappa_1 y^1, kappa_1^2 = s_t sigma/(S (sigma + S))   (zeta_t = c zhat^1 w_U/|w_U| + ..., c^2 = 1 + S/sigma),
  zhat^2 = kappa_2 y^2, kappa_2^2 = s_t S1 S2/S                 (the unit invariant direction orthogonal to w_U),
  standard parts: zhat = w_t om_L |f| z (orthonormal).
Conditions (H = 1/2 by default), for every pair t != t' and every t:
  (A_inv)  H^2 s_t s_t'/sigma^2 - ||Ghat^inv_tt'||_F^2 > 0,  (A_std) H^2 s_t s_t'/sigma^2 - (Ghat^std_tt')^2 > 0,
  (B_inv)  Ghat^inv_tt - H (1 - theta_t) I > 0 (three LDL pivots),  (B_std) Ghat^std_tt - H (1 - theta_t) > 0.
By Lemma W these give G > 0 for every m >= 2 and every T whose weights lie in {s, 1}.  With T-weights in {s, 1} the pair
(s_t, s_t') is (s, s), (s, 1), (1, s) or (1, 1) (the norm of (1, s) is that of (s, 1)), and sigma = s_t + s_t' + r with
r >= 0 (A), sigma = s_t + r with r > 0 (B).
Usage: python prove_gersh.py one|gen [Hnum Hden]"""
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from flintfield import Field, F as FE

Q = flint.fmpq
CASE = sys.argv[1] if __name__ == '__main__' else None
H = Q(int(sys.argv[2]), int(sys.argv[3])) if (__name__ == '__main__' and len(sys.argv) > 3) else Q(1, 2)
COLLECT = None          # when a list, conditions are collected (name, value, free) instead of certified
# optional: only=<tag> runs the single group of conditions with that tag (A_ss, A_s1, A_11, B_s, B_1, std), one per
# process, to keep each process under the memory limit
ONLY = next((a[5:] for a in sys.argv if a.startswith('only=')), None) if __name__ == '__main__' else None
import gc
BZ = 'bz' in sys.argv          # B_inv pivots in the z-coordinates (rho, z^1, z^2) instead of (rho, y^1, y^2)


def want(tag):
    return ONLY is None or ONLY == tag


def load(fn):
    D = pickle.load(open(fn, 'rb'))
    names = D['names']
    ctx = flint.fmpq_mpoly_ctx.get(tuple(names), 'lex')
    P = lambda d: ctx.from_dict({m: Q(p, q) for m, (p, q) in d.items()})
    return names, D['keep'], {k: (P(n), P(d)) for k, (n, d) in D['G'].items()}


VN = ('x1', 'x2', 'r', 's')
KX = Field(list(VN))
cx = KX.ctx
GX = dict(zip(VN, cx.gens()))
ONE = cx.from_dict({(0,) * len(VN): 1})


def F(p):
    return FE(KX, p, ONE, reduce=False)


def comp(nd, names, mp):
    n, d = nd
    args = [mp[nm] for nm in names]
    return FE(KX, n.compose(*args, ctx=cx), d.compose(*args, ctx=cx))


def cert(name, val, free):
    if COLLECT is not None:
        COLLECT.append((name, val, free))
        return True
    t0 = time.time()
    num, den = val.n, val.d
    idx = [VN.index(v) for v in free]
    cn = list(num.coeffs())
    neg = sum(1 for c in cn if c < 0)
    pos = len(cn) - neg
    sgn_num = 1 if neg == 0 else (-1 if pos == 0 else 0)
    freem = any(all(m[i] == 0 for i in idx) for m in num.monoms())
    c, fl = den.factor()
    sgn = 1 if c > 0 else -1
    den_ok = True
    for f, e in fl:
        fc = list(f.coeffs())
        if not any(all(m[i] == 0 for i in idx) for m in f.monoms()):
            den_ok = False
        if all(z > 0 for z in fc):
            continue
        if all(z < 0 for z in fc):
            if e % 2:
                sgn = -sgn
            continue
        den_ok = False
    ok = (sgn_num * sgn == 1) and freem and den_ok and sgn_num != 0
    degs = [max((m[i] for m in num.monoms()), default=0) for i in range(len(VN))]
    print('    %-40s num %7d terms (+%d -%d), degrees %s, free monomial %s, den %d factors one-signed %s -> %s [%.1fs]'
          % (name, len(cn), pos, neg, degs, freem, len(fl), den_ok, 'CERTIFIED' if ok else 'NOT certified',
             time.time() - t0), flush=True)
    return ok


def orig(prefix, mode):
    """the closed form of one T-mode in its own variables (names), as rational functions FE over Field(names)."""
    names, keep, G = load('%s_%s.pkl' % (prefix, mode))
    K = Field(names)
    GK = {k: FE(K, n, d, reduce=False) for k, (n, d) in G.items()}
    gen = dict(zip(names, K.gens))
    return names, K, GK, gen


def late(val, names, mp):
    """compose the final rational function (over Field(names)) into the orthant variables."""
    return comp((val.n, val.d), names, mp)


def level_sums(K, gen, u1fix, u2fix):
    one = K(1)
    u1 = gen['u1'] if 'u1' in gen else one * u1fix
    u2 = gen['u2'] if 'u2' in gen else one * u2fix
    S1 = u1 * gen['s1']
    S2 = u2
    return S1, S2


def run_inv(prefix, umap, free_u, u1fix=None, u2fix=None, split=False, bonly=False):
    """(A_inv), (B_inv) for the invariant part.  All algebra in the variables of the closed form (u1, u2, sig, s1);
    only the final rational function is composed into the orthant variables (u_i = shift + x_i, sigma = s_t + s_t' + r
    or s_t + r, s1 = s)."""
    ok = True
    s, r = GX['s'], GX['r']
    wts_x = {'s': s, '1': ONE}
    # with split=True (used for u1, u2 >= 2) the Frobenius sum is certified in two parts, the entries in the rho row or
    # column and the 2 x 2 z-block, each against the budget H^2/2 (their sum is then < H^2 s_t s_t'/sigma^2)
    pieces = [('', [(i, j) for i in range(3) for j in range(3)], Q(1))] if not split else         [('_rho', [(i, j) for i in range(3) for j in range(3) if i == 0 or j == 0], Q(1, 2)),
         ('_zz', [(i, j) for i in range(1, 3) for j in range(1, 3)], Q(1, 2))]
    for mode, (ka, kb) in (() if bonly else (('T11', ('s', 's')), ('T12', ('s', '1')), ('T22', ('1', '1')))):
        for suffix, entries, budget in pieces:
            tag = 'A_' + ka + kb
            if not (want(tag) or want(tag + suffix)):
                continue
            t0 = time.time()
            names, K, GK, gen = orig(prefix, mode)
            one = K(1)
            wK = {'s': gen['s1'], '1': one}
            a, b, sig = wK[ka], wK[kb], gen['sig']
            S1, S2 = level_sums(K, gen, u1fix, u2fix)
            S = S1 + S2
            Nm = [[one, K(0), K(0)], [K(0), one / S, S2 / S], [K(0), one / S, -S1 / S]]
            Gz = [[GK[(i, 3 + j)] for j in range(3)] for i in range(3)]
            T = [[sum((Gz[k][l] * Nm[l][j] for l in range(3) if Nm[l][j] != 0), K(0)) for j in range(3)] for k in range(3)]
            Gy = {}
            for (i, j) in entries:
                Gy[(i, j)] = sum((Nm[k][i] * T[k][j] for k in range(3) if Nm[k][i] != 0), K(0))
            del T, Gz, GK
            gc.collect()
            kap = lambda t: [one, t * sig / (S * (sig + S)), t * S1 * S2 / S]
            ka_, kb_ = kap(a), kap(b)
            fro = K(0)
            for (i, j) in entries:
                fro = fro + Gy[(i, j)] * Gy[(i, j)] / (ka_[i] * kb_[j])
                Gy[(i, j)] = None
                gc.collect()
            val = one * H * H * budget * a * b / (sig * sig) - fro
            mp = dict(umap, sig=wts_x[ka] + wts_x[kb] + r, s1=s)
            del Gy, fro
            gc.collect()
            if COLLECT is None:
                print('    (A_inv (%s, %s)%s in the variables %s: numerator %d terms, denominator %d terms [%.1fs])'
                      % (ka, kb, suffix, names, len(val.n), len(val.d), time.time() - t0), flush=True)
            label = "%s A_inv (s_t, s_t') = (%s, %s)" % (prefix, ka, kb) + (' %s part, budget H^2/2' % suffix[1:] if suffix else '')
            ok &= cert(label, late(val, names, mp), tuple(free_u) + ('r',))
            del val
            gc.collect()
    for mode, ka in (('T11', 's'), ('T22', '1')):
        if not (want('B_' + ka) or want('B_' + ka + '12')):
            continue
        npiv = 2 if ONLY == 'B_' + ka + '12' else 3          # B_s12, B_112: only the first two pivots
        names, K, GK, gen = orig(prefix, mode)
        one = K(1)
        wK = {'s': gen['s1'], '1': one}
        a, sig = wK[ka], gen['sig']
        S1, S2 = level_sums(K, gen, u1fix, u2fix)
        S = S1 + S2
        Gz = [[GK[(min(i, j), max(i, j))] for j in range(3)] for i in range(3)]
        om_th = one * H * (sig - a) / sig
        if BZ:
            # the same condition in the z-coordinates (a congruence): G_z - H (1 - theta) W_z, W_z = diag(1, Omega),
            # Omega = (s_t/S) [sigma/(sigma + S) v v^T + S1 S2 e e^T], v = (S1, S2), e = (1, -1)
            c1 = a * sig / (S * (sig + S))
            c2 = a * S1 * S2 / S
            Om = [[c1 * S1 * S1 + c2, c1 * S1 * S2 - c2], [c1 * S1 * S2 - c2, c1 * S2 * S2 + c2]]
            W = [[one, K(0), K(0)], [K(0), Om[0][0], Om[0][1]], [K(0), Om[1][0], Om[1][1]]]
            A = [[Gz[i][j] - om_th * W[i][j] for j in range(3)] for i in range(3)]
        else:
            Nm = [[one, K(0), K(0)], [K(0), one / S, S2 / S], [K(0), one / S, -S1 / S]]
            T = [[sum((Gz[k][l] * Nm[l][j] for l in range(3) if Nm[l][j] != 0), K(0)) for j in range(3)] for k in range(3)]
            Gy = [[sum((Nm[k][i] * T[k][j] for k in range(3) if Nm[k][i] != 0), K(0)) for j in range(3)] for i in range(3)]
            kap = [one, a * sig / (S * (sig + S)), a * S1 * S2 / S]
            A = [[Gy[i][j] - (om_th * kap[i] if i == j else K(0)) for j in range(3)] for i in range(3)]
        mp = dict(umap, sig=wts_x[ka] + r, s1=s)
        for p in range(npiv):
            pv = A[p][p]
            ok &= cert('%s B_inv s_t = %s pivot %d%s' % (prefix, ka, p + 1, ' (z-coordinates)' if BZ else ''),
                       late(pv, names, mp), tuple(free_u) + ('r',))
            gc.collect()
            if p == npiv - 1:
                break           # no elimination step beyond the last pivot asked for
            for i in range(p + 1, npiv):
                f = A[i][p] / pv
                for j in range(p + 1, npiv):
                    A[i][j] = A[i][j] - f * A[p][j]
    return ok


def run_std(prefix, umap, free_u, level):
    """(A_std), (B_std) for the standard part of `level` (squared weight s1 for level 1, 1 for level 2)."""
    ok = True
    if not want('std'):
        return ok
    s, r = GX['s'], GX['r']
    wts_x = {'s': s, '1': ONE}
    for mode, (ka, kb) in (('T11', ('s', 's')), ('T12', ('s', '1')), ('T22', ('1', '1'))):
        names, K, GK, gen = orig(prefix, mode)
        one = K(1)
        wK = {'s': gen['s1'], '1': one}
        a, b, sig = wK[ka], wK[kb], gen['sig']
        sL = gen['s1'] if level == 1 else one
        g = GK[(0, 1)]
        val = one * H * H * a * b / (sig * sig) - g * g / (a * b * sL * sL)
        mp = dict(umap, sig=wts_x[ka] + wts_x[kb] + r, s1=s)
        ok &= cert("%s A_std (s_t, s_t') = (%s, %s)" % (prefix, ka, kb), late(val, names, mp), tuple(free_u) + ('r',))
    for mode, ka in (('T11', 's'), ('T22', '1')):
        names, K, GK, gen = orig(prefix, mode)
        one = K(1)
        wK = {'s': gen['s1'], '1': one}
        a, sig = wK[ka], gen['sig']
        sL = gen['s1'] if level == 1 else one
        g = GK[(0, 0)]
        val = g / (a * sL) - one * H * (sig - a) / sig
        mp = dict(umap, sig=wts_x[ka] + r, s1=s)
        ok &= cert('%s B_std s_t = %s' % (prefix, ka), late(val, names, mp), tuple(free_u) + ('r',))
    return ok


def conditions(case, h=Q(1, 2)):
    """the list of (name, rational function of (x1, x2, r, s), free variables) for the case, at Gershgorin constant h"""
    global COLLECT, H, CASE
    COLLECT, H, CASE = [], Q(h), case
    main(verbose=False)
    out, COLLECT = COLLECT, None
    return out


def main(verbose=True):
    t0 = time.time()
    x1, x2, s = GX['x1'], GX['x2'], GX['s']
    ok = True
    if CASE == 'one':
        u2 = 2 * ONE + x2
        verbose and print('part tl_inv_one (u1 = 1, u2 = 2 + x2)', flush=True)
        ok &= run_inv('tl_inv_one', dict(u2=u2), ('x2',), u1fix=1)
        verbose and print('part tl_std2_one (u2 = 3 + x2)', flush=True)
        ok &= run_std('tl_std2_one', dict(u2=3 * ONE + x2), ('x2',), 2)
        verbose and print('part tl_std2_onetwo (u2 = 2)', flush=True)
        ok &= run_std('tl_std2_onetwo', dict(), (), 2)
    elif CASE.startswith('f'):
        # (B_inv) only, with u1 = N fixed and u2 = 2 + x2 symbolic (the third pivot is out of memory reach with both
        # sizes symbolic); with the exchange of the levels this covers min(u1, u2) <= N
        n = int(CASE[1:])
        verbose and print('part tl_inv_f%d (u1 = %d, u2 = 2 + x2): (B_inv)' % (n, n), flush=True)
        ok &= run_inv('tl_inv_f%d' % n, dict(u2=2 * ONE + x2), ('x2',), u1fix=n, bonly=True)
    else:
        u1, u2 = 2 * ONE + x1, 2 * ONE + x2
        verbose and print('part tl_inv_gen (u1 = 2 + x1, u2 = 2 + x2)', flush=True)
        ok &= run_inv('tl_inv_gen', dict(u1=u1, u2=u2), ('x1', 'x2'), split=True)
        verbose and print('part tl_std1_gen (u1 = 3 + x1, u2 = 1 + x2)', flush=True)
        ok &= run_std('tl_std1_gen', dict(u1=3 * ONE + x1, u2=ONE + x2), ('x1', 'x2'), 1)
        verbose and print('part tl_std1_two (u1 = 2, u2 = 1 + x2)', flush=True)
        ok &= run_std('tl_std1_two', dict(u2=ONE + x2), ('x2',), 1)
    if verbose:
        print('H = %s: all certified: %s [%.1fs]' % (H, ok, time.time() - t0))
    return ok


if __name__ == '__main__':
    print('PASS' if main() else 'FAIL')
