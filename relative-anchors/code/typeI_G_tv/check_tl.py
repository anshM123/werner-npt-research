"""Exact validation of the two-level closed forms of tl_m.py against direct exact elimination of the whole type I block
(cf.TypeI over Q; gauge: the coordinates P_{r0,h} deleted, r0 the first copy of U; second gauge r1 at some points).

U = U1 (u1 copies of amplitude om1) + U2 (u2 copies of amplitude 1, so s2 = 1), T = m copies with amplitudes wT (all 1
in the eqT forms, so that sigma = m).  The exact G_z (zeta_{t,h} = w_t w_h z_{t,h}) is projected on
  inv block (t, t'):  3 x 3, [G_z[rho_t, rho_t'], sum_{h in U_j} G_z[rho_t, z_t'h], ...; sum_{h in U_i, h' in U_j} G_z[z_th, z_t'h']]
  std_L block (t, t'): G_z[z_th, z_t'h] - G_z[z_th, z_t'h']  (h != h' in U_L),
and compared with the closed forms at (u1, u2, sigma, s1, a = s_t, b = s_t') (diagonal blocks: a = s_t).  Also: the exact
G_z is invariant under a transposition inside each level (so the projections are complete), the diagonal closed forms
do not involve b.
Usage: python check_tl.py tl_inv_one.pkl [more pkl files]"""
import sys, os, random, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from emlib import Gz

Q = flint.fmpq


def load(fn):
    D = pickle.load(open(fn, 'rb'))
    names = D['names']
    ctx = flint.fmpq_mpoly_ctx.get(tuple(names), 'lex')
    P = lambda d: ctx.from_dict({m: Q(p, q) for m, (p, q) in d.items()})
    G = {k: (P(n), P(d)) for k, (n, d) in D['G'].items()}
    return D, names, G


def evalf(nd, vals):
    n, d = nd
    return n(*vals) / d(*vals)


def main(fn, rng, npts):
    D, names, G = load(fn)
    kind, case = D['kind'], D['case']
    mode = next((x for x in ('eqT', 'T11', 'T12', 'T22') if fn.endswith('_%s.pkl' % x)), None)
    eqT = mode is not None
    nk = len(D['keep'])
    half = nk // 2
    ok = True
    # diagonal blocks must not involve b
    if 'b' in names:
        ib = names.index('b')
        for (i, j), (n, d) in G.items():
            if j < half and any(m[ib] for p in (n, d) for m in p.monoms()):
                ok = False
                print('diagonal entry depends on b', i, j)
    u1c, u2c = {'inv': dict([('one', (1, None)), ('gen', (None, None))] + [('f%d' % n, (n, None)) for n in range(2, 61)]),
                'std1': {'gen': (None, None), 'two': (2, None)},
                'std2': {'one': (1, None), 'onetwo': (1, 2), 'gen': (None, None), 'two': (None, 2)}}[kind][case]
    amps = [Q(1, 7), Q(1, 3), Q(1, 2), Q(2, 3), Q(3, 2), Q(2), Q(5, 2), Q(4), Q(1, 20), Q(10), Q(3, 4)]
    t0 = time.time()
    nb = 0
    for it in range(npts):
        u1 = u1c if u1c is not None else rng.choice([2, 3, 4] if kind != 'std1' else [3, 4, 5])
        if kind == 'std1' and case == 'gen' and u1 < 3:
            u1 = 3
        u2 = u2c if u2c is not None else rng.choice([1, 2, 3, 4] if kind == 'std1' else ([3, 4, 5] if kind == 'std2' else [2, 3, 4]))
        if kind == 'inv' and u2 < 2:
            u2 = 2
        if u1 + u2 > 7:
            u2 = max(1 if kind != 'inv' else 2, 7 - u1)
        if case.startswith('f'):
            u2 = rng.choice([2, 3])
        m = rng.choice([2, 2, 3]) if not case.startswith('f') else 2
        om1 = rng.choice(amps)
        if mode in ('eqT', 'T22'):
            wT = [Q(1)] * m
        elif mode == 'T11':
            wT = [om1] * m
        elif mode == 'T12':
            wT = [om1, Q(1)] + [rng.choice([om1, Q(1), rng.choice(amps)]) for _ in range(m - 2)]
        else:
            wT = [rng.choice(amps) for _ in range(m)]
        ab = {None: None, 'eqT': (Q(1), Q(1)), 'T22': (Q(1), Q(1)), 'T11': (om1 * om1, om1 * om1),
              'T12': (om1 * om1, Q(1))}[mode]
        wU = [om1] * u1 + [Q(1)] * u2
        Gq, mm, uu = Gz(wT, wU)
        d = 1 + u1 + u2
        U1 = list(range(1, 1 + u1))
        U2 = list(range(1 + u1, 1 + u1 + u2))
        good = True
        # invariance under a transposition inside each level
        for lvl in (U1, U2):
            if len(lvl) >= 2:
                perm = list(range(d))
                perm[lvl[0]], perm[lvl[1]] = perm[lvl[1]], perm[lvl[0]]
                Pm = [t * d + perm[i] for t in range(m) for i in range(d)]
                if any(Gq[i, j] != Gq[Pm[i], Pm[j]] for i in range(m * d) for j in range(m * d)):
                    good = False
                    print('  not invariant')
        sT = [x * x for x in wT]
        sig = sum(sT, Q(0))
        s1 = om1 * om1
        for t in range(m):
            for tp in range(m):
                if ab is not None:
                    if t == tp and sT[t] != ab[0]:
                        continue
                    if t != tp and (sT[t], sT[tp]) != ab:
                        continue
                b0, b1 = t * d, tp * d
                if kind == 'inv':
                    idx = [[0], U1, U2]
                    direct = [[sum((Gq[b0 + i, b1 + j] for i in idx[p] for j in idx[q]), Q(0)) for q in range(3)]
                              for p in range(3)]
                else:
                    Lv = U1 if kind == 'std1' else U2
                    direct = [[Gq[b0 + Lv[0], b1 + Lv[0]] - Gq[b0 + Lv[0], b1 + Lv[1]]]]
                vals = {'u1': Q(u1), 'u2': Q(u2), 'sig': sig, 's1': s1, 'a': sT[t], 'b': sT[tp]}
                if t == tp:
                    vals['b'] = Q(7, 3)          # must not matter
                args = [vals[nm] for nm in names]
                for p in range(half):
                    for q in range(half):
                        if t == tp:
                            key = (min(p, q), max(p, q))
                        else:
                            key = (p, half + q)
                        val = evalf(G[key], args)
                        if val != direct[p][q]:
                            good = False
                            print('  MISMATCH', kind, case, u1, u2, m, t, tp, p, q)
                nb += 1
        if it % 3 == 0 and u1 + u2 >= 2:
            G2, _, _ = Gz(wT, wU, gauge_r0=m + u1 + u2 - 1)
            if G2 != Gq:
                good = False
                print('  gauge dependence')
        ok &= good
        print('%s: u1=%d u2=%d m=%d wT=%s om1=%s: %s [%.1fs]' % (fn, u1, u2, m, [str(x) for x in wT], om1,
                                                               'equal' if good else 'MISMATCH', time.time() - t0), flush=True)
    print('%s: %d configurations, %d (t, t\') blocks compared: %s' % (fn, npts, nb, 'all equal' if ok else 'MISMATCH'))
    return ok


if __name__ == '__main__':
    rng = random.Random(20261106)
    ok = True
    for fn in sys.argv[1:]:
        ok &= main(fn, rng, 4 if '_f' in fn else 8)
    print('PASS' if ok else 'FAIL')
