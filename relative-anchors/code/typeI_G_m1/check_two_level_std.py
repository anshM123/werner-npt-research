"""Two-level U-weights: the off-diagonal entries G_z[h,h'] with h, h' in the same level (the part of (Z) not covered by
check_two_level.py).  lam1 = G_z[h,h] - G_z[h,h'] (h != h' in U1) from two_level_std.py; Q3[1,1] = sum over U1 x U1
of G_z from two_level.py.  Then  -u1 G_z[h,h'] = u1 lam1 - Q3[1,1] ... more precisely G_z[h,h'] = (Q3[1,1]/u1 - lam1)/u1,
so (Z) on these entries <=> lam1 - Q3[1,1]/u1 > 0.
Validation: lam1 = exact elimination at rational points.  Certificates (positive coefficients after the shifts):
  (i)   u1 = 3 + x1, u2 = 2 + x2   (lam1 from case gen, Q3 from two_level gen);
  (ii)  u1 = 3 + x1, u2 = 1        (lam1 from case gen, Q3 from two_level one with the two levels exchanged);
  (iii) u1 = 2, u2 = 2 + x2        (lam1 from case u1is2, Q3 from two_level gen at u1 = 2).
(u1 = 2, u2 = 1 is u = 3, covered by check_m_small.py.)  The family is symmetric under exchanging the levels, so the
pairs inside U2 are covered by the same certificates with the levels relabelled."""
import sys, os, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from gzexact import Gz_exact

Q = flint.fmpq
ok = True


def load(fn, key):
    D = pickle.load(open(fn, 'rb'))
    ctx = flint.fmpq_mpoly_ctx.get(tuple(D['names']), 'lex')
    if key == 'lam1':
        n, d = D['lam1']
        return D['names'], ctx.from_dict(n), ctx.from_dict(d)
    G = D['Q3']
    return D['names'], [[(ctx.from_dict(G[i][j][0]), ctx.from_dict(G[i][j][1])) for j in range(3)] for i in range(3)]


TX = flint.fmpq_mpoly_ctx.get(('sigma', 's1', 's2', 'x1', 'x2'), 'lex')
g = TX.gens()
base = {'sigma': g[0], 's1': g[1], 's2': g[2], 'x1': g[3], 'x2': g[4]}


def transport(p, names, mp):
    """polynomial p (variables `names`) -> TX, substituting name -> mp[name] (a TX polynomial)."""
    out = TX.from_dict({})
    cache = {}
    for m, c in p.to_dict().items():
        term = TX.from_dict({(0,) * 5: c})
        for nm, e in zip(names, m):
            if e:
                key = (nm, e)
                if key not in cache:
                    cache[key] = mp[nm] ** e
                term = term * cache[key]
        out = out + term
    return out


def signs(p):
    c = list(p.coeffs())
    return 1 if all(x > 0 for x in c) else (-1 if all(x < 0 for x in c) else 0)


# ---- validation of lam1
t0 = time.time()
lam = {case: load(f'two_level_std_{case}.pkl', 'lam1') for case in ('gen', 'u1is2')}
nv = 0
for (a1, a2) in [(3, 1), (3, 2), (4, 1), (4, 3), (5, 2), (2, 1), (2, 2), (2, 4)]:
    case = 'u1is2' if a1 == 2 else 'gen'
    names, ln, ld = lam[case]
    for wt, om1, om2 in [(Q(3, 4), Q(1, 2), Q(5, 3)), (Q(2), Q(1, 30), Q(1, 3)), (Q(1, 7), Q(5, 2), Q(1, 100))]:
        G = Gz_exact(wt, [om1] * a1 + [om2] * a2)
        direct = G[1, 1] - G[1, 2]
        vals = {'sigma': wt * wt, 's1': om1 * om1, 's2': om2 * om2, 'u1': Q(a1), 'u2': Q(a2)}
        args = [vals[nm] for nm in names]
        good = ln(*args) / ld(*args) == direct
        ok &= good
        nv += 1
        if not good:
            print('   MISMATCH', a1, a2, wt, om1, om2)
print(f'validation: lam1 = exact elimination at {nv} points: {"all equal" if ok else "MISMATCH"} ({time.time() - t0:.1f} s)',
      flush=True)

# ---- certificates
qg_names, Qg = load('two_level_gen.pkl', 'Q3')
qo_names, Qo = load('two_level_one.pkl', 'Q3')
one = TX.from_dict({(0,) * 5: 1})
t0 = time.time()
regimes = [
    ('(i) u1 = 3 + x1, u2 = 2 + x2', 'gen', {'u1': base['x1'] + 3, 'u2': base['x2'] + 2},
     ('gen', {'u1': base['x1'] + 3, 'u2': base['x2'] + 2}), base['x1'] + 3),
    ('(ii) u1 = 3 + x1, u2 = 1', 'gen', {'u1': base['x1'] + 3, 'u2': one},
     ('one-swapped', {'u2': base['x1'] + 3}), base['x1'] + 3),
    ('(iii) u1 = 2, u2 = 2 + x2', 'u1is2', {'u2': base['x2'] + 2},
     ('gen', {'u1': 2 * one, 'u2': base['x2'] + 2}), 2 * one),
]
for name, lcase, lmap, (qcase, qmap), U1 in regimes:
    names, ln, ld = lam[lcase]
    mp = dict(base)
    mp.update(lmap)
    Ln, Ld = transport(ln, names, mp), transport(ld, names, mp)
    if qcase == 'gen':
        mq = dict(base)
        mq.update(qmap)
        qn, qd = Qg[1][1]
        Qn, Qd = transport(qn, qg_names, mq), transport(qd, qg_names, mq)
    else:
        # two_level 'one': its level 1 is the single copy; exchange the levels: s1 <-> s2, its u2 = our u1;
        # our sum over U1 x U1 is its Q3[2][2]
        mq = {'sigma': base['sigma'], 's1': base['s2'], 's2': base['s1']}
        mq.update(qmap)
        qn, qd = Qo[2][2]
        Qn, Qd = transport(qn, qo_names, mq), transport(qd, qo_names, mq)
    # lam1 - Q/U1 = (Ln Qd U1 - Qn Ld) / (Ld Qd U1)
    num = Ln * Qd * U1 - Qn * Ld
    den = Ld * Qd * U1
    sn, sd = signs(num), signs(den)
    good = sn * sd == 1
    ok &= good
    print(f'{name}: numerator {len(num.coeffs())} terms sign {sn:+d}, denominator {len(den.coeffs())} terms sign {sd:+d}: '
          f'{"G_z[h,h\'] < 0 certified" if good else "NOT certified"}', flush=True)
print(f'certificates done ({time.time() - t0:.1f} s)')
print('PASS' if ok else 'FAIL')
