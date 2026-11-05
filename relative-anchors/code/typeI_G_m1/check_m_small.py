"""Exact re-check of Conjecture M for u = 1, 2, 3 (m = 1) from the closed forms of closed_t.py (gauge x_t = 0, a
different gauge from the one of typeI_G_u3/closed_gen.py):
  (Z) every off-diagonal numerator of G_z has only negative coefficients, every diagonal numerator and the common
      denominator only positive ones;
  (D) every entry of den * G_z v has only positive coefficients, v = (4 s_t (s_t + S), 3 (2 s_t + S), ...).
Also: the closed form agrees exactly with gzexact.py (direct elimination over Q) at three rational points per u."""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from gzexact import Gz_exact

Q = flint.fmpq


def sg(p):
    c = list(p.to_dict().values())
    return '+' if c and all(x > 0 for x in c) else ('-' if c and all(x < 0 for x in c) else ('0' if not c else '?'))


ok = True
for u in (1, 2, 3):
    D = pickle.load(open(f'closed_t_u{u}.pkl', 'rb'))
    ctx = flint.fmpq_mpoly_ctx.get(tuple(D['names']), 'lex')
    G = D['G']
    n = len(G)
    nums = [[ctx.from_dict(G[i][j][0]) for j in range(n)] for i in range(n)]
    dens = [[ctx.from_dict(G[i][j][1]) for j in range(n)] for i in range(n)]
    # common denominator
    L = dens[0][0]
    for i in range(n):
        for j in range(n):
            g = L.gcd(dens[i][j])
            L = L * (dens[i][j] / g)
    Nt = [[nums[i][j] * (L / dens[i][j]) for j in range(n)] for i in range(n)]
    z = sg(L) == '+' and all(sg(Nt[i][j]) == ('+' if i == j else '-') for i in range(n) for j in range(n))
    s = ctx.gens()
    s0, S = s[0], sum(s[1:])
    v = [4 * s0 * (s0 + S)] + [3 * (2 * s0 + S)] * u
    rows = [sum((Nt[i][j] * v[j] for j in range(n)), 0 * L) for i in range(n)]
    dpos = all(sg(r) == '+' for r in rows)
    # agreement with direct elimination
    agree = True
    for (wt, wU) in [(Q(3, 4), [Q(1, 2), Q(5, 3), Q(2, 7)][:u]), (Q(2), [Q(1, 9), Q(4), Q(3)][:u]),
                     (Q(1, 5), [Q(7, 3), Q(1, 2), Q(1, 11)][:u])]:
        Gd = Gz_exact(wt, wU)
        vals = [wt * wt] + [x * x for x in wU]
        for i in range(n):
            for j in range(n):
                if nums[i][j](*vals) / dens[i][j](*vals) != Gd[i, j]:
                    agree = False
    ok &= z and dpos and agree
    print(f'u = {u}: (Z) pattern {z}; (D) every row of den*G_z v positive {dpos} '
          f'(terms {[len(r.to_dict()) for r in rows]}); denominator terms {len(L.to_dict())}; '
          f'closed form = direct elimination at 3 points {agree}', flush=True)
print('PASS' if ok else 'FAIL')
