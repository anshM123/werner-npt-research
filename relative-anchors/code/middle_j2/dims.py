"""Family 2: the numbers of test vectors (the bases of keep_for2 in f2sym.py) times the irrep dimensions add up to the
numbers of states, k(k - 1)/2 (R), k(k - 1) (P) and (k - 1)(k - 2)/2 (Q), for every size pair (u0, u1) with
u0 + u1 = k - 1 and 5 <= k <= 60 (one level included).  Parts on level 0 at (u0, u1); parts on level 1 are the parts on
level 0 at (u1, u0).  Usage: python dims.py"""
from f2sym import keep_for2

dimf = {'inv': lambda u0, u1: 1, 'std': lambda u0, u1: u0 - 1, 'two': lambda u0, u1: u0 * (u0 - 3) // 2,
        'anti': lambda u0, u1: (u0 - 1) * (u0 - 2) // 2, 'mix': lambda u0, u1: (u0 - 1) * (u1 - 1)}
need = {'inv': lambda a, b: True, 'std': lambda a, b: a >= 2, 'two': lambda a, b: a >= 4,
        'anti': lambda a, b: a >= 3, 'mix': lambda a, b: a >= 2 and b >= 2}
ok = True
n = 0
for k in range(5, 61):
    u = k - 1
    for u0 in range(0, u + 1):
        u1 = u - u0
        tot = {'R': 0, 'P': 0, 'Q': 0}
        for (a0, a1, sw) in ((u0, u1, False), (u1, u0, True)):
            for p in ('inv', 'std', 'two', 'anti', 'mix'):
                if p in ('inv', 'mix') and sw:
                    continue
                if not need[p](a0, a1):
                    continue
                kR, kP, kQ, g, ker = keep_for2(p, min(a0, 3), min(a1, 3))
                d = dimf[p](a0, a1)
                tot['R'] += d * len(kR)
                tot['P'] += d * len(kP)
                tot['Q'] += d * len(kQ)
        want = {'R': k * (k - 1) // 2, 'P': k * (k - 1), 'Q': (k - 1) * (k - 2) // 2}
        n += 1
        ok &= tot == want
        if tot != want:
            print('  dimension count wrong:', k, (u0, u1), tot, want)
print(f'family 2: dimension counts right for all {n} size pairs with 5 <= k <= 60: {ok}')
print('PASS' if ok else 'FAIL')
