"""Bookkeeping for the proof: every isotypic part of the letter-1 block of (R_inf)_1 at every pair of level sizes
(u0, u1) with u0, u1 >= 0 and k = u0 + u1 >= 4 is covered by one of the certified regimes of run_all.sh.

A part is (irrep type, level of the irrep): 'inv'; 'std', 'two', 'anti' on level 0 or on level 1; 'mix'.  A part on
level 1 at (u0, u1) is the part on level 0 at (u1, u0) with the two weights exchanged (the certificates hold for every
weight ratio q > 0).  A regime (part, U0, U1) covers sizes (u0, u1) if each U is either the integer u or 'o+' with
u >= o.  Also checked: the parts exist exactly when the irreps occur (std needs u >= 2, anti u >= 3, two u >= 4, mix
u0, u1 >= 2), and the multiplicities times the irrep dimensions add up to k^2 and k(k-1)/2 for every k <= 60.
Usage: python coverage.py"""
from l1sym import keep_for

REGIMES = [('inv', '3+', '1'), ('inv', '2', '2'), ('inv', '3+', '2'), ('inv', '3+', '3+'), ('inv', '4+', '0'),
           ('std', '3+', '1'), ('std', '2', '2'), ('std', '2', '3+'), ('std', '3+', '2'), ('std', '3+', '3+'),
           ('std', '4+', '0'),
           ('two', '4+', '1'), ('two', '4+', '2'), ('two', '4+', '3+'), ('two', '4+', '0'),
           ('anti', '3+', '1'), ('anti', '3+', '2'), ('anti', '3+', '3+'), ('anti', '4+', '0'),
           ('mix', '2', '2'), ('mix', '2', '3+'), ('mix', '3+', '3+')]


def fits(spec, u):
    return u >= int(spec[:-1]) if spec.endswith('+') else u == int(spec)


def covered(part, a, c):
    """is the part with irrep on the level of size a (other level c) covered, directly or with levels exchanged
    (the latter only matters for 'inv' and 'mix', which involve both levels)?"""
    for p, s0, s1 in REGIMES:
        if p != part:
            continue
        if fits(s0, a) and fits(s1, c):
            return True
        if part in ('inv', 'mix') and fits(s0, c) and fits(s1, a):
            return True
    return False


def dims(u):
    return dict(std=u - 1 if u >= 2 else 0, two=u * (u - 3) // 2 if u >= 4 else 0,
                anti=(u - 1) * (u - 2) // 2 if u >= 3 else 0)


ok = True
n = 0
for k in range(4, 61):
    for u0 in range(0, k + 1):
        u1 = k - u0
        if u0 < u1:
            continue          # (u0, u1) and (u1, u0) are the same family (levels exchanged)
        parts = [('inv', u0, u1, 1)]
        for a, c in ((u0, u1), (u1, u0)):
            d = dims(a)
            for p in ('std', 'two', 'anti'):
                if d[p]:
                    parts.append((p, a, c, d[p]))
        if u0 >= 2 and u1 >= 2:
            parts.append(('mix', u0, u1, (u0 - 1) * (u1 - 1)))
        totP = totQ = 0
        for p, a, c, dim in parts:
            good = covered(p, a, c)
            if not good:
                print(f'not covered: k = {k}, (u0, u1) = ({u0}, {u1}), part {p} on a level of size {a}')
            ok &= good
            kp, kq, g = keep_for(p, a, c)
            totP += dim * len(kp)
            totQ += dim * (len(kq) if kq else 0)
            n += 1
        good = totP == k * k and totQ == k * (k - 1) // 2
        if not good:
            print(f'dimension count fails at (u0, u1) = ({u0}, {u1}): {totP}, {totQ}')
        ok &= good
print(f'{n} (size pair, part) cases for 4 <= k <= 60, all covered and all dimension counts right: {ok}')
print('PASS' if ok else 'FAIL')
