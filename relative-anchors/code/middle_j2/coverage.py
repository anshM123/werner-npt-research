"""Coverage of the statements of Theorems F1b and F2b by the certified regimes of run_all.sh (bookkeeping), for
5 <= k <= 60.  The regimes are uniform in the sizes they leave symbolic, so the same bookkeeping holds for every k.

python coverage.py M2 [M1=n]

Two symmetries reduce the regimes:
  (S1) exchanging the names of the two levels: divide every squared weight by q (the rescaled blocks are homogeneous of
       degree 1, so positivity and kernels do not change) and rename; q becomes 1/q, the level sizes are exchanged,
       a copy of T with the weight of level i gets the weight of level 1 - i, and a part on level i becomes the same
       part on level 1 - i;
  (S2) family 1 only: exchanging t and tau (a permutation of the states R <-> R').
Family 1 (Theorem F1b): every two-value configuration with U meeting both levels, except those with t and tau of
different weights and both U-levels larger than M1, is in the orbit under (S1), (S2) of a certified regime: (la, lb) =
(0, 0) at (1, 2+), (2+, 1), (2+, 2+); (0, 1) at (1, 2+) and (m, 2+), 2 <= m <= M1.  (U in one level: Theorem F1a.)
Family 2 (Theorem F2b): every two-value configuration with U in one level, or with min(u0, u1) <= M2, is, for every
part, in the orbit under (S1) of a certified regime (part on level 0; t-weight equal to that of level lt; inv and mix
for lt = 0 only; the inv part at (3+, 3+) is not certified and is replaced by (m, 3+), (3+, m), 3 <= m <= M2)."""
import sys

M = int(sys.argv[1])
opts = dict(a.split('=', 1) for a in sys.argv[2:] if '=' in a)
M1 = int(opts.get("M1", "1"))
ok = True


def inreg(spec, u):
    if spec.endswith('+'):
        return u >= int(spec[:-1])
    return u == int(spec)


# ---------------- family 1
F1_TWO = {(0, 0): [('1', '2+'), ('2+', '1'), ('2+', '2+')],
          (0, 1): [('1', '2+')] + [(str(m), '2+') for m in range(2, M1 + 1)]}


def f1_cert(la, lb, u0, u1):
    return any(inreg(r0, u0) and inreg(r1, u1) for r0, r1 in F1_TWO.get((la, lb), []))


n1 = 0
nskip1 = 0
for k in range(5, 61):
    u = k - 2
    for la in (0, 1):
        for lb in (0, 1):
            for u0 in range(1, u):
                u1 = u - u0
                orbit = set()
                for (x, y, v0, v1) in [(la, lb, u0, u1), (1 - la, 1 - lb, u1, u0)]:
                    orbit.add((x, y, v0, v1))
                    orbit.add((y, x, v0, v1))
                cov = any(f1_cert(*c) for c in orbit)
                if not cov and la != lb and min(u0, u1) > M1:
                    nskip1 += 1          # outside the statement of Theorem F1b
                    continue
                ok &= cov
                n1 += 1
                if not cov:
                    print('  family 1 not covered:', k, la, lb, u0, u1)
print(f'family 1: {n1} two-level configurations with 5 <= k <= 60 covered: {ok} ({nskip1} configurations with '
      f'T-weights 1, q and both U-levels > {M1} are outside the statement)')

# ---------------- family 2
REG = {('inv', 0): [('1', '3+'), ('3+', '1'), ('2', '2'), ('2', '3+'), ('3+', '2')] +
       [(str(m), '3+') for m in range(3, M + 1)] + [('3+', str(m)) for m in range(3, M + 1)],
       ('mix', 0): [('2', '2'), ('2', '3+'), ('3+', '2'), ('3+', '3+')]}
for lt in (0, 1):
    REG[('std', lt)] = [('2', '2'), ('2', '3+'), ('3+', '1'), ('3+', '2'), ('3+', '3+')]
    REG[('two', lt)] = [('4+', '1'), ('4+', '2'), ('4+', '3+')]
    REG[('anti', lt)] = [('3+', '1'), ('3+', '2'), ('3+', '3+')]
ONE = {'inv': '3+', 'std': '3+', 'two': '4+', 'anti': '3+'}
need = {'inv': lambda a, b: True, 'std': lambda a, b: a >= 2, 'two': lambda a, b: a >= 4,
        'anti': lambda a, b: a >= 3, 'mix': lambda a, b: a >= 2 and b >= 2}


def f2_cert(part, lt, v0, v1):
    if v1 == 0:
        return part in ONE and inreg(ONE[part], v0)
    if v0 == 0:
        return False
    return any(inreg(r0, v0) and inreg(r1, v1) for r0, r1 in REG.get((part, lt), []))


n2 = 0
ok2 = True
for k in range(5, 61):
    u = k - 1
    for u0 in range(0, u + 1):
        u1 = u - u0
        if u0 > 0 and u1 > 0 and min(u0, u1) > M:
            continue                     # outside the statement of Theorem F2b
        for lt in (0, 1):
            confs = [('inv', None), ('mix', None)] + [(p, pl) for p in ('std', 'two', 'anti') for pl in (0, 1)]
            for part, pl in confs:
                sizes = (u0, u1) if pl in (None, 0) else (u1, u0)
                if not need[part](*sizes):
                    continue
                n2 += 1
                cands = []
                for (q_lt, q_pl, w0, w1) in [(lt, pl, u0, u1), (1 - lt, None if pl is None else 1 - pl, u1, u0)]:
                    if q_pl in (None, 0):
                        cands.append((part, q_lt, w0, w1))
                cov = any(f2_cert(*c) for c in cands)
                ok2 &= cov
                if not cov:
                    print('  family 2 not covered:', k, (u0, u1), lt, part, pl)
print(f'family 2: {n2} (k, sizes, t-weight, part) cases with 5 <= k <= 60 and U in one level or min(u0, u1) <= {M}, '
      f'all covered: {ok2}')
ok &= ok2
print('PASS' if ok else 'FAIL')
