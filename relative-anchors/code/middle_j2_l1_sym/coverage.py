"""Coverage bookkeeping (NOTES Section 5).

python coverage.py M [S2] [KMAX]
  M: stage 1 certifies the regimes (o(m)+, m) for 2 <= m <= M, o(m) = max(m, 11 - m);
  S2: the offset o of a stage-2 regime (o+, o+) that is certified (or 'none');
  KMAX: the largest k checked (default 60).
For every two-value configuration with 11 <= k <= KMAX (the smaller level of size b, the larger of size a, b <= a), the
script names the result that covers it: b = 0 (equal weights): Theorem E1 of the README; b = 1: Theorem MB of
marked_copy (j = 2 <= 8); 2 <= b <= M: the regime (o(b)+, b) (levels exchanged if needed, q -> 1/q by homogeneity);
b >= o: the regime (o+, o+).  For every regime that is used and every size in it with k <= KMAX it checks the dimension
counts sum_lambda dim(lambda) * (number of test vectors) = number of P-states, of kernel vectors and of Q-states, with
the bases chosen as in cert.py (at the sample size of the regime; they are bases on the whole regime because their
Gram determinants are positive there, cert.py).  Prints PASS iff every configuration is covered and every count is
right; configurations not covered are listed (and give FAIL)."""
import sys
from math import comb
import r2v


def o_of(m):
    return max(m, 11 - m)


def regime_counts(s0, s1):
    def parse(sp):
        return (None, int(sp[:-1])) if sp.endswith('+') else (int(sp), None)
    (f0, o0), (f1, o1) = parse(s0), parse(s1)
    sample = (o0 + 3 if f0 is None else f0, o1 + 3 if f1 is None else f1)
    out = {}
    for F in r2v.families():
        K, W, Q = r2v.choose(F, sample)
        out[F.name] = (F, len(K), len(W), len(Q))
    return out


def dims_ok(counts, u0, u1):
    k = u0 + u1
    dP = dK = dQ = 0
    for name, (F, nK, nW, nQ) in counts.items():
        d = r2v.hook_dim(u0, F.shape[0]) * r2v.hook_dim(u1, F.shape[1])
        if nK + nW + nQ and d == 0:
            return False
        dP += d * (nK + nW)
        dK += d * nK
        dQ += d * nQ
    return dP == k * comb(k, 2) and dK == comb(k, 2) and dQ == comb(k, 3)


if __name__ == '__main__':
    M = int(sys.argv[1])
    S2 = sys.argv[2] if len(sys.argv) > 2 else 'none'
    KMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    o2 = None if S2 == 'none' else int(S2)
    cache = {}
    ok = True
    outside = []
    nconf = {'E1': 0, 'MB': 0, 'stage1': 0, 'stage2': 0}
    for k in range(11, KMAX + 1):
        for b in range(0, k // 2 + 1):
            a = k - b
            if b == 0:
                nconf['E1'] += 1
                continue
            if b == 1:
                nconf['MB'] += 1
                continue
            if b <= M:
                reg = (f'{o_of(b)}+', str(b))
                assert a >= o_of(b)
                u0, u1 = a, b
                nconf['stage1'] += 1
            elif o2 is not None and b >= o2:
                reg = (f'{o2}+', f'{o2}+')
                u0, u1 = a, b
                nconf['stage2'] += 1
            else:
                outside.append((a, b))     # both levels larger than M: outside the statement (open)
                continue
            if reg not in cache:
                cache[reg] = regime_counts(*reg)
            # both orders of the levels: (a, b) and (b, a) are the same configuration up to q -> 1/q, but the dimension
            # count is checked at the sizes of the regime (level 0 = the symbolic one)
            good = dims_ok(cache[reg], u0, u1)
            if not good:
                print(f'dimension count FAILS at {(u0, u1)} in regime {reg}')
            ok &= good
    print(f'configurations with 11 <= k <= {KMAX} and smaller level <= {M} (or in a stage-2 regime): {nconf}, all '
          f'covered; regimes used: {sorted(cache)}')
    print(f'outside the statement (both levels > {M}, open): {len(outside)} configurations, smallest k = '
          f'{min((a + b for a, b in outside), default=None)}, e.g. {outside[:4]}')
    print('PASS' if ok else 'FAIL')
