"""The symbolic reduced forms used by the certificates (cert.py, over Q(q, x0, x1) in flint), evaluated at rational
points, against exact projections of the letter-1 block built by j1lib.py (Lemma B).

For every certified regime (part, U0, U1) and several points (rational weight ratio q = w1^2/w0^2, integer level
sizes in the regime, k <= 7): the Schur complement S of the vacuum states on the P-vectors, its kernel vectors, and the
leading principal minors of the gauge-fixed S (in the order used by the certificate) are compared, exactly, with the
same objects computed from explicit test vectors on the j1lib block.  The weights in j1lib are w0 = 1 on level 0 and
w1 on level 1; also checked with the levels exchanged in j1lib (level-0 copies placed last), so both levels of the
concrete model carry the irrep.
Usage: python check_sym_points.py"""
import time
from fractions import Fraction as Fr
import flint
from cert import build
from check_reduce import explicit
from l1sym import det, keep_for

t0 = time.time()


def ev(f, vals):
    def pe(p):
        tot = Fr(0)
        for m, c in p.to_dict().items():
            c = flint.fmpq(c)
            term = Fr(int(c.p), int(c.q))
            for v, e in zip(vals, m):
                if e:
                    term *= v ** int(e)
            tot += term
        return tot
    return pe(f.n) / pe(f.d)


REGIMES = [('inv', '3+', '1'), ('inv', '2', '2'), ('inv', '3+', '2'), ('inv', '3+', '3+'),
           ('std', '3+', '1'), ('std', '2', '2'), ('std', '2', '3+'), ('std', '3+', '2'), ('std', '3+', '3+'),
           ('two', '4+', '1'), ('two', '4+', '2'), ('two', '4+', '3+'),
           ('anti', '3+', '1'), ('anti', '3+', '2'), ('anti', '3+', '3+'),
           ('mix', '2', '2'), ('mix', '2', '3+'), ('mix', '3+', '3+')]


def sizes(spec, extra):
    if spec.endswith('+'):
        o = int(spec[:-1])
        return [(o + x, x) for x in extra]
    return [(int(spec), None)]


ok = True
cnt = 0
ws = [Fr(1, 3), Fr(2), Fr(5, 7)]
for part, s0, s1 in REGIMES:
    K, M, red, gauge = build(part, s0, s1)
    S = red['S']
    n = len(S)
    keep = [i for i in range(n) if i not in gauge]
    names = K.ctx.names() if hasattr(K.ctx, 'names') else None
    for (u0, x0) in sizes(s0, [0, 1]):
        for (u1, x1) in sizes(s1, [0, 1]):
            if u0 + u1 > 7:
                continue
            for w1 in ws:
                q = w1 * w1
                vals = [q] + ([Fr(x0)] if x0 is not None else []) + ([Fr(x1)] if x1 is not None else [])
                Ssym = [[ev(S[i][j], vals) for j in range(n)] for i in range(n)]
                ex = explicit(part, Fr(1), w1, u0, u1, Fr(1))
                ex2 = explicit(part, Fr(1), w1, u0, u1, Fr(1), swap=True)
                same = Ssym == ex['S'] and Ssym == ex2['S']
                kerok = all(sum(ex['S'][i][j] * kv[j] for j in range(n)) == 0 for kv in red['ker'] for i in range(n))
                Sg = [[ex['S'][i][j] for j in keep] for i in keep]
                mins_ex = [det([r[:m] for r in Sg[:m]], Fr(1)) for m in range(1, len(keep) + 1)]
                Sgs = [[Ssym[i][j] for j in keep] for i in keep]
                mins_sym = [det([r[:m] for r in Sgs[:m]], Fr(1)) for m in range(1, len(keep) + 1)]
                pos = all(v > 0 for v in mins_ex)
                good = same and kerok and mins_ex == mins_sym and pos
                ok &= good
                cnt += 1
                print(f'{part:4s} U0={s0:3s} U1={s1:3s} (u0,u1)=({u0},{u1}) w1={str(w1):4s}: S symbolic = S from j1lib {same}, '
                      f'kernel {kerok}, minors positive {pos} [{time.time() - t0:.1f}s]', flush=True)
# larger sizes: the symbolic forms evaluated at x up to 9 against the same reduction computed in exact rational
# arithmetic at those integer sizes (k up to 21; a consistency check of the symbolic path at sizes where every
# pattern of the sums is realised)
from l1sym import Model, reduce_part
for part, s0, s1 in REGIMES:
    K, M, red, gauge = build(part, s0, s1)
    S = red['S']
    n = len(S)
    for x0 in ([3, 9] if s0.endswith('+') else [None]):
        for x1 in ([2, 7] if s1.endswith('+') else [None]):
            u0 = int(s0[:-1]) + x0 if x0 is not None else int(s0)
            u1 = int(s1[:-1]) + x1 if x1 is not None else int(s1)
            q = Fr(3, 11)
            vals = [q] + ([Fr(x0)] if x0 is not None else []) + ([Fr(x1)] if x1 is not None else [])
            Ssym = [[ev(S[i][j], vals) for j in range(n)] for i in range(n)]
            kp, kq, g = keep_for(part, u0, u1)
            Sfr = reduce_part(Model(Fr(1), [Fr(1), q], [Fr(u0), Fr(u1)]), part, kp, kq)['S']
            same = Ssym == Sfr
            ok &= same
            cnt += 1
            print(f'{part:4s} U0={s0:3s} U1={s1:3s} (u0,u1)=({u0},{u1}) q=3/11: symbolic = rational-arithmetic reduction '
                  f'{same} [{time.time() - t0:.1f}s]', flush=True)
print(f'{cnt} comparisons; time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
