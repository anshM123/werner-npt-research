"""4-kernel of the homogenised row numerator N_H: build it with exact polynomials and test coefficient signs."""
import sys, time, itertools
import flint
t0 = time.time()
ctx = flint.fmpq_mpoly_ctx.get(('p', 'q', 'r', 'e', 's1', 's2', 's3', 's4'), 'lex')
p, q, r, e, s1, s2, s3, s4 = ctx.gens()
S = [s1, s2, s3, s4]
xT = p + q + r; c = xT + 2 * e
al1 = q + r + 2 * e; al2 = p + r + 2 * e; be1 = 2 * q + r + 2 * e; be2 = 2 * p + r + 2 * e
def slot(s):
    t = xT + s
    den = t * (t + p) * (t + q) * (t + p + q)
    Xh = 2 * den + e * (t + p) * (t + q) * (t + p + q)
    Zh = (2 * al1 + p) * den + e * al1 * t * (t + q) * (t + p + q)
    Wh = (2 * al2 + q) * den + e * al2 * t * (t + p) * (t + p + q)
    Yh = (2 * be1 * be2 + p * be2 + q * be1) * den + e * be1 * be2 * t * (t + p) * (t + q)
    Mh = den
    return {'X': Xh, 'Y': Yh, 'Z': Zh, 'W': Wh, 'M': Mh}
slots = [slot(s) for s in S]
print('slots built', [len(slots[0][k]) for k in 'XYZWM'], f'[{time.time()-t0:.1f}s]', flush=True)
words = [(2 * q, 'ZZWY'), (-q * c * c, 'XWYM'), (-(q + r) * c * c, 'XYZM'), ((q + r) * (be1 + p) * (be2 + q), 'ZZWM')]
P = ctx.from_dict({})
for coef, w in words:
    acc = ctx.from_dict({})
    for sg in itertools.permutations(range(4)):
        term = slots[sg[0]][w[0]] * slots[sg[1]][w[1]]
        term = term * slots[sg[2]][w[2]]
        term = term * slots[sg[3]][w[3]]
        acc = acc + term
    P = P + coef * acc
    print('word', w, 'done, P terms', len(P), f'[{time.time()-t0:.1f}s]', flush=True)
neg = [(ex, cf) for ex, cf in P.to_dict().items() if cf < 0]
print('4-kernel numerator: terms', len(P), 'total degree', P.total_degree(), 'negative coefficients', len(neg), f'[{time.time()-t0:.1f}s]')
print(neg[:20])
