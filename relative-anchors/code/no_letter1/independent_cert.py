"""A separate implementation (from the formulas of NOTES.md, Sections 3 to 6, not from certlib.py) of the two kernel
certificates of Theorem DJ, with python-flint polynomials.
(B) row inequality: P(t1..t4) with t_m = x_T + s_m: all coefficients >= 0 and a positive monomial in (p, q, e).
(A) pair inequality: numerator of Dk(t, t') with t = x_T + 2e s/(1+s): coefficient signs; AM-GM repair."""
import itertools, time
import flint
t0 = time.time()
R = flint.fmpq_mpoly_ctx.get(('p', 'q', 'r', 'e', 's1', 's2', 's3', 's4'), 'lex')
p, q, r, e, s1, s2, s3, s4 = R.gens()
xT = p + q + r
c = xT + 2 * e
al1, al2 = c - p, c - q
be1, be2 = c - p + q, c - q + p

def den(t):
    return t * (t + p) * (t + q) * (t + p + q)

def Xh(t):   # e den (2/e + 1/t)
    return 2 * den(t) + e * (t + p) * (t + q) * (t + p + q)

def Zh(t):   # e den al1 ((2 + p/al1)/e + 1/(t + p))
    return den(t) * (2 * al1 + p) + e * al1 * t * (t + q) * (t + p + q)

def Wh(t):
    return den(t) * (2 * al2 + q) + e * al2 * t * (t + p) * (t + p + q)

def Yh(t):   # e den be1 be2 ((2 + p/be1 + q/be2)/e + 1/(t + p + q))
    return den(t) * (2 * be1 * be2 + p * be2 + q * be1) + e * be1 * be2 * t * (t + p) * (t + q)

def Mh(t):
    return den(t)

ts = [xT + s for s in (s1, s2, s3, s4)]
cache = {}
def slot(f, m):
    key = (f.__name__, m)
    if key not in cache:
        cache[key] = f(ts[m])
    return cache[key]

P = R.from_dict({})
for perm in itertools.permutations(range(4)):
    a, b, cc, d = perm
    P += 2 * q * slot(Zh, a) * slot(Zh, b) * slot(Wh, cc) * slot(Yh, d)
    P -= q * c ** 2 * slot(Xh, a) * slot(Wh, b) * slot(Yh, cc) * slot(Mh, d)
    P -= (q + r) * c ** 2 * slot(Xh, a) * slot(Yh, b) * slot(Zh, cc) * slot(Mh, d)
    P += (q + r) * (be1 + p) * (be2 + q) * slot(Zh, a) * slot(Zh, b) * slot(Wh, cc) * slot(Mh, d)
coeffs = list(P.coeffs())
neg = sum(1 for v in coeffs if v < 0)
pure = any(m[2] == 0 and all(m[i] == 0 for i in range(4, 8)) and v > 0 for m, v in zip(P.monoms(), coeffs))
print(f'(B) P: {len(P)} terms, total degree {P.total_degree()}, degrees {P.degrees()}, negative coefficients {neg}, '
      f'positive monomial in (p, q, e) alone: {pure} [{time.time() - t0:.1f}s]', flush=True)
okB = neg == 0 and pure

# (A) the pair kernel, with t = xT + 2e s/(1 + s): multiply through by the denominators
R2 = flint.fmpq_mpoly_ctx.get(('p', 'q', 'r', 'e', 's1', 's2'), 'lex')
p, q, r, e, s1, s2 = R2.gens()
xT = p + q + r; c = xT + 2 * e
al1, al2 = c - p, c - q
be1, be2 = c - p + q, c - q + p
# per-atom functions as fractions N/D with t = (xT (1 + s) + 2 e s)/(1 + s)
def tnum(s):
    return xT * (1 + s) + 2 * e * s
# X(t) = 2/e + 1/t = (2t + e)/(e t);  with t = tn/(1+s): X = (2 tn + e (1+s))/(e tn)
# Z(t) = (2 + g_i)/e + 1/(t + p), g_i = p/al1: Z = ((2 al1 + p)(t + p) + e al1)/(e al1 (t + p))
# W(t) = ((2 al2 + q)(t + q) + e al2)/(e al2 (t + q))
# Y(t) = ((2 be1 be2 + p be2 + q be1)(t + p + q) + e be1 be2)/(e be1 be2 (t + p + q))
# delta = (1 + g_i)(1 + g_l) = c^2/(al1 al2), eps = (1 + G_i)(1 + G_l) = (be1 + p)(be2 + q)/(be1 be2)
def parts(s):
    tn = tnum(s); o = 1 + s
    X = (2 * tn + e * o, e * tn)                                   # numerator, denominator (common factor 1/o cancels)
    Z = ((2 * al1 + p) * (tn + p * o) + e * al1 * o, e * al1 * (tn + p * o))
    W = ((2 * al2 + q) * (tn + q * o) + e * al2 * o, e * al2 * (tn + q * o))
    Y = ((2 * be1 * be2 + p * be2 + q * be1) * (tn + (p + q) * o) + e * be1 * be2 * o, e * be1 * be2 * (tn + (p + q) * o))
    return X, Y, Z, W
X1, Y1, Z1, W1 = parts(s1); X2, Y2, Z2, W2 = parts(s2)
# Dk = delta (X1 Y2 + X2 Y1) - eps (Z1 W2 + Z2 W1); clear: delta = c^2/(al1 al2), eps = (be1+p)(be2+q)/(be1 be2)
# common denominator: al1 al2 be1 be2 * [all per-atom denominators]
Dall = X1[1] * Y1[1] * Z1[1] * W1[1] * X2[1] * Y2[1] * Z2[1] * W2[1]
def frac_term(A, B):   # A[0]/A[1] * B[0]/B[1] times Dall
    return A[0] * B[0] * (Dall // (A[1] * B[1])) if True else None
# exact division check
for A, B in ((X1, Y2), (X2, Y1), (Z1, W2), (Z2, W1)):
    qq, rr = divmod(Dall, A[1] * B[1])
    assert rr == 0
ND = (c ** 2 * be1 * be2) * (X1[0] * Y2[0] * (Dall // (X1[1] * Y2[1])) + X2[0] * Y1[0] * (Dall // (X2[1] * Y1[1]))) \
     - ((be1 + p) * (be2 + q) * al1 * al2) * (Z1[0] * W2[0] * (Dall // (Z1[1] * W2[1])) + Z2[0] * W1[0] * (Dall // (Z2[1] * W1[1])))
# ND / (al1 al2 be1 be2 Dall) = Dk; the per-atom denominators share factors, so remove the common content first
g = ND.gcd(Dall)
NDr = ND // g
cs = list(NDr.coeffs()); ms = list(NDr.monoms())
negs = [(m, v) for m, v in zip(ms, cs) if v < 0]
print(f'(A) numerator of Dk after removing gcd with the denominator: {len(NDr)} terms, total degree {NDr.total_degree()}, '
      f'negative coefficients: {len(negs)} {[(m, str(v)) for m, v in negs[:4]]} [{time.time() - t0:.1f}s]', flush=True)
okA = False
if len(negs) <= 1:
    if len(negs) == 0:
        rep = NDr
    else:
        (m, v) = negs[0]
        # expect v * p q e^a s1^2 s2^2; repair with -v/2 * p q e^a s1 s2 (s1 - s2)^2
        mon = R2.from_dict({m: 1})
        base = R2.from_dict({tuple(m[:4]) + (1, 1): 1})   # same (p, q, r, e) part, s1 s2
        rep = NDr - (-v / 2) * base * (s1 - s2) ** 2
    cs2 = list(rep.coeffs())
    pure = any(mm[2] == 0 and mm[4] == 0 and mm[5] == 0 and vv > 0 for mm, vv in zip(rep.monoms(), cs2))
    okA = all(vv >= 0 for vv in cs2) and pure
    print(f'    after the repair: negative coefficients {sum(1 for vv in cs2 if vv < 0)}, positive monomial in (p, q, e) '
          f'alone: {pure}', flush=True)
# the remaining denominator g * ... : check the factors of Dall/g and al, be are positive on the domain
Dr = (Dall // g) * al1 * al2 * be1 * be2
cst, fac = Dr.factor()
okD = cst > 0 and all(all(v > 0 for v in f.coeffs()) for f, _ in fac)
print(f'    denominator factors all positive coefficients: {okD}', flush=True)
print('PASS' if (okA and okB and okD) else 'FAIL')
