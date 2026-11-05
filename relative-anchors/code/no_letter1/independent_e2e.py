"""(1) U empty: D and N as rational functions of (p, q, r) (c = x_T, Phi = 0); sign check.
(2) End to end, exact, a separate implementation of the closed form of Lemma 1: block PD by LDL, F_il < 0 and strict
diagonal dominance of F, at random rational configurations with |T| up to 12 and |U| up to 15."""
import sys, random, itertools
import flint
Q = flint.fmpq
# (1)
K = flint.fmpq_mpoly_ctx.get(('p', 'q', 'r'), 'lex')
p, q, r = K.gens()
xT = p + q + r; c = xT
al1, al2, be1, be2 = c - p, c - q, c - p + q, c - q + p
# work with numerators over the common denominator al1^2 al2 be1 be2
# X = 2, Y = 2 + p/be1 + q/be2, Z = 2 + p/al1, W = 2 + q/al2, delta = c^2/(al1 al2), eps = (be1+p)(be2+q)/(be1 be2)
Yn = 2 * be1 * be2 + p * be2 + q * be1          # Y = Yn/(be1 be2)
Zn = 2 * al1 + p; Wn = 2 * al2 + q              # Z = Zn/al1, W = Wn/al2
# D = delta X Y - eps Z W = c^2 * 2 * Yn/(al1 al2 be1 be2) - (be1+p)(be2+q) Zn Wn/(be1 be2 al1 al2)
Dn = 2 * c ** 2 * Yn - (be1 + p) * (be2 + q) * Zn * Wn
# N = 2qZ^2WY - q(1+g_i)^2 XWY - (q+r) delta X Y Z + (q+r) eps Z^2 W, times al1^2 al2 be1 be2:
Nn = 2 * q * Zn ** 2 * Wn * Yn - q * c ** 2 * 2 * Wn * Yn - (q + r) * c ** 2 * 2 * Yn * Zn + (q + r) * (be1 + p) * (be2 + q) * Zn ** 2 * Wn
for name, poly in (('D (U empty)', Dn), ('N (U empty)', Nn)):
    cs = list(poly.coeffs()); ms = list(poly.monoms())
    print(f'{name}: {len(poly)} terms, negative coefficients {sum(1 for v in cs if v < 0)}, '
          f'monomial in (p, q) alone: {any(m[2] == 0 and v > 0 for m, v in zip(ms, cs))}')
ok1 = all(all(v > 0 for v in P.coeffs()) for P in (Dn, Nn))
# (2)
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 99)
def closed(a, xU):
    n = len(a); xTv = sum(a, Q(0))
    Phi = lambda d: sum((x / (2 - 2 * x - xTv + d) for x in xU), Q(0))
    M = [[Q(0)] * (n + 1) for _ in range(n + 1)]
    for i in range(n):
        M[i][i] = Q(2)
    for i, l in itertools.combinations(range(n), 2):
        xS = xTv - a[i] - a[l]
        Gi = a[i] / (2 - 2 * a[i] - xS); Gl = a[l] / (2 - 2 * a[l] - xS)
        M[i][l] = M[l][i] = (1 + Gi) * (1 + Gl) / (2 + Gi + Gl + Phi(a[i] + a[l]))
    for i in range(n):
        gi = a[i] / (2 - a[i] - xTv)
        M[i][n] = M[n][i] = -2 * (1 + gi) / (2 + gi + Phi(a[i]))
    M[n][n] = 2 / (1 + Phi(Q(0)) / 2)
    return M
def ldl(Mx):
    A = [rw[:] for rw in Mx]; n = len(A)
    for cc in range(n):
        if A[cc][cc] <= 0:
            return False
        for i in range(cc + 1, n):
            f = A[i][cc] / A[cc][cc]
            for j in range(cc + 1, n):
                A[i][j] -= f * A[cc][j]
    return True
ok2 = True; cnt = 0
for trial in range(60):
    n = random.randint(2, 12); u = random.randint(0, 15)
    def rv():
        m = random.choice([1, 1, 1, 100, 10 ** 4, 10 ** 6])
        v = Q(random.randint(1, 50), random.randint(1, 50))
        return v * m if random.random() < 0.5 else v / m
    w = [rv() for _ in range(n + u)]
    S = sum(w, Q(0)); x = [v / S for v in w]
    a, xU = x[:n], x[n:]
    M = closed(a, xU)
    pd = ldl(M)
    # F = C - (Gamma/2) A A^T
    g = M[n][n]
    F = [[M[i][l] - M[i][n] * M[n][l] / g for l in range(n)] for i in range(n)]
    offneg = all(F[i][l] < 0 for i in range(n) for l in range(n) if i != l)
    dd = all(F[i][i] > sum((-F[i][l] for l in range(n) if l != i), Q(0)) for i in range(n))
    cnt += 1; ok2 &= pd and offneg and dd
    if not (pd and offneg and dd):
        print('FAILED at', n, u, [str(v) for v in w])
print(f'{cnt} random configurations (|T| <= 12, |U| <= 15): block PD, F_il < 0, strict diagonal dominance: {ok2}')
print('PASS' if (ok1 and ok2) else 'FAIL')
