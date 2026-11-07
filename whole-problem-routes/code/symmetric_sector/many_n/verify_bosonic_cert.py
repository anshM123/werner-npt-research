"""[Adapted to arbitrary certificates: usage python verify_bosonic_cert.py n FILE.pkl [maxden]]
Independent exact check of the greedy certificate on the permutation-symmetric sector (the statement checked by
exact_greedy.py, by a different construction).

Statement checked, for a given n: with A = (I - F)/2, B = (I - Omega)/2 on each pair (A_s, B_s), C_{3/4} = I/4 + A/2 + B/2, the
operator Q_n = 4^{-n} sum_{k+s+t=n} multinom(n;k,s,t) g(s,t) (Pi (x) Pi)(I^k (x) A^s (x) B^t)(Pi (x) Pi),
g = 1 (s > t), 1/2 (s = t >= 1), 0 otherwise, is positive semidefinite on Sym^n(R^3) (x) Sym^n(R^3), with kernel of the size of
Sym^{2n}; and Q_n + Gamma(Q_n) = (Pi (x) Pi) C_{3/4}^{(x)n} (Pi (x) Pi) - 4^{-n} Pi (x) Pi. Hence for symmetric x, y:
<y|P_{3/4}^{(x)n}(x x^T)|y> - 4^{-n}|x|^2|y|^2 = 2 <x y|Q_n|x y> >= 0.

Method (shares no code with exact_greedy.py and its modules): coherent-state symbols.  For real u, w, u', w' the matrix element
<u^n (x) w^n| I^k A^s B^t |u'^n (x) w'^n> = X^k (alpha/2)^s (beta/2)^t, X = (u.u')(w.w'), alpha = X - (u.w')(w.u'),
beta = X - (u.w)(u'.w').  Since u^{(x)n} = sum_c u^c f_c (f_c = sum of e_word over words of content c), the coefficient of
u^c w^d u'^c' w'^d' in the symbol is the matrix element <f_c (x) f_d| O |f_c' (x) f_d'>.  The symbols are expanded exactly
(python-flint fmpz_mpoly).  O is SO(3)-invariant, so O >= 0 iff its restriction to ker J_z is >= 0; ker J_z is computed as
the exact null space of the rotation generator (in coefficient coordinates), and the restricted matrix is tested by
fraction-free symmetric elimination with diagonal pivoting.
usage: python verify_bosonic.py n"""
import sys, time, itertools, math
import flint

n = int(sys.argv[1])
T0 = time.time()
def log(*a):
    print(*a, flush=True)

names = ['u1', 'u2', 'u3', 'w1', 'w2', 'w3', 'p1', 'p2', 'p3', 'q1', 'q2', 'q3']   # p = u', q = w'
ctx = flint.fmpz_mpoly_ctx.get(tuple(names), 'lex')
V = ctx.gens()
u, w, p, q = V[0:3], V[3:6], V[6:9], V[9:12]
dot = lambda a, b: a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
X = dot(u, p) * dot(w, q)
alpha = X - dot(u, q) * dot(w, p)
beta = X - dot(u, w) * dot(p, q)

# identity check: (X + alpha + beta)^n - X^n = (T + T|_{w<->w'})/2 with T = sum multinom * 2g * X^k alpha^s beta^t
import pickle
from fractions import Fraction
from verify_exact import certificate_poly
d = pickle.load(open(sys.argv[2], 'rb'))
maxden = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
yfree = {p: Fraction(float(v)).limit_denominator(maxden) for p, v in zip(d['free'], d['y'])}
P = certificate_poly(n, yfree)          # Q (times 4^n) as {(a,b,c): coeff} in p = X, q = X - alpha, r = X - beta
den = 1
for v in P.values():
    den = math.lcm(den, Fraction(v).denominator)
# T := 2 * den * 4^n Q  as integer polynomial in X, alpha, beta (p = X, q = X - alpha, r = X - beta)
qv = X - alpha; rv = X - beta
T = ctx.from_dict({})
for (a, b, c), v in P.items():
    coef = int(Fraction(v) * den * 2)
    T += coef * (X ** a) * (qv ** b) * (rv ** c)
log(f'certificate {sys.argv[2]} with denominators lcm {den}; non-greedy pairs: ' + str({k: str(v) for k, v in yfree.items() if v != 1}))
log(f'n = {n}: symbol T = 2 4^n Q_n expanded, {len(T)} terms ({time.time() - T0:.1f}s)')
# w <-> w' substitution
Tsw = T.compose(*(list(u) + list(q) + list(p) + list(w)))
lhs = (X + alpha + beta) ** n - X ** n
assert 2 * den * lhs == T + Tsw, 'certificate identity fails'
log(f'identity (X + alpha + beta)^n - X^n = (T + T(w <-> w\'))/2 holds exactly ({time.time() - T0:.1f}s)')

# compositions of n into 3 parts
comps = [c for c in itertools.product(range(n + 1), repeat=3) if sum(c) == n]
cidx = {c: i for i, c in enumerate(comps)}
N = len(comps)
dim = N * N
# matrix of 2 4^n Q_n in the basis f_c (x) f_d (row: bra (u, w), column: ket (u', w'))
M = [[0] * dim for _ in range(dim)]
for exps, coeff in zip(T.monoms(), T.coeffs()):
    r = cidx[tuple(exps[0:3])] * N + cidx[tuple(exps[3:6])]
    c = cidx[tuple(exps[6:9])] * N + cidx[tuple(exps[9:12])]
    M[r][c] = int(coeff)
for r in range(dim):
    for c in range(r):
        assert M[r][c] == M[c][r], 'matrix not symmetric'
log(f'matrix of 2 4^n Q_n on Sym^n (x) Sym^n built: dimension {dim}, symmetric ({time.time() - T0:.1f}s)')

# ker J_z: vectors v = sum c_{cd} f_c (x) f_d whose polynomial sum c_{cd} m(c) m(d) u^c w^d is invariant under rotations about z
mult = [math.factorial(n) // (math.factorial(c[0]) * math.factorial(c[1]) * math.factorial(c[2])) for c in comps]
def deriv_rows():
    # generator L = u1 d/du2 - u2 d/du1 + w1 d/dw2 - w2 d/dw1 acting on monomials u^c w^d (polynomial coordinates)
    D = [[0] * dim for _ in range(dim)]
    for i, c in enumerate(comps):
        for j, d in enumerate(comps):
            col = i * N + j
            for (vec, which) in ((c, 0), (d, 1)):
                a1, a2, a3 = vec
                if a2 > 0:   # u1 d/du2: exponent (a1+1, a2-1, a3), factor a2
                    nv = (a1 + 1, a2 - 1, a3)
                    row = (cidx[nv] * N + j) if which == 0 else (i * N + cidx[nv])
                    D[row][col] += a2
                if a1 > 0:   # - u2 d/du1
                    nv = (a1 - 1, a2 + 1, a3)
                    row = (cidx[nv] * N + j) if which == 0 else (i * N + cidx[nv])
                    D[row][col] -= a1
    return D
Dm = flint.fmpz_mat(deriv_rows())
Kp, nullity = Dm.nullspace()      # columns: kernel in polynomial coordinates
log(f'ker J_z has dimension {nullity} ({time.time() - T0:.1f}s)')
# convert to f-coordinates: c_{cd} = pi_{cd} / (m(c) m(d)); scale each column to integers
cols = []
for kcol in range(nullity):
    vec = [flint.fmpq(int(Kp[r, kcol]), mult[r // N] * mult[r % N]) for r in range(dim)]
    den = 1
    for x in vec:
        den = math.lcm(den, int(x.q))
    cols.append([int(x * den) for x in vec])
Bm = flint.fmpz_mat([[cols[k][r] for k in range(nullity)] for r in range(dim)])
Mm = flint.fmpz_mat(M)
Kmat = Bm.transpose() * Mm * Bm
log(f'restriction to ker J_z: {nullity} x {nullity} ({time.time() - T0:.1f}s)')

def psd_bareiss(A):
    m = len(A); A = [row[:] for row in A]; act = list(range(m)); prev = 1
    while act:
        pv = max(act, key=lambda i: A[i][i]); d = A[pv][pv]
        if d < 0:
            return False, None
        if d == 0:
            for i in act:
                for j in act:
                    if A[i][j] != 0:
                        return False, None
            return True, len(act)
        act.remove(pv); rp = A[pv]
        for i in act:
            Ai = A[i]; aip = Ai[pv]
            for j in act:
                qq, rr = divmod(d * Ai[j] - aip * rp[j], prev)
                assert rr == 0
                Ai[j] = qq
        prev = d
    return True, 0

Kl = [[int(Kmat[i, j]) for j in range(nullity)] for i in range(nullity)]
ok, kern = psd_bareiss(Kl)
log(f'n = {n}: certificate Q restricted to ker J_z is {"POSITIVE SEMIDEFINITE" if ok else "NOT positive semidefinite"}'
    + (f', kernel dimension {kern} (expected n + 1 = {n + 1}, the M = 0 vectors of Sym^{2 * n})' if ok else '')
    + f'  [{time.time() - T0:.1f}s]')
