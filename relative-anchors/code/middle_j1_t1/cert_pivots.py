"""Exact certificate for the |T| = 1 block of (R_inf)_1 at a fixed k and ALL positive weight vectors, through the
inverse-free congruent form Phi' = diag(A1, A2) Phi diag(A1, A2) of the P-first reduction (Lemma P; A1 = I + N1/2,
A2 = I + N'/2 are positive definite, so Phi' > 0 iff Phi > 0).

python cert_pivots.py k [order=...]
lev.py with l = k - 1 levels of size one (fixed), p = s_t = 1 (every rescaled entry of Phi' is homogeneous of degree 0
in the squared weights, so this loses nothing), variables q_1..q_{k-1} (the squared weights of the copies of U).  The
matrix is D Phi' D in the rescaled frame, of size 2k - 1, coordinates (rho_t, rho_1..rho_u, q_1..q_u).
 (1) validation: at two rational weight vectors with |w| = 1, D Phi' D computed from j1lib (Lemma P) and the closed
     A1, A2 of tlib.py (Lemma T1.3) equals the symbolic matrix evaluated there;
 (2) certificate: symmetric Gaussian elimination over Q(q_1..q_u); every pivot (ratio of consecutive leading principal
     minors) is checked to be a ratio of polynomials with coefficients of one sign, positive; then every leading
     principal minor is positive on the open orthant and Phi' > 0 (Sylvester).
Prints PASS if (1) and (2) hold."""
import sys, time, random
from fractions import Fraction as Fr
import flint
from lev import build, signcheck
from j1lib import J1, block_states
from tlib import closed

t0 = time.time()
k = int(sys.argv[1])
opts = dict(a.split('=', 1) for a in sys.argv[2:] if '=' in a)
u = k - 1
S = build(u, [1] * u, 1, want_phi=False, fixed=[True] * u, p_one=True, want_phic=True)
K = S['K']
A = S['Cinv']
n = len(A)
order = list(range(n))
if 'order' in opts:
    order = [int(v) for v in opts['order'].split(',')]
    A = [[A[i][j] for j in order] for i in order]
print(f'k = {k}: Phi\' of size {n} over Q({",".join(S["names"][1:])}), order {order} [{time.time() - t0:.1f}s]', flush=True)
ok = all(A[i][j] == A[j][i] for i in range(n) for j in range(n))


def ev(f, vals):
    def pr(m):
        out = Fr(1)
        for v, e in zip(vals, m):
            out *= v ** int(e)
        return out
    num = sum((Fr(int(flint.fmpq(c).p), int(flint.fmpq(c).q)) * pr(m) for m, c in f.n.to_dict().items()), Fr(0))
    den = sum((Fr(int(flint.fmpq(c).p), int(flint.fmpq(c).q)) * pr(m) for m, c in f.d.to_dict().items()), Fr(0))
    return num / den


def unit_weights(kk):
    v = [Fr(random.randint(1, 9), random.randint(1, 9)) for _ in range(kk - 1)]
    n2 = sum(a * a for a in v)
    w = [abs(2 * a / (1 + n2)) for a in v] + [abs((1 - n2) / (1 + n2))]
    return w if all(a != 0 for a in w) else unit_weights(kk)


def mul(X, Y):
    return [[sum(X[i][l] * Y[l][j] for l in range(len(Y))) for j in range(len(Y[0]))] for i in range(len(X))]


# (1) validation against j1lib + tlib
random.seed(29)
for trial in range(2):
    w = unit_weights(k)
    x = [a * a for a in w]
    st = block_states(k, (0,))
    M = J1(w).matrix(st)
    R = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 2]
    P = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 1]
    Q = [i for i, (a, y) in enumerate(st) if a is None]
    E = R + Q
    N = [Fr(1)] + [Fr(2)] * u
    Phi = [[M[E[a]][E[b]] - sum(M[P[r]][E[a]] * M[P[r]][E[b]] / N[r] for r in range(k)) for b in range(n)] for a in range(n)]
    Cf = closed(w, 0)
    Dg = [[Fr(0)] * n for _ in range(n)]
    for i in range(k):
        for j in range(k):
            Dg[i][j] = Cf['A1'][i][j]
    for i in range(u):
        for j in range(u):
            Dg[k + i][k + j] = Cf['A2'][i][j]
    Php = mul(mul(Dg, Phi), Dg)
    mu = [w[c] for c in range(k)] + [w[c] for c in range(1, k)]     # R_c: w_c; Q_c: w_c
    resc = [[mu[i] * Php[i][j] * mu[j] for j in range(n)] for i in range(n)]
    # the symbolic entries use the normalised weights x = (p, q)/L, so they are homogeneous of degree 0 in (p, q):
    # evaluate at p = 1, q_i = x_i / x_t
    vals = [Fr(1)] + [x[c] / x[0] for c in range(1, k)]
    Sym = [[ev(S['Cinv'][i][j], vals) for j in range(n)] for i in range(n)]
    good = all(Sym[i][j] == resc[i][j] for i in range(n) for j in range(n))
    ok &= good
    print(f'  (1) validation at w = {[str(a) for a in w]}: {good} [{time.time() - t0:.1f}s]', flush=True)

# (2) pivots
Mx = [row[:] for row in A]
for c in range(n):
    pv = Mx[c][c]
    sg, desc = signcheck(pv)
    good = sg == 1
    ok &= good
    print(f'  (2) pivot {c + 1}/{n}: numerator deg {pv.n.total_degree()} ({len(pv.n)} terms), denominator deg '
          f'{pv.d.total_degree()} ({len(pv.d)} terms); sign {sg} {desc[:80]} [{time.time() - t0:.1f}s]', flush=True)
    if not good:
        break
    for i in range(c + 1, n):
        if Mx[i][c] != 0:
            f = Mx[i][c] / pv
            for j in range(c + 1, n):
                if Mx[c][j] != 0:
                    Mx[i][j] = Mx[i][j] - f * Mx[c][j]
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
