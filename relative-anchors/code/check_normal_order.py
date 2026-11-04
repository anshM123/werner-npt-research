"""EXACT checks of Lemma N and Proposition K2 (relative-anchors/README.md).  Rational weights w, exact fmpq arithmetic.
Lemma N (normal-ordered form), on H_c for every c and all letters:
   T_(r,1) = 1/2 sigma_(r,1)^- J^+ + w_r N_r,   N_r = 3/2 |0><0|_r + 1/2 |1><1|_r + |2><2|_r,
   T_(t,2) = 1/2 sigma_(t,2)^- J^+ - 1/2 w_t |1><2|_t,
   E_c = Nn + 1/2 J^- J^+,   Nn = sum_t s_t N_t   (diagonal, >= s_min/2 > 0).
   Consequence: on ker J^+ every T_a is local (no hopping).
Proposition K2: for every k, the vector vac (x) J^-|1...1> (row vacuum, column J^-|1..1> in H_{k-1}) lies in the kernel
   of (R_inf)_{k-2}; checked as Y v = 0 exactly on its torus block (letter-1 block), k = 3..6.
Prints PASS/FAIL."""
import itertools, random
import flint
from ylib import Builder, point_ring, rinf_states, Nv
Q_ = flint.fmpq
random.seed(13)
ok = True


def zeros(n, m):
    return [[Q_(0)] * m for _ in range(n)]


def matmul(X, Yb):
    n, m = len(X), len(Yb)
    p = len(Yb[0]) if m else 0
    out = zeros(n, p)
    for i in range(n):
        for l in range(m):
            if X[i][l] != 0:
                xi = X[i][l]
                for jj in range(p):
                    if Yb[l][jj] != 0:
                        out[i][jj] += xi * Yb[l][jj]
    return out


for k in (2, 3, 4):
    w = [Q_(random.randint(1, 9), random.randint(1, 9)) for _ in range(k)]
    s = [x * x for x in w]
    cols = list(itertools.product(range(3), repeat=k))
    q = lambda y: sum(x != 0 for x in y)
    H = {c: [y for y in cols if q(y) == c] for c in range(k + 1)}
    idx = {c: {y: i for i, y in enumerate(H[c])} for c in H}

    def setl(y, r, l):
        y = list(y); y[r] = l; return tuple(y)

    def op(c_from, c_to, f):
        Mx = zeros(len(H[c_to]), len(H[c_from]))
        for y in H[c_from]:
            for coef, y2 in f(y):
                Mx[idx[c_to][y2]][idx[c_from][y]] += coef
        return Mx
    Jp = lambda c: op(c, c + 1, lambda y: [(w[t], setl(y, t, 1)) for t in range(k) if y[t] == 0])
    Jm = lambda c: op(c, c - 1, lambda y: [(w[t], setl(y, t, 0)) for t in range(k) if y[t] == 1])
    good = True
    for c in range(0, k + 1):
        n = len(H[c])
        Cd = op(c, c, lambda y: [(sum((s[t] * (2 * (y[t] == 0) + (y[t] == 2)) for t in range(k)), Q_(0)), y)])
        JpJm = matmul(Jp(c - 1), Jm(c)) if c >= 1 else zeros(n, n)
        E = [[Cd[i][l] + JpJm[i][l] / 2 for l in range(n)] for i in range(n)]
        Nn = op(c, c, lambda y: [(sum((s[t] * (Q_(3, 2) * (y[t] == 0) + Q_(1, 2) * (y[t] == 1) + (y[t] == 2)) for t in range(k)), Q_(0)), y)])
        JmJp = matmul(Jm(c + 1), Jp(c)) if c + 1 <= k else zeros(n, n)
        good &= all(E[i][l] == Nn[i][l] + JmJp[i][l] / 2 for i in range(n) for l in range(n))
        for r in range(k):
            for l in (1, 2):
                # T_a as in Lemma B
                sm = op(c, c - 1, lambda y: [(Q_(1), setl(y, r, 0))] if y[r] == l else []) if c >= 1 else None
                M1 = matmul(Jp(c - 1), sm) if c >= 1 else zeros(n, n)
                if l == 1:
                    M2 = op(c, c, lambda y: [(w[r] * (2 * (y[r] == 0) + (y[r] == 2)), y)])
                else:
                    M2 = op(c, c, lambda y: [(-w[r], setl(y, r, 1))] if y[r] == 2 else [])
                Tm = [[M1[i][j_] / 2 + M2[i][j_] for j_ in range(n)] for i in range(n)]
                # normal-ordered prediction
                sm_up = op(c + 1, c, lambda y: [(Q_(1), setl(y, r, 0))] if y[r] == l else []) if c + 1 <= k else None
                P1 = matmul(sm_up, Jp(c)) if c + 1 <= k else zeros(n, n)
                if l == 1:
                    P2 = op(c, c, lambda y: [(w[r] * (Q_(3, 2) * (y[r] == 0) + Q_(1, 2) * (y[r] == 1) + (y[r] == 2)), y)])
                else:
                    P2 = op(c, c, lambda y: [(-w[r] / 2, setl(y, r, 1))] if y[r] == 2 else [])
                pred = [[P1[i][j_] / 2 + P2[i][j_] for j_ in range(n)] for i in range(n)]
                good &= all(Tm[i][j_] == pred[i][j_] for i in range(n) for j_ in range(n))
    ok &= good
    print(f'k={k}: Lemma N (T_a normal-ordered, E = Nn + J^-J^+/2) on every H_c, all letters: {good}', flush=True)

for k in (3, 4, 5, 6):
    w = [Q_(random.randint(1, 9), random.randint(1, 9)) for _ in range(k)]
    s = [x * x for x in w]
    Bd = Builder(point_ring(s), letters=(0, 1))
    S = rinf_states(k, k - 2, letters=(0, 1))
    vac = (0,) * k
    v = {}
    for h in range(k):
        y = tuple(0 if r == h else 1 for r in range(k))
        v[(vac, y)] = w[h]                       # J^- |1..1> = sum_h w_h |hole at h>

    def mu(X):
        val = Q_(1)
        for t, nn in enumerate(Nv(*X)):
            val *= w[t] ** nn
        return val
    zero = True
    for X in S:
        tot = sum((Bd.Y(X, X2) / (mu(X) * mu(X2)) * c2 for X2, c2 in v.items()), Q_(0))
        zero &= (tot == 0)
    ok &= zero
    print(f'k={k}: Y (vac (x) J^-|1..1>) = 0 on the letter-1 block of (R_inf)_{k-2} ({len(S)} states): {zero}', flush=True)
print('PASS' if ok else 'FAIL')
