"""EXACT check of Lemma B (relative-anchors/README.md): the hard-core-boson form of the whole (R_inf)_j matrix, every j, all letters 0,1,2.
Unrescaled, rational weights w (exact fmpq).  Column space H_c (strings with c nonzero letters), Z = k - j.
For single excitations a = (letter l at copy r), b = (l' at r'):  sigma_a^- = |0><l| at r, sigma_b^+ = |l'><0| at r',
n_r = 1 - |0><0|_r, J^+ = sum_t w_t |1><0|_t, J^- = (J^+)^T,
C = sum_t s_t (2|0><0|_t + |2><2|_t),  E_c = C + 1/2 J^+ J^-  on H_c  (2^{k-c} D_c = E_c, Lemma 1 of relative-anchors/README.md),
T_a = 1/2 J^+ sigma_a^- - w_r |1><l|_r [l = 2] + w_r (2|0><0|_r + |2><2|_r) [l = 1]   (on H_c).
Claim:
  2^Z Y[(a,y'),(b,y)]   = <y'| (1/2 - [r = r']) sigma_b^+ sigma_a^- + [a = b](2 - n_r) - T_b^T E_j^{-1} T_a |y>,
  2^Z Y[(a,y'),(vac,y)] = -<y'| sigma_a^- - J^- E_{j+1}^{-1} T_a |y>                       (y' in H_j, y in H_{j+1}),
  2^Z Y[(vac,y'),(vac,y)] = 2 <y'| I - 1/2 J^- E_{j+2}^{-1} J^+ |y>                        (y', y in H_{j+1}).
k = 2, 3, 4, every j.  Prints PASS/FAIL."""
import itertools, random
import flint
from ylib import Builder, point_ring, rinf_states, Nv, inv
Q_ = flint.fmpq
random.seed(11)
ok = True
for k in (2, 3, 4):
    w = [Q_(random.randint(1, 9), random.randint(1, 9)) for _ in range(k)]
    s = [x * x for x in w]
    Bd = Builder(point_ring(s))
    cols = list(itertools.product(range(3), repeat=k))
    q = lambda y: sum(x != 0 for x in y)
    H = {c: [y for y in cols if q(y) == c] for c in range(k + 1)}
    idx = {c: {y: i for i, y in enumerate(H[c])} for c in H}

    def zeros(n, m):
        return [[Q_(0)] * m for _ in range(n)]

    def matmul(X, Yb):
        n, m, p = len(X), len(Yb), len(Yb[0]) if Yb else 0
        out = zeros(n, p)
        for i in range(n):
            for l in range(m):
                if X[i][l] != 0:
                    xi = X[i][l]
                    for jj in range(p):
                        if Yb[l][jj] != 0:
                            out[i][jj] += xi * Yb[l][jj]
        return out

    def tr(X):
        return [list(r) for r in zip(*X)] if X else []

    def op(c_from, c_to, f):
        """matrix of a map H_{c_from} -> H_{c_to}; f(y) returns list of (coef, y2)."""
        Mx = zeros(len(H[c_to]), len(H[c_from]))
        for y in H[c_from]:
            for coef, y2 in f(y):
                Mx[idx[c_to][y2]][idx[c_from][y]] += coef
        return Mx

    def setl(y, r, l):
        y = list(y); y[r] = l; return tuple(y)

    def Jp(c):   # H_c -> H_{c+1}
        return op(c, c + 1, lambda y: [(w[t], setl(y, t, 1)) for t in range(k) if y[t] == 0])

    def Jm(c):   # H_c -> H_{c-1}
        return op(c, c - 1, lambda y: [(w[t], setl(y, t, 0)) for t in range(k) if y[t] == 1])

    def E(c):
        if c > k or c < 0:
            return None
        Cd = op(c, c, lambda y: [(sum((s[t] * (2 * (y[t] == 0) + (y[t] == 2)) for t in range(k)), Q_(0)), y)])
        if c == 0:
            return Cd
        JJ = matmul(Jp(c - 1), Jm(c))
        return [[Cd[i][l] + JJ[i][l] / 2 for l in range(len(Cd))] for i in range(len(Cd))]

    def sig_minus(a, c):   # |0><l| at r : H_c -> H_{c-1}
        r, l = a
        return op(c, c - 1, lambda y: [(Q_(1), setl(y, r, 0))] if y[r] == l else [])

    def sig_plus(b, c):    # |l'><0| at r' : H_c -> H_{c+1}
        r, l = b
        return op(c, c + 1, lambda y: [(Q_(1), setl(y, r, l))] if y[r] == 0 else [])

    def Tmat(a, c):
        r, l = a
        M1 = matmul(Jp(c - 1), sig_minus(a, c)) if c >= 1 else zeros(len(H[c]), len(H[c]))
        M1 = [[x / 2 for x in row] for row in M1]
        if l == 2:
            M2 = op(c, c, lambda y: [(-w[r], setl(y, r, 1))] if y[r] == 2 else [])
        else:
            M2 = op(c, c, lambda y: [(w[r] * (2 * (y[r] == 0) + (y[r] == 2)), y)])
        return [[M1[i][l_] + M2[i][l_] for l_ in range(len(M1))] for i in range(len(M1))]

    def nmat(r, c):
        return op(c, c, lambda y: [(Q_(int(y[r] != 0)), y)])

    def single(a):
        r = next(i for i in range(k) if a[i]); return (r, a[r])

    for j in range(-1, k + 1):
        S = rinf_states(k, j)
        if not S:
            continue
        Z = k - j
        Ej_inv = inv(E(j), Q_(0), Q_(1)) if 0 <= j <= k else None
        Ej1_inv = inv(E(j + 1), Q_(0), Q_(1)) if j + 1 <= k else None
        Ej2_inv = inv(E(j + 2), Q_(0), Q_(1)) if j + 2 <= k else None
        cache = {}
        good = True
        for X1 in S:
            for X2 in S:
                (a, yp), (b, y) = X1, X2
                va, vb = sum(a) == 0, sum(b) == 0
                if not va and not vb:
                    A_, B_ = single(a), single(b)
                    key = ('11', A_, B_)
                    if key not in cache:
                        c = j
                        P1 = matmul(sig_plus(B_, c - 1), sig_minus(A_, c)) if c >= 1 else zeros(len(H[c]), len(H[c]))
                        coef = Q_(1, 2) - (A_[0] == B_[0])
                        Mx = [[coef * x for x in row] for row in P1]
                        if A_ == B_:
                            nr = nmat(A_[0], c)
                            for i in range(len(Mx)):
                                Mx[i][i] += 2 - nr[i][i]
                        TT = matmul(tr(Tmat(B_, c)), matmul(Ej_inv, Tmat(A_, c)))
                        cache[key] = [[Mx[i][l] - TT[i][l] for l in range(len(Mx))] for i in range(len(Mx))]
                    pred = cache[key][idx[j][yp]][idx[j][y]] / 2 ** Z
                elif not va and vb:
                    A_ = single(a)
                    key = ('10', A_)
                    if key not in cache:
                        sm = sig_minus(A_, j + 1)                       # H_{j+1} -> H_j
                        JEt = matmul(Jm(j + 1), matmul(Ej1_inv, Tmat(A_, j + 1)))
                        cache[key] = [[-(sm[i][l] - JEt[i][l]) for l in range(len(sm[0]))] for i in range(len(sm))]
                    pred = cache[key][idx[j][yp]][idx[j + 1][y]] / 2 ** Z
                elif va and not vb:
                    B_ = single(b)
                    key = ('10', B_)
                    if key not in cache:
                        sm = sig_minus(B_, j + 1)
                        JEt = matmul(Jm(j + 1), matmul(Ej1_inv, Tmat(B_, j + 1)))
                        cache[key] = [[-(sm[i][l] - JEt[i][l]) for l in range(len(sm[0]))] for i in range(len(sm))]
                    pred = cache[key][idx[j][y]][idx[j + 1][yp]] / 2 ** Z
                else:
                    key = ('00',)
                    if key not in cache:
                        n1 = len(H[j + 1])
                        if Ej2_inv is not None:
                            JEJ = matmul(Jm(j + 2), matmul(Ej2_inv, Jp(j + 1)))
                        else:
                            JEJ = zeros(n1, n1)
                        cache[key] = [[2 * ((Q_(1) if i == l else Q_(0)) - JEJ[i][l] / 2) for l in range(n1)] for i in range(n1)]
                    pred = cache[key][idx[j + 1][yp]][idx[j + 1][y]] / 2 ** Z
                ex = Bd.Y(X1, X2)
                muv = Q_(1)
                for t, n in enumerate(Nv(*X1)):
                    muv *= w[t] ** n
                for t, n in enumerate(Nv(*X2)):
                    muv *= w[t] ** n
                if ex / muv != pred:
                    good = False
        ok &= good
        print(f'k={k} j={j}: {len(S)} states, boson form exact: {good}', flush=True)
print('PASS' if ok else 'FAIL')
