"""EXACT check (rational arithmetic, python-flint fmpq) of Theorem T (relative-anchors/README.md): the letter-1 torus block of (R_inf)_{k-1}.
States P_{r,h} = (row e_r, column 1^k - e_h) and Q = (row vac, column 1^k); UNRESCALED entries.  Rational weights w_t
(so s_t = w_t^2 are squares and every entry is rational).  Checks, for k = 2..6 and several weight vectors:
 (1) the closed form of the block:
       A[P_rh, P_r'h'] = 1/4 [(r,h)=(r',h'), r != h] + [r=h, r'=h'] (3/4 [r=r'] + 1/4),  A[Q,Q] = 1,  A[Q,P_rh] = -1/2 [r=h]
       B[P_rh, P_r'h'] = m_{r',h}^T D^{-1} m_{r,h'},  B[Q,P_rh] = -1/2 w_r w_h/|s|,  B[Q,Q] = 0,
       D = 1/2 (diag(d) + 1/2 w w^T),  d_h = s_h + |s|/2,  m_{r,h} = 1/4 (w_r f_h + w_h f_r) + 1/2 [r=h] w_h f_h;
 (2) the radical: Y (p (x) f_h) = 0 and Y (u (x) J^- |1^k>) = 0 for u ⊥ w (J^-|1^k> = sum_h w_h f_h);
 (3) the reduction: for random rational X with w^T X = 0 and X w = 0,
       min_q <Psi, Y Psi> = F(X) = 1/4 sum_ab X_ab^2 + 1/2 sum_a X_aa^2 - 1/2 sum_a rho_a X_aa^2
                                   + 1/(2(2+G)) sum_ab rho_a rho_b X_ab X_ba,   rho_a = s_a/d_a, G = sum rho_a;
 (4) F(X) >= c ||X||^2 with the explicit c of the proof (pointwise check at the random X).
Prints PASS/FAIL."""
import random
import flint
from ylib import Builder, point_ring, Nv

Q_ = flint.fmpq
random.seed(5)
ok = True


def rat(x):
    return x if isinstance(x, flint.fmpq) else Q_(x)


for k in range(2, 7):
    for trial in range(3 if k <= 5 else 1):
        w = [Q_(random.randint(1, 9), random.randint(1, 9)) for _ in range(k)]
        if trial == 0 and k >= 3:
            w[1] = Q_(1, 50)          # one small weight
        s = [x * x for x in w]
        S = sum(s, Q_(0))
        Bd = Builder(point_ring(s), letters=(0, 1))
        e = [tuple(1 if r == t else 0 for r in range(k)) for t in range(k)]
        hole = lambda h: tuple(0 if r == h else 1 for r in range(k))
        Pst = {(r, h): (e[r], hole(h)) for r in range(k) for h in range(k)}
        Qst = ((0,) * k, (1,) * k)
        labels = [(r, h) for r in range(k) for h in range(k)] + ['Q']
        st = [Pst[l] if l != 'Q' else Qst for l in labels]

        def mu(X):
            v = Q_(1)
            for t, n in enumerate(Nv(*X)):
                v *= w[t] ** n
            return v
        Yb = [[Bd.Y(X1, X2) / (mu(X1) * mu(X2)) for X2 in st] for X1 in st]     # unrescaled, exact
        # closed form
        d = [x + S / 2 for x in s]
        Dm = [[(d[a] if a == b else Q_(0)) / 2 + w[a] * w[b] / 4 for b in range(k)] for a in range(k)]
        # invert Dm exactly
        n = k
        Aug = [row[:] + [Q_(1) if i == j else Q_(0) for j in range(n)] for i, row in enumerate(Dm)]
        Aug = [[rat(x) for x in row] for row in Aug]
        for c in range(n):
            piv = Aug[c][c]
            Aug[c] = [x / piv for x in Aug[c]]
            for i in range(n):
                if i != c and Aug[i][c] != 0:
                    f = Aug[i][c]
                    Aug[i] = [x - f * y for x, y in zip(Aug[i], Aug[c])]
        Dinv = [row[n:] for row in Aug]

        def m(r, h):
            v = [Q_(0)] * k
            v[h] += w[r] / 4
            v[r] += w[h] / 4
            if r == h:
                v[h] += w[h] / 2
            return v

        def Bp(r, h, r2, h2):
            a_, b_ = m(r2, h), m(r, h2)
            return sum((a_[i] * Dinv[i][j] * b_[j] for i in range(k) for j in range(k)), Q_(0))

        def Ap(r, h, r2, h2):
            v = Q_(0)
            if (r, h) == (r2, h2) and r != h:
                v += Q_(1, 4)
            if r == h and r2 == h2:
                v += Q_(3, 4) * (r == r2) + Q_(1, 4)
            return v
        Yc = []
        for l1 in labels:
            row = []
            for l2 in labels:
                if l1 == 'Q' and l2 == 'Q':
                    row.append(Q_(1))
                elif l1 == 'Q' or l2 == 'Q':
                    r, h = l2 if l1 == 'Q' else l1
                    row.append(-Q_(1, 2) * (r == h) + w[r] * w[h] / (2 * S))
                else:
                    row.append(Ap(*l1, *l2) - Bp(*l1, *l2))
            Yc.append(row)
        same = all(Yb[i][j] == Yc[i][j] for i in range(len(st)) for j in range(len(st)))
        # radical vectors
        N = len(st)
        idx = {l: i for i, l in enumerate(labels)}
        rad_ok = True
        vecs = []
        for h0 in range(k):
            v = [Q_(0)] * N
            for r in range(k):
                v[idx[(r, h0)]] = w[r]
            vecs.append(v)
        for t in range(1, k):
            u = [Q_(0)] * k
            u[0], u[t] = w[t], -w[0]
            v = [Q_(0)] * N
            for r in range(k):
                for h in range(k):
                    v[idx[(r, h)]] = u[r] * w[h]
            vecs.append(v)
        for v in vecs:
            rad_ok &= all(sum((Yb[i][j] * v[j] for j in range(N)), Q_(0)) == 0 for i in range(N))
        # reduction (3) and bound (4)
        rho = [s[a] / d[a] for a in range(k)]
        G = sum(rho, Q_(0))
        red_ok = True
        for _ in range(4):
            X = [[Q_(random.randint(-9, 9), random.randint(1, 5)) for _ in range(k)] for _ in range(k)]
            # project: X <- (I - w w^T/S) X (I - w w^T/S)  (exact)
            P = [[(Q_(1) if a == b else Q_(0)) - w[a] * w[b] / S for b in range(k)] for a in range(k)]
            X = [[sum((P[a][c] * X[c][e_] * P[e_][b] for c in range(k) for e_ in range(k)), Q_(0)) for b in range(k)]
                 for a in range(k)]
            x = [X[r][h] for r in range(k) for h in range(k)]
            quad = sum((x[i] * Yb[i][j] * x[j] for i in range(k * k) for j in range(k * k)), Q_(0))
            cQ = sum((Yb[N - 1][j] * x[j] for j in range(k * k)), Q_(0))
            val = quad - cQ * cQ / Yb[N - 1][N - 1]
            F = sum((X[a][b] ** 2 for a in range(k) for b in range(k)), Q_(0)) / 4 \
                + sum((X[a][a] ** 2 * (1 - rho[a]) for a in range(k)), Q_(0)) / 2 \
                + sum((rho[a] * rho[b] * X[a][b] * X[b][a] for a in range(k) for b in range(k)), Q_(0)) / (2 * (2 + G))
            red_ok &= (val == F)
            # bound: F >= sum_{a<b} (1/4 - rho_a rho_b/(2(2+G))) (X_ab^2 + X_ba^2) + sum_a (3/4 - rho_a/2) X_aa^2 >= 0
            lb = sum(((Q_(1, 4) - rho[a] * rho[b] / (2 * (2 + G))) * (X[a][b] ** 2 + X[b][a] ** 2)
                      for a in range(k) for b in range(a + 1, k)), Q_(0)) \
                + sum(((Q_(3, 4) - rho[a] / 2) * X[a][a] ** 2 for a in range(k)), Q_(0))
            red_ok &= (F >= lb) and all(Q_(1, 4) - rho[a] * rho[b] / (2 * (2 + G)) > 0 for a in range(k) for b in range(k))
        good = same and rad_ok and red_ok
        ok &= good
        print(f'k={k} w={[str(x) for x in w]}: closed form exact {same}, radical {rad_ok}, reduction+bound {red_ok}')
print('PASS' if ok else 'FAIL')
