"""EXACT check of Proposition T2 (relative-anchors/README.md): the torus blocks of (R_inf)_{k-1} in which one copy t carries the letter 2
in row AND column.  All states are (e_t^{(2)}, y_h), y_h = (2 at t and on T', hole at h, 1 on U - h), h in U, where
U' = U + {t}.  Claim (unrescaled):
    Y = 1/4 I - 1/2 Tm^T E^{-1} Tm,   Tm y_h = 1/2 (w_h u_t - w_t u_h),
    E = diag(s_x + s_{T'} + s_{U'}/2) + 1/2 w_{U'} w_{U'}^T   on span{u_x : x in U'}  (u_x = 2 on T', hole at x, 1 on U'-x),
and Y is positive definite (proof in relative-anchors/README.md: <g,E^{-1}g> <= 2|c|^2).  Checked for k = 2..5, three rational weight
vectors each, every such block: closed form exactly, and exact LDL pivots > 0.  Prints PASS/FAIL."""
import itertools, random
import flint
from ylib import Builder, point_ring, rinf_states, wstate, Nv, inv
Q_ = flint.fmpq
random.seed(7)
ok = True
for k in range(2, 6):
    for trial in range(3):
        w = [Q_(random.randint(1, 12), random.randint(1, 12)) for _ in range(k)]
        if trial == 1:
            w[0] = Q_(40)          # one dominant weight
        s = [x * x for x in w]
        Bd = Builder(point_ring(s))
        S = rinf_states(k, k - 1)
        blocks = {}
        for X in S:
            blocks.setdefault(wstate(*X), []).append(X)
        nblk = 0
        for key, sts in blocks.items():
            twos = key[1:]
            if 2 not in twos:
                continue
            nblk += 1
            t = twos.index(2)
            Tp = [x for x in range(k) if twos[x] == 1]
            U = [x for x in range(k) if twos[x] == 0]
            Up = sorted(U + [t])
            assert all(a == tuple(2 if r == t else 0 for r in range(k)) for (a, y) in sts)
            hole = lambda y: next(i for i in range(k) if y[i] == 0)
            sts = sorted(sts, key=lambda X: hole(X[1]))
            hs = [hole(y) for (a, y) in sts]
            assert hs == U

            def mu(X):
                v = Q_(1)
                for tt, n in enumerate(Nv(*X)):
                    v *= w[tt] ** n
                return v
            Yb = [[Bd.Y(X1, X2) / (mu(X1) * mu(X2)) for X2 in sts] for X1 in sts]
            sTp = sum((s[x] for x in Tp), Q_(0))
            sUp = sum((s[x] for x in Up), Q_(0))
            n_ = len(Up)
            ix = {x: i for i, x in enumerate(Up)}
            E = [[(s[Up[a]] + sTp + sUp / 2 if a == b else Q_(0)) + w[Up[a]] * w[Up[b]] / 2 for b in range(n_)] for a in range(n_)]
            Ei = inv(E, Q_(0), Q_(1))
            Tm = []
            for h in U:
                v = [Q_(0)] * n_
                v[ix[t]] += w[h] / 2
                v[ix[h]] -= w[t] / 2
                Tm.append(v)
            pred = [[(Q_(1, 4) if i == j else Q_(0)) - sum((Tm[i][a] * Ei[a][b] * Tm[j][b] for a in range(n_) for b in range(n_)), Q_(0)) / 2
                     for j in range(len(U))] for i in range(len(U))]
            same = all(Yb[i][j] == pred[i][j] for i in range(len(U)) for j in range(len(U)))
            # exact LDL pivots
            Mx = [row[:] for row in Yb]
            pd = True
            for c in range(len(Mx)):
                if not Mx[c][c] > 0:
                    pd = False
                    break
                for i in range(c + 1, len(Mx)):
                    f = Mx[i][c] / Mx[c][c]
                    for j in range(c, len(Mx)):
                        Mx[i][j] -= f * Mx[c][j]
            ok &= same and pd
            if not (same and pd):
                print('FAIL at', k, key, same, pd)
        print(f'k={k} w={[str(x) for x in w]}: {nblk} blocks with a doubly occupied letter 2: closed form + PD checked')
print('PASS' if ok else 'FAIL')
