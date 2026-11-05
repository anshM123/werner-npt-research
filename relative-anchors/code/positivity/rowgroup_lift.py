"""Exact check of the explicit lifts of Proposition R, at rational points with s_t = w_t^2 (w rational), k = 2, 3, 4,
every j, every torus sector.

Letter 1 at r (row a = 1@r), j = 0..k, coordinates psi'_y = s^{N(y)} x_(a,y) (x = ylib rescaled coordinates) and z'_x
(x in H_j):
  B'_a = sum_y (Delta(y) s_r / s^{N(y)}) f_y f_y^T + sum_x (C^(r)_x / s^{N(x)}) e_x e_x^T
         + (1/2) sum_{u in H_{j-1}, u_r != 0} s^{-N(u)} g_u g_u^T,
  f_y = e_{psi', y} + e_{z', y} + [y_r = 1] sum_{t: y_t = 0} e_{z', y - 1_r + 1_t},   g_u = sum_{t != r, u_t = 0} e_{z', u + 1_t},
  Delta(y) = 2, 1/2, 1 for y_r = 0, 1, 2, and C^(r)_x = sum_{t != r} s_t (2[x_t = 0] + [x_t = 2]).
  Claim: the Schur complement of B'_a onto the psi' coordinates equals 2^{k-j} Y_ylib((a,y),(a,y')) / (s^{N(y)} s^{N(y')}).
Vacuum row, j = -1..k-1 (columns y in H_{j+1}), coordinates phi'_y = s^{N(y)} x_(vac,y), z'_x (x in H_{j+2}):
  B'_vac = sum_y s^{-N(y)} f_y f_y^T + sum_x (2 C_x / s^{N(x)}) e_x e_x^T,   f_y = e_{phi', y} + sum_{t: y_t = 0} e_{z', y + 1_t}.
  Claim: the Schur complement onto phi' equals 2^{k-j-1} Y_ylib((vac,y),(vac,y')) / (s^{N(y)} s^{N(y')}).
Every weight is a Laurent monomial or a linear form with nonnegative coefficients, and every vector has entries 0 or 1,
so by Cauchy-Binet every principal minor of the row block is subtraction-free (Proposition R).
Usage: python rowgroup_lift.py"""
import itertools, random
import flint
from ylib import Builder, point_ring

Q = flint.fmpq


def strings(k, c):
    return [z for z in itertools.product((0, 1, 2), repeat=k) if sum(x != 0 for x in z) == c]


def schur_onto(Bm, keep, elim):
    nk, ne = len(keep), len(elim)
    if ne == 0:
        return [[Bm[a][b] for b in keep] for a in keep]
    E = flint.fmpq_mat(ne, ne, [Bm[a][b] for a in elim for b in elim])
    W = flint.fmpq_mat(ne, nk, [Bm[a][b] for a in elim for b in keep])
    A = flint.fmpq_mat(nk, nk, [Bm[a][b] for a in keep for b in keep])
    S = A - W.transpose() * E.solve(W)
    return [[S[i, j] for j in range(nk)] for i in range(nk)]


rng = random.Random(17)
ok = True
nchk = 0
for k in (2, 3, 4):
    for trial in range(2):
        w = [Q(rng.randint(1, 7), rng.randint(1, 7)) for _ in range(k)]
        s = [x * x for x in w]
        Bp = Builder(point_ring(s))

        def sN(z):
            v = Q(1)
            for t in range(k):
                if z[t] == 1:
                    v *= s[t]
            return v
        # letter-1 rows
        for j in range(0, k + 1):
            Hj, Hjm = strings(k, j), strings(k, j - 1)
            for r in range(k):
                a = tuple(1 if t == r else 0 for t in range(k))
                idx = {('p', y): i for i, y in enumerate(Hj)}
                for i, x in enumerate(Hj):
                    idx[('z', x)] = len(Hj) + i
                n = 2 * len(Hj)
                Bm = [[Q(0)] * n for _ in range(n)]

                def add(weight, vec):
                    for p, cp in vec.items():
                        for q_, cq in vec.items():
                            Bm[p][q_] += weight * cp * cq
                for y in Hj:
                    D = {0: Q(2), 1: Q(1, 2), 2: Q(1)}[y[r]]
                    f = {idx[('p', y)]: 1}
                    f[idx[('z', y)]] = f.get(idx[('z', y)], 0) + 1
                    if y[r] == 1:
                        for t in range(k):
                            if y[t] == 0:
                                y2 = list(y); y2[r] = 0; y2[t] = 1; y2 = tuple(y2)
                                f[idx[('z', y2)]] = f.get(idx[('z', y2)], 0) + 1
                    add(D * s[r] / sN(y), f)
                for x in Hj:
                    Cr = sum((s[t] * (2 * (x[t] == 0) + (x[t] == 2)) for t in range(k) if t != r), Q(0))
                    add(Cr / sN(x), {idx[('z', x)]: 1})
                for u in Hjm:
                    if u[r] != 0:
                        g = {}
                        for t in range(k):
                            if t != r and u[t] == 0:
                                u2 = u[:t] + (1,) + u[t + 1:]
                                g[idx[('z', u2)]] = 1
                        if g:
                            add(Q(1, 2) / sN(u), g)
                keep = [idx[('p', y)] for y in Hj]
                elim = [idx[('z', x)] for x in Hj]
                S = schur_onto(Bm, keep, elim)
                for i, y in enumerate(Hj):
                    for l, y2 in enumerate(Hj):
                        rhs = Q(2) ** (k - j) * Bp.Y((a, y), (a, y2)) / (sN(y) * sN(y2))
                        good = S[i][l] == rhs
                        ok &= good
                        nchk += 1
                        if not good:
                            print('letter-1 mismatch', k, j, r, y, y2, S[i][l], rhs)
                            break
        # vacuum row
        for j in range(-1, k):
            Hy, Hx = strings(k, j + 1), strings(k, j + 2) if j + 2 <= k else []
            idx = {('p', y): i for i, y in enumerate(Hy)}
            for i, x in enumerate(Hx):
                idx[('z', x)] = len(Hy) + i
            n = len(Hy) + len(Hx)
            Bm = [[Q(0)] * n for _ in range(n)]

            def add(weight, vec):
                for p, cp in vec.items():
                    for q_, cq in vec.items():
                        Bm[p][q_] += weight * cp * cq
            for y in Hy:
                f = {idx[('p', y)]: 1}
                for t in range(k):
                    if y[t] == 0 and j + 2 <= k:
                        y2 = y[:t] + (1,) + y[t + 1:]
                        f[idx[('z', y2)]] = 1
                add(1 / sN(y), f)
            for x in Hx:
                Cx = sum((s[t] * (2 * (x[t] == 0) + (x[t] == 2)) for t in range(k)), Q(0))
                add(2 * Cx / sN(x), {idx[('z', x)]: 1})
            keep = [idx[('p', y)] for y in Hy]
            elim = [idx[('z', x)] for x in Hx]
            vac = (0,) * k
            if any(Bm[e][e] == 0 for e in elim):
                # z-block singular (x = 1...1 with C = 0 and no f touching it is impossible; keep as a check)
                pass
            S = schur_onto(Bm, keep, elim)
            for i, y in enumerate(Hy):
                for l, y2 in enumerate(Hy):
                    rhs = Q(2) ** (k - j - 1) * Bp.Y((vac, y), (vac, y2)) / (sN(y) * sN(y2))
                    good = S[i][l] == rhs
                    ok &= good
                    nchk += 1
                    if not good:
                        print('vacuum mismatch', k, j, y, y2, S[i][l], rhs)
                        break
    print(f'k = {k}: {nchk} entries checked so far, all equal: {ok}', flush=True)
print('PASS' if ok else 'FAIL')
