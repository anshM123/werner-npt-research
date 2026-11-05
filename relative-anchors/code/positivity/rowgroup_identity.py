"""Exact checks of the algebra behind Proposition R (row groups), at rational points with s_t = w_t^2 (w rational).

Operators on column strings (README of relative-anchors, Section 1): J^+, J^-, C, E_c = C + J^+J^-/2, sigma_a, T_a,
and Delta_a = (2 - n_r) - (1/2) sigma_a^+ sigma_a^-  (the diagonal part of Lemma B for a = b).  Checks, for k = 3, 4,
every j = 0..k and every row a (letter 1 or 2 at copy r):
 (i)  Lemma B on the diagonal row block: 2^Z Y[(a,y'),(a,y)] = <y'| Delta_a - T_a^T E_j^{-1} T_a |y>, compared with
      the ylib point builder (rescaled entries divided by mu_X mu_X');
 (ii) letter 1 at r:  E_j - T_a Delta_a^{-1} T_a^T = C^(r) + (1/2) J^+_(r) (1 - P_{r,0}) J^-_(r),
      C^(r) = sum_{t != r} s_t (2|0><0|_t + |2><2|_t), J_(r) = J without copy r, P_{r,0} = projector on z_r = 0;
 (iii) letter 2 at r:  E_j - T_a Delta_a^{-1} T_a^T = C + (1/2) J^+_(r) (1 - P_{r,0}) J^-_(r)
      + sum_{t != r} w_r w_t (h_{t<-r} + h_{r<-t}),  h_{t<-r} moving a letter 1 from r to t (z_t = 0),
      and this matrix is positive semidefinite (exact LDL with pivoting) at the point.
Usage: python rowgroup_identity.py"""
import itertools, random
import flint
from ylib import Builder, point_ring, rinf_states, wstate

Q = flint.fmpq


def strings(k, c):
    return [z for z in itertools.product((0, 1, 2), repeat=k) if sum(x != 0 for x in z) == c]


def zeros(m, n):
    return [[Q(0)] * n for _ in range(m)]


def mat_mul(A, B):
    n, m, p = len(A), len(B), len(B[0]) if B else 0
    out = zeros(n, p)
    for i in range(n):
        Ai = A[i]
        for l in range(m):
            if Ai[l] != 0:
                Bl = B[l]
                for j in range(p):
                    if Bl[j] != 0:
                        out[i][j] += Ai[l] * Bl[j]
    return out


def T(A):
    return [list(r) for r in zip(*A)] if A else A


def inv(A):
    n = len(A)
    M = flint.fmpq_mat(n, n, [x for r in A for x in r])
    Mi = M.inv()
    return [[Mi[i, j] for j in range(n)] for i in range(n)]


def psd(A):
    n = len(A)
    A = [r[:] for r in A]
    act = list(range(n))
    while True:
        cand = [i for i in act if A[i][i] > 0]
        if not cand:
            break
        i = cand[0]
        rest = [r for r in act if r != i]
        for r in rest:
            if A[r][i] != 0:
                f = A[r][i] / A[i][i]
                for l in rest:
                    A[r][l] -= f * A[i][l]
        act = rest
    return all(A[i][i] == 0 for i in act) and all(A[i][j] == 0 for i in act for j in act)


rng = random.Random(5)
ok = True
count = dict(i=0, ii=0, iii=0)
for k in (3, 4):
    for trial in range(2):
        w = [Q(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(k)]
        s = [x * x for x in w]
        Bp = Builder(point_ring(s))
        H = {c: strings(k, c) for c in range(-1, k + 2)}
        idx = {c: {z: i for i, z in enumerate(H[c])} for c in H}

        def Jp(c):  # H_c -> H_{c+1}
            A = zeros(len(H[c + 1]), len(H[c]))
            for z in H[c]:
                for t in range(k):
                    if z[t] == 0:
                        z2 = z[:t] + (1,) + z[t + 1:]
                        A[idx[c + 1][z2]][idx[c][z]] += w[t]
            return A

        def Jp_excl(c, r):
            A = zeros(len(H[c + 1]), len(H[c]))
            for z in H[c]:
                for t in range(k):
                    if t != r and z[t] == 0:
                        z2 = z[:t] + (1,) + z[t + 1:]
                        A[idx[c + 1][z2]][idx[c][z]] += w[t]
            return A

        def Cdiag(c, skip=None):
            A = zeros(len(H[c]), len(H[c]))
            for z in H[c]:
                v = sum((s[t] * (2 * (z[t] == 0) + (z[t] == 2)) for t in range(k) if t != skip), Q(0))
                A[idx[c][z]][idx[c][z]] = v
            return A

        def E(c):
            A = Cdiag(c)
            if c >= 1:
                JJ = mat_mul(Jp(c - 1), T(Jp(c - 1)))
                for i in range(len(A)):
                    for j in range(len(A)):
                        A[i][j] += JJ[i][j] / 2
            return A

        def sigma_minus(c, ell, r):  # H_c -> H_{c-1}
            A = zeros(len(H[c - 1]), len(H[c]))
            for z in H[c]:
                if z[r] == ell:
                    z2 = z[:r] + (0,) + z[r + 1:]
                    A[idx[c - 1][z2]][idx[c][z]] = Q(1)
            return A

        def Top(c, ell, r):
            A = mat_mul(Jp(c - 1), sigma_minus(c, ell, r))
            A = [[x / 2 for x in row] for row in A]
            for z in H[c]:
                i = idx[c][z]
                if ell == 1:
                    A[i][i] += w[r] * (2 * (z[r] == 0) + (z[r] == 2))
                elif z[r] == 2:
                    z2 = z[:r] + (1,) + z[r + 1:]
                    A[idx[c][z2]][i] += -w[r]
            return A

        def Delta(c, ell, r):
            A = zeros(len(H[c]), len(H[c]))
            for z in H[c]:
                i = idx[c][z]
                nr = 1 if z[r] != 0 else 0
                A[i][i] = Q(2 - nr) - (Q(1, 2) if z[r] == ell else 0)
            return A

        for j in range(1, k + 1):
            Ej = E(j)
            Eji = inv(Ej)
            Z = k - j
            for ell in (1, 2):
                for r in range(k):
                    Ta = Top(j, ell, r)
                    Da = Delta(j, ell, r)
                    # (i) Lemma B diagonal row block against ylib
                    Ma = mat_mul(mat_mul(T(Ta), Eji), Ta)
                    a = tuple(ell if t == r else 0 for t in range(k))
                    for yi, y1 in enumerate(H[j]):
                        for y2 in H[j][yi:yi + 6]:
                            if wstate(a, y1) != wstate(a, y2):
                                continue
                            lhs = (Da[idx[j][y1]][idx[j][y2]] - Ma[idx[j][y1]][idx[j][y2]]) / 2 ** Z
                            mu = Q(1)
                            for t in range(k):
                                mu *= w[t] ** (int(a[t] == 1) * 2 + int(y1[t] == 1) + int(y2[t] == 1))
                            rhs = Bp.Y((a, y1), (a, y2)) / mu
                            ok &= lhs == rhs
                            count['i'] += 1
                    # (ii)/(iii) C_a
                    Ca = mat_mul(mat_mul(Ta, inv(Da)), T(Ta))
                    Ca = [[Ej[i][l] - Ca[i][l] for l in range(len(Ej))] for i in range(len(Ej))]
                    if j >= 1:
                        Jx = Jp_excl(j - 1, r)
                        P = zeros(len(H[j - 1]), len(H[j - 1]))
                        for z in H[j - 1]:
                            P[idx[j - 1][z]][idx[j - 1][z]] = Q(1) if z[r] != 0 else Q(0)
                        hop = mat_mul(mat_mul(Jx, P), T(Jx))
                    if ell == 1:
                        pred = Cdiag(j, skip=r)
                        for i in range(len(pred)):
                            for l in range(len(pred)):
                                pred[i][l] += hop[i][l] / 2
                        good = pred == Ca
                        ok &= good
                        count['ii'] += 1
                    else:
                        pred = Cdiag(j)
                        for i in range(len(pred)):
                            for l in range(len(pred)):
                                pred[i][l] += hop[i][l] / 2
                        for z in H[j]:
                            if z[r] == 1:
                                for t in range(k):
                                    if t != r and z[t] == 0:
                                        z2 = list(z); z2[r] = 0; z2[t] = 1; z2 = tuple(z2)
                                        pred[idx[j][z2]][idx[j][z]] += w[r] * w[t]
                                        pred[idx[j][z]][idx[j][z2]] += w[r] * w[t]
                        good = pred == Ca and psd(Ca)
                        ok &= good
                        count['iii'] += 1
                    if not ok:
                        print('FAILED at', k, trial, j, ell, r)
                        break
    print(f'k = {k}: checks so far {count}, all equal: {ok}', flush=True)
print('PASS' if ok else 'FAIL')
