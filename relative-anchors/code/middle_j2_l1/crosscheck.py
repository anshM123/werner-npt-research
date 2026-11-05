"""A second check of Theorem L (the letter-1 block of (R_inf)_2 at two weight values, 5 <= k <= 10), written
separately from the certificate.

Uses only the builders j1lib.py (Lemma B) and ylib.py (from the definitions), r2lib.py (states and kernel vectors),
python-flint and numpy.  None of l2block, young, fixedred, ufield, cert_fixed or saved is used.

  python crosscheck.py N             negative controls of the positive-definiteness test pd_cert (below)
  python crosscheck.py A KMIN KMAX   end to end from j1lib, every size pair (u0, u1), u0 >= u1, u0 + u1 = k, at a
                                     fresh two-value point (levels at random positions, fresh weight ratio, some with
                                     an overall scale): M v = 0 exactly for the C(k,2) vectors p (x) |1@X>, and the
                                     principal submatrix without the states (1@0, X) is positive definite (pd_cert).
                                     Then M >= 0 with kernel exactly span{p (x) |1@X>}.
  python crosscheck.py B             ylib at k = 5, 6 (non-square weight ratios q), same test in the rescaled frame;
                                     and j1lib (rescaled by mu) = 2^(k-2) ylib entry by entry at one point each.
  python crosscheck.py C DIR         the saved certified forms DIR/cert_k*_*_*.pkl: symmetric, stored irrep dimensions
                                     = hook-length dimensions, sum of dim * size = (k - 1) C(k, 2) (= dim of the gauge
                                     complement W), and positive definite (exact elimination) at 31 values of q in
                                     [10^-8, 10^8].

pd_cert (exact): A = a rational symmetric matrix.  With D the lcm of the denominators, A_int = D A (integers).  Z is an
integer lower-triangular matrix with nonzero diagonal (a per-row dyadic rounding of the inverse Cholesky factor
computed in floating point), so B = Z A_int Z^T is congruent to A_int.  If diag(2^-e) B diag(2^-e) is strictly
diagonally dominant with positive diagonal (checked in exact integer arithmetic), it is positive definite (Gershgorin),
hence so is A.  Floating point only chooses Z."""
import sys, time, random, math, os, pickle, itertools
from fractions import Fraction as Fr
from math import comb
import numpy as np
import flint

T0 = time.time()


def log(*a):
    print(*a, flush=True)


def p2(e):
    return Fr(2) ** e


def pd_cert(A, rounds=3):
    """A: symmetric list of lists of Fractions.  Returns (proved, info)."""
    n = len(A)
    D = 1
    for row in A:
        for x in row:
            D = D * x.denominator // math.gcd(D, x.denominator)
    cur = flint.fmpz_mat(n, n, [int(x * D) for row in A for x in row])
    f = [0] * n            # the current matrix is diag(2^-f) cur diag(2^-f) / D
    for rnd in range(rounds):
        G = max(f)
        # exact Gershgorin test on the current matrix
        ok = True
        worst = None
        for i in range(n):
            di = int(cur[i, i]) * (1 << (G - f[i]))
            off = sum(abs(int(cur[i, j])) * (1 << (G - f[j])) for j in range(n) if j != i)
            if di <= off:
                ok = False
            m = Fr(di - off, di) if di != 0 else Fr(-1)
            worst = m if worst is None or m < worst else worst
        if ok:
            return True, f'round {rnd}: strictly diagonally dominant, min relative margin {float(worst):.3g}'
        # one congruence step
        Cf = np.empty((n, n))
        for i in range(n):
            for j in range(n):
                Cf[i, j] = math.ldexp(float(Fr(int(cur[i, j]), D)), -f[i] - f[j])
        try:
            L = np.linalg.cholesky(Cf)
        except np.linalg.LinAlgError:
            ev = np.linalg.eigvalsh(Cf)
            return False, f'round {rnd}: floating-point Cholesky failed; smallest eigenvalue {ev[0]:.3g}, largest {ev[-1]:.3g}'
        Zf = np.linalg.inv(L)
        e = []
        rows = []
        for i in range(n):
            m = float(np.max(np.abs(Zf[i, :i + 1])))
            ei = 50 - math.frexp(m)[1]
            e.append(ei)
            row = [int(round(math.ldexp(float(Zf[i, j]), ei))) if j <= i else 0 for j in range(n)]
            assert row[i] != 0
            rows.append(row)
        # Z_r = diag(2^-e) Zint approximates the inverse factor of diag(2^-f) cur diag(2^-f); fold diag(2^-f) into Z
        F = max(f)
        Zp = flint.fmpz_mat(n, n, [rows[i][j] * (1 << (F - f[j])) for i in range(n) for j in range(n)])
        cur = Zp * cur * Zp.transpose()
        f = [e[i] + F for i in range(n)]
    return False, 'no certificate after the given rounds'


# ------------------------------------------------------------------------------------------------------------- (A)
def part_A(kmin, kmax):
    from j1lib import J1
    from r2lib import states_L1, kernel_vectors
    rng = random.Random(20261105)
    ratios = [Fr(3, 4), Fr(5, 3), Fr(2, 5), Fr(7, 3), Fr(1, 6), Fr(9, 2), Fr(4, 7), Fr(11, 5), Fr(1, 12), Fr(8),
              Fr(6, 5), Fr(3, 10)]
    scales = [Fr(1), Fr(1), Fr(2, 3), Fr(1), Fr(5, 4)]
    allok = True
    cnt = 0
    for k in range(5, kmax + 1):
        for u1 in range(0, k // 2 + 1):
            u0 = k - u1
            rho = ratios[cnt % len(ratios)]
            c0 = scales[cnt % len(scales)]
            cnt += 1
            lev = [0] * u0 + [1] * u1
            rng.shuffle(lev)
            if k < kmin:
                continue
            w = [c0 if l == 0 else c0 * rho for l in lev]
            if u1 == 0:
                w = [c0] * k
            t0 = time.time()
            st = states_L1(k)
            M = J1(w).matrix(st)
            n = len(st)
            cols, vecs = kernel_vectors(st, w)
            Mq = flint.fmpq_mat(n, n, [flint.fmpq(x.numerator, x.denominator) for row in M for x in row])
            V = flint.fmpq_mat(n, len(vecs), [flint.fmpq(vecs[c][i].numerator, vecs[c][i].denominator)
                                               if isinstance(vecs[c][i], Fr) else flint.fmpq(vecs[c][i])
                                               for i in range(n) for c in range(len(vecs))])
            MV = Mq * V
            kern = all(MV[i, c] == 0 for i in range(n) for c in range(len(vecs)))
            ngroups = len(cols) == comb(k, 2)
            sym = all(M[i][j] == M[j][i] for i in range(n) for j in range(i))
            drop = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 1 and a[1] == 0]
            assert len(drop) == comb(k, 2)
            keep = [i for i in range(n) if i not in set(drop)]
            A = [[M[i][j] for j in keep] for i in keep]
            proved, info = pd_cert(A)
            good = kern and ngroups and sym and proved
            allok &= good
            q = (rho * rho) if u1 else Fr(1)
            log(f'(A) k={k} (u0,u1)=({u0},{u1}) levels {"".join(map(str, lev))} w0={c0} q={q}: n={n}, symmetric {sym}, '
                f'{len(cols)} kernel vectors M v = 0: {kern}; gauge-fixed size {len(A)} PD: {proved} ({info}) '
                f'[{time.time() - t0:.1f}s]')
    return allok


# ------------------------------------------------------------------------------------------------------------- (B)
def part_B():
    from ylib import Builder, point_ring
    from j1lib import J1, to_ylib_state, mu
    from r2lib import states_L1
    rng = random.Random(7)
    allok = True

    def ylib_block(s):
        k = len(s)
        Bd = Builder(point_ring([flint.fmpq(x.numerator, x.denominator) for x in s]))
        st = states_L1(k)
        ys = [to_ylib_state(x, k) for x in st]
        n = len(st)
        Y = [[None] * n for _ in range(n)]
        for i in range(n):
            for j in range(i, n):
                v = Bd.Y(ys[i], ys[j])
                Y[i][j] = Y[j][i] = Fr(int(v.p), int(v.q)) * 2 ** (k - 2)
        return st, Y

    for (u0, u1, qv) in [(4, 1, Fr(5, 2)), (3, 2, Fr(3, 7)), (5, 1, Fr(2, 9)), (4, 2, Fr(7, 11)), (3, 3, Fr(13, 3))]:
        k = u0 + u1
        lev = [0] * u0 + [1] * u1
        rng.shuffle(lev)
        s = [Fr(1) if l == 0 else qv for l in lev]
        t0 = time.time()
        st, Y = ylib_block(s)
        n = len(st)
        groups = {}
        for i, (a, y) in enumerate(st):
            if a is not None and a[0] == 1:
                groups.setdefault(y, []).append(i)
        kern = len(groups) == comb(k, 2) and all(sum(Y[i][j] for j in g) == 0 for g in groups.values() for i in range(n))
        drop = set(i for i, (a, y) in enumerate(st) if a is not None and a[0] == 1 and a[1] == 0)
        keep = [i for i in range(n) if i not in drop]
        proved, info = pd_cert([[Y[i][j] for j in keep] for i in keep])
        good = kern and proved
        allok &= good
        log(f'(B) ylib k={k} (u0,u1)=({u0},{u1}) levels {"".join(map(str, lev))} q={qv}: column-group indicators in '
            f'the kernel: {kern}; gauge-fixed PD: {proved} ({info}) [{time.time() - t0:.1f}s]')
    # j1lib (rescaled by mu) against ylib, entry by entry, at one point each for k = 5, 6
    for w in ([Fr(1), Fr(3, 2), Fr(1), Fr(1), Fr(3, 2)], [Fr(2, 3), Fr(2, 3), Fr(5, 4), Fr(2, 3), Fr(5, 4), Fr(2, 3)]):
        k = len(w)
        t0 = time.time()
        st, Y = ylib_block([x * x for x in w])
        M = J1(w).matrix(st)
        mus = [mu(x, w) for x in st]
        eq = all(M[i][j] * mus[i] * mus[j] == Y[i][j] for i in range(len(st)) for j in range(len(st)))
        allok &= eq
        log(f'(B) k={k} w={[str(x) for x in w]}: j1lib rescaled by mu = 2^(k-2) ylib entry by entry: {eq} '
            f'[{time.time() - t0:.1f}s]')
    return allok


# ------------------------------------------------------------------------------------------------------------- (C)
def hook(lam):
    n = sum(lam)
    if n == 0:
        return 1
    conj = [sum(1 for x in lam if x > j) for j in range(lam[0])]
    h = 1
    for i, li in enumerate(lam):
        for j in range(li):
            h *= (li - j - 1) + (conj[j] - i - 1) + 1
    return math.factorial(n) // h


def ldl_pos(A):
    A = [list(r) for r in A]
    n = len(A)
    for c in range(n):
        p = A[c][c]
        if p <= 0:
            return False
        for i in range(c + 1, n):
            if A[i][c] != 0:
                fct = A[i][c] / p
                for j in range(c + 1, n):
                    A[i][j] -= fct * A[c][j]
    return True


def part_C(d):
    QS = [Fr(1, 10 ** 8), Fr(1, 10 ** 6), Fr(1, 10 ** 4), Fr(1, 1000), Fr(1, 100), Fr(1, 30), Fr(1, 10), Fr(1, 5),
          Fr(1, 3), Fr(7, 13), Fr(1, 2), Fr(2, 3), Fr(17, 19), Fr(99, 100), Fr(1), Fr(101, 100), Fr(19, 17), Fr(3, 2),
          Fr(2), Fr(13, 7), Fr(3), Fr(5), Fr(10), Fr(30), Fr(100), Fr(1000), Fr(10 ** 4), Fr(10 ** 6), Fr(10 ** 8),
          Fr(355, 113), Fr(1, 7)]
    assert len(QS) == 31
    allok = True
    files = sorted(f for f in os.listdir(d) if f.startswith('cert_k') and f.endswith('.pkl'))
    seen = set()
    for fn in files:
        k, u0, u1 = (int(x) for x in fn[len('cert_k'):-len('.pkl')].split('_'))
        seen.add((k, u0, u1))
        with open(os.path.join(d, fn), 'rb') as fh:
            dat = pickle.load(fh)
        good = (dat['u0'], dat['u1'], dat['b']) == (u0, u1, '1') and u0 + u1 == k
        tot = 0
        for lam, dim, rows in dat['parts']:
            l0, l1 = lam
            good &= dim == hook(tuple(l0)) * hook(tuple(l1))
            tot += dim * len(rows)
            n = len(rows)
            S = [[(tuple(Fr(c) for c in nc), tuple(Fr(c) for c in dc)) for (nc, dc) in row] for row in rows]
            good &= all(S[i][j] == S[j][i] for i in range(n) for j in range(n))
            if n == 0:
                continue

            def ev(pair, x):
                nc, dc = pair
                num = sum(c * x ** i for i, c in enumerate(nc))
                den = sum(c * x ** i for i, c in enumerate(dc))
                return num / den

            for x in (QS if u1 > 0 else [Fr(1)]):
                if not ldl_pos([[ev(S[i][j], x) for j in range(n)] for i in range(n)]):
                    good = False
                    log(f'   not PD: {fn} part {lam} at q = {x}')
        good &= tot == (k - 1) * comb(k, 2)
        allok &= good
        log(f'(C) {fn}: {len(dat["parts"])} parts, sum dim * size = {tot} = (k-1) C(k,2) = {(k - 1) * comb(k, 2)}; '
            f'symmetric, hook dimensions, PD at {31 if u1 else 1} values of q: {good}')
    want = set((k, k - u1, u1) for k in range(5, 11) for u1 in range(0, k // 2 + 1))
    allok &= want <= seen
    log(f'(C) all 27 size pairs with 5 <= k <= 10 present: {want <= seen}')
    return allok


def part_N():
    """negative controls of pd_cert: it must refuse (1) the full block at k = 5, 6 (singular, C(k,2) zero
    eigenvalues), (2) the gauge-fixed block minus 1.01 lambda_min v v^T/|v|^2 (v a rationalized eigenvector of its
    smallest eigenvalue), which has a negative direction, and (3) the gauge-fixed j1lib block with B -> (1 + 1/1000)B
    at the two-value points of neg_points.py style where it is not PD (searched over a few ratios; at least one must
    be refused)."""
    from j1lib import J1
    from r2lib import states_L1
    allok = True
    for (k, w) in [(5, [Fr(1), Fr(5, 3), Fr(1), Fr(1), Fr(5, 3)]), (6, [Fr(3, 2)] * 3 + [Fr(1)] * 3)]:
        st = states_L1(k)
        M = J1(w).matrix(st)
        n = len(st)
        r1, i1 = pd_cert(M)
        drop = set(i for i, (a, y) in enumerate(st) if a is not None and a[0] == 1 and a[1] == 0)
        keep = [i for i in range(n) if i not in drop]
        A = [[M[i][j] for j in keep] for i in keep]
        Af = np.array([[float(x) for x in row] for row in A])
        ev, U = np.linalg.eigh(Af)
        v = [Fr(x).limit_denominator(10 ** 6) for x in U[:, 0]]
        nv = sum(x * x for x in v)
        lam = Fr(float(ev[0])).limit_denominator(10 ** 12)
        c = Fr(101, 100) * lam / nv
        A2 = [[A[i][j] - c * v[i] * v[j] for j in range(len(A))] for i in range(len(A))]
        vAv = sum(v[i] * sum(A2[i][j] * v[j] for j in range(len(A))) for i in range(len(A)))
        r2, i2 = pd_cert(A2)
        good = (not r1) and (not r2) and vAv < 0
        allok &= good
        log(f'(N) k={k}: full singular block refused: {not r1} ({i1}); perturbed gauge-fixed block (v^T A v = '
            f'{float(vAv):.3g} < 0) refused: {not r2} ({i2})')
    refused = 0
    for rho in (Fr(1, 3), Fr(1, 2), Fr(2), Fr(3), Fr(1, 6), Fr(6)):
        w = [Fr(1), Fr(1), rho, Fr(1), rho]
        st = states_L1(5)
        M = J1(w, bscale=Fr(1001, 1000)).matrix(st)
        n = len(st)
        drop = set(i for i, (a, y) in enumerate(st) if a is not None and a[0] == 1 and a[1] == 0)
        keep = [i for i in range(n) if i not in drop]
        r3, i3 = pd_cert([[M[i][j] for j in keep] for i in keep])
        refused += (not r3)
        log(f'(N) k=5, B scaled by 1001/1000, w ratio {rho}: certified PD {r3} ({i3})')
    allok &= refused >= 1
    log(f'(N) scaled blocks refused at {refused} of 6 ratios (at least one needed)')
    return allok


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'N':
        ok = part_N()
    elif mode == 'A':
        ok = part_A(int(sys.argv[2]), int(sys.argv[3]))
    elif mode == 'B':
        ok = part_B()
    else:
        ok = part_C(sys.argv[2])
    log(f'time {time.time() - T0:.1f}s')
    log('PASS' if ok else 'FAIL')
