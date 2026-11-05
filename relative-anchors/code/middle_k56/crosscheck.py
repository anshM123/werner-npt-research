"""A second check of Theorems 5_2 and 6_2 (all of (R_inf) at k = 5, 6 for weights with at most two distinct values),
written separately from the certificates.

Uses only the builder j1lib.py (Lemma B, every j, units 2^(k-j)), python-flint and numpy.  None of k5lib, cert5,
cert_sym, symred, signtest, saved or the check scripts of this folder is used; the torus blocks, their known kernel
vectors and the gauge are set up here from the definitions.

  python crosscheck.py N                    negative controls of pd_cert
  python crosscheck.py A K NPTS SEED [JS]   all of (R_inf)_j at k = K (j in JS, default -1..K) at NPTS random two-value
                                            points: the states of G1_j + G0_(j+1) are grouped by the number of letters 2
                                            at each copy (row and column together), which the torus symmetry conserves
                                            (a random sample of entries between different groups is checked to be 0);
                                            in every group the known kernel vectors are annihilated exactly and the
                                            principal submatrix without a set of pivot coordinates of them is positive
                                            definite (pd_cert).  Then the group is >= 0 with kernel exactly their span.
  python crosscheck.py C DIR                the blocks saved by the certificates (DIR/cert_k*.pkl, the rescaled block over
                                            Q(q)): at 31 values of q in [1e-8, 1e8], the known kernel vectors (rescaled
                                            frame, set up here from the saved states) are annihilated and the block
                                            without pivot coordinates of them is positive definite (exact LDL, flint).

Known kernel vectors of a group (unrescaled frame): p (x) y = sum_r w_r (1@r, y) for every column y carried by the
letter-1 rows of the group; in the letter-1 group of (R_inf)_(k-2) also vac (x) J^-|1...1> (Proposition K2); in the
letter-1 group of (R_inf)_(k-1) also u (x) J^-|1...1>, u orthogonal to p (Theorem T).

pd_cert (exact): A rational symmetric.  A_int = D A with D the lcm of the denominators.  Z = an integer lower-triangular
matrix with nonzero diagonal (per-row dyadic rounding of the inverse of a floating-point Cholesky factor), so
B = Z A_int Z^T is congruent to A_int; if diag(2^-e) B diag(2^-e) is strictly diagonally dominant with positive
diagonal (exact integer arithmetic), A is positive definite.  Floating point only chooses Z."""
import sys, time, random, math, os, pickle, itertools
from fractions import Fraction as Fr
from math import comb
import numpy as np
import flint

T0 = time.time()


def log(*a):
    print(*a, flush=True)


# ------------------------------------------------------------------------------------------- exact PD certificate
def pd_cert(A, rounds=3):
    n = len(A)
    if n == 0:
        return True, 'empty'
    D = 1
    for row in A:
        for x in row:
            D = D * x.denominator // math.gcd(D, x.denominator)
    cur = flint.fmpz_mat(n, n, [int(x * D) for row in A for x in row])
    f = [0] * n
    for rnd in range(rounds):
        G = max(f)
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
            return True, f'round {rnd}, margin {float(worst):.3g}'
        Cf = np.empty((n, n))
        for i in range(n):
            for j in range(n):
                Cf[i, j] = math.ldexp(float(Fr(int(cur[i, j]), D)), -f[i] - f[j])
        try:
            L = np.linalg.cholesky(Cf)
        except np.linalg.LinAlgError:
            ev = np.linalg.eigvalsh(Cf)
            return False, f'round {rnd}: floating Cholesky failed; eigenvalues {ev[0]:.3g} .. {ev[-1]:.3g}'
        Zf = np.linalg.inv(L)
        e, rows = [], []
        for i in range(n):
            m = float(np.max(np.abs(Zf[i, :i + 1])))
            ei = 50 - math.frexp(m)[1]
            e.append(ei)
            row = [int(round(math.ldexp(float(Zf[i, j]), ei))) if j <= i else 0 for j in range(n)]
            assert row[i] != 0
            rows.append(row)
        F = max(f)
        Zp = flint.fmpz_mat(n, n, [rows[i][j] * (1 << (F - f[j])) for i in range(n) for j in range(n)])
        cur = Zp * cur * Zp.transpose()
        f = [e[i] + F for i in range(n)]
    return False, 'no certificate after the given rounds'


def ldl_pd_flint(A):
    """exact symmetric elimination over Q (flint fmpq); True iff positive definite"""
    n = len(A)
    M = [[flint.fmpq(x.numerator, x.denominator) for x in row] for row in A]
    for c in range(n):
        p = M[c][c]
        if p <= 0:
            return False
        for i in range(c + 1, n):
            if M[i][c] != 0:
                fct = M[i][c] / p
                Mi, Mc = M[i], M[c]
                for jj in range(c + 1, n):
                    if Mc[jj] != 0:
                        Mi[jj] -= fct * Mc[jj]
    return True


def pivots_of(vecs):
    """pivot coordinates of the row space of vecs (exact RREF); the restriction of span(vecs) to them is bijective"""
    rows = [[Fr(x) for x in v] for v in vecs]
    piv = []
    r = 0
    n = len(rows[0]) if rows else 0
    for c in range(n):
        p = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        pv = rows[r][c]
        rows[r] = [x / pv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                fc = rows[i][c]
                rows[i] = [x - fc * y for x, y in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
        if r == len(rows):
            break
    assert r == len(rows), 'kernel vectors not independent'
    return piv


# ------------------------------------------------------------------------------------------- states and groups
def strings(k, nnz):
    out = []
    for pos in itertools.combinations(range(k), nnz):
        for lets in itertools.product((1, 2), repeat=nnz):
            y = [0] * k
            for p_, l in zip(pos, lets):
                y[p_] = l
            out.append(tuple(y))
    return out


def rinf_states(k, j):
    st = []
    if 0 <= j <= k:
        for y in strings(k, j):
            for l in (1, 2):
                for r in range(k):
                    st.append(((l, r), y))
    if 0 <= j + 1 <= k:
        for y in strings(k, j + 1):
            st.append((None, y))
    return st


def label(state, k):
    a, y = state
    return tuple((1 if (a is not None and a[0] == 2 and a[1] == t) else 0) + (1 if y[t] == 2 else 0) for t in range(k))


def known_kernel(states, k, j, w):
    """known kernel vectors of a group, unrescaled frame (j1lib)"""
    idx = {s: i for i, s in enumerate(states)}
    n = len(states)
    vecs = []
    cols = []
    for (a, y) in states:
        if a is not None and a[0] == 1 and y not in cols:
            cols.append(y)
    for y in cols:
        v = [Fr(0)] * n
        for r in range(k):
            v[idx[((1, r), y)]] = Fr(w[r])
        vecs.append(v)
    lab = label(states[0], k)
    if all(x == 0 for x in lab):
        ones = tuple([1] * k)
        hole = lambda t: tuple(0 if s == t else 1 for s in range(k))
        if j == k - 2:
            v = [Fr(0)] * n
            for t in range(k):
                v[idx[(None, hole(t))]] = Fr(w[t])
            vecs.append(v)
        if j == k - 1:
            for i_ in range(1, k):
                u = [Fr(0)] * k
                u[i_] = Fr(w[0])
                u[0] = -Fr(w[i_])
                v = [Fr(0)] * n
                for r in range(k):
                    for t in range(k):
                        v[idx[((1, r), hole(t))]] += u[r] * w[t]
                vecs.append(v)
    return vecs


def check_group(M, vecs):
    n = len(M)
    kern = all(sum(M[i][c] * v[c] for c in range(n) if v[c] != 0) == 0 for v in vecs for i in range(n))
    piv = set(pivots_of(vecs)) if vecs else set()
    keep = [i for i in range(n) if i not in piv]
    proved, info = pd_cert([[M[i][jj] for jj in keep] for i in keep])
    return kern, proved, info


# ------------------------------------------------------------------------------------------- (A)
def part_A(k, npts, seed, js):
    from j1lib import J1
    rng = random.Random(seed)
    amps = [Fr(1, 3), Fr(2, 5), Fr(3, 4), Fr(5, 6), Fr(4, 3), Fr(3, 2), Fr(7, 4), Fr(5, 2), Fr(3)]
    allok = True
    for pt in range(npts):
        u1 = rng.randint(1, k - 1)
        lev = [0] * (k - u1) + [1] * u1
        rng.shuffle(lev)
        a0, a1 = rng.sample(amps, 2)
        w = [a0 if l == 0 else a1 for l in lev]
        J = J1(w)
        log(f'(A) k={k} point {pt}: levels {"".join(map(str, lev))}, weights {a0}, {a1} (q = {(a1 / a0) ** 2})')
        for j in js:
            t0 = time.time()
            st = rinf_states(k, j)
            groups = {}
            for s in st:
                groups.setdefault(label(s, k), []).append(s)
            # entries between different groups: a random sample must vanish
            labs = list(groups)
            zero_ok = True
            if len(labs) > 1:
                for _ in range(300):
                    g1, g2 = rng.sample(labs, 2)
                    s1, s2 = rng.choice(groups[g1]), rng.choice(groups[g2])
                    if J.entry(s1, s2) != 0:
                        zero_ok = False
            ok_j = zero_ok
            nker = 0
            worst = None
            for lab_, gs in groups.items():
                M = J.matrix(gs)
                vecs = known_kernel(gs, k, j, w)
                kern, proved, info = check_group(M, vecs)
                nker += len(vecs)
                if not (kern and proved):
                    ok_j = False
                    log(f'   FAILED group {lab_} (size {len(gs)}): kernel {kern}, PD {proved} ({info})')
            expected = (len(strings(k, j)) if 0 <= j <= k else 0) + (1 if j == k - 2 else 0) + (k - 1 if j == k - 1 else 0)
            ok_j &= nker == expected
            allok &= ok_j
            log(f'   (R_inf)_{j}: {len(st)} states in {len(groups)} torus groups, cross entries zero (sample) {zero_ok}; '
                f'every group >= 0 with kernel exactly its known vectors: {ok_j}; kernel dimension {nker} '
                f'(= |H_j| + K2 + T = {expected}) [{time.time() - t0:.1f}s]')
    return allok


# ------------------------------------------------------------------------------------------- (C)
QS = [Fr(1, 10 ** 8), Fr(1, 10 ** 6), Fr(1, 10 ** 4), Fr(1, 1000), Fr(1, 100), Fr(1, 30), Fr(1, 10), Fr(1, 5),
      Fr(1, 3), Fr(7, 13), Fr(1, 2), Fr(2, 3), Fr(17, 19), Fr(99, 100), Fr(1), Fr(101, 100), Fr(19, 17), Fr(3, 2),
      Fr(2), Fr(13, 7), Fr(3), Fr(5), Fr(10), Fr(30), Fr(100), Fr(1000), Fr(10 ** 4), Fr(10 ** 6), Fr(10 ** 8),
      Fr(355, 113), Fr(1, 7)]


def ev(cs, x):
    v = Fr(0)
    for c in reversed(cs):
        v = v * x + Fr(c)
    return v


def part_C(d, maxq):
    allok = True
    files = sorted(f for f in os.listdir(d) if f.startswith('cert_k') and f.endswith('.pkl'))
    for fn in files:
        t0 = time.time()
        with open(os.path.join(d, fn), 'rb') as fh:
            dat = pickle.load(fh)
        k, j = dat['k'], dat['j']
        # saved states are in the format of ylib.py: (row string, column string); convert to (None | (letter, copy), y)
        states = []
        for (a, y) in dat['states']:
            a, y = tuple(a), tuple(y)
            nz = [r for r in range(k) if a[r] != 0]
            assert len(nz) <= 1
            states.append(((a[nz[0]], nz[0]) if nz else None, y))
        A, B = dat['config']
        m = len(dat['T'])
        assert sorted(dat['T']) == list(range(m))
        ones = set(range(A)) | set(range(m, m + B))      # configuration (A, B): squared weight 1 there, q elsewhere
        n = len(states)
        # rescaled frame: kernel vectors are column-group indicators, and for the letter-1 block of (R_inf)_(k-2) the
        # vector (s_h) on the states (vac, 1..1 - e_h); s = squared weights of the configuration as functions of q
        cols = []
        for (a, y) in states:
            if a is not None and a[0] == 1 and y not in cols:
                cols.append(y)
        lab0 = label(states[0], k)
        letter1 = all(x == 0 for x in lab0)
        good = True
        qs = QS if len(QS) <= maxq else QS[::max(1, len(QS) // maxq)]
        for x in qs:
            sq = [Fr(1) if t in ones else x for t in range(k)]
            vecs = []
            for y in cols:
                v = [Fr(0)] * n
                for i, (a, yy) in enumerate(states):
                    if a is not None and a[0] == 1 and yy == y:
                        v[i] = Fr(1)
                vecs.append(v)
            if letter1 and j == k - 2:
                v = [Fr(0)] * n
                for i, (a, yy) in enumerate(states):
                    if a is None and sum(1 for c in yy if c == 1) == k - 1:
                        h = yy.index(0)
                        v[i] = sq[h]
                vecs.append(v)
            M = [[ev(nc, x) / ev(dc, x) for (nc, dc) in row] for row in dat['block']]
            kern = all(sum(M[i][c] * v[c] for c in range(n) if v[c] != 0) == 0 for v in vecs for i in range(n))
            piv = set(pivots_of(vecs)) if vecs else set()
            keep = [i for i in range(n) if i not in piv]
            pd = ldl_pd_flint([[M[i][jj] for jj in keep] for i in keep])
            if not (kern and pd):
                good = False
                log(f'   {fn}: q = {x}: kernel {kern}, PD {pd}')
        allok &= good
        log(f'(C) {fn}: {n} states, {len(qs)} values of q: known kernel annihilated and gauge-fixed block PD: {good} '
            f'[{time.time() - t0:.1f}s]')
    return allok


# ------------------------------------------------------------------------------------------- (N)
def part_N():
    from j1lib import J1
    allok = True
    k, j = 5, 3
    w = [Fr(1), Fr(1), Fr(5, 2), Fr(1), Fr(5, 2)]
    st = rinf_states(k, j)
    gs = [s for s in st if all(x == 0 for x in label(s, k))]
    M = J1(w).matrix(gs)
    r1, i1 = pd_cert(M)
    vecs = known_kernel(gs, k, j, w)
    piv = set(pivots_of(vecs))
    keep = [i for i in range(len(gs)) if i not in piv]
    A = [[M[i][jj] for jj in keep] for i in keep]
    Af = np.array([[float(x) for x in row] for row in A])
    evs, U = np.linalg.eigh(Af)
    v = [Fr(x).limit_denominator(10 ** 6) for x in U[:, 0]]
    nv = sum(x * x for x in v)
    c = Fr(101, 100) * Fr(float(evs[0])).limit_denominator(10 ** 12) / nv
    A2 = [[A[i][jj] - c * v[i] * v[jj] for jj in range(len(A))] for i in range(len(A))]
    vAv = sum(v[i] * sum(A2[i][jj] * v[jj] for jj in range(len(A))) for i in range(len(A)))
    r2, i2 = pd_cert(A2)
    # one kernel vector left out of the gauge: the block keeps a zero eigenvalue
    piv3 = set(pivots_of(vecs[:-1]))
    keep3 = [i for i in range(len(gs)) if i not in piv3]
    r3, i3 = pd_cert([[M[i][jj] for jj in keep3] for i in keep3])
    good = (not r1) and (not r2) and (not r3) and vAv < 0
    allok &= good
    log(f'(N) letter-1 block of (R_inf)_3, k = 5: full block (singular) refused {not r1}; perturbed gauge-fixed block '
        f'(v^T A v = {float(vAv):.3g}) refused {not r2}; gauge missing one kernel vector refused {not r3}')
    refused = 0
    for rho in (Fr(1, 10), Fr(1, 4), Fr(4), Fr(10)):
        w = [Fr(1), rho, Fr(1), Fr(1), rho]
        M = J1(w, bscale=Fr(1001, 1000)).matrix(gs)
        vecs = known_kernel(gs, k, j, w)
        piv = set(pivots_of(vecs))
        keep = [i for i in range(len(gs)) if i not in piv]
        r, info = pd_cert([[M[i][jj] for jj in keep] for i in keep])
        refused += (not r)
        log(f'(N) B scaled by 1001/1000, weight ratio {rho}: certified PD {r} ({info})')
    allok &= refused >= 1
    log(f'(N) scaled blocks refused at {refused} of 4 ratios (at least one needed)')
    return allok


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'N':
        ok = part_N()
    elif mode == 'A':
        k, npts, seed = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
        js = [int(x) for x in sys.argv[5].split(',')] if len(sys.argv) > 5 else list(range(-1, k + 1))
        ok = part_A(k, npts, seed, js)
    else:
        ok = part_C(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 31)
    log(f'time {time.time() - T0:.1f}s')
    log('PASS' if ok else 'FAIL')
