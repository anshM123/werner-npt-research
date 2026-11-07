"""Coset eigenvalues of the cyclic code at the sizes n = 14, ..., 17 and 27, 28, 29 (README, Theorem B).

Ring matrices over R = Z[t]/(t^3 - 1) are stored through their regular representation: a 9 x 9 matrix A over R becomes
the 27 x 27 integer matrix Ah[3i+p, 3j+q] = [t^((q-p) mod 3)] A_ij, so products of ring matrices are products of these
integer matrices and [t^a] Tr A = sum_i Ah[3i, 3i+a].  T_b has (T_b)_{(x,y),(y,z)} = f(y, x+z+b) t^y with f = -1 if
y = 0 and x + z + b = 0 (mod 3), else 2.  Z_{a,d} = [t^a] Tr(T_{d_1} ... T_{d_n}).

Modes:
  python gap_certificates.py control       least Z over all words for n = 4..9 by direct products and by the split
                                           method, and for n = 10..13 by the split method against the values printed
                                           by transfer_certificates.js
  python gap_certificates.py split n L1    every word of length n = L1 + L2 (all 3^n): Z exact through float64
                                           matrix products (every integer involved stays below 2^53, checked);
                                           PASS if every Z_{a,d} > 0
  python gap_certificates.py rows L        row bounds m_i = min, M_i = max of the coefficients of row i over all
                                           products of length L
  python gap_certificates.py block L1 L2   with the row bounds of length L1: for every word Q of length L2 checks
                                           beta(Q) = sum_i (m_i C+_i(Q) - M_i C-_i(Q)) > 0, which gives Z_{a,d} >= beta > 0
                                           for every word of length L1 + L2 (exact int64 and Python integers)
"""
import sys, time
import numpy as np

def T(b):
    R = np.zeros((27, 27), dtype=np.int64)
    for x in range(3):
        for y in range(3):
            for z in range(3):
                c = -1 if (y == 0 and (x + z + b) % 3 == 0) else 2
                i, j = 3 * x + y, 3 * y + z
                for p in range(3):
                    R[3 * i + p, 3 * j + (p + y) % 3] += c
    return R

TB = [T(b) for b in range(3)]


def products(L):
    """All 3^L products, word index = base-3 number with the first letter most significant; shape (3^L, 27, 27)."""
    P = np.eye(27, dtype=np.int64)[None]
    for _ in range(L):
        P = np.concatenate([P @ TB[b] for b in range(3)], axis=0).reshape(3, -1, 27, 27).transpose(1, 0, 2, 3).reshape(-1, 27, 27)
    return P


def zvals(P):
    """Z[w, a] = [t^a] Tr P_w."""
    idx = np.arange(9) * 3
    return np.stack([P[:, idx, idx + a].sum(axis=1) for a in range(3)], axis=1)


def split_min(n, L1, chunk=2000, report=True):
    L2 = n - L1
    P = products(L1).astype(np.float64)
    Q = products(L2).astype(np.float64)
    idx = np.arange(9) * 3
    Pv = P[:, idx, :].reshape(P.shape[0], 9 * 27)                 # P[3i, k] over (i, k)
    zmin, bad, smax = None, 0, 0.0
    for a in range(3):
        Qv = Q[:, :, idx + a].transpose(0, 2, 1).reshape(Q.shape[0], 9 * 27)   # Q[k, 3i+a] over (i, k)
        aP, aQ = np.abs(Pv), np.abs(Qv)
        for s in range(0, Pv.shape[0], chunk):
            G = Pv[s:s + chunk] @ Qv.T
            S = aP[s:s + chunk] @ aQ.T
            smax = max(smax, float(S.max()))
            m = float(G.min())
            zmin = m if zmin is None else min(zmin, m)
            bad += int((G <= 0).sum())
    exact = smax < 2.0 ** 53
    if report:
        print(f'n = {n} (split {L1} + {L2}): {3 ** n} words x 3 coefficients; least Z = {zmin:.0f}; nonpositive: {bad}; '
              f'largest sum of |terms| {smax:.3e} {"< 2^53, so every value is exact" if exact else ">= 2^53: NOT exact"}')
    return zmin, bad, exact


def control():
    ok = True
    for n in range(4, 10):
        P = products(n)
        z = zvals(P)
        direct = int(z.min())
        sm, bad, exact = split_min(n, n // 2, report=False)
        same = exact and int(round(sm)) == direct
        ok &= same and bad == 0
        print(f'n = {n}: least Z over all {3 ** n} words: direct {direct}, split {sm:.0f}, '
              f'eigenvalue {direct}/{3 ** (n - 1) * 2 ** n}{"" if same else "  MISMATCH"}')
    # n = 10..13 by the split method against the minima printed by transfer_certificates.js
    published = {10: 2717712, 11: 13049343, 12: 68934240, 13: 332300799}
    for n, val in published.items():
        sm, bad, exact = split_min(n, (n + 1) // 2, report=False)
        same = exact and bad == 0 and int(round(sm)) == val
        ok &= same
        print(f'n = {n}: least Z over all {3 ** n} words by the split method {sm:.0f}; transfer_certificates.js: {val}'
              f'{"" if same else "  MISMATCH"}')
    print('PASS' if ok else 'FAIL')


def rows(L):
    """Row bounds over all products of length L, computed in blocks (A of length L//2 times all B)."""
    L1 = L // 2
    A = products(L1)
    B = products(L - L1)
    idx = np.arange(9) * 3
    mn = np.full(9, np.iinfo(np.int64).max, dtype=np.int64)
    mx = np.full(9, np.iinfo(np.int64).min, dtype=np.int64)
    for w in range(A.shape[0]):
        Pw = A[w][None] @ B                                         # (3^(L-L1), 27, 27)
        rowsv = Pw[:, idx, :]                                       # (., 9, 27): row i, all (j, y)
        mn = np.minimum(mn, rowsv.min(axis=(0, 2)))
        mx = np.maximum(mx, rowsv.max(axis=(0, 2)))
    print(f'rows L = {L}: m_i = {mn.tolist()}')
    print(f'            M_i = {mx.tolist()}')
    print(f'            global min {mn.min()}, global max {mx.max()}, all m_i > 0: {bool((mn > 0).all())}')
    return mn, mx


def block(L1, L2):
    mn_, mx_ = rows(L1)
    mn, mx = [int(v) for v in mn_], [int(v) for v in mx_]
    assert all(v > 0 for v in mn)
    h = L2 // 2
    A = products(h)
    B = products(L2 - h)
    idx = np.arange(9) * 3
    worst, nfail, count, bmin = None, 0, 0, None
    t0 = time.time()
    for w in range(A.shape[0]):
        Qw = A[w][None] @ B                                         # products of length L2 starting with word w
        Qr = Qw[:, idx, :]                                          # rows 3j: Q[j, :, :] as (., 9 (j), 27 (3i+z))
        ent = Qr.reshape(Qr.shape[0], 9, 9, 3)                      # (., j, i, z)
        Cp = np.where(ent > 0, ent, 0).sum(axis=(1, 3))            # (., i)
        Cm = np.where(ent < 0, -ent, 0).sum(axis=(1, 3))
        # beta = sum_i (m_i C+_i - M_i C-_i), exact: float screen then Python integers for the closest cases
        betaf = (Cp * np.array(mn, dtype=np.float64)).sum(axis=1) - (Cm * np.array(mx, dtype=np.float64)).sum(axis=1)
        scale = (Cp * np.array(mn, dtype=np.float64)).sum(axis=1) + (Cm * np.array(mx, dtype=np.float64)).sum(axis=1)
        close = np.nonzero(betaf <= 1e-9 * scale)[0]
        for q in close:
            be = sum(mn[i] * int(Cp[q, i]) - mx[i] * int(Cm[q, i]) for i in range(9))
            if be <= 0:
                nfail += 1
        r = float((betaf / scale).min())
        worst = r if worst is None else min(worst, r)
        q = int(np.argmin(betaf))
        bq = sum(mn[i] * int(Cp[q, i]) - mx[i] * int(Cm[q, i]) for i in range(9))
        if float(bq) <= 1.0001 * float(betaf.min()):
            bmin = bq if bmin is None else min(bmin, bq)
        else:
            raise RuntimeError('float screen disagrees with the exact value')
        count += Qw.shape[0]
    n = L1 + L2
    print(f'block {L1} + {L2}: {count} tails Q (all 3^{L2}); words of length {n} covered; '
          f'least beta/scale {worst:.4f}; failures {nfail}; least beta {bmin} (so every coset eigenvalue is at least '
          f'{bmin}/{3 ** (n - 1) * 2 ** n} = {bmin / (3 ** (n - 1) * 2 ** n):.4e}); {time.time() - t0:.0f} s')
    print('PASS' if nfail == 0 and count == 3 ** L2 else 'FAIL')


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'control':
        control()
    elif mode == 'split':
        n, L1 = int(sys.argv[2]), int(sys.argv[3])
        zmin, bad, exact = split_min(n, L1)
        print('PASS' if bad == 0 and exact and zmin > 0 else 'FAIL')
    elif mode == 'rows':
        L = int(sys.argv[2])
        mn, mx = rows(L)
        good = bool((mn > 0).all())
        if L == 12:      # the global minimum and maximum printed by tail_certificates.js
            good &= int(mn.min()) == 365472 and int(mx.max()) == 18258048
        print('PASS' if good else 'FAIL')
    elif mode == 'block':
        block(int(sys.argv[2]), int(sys.argv[3]))
