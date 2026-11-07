"""Exact check (python-flint) of condition (O),  ||Rhat_raw,tt'||_N <= kappa sqrt(theta_t theta_t'),  on the FORMAL
domain: free measure with rational atoms x_i in [0, 1] and rational masses m_i > 0 (sum 1) unrelated to the
locations, and also on physical measures (m_i = x_i).  The raw pair is the closed form rawcf.raw_pair (checked
exactly against the whole block in check_exact.py for physical measures).  Test: kappa^2 N - Rhat N^{-1} Rhat^T is
positive semidefinite (exact LDL), Rhat = sigma R/(wa wb), N the metric rho^2 + E Z^2 - (E Z)^2/(sigma + 1).
Usage: python exact_formal.py seed npoints kappa"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from exact_offdiag import psd_exact, frac
from rawcf import raw_pair, metric

Q = flint.fmpq


def mat(rows):
    return flint.fmpq_mat(rows)


def check_point(sig, wa, wb, x, m, kap):
    M, st = raw_pair(x, m, wa * wa, wb * wb, sig, wa, wb, F=Q)
    n = len(M)
    Rh = mat([[M[i][j] * sig / (wa * wb) for j in range(n)] for i in range(n)])
    N = mat(metric(x, m, sig, F=Q))
    A = kap * kap * N - Rh * N.inv() * Rh.transpose()
    return psd_exact(A)


if __name__ == '__main__':
    seed, npts, kap = int(sys.argv[1]), int(sys.argv[2]), frac(sys.argv[3])
    rng = random.Random(seed)
    fails = 0
    for it in range(npts):
        n = rng.randint(1, 4)
        formal = it % 2 == 0
        if formal:
            x = [Q(rng.randint(0, 20), 20) for _ in range(n)]
            mw = [Q(rng.randint(1, 30)) for _ in range(n)]
            tot = sum(mw, Q(0)); m = [v / tot for v in mw]
        else:
            mw = [Q(rng.randint(1, 30)) for _ in range(n)]
            tot = sum(mw, Q(0)); m = [v / tot for v in mw]; x = list(m)
        e = rng.randint(-8, 10)
        base = Q(2) ** e
        wa = Q(rng.randint(1, 12), 8); wb = Q(rng.randint(1, 12), 8)
        if it % 5 == 0:
            wb = wa * Q(6, 7)                      # amplitude ratio 7:6
        spect = Q(rng.randint(0, 3), 4) * (wa * wa + wb * wb)
        scale = base / (wa * wa + wb * wb + spect)  # sigma = base; a, b scaled; need wa, wb rational: scale amplitudes
        # choose amplitudes so that sigma = (wa^2 + wb^2 + spect) * c^2 with c = 2^k rational
        k = e // 2
        c = Q(2) ** k
        wa2, wb2 = wa * c, wb * c
        sig = (wa * wa + wb * wb + spect) * c * c
        ok = check_point(sig, wa2, wb2, x, m, kap)
        if not ok:
            fails += 1
        print(f'point {it}: {"formal" if formal else "physical"} n={n} sigma={sig} wa={wa2} wb={wb2} '
              f'x={[str(v) for v in x]} m={[str(v) for v in m]}: {"pass" if ok else "FAIL"}', flush=True)
    print(f'kappa = {kap}: {npts} points, exact failures: {fails}')
