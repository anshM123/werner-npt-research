"""Exact pointwise certificates (python-flint fmpq): at rational weight vectors (w_t, w_U), all LDL^T pivots of the
relaxed form Psi (psi_exact.py) are positive.  By Lemma R1 (NOTES.md) this proves G > 0 at that weight vector, hence
(Theorem C) the type I block of (R_inf)_{k-1} with T = {t} is PSD with kernel exactly span{p (x) f_h : h in U}.

usage: python certify_points.py u_list npts_list [seed]   (npts_list: one count, or one count per u)
Families (cycled): log-uniform; one U-weight tiny; all U-weights shrunk (S << sig); all enlarged (S >> sig);
one U-weight dominant; geometric spread; two tiny U-weights; sig tiny with spread U-weights."""
import sys, time, random
import flint
from psi_exact import PsiExact, ldl_pivots

Q = flint.fmpq


def leading_minors_bareiss(Mx):
    """All leading principal minors of the rational symmetric matrix Mx, times den^k (den = lcm of the
    denominators, so the signs are those of the minors), by fraction-free Bareiss elimination over the integers."""
    n = Mx.nrows()
    L = flint.fmpz(1)
    for i in range(n):
        for j in range(n):
            q = Mx[i, j].q
            L = L * q // L.gcd(q)
    A = [[(Mx[i, j] * L).p for j in range(n)] for i in range(n)]
    mins = []
    prev = flint.fmpz(1)
    for k in range(n):
        piv = A[k][k]
        mins.append(piv)
        if piv == 0:
            return mins, L
        for i in range(k + 1, n):
            Arow = A[i]
            arow_k = Arow[k]
            Ak = A[k]
            for j in range(k + 1, n):
                Arow[j] = (piv * Arow[j] - arow_k * Ak[j]) // prev
        prev = piv
    return mins, L


def ratw(x, den=10 ** 3):
    """a positive rational close to x > 0 (x may be very small or large)."""
    import math
    if x <= 0:
        raise ValueError
    e = math.floor(math.log10(x))
    m = x / 10 ** e
    q = Q(int(round(m * den)), den)
    return q * (Q(10) ** e if e >= 0 else Q(1) / (Q(10) ** (-e)))


def family(u, i, rng):
    import math
    fam = i % 8
    wt = 1.0
    if fam == 0:
        wU = [math.exp(rng.uniform(-2, 2)) for _ in range(u)]
        wt = math.exp(rng.uniform(-2, 2))
    elif fam == 1:
        wU = [math.exp(rng.uniform(-1, 1)) for _ in range(u)]
        wU[0] = 10 ** rng.uniform(-4, -2)
    elif fam == 2:
        f = 10 ** rng.uniform(-4, -1)
        wU = [f * math.exp(rng.uniform(-1, 1)) for _ in range(u)]
    elif fam == 3:
        f = 10 ** rng.uniform(1, 4)
        wU = [f * math.exp(rng.uniform(-1, 1)) for _ in range(u)]
    elif fam == 4:
        wU = [math.exp(rng.uniform(-1, 1)) for _ in range(u)]
        wU[0] = 10 ** rng.uniform(1, 3)
    elif fam == 5:
        r = 10 ** rng.uniform(-1.5 / max(u - 1, 1), 1.5 / max(u - 1, 1))
        wU = [r ** j for j in range(u)]
    elif fam == 6:
        wU = [math.exp(rng.uniform(-1, 1)) for _ in range(u)]
        wU[0] = 10 ** rng.uniform(-4, -2)
        wU[1] = 10 ** rng.uniform(-4, -2)
    else:
        wU = [math.exp(rng.uniform(-2, 2)) for _ in range(u)]
        wt = 10 ** rng.uniform(-4, -2)
    return ratw(wt), [ratw(x) for x in wU]


if __name__ == '__main__':
    us = [int(x) for x in sys.argv[1].split(',')]
    nl = [int(x) for x in sys.argv[2].split(',')]
    nl = nl * len(us) if len(nl) == 1 else nl
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    rng = random.Random(seed)
    allok = True
    for u, npts in zip(us, nl):
        t0 = time.time()
        okc = 0
        minr = None
        for i in range(npts):
            wt, wU = family(u, i, rng)
            P = PsiExact(wt, wU)
            Mx, iy = P.matrix()
            mins, Lden = leading_minors_bareiss(Mx)
            ok = len(mins) == Mx.nrows() and all(x > 0 for x in mins)
            allok &= ok
            okc += ok
            # smallest pivot of the letter-2 part relative to its diagonal entry (a crude margin indicator)
            import math
            lg = lambda z: math.log2(abs(int(z))) if z != 0 else -1e9
            # pivot_k = mins[k] / (mins[k-1] * Lden); ratio to the diagonal entry Mx[k, k], in log2
            r = min(lg(mins[iy + j]) - lg(mins[iy + j - 1]) - lg(Lden) - lg(Mx[iy + j, iy + j].p) + lg(Mx[iy + j, iy + j].q)
                    for j in range(u + 1)) if ok else float('nan')
            minr = r if minr is None else min(minr, r)
            if not ok:
                print(f'  FAILED at u = {u}, family {i % 8}: w_t = {wt}, w_U = {wU}', flush=True)
        print(f'u = {u}: {okc}/{npts} weight vectors with all pivots of Psi positive; log2 of the smallest letter-2 pivot ratio '
              f'{minr:.1f}; {time.time() - t0:.1f} s', flush=True)
    print('PASS' if allok else 'FAIL')
