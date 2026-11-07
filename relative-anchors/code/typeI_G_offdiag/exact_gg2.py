"""Exact checks (python-flint) at rational weight vectors of the two conditions of Lemma GG' with one constant kappa:
  (O)  ||Rhat_raw,tt'||_op <= kappa sqrt(theta_t theta_t')          for t != t',
  (D)  Ghat_tt - ((1 - theta_t)/theta_t) Shat_free,tt - kappa (1 - theta_t) I  > 0    (strict) for every t,
with R_raw, S_free = R_raw - G as in exact_offdiag.py.  Lemma GG' says (O) and (D) for all t, t' give G > 0.
Also checks the conclusion G > 0 itself exactly (all LDL pivots of G positive).
Usage: python exact_gg2.py seed npoints kappa"""
import sys, random
import flint
from exact_offdiag import raw_split, Nmat, psd_exact, frac
from lemmaq_exact import direct_G, mat

Q = flint.fmpq


def pd_exact(A):
    return psd_exact(A) and A.det() != 0


seed, npts, kap = int(sys.argv[1]), int(sys.argv[2]), frac(sys.argv[3])
rng = random.Random(seed)
fails = {'O': 0, 'D': 0, 'G': 0}
nb = nd = 0
for it in range(npts):
    m = rng.randint(2, 4); u = rng.randint(2, 6); mode = it % 4
    def amp(lo, hi):
        x = 2.0 ** rng.uniform(lo, hi)
        return Q(int(round(x * 64)) or 1, 64)
    if mode == 0:
        wT = [amp(-3, 3) for _ in range(m)]; wU = [amp(-3, 3) for _ in range(u)]
    elif mode == 1:
        wT = [amp(2, 5) for _ in range(m)]; wU = [amp(-3, 0) for _ in range(u)]
    elif mode == 2:
        wT = [amp(-3, 0) for _ in range(m)]; wU = [amp(1, 4) for _ in range(u)]
    else:
        base = amp(1, 3); wT = [base, base * Q(rng.randint(5, 9), 7)] + [amp(-2, 2) for _ in range(m - 2)]
        wU = [amp(-3, 3) for _ in range(u)]
    G, R, info = raw_split(wT, wU)
    Gd, _, _ = direct_G(wT, wU)
    d, s, sig = info['d'], info['s'], info['sig']
    N = Nmat(info); Ninv = N.inv()
    if not pd_exact(Gd):
        fails['G'] += 1
    for t in range(m):
        th = s[t] / sig
        S = R[(t, t)] - G[(t, t)]
        nd += 1
        if not pd_exact(G[(t, t)] - (1 - th) / th * S - kap * (1 - th) * N):
            fails['D'] += 1
        for tp in range(m):
            if tp != t:
                nb += 1
                k2 = kap * kap * s[t] * s[tp] / sig ** 2
                if not psd_exact(k2 * N - R[(t, tp)] * Ninv * R[(tp, t)]):
                    fails['O'] += 1
    print(f'point {it} (m={m}, u={u}, mode {mode}): cumulative failures {fails}', flush=True)
print(f'kappa = {kap}: {nb} off-diagonal blocks, {nd} diagonal blocks at {npts} rational points; exact failures {fails}')
