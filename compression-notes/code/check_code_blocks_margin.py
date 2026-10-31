"""Exact separability certificates by the margin method (no exact LP needed).
   For a sector with weights p (exact Fractions) over the d^2 classes (d = 3^k):
   1. LP (floating): lambda >= 0 with sum_v lambda_v u_v = p - delta * 1, u_v uniform (1/d) on a Lagrangian coset.
   2. Round lambda to rationals (>= 0); r = p - sum lambda u computed EXACTLY.
   3. Isotropic lemma: averaging cosets through one class over all Lagrangians gives (1 + d e_i) / (d (d + 1)), so any
      r >= 0 with min(r) * (d^2 + d) >= sum(r) is in the cone. Check exactly."""
import sys, numpy as np
from fractions import Fraction
from scipy.optimize import linprog
from scipy.sparse import csc_matrix
import stab_fast as sf

def certify(name, gens, b, delta=Fraction(1, 1000)):
    import itertools
    stabs, destabs, logx, logz = sf.adapted_basis(gens, b)
    k = len(logx); m = len(stabs); d = 3 ** k
    allP = np.array(list(itertools.product(range(3), repeat=2 * b)))
    wt = ((allP[:, :b] != 0) | (allP[:, b:] != 0)).sum(1)
    Jb = sf.J(b)
    alpha = (allP @ Jb @ np.array(stabs).T) % 3
    bx = (allP @ Jb @ np.array(logz).T) % 3; bz = (-(allP @ Jb @ np.array(logx).T)) % 3
    beta = np.concatenate([bx, bz], 1)
    acode = alpha @ (3 ** np.arange(m)[::-1]); bcode = beta @ (3 ** np.arange(2 * k)[::-1])
    val = [Fraction(-1, 2) ** int(b - w) for w in range(b + 1)]
    p = {}
    for a, bb, w in zip(acode, bcode, wt): p[(a, bb)] = p.get((a, bb), Fraction(0)) + val[w]
    p = {key: v / 3 ** m for key, v in p.items()}
    lags, vecs = sf.lagrangians(k)
    pw = 3 ** np.arange(2 * k)[::-1]
    cols = sorted({tuple(sorted(int(c) for c in ((vecs[L] + vecs[t]) % 3) @ pw)) for L in lags for t in range(d * d)})
    rows, cc = [], []
    for j, col in enumerate(cols): rows.extend(col); cc.extend([j] * len(col))
    A = csc_matrix((np.full(len(rows), 1.0 / d), (rows, cc)), shape=(d * d, len(cols)))
    allok = True
    for a in sorted(set(acode)):
        pv = [p[(a, t)] for t in range(d * d)]
        target = np.array([float(x - delta) for x in pv])
        res = linprog(np.zeros(len(cols)), A_eq=A, b_eq=target, bounds=(0, None), method="highs-ds")
        if res.status != 0:
            print(f"  sector {a}: LP with margin {delta} infeasible"); allok = False; continue
        lam = {j: Fraction(float(x)).limit_denominator(10 ** 9) for j, x in enumerate(res.x) if x > 1e-13}
        r = list(pv)
        for j, l in lam.items():
            for c in cols[j]: r[c] -= l / d
        ok = min(r) >= 0 and min(r) * (d * d + d) >= sum(r)
        allok &= ok
        print(f"  sector {a}: {len(lam)} cosets; exact remainder min {float(min(r)):.3e}, needed >= "
              f"{float(sum(r) / (d * d + d)):.3e}: {'VALID' if ok else 'FAILED'}")
    print(f"{name}: exact certificates (margin method) in every sector: {allok}")
    return allok

if __name__ == "__main__":
    def Xall(n): v = np.zeros(2 * n, int); v[:n] = 1; return v
    def Zc(n, c): v = np.zeros(2 * n, int); v[n:] = c; return v
    certify("[[4,2,2]]", [Xall(4), Zc(4, [1, 2, 1, 2])], 4)
    certify("[[5,3,2]]", [Xall(5), Zc(5, [1, 1, 1, 1, 2])], 5)
    g = [np.array([1, 0, 0, 2, 0, 0, 1, 2, 0, 0])]
    for _ in range(3):
        g.append(np.concatenate([np.roll(g[-1][:5], 1), np.roll(g[-1][5:], 1)]))
    certify("[[5,1,3]]", g, 5)
