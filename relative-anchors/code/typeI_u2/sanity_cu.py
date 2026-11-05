"""Float sanity check (no claim rests on it): for random weights, u = 2..8, sigma/|s| over 14 decades:
 - smallest eigenvalue of the (C_u) principal block (cu_form.CU.full_matrix) relative to sigma/(sigma + S);
 - smallest eigenvalue of Phi1 (normalised block formulas) relative to 1 - t;
 - the certified per-coordinate margins min_h (LB_h - (4/5) ROW_h) and max_h COL_h at the actual parameters."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cu_form import CU
from phi1 import blocks_norm
from rowbounds2 import bounds

rng = np.random.default_rng(21)
t0 = time.time()
wC = np.inf; wP = np.inf; wI = np.inf; wII = 0
for u in (2, 3, 4, 6, 8):
    for trial in range(60 if u < 6 else 25):
        spread = rng.choice([0.3, 1, 3, 6])
        w = np.exp(rng.uniform(-spread, spread, u))
        if trial % 4 == 1:
            w[0] *= np.exp(rng.uniform(2, 6))
        S = (w ** 2).sum()
        sig = S * np.exp(rng.uniform(-14, 14) * np.log(10) / 2)
        cu = CU(list(w), sig, F=float)
        M = np.array(cu.full_matrix())
        lam = np.linalg.eigvalsh(M)[0]
        wC = min(wC, lam / (sig / (sig + S)))
        K = sig + S / 2; x = w ** 2 / K; t = sig / K
        ww = np.outer(w, w) / K
        Ma, Mc, Mca = [np.array(Z, float) for Z in blocks_norm(list(x), t, ww.tolist(), F=float)]
        Phi = np.block([[Ma, Mca.T], [Mca, Mc]])
        wP = min(wP, np.linalg.eigvalsh(Phi)[0] / (1 - t))
        rho = x / (1 + x); g = rho.sum(); eta2 = x / (1 + x) ** 2; N = eta2.sum(); Sx = x.sum()
        for h in range(u):
            Sp = Sx - x[h]; rSp = Sp / (1 + Sp); gp = g - rho[h]
            th = (gp - rSp) / (Sp - rSp) if Sp - rSp > 1e-300 else 0.0
            Dp = Sp - gp; Nub_p = gp - Dp / (1 + Sp)
            nu = (N - eta2[h]) / Nub_p if Nub_p > 0 else 0.0
            bb = bounds(Sx, x[h] / Sx, min(max(th, 0), 1), min(max(nu, 0), 1), one=1.0)
            wI = min(wI, bb['LB'] - 0.8 * max(bb['ROWp'], bb['ROWm']))
            wII = max(wII, max(bb['COLp'], bb['COLm']))
print(f"min lam(C_u block)/(sigma/(sigma+S)) = {wC:.4f};  min lam(Phi1)/(1-t) = {wP:.4f};  "
      f"min (LB - 4/5 ROW) = {wI:.4f};  max COL = {wII:.4f}  ({time.time() - t0:.0f} s)")
ok = wC > 0 and wP > 0 and wI > 0 and wII < 0.8
print("PASS" if ok else "FAIL")
