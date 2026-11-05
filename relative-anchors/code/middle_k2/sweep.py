"""Float sweep (no claim rests on it): letter-1 block of (R_inf)_{k-2}, k = 4..7.
For each weight vector:
 - the full block (Y units) is PSD with kernel dimension exactly C(k,2) + 1 (relative tolerance);
 - gap := smallest eigenvalue on the complement of the kernel, compared with the proved upper bound
   ub := min_h s_h/(2(|s| - s_h)) (Corollary S) and with x_min = min s/|s|;
 - hard part: lambda_min of F (units 4Y) on R = orthogonal complement of (S + p (x) H_{k-2}) in G1,
   S = span of the soft vectors Psi(a,h).
Also a direct comparison with the construction from M, V_p, D for k = 4, 5 (a few vectors)."""
import sys, time
import numpy as np
from math import comb
from mid import Mid
from ex1 import kernel_basis
from ex2 import soft_basis
from direct import direct_block

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 2026)
NS = {4: 300, 5: 200, 6: 80, 7: 30}
ok = True
t0 = time.time()
for k in (4, 5, 6, 7):
    worst = dict(neg=0.0, kerdev=0, ratio_ub=np.inf, ratio_x=np.inf, lamR=np.inf, ratio_ub_max=0)
    for trial in range(NS[k]):
        mode = trial % 6
        if mode == 0:
            w = rng.random(k) + 0.2
        elif mode == 1:
            w = np.exp(2 * rng.normal(size=k))
        elif mode == 2:
            w = np.exp(rng.normal(size=k)); w[rng.integers(k)] *= 10 ** rng.uniform(-3, -1)
        elif mode == 3:
            w = np.exp(rng.normal(size=k)); w[:2] *= 10 ** rng.uniform(-3, -1)
        elif mode == 4:
            w = np.exp(rng.normal(size=k)); w[rng.integers(k)] *= 10 ** rng.uniform(1, 2.5)
        else:
            w = np.exp(rng.normal(size=k)); w[:2] *= 10 ** rng.uniform(1, 2.5)
        m = Mid(list(w), float)
        Y = np.array(m.full(), float) / 4          # Y units
        n = m.np; N1 = k * n
        s = w ** 2; S = s.sum(); x = s / S
        # kernel basis K (p (x) H_{k-2} and vac (x) J^-|1..1>)
        K = np.zeros((N1 + k, n + 1)); K[:N1, :n] = kernel_basis(m); K[N1:, n] = w
        Qk, _ = np.linalg.qr(K)
        P = np.eye(N1 + k) - Qk @ Qk.T
        ev = np.linalg.eigvalsh(P @ Y @ P)
        scale = np.abs(ev).max()
        res = np.linalg.norm(Y @ K) / (np.linalg.norm(K) * scale)
        evc = np.linalg.eigvalsh(P @ Y @ P + scale * Qk @ Qk.T)
        gap = evc[0]
        ub = min(s[h] / (2 * (S - s[h])) for h in range(k))
        worst['neg'] = min(worst['neg'], np.linalg.eigvalsh(Y)[0] / scale)
        worst['kerdev'] = max(worst['kerdev'], res)
        worst['ratio_ub'] = min(worst['ratio_ub'], gap / ub)
        worst['ratio_ub_max'] = max(worst['ratio_ub_max'], gap / ub)
        worst['ratio_x'] = min(worst['ratio_x'], gap / x.min())
        if gap <= 1e-12 * scale or gap > ub * (1 + 1e-9):
            ok = False
            print("  violation:", k, w, gap, ub)
        # hard part
        Yf = 4 * Y
        A = Yf[:N1, :N1]; Bc = Yf[:N1, N1:]
        F = A - 0.5 * Bc @ Bc.T
        Sb = np.hstack(soft_basis(m))
        Q, _ = np.linalg.qr(np.hstack([kernel_basis(m), Sb]))
        Rb = np.linalg.svd(np.eye(N1) - Q @ Q.T)[0][:, :N1 - Q.shape[1]]
        worst['lamR'] = min(worst['lamR'], np.linalg.eigvalsh(Rb.T @ F @ Rb)[0])
    print(f"k={k} ({NS[k]} weight vectors): min eig/scale {worst['neg']:+.1e}; max |Y kappa| rel {worst['kerdev']:.1e}; "
          f"gap/ub in [{worst['ratio_ub']:.4f}, {worst['ratio_ub_max']:.4f}]; min gap/x_min {worst['ratio_x']:.4f}; "
          f"min lambda(F|R) {worst['lamR']:.4f}  ({time.time() - t0:.0f}s)", flush=True)
    ok &= worst['neg'] > -1e-12 and worst['kerdev'] < 1e-10
# direct construction spot checks
for k in (4, 5):
    for trial in range(2):
        w = np.exp(rng.normal(size=k)); w[0] *= 0.05
        Yd = direct_block(list(w))
        Yc = np.array(Mid(list(w), float).full(), float) / 4
        err = np.abs(Yd - Yc).max()
        print(f"direct vs closed form, k={k}: {err:.1e}")
        ok &= err < 1e-10
print("PASS" if ok else "FAIL")
