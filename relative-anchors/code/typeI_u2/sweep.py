"""Numerical sweep (floats; no claim rests on it beyond 'evidence') of the three remaining inequalities:
 (C_u)  lambda_min of the U-block (X_U modulo the kernel, q) after eliminating X_T  > 0,
 (A_u)  sup ||G^_tt'||_2 S_T / sqrt(s_t s_t'),
 (B_u)  inf lambda_min(G^_tt) / (1 - s_t/S_T)   (m >= 2),
 (W)    the weighted block dominance ratio max_t sum_t' ||G^_tt'|| sqrt(s_t'/s_t) / lambda_min(G^_tt),
and the full block: smallest eigenvalue off the kernel and |Y kernel|.  G^ = Schur complement on the T-blocks with
the w_U-component of zeta_t rescaled by c = sqrt(1 + s_U/S_T).  Weights: log-uniform with random spreads, and
families with s_T or s_U small."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI
from ex5 import null_basis


def analyse(w, T, rescale=True):
    ti = TypeI(w, T, F=float)
    M = np.array(ti.build())
    U = ti.U; u = len(U); m = len(T); w = np.asarray(w, float); s = w ** 2
    sT, sU = s[T].sum(), s[U].sum()
    n = ti.n
    wU = w[U]; nw = np.linalg.norm(wU); Bw = null_basis(wU[None, :])
    # kernel check
    Kv = np.zeros((n, u))
    for j, h in enumerate(U):
        for r in range(len(w)):
            Kv[ti.idx[('P', r, h)], j] = w[r]
    kres = np.abs(M @ Kv).max() / np.abs(M).max()
    Qk, _ = np.linalg.qr(Kv)
    Pc = np.eye(n) - Qk @ Qk.T
    evf = np.linalg.eigvalsh(Pc @ M @ Pc + 10 * np.abs(M).max() * Qk @ Qk.T)
    lam_full = evf[0]
    # U-block after eliminating X_T
    def XU(X):
        v = np.zeros(n)
        for a, r in enumerate(U):
            for b, h in enumerate(U):
                v[ti.idx[('P', r, h)]] = X[a, b]
        return v
    ub = [XU(np.outer(Bw[:, i], Bw[:, j])) for i in range(u - 1) for j in range(u - 1)]
    ub += [XU(np.outer(Bw[:, i], wU / nw)) for i in range(u - 1)]
    v = np.zeros(n); v[ti.idx[('Q',)]] = 1; ub.append(v)
    xt = []
    for t in T:
        for h in U:
            v = np.zeros(n); v[ti.idx[('P', t, h)]] = 1; xt.append(v)
    UB = np.array(ub).T; XT = np.array(xt).T
    Muu = UB.T @ M @ UB; Mux = UB.T @ M @ XT; Mxx = XT.T @ M @ XT
    SU = Muu - Mux @ np.linalg.solve(Mxx, Mux.T)
    lamC = np.linalg.eigvalsh((SU + SU.T) / 2)[0]
    # T-blocks
    c = np.sqrt(1 + sU / sT) if rescale else 1.0
    tb = []
    for t in T:
        v = np.zeros(n); v[ti.idx[('R', t)]] = 1; tb.append(v)
        v = np.zeros(n)
        for j, h in enumerate(U):
            v[ti.idx[('Z', t, h)]] = c * wU[j] / nw
        tb.append(v)
        for i in range(u - 1):
            v = np.zeros(n)
            for j, h in enumerate(U):
                v[ti.idx[('Z', t, h)]] = Bw[j, i]
            tb.append(v)
    TB = np.array(tb).T
    EL = np.hstack([UB, XT])
    A11 = EL.T @ M @ EL; A12 = EL.T @ M @ TB; A22 = TB.T @ M @ TB
    G = A22 - A12.T @ np.linalg.solve(A11, A12); G = (G + G.T) / 2
    d = u + 1
    lam = [np.linalg.eigvalsh(G[i * d:(i + 1) * d, i * d:(i + 1) * d])[0] for i in range(m)]
    A = 0.0; B = np.inf; R = 0.0
    for i, t in enumerate(T):
        if m > 1:
            B = min(B, lam[i] / (1 - s[t] / sT))
        acc = 0.0
        for j, tp in enumerate(T):
            if i == j:
                continue
            nr = np.linalg.norm(G[i * d:(i + 1) * d, j * d:(j + 1) * d], 2)
            A = max(A, nr * sT / np.sqrt(s[t] * s[tp]))
            acc += nr * np.sqrt(s[tp] / s[t])
        R = max(R, acc / lam[i])
    return dict(kres=kres, lam_full=lam_full, lamC=lamC, A=A, B=B, R=R, lamG=min(lam), relC=(lamC / (sT / (sT + sU)) if sT / (sT + sU) > 1e-8 else np.inf))


if __name__ == "__main__":
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    ntr = int(sys.argv[2]) if len(sys.argv) > 2 else 120
    t0 = time.time()
    ok = True
    print("u m | max|Y ker| | min eig full (off ker) | min lam(C_u) | min lam(C_u)/(s_T/|s|) [s_T/|s| > 1e-8] | sup A | inf B | max W | min lam(G^_tt) | control")
    for u in (2, 3, 4, 5):
        for m in (1, 2, 3, 4):
            if u + m > 8:
                continue
            k = u + m; T = list(range(m))
            agg = dict(kres=0, lam_full=np.inf, lamC=np.inf, relC=np.inf, A=0, B=np.inf, R=0, lamG=np.inf)
            for trial in range(ntr):
                spread = rng.choice([0.5, 1.5, 3, 6])
                w = np.exp(rng.uniform(-spread, spread, k))
                if trial % 6 == 1: w[:m] *= np.exp(rng.uniform(-6, -1))
                if trial % 6 == 2: w[m:] *= np.exp(rng.uniform(-6, -1))
                if trial % 6 == 3: w[rng.integers(0, k)] *= np.exp(rng.uniform(2, 6))
                r = analyse(w, T)
                for key in ('kres', 'A', 'R'):
                    agg[key] = max(agg[key], r[key])
                for key in ('lam_full', 'lamC', 'relC', 'B', 'lamG'):
                    agg[key] = min(agg[key], r[key])
            ctrl = 0.0
            if m >= 2:
                for trial in range(10):
                    w = np.exp(rng.uniform(-1, 1, k)); w[:m] *= 1e-3
                    ctrl = max(ctrl, analyse(w, T, rescale=False)['R'])
            agg['ctrl'] = ctrl
            good = (m == 1 or ctrl > 1) and agg['kres'] < 1e-10 and agg['lam_full'] > -1e-12 and agg['relC'] > 0.49 and agg['R'] < 1 and agg['lamG'] > 0
            ok &= good
            print(f"{u} {m} | {agg['kres']:.1e} | {agg['lam_full']:.2e} | {agg['lamC']:.2e} | {agg['relC']:.3f} | {agg['A']:.4f} | {agg['B']:.4f} | {agg['R']:.4f} | {agg['lamG']:.2e} | control (no rescaling) max W = {agg['ctrl']:.1f}", flush=True)
    print(f"time {time.time() - t0:.0f}s")
    print("PASS" if ok else "FAIL")
