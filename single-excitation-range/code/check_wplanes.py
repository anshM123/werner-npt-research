"""Checks for Theorem 2 (W-planes, arbitrary co-range).
(a) Equal-weight central spin: lambda_max(sum_{s=1}^m P^singlet_{sV}) = (m+1)/2, m = 1..8.
(b) Weighted central spin: lambda_max(sum_s w_s P^singlet_{sV}) <= (sum w + max w)/2 (random + adversarial), m = 1..6.
(c) Hole form of the sector inequality, built from the sector formula, against direct partial traces (n = 3, 4, all h).
(d) For range = W-plane and weighted W-planes inside the weight condition: the exact minimum over ALL co-ranges of
    q_n/||C||^2 (generalized eigenvalue problem), n = 2..5, with random rotations W_s of the excitation basis."""
import numpy as np, itertools
from scipy.optimize import minimize
from common import q_direct, strings, weight, single_excitations, sector_form, embed

rng = np.random.default_rng(17)
d = 3
sing = np.array([0, 1, -1, 0]) / np.sqrt(2); Ps = np.outer(sing, sing)


def embed2(P2, nq, a, b):
    op = np.kron(P2, np.eye(2 ** (nq - 2))).reshape([2] * (2 * nq))
    rest = [i for i in range(nq) if i not in (a, b)]
    perm = np.argsort([a, b] + rest)
    return op.transpose(list(perm) + [p + nq for p in perm]).reshape(2 ** nq, 2 ** nq)


ok = True
for m in range(1, 9):
    lm = np.linalg.eigvalsh(sum(embed2(Ps, m + 1, s, m) for s in range(m)))[-1]
    print(f"(a) m={m}: lambda_max = {lm:.10f}, (m+1)/2 = {(m + 1) / 2}")
    ok &= abs(lm - (m + 1) / 2) < 1e-9
worst_b = -np.inf
for m in range(1, 7):
    ops = [embed2(Ps, m + 1, s, m) for s in range(m)]
    # the inequality is homogeneous in w, so normalise w to sum 1 (softmax) to avoid overflow
    gap = lambda x: (lambda w: np.linalg.eigvalsh(sum(wi * o for wi, o in zip(w, ops)))[-1] - (w.sum() + w.max()) / 2)(
        np.exp(np.clip(x, -50, 50) - np.max(np.clip(x, -50, 50))) / np.sum(np.exp(np.clip(x, -50, 50) - np.max(np.clip(x, -50, 50)))))
    g = max([gap(rng.normal(size=m) * 2) for _ in range(1000)] +
            [-minimize(lambda x: -gap(x), rng.normal(size=m), method='Nelder-Mead', options={'maxiter': 2000}).fun for _ in range(5)])
    worst_b = max(worst_b, g)
    print(f"(b) m={m}: max of lambda_max - (sum w + max w)/2 = {g:+.2e}")
ok &= worst_b < 1e-9


def hole_form_matches(n, h, trials=3):
    """Compare RHS - LHS of the hole form (with the general gauge U^*U = I) against the sector formula."""
    rows = single_excitations(n); c = n - h
    cols = [st for st in strings(n) if weight(st) == c]
    err = 0
    for _ in range(trials):
        U = rng.normal(size=(2 * n, 2)) + 1j * rng.normal(size=(2 * n, 2))
        w, Q = np.linalg.eigh(U.conj().T @ U); U = U @ Q @ np.diag(w ** -0.5) @ Q.conj().T   # U^*U = I
        V = rng.normal(size=(len(cols), 2)) + 1j * rng.normal(size=(len(cols), 2))
        Xs = U @ V.conj().T
        F = {sig: V[j].conj() for j, sig in enumerate(cols)}                                 # X[(s,z);sig] = (U_s F_sig)_z
        A = [U[2 * s:2 * s + 2].conj().T @ U[2 * s:2 * s + 2] for s in range(n)]
        rhs = sum((F[sig].conj() @ (2 * np.eye(2) + 2 * sum((A[u] for u in range(n) if sig[u] == 0), np.zeros((2, 2)))) @ F[sig]).real for sig in cols)
        lhs = 0
        for rho in strings(n):
            if weight(rho) != c - 1:
                continue
            taus = []
            for s in range(n):
                if rho[s] != 0:
                    continue
                tau = 0
                for z in (1, 2):
                    sig = list(rho); sig[s] = z
                    tau += (U[2 * s:2 * s + 2] @ F[tuple(sig)])[z - 1]
                taus.append(tau)
            taus = np.array(taus)
            lhs += 2 * np.sum(abs(taus) ** 2) - abs(taus.sum()) ** 2
        direct = 2 ** (n - c + 1) * q_direct(embed(Xs, n, rows, cols), n)
        err = max(err, abs((rhs - lhs) - direct) / max(1, abs(direct)))
    return err


worst_c = 0
for n in (3, 4):
    for h in range(0, n):
        e = hole_form_matches(n, h)
        worst_c = max(worst_c, e)
        print(f"(c) n={n} h={h}: hole form vs direct partial traces, relative error {e:.1e}")
ok &= worst_c < 1e-9


def min_ratio(n, lam, Ws):
    D = d ** n; st_list = strings(n); idx = {s: i for i, s in enumerate(st_list)}
    U = np.zeros((D, 2), complex)
    for s in range(n):
        for z in range(2):
            st = [0] * n; st[s] = z + 1
            for k in range(2):
                U[idx[tuple(st)], k] += lam[s] * Ws[s][z, k]
    basis = []
    for j in range(D):
        for k in range(2):
            V = np.zeros((D, 2), complex); V[j, k] = 1; basis.append(U @ V.conj().T)
    T = [None] * len(basis)
    from common import q_direct as _q
    # Gram of the polarised form: q(a,b) via polarisation of the quadratic form
    def Bmap(X):
        Tt = X.reshape([d] * (2 * n))
        for s in range(n):
            tr = np.trace(Tt, axis1=s, axis2=n + s)
            full = np.moveaxis(np.einsum('ab,...->...ab', np.eye(d), tr), [2 * n - 2, 2 * n - 1], [s, n + s])
            Tt = Tt - 0.5 * full
        return Tt.reshape(D, D)
    imgs = [Bmap(X) for X in basis]
    Qm = np.array([[np.vdot(a, b) for b in imgs] for a in basis]); Nm = np.array([[np.vdot(a, b) for b in basis] for a in basis])
    Qm = (Qm + Qm.conj().T) / 2; Nm = (Nm + Nm.conj().T) / 2
    w, P = np.linalg.eigh(Nm); keep = w > 1e-10; Pk = P[:, keep] / np.sqrt(w[keep])
    return np.linalg.eigvalsh(Pk.conj().T @ Qm @ Pk)[0]


worst_d = np.inf
for n in range(2, 6):
    for kind in ("W-plane", "weighted"):
        vals = []
        for trial in range(3 if n < 5 else 1):
            lam2 = np.ones(n) if kind == "W-plane" else rng.uniform(1, 2, size=n)
            lam = np.sqrt(lam2 / lam2.sum())
            Ws = [np.eye(2)] * n if trial == 0 else [np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))[0] for _ in range(n)]
            vals.append(min_ratio(n, lam, Ws))
        worst_d = min(worst_d, min(vals))
        print(f"(d) n={n} {kind}: min over all co-ranges of q/||C||^2 = {min(vals):+.6e}", flush=True)
ok &= worst_d > -1e-10
print("RESULT: central spin exact; weighted bound gap", f"{worst_b:+.1e}", "; hole form error", f"{worst_c:.1e}",
      "; W-plane minimum", f"{worst_d:+.3e}")
assert ok
print("PASS")
