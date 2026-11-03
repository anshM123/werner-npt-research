"""Section 1: the compression conjecture QC.

   For a copy s and a unit vector w, C_{w perp} = P C P with P = I - |w><w| on copy s. The script
   - checks the exact formula
         Delta_s(w) = q_n(C) - q_n(C_{w perp})
                    = row(w) + col(w) - (1 - t) q_{n-1}(X_w) - 2t Re q_{n-1}(Tr_s C, X_w)
     and the trace identity Tr M_s = 2 q_n(C) for the quadratic part M_s, on random operators;
   - evaluates two canonical choices of w on the stored operators where they fail;
   - evaluates max_{s,w} Delta_s(w) (many local starts) on all stored operators."""
import glob
import numpy as np
from scipy.optimize import minimize
from werner import d, q, applyB, contract

rng = np.random.default_rng(0)
I3 = np.eye(d)


def blocks(C, n, s):
    return [[contract(C, n, s, I3[i], I3[j]) for j in range(d)] for i in range(d)]


def delta(C, n, s, w, t, X=None):
    w = w / np.linalg.norm(w)
    X = blocks(C, n, s) if X is None else X
    T = sum(X[a][a] for a in range(d))
    row = sum(q(sum(np.conj(w[i]) * X[i][l] for i in range(d)), n - 1, t) for l in range(d))
    col = sum(q(sum(w[j] * X[l][j] for j in range(d)), n - 1, t) for l in range(d))
    Xw = sum(np.conj(w[i]) * w[j] * X[i][j] for i in range(d) for j in range(d))
    return row + col - (1 - t) * q(Xw, n - 1, t) - 2 * t * np.vdot(applyB(T, n - 1, t), Xw).real


def quad_part(C, n, s, t):
    """Hermitian M_s with w^* M_s w = row(w) + col(w) - 2t Re q(Tr_s C, X_w)"""
    X = blocks(C, n, s)
    T = sum(X[a][a] for a in range(d))

    def Qd(w):
        row = sum(q(sum(np.conj(w[i]) * X[i][l] for i in range(d)), n - 1, t) for l in range(d))
        col = sum(q(sum(w[j] * X[l][j] for j in range(d)), n - 1, t) for l in range(d))
        Xw = sum(np.conj(w[i]) * w[j] * X[i][j] for i in range(d) for j in range(d))
        return row + col - 2 * t * np.vdot(applyB(T, n - 1, t), Xw).real

    M = np.zeros((d, d), complex)
    for i in range(d):
        M[i, i] = Qd(I3[i])
    for i in range(d):
        for j in range(i + 1, d):
            a = Qd(I3[i] + I3[j]) - M[i, i] - M[j, j]
            b = Qd(I3[i] + 1j * I3[j]) - M[i, i] - M[j, j]
            # Qd(e_i + i e_j) - Qd(e_i) - Qd(e_j) = -2 Im M_ij for Qd(w) = w^* M w
            M[i, j] = (a - 1j * b) / 2
            M[j, i] = np.conj(M[i, j])
    return M


def max_delta(C, n, t, starts=12):
    best = -np.inf
    for s in range(n):
        X = blocks(C, n, s)
        for _ in range(starts):
            r = minimize(lambda x: -delta(C, n, s, x[:3] + 1j * x[3:], t, X), rng.standard_normal(6),
                         method='L-BFGS-B', options={'maxiter': 2000})
            best = max(best, -r.fun)
    return best / np.vdot(C, C).real


def marginal(M, n, s):
    T = M.reshape([d] * (2 * n))
    T = np.moveaxis(T, [s, n + s], [0, n]).reshape(d, d ** (n - 1), d, d ** (n - 1))
    return np.einsum('aibi->ab', T)


def rule_least_weight(C, n, s, t):
    ev, V = np.linalg.eigh(marginal(C @ C.conj().T + C.conj().T @ C, n, s))
    return V[:, 0]


def rule_top_eigenvector(C, n, s, t):
    ev, V = np.linalg.eigh(quad_part(C, n, s, t))
    return V[:, -1]


print("exact formula and trace identity, random rank-2 C:")
for n in (2, 3, 4):
    D = d ** n
    for t in (0.3, 0.5):
        C = (rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2))) @ \
            (rng.standard_normal((2, D)) + 1j * rng.standard_normal((2, D)))
        s = int(rng.integers(n))
        w = rng.standard_normal(d) + 1j * rng.standard_normal(d); w /= np.linalg.norm(w)
        P = np.eye(1)
        for r_ in range(n):
            P = np.kron(P, I3 - np.outer(w, w.conj()) if r_ == s else I3)
        lhs = q(C, n, t) - q(P @ C @ P, n, t)
        rhs = delta(C, n, s, w, t)
        M = quad_part(C, n, s, t)
        tr = np.trace(M).real
        # w^* M_s w must equal Delta_s(w) + (1 - t) q_{n-1}(X_w) at complex w
        residual = abs(np.vdot(w, M @ w).real - rhs - (1 - t) * q(contract(C, n, s, w, w), n - 1, t))
        assert residual < 1e-11 * np.vdot(C, C).real
        print(f"  n={n} t={t} site {s}: q_n(C) - q_n(C_wperp) = {lhs:.10f}, formula = {rhs:.10f};  "
              f"Tr M_s = {tr:.8f}, 2 q_n(C) = {2 * q(C, n, t):.8f};  quadratic residual = {residual:.1e}", flush=True)

t = 0.5
print("\ncanonical rules on the stored two-copy operators (value = max over copies of Delta_s(w_s)/||C||^2):")
for name, rule in [('rule_least_weight_fails_n2', rule_least_weight),
                   ('rule_top_eigenvector_fails_n2', rule_top_eigenvector)]:
    C = np.load(f"../data/{name}.npy"); n = 2
    val = max(delta(C, n, s, rule(C, n, s, t), t) for s in range(n)) / np.vdot(C, C).real
    print(f"  {name}: rule value {val:+.4e};  max over all w: {max_delta(C, n, t):+.4e};  "
          f"q_2/||C||^2 = {q(C, n, t) / np.vdot(C, C).real:+.4e}", flush=True)

print("\nmax Delta on the stored adversarial operators:")
for name in sorted(glob.glob("../data/qc_adversarial_n*.npy")):
    C = np.load(name); n = int(round(np.log(C.shape[0]) / np.log(d)))
    print(f"  {name.split('/')[-1]}: max Delta/||C||^2 = {max_delta(C, n, t, starts=8):+.4e};  "
          f"q_n/||C||^2 = {q(C, n, t) / np.vdot(C, C).real:+.4e}", flush=True)
