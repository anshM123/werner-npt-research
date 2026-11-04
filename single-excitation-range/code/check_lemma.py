"""The rank-two diagonal-trace lemma: for Y of rank <= 2 made of m x m square blocks Y_st of size b,
    sum_s |Tr Y_ss|^2 <= ||Y||_F^2 + |Tr Y|^2 / 2.
(a) adversarial minimisation of the slack / ||Y||^2 over Y = U V^* (block sizes b = 1, 2, 3; m = 1..8);
(b) the equality cases named in the README;
(c) the matrix step: I - diag(u^2) + (1/2) u u^T is PSD whenever sum u^2 = 2 (random and adversarial u);
(d) the explicit block matrix showing that the uniform split fails for c >= 3: F_c(Y)/||Y||^2 = 2/c - 1."""
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(11)


def slack(Y, m, b):
    B = Y.reshape(m, b, m, b)
    trs = np.array([np.trace(B[s, :, s, :]) for s in range(m)])
    return np.linalg.norm(Y) ** 2 + abs(trs.sum()) ** 2 / 2 - np.sum(abs(trs) ** 2)


worst = np.inf
for b in (1, 2, 3):
    for m in range(1, 9):
        N = m * b
        def obj(x):
            U = (x[:2 * N] + 1j * x[2 * N:4 * N]).reshape(N, 2)
            V = (x[4 * N:6 * N] + 1j * x[6 * N:]).reshape(N, 2)
            Y = U @ V.conj().T
            return slack(Y, m, b) / np.linalg.norm(Y) ** 2
        best = min(minimize(obj, rng.normal(size=8 * N), method='L-BFGS-B').fun for _ in range(6))
        worst = min(worst, best)
        print(f"(a) block size {b}, m = {m}: adversarial min slack/||Y||^2 = {best:+.2e}", flush=True)
# (b) equality cases
Y1 = np.kron(np.diag([1.0, 0, 0]), np.eye(2))                    # one block equal to the identity
print("(b) Y = |1><1| (x) I_2: slack =", slack(Y1, 3, 2))
Y2 = np.kron(np.array([[1, 1], [-1, -1]]), np.eye(2))           # rank 2, blocks I, I, -I, -I
print("(b) Y = [[I, I], [-I, -I]]: slack =", slack(Y2, 2, 2))
Y3 = np.kron(np.array([[1, 1], [1, 1]]), np.eye(2))              # rank 2, all blocks I: strict
print("(b) Y = [[I, I], [I, I]] (strict case): slack =", slack(Y3, 2, 2))
# (c) the matrix step
mins = []
for trial in range(20000):
    m = rng.integers(1, 9)
    u2 = rng.dirichlet(np.ones(m) * rng.choice([0.2, 1.0, 5.0])) * 2.0
    u = np.sqrt(u2)
    M = np.eye(m) - np.diag(u2) + 0.5 * np.outer(u, u)
    mins.append(np.linalg.eigvalsh(M)[0])
print(f"(c) min eigenvalue of I - diag(u^2) + uu^T/2 over 20000 random u with sum u^2 = 2: {min(mins):+.2e}")
def mobj(x):
    u2 = 2 * np.exp(x) / np.exp(x).sum(); u = np.sqrt(u2)
    return np.linalg.eigvalsh(np.eye(len(x)) - np.diag(u2) + 0.5 * np.outer(u, u))[0]
madv = min(minimize(mobj, rng.normal(size=m), method='Nelder-Mead').fun for m in range(2, 9) for _ in range(5))
print(f"(c) adversarial min eigenvalue: {madv:+.2e}")
# (d) uniform split fails for c >= 3
a = 1.0
Y = np.kron(np.array([[1, 1], [-1, -1]]), np.eye(2)) * a          # Y_ss = aI, Y_tt = -aI, Y_st = aI, Y_ts = -aI; rank 2
B = Y.reshape(2, 2, 2, 2)
off = np.linalg.norm(B[0, :, 1, :]) ** 2 + np.linalg.norm(B[1, :, 0, :]) ** 2
dg = np.linalg.norm(B[0, :, 0, :]) ** 2 + np.linalg.norm(B[1, :, 1, :]) ** 2
trs = np.array([np.trace(B[0, :, 0, :]), np.trace(B[1, :, 1, :])])
for c in (1, 2, 3, 4, 6):
    F = (4 / c) * off + 2 * dg + abs(trs.sum()) ** 2 - 2 * np.sum(abs(trs) ** 2)
    print(f"(d) rank {np.linalg.matrix_rank(Y)}, c = {c}: F_c(Y)/||Y||^2 = {F / np.linalg.norm(Y) ** 2:+.6f} (2/c - 1 = {2 / c - 1:+.6f})")
print(f"RESULT: lemma adversarial min {worst:+.2e}; matrix step min {min(min(mins), madv):+.2e}")
assert worst > -1e-7 and min(mins) > -1e-12 and madv > -1e-9
print("PASS")
