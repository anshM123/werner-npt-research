"""Universal anchors (numerical, n = 2, 3). For a unit vector p on n copies let Y_p = (M - R_p)^Gamma, where
M = (I - |Omega><Omega|/2)^{(x)n} on rows (x) cols, R_p = M V D^{-1} V^* M, V = p (x) I, D = V^* M V, and Gamma is the
transpose on the rows. If Y_p >= 0, then q_n >= 0 on every operator whose range plane contains p (Section 1 is the
case of product p). Y_p always vanishes on pbar (x) cols. Prints the smallest eigenvalue and the gap above the kernel."""
import numpy as np, sys
d = 3
rng = np.random.default_rng(3)


def M_matrix(n):
    om = np.zeros(d * d); om[[a * d + a for a in range(d)]] = 1
    M1 = (np.eye(d * d) - 0.5 * np.outer(om, om)).reshape(d, d, d, d)
    M = M1
    for _ in range(1, n):
        M = np.multiply.outer(M, M1)
    out = [4 * s for s in range(n)] + [4 * s + 1 for s in range(n)]
    inn = [4 * s + 2 for s in range(n)] + [4 * s + 3 for s in range(n)]
    return M.transpose(out + inn).reshape(d ** (2 * n), d ** (2 * n))


def anchor(p, n, M):
    D = d ** n
    p = p / np.linalg.norm(p)
    V = np.kron(p.reshape(D, 1), np.eye(D))
    R = M @ V @ np.linalg.solve(V.conj().T @ M @ V, V.conj().T @ M)
    Y = (M - R).reshape(D, D, D, D).transpose(2, 1, 0, 3).reshape(D * D, D * D)
    ev = np.linalg.eigvalsh((Y + Y.conj().T) / 2)
    k = int((abs(ev) < 1e-9).sum())
    return ev[0], k, ev[k] if ev[0] > -1e-9 else float('nan')


def kron(*v):
    out = np.ones(1)
    for x in v:
        out = np.kron(out, x)
    return out


e = np.eye(3)
for n in (2, 3):
    M = M_matrix(n)
    tests = {'product |0...0>': kron(*[e[0]] * n),
             'GHZ |0..0> + |1..1>': kron(*[e[0]] * n) + kron(*[e[1]] * n),
             'W = sum_s |0..1_s..0>': sum(kron(*[e[1] if r == s else e[0] for r in range(n)]) for s in range(n))}
    for k in range(3):
        tests[f'random {k}'] = rng.normal(size=d ** n) + 1j * rng.normal(size=d ** n)
    for name, p in tests.items():
        lo, ker, gap = anchor(p, n, M)
        print(f"n={n} {name:24s}: smallest eigenvalue {lo:+.4e}, kernel {ker}, gap {gap:.4f}", flush=True)
