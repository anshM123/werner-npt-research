"""Fast check of the mixed-level candidate J at n = 4 (numpy, real+complex random and structured pairs, then
L-BFGS on the violation).  Prints the worst violation found and the components."""
import sys, itertools, numpy as np
from scipy.optimize import minimize
n, d = int(sys.argv[1]), 3
NS = int(sys.argv[2]); seed = int(sys.argv[3])
rng = np.random.default_rng(seed)
subsets = [T for k in range(n + 1) for T in itertools.combinations(range(n), k)]
D = d ** n
def red(x, T):
    X = x.reshape((d,) * n); rest = [i for i in range(n) if i not in T]
    X = np.transpose(X, list(T) + rest).reshape(d ** len(T), d ** len(rest)); return X @ X.conj().T
def tvec(x, y):
    return np.array([np.real(np.vdot(x, x) * np.vdot(y, y)) if len(T) == 0 else np.real(np.trace(red(x, T) @ red(y, T))) for T in subsets])
Wk = np.zeros((n + 1, len(subsets)))
for K in subsets:
    for j, T in enumerate(subsets): Wk[len(K), j] += (-1) ** len(set(K) & set(T)) / 2 ** n
def comps(z):
    a = z[:D] + 1j * z[D:2 * D]; b = z[2 * D:3 * D] + 1j * z[3 * D:]
    u = a / np.linalg.norm(a); b = b - np.vdot(u, b) * u; v = b / np.linalg.norm(b)
    pk = Wk @ tvec(u, v); Au = Wk @ tvec(u, u); Av = Wk @ tvec(v, v)
    lhs = sum((9 ** j - 1) * (3 * (pk[2 * j + 1] if 2 * j + 1 <= n else 0) - pk[2 * j]) for j in range(1, n // 2 + 1))
    rhs = sum((9 ** j - 1) * np.sqrt(max(Au[2 * j] * Av[2 * j], 0)) for j in range(1, n // 2 + 1))
    gfull = np.sqrt((1 + sum((9 ** j - 1) * Au[2 * j] for j in range(1, n // 2 + 1))) * (1 + sum((9 ** j - 1) * Av[2 * j] for j in range(1, n // 2 + 1))))
    return lhs, rhs, gfull, pk, Au, Av
def f(z):
    lhs, rhs, *_ = comps(z); return -(lhs - rhs)
best = (-1e9, None)
for s in range(NS):
    z0 = rng.standard_normal(4 * D) * (1.0 if s % 3 else 0.2)
    r = minimize(f, z0, method='L-BFGS-B', options={'maxiter': 3000})
    lhs, rhs, gfull, pk, Au, Av = comps(r.x)
    viol = lhs - rhs; psd = (1 + lhs) - gfull   # PSD endpoint violation (-g - sqrt(g g'))
    print(f"n={n} start {s}: J-viol {viol:+.4e} (lhs {lhs:.4f}, rhs {rhs:.4f});  endpoint margin {-psd:+.4e};  A_u={np.round(Au[2::2],4)} A_v={np.round(Av[2::2],4)}", flush=True)
    if viol > best[0]: best = (viol, r.x)
print(f"n={n}: MAX J violation {best[0]:+.4e}", flush=True)
np.save(f"jquick_n{n}_s{seed}.npy", best[1])
