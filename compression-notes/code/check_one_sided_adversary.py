"""Section 11, one-sided separability: adversary over planes V.  Maximises eps*(sigma_V) (the inner test of
check_one_sided.py, net of P qubit states) by ascent on V, with the envelope gradient d eps* = <-Y, d sigma> (Y the dual
of the equality constraint; a finite-difference check is printed first).  eps* > 0 beyond the net error would signal a
PPT-entangled compression.
usage: python check_one_sided_adversary.py n steps seed [P]   (e.g. 2 120 0 60; about half an hour)
"""
import sys, numpy as np, cvxpy as cp, torch
torch.set_default_dtype(torch.float64)
n, STEPS, SEED = (int(a) for a in sys.argv[1:4]); P = int(sys.argv[4]) if len(sys.argv) > 4 else 120
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2
def power(m):
    R = np.eye(1)
    for _ in range(m): R = np.kron(R, rho)
    perm = [2 * i for i in range(m)] + [2 * i + 1 for i in range(m)]
    return R.reshape([3] * (4 * m)).transpose(perm + [2 * m + q for q in perm]).reshape(9 ** m, 9 ** m)
N = 3 ** n; d = 2 * N
Rt = torch.tensor(power(n), dtype=torch.complex128)
def sigma_t(M):
    Q, _ = torch.linalg.qr(M)
    K = torch.kron(Q.conj().T, torch.eye(N, dtype=torch.complex128))
    S = K @ Rt @ K.conj().T; S = (S + S.conj().T) / 2
    return S / torch.trace(S).real
def net(P):
    i = np.arange(P) + 0.5; th = np.arccos(1 - 2 * i / P); ph = np.pi * (1 + 5 ** 0.5) * i
    return [np.array([np.cos(t / 2), np.exp(1j * f) * np.sin(t / 2)]) for t, f in zip(th, ph)]
A = net(P)
Sp = cp.Parameter((d, d), hermitian=True)
Bs = [cp.Variable((N, N), hermitian=True) for _ in range(P)]
eps = cp.Variable()
tot = sum(cp.kron(np.outer(a, a.conj()), B) for a, B in zip(A, Bs))
eq = tot == Sp + eps * np.eye(d) / d
prob = cp.Problem(cp.Minimize(eps), [B >> 0 for B in Bs] + [eq])
def solve(S):
    Sp.value = (S + S.conj().T) / 2
    prob.solve(solver=cp.CLARABEL)
    return eps.value, eq.dual_value
def ptmin(S):
    T = S.reshape(2, N, 2, N).transpose(0, 3, 2, 1).reshape(d, d)
    return np.linalg.eigvalsh(T)[0]
rng = np.random.default_rng(SEED)
M = torch.tensor(rng.standard_normal((N, 2)) + 1j * rng.standard_normal((N, 2)), requires_grad=True)
# finite-difference check of the envelope gradient (sign convention)
with torch.no_grad(): S0 = sigma_t(M).numpy()
e0, Y = solve(S0)
D = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d)); D = (D + D.conj().T) / 2; D -= np.trace(D) / d * np.eye(d)
h = 1e-4
e1, _ = solve(S0 + h * D)
pred = np.real(np.vdot(Y, D)) if Y is not None else np.nan
sgn = 1.0 if abs((e1 - e0) / h - pred) < abs((e1 - e0) / h + pred) else -1.0
print(f"FD check: (e1-e0)/h = {(e1 - e0) / h:+.4e}, <Y,D> = {pred:+.4e} -> sign {sgn:+.0f}", flush=True)
opt = torch.optim.Adam([M], lr=0.02)
best = (-np.inf, None)
for it in range(STEPS):
    S = sigma_t(M)
    e, Y = solve(S.detach().numpy())
    if e > best[0]: best = (e, M.detach().clone())
    if it % 5 == 0 or it == STEPS - 1:
        print(f"n = {n} seed {SEED} step {it}: eps* {e:+.3e} (PT-min {ptmin(S.detach().numpy()):+.2e}); best {best[0]:+.3e}", flush=True)
    Yt = torch.tensor(sgn * Y, dtype=torch.complex128)
    obj = torch.real(torch.sum(Yt.conj() * S))
    opt.zero_grad(); (-obj).backward(); opt.step()
np.save(f"adv_n{n}_s{SEED}.npy", best[1].numpy())
