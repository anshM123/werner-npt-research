"""Float search over splits q(k,s,t) + q(k,t,s) = c(k,s,t) at a = 3/4 for given n, on the M=0 slice.
Z_n(s,t) = f-matrix of the sum over words with s letters (I-F), t letters (I-Om), rest I  (= 4^n c(k,s,t) Y(k,s,t)).
Level n-2 classes stored; level n-1 classes generated one at a time and projected to level n via shifted bases.
"""
import sys, time
import numpy as np, cvxpy as cp
from bos import mindices, index_of
from exact_t2 import one_copy_int, trf
import exact_greedy as eg

n = int(sys.argv[1]); t0 = time.time()
I, F, Om = [x.astype(float) for x in one_copy_int()]
XA = I - F; XB = I - Om
# classes up to level n-2
Z = {(0, 0): np.array([[1.0]])}
for m in range(1, n - 1):
    newZ = {}
    for s in range(m + 1):
        for t in range(m + 1 - s):
            acc = 0
            if s + t <= m - 1: acc = acc + trf(I, Z[(s, t)], m)
            if s >= 1: acc = acc + trf(XA, Z[(s - 1, t)], m)
            if t >= 1: acc = acc + trf(XB, Z[(s, t - 1)], m)
            newZ[(s, t)] = acc
    Z = newZ
print(f"level {n-2} classes built ({len(Z)}), {time.time()-t0:.1f}s", flush=True)
cols, labels = eg.basis_cols(n); dB = len(cols)
Np = len(mindices(n - 1)); idxn = index_of(n); alp = mindices(n - 1)
# float shifted bases with column normalisation
Bfull = np.zeros((len(mindices(n)) ** 2, dB)); Nn = len(mindices(n))
for r, v in enumerate(cols):
    for (ia, ib), val in v.items():
        Bfull[ia * Nn + ib, r] = val
scale = np.linalg.norm(Bfull, axis=0); Bfull /= scale
Bs = {}
for i in range(3):
    for j in range(3):
        up_i = np.array([idxn[tuple(a[t] + (1 if t == i else 0) for t in range(3))] for a in alp])
        up_j = np.array([idxn[tuple(a[t] + (1 if t == j else 0) for t in range(3))] for a in alp])
        Bs[i, j] = Bfull[(up_i[:, None] * Nn + up_j[None, :]).reshape(-1), :]
Bd = Bs[0, 0] + Bs[1, 1] + Bs[2, 2]
def proj(Y, kind):
    M = np.zeros((dB, dB))
    for i in range(3):
        for j in range(3):
            W = Bs[i, j] if kind == 'I' else (Bs[i, j] - Bs[j, i] if kind == 'A' else Bs[i, j] - (Bd if i == j else 0))
            M += Bs[i, j].T @ (Y @ W)
    return M
Mz = {}
m = n - 1
for s in range(m + 1):
    for t in range(m + 1 - s):
        acc = 0
        if s + t <= m - 1: acc = acc + trf(I, Z[(s, t)], m)
        if s >= 1: acc = acc + trf(XA, Z[(s - 1, t)], m)
        if t >= 1: acc = acc + trf(XB, Z[(s, t - 1)], m)
        Y = acc
        for (ds, dt, kind) in ((0, 0, 'I'), (1, 0, 'A'), (0, 1, 'B')):
            key = (s + ds, t + dt)
            if key[0] + key[1] <= n:
                Mz[key] = Mz.get(key, 0) + proj(Y, kind)
        del acc, Y
print(f"level {n} classes projected ({len(Mz)}), {time.time()-t0:.1f}s", flush=True)
del Z
Mz = {k: (v + v.T) / 2 for k, v in Mz.items()}
np.savez(f'mz{n}.npz', **{f'{k[0]}_{k[1]}': v for k, v in Mz.items()})
if len(sys.argv) > 2 and sys.argv[2] == 'saveonly':
    print('saved', flush=True); sys.exit(0)
# greedy and forced kernel
fixed = sum(Mz[(s, 0)] for s in range(1, n + 1)) + 0.5 * sum(Mz[(s, s)] for s in range(1, n // 2 + 1))
pairs = [(s, t) for t in range(1, n) for s in range(t + 1, n + 1) if s + t <= n]
greedy = fixed + sum(Mz[p] for p in pairs)
print("greedy min eig (column-normalised slice):", np.linalg.eigvalsh(greedy)[:3], flush=True)
stack = np.vstack([greedy] + [Mz[p] - Mz[(p[1], p[0])] for p in pairs])
u, sv, vt = np.linalg.svd(stack)
ker = vt[sv < 1e-9 * sv[0]] if np.any(sv < 1e-9 * sv[0]) else np.zeros((0, dB))
nk = dB - np.sum(sv > 1e-9 * sv[0])
Vperp = vt[:dB - nk].T
print(f"forced common kernel dim {nk} (expected n+1 = {n+1})", flush=True)
phi = cp.Variable(len(pairs)); tt = cp.Variable()
expr = fixed + sum(phi[i] * Mz[p] + (1 - phi[i]) * Mz[(p[1], p[0])] for i, p in enumerate(pairs))
Mr = Vperp.T @ expr @ Vperp
prob = cp.Problem(cp.Maximize(tt), [(Mr + Mr.T) / 2 >> tt * np.eye(Vperp.shape[1]), cp.abs(phi) <= 20])
try:
    prob.solve(solver='CLARABEL')
except Exception as e:
    print('clarabel failed, SCS', e, flush=True)
    prob.solve(solver='SCS', eps=1e-7, max_iters=50000)
print(f"n={n}: best min eigenvalue over splits (on complement of forced kernel) = {tt.value:.4e} ({prob.status}); "
      f"greedy value = {np.linalg.eigvalsh(Vperp.T @ greedy @ Vperp)[0]:.4e}", flush=True)
print({p: round(float(phi.value[i]), 4) for i, p in enumerate(pairs)})
