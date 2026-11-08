"""Section 9, conjecture TN_n (floating point).

q(t, C) = sum_J (-t)^|J| ||Tr_J C||_F^2 on C^{d_1} (x) ... (x) C^{d_n}.  TN_n says that for rank-two C
    q(t, C) >= (1 - t)^(n-1) (||C||_F^2 - t ||C||_1^2),
with ||C||_1 = s_1 + s_2 the trace norm.  From random complex factors A, B (D x 2) this minimises
    gap = [q(t, AB^*) - (1 - t)^(n-1)(||C||_F^2 - t ||C||_1^2)] / ||C||_F^2
with L-BFGS.  s_1 s_2 = (det(A^*A) det(B^*B))^{1/2}, so the objective is smooth on rank two.
usage: python tn_bound.py t runs seed d_1 d_2 ... d_n      (qutrits: d_s = 3)"""
import sys, itertools
import numpy as np
import torch

torch.set_default_dtype(torch.float64)
t, runs, seed = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
dims = [int(x) for x in sys.argv[4:]]
n, D = len(dims), int(np.prod(dims))
subsets = [J for r in range(n + 1) for J in itertools.combinations(range(n), r)]


def q(C):
    T = C.reshape(dims + dims)
    tot = 0.
    for J in subsets:
        X, rem = T, list(range(n))
        for s in sorted(J, reverse=True):
            i = rem.index(s)
            X = torch.diagonal(X, dim1=i, dim2=len(rem) + i).sum(-1)
            rem.pop(i)
        tot = tot + (-t) ** len(J) * (X.real ** 2 + X.imag ** 2).sum()
    return tot


def gap(params):
    A = torch.complex(params[0], params[1])
    B = torch.complex(params[2], params[3])
    C = A @ B.conj().T
    n2 = (C.real ** 2 + C.imag ** 2).sum()
    GA, GB = A.conj().T @ A, B.conj().T @ B
    dA = (GA[0, 0] * GA[1, 1] - GA[0, 1] * GA[1, 0]).real
    dB = (GB[0, 0] * GB[1, 1] - GB[0, 1] * GB[1, 0]).real
    s1s2 = torch.sqrt(torch.clamp(dA * dB, min=1e-24))   # floor 1e-12 on s_1 s_2: keeps gradients finite near rank one
    rhs = (1 - t) ** (n - 1) * (n2 - t * (n2 + 2 * s1s2))
    return (q(C) - rhs) / n2, C


rng = np.random.default_rng(seed)
print(f"dims {dims}, t = {t}, {runs} runs, seed {seed}")
best = []
for r in range(runs):
    params = [torch.tensor(rng.standard_normal((D, 2)), requires_grad=True) for _ in range(4)]
    opt = torch.optim.LBFGS(params, lr=1, max_iter=3000, tolerance_grad=1e-13, tolerance_change=1e-16,
                            history_size=50, line_search_fn='strong_wolfe')

    def closure():
        opt.zero_grad()
        v, _ = gap(params)
        v.backward()
        return v
    for _ in range(3):
        opt.step(closure)
    v, C = gap(params)
    Cn = C.detach().numpy()
    if not np.isfinite(v.item()) or not np.all(np.isfinite(Cn)):
        print(f"run {r}: optimiser left the finite range, skipped", flush=True)
        continue
    sv = np.linalg.svd(Cn, compute_uv=False)[:3]
    qv = (q(C) / (C.real ** 2 + C.imag ** 2).sum()).item()
    best.append(v.item())
    print(f"run {r}: gap {v.item():.3e}   s_2/s_1 = {sv[1] / sv[0]:.4f}   q/||C||^2 = {qv:.3e}", flush=True)
print(f"smallest gap {min(best):.3e}")
