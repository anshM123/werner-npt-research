"""Census of the zeros of the endpoint form on rank-two operators (Section 8).

q_n(C) = sum_J (-1/2)^{|J|} ||Tr_J C||^2 on (C^3)^{(x) n}.  From random complex factors A, B (3^n x 2) minimise
q_n(AB^*)/||AB^*||^2 with L-BFGS and report, for the final point: the value, the singular values, the rank of the
one-copy marginal of range + co-range at each copy ("joint"), and the smallest ratio s_2/s_1 of operator Schmidt
coefficients over all cuts I | I^c (a ratio near 0 means C = R_I (x) K_{I^c}).
A zero is of the predicted form R_I (x) K_J when the joint rank is 2 at every copy (C lives on a product of local
qubits; Section 8) or when C factorises across a cut with the joint rank 2 on the J side.
usage: python zero_census.py n runs seed"""
import sys, itertools
import numpy as np
import torch

torch.set_default_dtype(torch.float64)
n, runs, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
d = 3 ** n
subsets = [K for r in range(n + 1) for K in itertools.combinations(range(n), r)]


def qn(C):
    T = C.reshape([3] * n + [3] * n)
    tot = 0.
    for K in subsets:
        X, rem = T, list(range(n))
        for s in sorted(K, reverse=True):
            i = rem.index(s)
            X = torch.diagonal(X, dim1=i, dim2=len(rem) + i).sum(-1)
            rem.pop(i)
        tot = tot + (-0.5) ** len(K) * (X.real ** 2 + X.imag ** 2).sum()
    return tot


def value(params):
    A = torch.complex(params[0], params[1])
    B = torch.complex(params[2], params[3])
    C = A @ B.conj().T
    return qn(C) / (C.real ** 2 + C.imag ** 2).sum(), C


def joint_ranks(C, tol=1e-6):
    U, S, Vh = np.linalg.svd(C)
    out = []
    for s in range(n):
        blocks = [np.moveaxis(M.reshape([3] * n + [2]), s, 0).reshape(3, -1) for M in (U[:, :2], Vh[:2].conj().T)]
        T = np.concatenate(blocks, axis=1)
        ev = np.linalg.eigvalsh(T @ T.conj().T)
        out.append(int((ev > tol * ev.max()).sum()))
    return out, S[:3] / S[0]


def cut_ratio(C, I):
    T = C.reshape([3] * n + [3] * n)
    J = [s for s in range(n) if s not in I]
    M = T.transpose(list(I) + [n + s for s in I] + J + [n + s for s in J]).reshape(9 ** len(I), 9 ** len(J))
    sv = np.linalg.svd(M, compute_uv=False)
    return sv[1] / sv[0]


rng = np.random.default_rng(seed)
vals, kinds = [], []
for r in range(runs):
    params = [torch.tensor(rng.standard_normal((d, 2)), requires_grad=True) for _ in range(4)]
    opt = torch.optim.LBFGS(params, lr=1, max_iter=4000, tolerance_grad=1e-13, tolerance_change=1e-16,
                            history_size=50, line_search_fn='strong_wolfe')

    def closure():
        opt.zero_grad()
        v, _ = value(params)
        v.backward()
        return v
    for _ in range(3):
        opt.step(closure)
    v, C = value(params)
    C = C.detach().numpy()
    jr, sv = joint_ranks(C)
    best = min(((I, cut_ratio(C, I)) for k in range(1, n) for I in itertools.combinations(range(n), k)),
               key=lambda z: z[1], default=((), 1.0))
    qubit = all(x == 2 for x in jr)
    kind = 'local-qubit' if qubit else 'OTHER'
    if not qubit and best[1] < 1e-4:
        for I in (best[0], tuple(s for s in range(n) if s not in best[0])):   # either side of the cut may be the tail
            if all(jr[s] == 2 for s in range(n) if s not in I):
                kind = 'rank-one tail on %s' % (I,)
                break
    vals.append(v.item()); kinds.append(kind)
    print(f"run {r:3d}: q/|C|^2 = {v.item(): .2e}  s/s1 = {np.round(sv, 4)}  joint local ranks {jr}  "
          f"best cut {best[0]} ratio {best[1]:.1e}  -> {kind}", flush=True)
print(f"n = {n}: {runs} runs, largest |q|/|C|^2 = {max(abs(x) for x in vals):.1e}; "
      f"local-qubit {kinds.count('local-qubit')}, rank-one tail {sum(k.startswith('rank-one') for k in kinds)}, "
      f"other {kinds.count('OTHER')}")
