"""Counterexample search on code subspaces with k logical qutrits (general stabilizer codes, random symplectic bases).
h_L = |S|^{-1} sum_{E in L+S} (-1/2)^{#identities(E)};  q_N(Enc C)/||C||^2 = sum_L h_L |Tr(W_L^dag C)|^2 / (3^k ||C||^2).
Annealing over codes (random transvections), each scored by the minimum over rank-two logical C (LBFGS restarts).
Checks before scoring: Parseval and rank-one positivity on the logical transform.
usage: python k_codes.py N k ncodes steps seed"""
import sys, itertools, time
import numpy as np, torch
torch.set_default_dtype(torch.float64)
N, k, ncodes, steps, seed = (int(a) for a in sys.argv[1:6])
rng = np.random.default_rng(seed); torch.manual_seed(seed)
D = 3 ** k
def transvect(M, v, c):
    s = (M[:, :N] @ v[N:] - M[:, N:] @ v[:N]) % 3
    return (M + c * np.outer(s, v)) % 3
def random_sp(nsteps, M=None):
    M = np.eye(2 * N, dtype=np.int64) if M is None else M
    for _ in range(nsteps):
        v = rng.integers(0, 3, 2 * N)
        if v.any(): M = transvect(M, v, int(rng.integers(1, 3)))
    return M
coeffs_S = np.array(list(itertools.product(range(3), repeat=N - k)), dtype=np.int64)
Lco = np.array(list(itertools.product(range(3), repeat=2 * k)), dtype=np.int64)    # (a_1..a_k, b_1..b_k) big-endian
def code_h(M):
    S = (coeffs_S @ M[[N + j for j in range(N - k)]]) % 3
    LX = M[[N - k + j for j in range(k)]]; LZ = M[[2 * N - k + j for j in range(k)]]
    h = np.zeros(len(Lco))
    for i, c in enumerate(Lco):
        L = (c[:k] @ LX + c[k:] @ LZ) % 3
        E = (S + L) % 3
        h[i] = ((-0.5) ** ((E[:, :N] == 0) & (E[:, N:] == 0)).sum(1)).sum() / len(S)
    return h
# logical Weyl transform: c[a, b] = sum_j w^{-b.j} C[j + a, j], indices little-endian base 3
digits = np.array(list(itertools.product(range(3), repeat=k)))[:, ::-1]
plus = np.array([(((digits + digits[a]) % 3) * 3 ** np.arange(k)).sum(1) for a in range(D)])
plus_t = torch.tensor(plus)
om = np.exp(2j * np.pi / 3)
D3t = torch.tensor(np.array([[om ** (-(a * b) % 3) for b in range(3)] for a in range(3)]))
def wcoeffs(C):
    V = torch.gather(C, 0, plus_t.T)                    # V[a, j] = C[j + a, j]
    T = V.reshape([D] + [3] * k)
    for ax in range(1, k + 1):
        T = torch.movedim(torch.tensordot(D3t, T, dims=([1], [ax])), 0, ax)
    return T.reshape(D, D)                              # [a_index (little-endian), b_index (C-order of digit axes)]
# map Lco (big-endian a, big-endian b) to the [a_le, b] layout used by wcoeffs
def to_grid(h):
    G = np.zeros((D, D))
    for i, c in enumerate(Lco):
        a = c[:k]; b = c[k:]
        ai = int((a * 3 ** np.arange(k)).sum())         # little-endian index of a (digit i <-> logical qutrit i)
        bi = int((b[::-1] * 3 ** np.arange(k)).sum()) if False else int((b * 3 ** np.arange(k)).sum())
        # wcoeffs' b index: C-order over the k digit axes of j, which hold j's digits k-1..0; b digit on axis ax pairs
        # with j digit (k-1-ax), so b index = sum_i b_i 3^i with b_i paired to j digit i (little-endian)
        G[ai, bi] = h[i]
    return torch.tensor(G)
Cr = torch.complex(torch.randn(D, D), torch.randn(D, D))
assert abs((wcoeffs(Cr).abs() ** 2).sum().item() / (D * (Cr.abs() ** 2).sum().item()) - 1) < 1e-9, "Parseval"
def logical_min(Hg, restarts):
    best = None
    for _ in range(restarts):
        p = torch.tensor(rng.normal(size=8 * D) / np.sqrt(D), requires_grad=True)
        def f(p):
            a = torch.complex(p[:2 * D].reshape(2, D), p[2 * D:4 * D].reshape(2, D))
            b = torch.complex(p[4 * D:6 * D].reshape(2, D), p[6 * D:].reshape(2, D))
            C = a.T @ b.conj()
            return (Hg * wcoeffs(C).abs() ** 2).sum() / (D * (C.abs() ** 2).sum())
        opt = torch.optim.LBFGS([p], max_iter=150, line_search_fn='strong_wolfe')
        def cl():
            opt.zero_grad(); v = f(p); v.backward(); return v
        try: opt.step(cl)
        except Exception: pass
        v = f(p).item(); best = v if best is None else min(best, v)
    return best
def rank_one_check(Hg):
    vals = []
    for _ in range(3):
        a = torch.complex(torch.randn(D), torch.randn(D)); C1 = torch.outer(a, a.conj())
        vals.append(((Hg * wcoeffs(C1).abs() ** 2).sum() / (D * (C1.abs() ** 2).sum())).item())
    return min(vals)
t0 = time.time(); overall = None
for c in range(ncodes):
    M = random_sp(8 * N); Hg = to_grid(code_h(M))
    r1 = rank_one_check(Hg)
    assert r1 > 0, f"rank-one value {r1} <= 0: transform or weights wrong"
    cur = (logical_min(Hg, 4), M)
    for s in range(steps):
        v = rng.integers(0, 3, 2 * N)
        if not v.any(): continue
        M2 = transvect(cur[1], v, int(rng.integers(1, 3)))
        val = logical_min(to_grid(code_h(M2)), 3)
        T = 0.02 * (1 - s / max(steps, 1))
        if val < cur[0] or rng.random() < np.exp(-(val - cur[0]) / max(T, 1e-9)): cur = (val, M2)
    fin = logical_min(to_grid(code_h(cur[1])), 8)
    overall = fin if overall is None else min(overall, fin)
    print(f"N={N} k={k} code {c}: annealed logical min {fin:+.6e}  ({time.time() - t0:.0f}s)", flush=True)
    if fin < -1e-8:
        np.save(f"candidate_N{N}_k{k}_{c}.npy", cur[1]); print("  NEGATIVE: saved code for exact re-verification", flush=True)
print(f"N={N} k={k}: smallest {overall:+.6e}  ({'NEGATIVE CANDIDATE' if overall < -1e-8 else 'no counterexample'})")
