"""Section 1: averaged compression rules.  For copy s let rho_s be the local support operator of C (reduced operator of
C C^* + C^* C on copy s) and draw w from the angular central Gaussian with shape rho_s^p.  The rule asks for
max_s E[Delta_s(w)] >= 0.  The adversary minimises this over unit rank-two C (expectation over a fixed sample set, soft
maximum over s), then rechecks the worst C with 20000 fresh samples.  For p = -0.5, -1, -2 it finds negative values at
two copies (see the log), so no rule of this form proves QC.
usage: python check_qc_averaged_rule.py n p samples steps seed   (e.g. 2 -1 600 500 5; needs PyTorch, a few minutes)
"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
n, p, NS, STEPS, SEED = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
d = 3; N = d ** n; t = 0.5
subsets = [S for r in range(n + 1) for S in itertools.combinations(range(n), r)]
letters = 'abcdefghijklmnopqrstuvwxyz'
def qn(C):
    """q_n(C) = sum_J (-t)^{|J|} ||Tr_J C||^2, C of shape (N, N); batched over leading dims."""
    lead = C.shape[:-2]; Ct = C.reshape(*lead, *([d] * (2 * n))); tot = 0
    nl = len(lead)
    for J in subsets:
        inp = [letters[i] for i in range(2 * n)]
        for j in J: inp[n + j] = inp[j]
        out = [inp[i] for i in range(2 * n) if not (i in J or (i >= n and i - n in J))]
        pre = 'XYZ'[:nl]
        T = torch.einsum('...' + ''.join(inp) + '->...' + ''.join(out), Ct)
        tot = tot + (-t) ** len(J) * (T.abs() ** 2).reshape(*lead, -1).sum(-1)
    return tot
def local_support(C, s):
    M = C @ C.conj().T + C.conj().T @ C
    T = M.reshape([d] * (2 * n)); inp = [letters[i] for i in range(2 * n)]
    for k in range(n):
        if k != s: inp[n + k] = inp[k]
    return torch.einsum(''.join(inp) + '->' + inp[s] + inp[n + s], T)
def compress(C, s, w):
    """C_{w perp} for a batch of unit w (K, d): P = I - w w^* on copy s, both sides."""
    K = w.shape[0]
    P = torch.eye(d, dtype=C.dtype).expand(K, d, d) - torch.einsum('ki,kj->kij', w, w.conj())
    Ct = C.reshape([d] * (2 * n))
    inp = [letters[i] for i in range(2 * n)]
    # left: P on output index s; right: P on input index n+s
    a, b = inp[s], inp[n + s]
    out = list(inp); out[s] = 'X'; out[n + s] = 'Y'
    expr = 'k' + 'X' + a + ',' + ''.join(inp) + ',' + 'k' + b + 'Y' + '->k' + ''.join(out)
    R = torch.einsum(expr, P, Ct, P)
    return R.reshape(K, N, N)
g_fixed = torch.randn(NS, d, dtype=torch.complex128, generator=torch.Generator().manual_seed(1000 + SEED))
def G(C, g, soft=True):
    C = C / torch.linalg.norm(C)
    q0 = qn(C)
    vals = []
    for s in range(n):
        ev, U = torch.linalg.eigh(local_support(C, s)); ev = torch.clamp(ev, min=1e-12 * ev.max().item())
        S = U @ torch.diag((ev ** (p / 2)).to(C.dtype)) @ U.conj().T
        w = g @ S.T; w = w / torch.linalg.norm(w, dim=1, keepdim=True)
        Cw = compress(C, s, w)
        vals.append((q0 - qn(Cw)).mean())
    vals = torch.stack(vals)
    return (torch.logsumexp(vals * 200, 0) / 200) if soft else vals.max(), q0
torch.manual_seed(SEED)
Yr, Yi, Zr, Zi = [torch.randn(N, 2, requires_grad=True) for _ in range(4)]
opt = torch.optim.Adam([Yr, Yi, Zr, Zi], lr=0.02)
best = (np.inf, None)
for it in range(STEPS):
    C = torch.complex(Yr, Yi) @ torch.complex(Zr, Zi).conj().T
    val, q0 = G(C, g_fixed)
    if val.item() < best[0]: best = (val.item(), C.detach().clone(), q0.item())
    opt.zero_grad(); val.backward(); opt.step()
    if it % 50 == 0 or it == STEPS - 1:
        print(f"n = {n}, p = {p}, seed {SEED}, step {it}: G_p {val.item():+.5f} (q_n {q0.item():+.5f}); best {best[0]:+.5f}", flush=True)
# recheck the worst C with fresh samples
Cw = best[1]
g2 = torch.randn(20000, d, dtype=torch.complex128, generator=torch.Generator().manual_seed(77))
with torch.no_grad():
    vals = []
    for chunk in torch.split(g2, 2000):
        v, q0 = G(Cw, chunk, soft=False); vals.append(v.item())
print(f"n = {n}, p = {p}, seed {SEED}: worst C rechecked with 20000 fresh samples: G_p = {np.mean(vals):+.5f} "
      f"(+- {np.std(vals) / np.sqrt(len(vals)):.5f}), q_n = {q0.item():+.5f}", flush=True)
np.save(f"qc_avg_worst_n{n}_p{p}_s{SEED}.npy", Cw.numpy())
