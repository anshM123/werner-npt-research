"""Five copies: exclusive antisymmetric-sector weights p_r(u,v) = sum_{|S|=r} ||Pi_S (u (x) v)||^2 and a_r(w) = p_r(w,w).
Level-wise Wu-Zou-type bounds:  L1 = (3 p3 - p2) - sqrt(a2(u) a2(v)),  L2 = (3 p5 - p4) - sqrt(a4(u) a4(v)).
Coupled (needed) bound: B5 - sqrt(A5(u) A5(v)) with B5 = (3p3 - p2) + 10 (3p5 - p4), A5 = a2 + 10 a4.
Maximise each over orthonormal u, v (complex) by gradient ascent from several starts."""
import itertools, sys, time
import numpy as np, torch
torch.set_default_dtype(torch.float64)
n = 5; d = 3; N = d ** n
subsets = [frozenset(S) for k in range(n + 1) for S in itertools.combinations(range(n), k)]
letters = 'abcdefghij'
def marg(P, K):
    keep = sorted(K)
    ins = ''.join(letters[s] for s in range(n))
    outs = ''.join(letters[s] if s not in K else letters[s].upper() for s in range(n))
    res = ''.join(letters[s] for s in keep) + ''.join(letters[s].upper() for s in keep)
    return torch.einsum(ins + outs + '->' + res, P.reshape([d] * (2 * n))).reshape(d ** len(K), d ** len(K))
def sector_weights(u, v):
    Pu, Pv = torch.outer(u, u.conj()), torch.outer(v, v.conj())
    tr = {}
    for K in subsets:
        if len(K) == 0:
            tr[K] = (u.conj() @ u).real * (v.conj() @ v).real
        else:
            tr[K] = torch.real(torch.trace(marg(Pu, K) @ marg(Pv, K)))
    incl = {R: sum(((-1) ** len(K)) * tr[K] for K in subsets if K <= R) / 2 ** len(R) for R in subsets}
    excl = {S: sum(((-1) ** len(R - S)) * incl[R] for R in subsets if R >= S) for S in subsets}
    return [sum(excl[S] for S in subsets if len(S) == r) for r in range(n + 1)]
def quantities(p):
    u = torch.complex(p[0], p[1]); v = torch.complex(p[2], p[3])
    u = u / torch.linalg.norm(u); v = v - (u.conj() @ v) * u; v = v / torch.linalg.norm(v)
    pr = sector_weights(u, v); au = sector_weights(u, u); av = sector_weights(v, v)
    L1 = (3 * pr[3] - pr[2]) - torch.sqrt(torch.clamp(au[2] * av[2], min=1e-30))
    L2 = (3 * pr[5] - pr[4]) - torch.sqrt(torch.clamp(au[4] * av[4], min=1e-30))
    B = (3 * pr[3] - pr[2]) + 10 * (3 * pr[5] - pr[4]); A_u = au[2] + 10 * au[4]; A_v = av[2] + 10 * av[4]
    Cpl = B - torch.sqrt(torch.clamp(A_u * A_v, min=1e-30))
    return L1, L2, Cpl, pr, au, av
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
t0 = time.time()
best = {'L1': -np.inf, 'L2': -np.inf, 'coupled': -np.inf}
for target in ('L1', 'L2', 'coupled'):
    for st in range(int(sys.argv[2]) if len(sys.argv) > 2 else 3):
        p = [torch.tensor(rng.standard_normal(N), requires_grad=True) for _ in range(4)]
        if st % 3 == 1:   # near the anti-state family: v = (i sigma_y)^{tensor} conj(u) on qubit planes
            pass
        opt = torch.optim.Adam(p, lr=0.02)
        for it in range(600):
            L1, L2, Cpl, *_ = quantities(p)
            f = {'L1': L1, 'L2': L2, 'coupled': Cpl}[target]
            opt.zero_grad(); (-f).backward(); opt.step()
        with torch.no_grad():
            L1, L2, Cpl, pr, au, av = quantities(p)
        val = {'L1': float(L1), 'L2': float(L2), 'coupled': float(Cpl)}[target]
        best[target] = max(best[target], val)
        print(f'maximise {target}, start {st}: {val:+.5f} (L1 {float(L1):+.4f}, L2 {float(L2):+.4f}, coupled {float(Cpl):+.4f}) [{time.time() - t0:.0f}s]', flush=True)
print('largest values found (positive = violation):', {k: round(v, 5) for k, v in best.items()})
