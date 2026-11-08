"""Sharp rank-two bound at odd n with autograd: minimise
   def = [2^n q_n(C) - 2 Tr C^2 + (Tr C)^2] / (Tr C)^2,  C = lam P_u + mu P_v (u, v orthonormal, complex),
   q_n(C) = sum_T (-1/2)^{|T|} ||Tr_T C||_2^2.  Negative would refute the sharp bound (not the endpoint)."""
import itertools, sys, time
import numpy as np, torch
torch.set_default_dtype(torch.float64)
n = int(sys.argv[1]); d = 3; N = d ** n
subsets = [S for k in range(n + 1) for S in itertools.combinations(range(n), k)]
letters = 'abcdefghijklmnop'
def ptr_norm2(C, T):
    # ||Tr_T C||^2 via einsum: C has indices (i_1..i_n, j_1..j_n); trace over T sets i_t = j_t
    keep = [s for s in range(n) if s not in T]
    ins = [letters[s] for s in range(n)]
    outs = [letters[s] if s in T else letters[s].upper() for s in range(n)]
    res = ''.join(letters[s] for s in keep) + ''.join(letters[s].upper() for s in keep)
    M = torch.einsum(''.join(ins) + ''.join(outs) + '->' + res, C.reshape([d] * (2 * n)))
    return (M.abs() ** 2).sum()
def defect(p):
    u = torch.complex(p[0], p[1]); v = torch.complex(p[2], p[3])
    u = u / torch.linalg.norm(u); v = v - (u.conj() @ v) * u; v = v / torch.linalg.norm(v)
    lam, mu = p[4][0] ** 2, p[4][1] ** 2; s = lam + mu
    C = lam * torch.outer(u, u.conj()) + mu * torch.outer(v, v.conj())
    q = sum((-0.5) ** len(T) * ptr_norm2(C, T) for T in subsets)
    return (2 ** n * q - 2 * (lam ** 2 + mu ** 2) + s ** 2) / s ** 2
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
t0 = time.time(); worst = np.inf
for st in range(int(sys.argv[3]) if len(sys.argv) > 3 else 6):
    init = [rng.standard_normal(N) for _ in range(4)] + [rng.standard_normal(2)]
    if st % 3 == 1:   # near qubit-supported vectors
        P = np.zeros(d); P[:2] = 1
        mask = np.array([1.0])
        for s in range(n):
            mask = np.kron(mask, P)
        init = [x * (mask + 0.05) for x in init[:4]] + [init[4]]
    if st % 3 == 2:   # near a product pair
        a = np.array([1.0]); b = np.array([1.0])
        for s in range(n):
            e = rng.standard_normal(3); f = rng.standard_normal(3); f -= (f @ e) / (e @ e) * e
            a = np.kron(a, e); b = np.kron(b, f)
        init = [a + 0.1 * init[0] / np.sqrt(N), 0.1 * init[1] / np.sqrt(N), b + 0.1 * init[2] / np.sqrt(N), 0.1 * init[3] / np.sqrt(N), init[4]]
    p = [torch.tensor(x, requires_grad=True) for x in init]
    opt = torch.optim.LBFGS(p, lr=1, max_iter=4000, history_size=50, line_search_fn='strong_wolfe', tolerance_grad=1e-13, tolerance_change=1e-16)
    def closure():
        opt.zero_grad(); f = defect(p); f.backward(); return f
    opt.step(closure)
    with torch.no_grad():
        val = float(defect(p))
    worst = min(worst, val)
    print(f'start {st}: min normalised defect {val:+.6e} [{time.time() - t0:.0f}s]', flush=True)
print(f'n = {n}: smallest sharp-bound defect found {worst:+.6e} (negative would refute the sharp bound)')
