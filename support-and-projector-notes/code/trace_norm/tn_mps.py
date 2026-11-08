"""Section 9: the trace-norm bound on matrix-product inputs (any number of copies).
alpha = sum_i |i>_R (x) MPS(l_i; A_1 .. A_n; r) and beta likewise, C = sum_i a_i b_i^* (rank <= 2);
q_n(t, C) = <alpha (x) beta| F_R (x) prod_s (I - t F_s)|alpha (x) beta> by transfer matrices; ||C||^2 and s_1 s_2 from the
Gram matrices of the reference marginals.  Checked against dense evaluation at small n.
usage: python tn_mps.py n chi t runs seed mode
mode 'norm': minimise [q - RHS]/||C||^2 ; mode 'ss': minimise [q - RHS]/(s1 s2) (keeps away from rank one)."""
import sys, numpy as np, torch
torch.set_default_dtype(torch.float64)
n, chi, t, runs, seed, mode = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
d = 3

def step(L, As, Bs, swap):
    T1 = torch.einsum('acbd,xae->xcbde', L, As)
    T2 = torch.einsum('xcbde,ybg->xycdeg', T1, Bs)
    if swap:
        T2 = T2 - t * T2.transpose(0, 1)
    T4 = torch.einsum('uvcdeg,ucf->vdegf', T2, As.conj())
    return torch.einsum('vdegf,vdh->efgh', T4, Bs.conj())

def contract(A, B, l, m, r, rp, swap):
    L = torch.einsum('ia,jc,jb,id->acbd', l, l.conj(), m, m.conj())
    for s in range(n):
        L = step(L, A[s], B[s], swap)
    return torch.einsum('efgh,e,f,g,h->', L, r, r.conj(), rp, rp.conj())

def gram(A, l, r):
    L = torch.einsum('ia,jc->ijac', l.conj(), l)
    for s in range(n):
        L = torch.einsum('ijac,xae,xcf->ijef', L, A[s].conj(), A[s])
    return torch.einsum('ijef,e,f->ij', L, r.conj(), r)

def unpack(ps):
    c = [torch.complex(ps[2 * i], ps[2 * i + 1]) for i in range(len(ps) // 2)]
    return c[:n], c[n:2 * n], c[2 * n], c[2 * n + 1], c[2 * n + 2], c[2 * n + 3]

def evaluate(ps):
    A, B, l, m, r, rp = unpack(ps)
    q = contract(A, B, l, m, r, rp, True).real
    n2 = contract(A, B, l, m, r, rp, False).real
    GA, GB = gram(A, l, r), gram(B, m, rp)
    dA = (GA[0, 0] * GA[1, 1] - GA[0, 1] * GA[1, 0]).real
    dB = (GB[0, 0] * GB[1, 1] - GB[0, 1] * GB[1, 0]).real
    s1s2 = torch.sqrt(torch.clamp(dA * dB, min=1e-30))
    rhs = (1 - t) ** (n - 1) * (n2 - t * (n2 + 2 * s1s2))
    return q - rhs, n2, s1s2, q

def objective(ps):
    g, n2, ss, q = evaluate(ps)
    return g / n2 if mode == 'norm' else g / ss

rng = np.random.default_rng(seed)
def init():
    ps = []
    for s in range(2 * n):
        for _ in range(2):
            ps.append(torch.tensor(rng.standard_normal((d, chi, chi)) / np.sqrt(chi), requires_grad=True))
    for shape in ((2, chi), (2, chi), (chi,), (chi,)):
        for _ in range(2):
            ps.append(torch.tensor(rng.standard_normal(shape), requires_grad=True))
    return ps

print(f"n={n} chi={chi} t={t} mode={mode}", flush=True)
best = 1e9
for rr in range(runs):
    ps = init()
    opt = torch.optim.LBFGS(ps, lr=1, max_iter=800, tolerance_grad=1e-12, tolerance_change=1e-15, history_size=50, line_search_fn='strong_wolfe')
    def cl():
        opt.zero_grad(); v = objective(ps); v.backward(); return v
    for _ in range(3):
        try:
            opt.step(cl)
        except Exception as e:
            print("  step failed:", e); break
    with torch.no_grad():
        g, n2, ss, q = evaluate(ps)
    v = (g / n2).item(); w = (g / ss).item()
    if not (np.isfinite(v) and np.isfinite(w)):
        print(f"run {rr}: non-finite, skipped", flush=True); continue
    ratio = (ss / n2).item()   # s1 s2 / (s1^2 + s2^2) <= 1/2
    best = min(best, v)
    print(f"run {rr}: gap/||C||^2 {v:.3e}  gap/(s1 s2) {w:.3e}  s1s2/||C||^2 {ratio:.3f}  q/||C||^2 {(q/n2).item():.3e}", flush=True)
print(f"smallest gap/||C||^2 {best:.3e}")
