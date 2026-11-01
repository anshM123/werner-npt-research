"""Section 11, remark on the k-family.  Validated test of R(k): min over unit psi, phi in C^k (x) (C^d)^{(x)n} of <phi| Omega_t(psi) |phi>,
Omega_t(psi) = sum_J (-t)^{|J|} psi_{qJ} (x) I_{J^c}.  Smooth bilinear objective (no eigvalsh), random and structured
starts (maximally entangled on one copy (x) products).  Control t = 1/k + 0.02 must come out negative.
usage: python check_k_family.py d k n restarts   (needs PyTorch)"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
d, k, n, RS = (int(a) for a in sys.argv[1:5])
N = d ** n
letters = 'abcdefghij'
def form(psi, phi, t):
    """<phi| Omega_t(psi) |phi> = sum_J (-t)^{|J|} Tr[psi_{qJ} phi_{qJ}]"""
    P = psi.reshape([k] + [d] * n); Q = phi.reshape([k] + [d] * n); tot = 0
    for m in range(n + 1):
        for J in itertools.combinations(range(n), m):
            # Tr[psi_{qJ} phi_{qJ}] = sum over (q, J) indices of (sum_{J^c} P Q^*) (sum_{J^c} P^* Q)
            ik = [letters[j] for j in range(n)]
            jk = [letters[j].upper() if j in J else letters[j] for j in range(n)]
            red_p = torch.einsum('q' + ''.join(ik) + ',Q' + ''.join(jk) + '->qQ' + ''.join(letters[j] + letters[j].upper() for j in J), P, P.conj())
            red_f = torch.einsum('q' + ''.join(ik) + ',Q' + ''.join(jk) + '->qQ' + ''.join(letters[j] + letters[j].upper() for j in J), Q, Q.conj())
            # Tr[A B] with A_{(q,J),(Q,J')} : contract A[q,Q,j,J] with B[Q,q,J,j]
            Ap = red_p.reshape(k, k, -1); Bf = red_f.reshape(k, k, -1)
            # need index swap within each (j, J) pair for B: build via permutation
            if m > 0:
                Bf = red_f.reshape([k, k] + [d, d] * m)
                perm = [1, 0] + sum([[2 + 2 * i + 1, 2 + 2 * i] for i in range(m)], [])
                Bf = Bf.permute(perm).reshape(k, k, -1)
                tot = tot + (-t) ** m * (Ap * Bf).sum().real
            else:
                tot = tot + (Ap * Bf.permute(1, 0, 2)).sum().real
    return tot
def run(t, seed):
    torch.manual_seed(seed); best = np.inf
    starts = []
    for s in range(RS):
        starts.append((torch.randn(k * N, dtype=torch.complex128), torch.randn(k * N, dtype=torch.complex128)))
    # structured: psi = max-entangled(k) on copy 0 (x) e_0 on the others; phi = psi
    v = torch.zeros(k, d, dtype=torch.complex128)
    for i in range(k): v[i, i] = 1
    rest = torch.zeros(d ** (n - 1), dtype=torch.complex128); rest[0] = 1
    ps = torch.einsum('qa,r->qar', v, rest).reshape(k * N)
    starts.append((ps + 0.05 * torch.randn(k * N, dtype=torch.complex128), ps + 0.05 * torch.randn(k * N, dtype=torch.complex128)))
    for (p0, f0) in starts:
        pr, pi = p0.real.clone().requires_grad_(), p0.imag.clone().requires_grad_()
        fr, fi = f0.real.clone().requires_grad_(), f0.imag.clone().requires_grad_()
        o = torch.optim.Adam([pr, pi, fr, fi], lr=0.02)
        for it in range(1500):
            psi = torch.complex(pr, pi); psi = psi / torch.linalg.norm(psi)
            phi = torch.complex(fr, fi); phi = phi / torch.linalg.norm(phi)
            val = form(psi, phi, t); o.zero_grad(); val.backward(); o.step()
        best = min(best, val.item())
    return best
for t in [1.0 / k + 0.02, 1.0 / k]:
    print(f"d = {d}, k = {k}, n = {n}, t = {t:.4f} ({'control, must be < 0' if t > 1.0 / k else 'claim'}): min <phi|Omega|phi> = {run(t, 0):+.3e}", flush=True)
