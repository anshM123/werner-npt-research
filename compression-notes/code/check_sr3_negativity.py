"""Section 11. Schmidt-rank-3 negativity of (rho^Gamma)^{(x)m}: minimise <w|(rho^Gamma)^{(x)m}|w> over unit w of Schmidt rank <= 3
(w = sum_{s=1..3} y_s (x) z_s, Alice/Bob = m qutrits each).  Report the local minima found and, for each, the per-copy
'Phi-content' p_i = <w| Pi_i |w> (Pi_i = projector onto the maximally entangled state of copy pair i) and the overlap with
the best 'single-copy' structure Phi_i (x) product.  Also: values along random SR-3 directions (how rare is negativity)."""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
m = int(sys.argv[1]); RS = int(sys.argv[2]) if len(sys.argv) > 2 else 20
N = 3 ** m
subsets = [S for r in range(m + 1) for S in itertools.combinations(range(m), r)]
def value(Y, Z):
    """<w|(rho^Gamma)^m|w> for w = sum_s y_s (x) z_s (columns of Y, Z: N x 3) = sum_J (-1/2)^{|J|} ||Tr_J C||^2, C = Y Z^T."""
    C = Y @ Z.T
    Ct = C.reshape([3] * m + [3] * m)
    tot = 0
    for J in subsets:
        idx = list(range(2 * m))
        # partial trace over copies in J: contract index j with m + j
        letters = 'abcdefghijklmnopqrstuvwxyz'
        inp = [letters[i] for i in range(2 * m)]
        for j in J: inp[m + j] = inp[j]
        out = [inp[i] for i in range(2 * m) if not (i in J or (i >= m and i - m in J))]
        T = torch.einsum(''.join(inp) + '->' + ''.join(out), Ct)
        tot = tot + (-0.5) ** len(J) * (T.abs() ** 2).sum()
    return tot
def phi_content(Y, Z):
    C = (Y @ Z.T).reshape([3] * m + [3] * m); nrm = (C.abs() ** 2).sum(); out = []
    for i in range(m):
        letters = 'abcdefghijklmnopqrstuvwxyz'
        inp = [letters[k] for k in range(2 * m)]; inp[m + i] = inp[i]
        o = [inp[k] for k in range(2 * m) if k not in (i, m + i)]
        T = torch.einsum(''.join(inp) + '->' + ''.join(o), C)
        out.append(float((T.abs() ** 2).sum() / nrm / 3))
    return out
torch.manual_seed(0); np.random.seed(0)
mins = []
for r in range(RS):
    Y = torch.randn(N, 3, dtype=torch.complex128, requires_grad=True); Z = torch.randn(N, 3, dtype=torch.complex128, requires_grad=True)
    opt = torch.optim.Adam([Y, Z], lr=0.02)
    for it in range(3000):
        C = Y @ Z.T; v = value(Y, Z) / (C.abs() ** 2).sum()
        opt.zero_grad(); v.backward(); opt.step()
    with torch.no_grad():
        C = Y @ Z.T; v = float(value(Y, Z) / (C.abs() ** 2).sum()); pc = phi_content(Y, Z)
    mins.append(v)
    print(f"m = {m}, restart {r}: min {v:+.5f}; per-copy Phi-content {np.round(pc, 4)}", flush=True)
# random SR-3 vectors: fraction negative
vals = []
for t in range(2000):
    Y = torch.randn(N, 3, dtype=torch.complex128); Z = torch.randn(N, 3, dtype=torch.complex128)
    with torch.no_grad():
        C = Y @ Z.T; vals.append(float(value(Y, Z) / (C.abs() ** 2).sum()))
vals = np.array(vals)
print(f"m = {m}: random SR-3 vectors: min {vals.min():+.4f}, fraction negative {np.mean(vals < 0):.4f}; "
      f"global minimum eigenvalue of (rho^Gamma)^m: -1/2 (Phi on one copy, product states orthogonal to Phi elsewhere)",
      flush=True)
