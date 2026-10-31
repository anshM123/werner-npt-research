"""Section 9: invariant cones and the reduction of 2m copies to m + 1.

   (a) Proposition 4 with site-dependent, rectangular tensors.  A random rank-two C on n copies is written as an MPS
       by successive SVDs (range and co-range families, common copy order); the bond dimensions are printed (for
       n = 2m they never exceed 3^m) and q_n(1/2, C) = <omega|X_0^G|omega> is checked against the dense q_n, with
       X_n = 1 and X_{s-1} = T_s(X_s).
   (b) Evidence for KL_1: least-squares fits of 3 x 3 local filterings of rho^{(x)2} by
       sum_j (a_j (x) b_j) rho (a_j (x) b_j)^* + sum_k |s_k t_k><s_k t_k|, against an in-cone and an out-of-cone control.
       Random filterings are separable (products alone fit them), so the targets are NPT filterings taken along a
       path that lowers lambda_min(Y^G)/Tr Y; for each, products alone and a single filtering of rho are also fitted.
   (c) Evidence for KL_2: the same for 9 x 9 local filterings of rho^{(x)3} by filterings of rho^{(x)2}, with a
       2-distillable control and the PPT entangled Tiles state (its range has no product vectors, so it is outside
       the cone by the range criterion).
   (d) The dimension bound is needed: rho (x) rho as a 9 x 9 operator against the cone of 9 x 9 filterings of one rho.
   (e) A Schmidt-number check for KL_1: every element of K_1 has Schmidt number <= 2, hence Y <= 2 Tr_B(Y) (x) I; the
       smallest eigenvalue of 2 Tr_B(Y) (x) I - Y over 3 x 3 filterings of rho^{(x)2} is minimised.
   (f) Calibration of the fits: a single generator of K_2 (in the cone exactly) against a large model started at
       random and against a one-term model. Large random-start models leave residuals of order 1e-4 on such extreme
       targets, so residuals of that size do not indicate a gap.
   (b)-(f) need PyTorch and take about half an hour; residuals are relative Frobenius norms. Fits are numerical
   evidence, not proofs."""
import numpy as np
from werner import d, q

rng = np.random.default_rng(9)
W4 = (np.eye(d * d) - 0.5 * np.eye(d * d)[[d * (k % d) + k // d for k in range(d * d)]]).reshape(d, d, d, d)


def cr(*s):
    return rng.standard_normal(s) + 1j * rng.standard_normal(s)


def mps(Y, tol=1e-10):
    """Y: 2 x 3^n (rows y_a). Returns l (2 x D_0), tensors A_s (3, D_{s-1}, D_s), r (D_n), bond dimensions."""
    n = round(np.log(Y.shape[1]) / np.log(d))
    U, S, Vh = np.linalg.svd(Y, full_matrices=False); D = int((S > tol * S[0]).sum())
    l, R, As, dims = U[:, :D] * S[:D], Vh[:D], [], [D]
    for s in range(n):
        U, S, Vh = np.linalg.svd(R.reshape(D * d, -1), full_matrices=False); Dn = int((S > tol * S[0]).sum())
        As.append(U[:, :Dn].reshape(D, d, Dn).transpose(1, 0, 2)); R = S[:Dn, None] * Vh[:Dn]; D = Dn; dims.append(D)
    return l, As, R[:, 0], dims


def transfer(A, B, X):
    """T(X) = sum W[(u v),(x y)] (A^u (x) B^v) X (A^x (x) B^y)^*, X on C^D (x) C^D' as a (D D') x (D D') matrix"""
    D, Dp = A.shape[2], B.shape[2]
    Xt = X.reshape(D, Dp, D, Dp)
    out = np.einsum('uvxy,uik,vjl,klnp,xmn,yop->ijmo', W4, A, B, Xt, A.conj(), B.conj(), optimize=True)
    return out.reshape(A.shape[1] * B.shape[1], A.shape[1] * B.shape[1])


def pt2(X, D, Dp):
    return X.reshape(D, Dp, D, Dp).transpose(0, 3, 2, 1).reshape(D * Dp, D * Dp)


print("(a) site-dependent transfer identity:")
for n in (2, 3, 4, 5):
    Y, Z = cr(2, d ** n), cr(2, d ** n)
    C = Y.T @ Z.conj()                          # C = sum_a |y_a><z_a|
    l, As, r, dA = mps(Y); m, Bs, rp, dB = mps(Z)
    X = np.outer(np.kron(r, rp), np.kron(r, rp).conj())
    for s in reversed(range(n)):
        X = transfer(As[s], Bs[s], X)
    omega = sum(np.kron(l[a].conj(), m[a]) for a in range(2))
    val = np.vdot(omega, pt2(X, l.shape[1], m.shape[1]) @ omega).real
    print(f"  n = {n}: bond dimensions {dA} (range), {dB} (co-range), max {max(dA + dB)}; "
          f"transfer {val:+.12e}, dense q_n {q(C, n, 0.5):+.12e}", flush=True)

# (b), (c): least-squares decompositions (PyTorch, L-BFGS). Model terms lie in the cone exactly, so a small residual
# is an explicit approximate decomposition; it is not a certificate of membership.
import torch
torch.manual_seed(0)
F9 = np.eye(d * d)[[d * (k % d) + k // d for k in range(d * d)]]
rho = np.eye(d * d) - 0.5 * F9


def tens(k):
    """rho^{(x)k} with the Alice copies first: (A_1..A_k)(B_1..B_k)"""
    R = np.array([[1.0]])
    for _ in range(k):
        R = np.kron(R, rho)
    perm = list(range(0, 2 * k, 2)) + list(range(1, 2 * k, 2))
    return R.reshape([d] * (4 * k)).transpose(perm + [p + 2 * k for p in perm]).reshape(d ** (2 * k), d ** (2 * k))


def fit(Y, base, D, mlf, mpr, iters):
    """min || sum_j (A_j (x) B_j) base (A_j (x) B_j)^* + sum_k |s_k t_k><s_k t_k| - Y ||, A_j, B_j : C^{dim} -> C^D"""
    din = round(np.sqrt(base.shape[0]))
    R0 = torch.tensor(base, dtype=torch.complex128); Yt = torch.tensor(Y, dtype=torch.complex128)
    P = [(0.3 * torch.randn(*s, dtype=torch.complex128)).requires_grad_()
         for s in ((max(mlf, 1), D, din), (max(mlf, 1), D, din), (max(mpr, 1), D), (max(mpr, 1), D))]

    def model():
        M = torch.zeros((D * D, D * D), dtype=torch.complex128)
        if mlf:
            K = torch.einsum('jab,jcd->jacbd', P[0], P[1]).reshape(mlf, D * D, din * din)
            M = M + (K @ R0 @ K.conj().transpose(1, 2)).sum(0)
        if mpr:
            V = torch.einsum('ka,kb->kab', P[2], P[3]).reshape(mpr, D * D)
            M = M + V.T @ V.conj()
        return M
    opt = torch.optim.LBFGS(P, lr=1, max_iter=iters, history_size=50, line_search_fn='strong_wolfe',
                            tolerance_grad=1e-14, tolerance_change=1e-16)

    def closure():
        opt.zero_grad(); loss = torch.linalg.norm(model() - Yt) ** 2; loss.backward(); return loss
    opt.step(closure)
    with torch.no_grad():
        return float(torch.linalg.norm(model() - Yt) / torch.linalg.norm(Yt))


def lf(base, D, k=1, prods=0):
    """random sum of k local filterings of base into D x D (plus prods product terms), unit trace"""
    din = round(np.sqrt(base.shape[0])); M = 0
    for _ in range(k):
        K = np.kron(cr(D, din), cr(D, din)); M = M + K @ base @ K.conj().T
    for _ in range(prods):
        v = np.kron(cr(D), cr(D)); M = M + 0.3 * np.outer(v, v.conj())
    return M / np.trace(M).real


def pt_t(Y, D):
    return Y.reshape(D, D, D, D).permute(0, 3, 2, 1).reshape(D * D, D * D)


def npt_targets(base, din, D, thresholds, seed, lr=0.01, steps=20000):
    """filterings of base into D x D along an Adam path that lowers lambda_min(Y^G)/Tr Y, kept at the thresholds"""
    torch.manual_seed(seed); Rt = torch.tensor(base, dtype=torch.complex128); out = []
    a = torch.randn(D, din, dtype=torch.complex128, requires_grad=True)
    b = torch.randn(D, din, dtype=torch.complex128, requires_grad=True)
    opt = torch.optim.Adam([a, b], lr=lr)
    for it in range(steps):
        K = torch.kron(a, b); Y = K @ Rt @ K.conj().T; Y = (Y + Y.conj().T) / 2; Y = Y / torch.trace(Y).real
        f = torch.linalg.eigvalsh(pt_t(Y, D))[0]
        if f.item() < thresholds[len(out)]:
            out.append((f.item(), Y.detach().numpy()))
            if len(out) == len(thresholds):
                break
        opt.zero_grad(); f.backward(); opt.step()
    return out


Pa = (np.eye(d * d) - F9) / 2
print("\n(b) KL_1 (two copies -> 3 x 3; KL_1 model = 6 filterings of rho + 12 products):")
r1, r2 = tens(1), tens(2)
print(f"  in-cone control (sum of 4 filterings of rho): {fit(lf(r1, 3, 4), r1, 3, 6, 12, 3000):.1e}")
print(f"  out-of-cone control (antisymmetric projector, 1-distillable): {fit(Pa / 3, r1, 3, 6, 12, 3000):.1e}")
print(f"  a random filtering of rho^(x)2 is separable: products only (20) {fit(lf(r2, 3), r1, 3, 0, 20, 3000):.1e}")
for f, Y in npt_targets(r2, 9, 3, [-0.01, -0.03], seed=100):
    print(f"  NPT filtering, lambda_min(Y^G)/Tr Y = {f:+.4f}: products only (20) {fit(Y, r1, 3, 0, 20, 3000):.1e}; "
          f"one filtering of rho {fit(Y, r1, 3, 1, 0, 3000):.1e}; KL_1 model {fit(Y, r1, 3, 6, 12, 3000):.1e}", flush=True)

print("\n(c) KL_2 (three copies -> 9 x 9; KL_2 model = 96 filterings of rho^(x)2 + 240 products):")
r3 = tens(3)
Yo = np.kron(Pa, np.diag([1.0, 1, 1, 0, 0, 0, 0, 0, 0])).reshape([3] * 8).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(81, 81)
print(f"  out-of-cone control (P_a (x) P, 2-distillable): {fit(Yo / np.trace(Yo), r2, 9, 96, 240, 8000):.1e}", flush=True)
ket = lambda i: np.eye(3)[i]
ups = [np.kron(ket(0), ket(0) - ket(1)), np.kron(ket(2), ket(1) - ket(2)), np.kron(ket(0) - ket(1), ket(2)),
       np.kron(ket(1) - ket(2), ket(0)), np.kron(np.ones(3), np.ones(3)) / 1.5]
tiles = np.eye(9) - sum(np.outer(u, u) / (u @ u) for u in ups)
e00 = np.zeros((3, 3)); e00[0, 0] = 1
Yt = np.kron(tiles, np.kron(e00, e00)).reshape([3] * 8).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(81, 81)
print(f"  Tiles (x) |00><00| (PPT entangled): {fit(Yt / np.trace(Yt), r2, 9, 96, 240, 8000):.1e}", flush=True)
print(f"  a random filtering of rho^(x)3 is separable: products only (700) {fit(lf(r3, 9), r2, 9, 0, 700, 20000):.1e}",
      flush=True)
for f, Y in npt_targets(r3, 27, 9, [-0.003, -0.015], seed=200):
    print(f"  NPT filtering, lambda_min(Y^G)/Tr Y = {f:+.4f}: products only (400) {fit(Y, r2, 9, 0, 400, 8000):.1e}; "
          f"one filtering of rho^(x)2 {fit(Y, r2, 9, 1, 0, 3000):.1e}; KL_2 model {fit(Y, r2, 9, 96, 240, 8000):.1e}",
          flush=True)

print("\n(d) without the dimension bound (9 x 9 filterings of one rho; 60 filterings + 120 products):")
print(f"  in-cone control: {fit(lf(r1, 9, 30, 60), r1, 9, 60, 120, 6000):.1e}", flush=True)
print(f"  rho (x) rho itself: {fit(r2 / np.trace(r2), r1, 9, 60, 120, 6000):.1e}", flush=True)

print("\n(e) Schmidt-number check for KL_1: min over filterings of lambda_min(2 Tr_B(Y) (x) I - Y) / Tr Y")
R2t = torch.tensor(r2, dtype=torch.complex128); I3 = torch.eye(3, dtype=torch.complex128); best = np.inf
for s in range(4):
    a = torch.randn(3, 9, dtype=torch.complex128, requires_grad=True); b = torch.randn(3, 9, dtype=torch.complex128, requires_grad=True)
    opt = torch.optim.Adam([a, b], lr=0.03)
    for it in range(2000):
        K = torch.kron(a, b); Y = K @ R2t @ K.conj().T; Y = (Y + Y.conj().T) / 2
        YA = torch.einsum('ijkj->ik', Y.reshape(3, 3, 3, 3))
        f = torch.linalg.eigvalsh(2 * torch.kron(YA, I3) - Y)[0] / torch.trace(Y).real
        opt.zero_grad(); f.backward(); opt.step()
    best = min(best, f.item())
print(f"  smallest value found: {best:+.2e} (a negative value would refute KL_1)")

print("\n(f) Calibration: the error floor of the fits on extreme targets. A single 6 x 6 filtering of rho^(x)2 is a generator")
print("    of K_2, so it is in the cone exactly. Large models started at random fit it poorly; a one-term model fits it.")
Yg = lf(r2, 6)
print(f"  48 filterings + 96 products, random start: {fit(Yg, r2, 6, 48, 96, 6000):.1e}; one filtering: {fit(Yg, r2, 6, 1, 0, 3000):.1e}",
      flush=True)
