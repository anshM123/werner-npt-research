"""Why two planar copies on one side cannot share a group.
Take f free copies and one group made of two planar copies s, t (A = A_s (x) A_t, dim 4), and an operator X on H_F (x) A.
The grouping needs 2^m q_n(X) = 2 q_{F+{A}}(X) (m = 2), where A counts as one copy of dimension 4.  Here
    G(X) = 2^m q_n(X) - 2 q_{F+{A}}(X) = 2 sum_{J in F} (-1/2)^|J| [ ||Y||^2 - ||Tr_s Y||^2 - ||Tr_t Y||^2 + ||Tr_st Y||^2 ],
with Y = Tr_J X.  This script checks that identity on random X, computes the smallest eigenvalue of G (it is -2), evaluates
the rank-two operator X = |0..0><1..1|_F (x) (1/sqrt2) 1_{A_s} (x) |a><a'| (value -2), and minimises G over rank-one X
(minimum 0, by the rank-one bound).   usage: python check_group_gap.py f starts seed"""
import sys, numpy as np, torch
torch.set_default_dtype(torch.float64); cplx = torch.complex128
f, ST, SEED = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); torch.manual_seed(SEED)
DF = 3 ** f
def q_dims(C, dims, coef=-0.5):  # sum_J coef^|J| ||Tr_J C||^2 on a product with the given factor dims
    k = len(dims); T = C.reshape(dims + dims); tot = 0.0
    for mask in range(2 ** k):
        J = [i for i in range(k) if mask >> i & 1]; keep = [i for i in range(k) if i not in J]
        letters = 'abcdefghijklmnop'; L = list(letters[:k]); R = [letters[k + i] if i in keep else letters[i] for i in range(k)]
        out = ''.join(L[i] for i in keep) + ''.join(R[i] for i in keep)
        tot = tot + coef ** len(J) * torch.sum(torch.abs(torch.einsum(''.join(L) + ''.join(R) + '->' + out, T)) ** 2)
    return tot
def G(X):  # 2^m q_n(X) - 2 q_{F+{A}}(X), m = 2
    return 4 * q_dims(X, [3] * f + [2, 2]) - 2 * q_dims(X, [3] * f + [4])
def G_closed(X):  # 2 x sum_{J in F} (-1/2)^|J| x (form with coefficient -1 on the two planar copies) of Tr_J X
    k = f + 2; T = X.reshape([3] * f + [2, 2] + [3] * f + [2, 2]); tot = 0.0; letters = 'abcdefghijklmnop'
    for mask in range(2 ** f):
        J = [i for i in range(f) if mask >> i & 1]; keep = [i for i in range(k) if i not in J]
        L = list(letters[:k]); R = [letters[k + i] if i in keep else letters[i] for i in range(k)]
        out = ''.join(L[i] for i in keep) + ''.join(R[i] for i in keep)
        Y = torch.einsum(''.join(L) + ''.join(R) + '->' + out, T); d = 3 ** (f - len(J))
        tot = tot + (-0.5) ** len(J) * q_dims_planar(Y.reshape(d * 4, d * 4), d)
    return 2 * tot
def q_dims_planar(Y, d):  # ||Y||^2 - ||Tr_s Y||^2 - ||Tr_t Y||^2 + ||Tr_st Y||^2 on [d, 2, 2], d-block never traced
    T = Y.reshape(d, 2, 2, d, 2, 2)
    a = torch.sum(torch.abs(T) ** 2)
    s = torch.sum(torch.abs(torch.einsum('xsiysj->xiyj', T)) ** 2)
    t = torch.sum(torch.abs(torch.einsum('xisyjs->xiyj', T)) ** 2)
    st = torch.sum(torch.abs(torch.einsum('xstyst->xy', T)) ** 2)
    return a - s - t + st
N = DF * 4
idmax = 0.0
for _ in range(5):
    X = torch.randn(N, N, dtype=cplx); X = X / torch.linalg.norm(X)
    idmax = max(idmax, abs(G(X) - G_closed(X)).item())
# smallest eigenvalue of the Hermitian form G on all of M_N (build the matrix of the linear map via polarization)
basis = torch.eye(N * N, dtype=cplx)
def Gpol(X, Y):  # sesquilinear form with G(X) = Gpol(X, X)
    return 0.25 * (G(X + Y) - G(X - Y)) + 0.25j * (G(X + 1j * Y) - G(X - 1j * Y))
if N * N <= 400:
    M = torch.zeros(N * N, N * N, dtype=cplx)
    for i in range(N * N):
        for j in range(i, N * N):
            v = Gpol(basis[i].reshape(N, N), basis[j].reshape(N, N)); M[i, j] = v; M[j, i] = v.conj()
    lam_min = torch.linalg.eigvalsh(M)[0].item(); lam_txt = f"{lam_min:+.6f}"
else:
    lam_txt = "(not built for this f; the closed form gives -2)"
# the rank-two example
e0 = torch.zeros(DF, dtype=cplx); e0[0] = 1; e1 = torch.zeros(DF, dtype=cplx); e1[-1] = 1
if f == 0: FX = torch.ones(1, 1, dtype=cplx)
else: FX = torch.outer(e0, e1.conj())
a = torch.tensor([1, 0], dtype=cplx); ap = torch.tensor([0, 1], dtype=cplx)
Xex = torch.kron(FX, torch.kron(torch.eye(2, dtype=cplx) / np.sqrt(2), torch.outer(a, ap.conj())))
ex = (G(Xex) / torch.linalg.norm(Xex) ** 2).item(); rk = torch.linalg.matrix_rank(Xex).item()
# adversarial minimum over rank-one and rank-two X
def adv(r):
    best = 1e9
    for _ in range(ST):
        p = torch.randn(4 * N * r, requires_grad=True)
        def mem(p):
            U = torch.complex(p[:N * r], p[N * r:2 * N * r]).reshape(N, r); V = torch.complex(p[2 * N * r:3 * N * r], p[3 * N * r:]).reshape(N, r)
            X = U @ V.conj().T; return X / torch.linalg.norm(X)
        opt = torch.optim.LBFGS([p], lr=1, max_iter=600, tolerance_grad=1e-13, tolerance_change=1e-17, history_size=30, line_search_fn='strong_wolfe')
        def closure():
            opt.zero_grad(); v = G(mem(p)); v.backward(); return v
        opt.step(closure)
        with torch.no_grad(): best = min(best, G(mem(p)).item())
    return best
print(f"f={f}: closed form error {idmax:.1e}; smallest eigenvalue of G/||X||^2 over all X = {lam_txt}; "
      f"example (rank {rk}) G/||X||^2 = {ex:+.6f}; adversarial min over rank one = {adv(1):+.2e}, over rank two = {adv(2):+.6f} ({ST} starts)")
