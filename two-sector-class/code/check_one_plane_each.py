"""Sector copies with at most one planar A_s and at most one planar B_s (types (2,1), (1,2), the rest (1,1)), f free copies.
Exact identity for EVERY C on H_F (x) (A + B):  2^m q_n(C) = 2 q_{F+K}(U C) + (2^m - 2)(qt_F(P) + qt_F(Q)),
U = P_A + (-1)^{m+1} P_B, K = A (+) B as one copy (Werner form).  Then adversarial min of q_n/|C|^2 over rank-2 C.
usage: python check_one_plane_each.py f types trials seed   (types e.g. '21,12,11')"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64); cplx = torch.complex128
f, TYPES, TR, SEED = int(sys.argv[1]), [tuple(int(c) for c in t) for t in sys.argv[2].split(',')], int(sys.argv[3]), int(sys.argv[4])
m = len(TYPES); n = f + m; D = 3 ** n; DF = 3 ** f
torch.manual_seed(SEED); rng = np.random.default_rng(SEED)
w9 = torch.eye(3, dtype=cplx).reshape(9) / np.sqrt(3)
def make_Phi(k):
    Dk = 3 ** k
    def Phi(C):
        if k == 0: return C
        X = C.reshape((3,) * (2 * k)).permute([i for s in range(k) for i in (s, k + s)]).reshape((9,) * k)
        for s in range(k):
            tr = torch.tensordot(w9.conj(), X, dims=([0], [s])); X = X - 1.5 * torch.moveaxis(torch.tensordot(w9, tr, dims=0), 0, s)
        return X.reshape([3, 3] * k).permute([2 * s for s in range(k)] + [2 * s + 1 for s in range(k)]).reshape(Dk, Dk)
    return Phi
Phi_n, Phi_F = make_Phi(n), make_Phi(f)
q = lambda C: torch.vdot(C.reshape(-1), Phi_n(C).reshape(-1)).real
qF = lambda Z: torch.vdot(Z.reshape(-1), Phi_F(Z).reshape(-1)).real
def kron_all(ms):
    out = torch.ones(1, 1, dtype=cplx)
    for M in ms: out = torch.kron(out, M)
    return out
Us = [torch.tensor(np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))[0]) for _ in range(m)]
EA = kron_all([torch.eye(DF, dtype=cplx)] + [Us[s][:, :TYPES[s][0]] for s in range(m)])
EB = kron_all([torch.eye(DF, dtype=cplx)] + [Us[s][:, TYPES[s][0]:sum(TYPES[s])] for s in range(m)])
dA, dB = EA.shape[1] // DF, EB.shape[1] // DF; dK = dA + dB
E = torch.cat([EA, EB], 1)                       # columns: (F index, A index) then (F index, B index)
def qFK(Ck):  # q_{F + K} for an operator given in coordinates (F (x) K), K = A (+) B
    C4 = Ck.reshape(DF, dK, DF, dK)
    tot = sum(qF(C4[:, i, :, j]) for i in range(dK) for j in range(dK))
    return tot - 0.5 * qF(torch.einsum('iaja->ij', C4))
def to_FK(Cs):  # Cs in the E-coordinates [(F,A) ; (F,B)] -> F (x) K ordering
    k = DF * dA
    XA = Cs.reshape(DF * dA + DF * dB, DF * dA + DF * dB)
    idx = []
    for i in range(DF):
        idx += [i * dA + a for a in range(dA)] + [k + i * dB + b for b in range(dB)]
    idx = torch.tensor(idx)
    return XA[idx][:, idx]
def qt(P, d1, d2):
    P4 = P.reshape(DF, d1, DF, d2); return sum(qF(P4[:, a, :, b]) for a in range(d1) for b in range(d2))
k = DF * dA
idmax = 0.0
for trial in range(3):
    Cs = torch.randn(DF * dK, DF * dK, dtype=cplx); C = E @ Cs @ E.conj().T
    Usgn = torch.ones(DF * dK, dtype=cplx); Usgn[k:] = (-1) ** (m + 1)
    lhs = 2 ** m * q(C)
    rhs = 2 * qFK(to_FK(torch.diag(Usgn) @ Cs)) + (2 ** m - 2) * (qt(Cs[:k, k:], dA, dB) + qt(Cs[k:, :k], dB, dA))
    idmax = max(idmax, (abs(lhs - rhs) / torch.linalg.norm(C) ** 2).item())
best = 1e9
for trial in range(TR):
    p = torch.randn(8 * DF * dK, requires_grad=True)
    def member(p):
        F_ = torch.complex(p[:2 * DF * dK], p[2 * DF * dK:4 * DF * dK]).reshape(DF * dK, 2)
        G_ = torch.complex(p[4 * DF * dK:6 * DF * dK], p[6 * DF * dK:]).reshape(DF * dK, 2)
        C = E @ F_ @ G_.conj().T @ E.conj().T; return C / torch.linalg.norm(C)
    opt = torch.optim.LBFGS([p], lr=1, max_iter=1500, tolerance_grad=1e-14, tolerance_change=1e-18, history_size=50, line_search_fn='strong_wolfe')
    def closure():
        opt.zero_grad(); v = q(member(p)); v.backward(); return v
    opt.step(closure)
    with torch.no_grad(): best = min(best, q(member(p)).item())
print(f"f={f} types={TYPES} (n={n}): identity max rel error {idmax:.1e}; adversarial min q_n/|C|^2 over rank-2 members = {best:+.3e} ({TR} starts)")
