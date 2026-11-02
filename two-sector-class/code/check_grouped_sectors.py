"""Grouped reduction: sector copies split into groups, each with at most one planar A_s and one planar B_s.
Identity: q_n(C) = 2^{g-m} q_{F+K_1..K_g}(U C) + (1 - 2^{g-m})(qt_F(P) + qt_F(Q)), U = P_A + (-1)^{m-g} P_B,
K_i = A_{G_i} (+) B_{G_i}.  Checks the identity for random C (any rank) and adversarial min of q_n over rank-2 members.
usage: python check_grouped_sectors.py f groups trials seed   e.g. groups '21+12,21+12'  (types within a group joined by +)"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64); cplx = torch.complex128
f = int(sys.argv[1]); GROUPS = [[tuple(int(c) for c in t) for t in g.split('+')] for g in sys.argv[2].split(',')]
TR, SEED = int(sys.argv[3]), int(sys.argv[4]); torch.manual_seed(SEED); rng = np.random.default_rng(SEED)
TYPES = [t for g in GROUPS for t in g]; m = len(TYPES); g = len(GROUPS); n = f + m; D = 3 ** n; DF = 3 ** f
w9 = torch.eye(3, dtype=cplx).reshape(9) / np.sqrt(3)
def q_dims(C, dims):  # endpoint form on a tensor product with arbitrary factor dims
    k = len(dims); T = C.reshape(dims + dims); tot = 0.0
    for mask in range(2 ** k):
        J = [i for i in range(k) if mask >> i & 1]; keep = [i for i in range(k) if i not in J]
        letters = 'abcdefghijklmnop'; L = list(letters[:k]); R_ = [letters[k + i] if i in keep else letters[i] for i in range(k)]
        out = ''.join(L[i] for i in keep) + ''.join(R_[i] for i in keep)
        X = torch.einsum(''.join(L) + ''.join(R_) + '->' + out, T)
        tot = tot + (-0.5) ** len(J) * torch.sum(torch.abs(X) ** 2)
    return tot
def kron_all(ms):
    out = torch.ones(1, 1, dtype=cplx)
    for M in ms: out = torch.kron(out, M)
    return out
Us = [torch.tensor(np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))[0]) for _ in range(m)]
Abas = [Us[s][:, :TYPES[s][0]] for s in range(m)]; Bbas = [Us[s][:, TYPES[s][0]:sum(TYPES[s])] for s in range(m)]
EA = kron_all([torch.eye(DF, dtype=cplx)] + Abas); EB = kron_all([torch.eye(DF, dtype=cplx)] + Bbas)
dA = EA.shape[1] // DF; dB = EB.shape[1] // DF
# group dims
sizes, gi = [], 0
for grp in GROUPS:
    dAg = int(np.prod([t[0] for t in grp])); dBg = int(np.prod([t[1] for t in grp])); sizes.append((dAg, dBg))
dK = [a + b for a, b in sizes]
def to_K(Cs):  # Cs in coordinates [(F,A);(F,B)] -> operator on H_F (x) K_1 (x) ... (x) K_g
    N = DF * int(np.prod(dK)); Mtx = torch.zeros(N, DF * dA + DF * dB, dtype=cplx)
    Adims = [t[0] for t in TYPES]; Bdims = [t[1] for t in TYPES]
    for i in range(DF):
        for ai, aidx in enumerate(itertools.product(*[range(d) for d in Adims])):
            # split the A multi-index by groups
            pos, kcoords = 0, []
            for gi_, grp in enumerate(GROUPS):
                sub = aidx[pos:pos + len(grp)]; pos += len(grp)
                kcoords.append(int(np.ravel_multi_index(sub, [t[0] for t in grp])) if len(sub) else 0)
            row = np.ravel_multi_index([i] + kcoords, [DF] + dK); Mtx[row, i * dA + ai] = 1
        for bi, bidx in enumerate(itertools.product(*[range(d) for d in Bdims])):
            pos, kcoords = 0, []
            for gi_, grp in enumerate(GROUPS):
                sub = bidx[pos:pos + len(grp)]; pos += len(grp)
                kcoords.append(sizes[gi_][0] + int(np.ravel_multi_index(sub, [t[1] for t in grp])))
            row = np.ravel_multi_index([i] + kcoords, [DF] + dK); Mtx[row, DF * dA + i * dB + bi] = 1
    return Mtx @ Cs @ Mtx.conj().T
def qtF(P, d1, d2):
    P4 = P.reshape(DF, d1, DF, d2); return sum(q_dims(P4[:, a, :, b], [3] * f) for a in range(d1) for b in range(d2))
E = torch.cat([EA, EB], 1); k = DF * dA
idmax = 0.0
for trial in range(3):
    Cs = torch.randn(DF * (dA + dB), DF * (dA + dB), dtype=cplx); C = E @ Cs @ E.conj().T
    U = torch.ones(DF * (dA + dB), dtype=cplx); U[k:] = (-1) ** (m - g)
    lhs = q_dims(C, [3] * n)
    rhs = 2.0 ** (g - m) * q_dims(to_K(torch.diag(U) @ Cs), [3] * f + dK) + (1 - 2.0 ** (g - m)) * (qtF(Cs[:k, k:], dA, dB) + qtF(Cs[k:, :k], dB, dA))
    idmax = max(idmax, (abs(lhs - rhs) / torch.linalg.norm(C) ** 2).item())
best = 1e9
dV = DF * (dA + dB)
for trial in range(TR):
    p = torch.randn(8 * dV, requires_grad=True)
    def member(p):
        F_ = torch.complex(p[:2 * dV], p[2 * dV:4 * dV]).reshape(dV, 2); G_ = torch.complex(p[4 * dV:6 * dV], p[6 * dV:]).reshape(dV, 2)
        C = E @ F_ @ G_.conj().T @ E.conj().T; return C / torch.linalg.norm(C)
    opt = torch.optim.LBFGS([p], lr=1, max_iter=800, tolerance_grad=1e-13, tolerance_change=1e-17, history_size=40, line_search_fn='strong_wolfe')
    def closure():
        opt.zero_grad(); v = q_dims(member(p), [3] * n); v.backward(); return v
    opt.step(closure)
    with torch.no_grad(): best = min(best, q_dims(member(p), [3] * n).item())
print(f"f={f} groups={GROUPS} (n={n}, g={g}): identity max rel error {idmax:.1e}; adversarial min q_n/|C|^2 over rank-2 = {best:+.3e} ({TR} starts)")
