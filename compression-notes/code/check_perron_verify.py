"""Section 4: finite-n check for a Perron-operator candidate (output A, B of check_perron_operator.py).
Builds the Perron operator X_*, minimises lambda_min((U (x) V)^* X_*^Gamma (U (x) V)) over isometries, takes the
eigenvector psi in C^2 (x) C^2, omega = (U (x) V) psi (Schmidt rank <= 2), splits omega = sum_a conj(l_a) (x) m_a,
picks boundary vectors r, r' (random, fixed seed) and evaluates the ratio q_n(1/2, C)/q_n(0, C) for the MPS operator
C = sum_a |y_a><z_a| (Proposition 4): transfer-channel evaluation in 50-digit arithmetic (mpmath) for n up to NMAX,
and a dense double-precision cross-check of q_n for n <= 6.   usage: python check_perron_verify.py file.npz NMAX [hp]   (needs PyTorch and mpmath; 'hp' = 50-digit transfer evaluation, small D only)"""
import sys, itertools, numpy as np, torch, mpmath as mp
torch.set_default_dtype(torch.float64)
mp.mp.dps = 50
FN, NMAX = sys.argv[1], int(sys.argv[2]); HP = len(sys.argv) > 3 and sys.argv[3] == 'hp'
P = np.load(FN); A, B = P["A"], P["B"]; D = A.shape[1]
def transfer_np(X, W):
    Y = np.zeros_like(X)
    for u in range(3):
        for v in range(3):
            Muv = np.kron(A[u], B[v])
            for x in range(3):
                for y in range(3):
                    w = W[x, y, u, v]
                    if w != 0: Y = Y + w * Muv @ X @ np.kron(A[x], B[y]).conj().T
    return Y
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
rho = (np.eye(9) - F / 2).reshape(3, 3, 3, 3); idt = np.eye(9).reshape(3, 3, 3, 3)
def pt(X): return X.reshape(D, D, D, D).transpose(0, 3, 2, 1).reshape(D * D, D * D)
X = np.eye(D * D, dtype=complex) / D ** 2
for k in range(2000):
    Xn = transfer_np(X, rho); Xn = (Xn + Xn.conj().T) / 2; Xn /= np.trace(Xn).real
    if np.linalg.norm(Xn - X) < 1e-15: break
    X = Xn
XG = torch.tensor(pt(X))
best = (np.inf, None)
for r in range(40):
    U = torch.randn(D, 2, dtype=torch.complex128, requires_grad=True); V = torch.randn(D, 2, dtype=torch.complex128, requires_grad=True)
    o = torch.optim.Adam([U, V], lr=0.05)
    for it in range(600):
        K = torch.kron(torch.linalg.qr(U)[0], torch.linalg.qr(V)[0]); e = torch.linalg.eigvalsh(K.conj().T @ XG @ K)[0]
        o.zero_grad(); e.backward(); o.step()
    if e.item() < best[0]: best = (e.item(), (torch.linalg.qr(U)[0].detach().numpy(), torch.linalg.qr(V)[0].detach().numpy()))
eta, (Uq, Vq) = best
K = np.kron(Uq, Vq); w, Z = np.linalg.eigh(K.conj().T @ pt(X) @ K); omega = K @ Z[:, 0]
print(f"D = {D}: eta(X_*) = {eta:+.6e} (Perron converged after {k} iterations)", flush=True)
Om = omega.reshape(D, D); Uo, s, Vh = np.linalg.svd(Om)
l = [np.conj(Uo[:, a] * np.sqrt(s[a])) for a in range(2)]; m = [Vh[a] * np.sqrt(s[a]) for a in range(2)]   # omega = sum conj(l_a) (x) m_a
rng = np.random.default_rng(0); r = rng.standard_normal(D) + 1j * rng.standard_normal(D); rp = rng.standard_normal(D) + 1j * rng.standard_normal(D)
# high-precision transfer evaluation
def mpmat(M): return mp.matrix([[mp.mpc(complex(M[i, j])) for j in range(M.shape[1])] for i in range(M.shape[0])])
Ms = {(u, v): mpmat(np.kron(A[u], B[v])) for u in range(3) for v in range(3)}
def T_mp(X, W):
    Y = mp.zeros(D * D, D * D)
    for u in range(3):
        for v in range(3):
            for x in range(3):
                for y in range(3):
                    w = W[x, y, u, v]
                    if w != 0: Y += mp.mpf(float(w)) * (Ms[(u, v)] * X * Ms[(x, y)].H)
    return Y
def pt_mp(X):
    Y = mp.zeros(D * D, D * D)
    for a in range(D):
        for b in range(D):
            for c in range(D):
                for d in range(D):
                    Y[a * D + b, c * D + d] = X[a * D + d, c * D + b]
    return Y
v0 = np.kron(r, rp)
if HP:
    X0 = mpmat(np.outer(v0, v0.conj())); Xr, Xi = X0, X0
    om = mp.matrix([mp.mpc(complex(z)) for z in omega])
    for nn in range(1, NMAX + 1):
        Xr = T_mp(Xr, rho); Xi = T_mp(Xi, idt)
        sc = mp.mpf(1) / mp.norm(Xi); Xr *= sc; Xi *= sc
        if nn in (1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64) or nn == NMAX:
            q = (om.H * pt_mp(Xr) * om)[0].real; q0 = (om.H * pt_mp(Xi) * om)[0].real
            print(f"   n = {nn} (50 digits): q_n(1/2, C)/q_n(0, C) = {mp.nstr(q / q0, 12)}", flush=True)
else:
    Xr = Xi = np.outer(v0, v0.conj())
    for nn in range(1, NMAX + 1):
        Xr = transfer_np(Xr, rho); Xi = transfer_np(Xi, idt)
        sc = 1 / np.linalg.norm(Xi); Xr = Xr * sc; Xi = Xi * sc
        if nn in (1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64) or nn == NMAX:
            q = (omega.conj() @ pt(Xr) @ omega).real; q0 = (omega.conj() @ pt(Xi) @ omega).real
            print(f"   n = {nn} (double): q_n(1/2, C)/q_n(0, C) = {q / q0:+.12e}", flush=True)
# dense cross-check for small n
def mps_vec(lv, T_, rv, n):
    v = lv.copy().reshape(1, D)
    for s_ in range(n): v = np.einsum('ka,xab->kxb', v.reshape(-1, D), T_).reshape(-1, D)
    return (v @ rv)
for nn in (2, 3, 4):
    Y = np.stack([mps_vec(l[a], A, r, nn) for a in range(2)], 1)   # y_a (3^n)
    Zm = np.stack([mps_vec(m[a], B, rp, nn) for a in range(2)], 1)
    C = Y @ Zm.conj().T
    Ct = C.reshape([3] * (2 * nn)); letters = 'abcdefghijklmnop'; tot = 0
    for rr in range(nn + 1):
        for J in itertools.combinations(range(nn), rr):
            inp = [letters[i] for i in range(2 * nn)]
            for j in J: inp[nn + j] = inp[j]
            out = [inp[i] for i in range(2 * nn) if not (i in J or (i >= nn and i - nn in J))]
            T_ = np.einsum(''.join(inp) + '->' + ''.join(out), Ct); tot += (-0.5) ** rr * np.sum(np.abs(T_) ** 2)
    print(f"   dense n = {nn}: q_n(1/2, C)/||C||^2 = {tot / np.sum(np.abs(C) ** 2):+.12f}", flush=True)
