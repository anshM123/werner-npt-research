"""Section 9, L_1(6): are 6 x 6 local filterings of rho (x) rho in K_1(6) = cone{products, (a(x)b) rho (a(x)b)^*}?
Fits by J filterings of rho plus P products (relative Frobenius residuals, L-BFGS from random starts):
  (a) controls: an element of K_1(6) (3 filterings + 10 products) and a Schmidt-rank-4 pure state;
  (b) a random filtering and an NPT filtering (partial-transpose minimum pushed to -0.01), at two term budgets;
  (c) the pair compressions tau_T of random rank-two operators on three copies (removal by compression at n = 3).
Needs PyTorch; about half an hour."""
import numpy as np, torch
torch.set_default_dtype(torch.float64)
D = 6
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2
rho2 = np.kron(rho, rho).reshape([3] * 8).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(81, 81)   # (A1 A2)(B1 B2)
R1 = torch.tensor(rho, dtype=torch.complex128)
R2 = torch.tensor(rho2, dtype=torch.complex128)


def pt(Y):
    return Y.reshape(D, D, D, D).permute(0, 3, 2, 1).reshape(D * D, D * D)


def fit(Y, J, P, IT=4000, RS=3):
    Yt = torch.tensor(Y, dtype=torch.complex128)
    best = np.inf
    for _ in range(RS):
        p = [(0.4 * torch.randn(*s, dtype=torch.complex128)).requires_grad_()
             for s in ((J, D, 3), (J, D, 3), (P, D), (P, D))]

        def model():
            V = torch.einsum('ka,kb->kab', p[2], p[3]).reshape(P, D * D)
            Kt = torch.einsum('jab,jcd->jacbd', p[0], p[1]).reshape(J, D * D, 9)
            return V.T @ V.conj() + (Kt @ R1 @ Kt.conj().transpose(1, 2)).sum(0)
        opt = torch.optim.LBFGS(p, lr=1, max_iter=IT, history_size=50, line_search_fn='strong_wolfe',
                                tolerance_grad=1e-15, tolerance_change=1e-18)

        def closure():
            opt.zero_grad()
            loss = torch.linalg.norm(model() - Yt) ** 2
            loss.backward()
            return loss
        opt.step(closure)
        with torch.no_grad():
            best = min(best, float(torch.linalg.norm(model() - Yt) / torch.linalg.norm(Yt)))
    return best


torch.manual_seed(0)
g = np.random.default_rng(7)
a = g.standard_normal((3, D, 3)) + 1j * g.standard_normal((3, D, 3))
b = g.standard_normal((3, D, 3)) + 1j * g.standard_normal((3, D, 3))
Yc = sum(np.kron(a[i], b[i]) @ rho @ np.kron(a[i], b[i]).conj().T for i in range(3))
xs = g.standard_normal((10, D)) + 1j * g.standard_normal((10, D))
ys = g.standard_normal((10, D)) + 1j * g.standard_normal((10, D))
Yc = Yc + sum(np.outer(np.kron(x, y), np.kron(x, y).conj()) for x, y in zip(xs, ys))
Yc /= np.trace(Yc).real
v = np.zeros((D, D), complex)
v[0, 1] = v[2, 3] = 0.5
v[1, 0] = v[3, 2] = -0.5
Yo = np.outer(v.reshape(-1), v.reshape(-1).conj())
print("(a) controls, 16 filterings + 80 products:")
print(f"    element of K_1(6): {fit(Yc, 16, 80):.2e};  Schmidt-rank-4 state: {fit(Yo, 16, 80, RS=2):.2e}", flush=True)

torch.manual_seed(300)
A = torch.randn(D, 9, dtype=torch.complex128)
B = torch.randn(D, 9, dtype=torch.complex128)
K = torch.kron(A, B)
Y0 = K @ R2 @ K.conj().T
Y0 = (Y0 + Y0.conj().T) / 2
Y0 = (Y0 / torch.trace(Y0).real).numpy()
A = A.clone().requires_grad_()
B = B.clone().requires_grad_()
opt = torch.optim.Adam([A, B], lr=0.005)
for it in range(4000):
    K = torch.kron(A, B)
    Yt = K @ R2 @ K.conj().T
    Yt = (Yt + Yt.conj().T) / 2
    Yt = Yt / torch.trace(Yt).real
    f = torch.linalg.eigvalsh(pt(Yt))[0]
    if f.item() < -0.01:
        break
    opt.zero_grad()
    f.backward()
    opt.step()
Y1 = Yt.detach().numpy()
print("(b) filterings of rho (x) rho:")
for (J, P) in [(16, 80), (32, 160)]:
    print(f"    {J} + {P}: element of K_1(6) {fit(Yc, J, P):.2e} | random filtering {fit(Y0, J, P):.2e} | "
          f"NPT filtering (PT minimum {f.item():+.4f}) {fit(Y1, J, P):.2e}", flush=True)

print("(c) pair compressions at n = 3 (range and co-range of C at copies 1, 2):")
rng = np.random.default_rng(11)
for c in range(3):
    y = rng.standard_normal((27, 2)) + 1j * rng.standard_normal((27, 2))
    z = rng.standard_normal((27, 2)) + 1j * rng.standard_normal((27, 2))
    Ys = y.reshape(9, 3, 2).transpose(0, 2, 1).reshape(9, 6)
    Zs = z.reshape(9, 3, 2).transpose(0, 2, 1).reshape(9, 6)
    US = np.linalg.qr(Ys)[0]
    UZ = np.linalg.qr(Zs)[0]
    Kmat = np.kron(US.conj().T, UZ.conj().T)
    tau = Kmat @ rho2 @ Kmat.conj().T
    tau = (tau + tau.conj().T) / 2
    tau /= np.trace(tau).real
    ptmin = np.linalg.eigvalsh(pt(torch.tensor(tau)).numpy())[0]
    print(f"    C {c}: PT minimum {ptmin:+.4f}; fit by 16 + 80: {fit(tau, 16, 80):.2e}", flush=True)
