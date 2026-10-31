"""Structured counterexample search for EP_b inside NPT code blocks ([[b, b-2, 2]] codes, b >= 7).

For C with range in the code space Q (3^k dims, k = b - 2) and co-range in one syndrome sector of the other side,
q_b(1/2, C) = sum_E p(E) |<Phi_E|psi_C>|^2 (Proposition 14(a)). In logical coordinates (Phi_{E0 + l} = (L(l) (x) I)
Phi_{E0}), this is   q(M) = sum_{u, v in F_3^k} p(E0 + u.X~ + v.Z~) |Tr((X^u Z^v)^dagger M)|^2,   M = rank-two
operator on C^(3^k). We minimise q(M) / sum |Tr(.)|^2 (a weighted average of p) over rank-two M.

The weight of a Pauli X^a Z^z (a, z in F_3^b) is (1/9) sum_{cells} (-1/2)^{N_cell}, where site i sits in cell
(a_i, z_i / c_i)."""
import itertools, sys, numpy as np, torch
torch.set_default_dtype(torch.float64)


def logical_basis(c):
    """X~_j = X^{a_j} (a_j . c = 0), Z~_j = Z^{e_j} (sum e_j = 0), with a_i . e_j = delta_ij (mod 3)."""
    b = len(c); c = np.array(c) % 3; one = np.ones(b, int)
    # candidates
    A = [np.array(v) for v in itertools.product(range(3), repeat=b) if (np.array(v) @ c) % 3 == 0]
    E = [np.array(v) for v in itertools.product(range(3), repeat=b) if np.array(v).sum() % 3 == 0]
    def rank(vs):
        M = [list(v) for v in vs]; r = 0
        for col in range(b):
            piv = next((i for i in range(r, len(M)) if M[i][col] % 3), None)
            if piv is None: continue
            M[r], M[piv] = M[piv], M[r]; inv = 1 if M[r][col] % 3 == 1 else 2
            M[r] = [(x * inv) % 3 for x in M[r]]
            for i in range(len(M)):
                if i != r and M[i][col] % 3:
                    f = M[i][col]; M[i] = [(x - f * y) % 3 for x, y in zip(M[i], M[r])]
            r += 1
        return r
    xs, zs = [], []
    k = b - 2
    for _ in range(k):
        # pick a with a.c = 0, independent of {1} + xs, and a . z = 0 for previous zs
        a = next(v for v in A if rank([one] + xs + [v]) > 1 + len(xs) and all((v @ z) % 3 == 0 for z in zs))
        # pick e with sum e = 0, independent of {c} + zs, e . x = 0 for previous xs, e . a = 1
        e = next(v for v in E if rank([c] + zs + [v]) > 1 + len(zs) and all((v @ x) % 3 == 0 for x in xs)
                 and (v @ a) % 3 == 1)
        xs.append(a); zs.append(e)
    return np.array(xs), np.array(zs)


def weights(c, sector):
    """p over (u, v) in F_3^k x F_3^k for the given sector (alpha1, alpha2): E0 = Z^{f} X^{g}."""
    b = len(c); c = np.array(c) % 3; k = b - 2
    xs, zs = logical_basis(c)
    cinv = np.array([1 if ci == 1 else 2 for ci in c])
    # E0: X-part g with g.c = alpha2, Z-part f with sum f = alpha1 (any representative)
    g = np.zeros(b, int); g[0] = (sector[1] * cinv[0]) % 3
    f = np.zeros(b, int); f[0] = sector[0] % 3
    U = np.array(list(itertools.product(range(3), repeat=k)))           # 3^k x k
    ax = (g[None, :] + U @ xs) % 3                                        # X-exponents for each u
    zz = (f[None, :] + U @ zs) % 3                                        # Z-exponents for each v
    # cell index per site: (a_i, z_i * c_i^{-1}) -> 3 a + z'
    cellx = ax                                                            # [u, i]
    cellz = (zz * cinv[None, :]) % 3                                      # [v, i]
    cell = 3 * cellx[:, None, :] + cellz[None, :, :]                      # [u, v, i]
    counts = np.stack([(cell == m).sum(-1) for m in range(9)], -1)
    p = ((-0.5) ** counts).sum(-1) / 9.0
    return p                                                              # shape (3^k, 3^k), [u, v]


def pauli_spectrum(Ur, Ui, Wr, Wi, k):
    """|Tr((X^u Z^v)^dag M)|^2 for M = U W^T (complex), via shifted diagonals and an FFT."""
    d = 3 ** k
    Uc = torch.complex(Ur, Ui); Wc = torch.complex(Wr, Wi)
    J = torch.tensor(np.array(list(itertools.product(range(3), repeat=k))))  # d x k
    pw = torch.tensor(3 ** np.arange(k)[::-1].copy())
    shift = ((J[:, None, :] + J[None, :, :]) % 3) @ pw                      # [u, j] -> index of j + u
    D = (Uc[shift] * Wc[None, :, :]).sum(-1)                                # M[j+u, j] = sum_a U[j+u,a] W[j,a]
    D = D.reshape((d,) + (3,) * k)
    C = torch.fft.fftn(D, dim=tuple(range(1, k + 1))).reshape(d, d)        # [u, v]: sum_j w^{-v.j} M[j+u, j]
    return C.real ** 2 + C.imag ** 2


def search(c, restarts=8, iters=1500, seed=0, sectors=None):
    b = len(c); k = b - 2; d = 3 ** k
    best_all = None
    for sector in (sectors or list(itertools.product(range(3), repeat=2))):
        p = torch.tensor(weights(c, sector))
        gen = torch.Generator().manual_seed(seed)
        best = 1e9
        for r in range(restarts):
            params = [torch.randn(d, 2, generator=gen, requires_grad=True) for _ in range(4)]
            opt = torch.optim.Adam(params, lr=0.05)
            for it in range(iters):
                S = pauli_spectrum(*params, k)
                val = (p * S).sum() / S.sum()
                opt.zero_grad(); val.backward(); opt.step()
            best = min(best, val.item())
        print(f"  [[{b},{k},2]] c = {c}, sector {sector}: p in [{p.min().item():+.4f}, {p.max().item():+.4f}], "
              f"negative fraction {(p < 0).double().mean().item():.4f}; best normalised q = {best:+.6f}", flush=True)
        best_all = best if best_all is None else min(best_all, best)
    return best_all


if __name__ == "__main__":
    b = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    cs = {7: [1, 1, 1, 1, 1, 2, 2], 8: [1, 1, 1, 1, 1, 1, 1, 2], 9: [1] * 9}
    R = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    I = int(sys.argv[3]) if len(sys.argv) > 3 else 1500
    secs = None if len(sys.argv) <= 4 else [tuple(int(x) for x in t) for t in sys.argv[4].split(",")]
    print(f"b = {b}: minimum over sectors {search(cs[b], restarts=R, iters=I, sectors=secs):+.6f}")
