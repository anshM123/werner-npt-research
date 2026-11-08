"""Section 9: the one-step inequality that slice-by-slice induction would need on local qubits (it fails).
M = A B^* rank two on C^D (x) C^2 (last factor = peeled qubit).  M_w = (I (x) w^T) M (I (x) w).
OSI-Q:  ||M||_1^2 >= ||M_01 - M_10||_F^2 + 12 E_w s1 s2 (M_w),  w uniform on the unit sphere of C^2.
E_w by Gauss-Legendre (cos theta) x uniform (phi) quadrature on the Bloch sphere (integrand is phase invariant).
usage: python one_step_qubit.py D runs seed nth nph"""
import sys
import numpy as np
import torch
torch.set_default_dtype(torch.float64)
D, runs, seed, nth, nph = map(int, sys.argv[1:6])

def grid(nth, nph):
    x, wx = np.polynomial.legendre.leggauss(nth)        # x = cos theta
    ph = np.arange(nph) * 2 * np.pi / nph
    th = np.arccos(x)
    W = []; Wt = []
    for i in range(nth):
        for j in range(nph):
            W.append([np.cos(th[i] / 2), np.exp(1j * ph[j]) * np.sin(th[i] / 2)])
            Wt.append(wx[i] / 2 / nph)
    return torch.tensor(np.array(W)), torch.tensor(np.array(Wt))

Wg, Wt = grid(nth, nph)

def parts(A, B):
    # A, B: (2D) x 2 complex, index = d*2 + a  (last factor is the qubit)
    A3 = A.reshape(D, 2, 2); B3 = B.reshape(D, 2, 2)
    Aw = torch.einsum('dak,ga->gdk', A3, Wg)              # (I x w^T) A
    Bw = torch.einsum('dak,ga->gdk', B3, Wg.conj())       # (I x w^dagger) B
    GA = Aw.conj().transpose(1, 2) @ Aw; GB = Bw.conj().transpose(1, 2) @ Bw
    dA = (GA[:, 0, 0] * GA[:, 1, 1] - GA[:, 0, 1] * GA[:, 1, 0]).real
    dB = (GB[:, 0, 0] * GB[:, 1, 1] - GB[:, 0, 1] * GB[:, 1, 0]).real
    Es = (Wt * torch.sqrt(torch.clamp(dA * dB, min=0) + 1e-300)).sum()
    M = A @ B.conj().T
    M4 = M.reshape(D, 2, D, 2)
    J = M4[:, 0, :, 1] - M4[:, 1, :, 0]
    nJ = (J.real ** 2 + J.imag ** 2).sum()
    n2 = (M.real ** 2 + M.imag ** 2).sum()
    G1 = A.conj().T @ A; G2 = B.conj().T @ B
    s1s2 = torch.sqrt(torch.clamp(((G1[0, 0] * G1[1, 1] - G1[0, 1] * G1[1, 0]) * (G2[0, 0] * G2[1, 1] - G2[0, 1] * G2[1, 0])).real, min=0) + 1e-300)
    tn2 = n2 + 2 * s1s2
    return tn2, nJ, Es, n2

def gap(params):
    A = torch.complex(params[0], params[1]); B = torch.complex(params[2], params[3])
    tn2, nJ, Es, n2 = parts(A, B)
    return (tn2 - nJ - 12 * Es) / n2

rng = np.random.default_rng(seed)
res = []
for r in range(runs):
    params = [torch.tensor(rng.standard_normal((2 * D, 2)), requires_grad=True) for _ in range(4)]
    g0 = gap(params).item()
    opt = torch.optim.LBFGS(params, lr=1, max_iter=2000, tolerance_grad=1e-12, tolerance_change=1e-15,
                            history_size=50, line_search_fn='strong_wolfe')
    def closure():
        opt.zero_grad(); v = gap(params); v.backward(); return v
    for _ in range(3):
        opt.step(closure)
    g = gap(params).item()
    res.append(g)
    print(f"run {r}: random start gap {g0:.4f} -> minimised gap {g:.3e}", flush=True)
print(f"D={D}: smallest OSI-Q gap {min(res):.3e} (quadrature {nth}x{nph})")
