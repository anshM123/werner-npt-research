"""Check the three operators in data/ (t = 1/2): two at three copies, one at four copies.

data/single_copy_min_n3.npy: thick, critical for squeezes at copy 0, and the copy-0 squeeze form F_0 is positive on all
    traceless directions, so C is a strict minimum of R over its whole copy-0 squeeze orbit (one copy is not enough).
data/copywise_min_n3.npy, data/copywise_min_n4.npy: thick, critical at every copy, and every single-copy form F_s is
    positive semidefinite, yet the joint squeeze Hessian has negative directions, already for each pair of copies.

For each file: R, thickness at every copy, the squeeze gradient at every copy, the eigenvalues of the exact 9 x 9 forms
F_s (block formula, F_s(I) = 0 gives one zero), the smallest eigenvalue of the full squeeze Hessian and of each pair block.
usage: python one_copy_counterexamples.py   (run from orbit-notes/)"""
import sys, itertools
import numpy as np
import torch
sys.path.insert(0, 'code')
HB = []
for a in range(3):
    E = np.zeros((3, 3), complex); E[a, a] = 1; HB.append(E)
for a, b in itertools.combinations(range(3), 2):
    E = np.zeros((3, 3), complex); E[a, b] = E[b, a] = 1 / np.sqrt(2); HB.append(E)
    E = np.zeros((3, 3), complex); E[a, b] = -1j / np.sqrt(2); E[b, a] = 1j / np.sqrt(2); HB.append(E)


def F_matrix(C, s, r):
    F = np.zeros((9, 9))
    f = lambda P: SH.F_form(C, s, r, P)
    for i in range(9):
        for j in range(i, 9):
            F[i, j] = F[j, i] = 0.5 * (f(HB[i] + HB[j]) - f(HB[i]) - f(HB[j]))
    return F


def thickness_each(C):
    U, S, Vh = np.linalg.svd(C); out = []
    for M in (U[:, :2], Vh[:2].conj().T):
        P = (M @ M.conj().T).reshape((3,) * (2 * n))
        for s in range(n):
            T = np.moveaxis(P, [s, n + s], [0, 1]).reshape(3, 3, 3 ** (n - 1), 3 ** (n - 1))
            out.append(np.linalg.eigvalsh(np.einsum('abxx->ab', T))[0])
    return out


for fn, n in (('data/single_copy_min_n3.npy', 3), ('data/copywise_min_n3.npy', 3), ('data/copywise_min_n4.npy', 4)):
    sys.argv = [sys.argv[0], str(n), '0.5', '0.1', '0', '0']
    if 'squeeze_hessian' in sys.modules:
        del sys.modules['squeeze_hessian']
    import squeeze_hessian as SH
    C = np.load(fn); C = C / np.linalg.norm(C); r = SH.q_np(C, n)
    Ct = torch.tensor(C); h = torch.zeros(8 * n, requires_grad=True)
    g = torch.autograd.grad(SH.squeeze_R(Ct, h), h)[0].numpy()
    H = torch.autograd.functional.hessian(lambda hh: SH.squeeze_R(Ct, hh), torch.zeros(8 * n)); H = ((H + H.T) / 2).numpy()
    print(f"{fn}: R = {r:+.6f}, rank = {np.sum(np.linalg.svd(C, compute_uv=False) > 1e-10)}")
    print(f"   thickness (range at each copy, then co-range at each copy): {np.round(thickness_each(C), 4)}")
    for s in range(n):
        w = np.linalg.eigvalsh(F_matrix(C, s, r))
        print(f"   copy {s}: |squeeze gradient| {np.linalg.norm(g[8 * s:8 * s + 8]):.1e}; F_{s} eigenvalues "
              f"{np.array2string(w, precision=6, max_line_width=200)}")
    pairs = [np.linalg.eigvalsh(H[np.ix_(list(range(8 * a, 8 * a + 8)) + list(range(8 * b, 8 * b + 8)),
                                          list(range(8 * a, 8 * a + 8)) + list(range(8 * b, 8 * b + 8)))])[0]
             for a, b in itertools.combinations(range(n), 2)]
    w = np.linalg.eigvalsh(H)
    print(f"   full squeeze Hessian: smallest eigenvalue {w[0]:+.4f}, {np.sum(w < -1e-6)} negative of {8 * n}; "
          f"smallest eigenvalue of each pair block {np.round(pairs, 4)}")
