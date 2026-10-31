"""Section 11: check of the dual certificate for the congruence question.  sigma = (7 Q_sss + 12 Q_sas + 12 Q_asa +
42 Q_aaa)/1998 (labels: qubit pair, copy 1, copy 2; s/a = symmetric/antisymmetric).  Prints <W_2, sigma> (negative)
and the minimum of Tr[W_1 (X (x) Y) sigma (X (x) Y)^*] over unit X, Y : C^18 -> C^6 found by 120 local searches
(complex and real starts, Adam then L-BFGS).  A negative minimum would invalidate the certificate.
usage: python check_congruence_certificate.py   (needs PyTorch; about twenty minutes)
"""
import sys, itertools, numpy as np, torch
from scipy.optimize import linprog
torch.set_default_dtype(torch.float64)
ROUNDS, RS, IT = 0, 0, 0
def swap(d):
    S = np.zeros((d * d, d * d))
    for i in range(d):
        for j in range(d): S[d * i + j, d * j + i] = 1
    return S
F2, F3 = swap(2), swap(3)
P2 = {'s': (np.eye(4) + F2) / 2, 'a': (np.eye(4) - F2) / 2}
P3 = {'s': (np.eye(9) + F3) / 2, 'a': (np.eye(9) - F3) / 2}
rho = np.eye(9) - F3 / 2
def arrange(ops, dims):
    """ops[k] acts on pair k (A_k B_k) with local dim dims[k]; return operator on (A_0 A_1 ..)(B_0 B_1 ..)."""
    m = len(ops); R = np.eye(1)
    for o in ops: R = np.kron(R, o)
    shape = []
    for d in dims: shape += [d, d]
    perm = [2 * i for i in range(m)] + [2 * i + 1 for i in range(m)]
    R = R.reshape(shape + shape).transpose(perm + [2 * m + q for q in perm])
    D = int(np.prod(dims)); return R.reshape(D * D, D * D)
labels = list(itertools.product('sa', repeat=3))
Q = [arrange([P2[e[0]], P3[e[1]], P3[e[2]]], [2, 3, 3]) for e in labels]
W1 = arrange([F2, rho], [2, 3]); W2 = arrange([F2, rho, rho], [2, 3, 3])
w = np.array([np.trace(W2 @ q) for q in Q]); trQ = np.array([np.trace(q) for q in Q])
print("labels", [''.join(e) for e in labels]); print("w =", np.round(w, 4)); print("Tr Q =", trQ)
Qt = torch.tensor(np.stack(Q), dtype=torch.complex128); W1t = torch.tensor(W1, dtype=torch.complex128)
def t_of(X, Y):
    K = torch.kron(X, Y)                                     # 36 x 324
    G = K.conj().T @ W1t @ K
    return torch.einsum('ij,eji->e', G, Qt).real
def oracle(c):
    sig = torch.tensor(np.einsum('e,eij->ij', c, np.stack(Q)), dtype=torch.complex128)
    best = (np.inf, None)
    for r in range(RS):
        torch.manual_seed(7919 * len(cuts) + r)
        X = torch.randn(6, 18, dtype=torch.complex128, requires_grad=True)
        Y = torch.randn(6, 18, dtype=torch.complex128, requires_grad=True)
        def fval():
            Xn = X / torch.linalg.norm(X); Yn = Y / torch.linalg.norm(Y)
            K = torch.kron(Xn, Yn)
            return torch.trace(W1t @ K @ sig @ K.conj().T).real
        opt = torch.optim.Adam([X, Y], lr=0.03)
        for it in range(IT):
            f = fval(); opt.zero_grad(); f.backward(); opt.step()
        opt2 = torch.optim.LBFGS([X, Y], lr=1, max_iter=400, line_search_fn='strong_wolfe', tolerance_grad=1e-14, tolerance_change=1e-16)
        def cl():
            opt2.zero_grad(); v = fval(); v.backward(); return v
        opt2.step(cl)
        f = fval()
        if f.item() < best[0]: best = (f.item(), ((X / torch.linalg.norm(X)).detach().clone(), (Y / torch.linalg.norm(Y)).detach().clone()))
    return best

c = np.zeros(8)
for lab, wt in (('sss', 7), ('sas', 12), ('asa', 12), ('aaa', 42)):
    c[[''.join(e) for e in labels].index(lab)] = wt / 1998
M = sum(ci * q for ci, q in zip(c, Q))
print("Tr sigma =", round(float(np.trace(M).real), 12), "| <W_2, sigma> =", float(np.trace(W2 @ M).real))
Mt = torch.tensor(M, dtype=torch.complex128)
def f(X, Y):
    K = torch.kron(X / torch.linalg.norm(X), Y / torch.linalg.norm(Y))
    return torch.trace(W1t @ K @ Mt @ K.conj().T).real
best = np.inf
for r in range(120):
    torch.manual_seed(55555 + r)
    if r % 2 == 0:
        X = torch.randn(6, 18, dtype=torch.complex128, requires_grad=True)
        Y = torch.randn(6, 18, dtype=torch.complex128, requires_grad=True)
    else:
        X = torch.randn(6, 18).to(torch.complex128).requires_grad_()
        Y = torch.randn(6, 18).to(torch.complex128).requires_grad_()
    opt = torch.optim.Adam([X, Y], lr=0.03)
    for it in range(1500):
        v = f(X, Y); opt.zero_grad(); v.backward(); opt.step()
    opt2 = torch.optim.LBFGS([X, Y], lr=1, max_iter=500, line_search_fn='strong_wolfe', tolerance_grad=1e-14,
                             tolerance_change=1e-16)
    def cl():
        opt2.zero_grad(); v = f(X, Y); v.backward(); return v
    opt2.step(cl)
    best = min(best, f(X, Y).item())
    if r % 20 == 19:
        print(f"after {r + 1} searches: minimum {best:+.3e}", flush=True)
