"""Section 9: checks of the identities used there, on random operators (n = 2, 3): the crossed form, the Hermitian
formula q_n(t, C) = Tr[(C (x) C) prod_s (F_s - t)], identity (9.3), the bookkeeping behind (9.2), and the odd-n
equality family lambda_1 |0^n><0^n| + lambda_2 |1^n><1^n| at the endpoint (n = 3, 5)."""
import itertools, numpy as np
rng = np.random.default_rng(7)
def ptr_norms(C, n, t):
    T = C.reshape([3]*n*2); tot = 0
    for r in range(n+1):
        for K in itertools.combinations(range(n), r):
            X = T; rem = list(range(n))
            for s in sorted(K, reverse=True):
                i = rem.index(s); X = np.trace(X, axis1=i, axis2=len(rem)+i); rem.pop(i)
            tot += (-t)**len(K) * np.sum(abs(X)**2)
    return tot
def perm_op(n, which):
    """operator on (C^3)^n (x) (C^3)^n swapping copy s of the two factors for s in which"""
    d = 3**n; M = np.zeros((d*d, d*d))
    for idx in itertools.product(range(3), repeat=2*n):
        a, b = list(idx[:n]), list(idx[n:])
        for s in which: a[s], b[s] = b[s], a[s]
        src = np.ravel_multi_index(idx, [3]*2*n); dst = np.ravel_multi_index(tuple(a+b), [3]*2*n)
        M[dst, src] = 1
    return M
def crossed(A, B, n, t):
    d = 3**n; tot = 0
    # N = (I - tF)^{(x) n}, F_s swaps copy s of the A-factor and the B-factor
    N = np.eye(d*d)
    Fs = [perm_op(n, [s]) for s in range(n)]
    N = np.eye(d*d)
    for s in range(n): N = N @ (np.eye(d*d) - t*Fs[s])
    for i in range(2):
        for j in range(2):
            tot += np.vdot(np.kron(A[:, i], B[:, j]), N @ np.kron(A[:, j], B[:, i]))
    return tot, Fs
for n in (2, 3):
    d = 3**n; t = 0.37
    A = rng.standard_normal((d, 2)) + 1j*rng.standard_normal((d, 2)); B = rng.standard_normal((d, 2)) + 1j*rng.standard_normal((d, 2))
    C = A @ B.conj().T
    cr, Fs = crossed(A, B, n, t)
    print(f"n={n}: crossed form {cr.real:.12f} {cr.imag:+.1e}i   direct {ptr_norms(C, n, t):.12f}")
    # Hermitian identity q = Tr[(C x C) prod (F_s - t)]
    H = C + C.conj().T
    W = np.eye(d*d)
    for s in range(n): W = W @ (Fs[s] - t*np.eye(d*d))
    print(f"   Hermitian: Tr[(C x C) prod(F_s - t)] = {np.trace(np.kron(H, H) @ W).real:.12f}   direct {ptr_norms(H, n, t):.12f}")
    # (9.3) at t = 1/2: 2^n q - (2||H||^2 - (TrH)^2) = sum_T c_|T| Tr[(H x H) Pi_T]
    PS = [(np.eye(d*d) + F)/2 for F in Fs]; PA = [(np.eye(d*d) - F)/2 for F in Fs]
    rhs = 0
    for r in range(n+1):
        for T in itertools.combinations(range(n), r):
            Pi = np.eye(d*d)
            for s in range(n): Pi = Pi @ (PA[s] if s in T else PS[s])
            c = (-3)**r - 2*(-1)**r + 1
            rhs += c * np.trace(np.kron(H, H) @ Pi).real
    lhs = 2**n*ptr_norms(H, n, 0.5) - (2*np.linalg.norm(H)**2 - np.trace(H).real**2)
    print(f"   (9.3): lhs {lhs:.10f}  rhs {rhs:.10f}")
    # (9.2) frame form: q - RHS = (1-t)^n tau [s1^2 G11 + s2^2 G22 + 2 s1 s2 (Re Gx + 1)]
    U, S, Vh = np.linalg.svd(C); u = U[:, :2]; v = Vh[:2].conj().T; s1, s2 = S[:2]
    tau = t/(1-t)
    G = np.eye(d*d)
    for s in range(n): G = G @ (np.eye(d*d) + 2*tau*PA[s])
    G = (G - np.eye(d*d))/tau
    G11 = np.vdot(np.kron(u[:, 0], v[:, 0]), G @ np.kron(u[:, 0], v[:, 0])).real
    G22 = np.vdot(np.kron(u[:, 1], v[:, 1]), G @ np.kron(u[:, 1], v[:, 1])).real
    Gx = np.vdot(np.kron(u[:, 0], v[:, 1]), G @ np.kron(u[:, 1], v[:, 0]))
    lhs2 = ptr_norms(C, n, t) - (1-t)**(n-1)*(np.linalg.norm(C)**2 - t*(s1+s2)**2)
    rhs2 = (1-t)**n*tau*(s1**2*G11 + s2**2*G22 + 2*s1*s2*(Gx.real + 1))
    print(f"   (9.2) bookkeeping: q - RHS = {lhs2:.10f}   frame expression {rhs2:.10f}")
# odd-n GHZ equality at t = 1/2
for n in (3, 5):
    d = 3**n; C = np.zeros((d, d)); l1, l2 = 0.8, 0.3
    C[0, 0] = l1; i1 = np.ravel_multi_index((1,)*n, (3,)*n); C[i1, i1] = l2
    print(f"GHZ n={n}: 2^n q_n = {2**n*ptr_norms(C, n, 0.5):.12f}   (l1-l2)^2 = {(l1-l2)**2:.12f}")
