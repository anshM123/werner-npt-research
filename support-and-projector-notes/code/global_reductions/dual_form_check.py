"""(a) Dual enumerator form: q_n(P) = 2 (5/6)^n sum_j (-1/5)^j B_j for rank-two projectors P (random, n = 2, 3).
(b) Stabilizer planes: for a stabilizer state with group M (Lagrangian in F_3^{2N}) and a Pauli x,
    q(P_c + P_xc) >= 0 and q(P_c - P_xc) >= 0  <=>  |W_{x+M}| <= W_M,  W_A = sum_{y in A} (-1/5)^{wt y}.
    Checked on every coset of random Lagrangians, N <= 8, through u_hat(chi) = sum_{h in M} (-1/2)^{Z(h)} chi(h)."""
import itertools, numpy as np
rng = np.random.default_rng(124)
w3 = np.exp(2j*np.pi/3); X = np.roll(np.eye(3), 1, axis=0); Zm = np.diag([1, w3, w3**2])
def weyl(a, b): return np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Zm, b)
def q(C, n):
    T = C.reshape([3]*(2*n)).astype(complex)
    for s in range(n):
        Tp = np.moveaxis(T, [s, n+s], [2*n-2, 2*n-1]); tr = np.trace(Tp, axis1=-2, axis2=-1)
        T = np.moveaxis(Tp - 0.5*tr[..., None, None]*np.eye(3), [2*n-2, 2*n-1], [s, n+s])
    return np.vdot(C, T.reshape(3**n, 3**n)).real
err = 0
for n in (2, 3):
    for _ in range(3):
        Q = np.linalg.qr(rng.normal(size=(3**n, 2)) + 1j*rng.normal(size=(3**n, 2)))[0]; P = Q @ Q.conj().T
        B = np.zeros(n+1)
        for e in itertools.product(range(9), repeat=n):
            E = np.ones((1, 1))
            for s in e: E = np.kron(E, weyl(s // 3, s % 3))
            B[sum(1 for s in e if s)] += 0.5*np.trace(E @ P @ E.conj().T @ P).real
        err = max(err, abs(q(P, n) - 2*(5/6)**n*sum((-0.2)**j*B[j] for j in range(n+1))))
print(f"(a) dual form, random rank-two projectors n = 2, 3: max error {err:.1e}")
def random_lagrangian(N):
    """rows (x|z) in F_3^{2N}: random symplectic Gram-Schmidt."""
    rows = []
    def sp(u, v): return int((u[:N] @ v[N:] - u[N:] @ v[:N]) % 3)
    while len(rows) < N:
        v = rng.integers(0, 3, size=2*N)
        for r in rows:                       # make v commute with r using a partner trick: project out
            pass
        # brute force: accept v if it commutes with all rows and is independent
        if all(sp(v, r) == 0 for r in rows):
            M = np.array(rows + [v]) % 3
            # rank mod 3
            A = M.copy(); rk = 0
            for c in range(2*N):
                piv = [i for i in range(rk, len(A)) if A[i, c] % 3]
                if not piv: continue
                A[[rk, piv[0]]] = A[[piv[0], rk]]; A[rk] = (A[rk]*pow(int(A[rk, c]), -1, 3)) % 3
                for i in range(len(A)):
                    if i != rk and A[i, c] % 3: A[i] = (A[i] - A[i, c]*A[rk]) % 3
                rk += 1
            if rk == len(rows)+1: rows.append(v)
    return np.array(rows)
worst = np.inf; ncodes = 0
for N in range(2, 9):
    for _ in range(40 if N <= 6 else 8):
        G = random_lagrangian(N)
        coef = np.array(list(itertools.product(range(3), repeat=N)))          # 3^N x N
        H = (coef @ G) % 3                                                    # all elements of M
        Zc = ((H[:, :N] == 0) & (H[:, N:] == 0)).sum(axis=1)
        u = ((-0.5)**Zc).reshape([3]*N)
        uh = np.fft.fftn(u)                                                   # all characters of M
        U1 = uh.flat[0].real
        worst = min(worst, (U1 - np.abs(uh).max())/U1); ncodes += 1
print(f"(b) {ncodes} random Lagrangians, N = 2..8, all 3^N characters each: min of (u_hat(1) - max|u_hat(chi)|)/u_hat(1) = {worst:.3e} (>= 0 means the inequality holds)")
