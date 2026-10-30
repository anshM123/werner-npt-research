"""Checks for the statements in the Remarks section of the paper.

1. The graded (Bernstein) identity. For every operator C on n sites and 0 <= t <= 1/2, with u = 1 - 2t,
       f_C(t) - (1-2t)(1-t)^(n-1)||C||^2 = sum_k u^(n-k) (1-u)^k [A_k(C) - C(n-1,k) 2^(-k) ||C||^2],
   where f_C(t) = sum_J (-t)^|J| ||Tr_J C||^2 and A_k(C) = sum_{m<=k} (-1/2)^m C(n-m,k-m) S_m(C).
   This is an algebraic identity, so random C of any rank are used.

2. The sharp rank-two bound for normal C:
       sum_i ||Tr_i C||^2 <= n ||C||^2 + |Tr C|^2 / 2.
   Random normal rank-two C are used, including product-type eigenvectors near the equality case.

3. The same bound for C = V M V^* (range of C = range of C^*), with V an isometry from C^2 and M an arbitrary
   2 x 2 matrix, and the additivity ||Tr_J C||^2 = ||Tr_J Z1||^2 + ||Tr_J Z2||^2 for C = Z1 + i Z2 (Z1, Z2 Hermitian).

4. The three-copy identity
       11 S_0 - 4 S_1 + S_2 = (5/6) 8 q_3 + (2/3) sum_i [2 S_0 + S_3/2 - ||Tr_i C||^2 - ||Tr_{jk} C||^2] + (1/6)(2 S_0 - S_3),
   {j,k} the complement of i (the bracket is the two-copy slack for the cut {i} | {j,k}); also that the three brackets
   are nonnegative for random rank-two C.

5. The sharp two-factor bound with a flatness term: for rank C <= r on H_1 (x) H_2,
       ||Tr_1 C||^2 + ||Tr_2 C||^2 <= ||C||^2 + (1 - 1/r) ||C||_1^2 + |Tr C|^2 / r
   (random and near-equality samples, r = 2, 3), and for r = 2 the sum-of-squares identity of the paper, with a
   balanced decomposition C = x1 y1^* + x2 y2^* (equal Gram matrices, <y_i, x_i> = Tr C / 2) built from the SVD.

6. The n-factor sharp bound with a flatness term for NORMAL rank-two C:
       sum_i ||Tr_i C||^2 <= (n-1) ||C||^2 + (||C||_1^2 + |Tr C|^2) / 2.

7. Its rank-one case for every n, as an operator inequality: sum_i rho_i (x) I <= (n-1) (Tr rho) I + rho for rho = xx^*,
   with rho_i the one-factor marginals (largest eigenvalue of the difference must be <= 0).

Run: python check_remarks.py
"""
import itertools
from math import comb
import numpy as np


def partial_trace_norms(C, d, n):
    """S[m] = sum over |J| = m of ||Tr_J C||_2^2."""
    T = C.reshape((d,) * (2 * n))
    letters = 'abcdefghijklmnop'
    S = [0.0] * (n + 1)
    for m in range(n + 1):
        for J in itertools.combinations(range(n), m):
            rows = list(letters[:n])
            cols = list(letters[n:2 * n])
            for k in J:
                cols[k] = rows[k]
            keep = [k for k in range(n) if k not in J]
            out = ''.join(rows[k] for k in keep) + ''.join(cols[k] for k in keep)
            R = np.einsum(''.join(rows) + ''.join(cols) + '->' + out, T)
            S[m] += np.sum(np.abs(R) ** 2)
    return S


def check_graded(rng, d=3, n=3, trials=20):
    worst = 0.0
    D = d ** n
    for _ in range(trials):
        C = rng.standard_normal((D, D)) + 1j * rng.standard_normal((D, D))
        S = partial_trace_norms(C, d, n)
        A = [sum((-0.5) ** m * comb(n - m, k - m) * S[m] for m in range(k + 1)) for k in range(n + 1)]
        for t in np.linspace(0, 0.5, 11):
            u = 1 - 2 * t
            lhs = sum((-t) ** m * S[m] for m in range(n + 1)) - (1 - 2 * t) * (1 - t) ** (n - 1) * S[0]
            rhs = sum(u ** (n - k) * (1 - u) ** k * (A[k] - comb(n - 1, k) * 2.0 ** (-k) * S[0])
                      for k in range(n + 1))
            worst = max(worst, abs(lhs - rhs) / S[0])
    return worst


def check_normal_T1(rng, d=3, n=3, trials=2000):
    D = d ** n
    worst = -np.inf
    for t in range(trials):
        if t % 3 == 0:
            c = rng.standard_normal(d) + 1j * rng.standard_normal(d)
            c /= np.linalg.norm(c)
            e = np.linalg.qr(rng.standard_normal((d, 2)) + 1j * rng.standard_normal((d, 2)))[0]
            Q = np.stack([np.kron(np.kron(e[:, k], c), c) for k in range(2)], 1)
        else:
            Q = np.linalg.qr(rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2)))[0]
        lam = rng.standard_normal(2) + 1j * rng.standard_normal(2)
        C = lam[0] * np.outer(Q[:, 0], Q[:, 0].conj()) + lam[1] * np.outer(Q[:, 1], Q[:, 1].conj())
        S = partial_trace_norms(C, d, n)
        excess = (S[1] - n * S[0] - abs(np.trace(C)) ** 2 / 2) / S[0]
        worst = max(worst, excess)
    return worst


def check_same_range(rng, d=3, n=3, trials=1000):
    D = d ** n
    worst, add_err = -np.inf, 0.0
    for t in range(trials):
        if t % 3 == 0:
            c = rng.standard_normal(d) + 1j * rng.standard_normal(d)
            c /= np.linalg.norm(c)
            e = np.linalg.qr(rng.standard_normal((d, 2)) + 1j * rng.standard_normal((d, 2)))[0]
            V = np.stack([np.kron(np.kron(e[:, k], c), c) for k in range(2)], 1)
            V = V + 0.005 * (rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2)))
            V = np.linalg.qr(V)[0]
        else:
            V = np.linalg.qr(rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2)))[0]
        M = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
        if t % 3 == 0:  # near the equality case P_U (x) cc^* (x) cc^*
            M = np.eye(2) + 0.05 * M
        C = V @ M @ V.conj().T
        S = partial_trace_norms(C, d, n)
        worst = max(worst, (S[1] - n * S[0] - abs(np.trace(C)) ** 2 / 2) / S[0])
        Z1, Z2 = (C + C.conj().T) / 2, (C - C.conj().T) / 2j
        S1, S2 = partial_trace_norms(Z1, d, n), partial_trace_norms(Z2, d, n)
        add_err = max(add_err, max(abs(S[m] - S1[m] - S2[m]) for m in range(n + 1)) / S[0])
    return worst, add_err


def cut_terms(C, d):
    """||Tr_i C||^2 and ||Tr_{jk} C||^2 for the three cuts of three factors."""
    T = C.reshape((d,) * 6)
    one = [np.sum(np.abs(np.trace(T, axis1=i, axis2=3 + i)) ** 2) for i in range(3)]
    two = []
    for i in range(3):
        j, k = [m for m in range(3) if m != i]
        R = np.trace(np.trace(T, axis1=j, axis2=3 + j), axis1=k - (1 if k > j else 0), axis2=3 + k - (1 if k > j else 0) - 1)
        two.append(np.sum(np.abs(R) ** 2))
    return one, two


def check_n3_identity(rng, d=3, trials=10):
    worst, worst_cut = 0.0, np.inf
    for t in range(2 * trials):
        if t % 2 == 0:
            C = rng.standard_normal((d ** 3, d ** 3)) + 1j * rng.standard_normal((d ** 3, d ** 3))
        else:
            X = rng.standard_normal((d ** 3, 2)) + 1j * rng.standard_normal((d ** 3, 2))
            Y = rng.standard_normal((d ** 3, 2)) + 1j * rng.standard_normal((d ** 3, 2))
            C = X @ Y.conj().T
        S = partial_trace_norms(C, d, 3)
        one, two = cut_terms(C, d)
        assert abs(sum(one) - S[1]) < 1e-9 * S[0] and abs(sum(two) - S[2]) < 1e-9 * S[0]
        q3 = S[0] - S[1] / 2 + S[2] / 4 - S[3] / 8
        cuts = [2 * S[0] + S[3] / 2 - one[i] - two[i] for i in range(3)]
        lhs = 11 * S[0] - 4 * S[1] + S[2]
        rhs = (5 / 6) * 8 * q3 + (2 / 3) * sum(cuts) + (1 / 6) * (2 * S[0] - S[3])
        worst = max(worst, abs(lhs - rhs) / S[0])
        if t % 2 == 1:
            worst_cut = min(worst_cut, min(cuts) / S[0])
    return worst, worst_cut


def two_factor_terms(C, d1, d2):
    T = C.reshape(d1, d2, d1, d2)
    return np.sum(np.abs(np.einsum('abac->bc', T)) ** 2) + np.sum(np.abs(np.einsum('abcb->ac', T)) ** 2)


def balanced_rank2(C, rng):
    U, s, Vh = np.linalg.svd(C)
    X = U[:, :2] * np.sqrt(s[:2]); Y = Vh[:2].conj().T * np.sqrt(s[:2])
    K = Y.conj().T @ X
    paulis = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
    K0 = K - np.trace(K) / 2 * np.eye(2)
    H1, H2 = (K0 + K0.conj().T) / 2, (K0 - K0.conj().T) / 2j
    h1 = np.real([np.trace(H1 @ P) / 2 for P in paulis]); h2 = np.real([np.trace(H2 @ P) / 2 for P in paulis])
    nv = np.cross(h1, h2)
    if np.linalg.norm(nv) < 1e-12:
        nv = np.cross(h1 + h2, rng.standard_normal(3))
    nv /= np.linalg.norm(nv)
    th, ph = np.arccos(np.clip(nv[2], -1, 1)), np.arctan2(nv[1], nv[0])
    q = np.array([np.cos(th / 2), np.exp(1j * ph) * np.sin(th / 2)])
    Q = np.stack([q, np.array([-np.conj(q[1]), np.conj(q[0])])], 1)
    return X @ Q, Y @ Q


def check_sharp_two_factor(rng, trials=600):
    worst, ident = -np.inf, 0.0
    for t in range(trials):
        d1, d2 = [(2, 2), (2, 3), (3, 3), (3, 2)][t % 4]
        r = 2 if t % 3 else 3
        D = d1 * d2
        if t % 5 == 0:  # near the equality case rho (x) |c><c| with rho >= 0 of rank r
            c = rng.standard_normal(d2) + 1j * rng.standard_normal(d2); c /= np.linalg.norm(c)
            A = rng.standard_normal((d1, min(r, d1))) + 1j * rng.standard_normal((d1, min(r, d1)))
            C = np.kron(A @ A.conj().T, np.outer(c, c.conj()))
            C = C + 1e-3 * ((rng.standard_normal((D, r)) + 1j * rng.standard_normal((D, r)))
                            @ (rng.standard_normal((D, r)) + 1j * rng.standard_normal((D, r))).conj().T)
            U, sv, Vh = np.linalg.svd(C); C = (U[:, :r] * sv[:r]) @ Vh[:r]
        else:
            C = (rng.standard_normal((D, r)) + 1j * rng.standard_normal((D, r))) @ (
                rng.standard_normal((D, r)) + 1j * rng.standard_normal((D, r))).conj().T
        n1 = np.linalg.svd(C, compute_uv=False).sum()
        N = np.linalg.norm(C) ** 2
        rhs = N + (1 - 1 / r) * n1 ** 2 + abs(np.trace(C)) ** 2 / r
        worst = max(worst, (two_factor_terms(C, d1, d2) - rhs) / N)
        if r == 2:
            X, Y = balanced_rank2(C, rng)
            x1, x2, y1, y2 = X[:, 0], X[:, 1], Y[:, 0], Y[:, 1]
            G = X.conj().T @ X
            def Pi(v):
                T = v.reshape(d1, d2, d1, d2)
                return (T - T.transpose(2, 1, 0, 3) - T.transpose(0, 3, 2, 1) + T.transpose(2, 3, 0, 1)) / 4
            sos = (0.5 * (G[0, 0].real - G[1, 1].real) ** 2
                   + 4 * sum(np.sum(np.abs(Pi(np.kron(a, b))) ** 2) for a, b in [(x1, y1), (x2, y2), (x1, x2), (y1, y2)])
                   + two_factor_terms(np.outer(x2, x1.conj()) - np.outer(y2, y1.conj()), d1, d2))
            ident = max(ident, abs(rhs - two_factor_terms(C, d1, d2) - sos) / N)
    return worst, ident


def check_normal_t1plus(rng, d=3, n=3, trials=1500):
    D = d ** n
    worst = -np.inf
    for t in range(trials):
        if t % 3 == 0:  # near the equality case rho (x) cc^* (x) cc^*, rho >= 0 of rank two
            c = rng.standard_normal(d) + 1j * rng.standard_normal(d); c /= np.linalg.norm(c)
            e = np.linalg.qr(rng.standard_normal((d, 2)) + 1j * rng.standard_normal((d, 2)))[0]
            Q = np.stack([np.kron(np.kron(e[:, k], c), c) for k in range(2)], 1)
            Q = np.linalg.qr(Q + 0.003 * (rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2))))[0]
            lam = rng.random(2) * np.exp(0.02j * rng.standard_normal(2))
        else:
            Q = np.linalg.qr(rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2)))[0]
            lam = rng.standard_normal(2) + 1j * rng.standard_normal(2)
        C = (Q * lam) @ Q.conj().T
        S = partial_trace_norms(C, d, n)
        rhs = (n - 1) * S[0] + 0.5 * (np.abs(lam).sum() ** 2 + abs(np.trace(C)) ** 2)
        worst = max(worst, (S[1] - rhs) / S[0])
    return worst


def check_rank_one_operator(rng, d=2, n=4, trials=300):
    D = d ** n
    worst = -np.inf
    for t in range(trials):
        if t % 3 == 0:  # product vectors: equality directions
            x = np.ones(1, complex)
            for _ in range(n):
                v = rng.standard_normal(d) + 1j * rng.standard_normal(d)
                x = np.kron(x, v / np.linalg.norm(v))
        else:
            x = rng.standard_normal(D) + 1j * rng.standard_normal(D)
            x /= np.linalg.norm(x)
        rho = np.outer(x, x.conj())
        T = rho.reshape((d,) * (2 * n))
        L = np.zeros((D, D), complex)
        for i in range(n):
            # one-factor marginal at i, tensored with identity on the other factors
            rows = list(range(n)); cols = list(range(n, 2 * n))
            letters = 'abcdefghijklmnop'
            r = [letters[k] for k in rows]; c = [letters[n + k] for k in rows]
            for k in range(n):
                if k != i:
                    c[k] = r[k]
            rho_i = np.einsum(''.join(r) + ''.join(c) + '->' + r[i] + c[i], T)
            L += np.kron(np.kron(np.eye(d ** i), rho_i), np.eye(d ** (n - i - 1)))
        M = L - (n - 1) * np.eye(D) - rho
        worst = max(worst, np.linalg.eigvalsh((M + M.conj().T) / 2)[-1])
    return worst


if __name__ == '__main__':
    rng = np.random.default_rng(2026)
    print('graded identity, max relative error (n = 3):', check_graded(rng))
    print('graded identity, max relative error (n = 4, d = 2):', check_graded(rng, d=2, n=4))
    w = check_normal_T1(rng)
    print('normal rank-two bound, max excess / ||C||^2 (must be <= 0):', w)
    assert w <= 1e-12
    w2, a2 = check_same_range(rng)
    print('C = V M V^*: max excess / ||C||^2 (must be <= 0):', w2, '  additivity error:', a2)
    assert w2 <= 1e-12 and a2 <= 1e-12
    e4, c4 = check_n3_identity(rng)
    print('three-copy identity, max relative error:', e4, '  smallest cut slack / ||C||^2 (rank two, must be >= 0):', c4)
    assert e4 <= 1e-12 and c4 >= -1e-12
    w5, i5 = check_sharp_two_factor(rng)
    print('sharp two-factor bound, max excess / ||C||^2 (must be <= 0):', w5, '  r = 2 identity error:', i5)
    assert w5 <= 1e-12 and i5 <= 1e-12
    w6 = check_normal_t1plus(rng)
    print('n-factor flatness bound for normal C, max excess / ||C||^2 (must be <= 0):', w6)
    assert w6 <= 1e-12
    w7 = check_rank_one_operator(rng)
    print('rank-one n-factor operator bound, largest eigenvalue of the difference (must be <= 0):', w7)
    assert w7 <= 1e-12
