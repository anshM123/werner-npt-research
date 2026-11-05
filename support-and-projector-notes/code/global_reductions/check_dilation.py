"""Projector dilation: q_{k+r}(P) = 2 q_k(C) + 2^{-r} q_k(D0 + (-1)^r D1), and the projector choice P = lambda Pi."""
import numpy as np
from qtools import q_dense, rand_c

rng = np.random.default_rng(7)
worst = 0.0
for k in (1, 2):
    D = 3 ** k
    for trial in range(4):
        # a random rank-two C, written with an arbitrary factorization
        ls = [rand_c(D, rng) for _ in range(2)]
        rs = [rand_c(D, rng) for _ in range(2)]
        C = sum(np.outer(l, r.conj()) for l, r in zip(ls, rs))
        for r in (1, 2, 3):
            f0 = np.zeros(3 ** r); f0[0] = 1.0                      # 0^r
            f1 = np.zeros(3 ** r); f1[sum(3 ** i for i in range(r))] = 1.0  # 1^r
            ws = [np.kron(l, f0) + np.kron(rr, f1) for l, rr in zip(ls, rs)]
            P = sum(np.outer(w, w.conj()) for w in ws)
            D0 = sum(np.outer(l, l.conj()) for l in ls)
            D1 = sum(np.outer(rr, rr.conj()) for rr in rs)
            lhs = q_dense(P, k + r)
            rhs = 2 * q_dense(C, k) + 2.0 ** (-r) * q_dense(D0 + (-1) ** r * D1, k)
            worst = max(worst, abs(lhs - rhs) / max(1.0, abs(lhs)))
        # projector choice from the singular value decomposition
        U, sig, Vh = np.linalg.svd(C)
        lam = 2 * sig[0]
        al = [(lam + np.sqrt(lam ** 2 - 4 * s ** 2)) / 2 for s in sig[:2]]
        be = [(lam - np.sqrt(lam ** 2 - 4 * s ** 2)) / 2 for s in sig[:2]]
        ls2 = [np.sqrt(al[i]) * U[:, i] for i in range(2)]
        rs2 = [np.sqrt(be[i]) * Vh[i].conj() for i in range(2)]
        C2 = sum(np.outer(l, r.conj()) for l, r in zip(ls2, rs2))
        assert np.allclose(C2, C)
        for r in (1, 2, 3):
            f0 = np.zeros(3 ** r); f0[0] = 1.0
            f1 = np.zeros(3 ** r); f1[sum(3 ** i for i in range(r))] = 1.0
            ws = [np.kron(l, f0) + np.kron(rr, f1) for l, rr in zip(ls2, rs2)]
            G = np.array([[np.vdot(a, b) for b in ws] for a in ws])
            assert np.allclose(G, lam * np.eye(2)), G
            Pi = sum(np.outer(w, w.conj()) for w in ws) / lam
            assert np.allclose(Pi @ Pi, Pi)
            D0 = sum(np.outer(l, l.conj()) for l in ls2)
            D1 = sum(np.outer(rr, rr.conj()) for rr in rs2)
            lhs = q_dense(Pi, k + r)
            rhs = (2 * q_dense(C, k) + 2.0 ** (-r) * q_dense(D0 + (-1) ** r * D1, k)) / lam ** 2
            worst = max(worst, abs(lhs - rhs))
            # error term bound |2^{-r} q_k(D0 +- D1)| / lam^2 <= 4 * 2^{-r}
            assert abs(2.0 ** (-r) * q_dense(D0 + (-1) ** r * D1, k)) / lam ** 2 <= 4 * 2.0 ** (-r) + 1e-12
            # and 2 ||C||^2 / lam^2 >= 1/2
            assert 2 * np.linalg.norm(C) ** 2 / lam ** 2 >= 0.5 - 1e-12
print('max deviation', worst)
print('PASS' if worst < 1e-9 else 'FAIL')
