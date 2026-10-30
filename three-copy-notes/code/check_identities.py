"""Exact identities for q_3(P_V), where P_V is the projector onto a plane V spanned by orthonormal u, v in
(C^d)^{(x)3} and q_3(C) = sum_J (-1/2)^{|J|} ||Tr_J C||_2^2.

  (a) q_3(P_V) = 3/2 - sum_l C_l,  C_l = (1/4)||alpha_l - beta_l||^2 + ||tau_l||^2      (any d)
      alpha_l, beta_l = one-site marginals of uu*, vv*;  tau_l = Tr_{other two sites}(u v*).
  (b) q_3(P_V) = 5/2 - sum_l Tr(psi_l^2) - ||u [x] v||^2                                   (d = 3)
      psi_l = one-site marginal of (uu* + vv*)/2,  (u [x] v)_{pqr} = sum eps_{pii'} eps_{qjj'} eps_{rkk'} u_{ijk} v_{i'j'k'}.
  (c) q_3(P_V) = (3/4) ( E_{x in V} sum_l S_L(x_l) - ||u [x] v||^2 )                         (d = 3)
      S_L(rho) = 1 - Tr rho^2, x Haar-random unit vector in V; the average is computed exactly from
      E[xx* (x) xx*] = (P_V (x) P_V)(I + F)/6.
Run: python check_identities.py"""
import itertools
import numpy as np

rng = np.random.default_rng(2026)


def q3(X, d):
    T = X.reshape([d] * 6)
    tot = 0.0
    for J in itertools.product([0, 1], repeat=3):
        A, m = T, 3
        for i in reversed(range(3)):
            if J[i]:
                A = np.trace(A, axis1=i, axis2=i + m)
                m -= 1
        M = A.reshape(d ** m, d ** m)
        tot += (-0.5) ** sum(J) * np.vdot(M, M).real
    return tot


def one_site(x, y, l, d):
    X, Y = x.reshape(d, d, d), y.reshape(d, d, d)
    ax = [i for i in range(3) if i != l]
    return np.tensordot(X, Y.conj(), axes=(ax, ax))


eps = np.zeros((3, 3, 3))
for p in itertools.permutations(range(3)):
    eps[p] = np.linalg.det(np.eye(3)[list(p)])


def random_plane(d):
    Q, _ = np.linalg.qr(rng.standard_normal((d ** 3, 2)) + 1j * rng.standard_normal((d ** 3, 2)))
    return Q[:, 0], Q[:, 1]


err_a = 0.0
for d in (2, 3, 4):
    for _ in range(30):
        u, v = random_plane(d)
        P = np.outer(u, u.conj()) + np.outer(v, v.conj())
        C = sum(0.25 * np.linalg.norm(one_site(u, u, l, d) - one_site(v, v, l, d)) ** 2
                + np.linalg.norm(one_site(u, v, l, d)) ** 2 for l in range(3))
        err_a = max(err_a, abs(q3(P, d) - (1.5 - C)))
print(f"(a) monogamy form, d = 2, 3, 4: max error {err_a:.1e}")

err_b = err_c = 0.0
F27 = None
for _ in range(50):
    u, v = random_plane(3)
    P = np.outer(u, u.conj()) + np.outer(v, v.conj())
    w = np.einsum('pai,qbj,rck,abc,ijk->pqr', eps, eps, eps, u.reshape(3, 3, 3), v.reshape(3, 3, 3))
    cross = np.vdot(w, w).real
    pur = sum(np.vdot(r, r).real for r in (0.5 * (one_site(u, u, l, 3) + one_site(v, v, l, 3)) for l in range(3)))
    err_b = max(err_b, abs(q3(P, 3) - (2.5 - pur - cross)))
    # exact Haar average of sum_l S_L(x_l): E[Tr rho_l(x)^2] = (||rho_l(P)||^2 + ||rho_{other two}(P)||^2) / 6
    avg = 0.0
    for l in range(3):
        r1 = one_site(u, u, l, 3) + one_site(v, v, l, 3)
        PP = P.reshape(3, 3, 3, 3, 3, 3)
        r2 = np.trace(PP, axis1=l, axis2=l + 3).reshape(9, 9)
        avg += 1 - (np.vdot(r1, r1).real + np.vdot(r2, r2).real) / 6
    err_c = max(err_c, abs(q3(P, 3) - 0.75 * (avg - cross)))
print(f"(b) cross-product form: max error {err_b:.1e}")
print(f"(c) entanglement form:  max error {err_c:.1e}")
