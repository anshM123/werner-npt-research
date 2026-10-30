"""Numerical checks for the thickness margin and the plane-class bounds in the remarks of the paper.

1. Qubit lemma: for A >= 0 on C^2 (x) K with Tr_K A = I_2,
   (a)  lam_max(A) <= 1 + ||A_00 - A_11||_F / sqrt(2);
   (b)  ||A_00 - A_11||_F^2 >= 2 eta (eta + eps) for every unit psi with <psi|A|psi> = 1 + eta, eta >= 0,
        where eps = Tr((I (x) (I - E)) A) and E projects onto the K-support of psi.
   We check the explicit certificate of the proof (P is a rank-two projector with P psi = psi), the identity
   ||Y||^2 = 4 l (1 - l) for the smallest eigenvalue l of Z, both parts on random A, and the equality family of (a).
2. Pair lemma: for orthonormal w1, w2 in H1 (x) H2,
       4 ||Pi(w1 (x) w2)||^2 >= min(tau_1, tau_2) / 2   (sharp: orthogonal unitaries),
   the identity 4||Pi(w1 (x) w2)||^2 = 1 - Tr(rho_2^2)/2 + ||A1 - A2||^2/2 used in step (i), and the refined bound
   (1/2)(2 - mu_1) tau_2 + (mu_1 - 1) tau_1 when mu_1 > 1.
3. Proposition: q_2(C) >= min(tau_1, tau_2) ||C||^2 / (48 - 8 sqrt 2) for rank-two C, tau computed for the range.
4. Range and co-range in one plane (remark in the paper): the formula for q_2, the bounds (a) with lambda_3 and 1 - 1/sqrt2, (b) with
   max tau and (3 - sqrt5)/4, (c) for C >= 0 with tau and 1/4,
   the facts and product identities used in the proofs, the Re/Im identity, operators C = VMV*, and the examples behind
   the sharpness claims, including the one showing
   that tau cannot replace lambda_3 in (a) for d >= 5.

Plain NumPy; runs in about a minute.
"""
import numpy as np

rng = np.random.default_rng(2026)
sx = np.array([[0., 1.], [1., 0.]])
sz = np.diag([1., -1.])


def certificate_matrix(p, th):
    c, s = np.cos(th / 2), np.sin(th / 2)
    r = np.sqrt(p * (1 - p))
    k = (2 * p - 1) * c * s
    P = np.array([[c * c, -2 * r * c * s, k, 0],
                  [-2 * r * c * s, s * s, 0, k],
                  [k, 0, s * s, 2 * r * c * s],
                  [0, k, 2 * r * c * s, c * c]])
    psi = np.array([np.sqrt(p) * c, -np.sqrt(1 - p) * s, np.sqrt(p) * s, np.sqrt(1 - p) * c])
    Z = np.eye(2) / 2 + (p - 0.5) * np.sin(th) * sx
    Y = (np.cos(th) * sz - 2 * r * np.sin(th) * sx) / np.sqrt(2)
    return P, psi, Z, Y


def random_A(d, k):
    G = rng.standard_normal((2 * d, k)) + 1j * rng.standard_normal((2 * d, k))
    A = G @ G.conj().T
    Q = np.einsum('iaja->ij', A.reshape(2, d, 2, d))
    w, V = np.linalg.eigh(Q)
    T = np.kron(V @ np.diag(w ** -0.5) @ V.conj().T, np.eye(d))
    return T @ A @ T.conj().T


def check_qubit_lemma():
    err = 0.0
    for p in np.linspace(0, 1, 41):
        for th in np.linspace(0, np.pi, 41):
            P, psi, Z, Y = certificate_matrix(p, th)
            R = np.kron(Z, np.eye(2)) + np.kron(sz, Y) / np.sqrt(2)
            l = np.linalg.eigvalsh(Z)[0]
            err = max(err, np.abs(R - P).max(), np.abs(P @ P - P).max(), np.abs(P @ psi - psi).max(),
                      abs(np.trace(Z) - 1), max(0.0, np.linalg.norm(Y) - 1), max(0.0, -l),
                      abs(np.linalg.norm(Y) ** 2 - 4 * l * (1 - l)))
    print(f"1a. certificate: P = Z(x)I + sz(x)Y/sqrt2, P^2 = P, P psi = psi, Tr Z = 1, ||Y|| <= 1, Z >= 0, "
          f"||Y||^2 = 4l(1-l): max error {err:.1e}")
    worst_a, worst_b, nb = -np.inf, np.inf, 0
    for d in (2, 3, 4):
        for t in range(3000):
            A = random_A(d, rng.integers(1, 2 * d + 1))
            D = A[:d, :d] - A[d:, d:]
            ev, U = np.linalg.eigh(A)
            worst_a = max(worst_a, ev[-1] - 1 - np.linalg.norm(D) / np.sqrt(2))
            # part (b): top eigenvector and a perturbed (non-eigen) vector
            for psi in (U[:, -1], U[:, -1] + 0.3 * (rng.standard_normal(2 * d) + 1j * rng.standard_normal(2 * d))):
                psi = psi / np.linalg.norm(psi)
                eta = (psi.conj() @ A @ psi).real - 1
                if eta <= 0:
                    continue
                Kred = psi.reshape(2, d).T @ psi.reshape(2, d).conj()        # reduced operator on K
                w, V = np.linalg.eigh(Kred)
                E = V[:, w > 1e-12] @ V[:, w > 1e-12].conj().T
                eps = np.trace((np.eye(d) - E) @ (A[:d, :d] + A[d:, d:])).real
                worst_b = min(worst_b, np.linalg.norm(D) ** 2 - 2 * eta * (eta + eps))
                nb += 1
    print(f"1b. qubit lemma (a) on 9000 random A (d = 2, 3, 4): max of lam_max - 1 - ||A00 - A11||/sqrt2 = {worst_a:+.1e}")
    print(f"1c. qubit lemma (b) on {nb} pairs (A, psi): min of ||A00 - A11||^2 - 2 eta (eta + eps) = {worst_b:+.1e}")
    a, b = np.array([1., 0, 0]), np.array([0., 1, 0])
    phi = (np.kron([1, 0], a) + np.kron([0, 1], b)) / np.sqrt(2)
    phi2 = (np.kron([1, 0], b) - np.kron([0, 1], a)) / np.sqrt(2)
    gaps = []
    for lam in (1.0, 1.3, 1.7, 2.0):
        A = lam * np.outer(phi, phi) + (2 - lam) * np.outer(phi2, phi2)
        gaps.append(np.linalg.eigvalsh(A)[-1] - 1 - np.linalg.norm(A[:3, :3] - A[3:, 3:]) / np.sqrt(2))
    print(f"1d. equality family of (a): gaps {np.round(gaps, 14)}")


def pair_quantities(W1, W2):
    A1, A2 = W1 @ W1.conj().T, W2 @ W2.conj().T
    B1, B2 = W1.conj().T @ W1, W2.conj().T @ W2
    delta = 1 - np.linalg.norm(W1.conj().T @ W2) ** 2 - np.linalg.norm(W1 @ W2.conj().T) ** 2
    ea = np.sort(np.linalg.eigvalsh(A1 + A2))[::-1]
    eb = np.sort(np.linalg.eigvalsh(B1 + B2))[::-1]
    ident = 1 - 0.5 * np.sum(eb ** 2) + 0.5 * np.linalg.norm(A1 - A2) ** 2
    return delta, ea[2:].sum(), eb[2:].sum(), eb[0], abs(delta - ident)


def check_pair_lemma():
    worst, worst_add, ierr = np.inf, np.inf, 0.0
    for d in (3, 4, 5):
        for t in range(3000):
            Zm = rng.standard_normal((2, d * d)) + 1j * rng.standard_normal((2, d * d))
            if t % 3 == 1:
                Zm = 0.05 * Zm
                Zm[:, [0, 1, d, d + 1]] += rng.standard_normal((2, 4))
            Q, _ = np.linalg.qr(Zm.T)
            for W1, W2 in ((Q[:, 0].reshape(d, d), Q[:, 1].reshape(d, d)),
                           (Q[:, 0].reshape(d, d).T, Q[:, 1].reshape(d, d).T)):
                delta, ta, tb, mu1, e = pair_quantities(W1, W2)
                ierr = max(ierr, e)
                if min(ta, tb) > 1e-12:
                    worst = min(worst, delta / min(ta, tb))
                if mu1 > 1:
                    worst_add = min(worst_add, delta - 0.5 * (2 - mu1) * tb - (mu1 - 1) * ta)
    U1 = np.linalg.qr(rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)))[0]
    U2 = U1 @ np.diag([1, np.exp(2j * np.pi / 3), np.exp(4j * np.pi / 3)])
    delta, ta, tb, _, _ = pair_quantities(U1 / np.sqrt(3), U2 / np.sqrt(3))
    print(f"2.  pair lemma on 18000 random pairs (d = 3, 4, 5): min Delta/min(tau) = {worst:.4f} (claim >= 1/2); "
          f"refined bound min slack {worst_add:+.1e}; identity (i) error {ierr:.1e}; "
          f"orthogonal unitaries (d = 3): Delta = {delta:.6f}, min(tau)/2 = {min(ta, tb) / 2:.6f}")


def q2(C, d):
    T = C.reshape(d, d, d, d)
    return (np.linalg.norm(C) ** 2
            - 0.5 * (np.sum(abs(np.einsum('abac->bc', T)) ** 2) + np.sum(abs(np.einsum('abcb->ac', T)) ** 2))
            + 0.25 * abs(np.trace(C)) ** 2)


def check_theorem():
    const = 1 / (48 - 8 * np.sqrt(2))
    worst = np.inf
    for d in (3, 4):
        for t in range(4000):
            X = rng.standard_normal((d * d, 2)) + 1j * rng.standard_normal((d * d, 2))
            Y = rng.standard_normal((d * d, 2)) + 1j * rng.standard_normal((d * d, 2))
            if t % 2:
                X, Y = 0.05 * X, 0.05 * Y
                X[[0, d + 1], [0, 1]] += 1
                Y[[0, d + 1], [0, 1]] += 1
            C = X @ Y.conj().T
            C /= np.linalg.norm(C)
            U = np.linalg.svd(C)[0][:, :2]
            P = (U @ U.conj().T).reshape(d, d, d, d)
            ea = np.sort(np.linalg.eigvalsh(np.einsum('abcb->ac', P)))[::-1]
            eb = np.sort(np.linalg.eigvalsh(np.einsum('abac->bc', P)))[::-1]
            m = min(ea[2:].sum(), eb[2:].sum())
            if m > 1e-10:
                worst = min(worst, q2(C, d) / m)
    print(f"3.  proposition on 8000 random rank-two C (d = 3, 4): min q_2/(min tau ||C||^2) = {worst:.3f} "
          f"(claim >= {const:.5f})")


def plane_spectra(V, d):
    P = (V @ V.conj().T).reshape(d, d, d, d)
    ea = np.sort(np.linalg.eigvalsh(np.einsum('abcb->ac', P)))[::-1]
    eb = np.sort(np.linalg.eigvalsh(np.einsum('abac->bc', P)))[::-1]
    return ea, eb


def check_normal():
    """Plane class: the formula for q_2; PSD bounds with tau; (b) normal bound with lambda_3; the d >= 5 example."""
    ferr, wa1, wa2, wb = 0.0, np.inf, np.inf, np.inf
    for d in (3, 4, 5):
        for t in range(3000):
            Z = rng.standard_normal((d * d, 2)) + 1j * rng.standard_normal((d * d, 2))
            if t % 3 == 1:
                Z = 0.05 * Z
                Z[[0, d + 1], [0, 1]] += 1
            V = np.linalg.qr(Z)[0]
            v1, v2 = V[:, 0], V[:, 1]
            s2, phi = rng.random(), rng.random() * 2 * np.pi
            lam1, lam2 = 1.0, s2 * np.exp(1j * phi)
            C = lam1 * np.outer(v1, v1.conj()) + lam2 * np.outer(v2, v2.conj())
            W1, W2 = v1.reshape(d, d), v2.reshape(d, d)
            A1, A2 = W1 @ W1.conj().T, W2 @ W2.conj().T
            B1, B2 = W1.T @ W1.conj(), W2.T @ W2.conj()
            E1, E2 = 1 - np.trace(A1 @ A1).real, 1 - np.trace(A2 @ A2).real
            Dl = 1 - np.trace(A1 @ A2).real - np.trace(B1 @ B2).real
            formula = (0.25 + E1) + s2 ** 2 * (0.25 + E2) + s2 * np.cos(phi) * (Dl - 0.5)
            ferr = max(ferr, abs(q2(C, d) - formula))
            ea, eb = plane_spectra(V, d)
            m, M, m3 = min(ea[2:].sum(), eb[2:].sum()), max(ea[2:].sum(), eb[2:].sum()), min(ea[2], eb[2])
            Cp = np.outer(v1, v1.conj()) + s2 * np.outer(v2, v2.conj())
            nC, nCp = np.linalg.norm(C) ** 2, np.linalg.norm(Cp) ** 2
            if m > 1e-10:
                wa1 = min(wa1, q2(Cp, d) / (m * nCp))
            if M > 1e-10:
                wa2 = min(wa2, q2(Cp, d) / (M * nCp))
            if m3 > 1e-10:
                wb = min(wb, q2(C, d) / (m3 * nC))
    kap = 1 - 1 / np.sqrt(2)
    gold = (3 - np.sqrt(5)) / 4
    print(f"4a. the formula for q_2 on 9000 normal C: max error {ferr:.1e}")
    print(f"4b. plane class (c), PSD: min q2/(min tau ||C||^2) = {wa1:.3f} (claim >= 1/4); "
          f"plane class (b) on PSD: min q2/(max tau ||C||^2) = {wa2:.3f} (claim >= (3 - sqrt5)/4 = {gold:.4f})")
    print(f"4c. plane class (a), normal: min q2/(min lambda_3 ||C||^2) = {wb:.3f} (claim >= 1 - 1/sqrt2 = {kap:.4f})")
    # 4e. the facts used in (a), (b), the two product identities, the Re/Im identity, and C = V M V^*, M arbitrary
    s_, e1_, e2_ = rng.random(20000) * 0.5, rng.random(20000) * 0.5, rng.random(20000) * 0.5
    s_ = e1_ + e2_
    idA = max(abs((0.25 + e1_ - kap * s_) * (0.25 + e2_ - kap * s_) - 1 / 16
                  - (np.sqrt(2) - 1) / 4 * s_ * (1 - 2 * s_) - e1_ * e2_))
    E1_, E2_ = rng.random(20000), rng.random(20000)
    S_ = E1_ + E2_
    idB = max(abs((0.25 + E1_ - gold * S_) * (0.25 + E2_ - gold * S_) - 1 / 16
                  - gold * (1 - gold) * S_ * (1 - S_) - E1_ * E2_))
    prodA = min((0.25 + e1_ - kap * s_) * (0.25 + e2_ - kap * s_)) - 1 / 16
    prodB = min((0.25 + E1_ - gold * S_) * (0.25 + E2_ - gold * S_)) - 1 / 16
    print(f"4f. product identities: errors {idA:.1e} ((a)), {idB:.1e} ((b)); min of product - 1/16: {prodA:+.2e} "
          f"(a, e_j in [0,1/2]), {prodB:+.2e} (b, E_j in [0,1))")
    fmin, wv, wvM, rerr = np.inf, np.inf, np.inf, 0.0
    for d in (3, 4, 5):
        for t in range(2000):
            V = np.linalg.qr(rng.standard_normal((d * d, 2)) + 1j * rng.standard_normal((d * d, 2)))[0]
            ea, eb = plane_spectra(V, d)
            m3 = min(ea[2], eb[2])
            xs, Es = [], []
            for j in range(2):
                r = np.sort(np.linalg.eigvalsh(V[:, j].reshape(d, d) @ V[:, j].reshape(d, d).conj().T))[::-1]
                xs.append(r[1]); Es.append(1 - (r ** 2).sum())
            fmin = min(fmin, Es[0] - 2 * xs[0] * (1 - xs[0]), Es[1] - 2 * xs[1] * (1 - xs[1]), xs[0] + xs[1] - m3)
            M = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
            C = V @ M @ V.conj().T
            X, Y = (C + C.conj().T) / 2, (C - C.conj().T) / 2j
            rerr = max(rerr, abs(q2(C, d) - q2(X, d) - q2(Y, d)))
            if m3 > 1e-10:
                wv = min(wv, q2(C, d) / (m3 * np.linalg.norm(C) ** 2))
            wvM = min(wvM, q2(C, d) / (max(ea[2:].sum(), eb[2:].sum()) * np.linalg.norm(C) ** 2))
    print(f"4e. facts E_j >= 2x_j(1-x_j), lambda_3 <= x1 + x2: min slack {fmin:+.1e}; Re/Im identity error {rerr:.1e}; "
          f"C = VMV*: min q2/(min lambda_3 ||C||^2) = {wv:.3f} (claim >= {kap:.4f}), "
          f"min q2/(max tau ||C||^2) = {wvM:.3f} (claim >= {gold:.4f})")
    out = []
    for d in (3, 4, 5, 8):
        e = np.eye(d)
        v1 = np.kron(e[0], e[0])
        v2 = sum(np.kron(e[i], e[i]) for i in range(1, d)) / np.sqrt(d - 1)
        ea, _ = plane_spectra(np.stack([v1, v2], 1), d)
        r = min(q2(np.outer(v1, v1) - t * np.outer(v2, v2), d) / (1 + t * t) for t in np.linspace(0, 1, 4001))
        e_ = (d - 2) / (d - 1)
        closed = 0.25 * (1 + 2 * e_ - np.sqrt(1 + 4 * e_ ** 2))
        Cp = np.outer(v1, v1) + 1e-5 * np.outer(v2, v2)
        rp = q2(Cp, d) / (ea[2:].sum() * np.linalg.norm(Cp) ** 2)
        out.append(f"d={d}: tau-ratio {r / ea[2:].sum():.4f} (closed form {closed / e_:.4f}), lambda_3-ratio "
                   f"{r / ea[2]:.4f}, PSD s->0 ratio {rp:.4f} vs (d-1)/(4(d-2)) = {(d - 1) / (4 * (d - 2)):.4f}")
    print("4d. examples v1 = |00>, v2 = sum_i |ii>/sqrt(d-1); C = v1v1* - t v2v2* and C = v1v1* + s v2v2*:\n    "
          + "\n    ".join(out) + f"\n    (limit of the tau-ratio as d -> infinity: (3 - sqrt5)/4 = {gold:.4f})")


if __name__ == '__main__':
    check_qubit_lemma()
    check_pair_lemma()
    check_theorem()
    check_normal()
