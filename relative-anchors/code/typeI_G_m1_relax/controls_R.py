"""Negative controls for Conjecture R (exact, python-flint), at explicit rational weight vectors.

(C1) the couplings between the soft variables (a, c) and the letter-2 variables (rho, zeta) in Psi multiplied by
     1001/1000: Psi must become indefinite somewhere (no slack in the degenerate regimes);
(C2) the projection Pi omitted (Pi = I, the unprojected representative): Psi indefinite;
(C3) Lemma H's diagonal weight Lambda_h = 3/2 - rho_h replaced by the off-diagonal one, 5/12 (a valid but weaker
     relaxation): reported, no expectation;
(C4) the coupling of the letter-2 states to Q (the vector m of Acal = 2 rho^2 + |zeta|^2 - m^2/2) multiplied by
     1001/1000: Psi must become indefinite where the U-copies are light (there Acal_rho_rho = 2 - kappa^2/2 -> 0);
(C0) the unmodified Psi at the same points (all leading minors positive)."""
import flint
from psi_exact import PsiExact
from certify_points import leading_minors_bareiss

Q = flint.fmpq


def posdef(Mx):
    mins, _ = leading_minors_bareiss(Mx)
    return len(mins) == Mx.nrows() and all(x > 0 for x in mins)


points = []
for u in (4, 6):
    points.append((Q(1), [Q(1, 100) * Q(k + 3, k + 2) for k in range(u)]))          # S/sigma about 1e-4 * u
    points.append((Q(1), [Q(1, 1000) * Q(k + 3, k + 2) for k in range(u)]))
    points.append((Q(1, 100), [Q(k + 2, k + 1) for k in range(u)]))                 # S/sigma about 1e4 * u
    points.append((Q(1), [Q(k + 2, k + 1) for k in range(u)]))
ok = True
for wt, wU in points:
    u = len(wU)
    P = PsiExact(wt, wU)
    M0, iy = P.matrix()
    c0 = posdef(M0)
    n = M0.nrows()
    M1 = flint.fmpq_mat(M0)
    for i in range(iy):
        for j in range(iy, n):
            M1[i, j] = M0[i, j] * Q(1001, 1000)
            M1[j, i] = M0[j, i] * Q(1001, 1000)
    c1 = posdef(M1)
    P2 = PsiExact(wt, wU)
    P2.Pi = flint.fmpq_mat(u, u, [Q(1) if i == j else Q(0) for i in range(u) for j in range(u)])
    M2, _ = P2.matrix()
    c2 = posdef(M2)
    P3 = PsiExact(wt, wU)
    P3.lam = [[Q(5, 12)] * u for _ in range(u)]
    M3, _ = P3.matrix()
    c3 = posdef(M3)
    mv = [P.kappa] + [P.wt * P.w[h] / (P.sig + P.S) for h in range(u)]
    lam2 = Q(1001, 1000) ** 2
    M4 = flint.fmpq_mat(M0)
    for i in range(u + 1):
        for j in range(u + 1):
            M4[iy + i, iy + j] -= (lam2 - 1) * mv[i] * mv[j] / 2
    c4 = posdef(M4)
    ok &= c0
    print(f'u = {u}, w_t = {wt}, w_U = {[str(x) for x in wU]}: (C0) Psi > 0: {c0}; (C1) couplings x 1.001: {c1}; '
          f'(C2) no projection: {c2}; (C3) diagonal weight 5/12: {c3}; (C4) m x 1.001: {c4}', flush=True)
print('C0 PASS' if ok else 'C0 FAIL')
