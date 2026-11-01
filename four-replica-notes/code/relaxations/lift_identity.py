"""Check (before twirling, so for every physical pair) that the covariance objects of a rank-2 operator C are
compressions of whole-replica partial transposes of the lifted moments Y_f of Tomar's lifting.
Lifting: ||C||_F^2 = 1/2, A isometry Q -> H with ran C in ran A, B = 2 C^* A,
a = sum_{x,r} A_{xr} |x>|r> / sqrt2, b = sum B_{xr} |x>|r> / sqrt2 in V = H (x) Q, so that C = Tr_Q |a><b|.
Claim: vec C = K^*(a (x) conj b) with K = I_H (x) |Omega_Q> (x) I_H on the Q factors of the pair, hence
   rho_x (x) rho_y (x) |vec C><vec C| = (Tr_{Q1 Q2} (x) K^*) [ |x y a conj(b)><x y a conj(b)| ] (... K)
and |x y a conj(b)><x y a conj(b)| = (Y_{x y a b})^{Gamma_4} is a whole-replica partial transpose of a moment family.
usage: python lift_identity.py [n]"""
import sys, itertools, numpy as np

n = int(sys.argv[1]) if len(sys.argv) > 1 else 1
rng = np.random.default_rng(3)
dH = 3 ** n; dQ = 2
# random rank-2 C with ||C||_F^2 = 1/2
C = (rng.standard_normal((dH, 2)) + 1j * rng.standard_normal((dH, 2))) @ (rng.standard_normal((2, dH)) + 1j * rng.standard_normal((2, dH)))
C *= np.sqrt(0.5) / np.linalg.norm(C)
Uc, s, Vh = np.linalg.svd(C)
A = Uc[:, :2]                                   # isometry Q -> H, ran C in ran A
B = 2 * C.conj().T @ A
a = (A / np.sqrt(2)).reshape(dH * dQ)           # index (x, r)
b = (B / np.sqrt(2)).reshape(dH * dQ)
Cr = np.einsum('xr,yr->xy', a.reshape(dH, dQ), b.reshape(dH, dQ).conj())
print(f"n={n}: ||a|| = {np.linalg.norm(a):.12f}, ||b|| = {np.linalg.norm(b):.12f}, |C - Tr_Q|a><b|| = {np.abs(Cr - C).max():.1e}, "
      f"|Tr_H |a><a| - I/2| = {np.abs(np.einsum('xr,xs->rs', a.reshape(dH, dQ), a.reshape(dH, dQ).conj()) - np.eye(2) / 2).max():.1e}")
vecC = C.reshape(-1)                            # vec C = sum C_xy |x>|y>
# K^*(a (x) conj b): contract Q of the third replica with Q of the fourth replica by <Omega|
kb = np.einsum('xr,yr->xy', a.reshape(dH, dQ), b.reshape(dH, dQ).conj()).reshape(-1)
print(f"   |vec C - K^*(a (x) conj b)| = {np.abs(kb - vecC).max():.1e}")
letters = {'a': a.reshape(dH, dQ), 'b': b.reshape(dH, dQ)}
for x, y in [('a', 'a'), ('a', 'b'), ('b', 'b')]:
    # left side: rho_x (x) rho_y (x) |vec C><vec C| with rho = Tr_Q |x><x|
    rx = np.einsum('xr,yr->xy', letters[x], letters[x].conj()); ry = np.einsum('xr,yr->xy', letters[y], letters[y].conj())
    lhs = np.kron(np.kron(rx, ry), np.outer(vecC, vecC.conj()))
    # right side: vector u = x (x) y (x) a (x) conj(b) in V^4, compressed: Q1, Q2 traced (as a sum over basis
    # vectors), Q3 Q4 contracted with Omega
    rhs = np.zeros_like(lhs)
    for q1 in range(dQ):
        for q2 in range(dQ):
            v = np.einsum('i,j,kl->ijkl', letters[x][:, q1], letters[y][:, q2],
                          np.einsum('kr,lr->kl', letters['a'], letters['b'].conj())).reshape(-1)
            rhs += np.outer(v, v.conj())
    fam = ''.join(sorted(x + y + 'ab'))
    print(f"   letters ({x},{y}): |lhs - compression of (Y_{{{x}{y}ab}})^Gamma_4| = {np.abs(lhs - rhs).max():.1e}  "
          f"[family {fam.count('a')}{fam.count('b')}, transposed replica carries letter b]")
