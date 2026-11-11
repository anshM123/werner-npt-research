"""Where the squeeze descent lives (orbit-notes/NOTES.md, section "Weak planes and the residual sandwich").

For a squeeze-critical rank-two C (||C|| = 1, r = R_t(C), Y = Phi_t(C) - r C with Phi_t = (id - t D)^(x)n):
  * cross second variation for H = A at copy s plus B at copy s' (G = exp(eps H)):
        X(A, B) = 4 Re B_r({A,C},{B,C}) + 4 Re<Y, {AB, C}> + 4 Re<Y, ACB + BCA>,   B_r(U,V) = <U, Phi_t(V) - r V>,
    the Hessian cross entry being X/2 and the diagonal entry S(A) = 2 q_r({A,C}) + 2 Re<Y, A^2 C + 2ACA + CA^2>;
  * weak plane K_s = span of the two smallest eigenvectors of mu_s = Tr_other(CC* + C*C);
  * paired directions H^mu = sigma_mu^(s) + (O sigma)_mu^(s'), sigma_mu the Pauli matrices of K_s in the mu_s
    eigenbasis; PT = sum_mu [S(sigma_mu^(s)) + S((O sigma)_mu^(s')) + X(sigma_mu, (O sigma)_mu)];
  * 'discrete': min over the 16 identifications O = diag(+-1,+-1,+-1) x {1, swap of sigma_x and sigma_y} and over pairs;
    'polar': O from the polar decomposition of the sandwich matrix N_mu,nu = 4 Re<Y, s_mu C s'_nu + s'_nu C s_mu>;
  * control: the same with 20 random planes instead of the weak planes (fraction with PT < 0).
Checks the cross formula against finite differences of R along exp(eps H) C exp(eps H) at the first file.
usage: python weak_plane_descent.py t file1.npy [file2.npy ...]"""
import sys, itertools
import numpy as np
from scipy.linalg import expm

T = float(sys.argv[1]); FILES = sys.argv[2:]
rng = np.random.default_rng(7)


def phi(C, n):
    Tn = C.reshape((3,) * (2 * n))
    for s in range(n):
        Tp = np.moveaxis(Tn, [s, n + s], [2 * n - 2, 2 * n - 1])
        Tp = Tp - T * np.trace(Tp, axis1=-2, axis2=-1)[..., None, None] * np.eye(3)
        Tn = np.moveaxis(Tp, [2 * n - 2, 2 * n - 1], [s, n + s])
    return Tn.reshape(C.shape)


def marg(M, s, n):
    Tn = M.reshape((3,) * (2 * n)); L = 'abcdefghij'; inp = [L[i] for i in range(2 * n)]
    for a in range(n):
        if a != s:
            inp[n + a] = inp[a]
    return np.einsum(''.join(inp) + '->' + L[s] + inp[n + s], Tn)


def emb(M, s, n):
    out = np.ones((1, 1), complex)
    for k in range(n):
        out = np.kron(out, M if k == s else np.eye(3))
    return out


def paulis(U):
    E = lambda i, j: np.outer(U[:, i], U[:, j].conj())
    return [E(0, 0) - E(1, 1), E(0, 1) + E(1, 0), 1j * (E(0, 1) - E(1, 0))]


OPTS = []
for sg in itertools.product((1, -1), repeat=3):
    for sw in (False, True):
        O = np.diag(sg).astype(float)
        OPTS.append(O @ np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]]) if sw else O)


def analyse(C, n, planes):
    r = np.vdot(C, phi(C, n)).real
    Pr = lambda Z: phi(Z, n) - r * Z
    Y = Pr(C); qr = lambda Z: np.vdot(Z, Pr(Z)).real; Br = lambda U, V: np.vdot(U, Pr(V))
    S = lambda a: 2 * qr(a @ C + C @ a) + 2 * np.vdot(Y, a @ a @ C + 2 * a @ C @ a + C @ a @ a).real
    X = lambda a, b: (4 * Br(a @ C + C @ a, b @ C + C @ b).real + 4 * np.vdot(Y, a @ b @ C + C @ a @ b).real
                      + 4 * np.vdot(Y, a @ C @ b + b @ C @ a).real)
    best_d = best_p = np.inf
    for s, u in itertools.combinations(range(n), 2):
        A = [emb(M, s, n) for M in paulis(planes[s])]; B = [emb(M, u, n) for M in paulis(planes[u])]
        SA = sum(S(a) for a in A); Xm = np.array([[X(a, b) for b in B] for a in A])
        for O in OPTS:
            Bp = [sum(O[nu, mu] * B[nu] for nu in range(3)) for mu in range(3)]
            best_d = min(best_d, SA + sum(S(b) for b in Bp) + sum(Xm[mu] @ O[:, mu] for mu in range(3)))
        Ns = np.array([[4 * np.vdot(Y, a @ C @ b + b @ C @ a).real for b in B] for a in A])
        Uo, _, Vt = np.linalg.svd(Ns); O = -(Uo @ Vt).T
        Bp = [sum(O[nu, mu] * B[nu] for nu in range(3)) for mu in range(3)]
        best_p = min(best_p, SA + sum(S(b) for b in Bp) + sum(Xm[mu] @ O[:, mu] for mu in range(3)))
    return r, best_d, best_p, X


for i, f in enumerate(FILES):
    C = np.load(f); C = C / np.linalg.norm(C); n = round(np.log(C.shape[0]) / np.log(3))
    weak = [np.linalg.eigh(marg(C @ C.conj().T + C.conj().T @ C, s, n))[1] for s in range(n)]
    r, pd, pp, X = analyse(C, n, weak)
    rand = [analyse(C, n, [np.linalg.qr(rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)))[0]
                           for _ in range(n)])[1] for _ in range(20)]
    line = (f"{f.replace(chr(92), '/').split('/')[-1]:34s} t={T} n={n} R={r:+.4f}  weak-plane PT: discrete {pd:+.4f},"
            f" polar {pp:+.4f};  random planes: {np.mean(np.array(rand) < 0):.2f} with PT<0")
    if i == 0:   # finite-difference check of the cross formula
        a8, b8 = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)), rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
        A = emb((a8 + a8.conj().T) / 2, 0, n); B = emb((b8 + b8.conj().T) / 2, 1, n)
        R = lambda Z: np.vdot(Z, phi(Z, n)).real / np.vdot(Z, Z).real
        g = lambda x, y: R(expm(x * A + y * B) @ C @ expm(x * A + y * B))
        h = 1e-4
        fd = (g(h, h) - g(h, -h) - g(-h, h) + g(-h, -h)) / (4 * h * h)
        line += f"\n    check: cross second derivative {fd:+.6f} vs X/2 = {X(A, B) / 2:+.6f}"
    print(line, flush=True)
