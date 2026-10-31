"""Section 10: antisymmetric and code blocks.

   (a) Exact (SymPy) partial-transpose spectrum of the one-sided compression (P_S (x) I) rho^{(x)b} (P_S (x) I) for
       S = Lambda (the antisymmetric subspace of two copies, b = 2) and S = Q (the qutrit [[3,1,2]] code, b = 3).
   (b) The same smallest eigenvalue in floating point for a few larger candidate blocks (Section 10, "Larger blocks").
   (c) Dense checks of the product bound for random operators with ranges in products of local-unitary images of the
       blocks (n = 2, 4, 5).
   (d) Proposition 11: the identity for the compression to Lambda, the 2-design decomposition of I + |Phi'><Phi'| over
       the 12 vectors of four mutually unbiased bases, the resulting explicit separable decomposition, and a dense
       check at n = 4.
   (e) Proposition 12: the compression to Q is block diagonal over the nine syndrome sectors of the second side; each
       block is the identity plus explicit product terms (MUB bases), so it is separable. Dense checks of removal."""
import itertools
import numpy as np
import sympy as sp
from werner import q

R1 = sp.zeros(9, 9)
for i in range(3):
    for j in range(3):
        R1[3 * i + j, 3 * i + j] += 1
        R1[3 * i + j, 3 * j + i] -= sp.Rational(1, 2)


def onesided_pt_exact(V, b):
    """V: list of k orthonormal vectors in (C^3)^{(x)b} (exact). Returns the partial transpose (on Bob) of the one-sided
       compression of rho^{(x)b} to span(V), as a (k 3^b) x (k 3^b) exact matrix."""
    N = 3 ** b; k = len(V); idx = list(itertools.product(range(3), repeat=b))

    def rb(a, bb, a2, b2):
        val = sp.Integer(1)
        for s in range(b):
            val *= R1[3 * a[s] + bb[s], 3 * a2[s] + b2[s]]
            if val == 0:
                return 0
        return val
    M = sp.zeros(k * N, k * N)
    for i, j in itertools.product(range(k), repeat=2):
        for (qi, qq), (pi, pp) in itertools.product(enumerate(idx), repeat=2):
            tot = 0
            for ai, a in enumerate(idx):
                if V[i][ai] == 0:
                    continue
                for a2i, a2 in enumerate(idx):
                    if V[j][a2i] != 0:
                        tot += sp.conjugate(V[i][ai]) * V[j][a2i] * rb(a, qq, a2, pp)
            M[i * N + pi, j * N + qi] = tot          # partial transpose: Bob's q and p exchanged
    return M


def vec(b, terms, norm):
    idx = list(itertools.product(range(3), repeat=b)); v = [sp.Integer(0)] * len(idx)
    for coef, ks in terms:
        v[idx.index(tuple(ks))] += sp.Integer(coef) / norm
    return v


Lam = [vec(2, [(1, (i, j)), (-1, (j, i))], sp.sqrt(2)) for i, j in ((0, 1), (0, 2), (1, 2))]
Q = [vec(3, [(1, (k, (k + j) % 3, (k + 2 * j) % 3)) for k in range(3)], sp.sqrt(3)) for j in range(3)]
print("(a) exact partial-transpose spectra of one-sided compressions:")
for name, V, b in (("Lambda (antisymmetric subspace of two copies)", Lam, 2), ("Q (qutrit [[3,1,2]] code, three copies)", Q, 3)):
    ev = onesided_pt_exact(V, b).eigenvals()
    print(f"  {name}: " + ", ".join(f"{sp.nsimplify(e)} (x{m})" for e, m in sorted(ev.items(), key=lambda t: float(t[0]))), flush=True)

# (b) floating point, larger candidates
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2


def tens(c):
    R = np.array([[1.0]])
    for _ in range(c):
        R = np.kron(R, rho)
    perm = list(range(0, 2 * c, 2)) + list(range(1, 2 * c, 2))
    return R.reshape([3] * (4 * c)).transpose(perm + [p + 2 * c for p in perm]).reshape(9 ** c, 9 ** c)


R2, R3 = tens(2), tens(3)


def ptmin(V, R, N):
    V = np.linalg.qr(V)[0]; k = V.shape[1]
    W = np.kron(V, np.eye(N)); s = W.conj().T @ R @ W
    s = s.reshape(k, N, k, N).transpose(0, 3, 2, 1).reshape(k * N, k * N)
    return np.linalg.eigvalsh((s + s.conj().T) / 2)[0]


e = np.eye(3)
k2 = lambda i, j: np.kron(e[i], e[j])
k3 = lambda i, j, k: np.kron(k2(i, j), e[k])
L2 = np.array([(k2(i, j) - k2(j, i)) / 2 ** .5 for i, j in ((0, 1), (0, 2), (1, 2))]).T
Qn = np.array([sum(k3(k, (k + j) % 3, (k + 2 * j) % 3) for k in range(3)) / 3 ** .5 for j in range(3)]).T
eps = sum(np.linalg.det(np.eye(3)[list(p)]) * k3(*p) for p in itertools.permutations(range(3))) / 6 ** .5
print("\n(b) smallest partial-transpose eigenvalue of other candidate blocks:")
for name, V, R, N in (("Lambda + |00>", np.column_stack([L2, k2(0, 0)]), R2, 9),
                      ("Lambda + sum_i |ii>/sqrt3", np.column_stack([L2, sum(k2(i, i) for i in range(3)) / 3 ** .5]), R2, 9),
                      ("Lambda + (|01> + |10>)/sqrt2", np.column_stack([L2, (k2(0, 1) + k2(1, 0)) / 2 ** .5]), R2, 9),
                      ("Sym^2(C^3)", np.array([k2(i, i) for i in range(3)] + [(k2(i, j) + k2(j, i)) / 2 ** .5 for i, j in ((0, 1), (0, 2), (1, 2))]).T, R2, 9),
                      ("Q + totally antisymmetric vector", np.column_stack([Qn, eps]), R3, 27),
                      ("C^3 (x) |00> (three copies)", np.kron(np.eye(3), k2(0, 0)[:, None]), R3, 27)):
    print(f"  {name}: {ptmin(V, R, N):+.4f}", flush=True)

# (c) dense checks
rng = np.random.default_rng(12)
cr = lambda *s: rng.standard_normal(s) + 1j * rng.standard_normal(s)


def haar(d):
    Qm, Rm = np.linalg.qr(cr(d, d)); return Qm * (np.diag(Rm) / abs(np.diag(Rm)))


print("\n(c) dense q_n / ||C||^2 for random C with range in a product of blocks (any rank, arbitrary co-range):")
for n, blocks, bound in ((2, ["L"], 1 / 4), (4, ["L", "L"], 1 / 16), (5, ["L", "Q"], 1 / 8)):
    worst = np.inf
    for t in range(40 if n < 5 else 15):
        S = np.eye(1)
        for bl in blocks:
            S = np.kron(S, np.kron(haar(3), haar(3)) @ L2 if bl == "L" else np.kron(np.kron(haar(3), haar(3)), haar(3)) @ Qn)
        r = rng.integers(1, S.shape[1] + 1)
        C = S @ cr(S.shape[1], r) @ cr(r, 3 ** n)
        worst = min(worst, q(C, n, 0.5) / np.vdot(C, C).real)
    print(f"  n = {n}, blocks {blocks}: smallest value {worst:.4f} (proved bound {bound:.4f})", flush=True)

# (d) Proposition 11: the compression to Lambda is separable, with an explicit decomposition
print("\n(d) separability of the compression to Lambda (Proposition 11):")
cyc = {0: (1, 2), 1: (2, 0), 2: (0, 1)}
A = np.array([(k2(*cyc[z]) - k2(*cyc[z][::-1])) / 2 ** .5 for z in range(3)]).T          # a_z = sum eps_{zxy}|xy>/sqrt2
Wl = np.kron(A, np.eye(9)); sig = Wl.T @ R2 @ Wl                                          # on Lambda (x) (B1 B2)
Sw = np.eye(9)[[3 * (k % 3) + k // 3 for k in range(9)]]
Psym, Panti = (np.eye(9) + Sw) / 2, (np.eye(9) - Sw) / 2
Phi = np.eye(3).reshape(-1)
D = np.kron(np.outer(Phi, Phi), np.eye(3))
form = 0.75 * np.kron(np.eye(3), Panti) + 0.5 * np.kron(np.eye(3), Psym) @ (np.eye(27) + D) @ np.kron(np.eye(3), Psym)
print(f"  sigma_Lambda = 3/4 I (x) P_anti + 1/2 (I (x) P_sym)(I + D)(I (x) P_sym): max deviation {np.abs(sig - form).max():.1e}")
w = np.exp(2j * np.pi / 3)
mub = [np.eye(3)[b] for b in range(3)] + [np.array([w ** (a * k * k + b * k) for k in range(3)]) / 3 ** .5 for a in range(3) for b in range(3)]
Ssum = sum(np.kron(np.outer(u, u.conj()), np.outer(u.conj(), u)) for u in mub)
print(f"  sum over the 12 MUB vectors of |u><u| (x) |conj u><conj u| = I + |Phi'><Phi'|: max deviation {np.abs(Ssum - np.eye(9) - np.outer(Phi, Phi)).max():.1e}")
sep = 0.75 * np.kron(np.eye(3), Panti) + 0.5 * sum(np.kron(np.outer(u, u.conj()), Psym @ np.kron(np.outer(u.conj(), u), np.eye(3)) @ Psym) for u in mub)
print(f"  explicit separable decomposition reproduces sigma_Lambda: max deviation {np.abs(sep - sig).max():.1e}")
worst = np.inf
for t in range(100):
    B = np.kron(np.kron(haar(3), haar(3)) @ L2, np.eye(9))
    C = B @ cr(27, 2) @ cr(2, 81)
    worst = min(worst, q(C, 4, 0.5) / np.vdot(C, C).real)
print(f"  n = 4, rank two, range at copies 1, 2 in a local-unitary image of Lambda: smallest q_4/||C||^2 = {worst:.4f} (>= 0)")

# (e) Proposition 12: the compression to Q is separable (sector-wise explicit decomposition), triple removal
print("\n(e) separability of the compression to Q (Proposition 12):")
sQ = np.kron(Qn, np.eye(27)).T @ R3 @ np.kron(Qn, np.eye(27))
Xp = np.roll(np.eye(3), 1, axis=0); Zp = np.diag([1, w, w ** 2])
Pauli = lambda a, b: np.linalg.matrix_power(Xp, a) @ np.linalg.matrix_power(Zp, b)
mubs = [np.eye(3)] + [np.array([[w ** (a * k * k + b * k) for b in range(3)] for k in range(3)]) / 3 ** .5 for a in range(3)]
F9 = np.eye(9)[[3 * (k % 3) + k // 3 for k in range(9)]]
total = np.zeros((81, 81), complex); nterms = 0
for a, b in itertools.product(range(3), repeat=2):
    Bs = np.kron(np.kron(Pauli(a, b), np.eye(3)), np.eye(3)) @ Qn           # sector a,b of the second side: (X^a Z^b (x) I (x) I) Q
    emb = np.kron(np.eye(3), Bs); M = emb.conj().T @ sQ @ emb
    if (a, b) == (0, 0):
        dec = 3 / 8 * np.eye(9) + 1 / 8 * sum(np.kron(np.outer(u, u.conj()), np.outer(u, u.conj())) for mb in mubs for u in mb.T)
        print(f"  sector (0,0): M = I/2 + F/8: {np.allclose(M, 0.5 * np.eye(9) + F9 / 8)}; = 3/8 I + 1/8 sum over 12 MUB vectors |uu><uu|: deviation {np.abs(dec - M).max():.1e}")
        nterms += 12
    else:
        Pi = 4 * (M - 0.5 * np.eye(9)); dec = None
        for Bl, Bb in itertools.product(mubs, [m for mb in mubs for m in (mb, mb.conj())]):
            G = np.kron(Bl, Bb).conj().T @ Pi @ np.kron(Bl, Bb)
            if np.abs(G - np.diag(np.diag(G))).max() < 1e-10:
                d = np.round(np.diag(G).real).astype(int)
                dec = 0.5 * np.eye(9) + 0.25 * sum(np.kron(np.outer(Bl[:, x], Bl[:, x].conj()), np.outer(Bb[:, y], Bb[:, y].conj()))
                                                   for x in range(3) for y in range(3) if d[3 * x + y] == 1)
                break
        nterms += 3
    total += emb @ dec @ emb.conj().T
print(f"  sigma_Q = sum over the 9 sectors of (identity + explicit product terms): deviation {np.abs(total - sQ).max():.1e} ({nterms} product terms besides identities)")
for n, blocks, extra, why in ((4, ["Q"], 1, "EP_1"), (5, ["Q"], 2, "EP_2"), (5, ["L"], 3, "EP_3")):
    worst = np.inf
    for t in range(25):
        S = np.eye(1)
        for bl in blocks:
            S = np.kron(S, np.kron(haar(3), haar(3)) @ L2 if bl == "L" else np.kron(np.kron(haar(3), haar(3)), haar(3)) @ Qn)
        S = np.kron(S, np.eye(3 ** extra))
        C = S @ cr(S.shape[1], 2) @ cr(2, 3 ** n)
        worst = min(worst, q(C, n, 0.5) / np.vdot(C, C).real)
    print(f"  n = {n}, rank two, range in {blocks} (x) (C^3)^(x){extra} (free copies may be entangled with the block): "
          f"smallest q_n/||C||^2 = {worst:.4f} (>= 0 by removal and {why})", flush=True)
