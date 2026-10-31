"""Section 2 checks (occupation-sector positivity).
   (a) The structure used in the proof. Pair space C^9 with basis f = (|00>+|11>)/sqrt2, g = |22>, s1..s7.
       On Sym^n(C^9), restricted to a fixed Alice z-number, B^{(x)n} (B = I - 3t|phi><phi|) is diagonal in the
       occupation basis with entries d(a,b) = sum_j C(a,j) C(b,j) (1-2t)^{a-j} (2t^2)^j (1-t)^{b-j}  (n <= 3).
   (b) q_n(t,C) = <<C| B^{(x)n} |C>> for a random C on (C^3)^{(x)3}.
   (c) Minimum eigenvalue of the q-form on Hom(Sym^n, E_l) (range in one occupation sector), n <= 5.
   (d) Block patterns: at most one block per diagonal l - m gives a PSD form; a repeated diagonal can be negative."""
import itertools
import numpy as np
from math import comb, sqrt
from bosonic import basis, q_and_grad, q_direct

def d(a, b, t):
    return sum(comb(a, j) * comb(b, j) * (1 - 2 * t) ** (a - j) * (2 * t * t) ** j * (1 - t) ** (b - j)
               for j in range(min(a, b) + 1))

def ket(i, j):
    v = np.zeros(9); v[3 * i + j] = 1; return v

s2 = 1 / sqrt(2)
modes = [s2 * (ket(0, 0) + ket(1, 1)), ket(2, 2), s2 * (ket(0, 0) - ket(1, 1)),
         ket(0, 1), ket(1, 0), ket(0, 2), ket(1, 2), ket(2, 0), ket(2, 1)]
U = np.array(modes).T
phi = sum(ket(i, i) for i in range(3)) / sqrt(3)

print("(a) diagonal structure on fixed Alice z-number")
for t in (0.5, 0.41):
    B = np.eye(9) - 3 * t * np.outer(phi, phi)
    worst_off = worst_diag = 0.0; min_d = np.inf
    for n in (1, 2, 3):
        Bn = B
        for _ in range(n - 1): Bn = np.kron(Bn, B)
        occs = [o for o in itertools.product(range(n + 1), repeat=9) if sum(o) == n]
        vecs = {}
        for o in occs:
            word = [m for m in range(9) for _ in range(o[m])]
            v = np.zeros(9 ** n)
            for p in set(itertools.permutations(word)):
                w = U[:, p[0]]
                for k in p[1:]: w = np.kron(w, U[:, k])
                v += w
            vecs[o] = v / np.linalg.norm(v)
        for l in range(n + 1):
            sec = [o for o in occs if o[1] + o[7] + o[8] == l]       # g, s6, s7 carry Alice's z
            M = np.array([[vecs[o1] @ Bn @ vecs[o2] for o2 in sec] for o1 in sec])
            worst_off = max(worst_off, np.abs(M - np.diag(np.diag(M))).max())
            worst_diag = max(worst_diag, max(abs(M[i, i] - d(o[0], o[1], t)) for i, o in enumerate(sec)))
            min_d = min(min_d, np.diag(M).min())
    print(f"  t={t}: max off-diagonal {worst_off:.1e}, max |diagonal - d(a,b)| {worst_diag:.1e}, min diagonal {min_d:+.1e}")

print("(b) pair-vectorised form vs partial traces (n = 3, random C)")
rng = np.random.default_rng(0); t = 0.5
B = np.eye(9) - 3 * t * np.outer(phi, phi)
C = rng.standard_normal((27, 27)) + 1j * rng.standard_normal((27, 27))
vecC = np.transpose(C.reshape([3] * 6), (0, 3, 1, 4, 2, 5)).reshape(-1)
print(f"  partial traces {q_direct(C, 3, t):.10f}   <<C|B^(x)3|C>> {np.vdot(vecC, np.kron(np.kron(B, B), B) @ vecC).real:.10f}")

def qform(n, t=0.5):
    D = len(basis(n)); Q = np.zeros((D * D, D * D))
    for e in range(D * D):
        X = np.zeros(D * D); X[e] = 1
        Q[:, e] = q_and_grad(X.reshape(D, D), n, t)[1].ravel()
    return (Q + Q.T) / 2, D, np.array([b[2] for b in basis(n)])

print("(c) range in one occupation sector E_l (z = e2): minimum eigenvalue of the form")
for n in range(1, 6):
    Q, D, occ = qform(n)
    row = []
    for l in range(n + 1):
        idx = [r * D + c for r in range(D) for c in range(D) if occ[r] == l]
        row.append(round(float(np.linalg.eigvalsh(Q[np.ix_(idx, idx)])[0]), 6))
    print(f"  n={n}: l = 0..n -> {row}")

print("(d) block patterns {(l,m)}: one block per diagonal vs a repeated diagonal")
for n in (2, 3, 4):
    Q, D, occ = qform(n)
    allb = [(l, m) for l in range(n + 1) for m in range(n + 1)]
    good = bad = np.inf
    for k in (2, 3):
        for blocks in itertools.combinations(allb, k):
            idx = [r * D + c for r in range(D) for c in range(D) if (occ[r], occ[c]) in blocks]
            v = np.linalg.eigvalsh(Q[np.ix_(idx, idx)])[0]
            diag = [l - m for l, m in blocks]
            if len(set(diag)) == len(diag): good = min(good, v)
            else: bad = min(bad, v)
    print(f"  n={n}: one block per diagonal: min eig {good:+.1e};  repeated diagonal: min eig {bad:+.4f}")
