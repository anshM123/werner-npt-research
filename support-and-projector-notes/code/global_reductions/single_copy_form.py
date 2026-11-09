"""Single-copy form of the stabilizer-domination sums (Section 7 of NOTES.md).

In the coordinates sigma = a + b, delta = b - a of two qutrits (a bijection of F_3^2, so a permutation of the
computational basis), the swap F becomes I_sigma (x) P_delta with P|x> = |-x>.  For a graph state g (amplitudes
3^{-N/2} w^{Q(x)}, Q(x) = sum_{i<j} Gamma_ij x_i x_j) the parallelogram law Q(s - d) + Q(s + d) = 2Q(s) + 2Q(d) makes
g (x) g a product g-bar (x) g-bar in these coordinates, and Z^x (x) Z^y becomes Z^{-(x+y)} (x) Z^{x-y}.  The result is

    S(d, v) = 3^N <Z^v g| B_d |Z^{-v} g>,    B_d = (x)_s (I - t Z^{d_s} P_s)  (positive definite, >= (1-t)^N I),

a statement about one copy.  This script checks:
  (1) the identity against the defining sum S(d, v) = sum_a w^{v.a} (1-3t)^{N - w_d(a)} on random graphs;
  (1') more generally <Z^x g|B_d|Z^y g> = S(d - x - y, y - x)/3^N for all x, y;
  (2) that the three vectors Z^{jv} g have the PSD Gram matrix G_jl = S(d - (j+l)v, (l-j)v) / 3^N under B_d;
  (3) that the inequality |<Z^v h|B_0|Z^{-v} h>| <= tau <h|B_0|h> fails for general states h (so SD_t uses the
      stabilizer structure), while it holds for Pauli images h = Z^d g of graph states.
usage: python single_copy_form.py N t trials seed"""
import sys, itertools
import numpy as np

N, t, trials, seed = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(seed)
lam = 1 - 3 * t; tau = t / (1 - t); om = np.exp(2j * np.pi / 3)
X = np.array(list(itertools.product(range(3), repeat=N)), dtype=np.int64)
D = 3 ** N
Z1 = np.diag([1, om, om ** 2]).astype(complex)
P1 = np.zeros((3, 3), complex); P1[0, 0] = 1; P1[1, 2] = P1[2, 1] = 1


def kron_list(ms):
    out = np.array([[1.0 + 0j]])
    for m in ms:
        out = np.kron(out, m)
    return out


def Zv(v):
    return kron_list([np.linalg.matrix_power(Z1, int(x) % 3) for x in v])


def Bd(d):
    return kron_list([np.eye(3) - t * np.linalg.matrix_power(Z1, int(x) % 3) @ P1 for x in d])


def random_graph():
    M = np.triu(rng.integers(0, 3, size=(N, N)), 1)
    return (M + M.T) % 3


def graph_state(G):
    ph = np.zeros(D, dtype=np.int64)
    for i in range(N):
        for j in range(i + 1, N):
            ph += G[i, j] * X[:, i] * X[:, j]
    return om ** (ph % 3) / np.sqrt(D)


def S_direct(G, d, v):
    w = ((X != 0) | ((X @ G.T + d) % 3 != 0)).sum(1)
    return np.sum(om ** ((X @ v) % 3) * lam ** (N - w))


err1 = err2 = 0.0; mineig = np.inf
for _ in range(trials):
    G = random_graph(); g = graph_state(G)
    d = rng.integers(0, 3, N); v = rng.integers(0, 3, N)
    B = Bd(d)
    val = 3 ** N * (Zv(v) @ g).conj() @ B @ (Zv(-v) @ g)
    err1 = max(err1, abs(val - S_direct(G, d, v)))
    vecs = [Zv(j * v) @ g for j in range(3)]
    Gram = np.array([[vecs[j].conj() @ B @ vecs[l] for l in range(3)] for j in range(3)])
    pred = np.array([[S_direct(G, (d - (j + l) * v) % 3, ((l - j) * v) % 3) / 3 ** N for l in range(3)]
                     for j in range(3)])
    err2 = max(err2, np.abs(Gram - pred).max())
    mineig = min(mineig, np.linalg.eigvalsh((Gram + Gram.conj().T) / 2).min())
err4 = 0.0
for _ in range(trials):
    G = random_graph(); g = graph_state(G)
    d, x, y = (rng.integers(0, 3, N) for _ in range(3))
    val = (Zv(x) @ g).conj() @ Bd(d) @ (Zv(y) @ g)
    err4 = max(err4, abs(val - S_direct(G, (d - x - y) % 3, (y - x) % 3) / 3 ** N))
print(f"N={N} t={t}: (1) max |3^N <Z^v g|B_d|Z^-v g> - S(d,v)| = {err1:.2e}")
print(f"           (1') max |<Z^x g|B_d|Z^y g> - S(d-x-y, y-x)/3^N| over random x, y = {err4:.2e}")
print(f"           (2) max |Gram - S-array| = {err2:.2e}, least Gram eigenvalue = {mineig:.4f} (PSD)")

B0 = Bd(np.zeros(N, dtype=np.int64)); Pn = kron_list([P1] * N)
nonzero_v = [np.array(v) for v in itertools.product(range(3), repeat=N) if any(v)]


def worst_ratio(h):
    b = (h.conj() @ B0 @ h).real
    return max(abs((Zv(v) @ h).conj() @ B0 @ (Zv(-v) @ h)) / b for v in nonzero_v)


for name in ['Pauli images Z^d g', 'random P-symmetric states', 'random states']:
    worst = 0.0
    for _ in range(trials):
        if name.startswith('Pauli'):
            h = Zv(rng.integers(0, 3, N)) @ graph_state(random_graph())
        else:
            h = rng.normal(size=D) + 1j * rng.normal(size=D)
            if 'symmetric' in name:
                h = h + Pn @ h
            h /= np.linalg.norm(h)
        worst = max(worst, worst_ratio(h))
    print(f"           (3) {name:26s}: max |<Z^v h|B_0|Z^-v h>| / <h|B_0|h> = {worst:.6f}   (tau = {tau:.6f})")
