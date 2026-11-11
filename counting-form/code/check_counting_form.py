"""Numerical checks for counting-form/README.md (floating point; the identities are exact algebra).

1. Counting identity: for disjoint U, V in {0,1,2}^N and C = 1_U 1_V^T - 1_V 1_U^T,
       2^N q_N(C) = sum_S (-1)^{|S|} 2^{N-|S|} sum_{i,j} (N_S^{UV}(i,j) - N_S^{VU}(i,j))^2,
   N_S^{UV}(i,j) = #{s in {0,1,2}^S : (i,s) in U, (j,s) in V}; compared with q from partial traces.
2. Flag constants on a binary code: x = |+>^m, y(z) = 2^{-m/2} (-1)^{sum_{i<j} z_i z_j}:
       d(x,x) = 2^-m, d(x,y) = d(y,x) = (3/4)^m, d(y,y) = 1/2 + 2^-m   (m even), d(a,b) = q_m(a b^T).
3. |<ab, H cd>| <= sqrt(d(a,d) d(c,b)) for real vectors (random check).
4. Disjoint binary flags a = (00-11)/sqrt2, b = (01+10)/sqrt2: <ab,H ab> = <ab,H ba> = 1/2, and
   q_{n+2}((u x a)(v x b)^T - (v x b)(u x a)^T) = q_n(uv^T - vu^T)/2.
5. Window flags a+ = (00+11)/sqrt2, b: direct 1/2, crossed -1/2; two in a row: 1/4 and 1/4; the wedge score is
   multiplied by exactly 1/4.
6. Near-binary tags y1 = (0^L+1^L)/sqrt2, y2 = (w+wbar)/sqrt2, |w| = L/2 (L = 4): the Gram of y_i y_j^T under H_L,
   and the decomposition of the score of u_L ^ v_L with error at most 20 b t.
usage: python check_counting_form.py"""
import itertools
import numpy as np


def q(C, N):
    T = C.reshape((3,) * (2 * N)); tot = 0.0
    for mask in range(2 ** N):
        K = [s for s in range(N) if mask >> s & 1]; Y = T
        for s in sorted(K, reverse=True):
            Y = np.trace(Y, axis1=s, axis2=Y.ndim // 2 + s)
        tot += (-0.5) ** len(K) * np.sum(np.abs(Y) ** 2)
    return tot


def counting(U, V, N):
    tot = 0
    for mask in range(2 ** N):
        S = [s for s in range(N) if mask >> s & 1]; Sc = [s for s in range(N) if s not in S]
        a = np.transpose(U, S + Sc).reshape(3 ** len(S), 3 ** len(Sc)).astype(np.int64)
        b = np.transpose(V, S + Sc).reshape(3 ** len(S), 3 ** len(Sc)).astype(np.int64)
        M = a.T @ b                                  # N_S^{UV}(i, j)
        tot += (-1) ** len(S) * 2 ** (N - len(S)) * int(np.sum((M - M.T) ** 2))
    return tot


def H_inner(X, Y, m):
    return (q(X + Y, m) - q(X - Y, m)) / 4


def word(w, m):
    out = np.zeros(3 ** m); out[np.ravel_multi_index(tuple(w), (3,) * m)] = 1; return out


wedge = lambda a, c: np.outer(a, c) - np.outer(c, a)
rng = np.random.default_rng(0)

print("1. counting identity")
for N in (2, 3, 4):
    for _ in range(3):
        L = rng.integers(0, 3, size=(3,) * N); U, V = L == 1, L == 2
        u, v = U.reshape(-1).astype(float), V.reshape(-1).astype(float)
        print(f"   N={N} |U|={U.sum():3d} |V|={V.sum():3d}: counting {counting(U, V, N):8d}   2^N q = {2 ** N * q(wedge(u, v), N):12.4f}")
for Us, Vs in ((['02'], ['10', '20', '01', '11']), (['20', '02', '22'], ['00', '21'])):
    L = np.zeros((3, 3), int)
    for w in Us: L[tuple(int(c) for c in w)] = 1
    for w in Vs: L[tuple(int(c) for c in w)] = 2
    print(f"   U={Us} V={Vs}: counting {counting(L == 1, L == 2, 2)}")

print("2. flag constants on a binary code")
for m in (2, 4, 6):
    xb = np.full((2,) * m, 2 ** (-m / 2)); yb = np.zeros((2,) * m)
    for z in itertools.product((0, 1), repeat=m):
        yb[z] = 2 ** (-m / 2) * (-1) ** sum(z[i] * z[j] for i in range(m) for j in range(i + 1, m))
    emb = lambda vb: sum(vb[z] * word(z, m) for z in itertools.product((0, 1), repeat=m))
    x, y = emb(xb), emb(yb)
    print(f"   m={m}: d(x,x) = {q(np.outer(x, x), m):.6f} [{2 ** -m:.6f}]  d(x,y) = {q(np.outer(x, y), m):.6f}, d(y,x) = {q(np.outer(y, x), m):.6f}"
          f" [{0.75 ** m:.6f}]  d(y,y) = {q(np.outer(y, y), m):.6f} [{0.5 + 2 ** -m:.6f}]")

print("3. |<ab,H cd>| <= sqrt(d(a,d) d(c,b))")
worst = 0
for _ in range(200):
    a, b, c, d = (rng.standard_normal(9) for _ in range(4))
    worst = max(worst, abs(H_inner(np.outer(a, b), np.outer(c, d), 2)) / np.sqrt(q(np.outer(a, d), 2) * q(np.outer(c, b), 2)))
print(f"   largest ratio over 200 random real quadruples (m = 2): {worst:.4f}")

print("4. disjoint binary flags")
a = (word([0, 0], 2) - word([1, 1], 2)) / np.sqrt(2); b = (word([0, 1], 2) + word([1, 0], 2)) / np.sqrt(2)
print(f"   <ab,H ab> = {H_inner(np.outer(a, b), np.outer(a, b), 2):.6f}, <ab,H ba> = {H_inner(np.outer(a, b), np.outer(b, a), 2):.6f}")
u, v = rng.standard_normal(3), rng.standard_normal(3)
print(f"   wedge with the flags: {q(wedge(np.kron(u, a), np.kron(v, b)), 3):.6f} vs q_1/2 = {q(wedge(u, v), 1) / 2:.6f}")

print("5. window flags")
ap = (word([0, 0], 2) + word([1, 1], 2)) / np.sqrt(2)
print(f"   one flag: direct {H_inner(np.outer(ap, b), np.outer(ap, b), 2):+.6f}, crossed {H_inner(np.outer(ap, b), np.outer(b, ap), 2):+.6f}")
A2, B2 = np.kron(ap, ap), np.kron(b, b)
print(f"   two flags: direct {H_inner(np.outer(A2, B2), np.outer(A2, B2), 4):+.6f}, crossed {H_inner(np.outer(A2, B2), np.outer(B2, A2), 4):+.6f}")
print(f"   wedge with two flags: {q(wedge(np.kron(u, A2), np.kron(v, B2)), 5):.6f} vs q_1/4 = {q(wedge(u, v), 1) / 4:.6f}")

print("6. near-binary tags (L = 4)")
Lt = 4; w = [0, 0, 1, 1]; wb = [1, 1, 0, 0]
y1 = (word([0] * Lt, Lt) + word([1] * Lt, Lt)) / np.sqrt(2); y2 = (word(w, Lt) + word(wb, Lt)) / np.sqrt(2)
ys = [y1, y2]; tags = [(0, 0), (0, 1), (1, 0), (1, 1)]
G = np.array([[H_inner(np.outer(ys[i], ys[j]), np.outer(ys[k], ys[l]), Lt) for (k, l) in tags] for (i, j) in tags])
bb = 2 ** (-Lt / 2); dd = 0.5 + bb ** 2
print(f"   Gram in the order (11,12,21,22):\n{np.round(G, 6)}\n   expected b = {bb}, d = {dd}, b^2 = {bb ** 2}")
M = np.linalg.qr(rng.standard_normal((9, 2)))[0]; u2, v2 = M[:, 0], M[:, 1]
qpq, quv = q(wedge(a, b), 2), q(wedge(u2, v2), 2)
for t in (0.05, 0.2, 0.5):
    uL = np.sqrt(1 - t) * np.kron(a, y1) + np.sqrt(t) * np.kron(u2, y2); vL = np.sqrt(1 - t) * np.kron(b, y1) + np.sqrt(t) * np.kron(v2, y2)
    val = q(wedge(uL, vL), 2 + Lt); main = (1 - t) ** 2 * dd * qpq + t ** 2 * dd * quv
    print(f"   t={t}: score {val:+.6f}, main terms {main:+.6f}, |difference| {abs(val - main):.6f} <= 20 b t = {20 * bb * t:.6f}")
