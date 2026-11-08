"""Section 9, Propositions 28 and 29 and the crossed form of the merging budget (floating point checks).
(a) (9.4): ||Tr_1 C||^2 + ||Tr_2 C||^2 <= ||C||^2 + (1 - 1/r)||C||_1^2 + |Tr C|^2 / r for rank(C) <= r, random complex C,
    local dimensions 2..4, and its two equality families.
(b) (9.5) at rank two: q_2(t,C) - (1-t)(||C||^2 - t||C||_1^2) >= t(1/2 - t)(||C||_1^2 - |Tr C|^2).
(c) the crossed form of the merging budget: B >= (1-t)^(m-1)/2 and
    Delta_d - gamma Delta_d' = s1^2 D1 + s2^2 D2 + 2 s1 s2 (Re c + eta).
(d) Proposition 29: the merging budget on C = s1 e^{i phi}|a><a| + s2 u v^* with a a product vector and u, v orthogonal to a.
usage: python refined_partial_trace.py [seed]"""
import sys, itertools, numpy as np
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)

def crand(*s): return rng.standard_normal(s) + 1j * rng.standard_normal(s)

def ptraces(C, dims, K):
    n = len(dims); X, rem = C.reshape(list(dims) * 2), list(range(n))
    for s in sorted(K, reverse=True):
        i = rem.index(s); X = np.trace(X, axis1=i, axis2=len(rem) + i); rem.pop(i)
    return X

def q(C, dims, t):
    n = len(dims)
    return sum((-t) ** r * np.sum(np.abs(ptraces(C, dims, K)) ** 2) for r in range(n + 1) for K in itertools.combinations(range(n), r))

def NST(C):
    sv = np.linalg.svd(C, compute_uv=False)
    return np.sum(sv ** 2), sv.sum(), abs(np.trace(C)) ** 2

def gap(C, dims, t):
    N, S, T = NST(C); return q(C, dims, t) - (1 - t) ** (len(dims) - 1) * (N - t * S * S)

def fb(C, dims, t):
    return gap(C, dims, t) - (1 - 2 * t) / 2 * gap(C, [dims[0] * dims[1]] + list(dims[2:]), t)

# (a)
worst, cnt = -np.inf, 0
for _ in range(4000):
    d1, d2 = rng.integers(2, 5, size=2); D = d1 * d2; r = int(rng.integers(1, min(D, 6) + 1))
    C = crand(D, r) @ crand(D, r).conj().T
    N, S, T = NST(C); P = np.sum(np.abs(ptraces(C, [d1, d2], [0])) ** 2) + np.sum(np.abs(ptraces(C, [d1, d2], [1])) ** 2)
    worst = max(worst, (P - (N + (1 - 1 / r) * S * S + T / r)) / N); cnt += 1
print(f"(a) {cnt} random operators of rank r <= 6: max [LHS - RHS of (9.4)]/||C||^2 = {worst:.3e}")
eq = 0.0
for r in (1, 2, 3):
    d1, d2 = 4, 3; a = np.eye(d2)[0]; b = np.eye(d2)[1]; Pr = np.diag([1.0] * r + [0.0] * (d1 - r))
    for C in (np.kron(Pr, np.outer(a, b)), np.kron((lambda M: M @ M.conj().T)(crand(d1, r)), np.outer(a, a))):
        N, S, T = NST(C); P = np.sum(np.abs(ptraces(C, [d1, d2], [0])) ** 2) + np.sum(np.abs(ptraces(C, [d1, d2], [1])) ** 2)
        eq = max(eq, abs(P - (N + (1 - 1 / r) * S * S + T / r)) / N)
print(f"    equality families P (x) |a><b| (a orthogonal to b) and X (x) |a><a| (X >= 0), r = 1, 2, 3: max relative defect {eq:.1e}")
# (b)
worst = np.inf
for _ in range(1500):
    d1, d2 = rng.integers(2, 5, size=2); D = d1 * d2; t = rng.uniform(0, 0.5)
    C = crand(D, 2) @ crand(D, 2).conj().T; N, S, T = NST(C)
    worst = min(worst, (gap(C, [d1, d2], t) - t * (0.5 - t) * (S * S - T)) / N)
print(f"(b) 1500 random rank-two operators, random t in [0, 1/2]: min [Delta_2 - t(1/2 - t)(S^2 - T)]/||C||^2 = {worst:.3e}")
# (c)
def swap_op(dims, sites):
    D = int(np.prod(dims)); n = len(dims); idx = np.arange(D * D).reshape(list(dims) * 2); perm = list(range(2 * n))
    for s in sites: perm[s], perm[n + s] = perm[n + s], perm[s]
    F = np.zeros((D * D, D * D)); F[np.arange(D * D), np.transpose(idx, perm).reshape(-1)] = 1.0
    return F
k = np.kron; worst_id = worst_eig = 0.0
for dims in ([2, 2, 2], [3, 2, 2], [2, 3, 2], [2, 2, 3], [2, 2, 2, 2]):
    n = len(dims); D = int(np.prod(dims)); I = np.eye(D * D)
    for t in (0.2, 0.4, 0.5):
        m, g = 1 - t, (1 - 2 * t) / 2; F1, F2 = swap_op(dims, [0]), swap_op(dims, [1])
        B = (I - t * F1) @ (I - t * F2) - g * (I - t * F1 @ F2)
        for s in range(2, n): B = B @ (I - t * swap_op(dims, [s]))
        lam, eta = m ** (n - 1) / 2, t * m ** (n - 2) / 2
        worst_eig = max(worst_eig, abs(np.linalg.eigvalsh((B + B.T) / 2).min() - lam))
        for _ in range(3):
            U = np.linalg.qr(crand(D, 2))[0]; V = np.linalg.qr(crand(D, 2))[0]; s1, s2 = rng.random(2) + 0.1
            C = s1 * np.outer(U[:, 0], V[:, 0].conj()) + s2 * np.outer(U[:, 1], V[:, 1].conj())
            Dg = [np.vdot(k(U[:, i], V[:, i]), B @ k(U[:, i], V[:, i])).real - lam for i in (0, 1)]
            c = np.vdot(k(U[:, 0], V[:, 1]), B @ k(U[:, 1], V[:, 0]))
            worst_id = max(worst_id, abs(fb(C, dims, t) - (s1 ** 2 * Dg[0] + s2 ** 2 * Dg[1] + 2 * s1 * s2 * (c.real + eta))))
print(f"(c) smallest eigenvalue of B minus (1-t)^(m-1)/2: max |.| = {worst_eig:.1e};  crossed identity: max error {worst_id:.1e}")
# (d)
worst = np.inf; cnt = 0
while cnt < 300:
    n = int(rng.integers(2, 5)); dims = [int(x) for x in rng.integers(2, 4, size=n)]; D = int(np.prod(dims))
    if D > 100: continue
    a = np.array([1.0 + 0j])
    for d in dims:
        v = crand(d); a = np.kron(a, v / np.linalg.norm(v))
    P = np.eye(D) - np.outer(a, a.conj()); u = P @ crand(D); v = P @ crand(D)
    u /= np.linalg.norm(u); v /= np.linalg.norm(v); s1, s2 = rng.random(2) * 2
    C = s1 * np.exp(2j * np.pi * rng.random()) * np.outer(a, a.conj()) + s2 * np.outer(u, v.conj())
    worst = min(worst, fb(C, dims, rng.uniform(0, 0.5)) / np.linalg.norm(C) ** 2); cnt += 1
print(f"(d) {cnt} operators of Proposition 29 (2 to 4 sites, local dimensions 2 and 3): min (Delta_d - gamma Delta_d')/||C||^2 = {worst:.3e}")
lam1, lam2, t = 0.8, 0.5, 0.4
e0 = np.zeros((3, 3)); e0[0, 0] = 1; e1 = np.zeros((3, 3)); e1[1, 1] = 1
C = k(lam1 * k(e0, e0) - lam2 * k(e1, e1), e0)
print(f"    sharpness: Delta_d / Delta_d' on (l1|00><00| - l2|11><11|) (x) |0><0| at t = 0.4: "
      f"{gap(C, [3, 3, 3], t) / gap(C, [9, 3], t):.12f}  (gamma = {(1 - 2 * t) / 2})")
