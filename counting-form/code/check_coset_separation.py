"""Checks for Section 4 of counting-form/README.md (floating point; the identities are exact algebra).

Claim. Let a in F_3^N and alpha != beta. If u is supported on {x : a.x = alpha} and v on {x : a.x = beta}
(any real or complex entries), then the cross term <vec(uv^T), H_N vec(vu^T)> vanishes and

    2^N q_N(uv^T - vu^T) = 2 sum_S 2^{|S|} ||P_-^{(S)} (u (x) v)||^2 >= 2 ||u||^2 ||v||^2,

where P_-^{(S)} antisymmetrises the two factors at the copies in S. Indicator vectors of two sets lying in distinct
cosets of a ternary linear code are the special case used in the counting form.

1. Random codes L (random generator matrices, dimension 0..N-1), random distinct cosets, random subsets
   (indicator vectors) and random real vectors on the cosets, N = 2..5.
2. Controls: N = 2, U = {00, 11}, V = {01, 10} (separated only by a binary parity, not by a ternary hyperplane):
   the cross term is nonzero and the identity fails.
usage: python check_coset_separation.py"""
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


def cross(A, B, N):
    """<vec A, H_N vec B> = sum_J (-1/2)^{|J|} <Tr_J A, Tr_J B>."""
    TA = A.reshape((3,) * (2 * N)); TB = B.reshape((3,) * (2 * N)); tot = 0.0
    for mask in range(2 ** N):
        K = [s for s in range(N) if mask >> s & 1]; X, Y = TA, TB
        for s in sorted(K, reverse=True):
            X = np.trace(X, axis1=s, axis2=X.ndim // 2 + s); Y = np.trace(Y, axis1=s, axis2=Y.ndim // 2 + s)
        tot += (-0.5) ** len(K) * np.sum(np.conj(X) * Y).real
    return tot


def antisym_sum(u, v, N):
    W = np.einsum('i,j->ij', u, v).reshape((3,) * (2 * N)); tot = 0.0
    for mask in range(2 ** N):
        X = W
        for s in range(N):
            if mask >> s & 1:
                X = 0.5 * (X - np.swapaxes(X, s, N + s))
        tot += 2 ** bin(mask).count('1') * np.sum(np.abs(X) ** 2)
    return tot


def index(x):
    return int(sum(int(c) * 3 ** (len(x) - 1 - i) for i, c in enumerate(x)))


rng = np.random.default_rng(20261112)
stats = {'cases': 0, 'cross': 0.0, 'ident': 0.0, 'ratio': np.inf}
for N in (2, 3, 4, 5):
    for trial in range(60 if N < 5 else 16):
        k = int(rng.integers(0, N))
        G = rng.integers(0, 3, size=(k, N))
        L = np.array(sorted({tuple(np.dot(c, G) % 3) for c in itertools.product(range(3), repeat=k)}) if k
                     else [(0,) * N])
        Lset = {tuple(l) for l in L}
        p = rng.integers(0, 3, N); qv = rng.integers(0, 3, N)
        if tuple((qv - p) % 3) in Lset:
            continue
        cp = [(p + l) % 3 for l in L]; cq = [(qv + l) % 3 for l in L]
        for kind in ('indicator', 'real'):
            u = np.zeros(3 ** N); v = np.zeros(3 ** N)
            for x in cp:
                if rng.random() < 0.6 or x is cp[0]:
                    u[index(x)] = 1.0 if kind == 'indicator' else rng.normal()
            for x in cq:
                if rng.random() < 0.6 or x is cq[0]:
                    v[index(x)] = 1.0 if kind == 'indicator' else rng.normal()
            A = np.outer(u, v); B = np.outer(v, u)
            lhs = 2 ** N * q(A - B, N); rhs = 2 * antisym_sum(u, v, N)
            stats['cases'] += 1
            stats['cross'] = max(stats['cross'], abs(cross(A, B, N)))
            stats['ident'] = max(stats['ident'], abs(lhs - rhs) / max(1.0, abs(lhs)))
            stats['ratio'] = min(stats['ratio'], lhs / (2 * np.dot(u, u) * np.dot(v, v)))
print(f"1. coset pairs: {stats['cases']} cases, N = 2..5")
print(f"   max |cross term| = {stats['cross']:.2e}, max rel. identity error = {stats['ident']:.2e}")
print(f"   min 2^N q / (2 |u|^2 |v|^2) = {stats['ratio']:.4f}  (claim: >= 1)")
u = np.zeros(9); v = np.zeros(9)
for x in ('00', '11'):
    u[index(x)] = 1
for x in ('01', '10'):
    v[index(x)] = 1
A = np.outer(u, v); B = np.outer(v, u)
print(f"2. control U={{00,11}}, V={{01,10}}: cross term = {cross(A, B, 2):.4f}, "
      f"4 q = {4 * q(A - B, 2):.4f}, 2 sum = {2 * antisym_sum(u, v, 2):.4f} (identity fails, as it should)")
