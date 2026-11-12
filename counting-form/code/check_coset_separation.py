"""Checks for Sections 4 and 5 of counting-form/README.md (floating point; the identities are exact algebra).

Claim. Let a in F_3^N and alpha != beta. If u is supported on {x : a.x = alpha} and v on {x : a.x = beta}
(any real or complex entries), then the cross term <vec(uv^T), H_N vec(vu^T)> vanishes and

    2^N q_N(uv^T - vu^T) = 2 sum_S 2^{|S|} ||P_-^{(S)} (u (x) v)||^2 >= 2 ||u||^2 ||v||^2,

where P_-^{(S)} antisymmetrises the two factors at the copies in S. Indicator vectors of two sets lying in distinct
cosets of a ternary linear code are the special case used in the counting form.

1. Random codes L (random generator matrices, dimension 0..N-1), random distinct cosets, random subsets
   (indicator vectors) and random real vectors on the cosets, N = 2..5.
2. Controls: N = 2, U = {00, 11}, V = {01, 10} (separated only by a binary parity, not by a ternary hyperplane):
   the cross term is nonzero and the identity fails.
3. Level sets of random integer additive functions F(x) = sum_s f_s(x_s), N = 2..4 (product phase unitary).
4. Section 5 peel: u = w (x) ut at one copy, q_N(u^v) = 1/2 q_{N-1}(ut^vt_0) + 2 sum_{k=1,2} q_{N-1}(ut vt_k^T).
5. Balanced rectangles: the cross term equals the signed rectangle sum; rectangle-free pairs give zero; the 3-copy example.
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

# 3. Level sets of an integer additive function F(x) = sum_s f_s(x_s) (phase unitary prod_s diag(e^{i t f_s})).
stats = {'cases': 0, 'cross': 0.0, 'ident': 0.0, 'ratio': np.inf}
for N in (2, 3, 4):
    for trial in range(40):
        f = rng.integers(-2, 3, size=(N, 3))
        pts = list(itertools.product(range(3), repeat=N))
        F = {x: int(sum(f[s, x[s]] for s in range(N))) for x in pts}
        levels = sorted(set(F.values()))
        if len(levels) < 2:
            continue
        c1, c2 = rng.choice(levels, 2, replace=False)
        u = np.zeros(3 ** N); v = np.zeros(3 ** N)
        for x in pts:
            if F[x] == c1 and rng.random() < 0.7:
                u[index(''.join(map(str, x)))] = rng.normal()
            if F[x] == c2 and rng.random() < 0.7:
                v[index(''.join(map(str, x)))] = rng.normal()
        if not u.any() or not v.any():
            continue
        A = np.outer(u, v); B = np.outer(v, u)
        lhs = 2 ** N * q(A - B, N); rhs = 2 * antisym_sum(u, v, N)
        stats['cases'] += 1
        stats['cross'] = max(stats['cross'], abs(cross(A, B, N)))
        stats['ident'] = max(stats['ident'], abs(lhs - rhs) / max(1.0, abs(lhs)))
        stats['ratio'] = min(stats['ratio'], lhs / (2 * np.dot(u, u) * np.dot(v, v)))
print(f"3. integer additive level sets: {stats['cases']} cases, N = 2..4, max |cross| = {stats['cross']:.2e}, "
      f"identity error {stats['ident']:.2e}, min ratio {stats['ratio']:.4f}")

# 4. One-copy product peel (Section 5): u = w (x) ut at copy s; rotate w to e_0 at that copy:
#    q_N(u^v) = 1/2 q_{N-1}(ut ^ vt_0) + 2 sum_{k=1,2} q_{N-1}(ut vt_k^T).
worst = 0.0; count = 0
for N in (2, 3, 4):
    for trial in range(20):
        s = int(rng.integers(0, N)); w = rng.normal(size=3); w /= np.linalg.norm(w); ut = rng.normal(size=3 ** (N - 1))
        u = np.moveaxis(np.einsum('a,b->ab', w, ut).reshape((3,) * N), 0, s).reshape(-1)
        v = rng.normal(size=3 ** N)
        Q, _ = np.linalg.qr(np.column_stack([w, rng.normal(size=(3, 2))])); R = Q.T
        R[0] = w
        Vt = np.tensordot(R, np.moveaxis(v.reshape((3,) * N), s, 0), axes=(1, 0)).reshape(3, -1)
        lhs = q(np.outer(u, v) - np.outer(v, u), N)
        rhs = 0.5 * q(np.outer(ut, Vt[0]) - np.outer(Vt[0], ut), N - 1) + 2 * sum(q(np.outer(ut, Vt[k]), N - 1) for k in (1, 2))
        worst = max(worst, abs(lhs - rhs) / max(1.0, abs(lhs))); count += 1
print(f"4. one-copy product peel: {count} cases, N = 2..4, max rel. error {worst:.2e}")

# 5. Balanced rectangles: the cross term equals the signed rectangle sum; it vanishes when there is none;
#    the 3-copy example U = {200,222,102,221}, V = {021,202,120,211}.
def rect_sum(u, v, N):
    pts = list(itertools.product(range(3), repeat=N)); tot = 0.0
    U = [x for x in pts if u[index(''.join(map(str, x)))] != 0]; V = [y for y in pts if v[index(''.join(map(str, y)))] != 0]
    val = lambda w, x: w[index(''.join(map(str, x)))]
    rects = 0
    for x in U:
        for xp in U:
            for y in V:
                for yp in V:
                    c = 1.0
                    for s in range(N):
                        a, b, cc, d = x[s], xp[s], y[s], yp[s]
                        c *= (1.0 if (a == d and b == cc) else 0.0) - 0.5 * (1.0 if (a == cc and b == d) else 0.0)
                        if c == 0.0:
                            break
                    if c != 0.0:
                        tot += val(u, x) * val(u, xp) * val(v, y) * val(v, yp) * c; rects += 1
    return tot, rects
worst = 0.0; empty = 0; empty_zero = 0; count = 0
for N in (2, 3):
    pts = list(itertools.product(range(3), repeat=N))
    for trial in range(40):
        k = int(rng.integers(1, 5)); sel = rng.choice(len(pts), 2 * k, replace=False)
        u = np.zeros(3 ** N); v = np.zeros(3 ** N)
        for i in sel[:k]:
            u[index(''.join(map(str, pts[i])))] = rng.normal()
        for i in sel[k:]:
            v[index(''.join(map(str, pts[i])))] = rng.normal()
        rs, nrect = rect_sum(u, v, N); cr = cross(np.outer(u, v), np.outer(v, u), N); count += 1
        worst = max(worst, abs(rs - cr) / max(1.0, abs(cr)))
        if nrect == 0:
            empty += 1; empty_zero += abs(cr) < 1e-12
U5 = ['200', '222', '102', '221']; V5 = ['021', '202', '120', '211']
u = np.zeros(27); v = np.zeros(27)
for x in U5:
    u[index(x)] = 1
for y in V5:
    v[index(y)] = 1
rs, nrect = rect_sum(u, v, 3)
print(f"5. rectangle formula vs cross term: {count} random pairs (N = 2, 3), max rel. difference {worst:.1e}; "
      f"{empty} pairs without rectangles, cross term zero in {empty_zero}")
print(f"   example: balanced rectangles {nrect}, cross term {cross(np.outer(u, v), np.outer(v, u), 3):.3f}, "
      f"q_3 = {q(np.outer(u, v) - np.outer(v, u), 3):.4f} (33/2 = 16.5)")
