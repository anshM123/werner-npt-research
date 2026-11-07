"""Numerics for Conjecture QT_m and Remark 7.21 (tensor products of the qutrit maps P_a with different parameters).

python multiparam_numerics.py real n seed
    For listed and random parameter vectors a = (a_1, ..., a_n), estimate the least value of
    <y| P_{a_1} (x) ... (x) P_{a_n} (x x^T) |y> over real unit x, y (alternating minimisation: y is the bottom eigenvector
    for fixed x and then x for fixed y; 10 to 12 starts, three of them near plane GHZ vectors on random subsets of the
    copies) and compare its sign with the rule prod_{s in S} (2 a_s - 1/2) >= 1/2 for every nonempty S.
python multiparam_numerics.py two
    n = 2: 3000 random parameter pairs, least eigenvalue over the reduced inputs sum_i s_i e_i (x) e_i (every real input
    up to rotations), against the same rule.
python multiparam_numerics.py complex
    (a) the explicit zero at n = 7 with every Re a_s = 0.76 (direct evaluation with 3^7-dimensional vectors);
    (a') the explicit zero at n = 5 with Re a_1 = 0.53 and Re a_2..a_5 = 0.76 (inside the region of Conjecture QT_m);
    (b) adversarial search at n = 3: the largest Re a_1 of a zero with Re a_2, Re a_3 > 3/4 (expected 1/2).
Floating point throughout; prints MISMATCH lines if the rule and the sign disagree."""
import itertools, sys
import numpy as np
from scipy.optimize import minimize

d = 3

def apply_map(a_list, X):
    n = len(a_list)
    T = X.reshape([d] * (2 * n)).astype(complex if np.iscomplexobj(np.array(a_list)) else float)
    for s, a in enumerate(a_list):
        T = np.moveaxis(T, [s, n + s], [0, 1]); sh = T.shape
        Mi = T.reshape(3, 3, -1); out = np.zeros_like(Mi); tr = np.einsum('iir->r', Mi)
        for k in range(3): out[k, k] += a * tr
        out -= 0.25 * (Mi + Mi.transpose(1, 0, 2))
        T = out.reshape(sh); T = np.moveaxis(T, [0, 1], [s, n + s])
    return T.reshape(d**n, d**n)

def rule(a_list):
    vals = [2 * a - 0.5 for a in a_list]; n = len(a_list)
    return min(np.prod([vals[s] for s in S]) for r in range(1, n + 1) for S in itertools.combinations(range(n), r))

def minval(a_list, rng, starts=12, iters=60):
    n = len(a_list); best = np.inf
    p = np.array([1, 1j, 0]) / np.sqrt(2)
    for st in range(starts):
        if st < 3:
            S = [s for s in range(n) if rng.random() < 0.7] or [0]
            g1 = np.ones(1, dtype=complex); g2 = np.ones(1, dtype=complex)
            for s in range(n):
                if s in S: g1 = np.kron(g1, p); g2 = np.kron(g2, p.conj())
                else: g1 = np.kron(g1, np.array([0, 0, 1.0])); g2 = np.kron(g2, np.array([0, 0, 1.0]))
            x = ((g1 + g2) / np.sqrt(2)).real + 0.05 * rng.standard_normal(d**n)
        else:
            x = rng.standard_normal(d**n)
        x /= np.linalg.norm(x)
        for it in range(iters):
            W = apply_map(a_list, np.outer(x, x)); w, V = np.linalg.eigh((W + W.T) / 2); y = V[:, 0]
            W2 = apply_map(a_list, np.outer(y, y)); w2, V2 = np.linalg.eigh((W2 + W2.T) / 2); x = V2[:, 0]
        best = min(best, w2[0])
    return best

mode = sys.argv[1]
if mode == 'real':
    n = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3]))
    if n == 3:
        tests = [(0.5, 0.75, 0.75), (0.5, 0.76, 0.76), (0.52, 0.74, 0.76), (0.6, 0.6, 1.0), (0.62, 0.62, 0.65),
                 (0.64, 0.64, 0.7), (0.646850, 0.646850, 0.646850), (0.6475, 0.6475, 0.6475), (0.55, 0.7, 0.72),
                 (0.55, 0.72, 0.72), (0.58, 0.66, 0.7), (0.6, 0.65, 0.66), (0.6, 0.66, 0.665)]
        lo, hi, extra = 0.5, 0.8, 8
    else:
        tests = [(0.5, 0.75, 0.75, 0.75), (0.52, 0.75, 0.76, 0.77), (0.6036, 0.6036, 0.75, 0.75), (0.60, 0.6036, 0.8, 0.8),
                 (0.6704, 0.6704, 0.6704, 0.6704), (0.6710, 0.6710, 0.6710, 0.6710), (0.62, 0.66, 0.7, 0.75),
                 (0.62, 0.66, 0.72, 0.75), (0.55, 0.72, 0.72, 0.74), (0.55, 0.73, 0.735, 0.74), (0.64, 0.65, 0.66, 0.9)]
        lo, hi, extra = 0.55, 0.8, 6
    for _ in range(extra): tests.append(tuple(float(z) for z in np.round(rng.uniform(lo, hi, n), 4)))
    mism = 0
    for a in tests:
        mv = minval(list(a), rng, starts=12 if n == 3 else 10, iters=60 if n == 3 else 50); g = rule(list(a))
        pred = 'pos' if g >= 0.5 - 1e-12 else 'NEG'; obs = 'pos' if mv >= -1e-9 else 'NEG'
        mism += pred != obs
        print(f"a = {a}: min over S of prod (2a_s - 1/2) = {g:.5f} (predict {pred}); least value found {mv:+.3e} ({obs})"
              f"{'' if pred == obs else '  MISMATCH'}", flush=True)
    print(f"n = {n}: {len(tests)} parameter vectors, {mism} mismatches")
elif mode == 'two':
    rng = np.random.default_rng(3)
    def P(a, X): return a * np.trace(X) * np.eye(3) - 0.25 * (X + X.T)
    def lam(a1, a2, s):
        x = np.zeros(9)
        for i in range(3): x[4 * i] = s[i]
        X = np.outer(x, x).reshape(3, 3, 3, 3); out = np.zeros((9, 9))
        for i1, j1, i2, j2 in itertools.product(range(3), repeat=4):
            c = X[i1, i2, j1, j2]
            if c == 0: continue
            E = np.zeros((3, 3)); E[i1, j1] = 1; F = np.zeros((3, 3)); F[i2, j2] = 1
            out += c * np.kron(P(a1, E), P(a2, F))
        return np.linalg.eigvalsh((out + out.T) / 2)[0]
    bad = 0; worst_in = np.inf; worst_out = -np.inf
    for _ in range(3000):
        a1, a2 = rng.uniform(0.45, 1.2, 2)
        ok = a1 >= 0.5 and a2 >= 0.5 and (2 * a1 - 0.5) * (2 * a2 - 0.5) >= 0.5
        cands = [rng.standard_normal(3) for _ in range(30)] + [np.array([1, 1, 0.0]), np.array([1, 0, 0.0]), np.array([1, 1, 1.0])]
        m = min(lam(a1, a2, c / np.linalg.norm(c)) for c in cands)
        if ok: worst_in = min(worst_in, m)
        else: worst_out = max(worst_out, m)
        bad += (m >= -1e-12) != ok
    print(f"n = 2: 3000 parameter pairs, {bad} mismatches; least value inside the region {worst_in:.2e}, "
          f"largest least value outside {worst_out:.2e}")
elif mode == 'complex':
    from scipy.optimize import brentq
    n = 7; eps = 0.01
    phase = lambda eta: np.angle(0.75 + eps + 1j * eta) - np.angle(0.25 + eps + 1j * eta)
    etas = np.linspace(0.01, 2, 2000); k = int(np.argmin([phase(e) for e in etas]))
    eta = brentq(lambda e: n * phase(e) + np.pi, 0.01, etas[k])
    a = 0.75 + eps + 1j * eta
    p = abs((a - 0.5)**n) / (abs(a**n) + abs((a - 0.5)**n))
    e1 = np.array([1.0, 0, 0]); e2 = np.array([0, 1.0, 0]); u = np.ones(1); v = np.ones(1)
    for _ in range(n): u = np.kron(u, e2); v = np.kron(v, e1)
    x = np.sqrt(p) * u + np.sqrt(1 - p) * v
    val = v @ apply_map([a] * n, np.outer(x, x)) @ v
    print(f"(a) n = 7, a_s = {a:.6f} for every s (Re a_s = {a.real:.2f}), p = {p:.6f}: "
          f"<y| (x) P_a (x x^T) |y> = {val:.2e}; formula p prod a + (1-p) prod (a - 1/2) = {p * a**n + (1 - p) * (a - 0.5)**n:.2e}")
    from scipy.optimize import fsolve
    n5 = 5; aa = 0.76 + 0.4445j
    def G(w):
        a1_ = 0.53 + 1j * w[0]; f = w[1] * a1_ * aa**4 + (1 - w[1]) * (a1_ - 0.5) * (aa - 0.5)**4
        return [f.real, f.imag]
    w = fsolve(G, [0.0995, 0.022]); a1 = 0.53 + 1j * w[0]; p5 = w[1]
    u5 = np.ones(1); v5 = np.ones(1)
    for _ in range(n5): u5 = np.kron(u5, e2); v5 = np.kron(v5, e1)
    x5 = np.sqrt(p5) * u5 + np.sqrt(1 - p5) * v5
    val5 = v5 @ apply_map([a1] + [aa] * 4, np.outer(x5, x5)) @ v5
    rule5 = min(np.prod([2 * z.real - 0.5 for z, keep in zip([a1] + [aa] * 4, S) if keep])
                for S in itertools.product([0, 1], repeat=5) if any(S))
    print(f"(a') n = 5, a_1 = {a1:.6f}, a_2..a_5 = {aa:.4f}, p = {p5:.6f}: <y| (x) P_a (x x^T) |y> = {val5:.2e}; "
          f"real parts satisfy the rule of QT_m (least subset product {rule5:.3f} >= 1/2)")
    # (b) adversarial, n = 3
    n = 3; rng = np.random.default_rng(2)
    subsets = [U for r in range(n + 1) for U in itertools.combinations(range(n), r)]
    def coeffs(x, y):
        X = x.reshape([d] * n); Y = y.reshape([d] * n); c = []
        for U in subsets:
            r = len(U); rest = [s for s in range(n) if s not in U]
            Xm = np.moveaxis(X, list(U) + rest, list(range(n))).reshape(d**r, -1)
            Ym = np.moveaxis(Y, list(U) + rest, list(range(n))).reshape(d**r, -1)
            rx = Xm @ Xm.T; ry = Ym @ Ym.T
            if r == 0: c.append(float(rx[0, 0] * ry[0, 0])); continue
            T = rx.reshape([d] * (2 * r)); acc = np.zeros_like(rx)
            for kk in range(r + 1):
                for W in itertools.combinations(range(r), kk):
                    perm = list(range(2 * r))
                    for w in W: perm[w], perm[r + w] = perm[r + w], perm[w]
                    acc = acc + np.transpose(T, perm).reshape(d**r, d**r)
            c.append((-0.5)**r * float(np.sum(ry * acc)) / 2**r)
        return np.array(c)
    def fe(c, a):
        tot = 0
        for cu, U in zip(c, subsets):
            pr = 1
            for s in range(n):
                if s not in U: pr *= a[s]
            tot += cu * pr
        return tot
    def obj(v):
        D = d**n; x = v[:D] / np.linalg.norm(v[:D]); y = v[D:2 * D] / np.linalg.norm(v[D:2 * D]); ar = v[2 * D:]
        rest = [0.75 + np.log1p(np.exp(ar[2 * j])) + 1j * ar[2 * j + 1] for j in range(n - 1)]
        c = coeffs(x, y); A0 = fe(c, [0] + rest); B = fe(c, [1] + rest) - A0
        r = -A0 / B
        return 10.0 if not np.isfinite(r) else -r.real
    best = -np.inf
    for st in range(6):
        D = d**n
        v0 = np.concatenate([rng.standard_normal(2 * D), rng.normal(size=2 * (n - 1)) * np.array([1, 3] * (n - 1))])
        if st % 2: v0[D:2 * D] = v0[:D] + 0.1 * rng.standard_normal(D)
        r = minimize(obj, v0, method='L-BFGS-B', options={'maxiter': 3000, 'maxfun': 200000})
        best = max(best, -r.fun)
        print(f"(b) n = 3, start {st}: largest Re a_1 of a zero with Re a_2, Re a_3 > 3/4: {-r.fun:.6f}", flush=True)
    print(f"(b) n = 3: largest Re a_1 found = {best:.6f}")
