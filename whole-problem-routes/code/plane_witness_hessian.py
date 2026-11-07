"""Checks for Proposition 7.22 (the plane witness is a second-order local minimum on the facet, every n).

On the facet prod_s (2 a_s - 1/2) = 1/2, a_s in [1/2, 3/4], take the vectors x0, y0 of Theorem 7.1 (on the real plane
spanned by e1, e2; y0 = x0 for odd n and (Y (x) I) x0 for even n) and f(x, y) = <y| (x)_s P_{a_s}(x x^T) |y>.
Writing M_x = (x)P(x0 x0^T), M_y = (x)P(y0 y0^T) and K_ab = <e_a| (x)P(x0 e_b^T + e_b x0^T) |y0>, the Hessian of f at
(x0, y0) is H = [[M_y, K^T], [K, M_x]].

For n = 2, ..., 6 at the centre (a_s = a_n), at the vertex (1/2, 3/4, ..., 3/4), at an edge point and at two random
facet points:
1. f(x0, y0) = 0 and M_x y0 = M_y x0 = 0 (first derivatives vanish);
2. H is positive semidefinite (least eigenvalue >= -1e-12);
3. the part of H on vectors with a component along e3 at some copy (the 'thick' part) has least eigenvalue 0 and least
   positive eigenvalue equal to the closed form min over T != {} and sigma of A_T - 4^-|T| |w_sigma| (where that is
   positive), with A_T = prod_{s in T} a_s prod_{s not in T} (a_s - 1/4) and w_sigma as in the proof;
4. away from the vertices, the kernel of H has dimension 2n + 3 (rotations at single copies, 2n + 1 of them, and the
   two scalings), at the centre and the random points.
Also: 5. the inequality prod_{s in T} 4 a_s >= 1 + 2 prod_{s in T} (2 a_s - 1/2) at 2000 random facet points and all T.
Ends with PASS or FAIL."""
import itertools
import numpy as np

d = 3
ok = True
rng = np.random.default_rng(22)


def apply_map(a_list, X):
    n = len(a_list)
    T = X.reshape([d] * (2 * n)).astype(float)
    for s, a in enumerate(a_list):
        T = np.moveaxis(T, [s, n + s], [0, 1]); sh = T.shape
        Mi = T.reshape(3, 3, -1); out = np.zeros_like(Mi); tr = np.einsum('iir->r', Mi)
        for k in range(3): out[k, k] += a * tr
        out -= 0.25 * (Mi + Mi.transpose(1, 0, 2))
        T = out.reshape(sh); T = np.moveaxis(T, [0, 1], [s, n + s])
    return T.reshape(d**n, d**n)


def witness(n):
    p = np.array([1, 1j, 0]) / np.sqrt(2)
    g1 = np.ones(1, dtype=complex); g2 = np.ones(1, dtype=complex)
    for _ in range(n): g1 = np.kron(g1, p); g2 = np.kron(g2, p.conj())
    x0 = ((g1 + g2) / np.sqrt(2)).real; x0 /= np.linalg.norm(x0)
    if n % 2 == 1:
        y0 = x0.copy()
    else:
        Ym = np.zeros((3, 3)); Ym[0, 1] = -1; Ym[1, 0] = 1
        Op = Ym
        for _ in range(n - 1): Op = np.kron(Op, np.eye(3))
        y0 = Op @ x0
    return x0, y0


def closed_form_gap(a_list):
    """least positive value of A_T - 4^-|T| |w_sigma| over T != {} and sigma."""
    n = len(a_list); u = [2 * a - 0.5 for a in a_list]; vals = []
    for t in range(1, n + 1):
        for T in itertools.combinations(range(n), t):
            Tc = [s for s in range(n) if s not in T]
            A = np.prod([a_list[s] for s in T]) * np.prod([a_list[s] - 0.25 for s in Tc])
            for K in itertools.product([0, 1], repeat=len(Tc)):      # 1 = p, 0 = conj p at each copy of T^c
                c = np.prod([(a_list[s] - 0.25) if k else -0.5 for s, k in zip(Tc, K)])
                cb = np.prod([-0.5 if k else (a_list[s] - 0.25) for s, k in zip(Tc, K)])
                w = c + cb if n % 2 == 1 else c - cb
                vals.append(A - 4.0**(-t) * abs(w))
    vals = np.array(vals)
    return vals.min(), vals[vals > 1e-12].min()


for n in range(2, 7):
    an = (1 + 2**(1 - 1 / n)) / 4
    pts = [('centre', [an] * n), ('vertex', [0.5] + [0.75] * (n - 1)),
           ('edge', [(2**-0.5 + 0.5) / 2] * 2 + [0.75] * (n - 2))]
    for k in range(2):
        while True:
            t = rng.dirichlet(np.ones(n)) * (-np.log(2))
            if all(-np.log(2) <= ts <= 0 for ts in t): break
        pts.append((f'random{k}', [(np.exp(ts) + 0.5) / 2 for ts in t]))
    x0, y0 = witness(n)
    D = d**n
    for name, a_list in pts:
        Mx = apply_map(a_list, np.outer(x0, x0)); My = apply_map(a_list, np.outer(y0, y0))
        Mx = (Mx + Mx.T) / 2; My = (My + My.T) / 2
        f0 = y0 @ Mx @ y0
        r1, r2 = np.linalg.norm(Mx @ y0), np.linalg.norm(My @ x0)
        K = np.zeros((D, D))
        for b in range(D):
            eb = np.zeros(D); eb[b] = 1
            K[:, b] = apply_map(a_list, np.outer(x0, eb) + np.outer(eb, x0)) @ y0
        H = np.block([[My, K.T], [K, Mx]]); H = (H + H.T) / 2
        w = np.linalg.eigvalsh(H)
        thick = [i for i in range(D) if 2 in np.unravel_index(i, [d] * n)]
        idx = thick + [D + i for i in thick]
        wt = np.linalg.eigvalsh(H[np.ix_(idx, idx)])
        cf_min, cf_gap = closed_form_gap(a_list)
        num_gap = wt[wt > 1e-10].min()
        kdim = int(np.sum(np.abs(w) < 1e-10))
        c1 = abs(f0) < 1e-12 and r1 < 1e-12 and r2 < 1e-12
        c2 = w[0] > -1e-12
        c3 = abs(wt[0]) < 1e-10 and abs(cf_min) < 1e-12 and abs(num_gap - cf_gap) < 1e-9
        c4 = (kdim == 2 * n + 3) if name in ('centre', 'random0', 'random1') else True
        good = c1 and c2 and c3 and c4
        ok &= good
        print(f"n = {n} {name:8s}: f0 {f0:+.0e}, first order {max(r1, r2):.0e}; least eig H {w[0]:+.1e}; thick part: least "
              f"{wt[0]:+.1e}, least positive {num_gap:.6e} (closed form {cf_gap:.6e}); kernel dim {kdim}"
              f"{'' if good else '  FAIL'}", flush=True)
# 5. the inequality
worst2 = np.inf; dev1 = 0.0
for _ in range(2000):
    n = int(rng.integers(2, 9))
    while True:
        t = rng.dirichlet(np.ones(n)) * (-np.log(2))
        if all(-np.log(2) <= ts <= 0 for ts in t): break
    a = [(np.exp(ts) + 0.5) / 2 for ts in t]
    for k in range(1, n + 1):
        for T in itertools.combinations(range(n), k):
            lhs = np.prod([4 * a[s] for s in T]); rhs = 1 + 2 * np.prod([2 * a[s] - 0.5 for s in T])
            if k == 1: dev1 = max(dev1, abs(lhs - rhs))
            else: worst2 = min(worst2, lhs - rhs)
i_ok = dev1 < 1e-12 and worst2 >= 1 - 1e-12
print(f"5. prod_T 4a_s - 1 - 2 prod_T (2a_s - 1/2) over 2000 random facet points (n <= 8): |T| = 1 identity (max "
      f"deviation {dev1:.1e}); |T| >= 2 least value {worst2:.4f} (>= 1): {i_ok}")
ok &= i_ok
print("PASS" if ok else "FAIL")
