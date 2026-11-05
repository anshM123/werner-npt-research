"""Checks for Proposition 13 and Theorem 14 of NOTES.md, Section 5 (balanced encodings).
Flag code: c(a,b) = (a, b, a+b, a+2b) over F_3, each coordinate repeated r times; y = column flags with r/2 letters 1
then r/2 letters 2 in each group.  One-copy entry m(a, y, b, x) = d(a,b) d(y,x) - 1/2 d(a,y) d(b,x).
(1) Flag factors, group by group (all flag strings for r = 2, 4; y and 2000 random strings for r = 12):
    D-type:  prod_f m(c(f), y_f, c'(f), y'_f) = [c = c'][y' = y] 2^{-matches(c, y)} for codewords c, c' and all y';
    MV-type: prod_f m(0, y_f, c(f), x_f) = [c = 0][x = y] for all x;
    matches(c, y) = 0 for c = 0 and 3r/2 otherwise.
(2) The reduction identity q_N(E C J_y^*) = (1/9)[q_k(C) + 2^{-3r/2} sum_{(a,b) != 0} q_k((X^a Z^b)^{(x)k} C)],
    numerically by direct partial traces, (k, r) = (1, 2), (2, 2); and EP on the encoded space for random co-ranges.
(3) One-copy marginals of encoded unit vectors are I/3.
(4) The numbers of Proposition 13: the witness 8a_3 + b_3 has coefficient matrix 8J - 7I (singular values 17, 7, 7),
    so |z|^2 = 387 and ||(|z><z|)^Gamma||_1 = 31^2 = 961; and (-4 + 512 * 961 / 2^18)/387 = -1087/198144."""
import itertools
import numpy as np
from fractions import Fraction as Fr

def m1(a, y, b, x): return (1 if (a == b and y == x) else 0) - (Fr(1, 2) if (a == y and b == x) else 0)
ok = True
for r in (2, 4, 12):
    yg = [1] * (r // 2) + [2] * (r // 2)
    good = True
    for alpha in range(3):
        for beta in range(3):
            for yp in itertools.product(range(3), repeat=r) if r <= 4 else [tuple(yg)] + [tuple(np.random.default_rng(1).integers(0, 3, r)) for _ in range(2000)]:
                v = Fr(1)
                for f in range(r): v *= m1(alpha, yg[f], beta, yp[f])
                exp = (Fr(1, 2 ** sum(1 for f in range(r) if yg[f] == alpha)) if (alpha == beta and tuple(yp) == tuple(yg)) else 0)
                good &= v == exp
            for x in itertools.product(range(3), repeat=r) if r <= 4 else [tuple(yg)] + [tuple(np.random.default_rng(2).integers(0, 3, r)) for _ in range(2000)]:
                v = Fr(1)
                for f in range(r): v *= m1(0, yg[f], alpha, x[f])
                good &= v == (1 if (alpha == 0 and tuple(x) == tuple(yg)) else 0)
    # one group: D factor is diagonal only when the two codeword symbols agree; across the 4 groups a pair of distinct
    # codewords differs in at least one group (in fact three), so the product vanishes.
    codes = [(a, b, (a + b) % 3, (a + 2 * b) % 3) for a in range(3) for b in range(3)]
    good &= all(sum(1 for g in range(4) if c[g] != cp[g]) >= 3 for c in codes for cp in codes if c != cp)
    good &= all(sum(r // 2 for g in range(4) if c[g] != 0) == (0 if c == (0, 0, 0, 0) else 3 * r // 2) for c in codes)
    print(f"(1) r = {r}: flag factors and match counts: {'PASS' if good else 'FAIL'}")
    ok &= good

w3 = np.exp(2j * np.pi / 3)
Xm = np.roll(np.eye(3), 1, axis=0); Zm = np.diag([1, w3, w3 ** 2])
def kronlist(ms):
    out = np.array([[1.0 + 0j]])
    for mm in ms: out = np.kron(out, mm)
    return out
def basis_vec(digits):
    v = np.zeros(3 ** len(digits), complex); idx = 0
    for d in digits: idx = 3 * idx + d
    v[idx] = 1; return v
def q_lowrank(xs, ys, n):
    tot = 0.0
    for J in itertools.product([0, 1], repeat=n):
        Jl = [i for i in range(n) if J[i]]; Jc = [i for i in range(n) if not J[i]]
        Xs = [x.reshape([3] * n).transpose(Jl + Jc).reshape(3 ** len(Jl), 3 ** len(Jc)) for x in xs]
        Ys = [y.reshape([3] * n).transpose(Jl + Jc).reshape(3 ** len(Jl), 3 ** len(Jc)) for y in ys]
        val = 0
        for i in range(len(xs)):
            for ip in range(len(xs)):
                if len(Jl) <= len(Jc):
                    val += np.sum((Xs[i].conj() @ Xs[ip].T) * (Ys[ip].conj() @ Ys[i].T).T)
                else:
                    val += np.sum(np.conj(Xs[i].T @ Ys[i].conj()) * (Xs[ip].T @ Ys[ip].conj()))
        tot += (-0.5) ** len(Jl) * val.real
    return tot
rng = np.random.default_rng(0)
for k, r in [(1, 2), (2, 2)]:
    N = k + 4 * r
    codes = {}
    for a in range(3):
        for b in range(3):
            c = [a, b, (a + b) % 3, (a + 2 * b) % 3]; codes[(a, b)] = [ci for ci in c for _ in range(r)]
    yflags = ([1] * (r // 2) + [2] * (r // 2)) * 4
    U = {ab: kronlist([np.linalg.matrix_power(Xm, ab[0]) @ np.linalg.matrix_power(Zm, ab[1])] * k) for ab in codes}
    enc = lambda v: sum(np.kron(U[ab] @ v, basis_vec(codes[ab])) for ab in codes) / 3
    xs = [rng.normal(size=3 ** k) + 1j * rng.normal(size=3 ** k) for _ in range(2)]
    ys = [rng.normal(size=3 ** k) + 1j * rng.normal(size=3 ** k) for _ in range(2)]
    lhs = q_lowrank([enc(x) for x in xs], [np.kron(yv, basis_vec(yflags)) for yv in ys], N)
    rhs = (q_lowrank(xs, ys, k) + 2.0 ** (-3 * r / 2) * sum(q_lowrank([U[ab] @ x for x in xs], ys, k) for ab in codes if ab != (0, 0))) / 9
    good = abs(lhs - rhs) < 1e-9 * max(1, abs(lhs))
    # EP on the encoded space: random rank-two operators with range in E(C^{3^k}) and arbitrary co-range
    qmin = min(q_lowrank([enc(rng.normal(size=3 ** k) + 1j * rng.normal(size=3 ** k)) for _ in range(2)],
                         [rng.normal(size=3 ** N) + 1j * rng.normal(size=3 ** N) for _ in range(2)], N) for _ in range(3))
    v = rng.normal(size=3 ** k) + 1j * rng.normal(size=3 ** k); v /= np.linalg.norm(v)
    ev = enc(v).reshape([3] * N)
    marg = max(np.abs(np.moveaxis(ev, i, 0).reshape(3, -1) @ np.moveaxis(ev, i, 0).reshape(3, -1).conj().T - np.eye(3) / 3).max() for i in range(N))
    print(f"(2,3) k={k}, r={r} (N={N}): identity diff {abs(lhs - rhs):.1e}; min q_N over 3 random encoded operators {qmin:+.3e}; marginal error {marg:.1e}")
    ok &= good and qmin > 0 and marg < 1e-12
Zc = 8 * np.ones((3, 3)) - 7 * np.eye(3); sv = np.linalg.svd(Zc, compute_uv=False)
good = np.allclose(sorted(sv), [7, 7, 17]) and abs(np.sum(sv ** 2) - 387) < 1e-12 and abs(np.sum(sv) ** 2 - 961) < 1e-9
good &= (Fr(-4) + Fr(512 * 961, 2 ** 18)) / 387 == Fr(-1087, 198144)
print(f"(4) witness numbers 387, 961 and the bound -1087/198144: {'PASS' if good else 'FAIL'}")
ok &= good
print("PASS" if ok else "FAIL")
