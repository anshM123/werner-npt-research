"""Checks for the reduction below the endpoint (README, "Below the endpoint").  The form at a fixed t is
b_n(t; X, Y) = <X, M_t^{(x)n}(Y)>, M_t(Z) = Z - t Tr(Z) I_3, q_n(t, X) = b_n(t; X, X), m = 1 - t.
(1) exact (sympy): the b_2(t)-Gram of the flag operators of (f, g) and of (h_0, h_1) is
    [[a,b,0,0],[b,a,0,0],[0,0,c,0],[0,0,0,c]] with a = 1 - 2t/3 + t^2, b = t^2 - 2t/3, c = 1 - 2t/3, as polynomials in t;
    and the real orthogonal design (1/12) sum O^T H O = Tr(H) I/3.
(2) the dilation identity q_{n+k}(t,P) = 2q_n(t,C) + m^k[q_n(t,D_0) + q_n(t,D_1)] + 2(-t)^k Re b_n(t;D_0,D_1),
    P/lambda a projector, and q_{n+k}(t, P/lambda) <= 2 q_n(t,C)/lambda^2 + 4 m^k.
(3) the realification identity q_{n+2r}(t,C_r) = 1/2[a^r q_n(t,C) + b^r Re b_n(t;C,conj C) + c^r q_n(t,T)].
(4) the balancing expansion on a small instance: one input copy, the two-element group {I, O}, one label bit,
    l = 2 repetitions (output on 5 copies): the exact four-index expansion, and
    |q - a^2 q_1/2| <= (a^2/2)(2^3 - 1) rho^2 |C|^2 with rho = c/a.
(5) the merging telescope on C = a1|0^N><0^N| + a2|1^N><1^N|:
    sum_k m^(k-1) X_k = 2 a1 a2 [m^(N-1) - (-t)^(N-1)], X_k the two-site term of Lemma 27 at the k-th merge.
usage: python below_endpoint_check.py"""
import itertools
import numpy as np
import sympy as sp

rng = np.random.default_rng(11)
om = np.exp(2j * np.pi / 3)


def bform(X, Y, n, t):
    T = Y.reshape((3,) * (2 * n))
    for i in range(n):
        tr = np.trace(T, axis1=i, axis2=n + i)
        T = T - t * np.moveaxis(np.multiply.outer(tr, np.eye(3)), [-2, -1], [i, n + i])
    return np.vdot(X.reshape(-1), T.reshape(-1))


ok = True

# (1) exact flag Gram, as polynomials in t
ts = sp.symbols('t')
w = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2          # exp(2 pi i / 3)


def svec(entries):
    v = sp.zeros(9, 1)
    for (i, j), val in entries.items():
        v[3 * i + j] = val
    return v / sp.sqrt(3)


def b2(X, Y):
    """b_2(t; X, Y) for 9x9 sympy matrices, by explicit partial traces."""
    def red(A, site):
        B = sp.zeros(9, 9)
        for a0, a1, b0, b1 in itertools.product(range(3), repeat=4):
            if site == 0 and a0 == b0:
                B[3 * a0 + a1, 3 * b0 + b1] = sum(A[3 * k + a1, 3 * k + b1] for k in range(3))
            if site == 1 and a1 == b1:
                B[3 * a0 + a1, 3 * b0 + b1] = sum(A[3 * a0 + k, 3 * b0 + k] for k in range(3))
        return B
    R = Y - ts * red(Y, 0)
    R = R - ts * red(R, 1)
    return sp.expand(sum(sp.conjugate(X[i, j]) * R[i, j] for i in range(9) for j in range(9)))


a_, b_, c_ = 1 - 2 * ts / 3 + ts ** 2, ts ** 2 - 2 * ts / 3, 1 - 2 * ts / 3
want = sp.Matrix([[a_, b_, 0, 0], [b_, a_, 0, 0], [0, 0, c_, 0], [0, 0, 0, c_]])
fS = svec({(k, k): w ** k for k in range(3)})
gS = fS.conjugate()
h0 = svec({(k, k): 1 for k in range(3)})
h1 = svec({(k, (k + 1) % 3): 1 for k in range(3)})
for name, ops in {'(f, g)': [fS * fS.H, gS * gS.H, fS * gS.H, gS * fS.H],
                  '(h0, h1)': [h0 * h0.T, h1 * h1.T, h0 * h1.T, h1 * h0.T]}.items():
    G = sp.Matrix(4, 4, lambda i, j: b2(ops[i], ops[j]))
    good = all(sp.simplify(sp.expand(e)) == 0 for e in (G - want))
    ok &= good
    print(f"(1) flag Gram of {name} equals [[a,b,0,0],[b,a,0,0],[0,0,c,0],[0,0,0,c]] as polynomials in t: {good}")
S = sp.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
Hs = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'H{i}{j}'))
avg = sp.zeros(3, 3)
for j in range(3):
    for e in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]:
        O = S ** j * sp.diag(*e)
        avg += O.T * Hs * O
good = sp.simplify(avg / 12 - Hs.trace() / 3 * sp.eye(3)) == sp.zeros(3, 3)
ok &= good
print(f"(1) design: (1/12) sum O^T H O = Tr(H) I/3 for a symbolic 3x3 H: {good}")

# (2) dilation below the endpoint
worst = 0.0
worst_b = -np.inf
for t in (0.2, 0.35, 0.4, 0.5):
    m = 1 - t
    for n, k in [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2)]:
        D = 3 ** n
        U = np.linalg.qr(rng.normal(size=(D, 2)) + 1j * rng.normal(size=(D, 2)))[0]
        V = np.linalg.qr(rng.normal(size=(D, 2)) + 1j * rng.normal(size=(D, 2)))[0]
        s = np.sort(rng.uniform(0.3, 2, 2))[::-1]
        C = U @ np.diag(s) @ V.conj().T
        lam = 2 * s[0]
        al = (lam + np.sqrt(lam ** 2 - 4 * s ** 2)) / 2
        be = lam - al
        l = U * np.sqrt(al)
        r = V * np.sqrt(be)
        f0 = np.zeros(3 ** k)
        f0[0] = 1
        f1 = np.zeros(3 ** k)
        f1[int('1' * k, 3)] = 1
        wv = [np.kron(l[:, i], f0) + np.kron(r[:, i], f1) for i in range(2)]
        P = sum(np.outer(x, x.conj()) for x in wv)
        D0 = l @ l.conj().T
        D1 = r @ r.conj().T
        lhs = bform(P, P, n + k, t).real
        rhs = (2 * bform(C, C, n, t).real + m ** k * (bform(D0, D0, n, t).real + bform(D1, D1, n, t).real)
               + 2 * (-t) ** k * bform(D0, D1, n, t).real)
        Pi = P / lam
        worst = max(worst, abs(lhs - rhs), np.abs(Pi @ Pi - Pi).max())
        worst_b = max(worst_b, bform(Pi, Pi, n + k, t).real - (2 * bform(C, C, n, t).real / lam ** 2 + 4 * m ** k))
print(f"(2) dilation identity and projector property: max error {worst:.1e};  max of q(Pi) - bound: {worst_b:.3e} (must be <= 0)")
ok &= worst < 1e-9 and worst_b <= 1e-12

# (3) realification below the endpoint
f = np.zeros(9, complex)
for k in range(3):
    f[4 * k] = om ** k / np.sqrt(3)
g = f.conj()
worst = 0.0
for t in (0.2, 0.35, 0.4, 0.5):
    a, b, c = 1 - 2 * t / 3 + t * t, t * t - 2 * t / 3, 1 - 2 * t / 3
    for n, rr in [(1, 1), (1, 2), (2, 1)]:
        D = 3 ** n
        U = np.linalg.qr(rng.normal(size=(D, 2)) + 1j * rng.normal(size=(D, 2)))[0]
        V = np.linalg.qr(rng.normal(size=(D, 2)) + 1j * rng.normal(size=(D, 2)))[0]
        s = rng.uniform(0.3, 2, 2)
        C = U @ np.diag(s) @ V.conj().T
        T = U @ np.diag(s) @ V.T
        F, G = f, g
        for _ in range(rr - 1):
            F, G = np.kron(F, f), np.kron(G, g)

        def R(z, F=F, G=G):
            return (np.kron(z, F) + np.kron(z.conj(), G)) / np.sqrt(2)
        Cr = sum(s[j] * np.outer(R(U[:, j]), R(V[:, j])) for j in range(2))
        lhs = bform(Cr, Cr, n + 2 * rr, t).real
        rhs = 0.5 * (a ** rr * bform(C, C, n, t).real + b ** rr * bform(C, C.conj(), n, t).real
                     + c ** rr * bform(T, T, n, t).real)
        worst = max(worst, abs(lhs - rhs), np.abs(Cr.imag).max(), abs(np.linalg.norm(Cr) - np.linalg.norm(C)))
print(f"(3) realification identity, reality and norm: max error {worst:.1e}")
ok &= worst < 1e-9

# (4) balancing expansion, small instance
h0n = np.zeros(9)
h1n = np.zeros(9)
for k in range(3):
    h0n[4 * k] = 1 / np.sqrt(3)
    h1n[3 * k + (k + 1) % 3] = 1 / np.sqrt(3)
hj = [h0n, h1n]
Sn = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], float)
O = Sn @ np.diag([1, -1, -1.])
group = [np.eye(3), O]
H = [np.kron(h0n, h0n), np.kron(h1n, h1n)]          # label bit j -> flag h_j, repeated l = 2 times
Vb = sum(np.kron(group[j], H[j].reshape(-1, 1)) for j in range(2)) / np.sqrt(2)
worst = np.abs(Vb.T @ Vb - np.eye(3)).max()
viol = -np.inf
for t in (0.35, 0.4, 0.5):
    a, b, c = 1 - 2 * t / 3 + t * t, t * t - 2 * t / 3, 1 - 2 * t / 3
    rho = c / a
    gram = {q4: bform(np.outer(hj[q4[0]], hj[q4[1]]), np.outer(hj[q4[2]], hj[q4[3]]), 2, t).real ** 2
            for q4 in itertools.product(range(2), repeat=4)}
    for trial in range(3):
        C = rng.normal(size=(3, 2)) @ rng.normal(size=(2, 3))            # real, rank two, one copy
        Ct = Vb @ C @ Vb.T
        lhs = bform(Ct, Ct, 5, t).real
        rhs = sum(bform(group[j] @ C @ group[k_].T, group[jj] @ C @ group[kk].T, 1, t).real * gram[(j, k_, jj, kk)]
                  for (j, k_, jj, kk) in itertools.product(range(2), repeat=4)) / 4
        main = a ** 2 * bform(C, C, 1, t).real / 2
        worst = max(worst, abs(lhs - rhs))
        viol = max(viol, abs(lhs - main) - (a ** 2 / 2) * (2 ** 3 - 1) * rho ** 2 * np.linalg.norm(C) ** 2)
print(f"(4) balancing: isometry and four-index expansion, max error {worst:.1e}; "
      f"max of |q - main| - bound: {viol:.3e} (must be <= 0)")
ok &= worst < 1e-9 and viol <= 1e-12


# (5) merging telescope on the GHZ-type pair
def qd(C, dims, t):
    n = len(dims)
    T0 = C.reshape(list(dims) * 2)
    tot = 0.0
    for rr in range(n + 1):
        for K in itertools.combinations(range(n), rr):
            X, rem = T0, list(range(n))
            for s_ in sorted(K, reverse=True):
                i = rem.index(s_)
                X = np.trace(X, axis1=i, axis2=len(rem) + i)
                rem.pop(i)
            tot += (-t) ** rr * np.sum(np.abs(X) ** 2)
    return tot


worst = 0.0
for N in (3, 4, 5):
    for t in (0.3, 0.4, 0.5):
        m = 1 - t
        a1, a2 = rng.uniform(0.2, 2, 2)
        e0 = np.zeros(3 ** N)
        e0[0] = 1
        e1 = np.zeros(3 ** N)
        e1[int('1' * N, 3)] = 1
        C = a1 * np.outer(e0, e0) + a2 * np.outer(e1, e1)
        Xs = []
        for k in range(1, N):     # q_{d_k} = m q_{d_{k+1}} + t X_k with d_k = (3^k, 3, ..., 3)
            dk = [3 ** k] + [3] * (N - k)
            dk1 = [3 ** (k + 1)] + [3] * (N - k - 1)
            Xs.append((qd(C, dk, t) - m * qd(C, dk1, t)) / t)
        wsum = sum(m ** (k - 1) * Xs[k - 1] for k in range(1, N))
        worst = max(worst, abs(wsum - 2 * a1 * a2 * (m ** (N - 1) - (-t) ** (N - 1))))
t, m = 0.4, 0.6
e0 = np.zeros(27)
e0[0] = 1
e1 = np.zeros(27)
e1[13] = 1
C1 = np.outer(e0, e0) + np.outer(e1, e1)
X1 = (qd(C1, [3, 3, 3], t) - m * qd(C1, [9, 3], t)) / t
X2 = (qd(C1, [9, 3], t) - m * qd(C1, [27], t)) / t
print(f"(5) |000><000| + |111><111| at t = 0.4: first two-site term {X1:+.6f} (-2t), second {X2:+.6f} (2), "
      f"weighted sum {X1 + m * X2:+.6f} (2(1 - 2t))")
print(f"(5) telescope on a1|0^N><0^N| + a2|1^N><1^N|, N = 3, 4, 5: max |sum - 2 a1 a2 [m^(N-1) - (-t)^(N-1)]| = {worst:.1e}")
ok &= worst < 1e-9
print("PASS" if ok else "FAIL")
