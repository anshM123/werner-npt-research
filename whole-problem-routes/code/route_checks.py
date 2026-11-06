"""Checks of the identities behind the whole-problem routes (README of this folder), written from the definitions.

q_{n,t}(C) = sum_J (-t)^{|J|} ||Tr_J C||_F^2 and b_{n,t}(X, Y) = sum_J (-t)^{|J|} <Tr_J X, Tr_J Y>, t = 1/2 unless said.
Floating point (agreement to 1e-10) unless a section says exact.  Prints PASS at the end if every check passed."""
import itertools
from fractions import Fraction as Fr
import numpy as np

rng = np.random.default_rng(2026)
FAIL = []


def ok(name, cond, detail=''):
    print(('ok   ' if cond else 'FAIL ') + name + ('  ' + detail if detail else ''))
    if not cond:
        FAIL.append(name)


def ptrace(X, n, J, d=3):
    """Partial trace of the operator X on (C^d)^{n} over the copies in J."""
    T = X.reshape([d] * (2 * n))
    keep = [i for i in range(n) if i not in J]
    letters = 'abcdefghijklmnopqrstuvwxyz'
    rows = [letters[i] for i in range(n)]
    cols = [letters[i] if i in J else letters[n + i] for i in range(n)]
    out = ''.join(rows[i] for i in keep) + ''.join(cols[i] for i in keep)
    Y = np.einsum(''.join(rows) + ''.join(cols) + '->' + out, T)
    m = d ** len(keep)
    return Y.reshape(m, m)


def subsets(n):
    for r in range(n + 1):
        for J in itertools.combinations(range(n), r):
            yield J


def bform(X, Y, n, t=0.5, d=3):
    return sum((-t) ** len(J) * np.vdot(ptrace(X, n, J, d), ptrace(Y, n, J, d)) for J in subsets(n))


def qform(C, n, t=0.5, d=3):
    return bform(C, C, n, t, d).real


def ptranspose(X, n, K, d=3):
    T = X.reshape([d] * (2 * n))
    axes = list(range(2 * n))
    for k in K:
        axes[k], axes[n + k] = axes[n + k], axes[k]
    return T.transpose(axes).reshape(d ** n, d ** n)


def Lnorm2(X, n):
    """<X, (Id - T/2)^{tensor n}(X)> for a real matrix X on n copies (T = transpose on one copy)."""
    return sum((-0.5) ** len(K) * np.sum(X * ptranspose(X, n, K)) for K in subsets(n))


def rand_on_pair(D):
    A = rng.normal(size=(D, 2))
    Q, _ = np.linalg.qr(A)
    return Q[:, 0], Q[:, 1]


def close(a, b, tol=1e-10):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


# ---------------------------------------------------------------- 1. ratio, real formulation and fusion
print('== 1. the cross ratio and fusion')
for (n1, n2) in [(1, 1), (2, 1), (2, 2)]:
    D1, D2 = 3 ** n1, 3 ** n2
    u1, v1 = rand_on_pair(D1)
    u2, v2 = rand_on_pair(D2)
    vals = []
    for (u, v, n) in [(u1, v1, n1), (u2, v2, n2)]:
        z = (u + 1j * v) / np.sqrt(2)
        Cz, Czb = np.outer(z, z.conj()), np.outer(z.conj(), z)
        q, b = qform(Cz, n), bform(Cz, Czb, n)
        P = np.outer(u, u) + np.outer(v, v)
        H1, H2, J = np.outer(u, u) - np.outer(v, v), np.outer(u, v) + np.outer(v, u), np.outer(u, v) - np.outer(v, u)
        Delta, eps = Lnorm2(H1, n) + Lnorm2(H2, n), Lnorm2(P, n) - Lnorm2(J, n)
        ok(f'real formulation n={n}: Delta = 4 q(C_z), eps = 4 b(C_z, C_zbar), q(P) = (Delta + eps)/2',
           close(Delta, 4 * q) and abs(b.imag) < 1e-12 and close(eps, 4 * b.real) and close(qform(P, n), (Delta + eps) / 2))
        ok(f'rank-one lower bound n={n}: q(C_z) >= 2^-n', q >= 2.0 ** -n - 1e-12, f'q = {q:.6f}')
        vals.append((Delta, eps))
    us = (np.kron(u1, u2) - np.kron(v1, v2)) / np.sqrt(2)
    vs = (np.kron(u1, v2) + np.kron(v1, u2)) / np.sqrt(2)
    n = n1 + n2
    z = (us + 1j * vs) / np.sqrt(2)
    Cz, Czb = np.outer(z, z.conj()), np.outer(z.conj(), z)
    Ds, es = 4 * qform(Cz, n), 4 * bform(Cz, Czb, n).real
    ok(f'fusion ({n1},{n2}): orthonormal, Delta* = Delta1 Delta2 / 4, eps* = eps1 eps2 / 4',
       close(us @ vs, 0) and close(us @ us, 1) and close(vs @ vs, 1)
       and close(Ds, vals[0][0] * vals[1][0] / 4) and close(es, vals[0][1] * vals[1][1] / 4))
    Ps = np.outer(us, us) + np.outer(vs, vs)
    P1 = np.outer(u1, u1) + np.outer(v1, v1); P2 = np.outer(u2, u2) + np.outer(v2, v2)
    J1 = np.outer(u1, v1) - np.outer(v1, u1); J2 = np.outer(u2, v2) - np.outer(v2, u2)
    ok(f'fusion ({n1},{n2}): P* = (P1 x P2 - J1 x J2)/2', np.allclose(Ps, (np.kron(P1, P2) - np.kron(J1, J2)) / 2))

# balanced three-copy code and its fusion with itself
e = np.eye(3)
g = sum(np.kron(np.kron(e[a], e[a]), e[a]) for a in range(3)) / np.sqrt(3)
h = (np.kron(np.kron(e[0], e[1]), e[2]) + np.kron(np.kron(e[1], e[2]), e[0]) + np.kron(np.kron(e[2], e[0]), e[1])) / np.sqrt(3)
gm = []
for i in range(3):
    for j in range(3):
        E = np.zeros((3, 3)); E[i, j] = 1
        if i == j:
            E = E - np.eye(3) / 3
        gm.append(E)


def balanced(u, v, n):
    Pi = np.outer(u, u) + np.outer(v, v)
    worst = 0.0
    for s in range(n):
        for E in gm:
            Es = np.kron(np.kron(np.eye(3 ** s), E), np.eye(3 ** (n - s - 1)))
            worst = max(worst, np.abs(Pi @ Es @ Pi).max())
    return worst


ok('three-copy code span{sum_a |aaa>, |012>+|120>+|201>} is balanced', balanced(g, h, 3) < 1e-12)
ok('its q_3(P)/2 = 3/4', close(qform(np.outer(g, g) + np.outer(h, h), 3) / 2, 0.75))
us = (np.kron(g, g) - np.kron(h, h)) / np.sqrt(2)
vs = (np.kron(g, h) + np.kron(h, g)) / np.sqrt(2)
ok('its fusion with itself (six copies) is balanced', balanced(us, vs, 6) < 1e-12)

# ---------------------------------------------------------------- 2. monotone minima and GHZ padding (exact)
print('== 2. padding')
D = np.zeros((3, 3)); D[0, 1] = 1
ok('q_1(|0><1|) = 1', close(qform(D, 1), 1.0))
C = rng.normal(size=(9, 2)) @ rng.normal(size=(2, 9)) + 1j * rng.normal(size=(9, 2)) @ rng.normal(size=(2, 9))
ok('q_3(C x |0><1|) = q_2(C) for a random rank-two C', close(qform(np.kron(C, D), 3), qform(C, 2)))
def kpow(v, r):
    out = np.ones(1)
    for _ in range(r):
        out = np.kron(out, v)
    return out


for r in range(2, 6):
    gr = sum(kpow(e[a], r) for a in range(3)) / np.sqrt(3)
    Dr = np.outer(gr, gr)
    for t in (0.5, 0.6):
        val = qform(Dr, r, t)
        ok(f'h_{r}({t}) = ((1-t)^r + 2(1 + (-t)^r))/3', close(val, ((1 - t) ** r + 2 * (1 + (-t) ** r)) / 3))
half = Fr(1, 2)
for r in range(2, 12):
    hr = ((1 - half) ** r + 2 * (1 + (-half) ** r)) / 3
    want = Fr(2, 3) + Fr(1, 2 ** r) if r % 2 == 0 else Fr(2, 3) - Fr(1, 3 * 2 ** r)
    ok(f'exact h_{r}(1/2) = {hr}', hr == want and hr >= half)
P3 = np.outer(g, g) + np.outer(h, h)
g2 = sum(np.kron(e[a], e[a]) for a in range(3)) / np.sqrt(3)
ok('q_5(P_code x D_2) = q_3(P_code) h_2(1/2) = (3/2)(11/12)', close(qform(np.kron(P3, np.outer(g2, g2)), 5), 1.5 * 11 / 12))

# ---------------------------------------------------------------- 3. the real skew lift
print('== 3. real skew lift')
for (n, r) in [(1, 1), (1, 2), (2, 1), (1, 3)]:
    Dn = 3 ** n
    U, _ = np.linalg.qr(rng.normal(size=(Dn, 2))); V, _ = np.linalg.qr(rng.normal(size=(Dn, 2)))
    s1, s2 = rng.uniform(0.3, 2, 2)
    u1, u2, v1, v2 = U[:, 0], U[:, 1], V[:, 0], V[:, 1]
    C = s1 * np.outer(u1, v1) + s2 * np.outer(u2, v2)
    s = s1 + s2
    f0 = np.zeros(3 ** r); f0[0] = 1
    f1 = np.zeros(3 ** r); f1[sum(3 ** i for i in range(r))] = 1
    a = (np.sqrt(s1) * np.kron(u1, f0) + np.sqrt(s2) * np.kron(v2, f1)) / np.sqrt(s)
    b = (np.sqrt(s2) * np.kron(u2, f0) - np.sqrt(s1) * np.kron(v1, f1)) / np.sqrt(s)
    A = np.outer(a, b) - np.outer(b, a)
    JU, JV = np.outer(u1, u2) - np.outer(u2, u1), np.outer(v1, v2) - np.outer(v2, v1)
    rhs = (2 * qform(C, n) + s1 * s2 * 2.0 ** -r * (qform(JU, n) + qform(JV, n) + 2 * (-1) ** r * bform(JU, JV, n).real)) / s ** 2
    ok(f'skew lift n={n}, r={r}: a, b orthonormal and q_(n+r)(A_r) equals the formula',
       close(a @ b, 0) and close(a @ a, 1) and close(b @ b, 1) and close(qform(A, n + r), rhs))

# ---------------------------------------------------------------- 4. harmonic formulation: symbol and pair bound
print('== 4. harmonic formulation')
for N in (2, 3):
    x, y = rng.normal(size=3 ** N), rng.normal(size=3 ** N)
    lhs = 0.0
    for S in subsets(N):
        # x_{a i}: a indexes the copies in S, i the others
        perm = list(S) + [i for i in range(N) if i not in S]
        X = x.reshape([3] * N).transpose(perm).reshape(3 ** len(S), -1)
        Y = y.reshape([3] * N).transpose(perm).reshape(3 ** len(S), -1)
        M = X.T @ Y - Y.T @ X
        lhs += (-0.5) ** len(S) * np.sum(np.triu(M, 1) ** 2)
    ok(f'symbol identity N={N}: sum_S (-1/2)^|S| sum_(i<j) [...]^2 = q_N(xy^T - yx^T)/2',
       close(lhs, qform(np.outer(x, y) - np.outer(y, x), N) / 2))
for N in range(2, 12):
    for m in (1, 2, 10, 26, 27, 40):
        lower = Fr(m, 2 ** N) * ((26 - m) * 5 ** (N - 2) + m - 2)
        direct = m * (Fr(5, 2) ** N - Fr(1, 2) ** N) + m * (m - 1) * (Fr(1, 2 ** N) - Fr(1, 4) * Fr(5, 2) ** (N - 2))
        if lower != direct:
            ok(f'pair-bound algebra N={N} m={m}', False)
ok('pair-bound algebra: m[(5/2)^N - 2^-N] + m(m-1)[2^-N - (5/2)^(N-2)/4] = (m/2^N)[(26-m)5^(N-2) + m - 2] (exact, N <= 11)', True)
ok('the bound is positive for m <= 26 at every N = 2..40, and for each m = 27..60 it is negative at some N <= 40 (exact)',
   all(Fr(m, 2 ** N) * ((26 - m) * 5 ** (N - 2) + m - 2) > 0 for N in range(2, 41) for m in range(1, 27))
   and all(any(Fr(m, 2 ** N) * ((26 - m) * 5 ** (N - 2) + m - 2) < 0 for N in range(2, 41)) for m in range(27, 61)))

# ---------------------------------------------------------------- 5. the map route
print('== 5. the candidate map P_a(X) = a Tr(X) I - (X + X^T)/4')
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        F[3 * i + j, 3 * j + i] = 1
Om = sum(np.kron(e[i], e[i]) for i in range(3))
for a in (0.5, 0.75, 0.9, 1.0):
    Ca = a * np.eye(9) - (F + np.outer(Om, Om)) / 4
    ev = np.sort(np.linalg.eigvalsh(Ca))
    want = np.sort([a - 1] + [a - 0.25] * 5 + [a + 0.25] * 3)
    choi = sum(np.kron(np.outer(e[i], e[j]), a * np.trace(np.outer(e[i], e[j])) * np.eye(3)
                       - (np.outer(e[i], e[j]) + np.outer(e[j], e[i])) / 4) for i in range(3) for j in range(3))
    ok(f'a={a}: Choi matrix a I - (F + Omega Omega^T)/4 with eigenvalues a-1, a-1/4 (x5), a+1/4 (x3); invariant under partial transpose',
       np.allclose(choi, Ca) and np.allclose(ev, want)
       and np.allclose(Ca.reshape(3, 3, 3, 3).transpose(0, 3, 2, 1).reshape(9, 9), Ca))
psi = (np.kron(e[1], e[1]) + np.kron(e[2], e[2])) / np.sqrt(2)
ok('two-positivity threshold: <psi|(F + Omega Omega^T)|psi> = 3 at psi = (|11> + |22>)/sqrt2', close(psi @ (F + np.outer(Om, Om)) @ psi, 3.0))
best = 0
for it in range(300):
    A = rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2)); B = rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2))
    w = (A @ B.T).reshape(9); w = w / np.linalg.norm(w)
    best = max(best, (w.conj() @ (F + np.outer(Om, Om)) @ w).real)
ok('random Schmidt-rank-two vectors never exceed 3', best <= 3 + 1e-12, f'max {best:.6f}')
# Kraus form of R(X) = (Tr X I - X^T)/4
X = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
Ks = [(np.outer(e[i], e[j]) - np.outer(e[j], e[i])) / 2 for i in range(3) for j in range(i + 1, 3)]
ok('R(X) = (Tr X I - X^T)/4 has Kraus operators (|i><j| - |j><i|)/2', np.allclose(sum(K @ X @ K.conj().T for K in Ks), (np.trace(X) * np.eye(3) - X.T) / 4))
# a = 1 measure-and-prepare formula, exact through the moments of the uniform measure on S^2
ok_mp = True
for trial in range(5):
    Xq = [[Fr(int(v)) for v in row] for row in rng.integers(-9, 10, size=(3, 3))]
    tr = sum(Xq[i][i] for i in range(3))
    for i in range(3):
        for j in range(3):
            # E[(I - uu^T)_{ij} u^T X u] = delta_ij Tr X / 3 - (delta_ij Tr X + X_ij + X_ji)/15
            mp = Fr(15, 4) * ((tr / 3 if i == j else 0) - ((tr if i == j else 0) + Xq[i][j] + Xq[j][i]) / 15)
            direct = (tr if i == j else 0) - (Xq[i][j] + Xq[j][i]) / 4
            ok_mp &= (mp == direct)
ok('a = 1: P_1(X) = (15/4) E_u[(I - uu^T) u^T X u] (exact, five integer matrices)', ok_mp)
# twirl and the conditional state
alpha = beta = Fr(1, 4)
ok('U(x)U twirl of Omega Omega^T is (I + F)/4: 9 alpha + 3 beta = 3 and 3 alpha + 9 beta = 3 (exact)', 9 * alpha + 3 * beta == 3 and 3 * alpha + 9 * beta == 3)
a = Fr(9, 10)
ok('twirl of C_(9/10) is (67 I - 25 F)/80 (exact)', a - Fr(1, 16) == Fr(67, 80) and Fr(5, 16) == Fr(25, 80))
rho = (67 * np.eye(9) - 25 * F) / 528
ev = np.sort(np.linalg.eigvalsh(rho)); evG = np.sort(np.linalg.eigvalsh(rho.reshape(3, 3, 3, 3).transpose(0, 3, 2, 1).reshape(9, 9)))
ok('rho = (67I - 25F)/528: trace 1, eigenvalues 7/88 (x6), 23/132 (x3); partial transpose -1/66 (x1), 67/528 (x8)',
   close(np.trace(rho), 1) and np.allclose(ev, sorted([7 / 88] * 6 + [23 / 132] * 3)) and np.allclose(evG, sorted([-1 / 66] + [67 / 528] * 8)))


# tensor powers of P_a as matrices on vec(X) (row-major vec), for the recursion and the flag identity
def Pa_matrix(a):
    M = np.zeros((9, 9))
    for i in range(3):
        for j in range(3):
            Eij = np.outer(e[i], e[j])
            M[:, 3 * i + j] = (a * np.trace(Eij) * np.eye(3) - (Eij + Eij.T) / 4).reshape(9)
    return M


def Phi(a, n, X):
    """P_a^{tensor n}(X) for X on n copies."""
    T = X.reshape([3] * (2 * n))
    M = Pa_matrix(a).reshape(3, 3, 3, 3)
    for c in range(n):
        T = np.moveaxis(np.tensordot(M, T, axes=([2, 3], [c, n + c])), [0, 1], [c, n + c])
    return T.reshape(3 ** n, 3 ** n)


for a in (0.6, 0.75, 0.9, 1.0):
    m_ = a - 0.5
    for n in (1, 2):
        x, y = rng.normal(size=3 ** (n + 1)), rng.normal(size=3 ** (n + 1))
        xs, ys = x.reshape(3, -1), y.reshape(3, -1)
        Bn = lambda u, v: v @ Phi(a, n, np.outer(u, u)) @ v
        rhs = m_ * sum(Bn(xs[j], ys[i]) for i in range(3) for j in range(3))
        for i in range(3):
            for j in range(i + 1, 3):
                H = Bn(xs[j], ys[i]) + Bn(xs[i], ys[j]) - 2 * ys[i] @ Phi(a, n, np.outer(xs[j], xs[i])) @ ys[j]
                rhs += H / 2
        ok(f'recursion a={a}, n={n}: B_(n+1)(x, y) = m sum B_n(x_j, y_i) + (1/2) sum_(i<j) H_ij', close(y @ Phi(a, n + 1, np.outer(x, x)) @ y, rhs))
    for n in (1, 2):
        for k in (1, 2):
            x, y = rng.normal(size=3 ** n), rng.normal(size=3 ** n)
            f0 = np.zeros(3 ** k); f0[0] = 1
            f1 = np.zeros(3 ** k); f1[sum(3 ** i for i in range(k))] = 1
            z = np.kron(x, f0) + np.kron(y, f1)
            s = lambda w, nn: w @ Phi(a, nn, np.outer(w, w)) @ w
            bn = y @ Phi(a, n, np.outer(x, x)) @ y
            kn = x @ Phi(a, n, np.outer(x, y)) @ y
            rhs = m_ ** k * (s(x, n) + s(y, n)) + 2 * a ** k * bn + 4 * (-0.25) ** k * kn
            ok(f'flag identity a={a}, n={n}, k={k}', close(s(z, n + k), rhs))
x = np.kron(e[0], e[0]) + np.kron(e[1], e[1]); y = np.kron(e[0], e[1]) - np.kron(e[1], e[0])
z = np.kron(x, e[0]) + np.kron(y, e[1])
for a in (0.6, 0.75, 0.9):
    ok(f'control a={a}: s_3(z) = ((4a - 1)^3 - 4)/4', close(z @ Phi(a, 3, np.outer(z, z)) @ z, ((4 * a - 1) ** 3 - 4) / 4))
ok('control value at a = 3/5 is -157/500 (exact)', ((4 * Fr(3, 5) - 1) ** 3 - 4) / 4 == Fr(-157, 500))
# product minima m^n at equal real unit vectors
for a in (0.75, 0.9):
    for n in (1, 2, 3):
        w = np.zeros(3 ** n); w[0] = 1
        ok(f'product value a={a}, n={n}: <e|P_a^n(ee^T)|e> = (a - 1/2)^n', close(w @ Phi(a, n, np.outer(w, w)) @ w, (a - 0.5) ** n))

print()
print('PASS' if not FAIL else 'FAIL: ' + '; '.join(FAIL))
