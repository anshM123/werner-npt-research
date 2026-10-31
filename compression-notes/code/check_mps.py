"""Section 4: translation-invariant MPS filters and the transfer channel.

   y_a = sum_x (l_a^T A^{x_1} ... A^{x_n} r) |x>,  z_a = sum_x (m_a^T B^{x_1} ... B^{x_n} r') |x>,  C = sum_a |y_a><z_a|.
   With the one-copy channel T(X) = sum W[(x x'),(u u')] (A^u (x) B^u') X (A^x (x) B^x')^*, W = I - t F (a local
   filtering of X (x) W, completely positive since W >= 0), the script checks
       q_n(t, C) = <omega| [T^n(|r r'><r r'|)]^G |omega>,   omega = sum_a conj(l_a) (x) m_a,
   G = partial transpose on the second bond factor, against the dense q_n for D = 2, 3 and n <= 4. It then shows the
   Perron behaviour for a random channel, and re-evaluates a stored D = 2 search output at 60 digits for n = 33."""
import numpy as np
from werner import d, q

rng = np.random.default_rng(4)


def channel(A, B, t):
    """T as a D^4 x D^4 matrix acting on vec(X), X indexed by (iA iB),(jA jB)"""
    D = A.shape[1]
    W = (np.eye(d * d) - t * np.eye(d * d)[[d * (k % d) + k // d for k in range(d * d)]]).reshape(d, d, d, d)
    T = np.einsum('uvxy,uik,vjl,xmn,yop->ijmoklnp', W, A, B, A.conj(), B.conj())
    return T.reshape(D ** 4, D ** 4)


def pt2(X, D):
    return X.reshape(D, D, D, D).transpose(0, 3, 2, 1).reshape(D * D, D * D)


def mps_vec(M, lv, rv, n):
    v = lv[None, :]
    for _ in range(n):
        v = np.einsum('pi,xij->pxj', v, M).reshape(-1, M.shape[2])
    return v @ rv


for D in (2, 3):
    A = rng.standard_normal((d, D, D)) + 1j * rng.standard_normal((d, D, D))
    B = rng.standard_normal((d, D, D)) + 1j * rng.standard_normal((d, D, D))
    l = rng.standard_normal((2, D)) + 1j * rng.standard_normal((2, D))
    m = rng.standard_normal((2, D)) + 1j * rng.standard_normal((2, D))
    r = rng.standard_normal(D) + 1j * rng.standard_normal(D)
    rp = rng.standard_normal(D) + 1j * rng.standard_normal(D)
    T = channel(A, B, 0.5)
    omega = sum(np.kron(l[a].conj(), m[a]) for a in range(2))
    X = np.outer(np.kron(r, rp), np.kron(r, rp).conj()).reshape(-1)
    for n in range(1, 5):
        X = T @ X
        val = np.vdot(omega, pt2(X.reshape(D * D, D * D), D) @ omega).real
        C = sum(np.outer(mps_vec(A, l[a], r, n), mps_vec(B, m[a], rp, n).conj()) for a in range(2))
        print(f"D={D} n={n}: transfer formula {val:+.10e}   dense q_n {q(C, n, 0.5):+.10e}", flush=True)

# Perron behaviour for one random channel (D = 2)
D = 2
A = rng.standard_normal((d, D, D)) + 1j * rng.standard_normal((d, D, D))
B = rng.standard_normal((d, D, D)) + 1j * rng.standard_normal((d, D, D))
T = channel(A, B, 0.5)
ev, V = np.linalg.eig(T); k = np.argmax(abs(ev)); mu = ev[k]
Xs = V[:, k].reshape(D * D, D * D); Xs = Xs / np.trace(Xs); Xs = (Xs + Xs.conj().T) / 2
print(f"\nrandom D=2 channel: Perron eigenvalue {mu:.6f} (next |ev| {np.sort(abs(ev))[-2]:.6f});"
      f" eigenvalues of the Perron operator {np.round(np.linalg.eigvalsh(Xs), 6)};"
      f" of its partial transpose {np.round(np.linalg.eigvalsh(pt2(Xs, D)), 6)}")

# a stored output of the D = 2 search, evaluated at 60 digits
try:
    import mpmath as mp
    mp.mp.dps = 60
    th = np.load('../data/mps_search_D2.npy')
    c = [mp.mpc(float(th[2 * i]), float(th[2 * i + 1])) for i in range(len(th) // 2)]
    k = 0

    def take(cnt):
        global k
        out = c[k:k + cnt]; k += cnt; return out
    A2 = np.array(take(3 * D * D), dtype=object).reshape(3, D, D); B2 = np.array(take(3 * D * D), dtype=object).reshape(3, D, D)
    l2 = np.array(take(2 * D), dtype=object).reshape(2, D); m2 = np.array(take(2 * D), dtype=object).reshape(2, D)
    r2 = np.array(take(2 * D), dtype=object).reshape(2, D); rp2 = np.array(take(2 * D), dtype=object).reshape(2, D)
    cj = np.vectorize(lambda z: mp.conj(z), otypes=[object])

    def Tmp(tt):
        N = D ** 4; M = mp.matrix(N, N)
        for x in range(3):
            for y in range(3):
                K0 = np.kron(np.kron(np.kron(cj(A2[x]), cj(B2[y])), A2[x]), B2[y])
                KF = np.kron(np.kron(np.kron(cj(A2[x]), cj(B2[y])), A2[y]), B2[x])
                for i in range(N):
                    for j in range(N):
                        M[i, j] += K0[i, j] - tt * KF[i, j]
        return M
    k4 = lambda a1, a2, a3, a4: np.kron(np.kron(np.kron(a1, a2), a3), a4)
    Th, T0 = Tmp(mp.mpf(1) / 2), Tmp(mp.mpf(0))
    Rs = [mp.matrix(list(k4(cj(r2[a]), cj(rp2[b]), r2[b], rp2[a]))) for a in range(2) for b in range(2)]
    Ls = [list(k4(cj(l2[a]), cj(m2[b]), l2[b], m2[a])) for a in range(2) for b in range(2)]
    V1 = [R.copy() for R in Rs]; V0 = [R.copy() for R in Rs]
    for _ in range(33):
        V1 = [Th * v for v in V1]; V0 = [T0 * v for v in V0]
    qv = sum(sum(Ls[i][j] * V1[i][j] for j in range(D ** 4)) for i in range(4))
    nv = sum(sum(Ls[i][j] * V0[i][j] for j in range(D ** 4)) for i in range(4))
    print(f"stored D=2 search output, n = 33, 60 digits: 2^n q_n/||C||^2 = {mp.nstr(mp.re(qv) / mp.re(nv) * 2 ** 33, 12)}")
except ImportError:
    print("mpmath not available; skipping the 60-digit evaluation")
