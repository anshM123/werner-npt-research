"""The completely positive half of P_{5/6} (whole-problem-routes/README.md, Section 5, remark after Proposition M).
F = (5D - id - 2T)/12 and T.F = (5D - 2id - T)/12, so F + T.F = P_{5/6} = (5/6)D - (id + T)/4.
 (a) Choi spectra of F and T.F; Choi(F) = (1/4)(I - P_Phi) + (1/3)P_anti, so F breaks Schmidt number three;
     the interval Phi_mu = P_0 - P_2/4 + mu P_1 (spin sectors of M_3 under SO(3)) is CP iff 1/12 <= mu <= 7/12.
 (b) one-copy minimum output eigenvalue 1/6 for F and T.F.
 (c) the split F = D/6 + G does not work termwise: G (x) G is not positive (least value -1/72).
 (d) alternating minimisation of <phi| T_S F^(x)n (psi psi^*) |phi> over unit psi, phi (random and structured starts):
     the least value found equals 6^-n.
 (e) (MP2): two-positivity of T_S F^(x)m (a qubit reference attached), m <= m2_max, every |S|; three-positivity fails
     already at m = 2.
usage: python completion_five_sixths.py [n_max (default 5)] [random starts (default 8)] [m2_max (default 3)]"""
import sys, itertools
import numpy as np

rng = np.random.default_rng(56)
I3 = np.eye(3)


def lin(a, b, c):
    """X -> a Tr(X) I + b X + c X^T on the last two axes."""
    def f(Tp):
        tr = np.trace(Tp, axis1=-2, axis2=-1)
        return b * Tp + c * np.swapaxes(Tp, -1, -2) + a * tr[..., None, None] * I3
    return f


def apply(X, n, maps):
    T = X.reshape([3] * (2 * n))
    for s in range(n):
        Tp = np.moveaxis(T, [s, n + s], [2 * n - 2, 2 * n - 1])
        T = np.moveaxis(maps[s](Tp), [2 * n - 2, 2 * n - 1], [s, n + s])
    return T.reshape(3 ** n, 3 ** n)


def choi(f):
    out = np.zeros((9, 9), dtype=complex)
    for a in range(3):
        for b in range(3):
            E = np.zeros((1, 1, 3, 3), dtype=complex); E[0, 0, a, b] = 1
            out += np.kron(np.outer(I3[a], I3[b]), f(E).reshape(3, 3))
    return out


F = lin(5 / 12, -1 / 12, -2 / 12); TF = lin(5 / 12, -2 / 12, -1 / 12); G = lin(3 / 12, -1 / 12, -2 / 12)
Fsw = np.zeros((9, 9)); Om = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        Fsw[3 * i + j, 3 * j + i] = 1; Om[4 * i, 4 * j] = 1
Panti = (np.eye(9) - Fsw) / 2; PPhi = Om / 3
ok = []
ev = np.round(12 * np.linalg.eigvalsh(choi(F)), 9); evt = np.round(12 * np.linalg.eigvalsh(choi(TF)), 9)
ok.append(np.allclose(choi(F), (np.eye(9) - PPhi) / 4 + Panti / 3))
print(f"(a) 12 x Choi spectrum of F: {sorted(set(ev))}, of T.F: {sorted(set(evt))};"
      f" Choi(F) = (I - P_Phi)/4 + P_anti/3: {ok[-1]}")
cp_ok = True
for mu in np.linspace(-0.2, 0.7, 91):
    Phi = lin(5 / 12, mu / 2 - 1 / 8, -(mu / 2 + 1 / 8))
    cp = np.linalg.eigvalsh(choi(Phi)).min() > -1e-12
    cp_ok &= cp == (1 / 12 - 1e-12 <= mu <= 7 / 12 + 1e-12)
ok.append(cp_ok)
print(f"(a) Phi_mu = (5/12)D + (mu/2 - 1/8)id - (mu/2 + 1/8)T is CP exactly for 1/12 <= mu <= 7/12: {cp_ok}"
      f" (F = Phi_(1/12), T.F = Phi_(-1/12))")


def least(n, maps, starts, structured=True):
    D = 3 ** n; best = np.inf
    seeds = []
    if structured:
        g = np.zeros(D)
        for i in range(3): g[sum(i * 3 ** s for s in range(n))] = 1
        seeds.append(g)
        p = np.array([1, 1j, 0]) / np.sqrt(2); m = p.conj(); a_, b_ = np.ones(1, complex), np.ones(1, complex)
        for _ in range(n): a_, b_ = np.kron(a_, p), np.kron(b_, m)
        seeds.append(a_ + b_)
        sing = np.kron(I3[0], I3[1]) - np.kron(I3[1], I3[0]); v = np.ones(1)
        for _ in range(n // 2): v = np.kron(v, sing)
        if n % 2: v = np.kron(v, I3[2])
        seeds.append(v)
        w = np.zeros(D)
        for s in range(n): w[3 ** s] = 1
        seeds.append(w)
    for _ in range(starts):
        seeds.append(rng.normal(size=D) + 1j * rng.normal(size=D))
    for psi in seeds:
        psi = np.asarray(psi, complex) + 1e-3 * (rng.normal(size=D) + 1j * rng.normal(size=D)); psi /= np.linalg.norm(psi)
        old = np.inf
        for _ in range(400):
            Y = apply(np.outer(psi, psi.conj()), n, maps); e1, U1 = np.linalg.eigh((Y + Y.conj().T) / 2); phi = U1[:, 0]
            Z = apply(np.outer(phi, phi.conj()), n, maps); e2, U2 = np.linalg.eigh((Z + Z.conj().T) / 2); psi = U2[:, 0]
            if abs(e2[0] - old) < 1e-15: break
            old = e2[0]
        best = min(best, e2[0])
    return best


b1 = [least(1, [F], 20, False), least(1, [TF], 20, False)]
ok.append(np.allclose(b1, [1 / 6, 1 / 6], atol=1e-10))
print(f"(b) one copy: least output eigenvalue of F {b1[0]:.12f}, of T.F {b1[1]:.12f} (1/6 = {1/6:.12f})")
g2 = least(2, [G, G], 20, True)
ok.append(abs(g2 + 1 / 72) < 1e-9)
print(f"(c) G = F - D/6 = (3D - id - 2T)/12: least output eigenvalue of G (x) G = {g2:.12f} (-1/72 = {-1/72:.12f})")
nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5; starts = int(sys.argv[2]) if len(sys.argv) > 2 else 8
for n in range(2, nmax + 1):
    for k in range(1, n // 2 + 1):
        val = least(n, [TF] * k + [F] * (n - k), starts)
        ok.append(abs(val * 6 ** n - 1) < 1e-8)
        print(f"(d) n = {n}, |S| = {k}: least value found {val:.6e}, times 6^n = {val * 6 ** n:.9f}", flush=True)


def apply_ref(X, m, maps, r):
    T = X.reshape([r] + [3] * m + [r] + [3] * m)
    for s_ in range(m):
        a_, b_ = 1 + s_, 2 + m + s_
        Tp = np.moveaxis(T, [a_, b_], [-2, -1]); T = np.moveaxis(maps[s_](Tp), [-2, -1], [a_, b_])
    return T.reshape(r * 3 ** m, r * 3 ** m)


def least_ref(m, maps, r, starts):
    D = r * 3 ** m; best = np.inf
    for _ in range(starts):
        psi = rng.normal(size=D) + 1j * rng.normal(size=D); psi /= np.linalg.norm(psi); old = np.inf
        for _ in range(500):
            Y = apply_ref(np.outer(psi, psi.conj()), m, maps, r); e1, U1 = np.linalg.eigh((Y + Y.conj().T) / 2); phi = U1[:, 0]
            Z = apply_ref(np.outer(phi, phi.conj()), m, maps, r); e2, U2 = np.linalg.eigh((Z + Z.conj().T) / 2); psi = U2[:, 0]
            if abs(e2[0] - old) < 1e-15: break
            old = e2[0]
        best = min(best, e2[0])
    return best


m2max = int(sys.argv[3]) if len(sys.argv) > 3 else 3
for m in range(1, m2max + 1):
    for k in range(1, m + 1):
        val = least_ref(m, [TF] * k + [F] * (m - k), 2, 12 if m < 4 else 8)
        ok.append(val > -1e-12)
        print(f"(e) two-positivity, m = {m}, |S| = {k}: least value with a qubit reference {val:+.3e}", flush=True)
v3 = least_ref(2, [TF, F], 3, 12)
ok.append(v3 < -1e-3)
print(f"(e) three-positivity fails at m = 2: least value with a qutrit reference {v3:+.6f}")
print(f"\n{sum(ok)} of {len(ok)} checks as stated")
