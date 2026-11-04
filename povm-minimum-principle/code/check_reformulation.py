"""Exact facts of the POVM picture, checked numerically at n = 2, 3, 4.
(a) m = lambda_min K agrees with the minimum over co-ranges of q_n(C)/||C||^2 computed from q directly
    (single-excitation and tangent ranges).
(b) m depends only on the POVM: unchanged under U_s -> W_s U_s (W_s unitary on E_s) and under simultaneous
    conjugation A_s -> V^* A_s V.
(c) The sector form: K preserves the excitation number c of co-range strings, and 2^(n-c) K_c equals
       Kt_c = I + sum_{s vacuum} Abar_s - 1/2 sum_{s excited} |Ut_s><Ut_s| + 1/2 sum_{t excited, s vacuum} Z_s^* Z_t;
    sector 0 is the constant 2^(1-n); sector n is >= 0; and Kt_c >= 2I - 3 sum_{t excited} Abar_t (weight bound).
(d) At extreme POVMs (every element of rank one) m >= 0."""
import numpy as np, itertools
from povm_form import strings, psi_n, lam_map, q, sqrt_psd, tangent_range, K_of, m_povm, povm_from

rng = np.random.default_rng(0); ok = True
def m_direct(u, n):
    D = 3 ** n; basis = []
    for r in range(2):
        for j in range(D):
            V = np.zeros((2, D), complex); V[r, j] = 1; basis.append(np.outer(u[r], V[r].conj()))
    B = np.array([b.ravel() for b in basis]); L = np.array([lam_map(b, n).ravel() for b in basis])
    Qm = B.conj() @ L.T; Qm = (Qm + Qm.conj().T) / 2; Nm = B.conj() @ B.T; Nm = (Nm + Nm.conj().T) / 2
    w, P = np.linalg.eigh(Nm); Pk = P[:, w > 1e-10 * w.max()] / np.sqrt(w[w > 1e-10 * w.max()])
    return np.linalg.eigvalsh(Pk.conj().T @ Qm @ Pk)[0]
for n in (2, 3, 4):
    D = 3 ** n
    # (a), single excitations and tangent
    A = povm_from(rng.normal(size=8 * n), n); u = tangent_range([sqrt_psd(x) for x in A], np.zeros(2), n)
    ea = abs(m_povm(A, n) - m_direct(u, n))
    aa = rng.normal(size=2) + 1j * rng.normal(size=2); aa *= 0.6 / np.linalg.norm(aa)
    S = np.outer(aa.conj(), aa); Ab = povm_from(rng.normal(size=8 * n), n); w, V = np.linalg.eigh(np.eye(2) - S)
    R = V @ np.diag(np.sqrt(w)) @ V.conj().T; At = [R @ x @ R for x in Ab]          # |abar><abar| + sum At = I
    ut = tangent_range([sqrt_psd(x) for x in At], aa, n)
    G = np.array([[np.vdot(ut[r], ut[rp]) for rp in range(2)] for r in range(2)])
    ea = max(ea, abs(m_povm(At, n, aa) - m_direct(ut, n)), np.abs(G - np.eye(2)).max())
    # (b)
    Wl = [np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))[0] for _ in range(n)]
    u2 = tangent_range([w_ @ sqrt_psd(x) for w_, x in zip(Wl, A)], np.zeros(2), n)
    Vc = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))[0]
    eb = max(abs(np.linalg.eigvalsh(K_of(u2, n))[0] - m_povm(A, n)), abs(m_povm([Vc.conj().T @ x @ Vc for x in A], n) - m_povm(A, n)))
    # (c)
    st = strings(n); idx = {s: i for i, s in enumerate(st)}; U = [sqrt_psd(x) for x in A]; K = K_of(u, n); ec = 0; csum = []
    for c in range(n + 1):
        strs = [s for s in st if sum(x != 0 for x in s) == c]; mm = len(strs); pos = {s: i for i, s in enumerate(strs)}
        sel = [r * D + idx[s] for r in range(2) for s in strs]; Kc = K[np.ix_(sel, sel)]
        Kt = np.zeros((2 * mm, 2 * mm), complex); Wb = np.zeros((2 * mm, 2 * mm), complex); at = lambda r, s: r * mm + pos[s]
        for sg in strs:
            exc = [t for t in range(n) if sg[t] != 0]; vac = [t for t in range(n) if sg[t] == 0]
            for r in range(2):
                for rp in range(2):
                    Kt[at(r, sg), at(rp, sg)] += (r == rp) + sum(A[t][rp, r] for t in vac)
                    Wb[at(r, sg), at(rp, sg)] += 2 * (r == rp) - 3 * sum(A[t][rp, r] for t in exc)
            for t in exc:
                for zp in (1, 2):
                    sp = list(sg); sp[t] = zp; sp = tuple(sp)
                    for r in range(2):
                        for rp in range(2):
                            Kt[at(r, sg), at(rp, sp)] += -0.5 * U[t][sg[t] - 1, r] * np.conj(U[t][zp - 1, rp])
                for s in vac:
                    for z in (1, 2):
                        sp = list(sg); sp[t] = 0; sp[s] = z; sp = tuple(sp)
                        for r in range(2):
                            for rp in range(2):
                                Kt[at(r, sp), at(rp, sg)] += 0.5 * U[s][z - 1, r] * np.conj(U[t][sg[t] - 1, rp])
        ec = max(ec, np.abs(2 ** (n - c) * Kc - Kt).max())
        lm = np.linalg.eigvalsh((Kt + Kt.conj().T) / 2)[0]; lw = np.linalg.eigvalsh((Kt - Wb + (Kt - Wb).conj().T) / 2)[0]
        csum.append((c, lm, lw))
        ok &= lw > -1e-10 and (c != 0 or abs(lm - 2) < 1e-10) and (c != n or lm > -1e-10)
    # cross-sector blocks of K vanish
    exc_of = np.array([sum(x != 0 for x in s) for s in st]); exc2 = np.concatenate([exc_of, exc_of])
    cross = max(np.abs(K[np.ix_(exc2 == c1, exc2 == c2)]).max() for c1 in range(n + 1) for c2 in range(n + 1) if c1 != c2)
    # (d) extreme POVMs: rank-one elements (random n rank-one vectors completed to a POVM is not always possible; use
    #     two random orthonormal bases split among the copies, which gives rank-one elements summing to I)
    md = np.inf
    for _ in range(5):
        Qa = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))[0]; Qb = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))[0]
        vecs = [Qa[:, 0], Qa[:, 1], Qb[:, 0], Qb[:, 1]]; wts = [0.5, 0.5, 0.5, 0.5]
        E = [np.zeros((2, 2), complex) for _ in range(n)]
        for k, (v, wt) in enumerate(zip(vecs, wts)): E[k % n] += wt * np.outer(v, v.conj())
        if all(np.linalg.matrix_rank(e, 1e-9) <= 1 for e in E): md = min(md, m_povm(E, n))
    md = min(md, m_povm([np.eye(2)] + [np.zeros((2, 2))] * (n - 1), n))
    print(f"n={n}: (a) |m - direct| {ea:.1e}; (b) invariance {eb:.1e}; (c) sector formula {ec:.1e}, cross-sector blocks {cross:.1e}, "
          + "; ".join(f"c={c}: lambda_min {lm:+.3f}, weight-bound slack {lw:+.1e}" for c, lm, lw in csum) + f"; (d) min m at extreme POVMs {md:+.2e}", flush=True)
    ok &= ea < 1e-10 and eb < 1e-10 and ec < 1e-10 and cross < 1e-12 and md > -1e-12
assert ok
print("PASS")
