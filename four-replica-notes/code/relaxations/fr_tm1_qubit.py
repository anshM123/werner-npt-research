"""Facial reduction of the two-letter relaxation (herm4) at t = -1, n = 4, with the qubit spin-flip zero family.

Zero families used (all checked directly via q_4):
  k2, k4      u, v orthogonal products on 2 (4) sites sharing psi on the rest           (fr2_herm4.zero_point_k)
  inner       u = a (x) phi, v = a' (x) chi, (phi, chi) an orthonormal basis of span(b (x) c (x) d, b' (x) c' (x) d')
  spinflip    u in Q_1 (x) ... (x) Q_4 (Q_i 2-dim), v = Y conj(u), Y = sigma_y^(x)4 in the Q_i bases, <u, Y conj(u)> = 0.
              In the qubit sector q_n(C) = 2^n * (weight of C on full-weight Pauli strings), and conjugation by the spin
              flip keeps every full-weight correlation for even n, so P_u - P_v has none: q_4 = 0.
Runs the reduced problem for several constraint filters, with real zeros only and with complex zeros added, and
measures how far complex zeros stick out of the face of the real ones.
usage: python fr_tm1_qubit.py [n_spinflip_real] [n_spinflip_complex] [filters]"""
import sys, os, time, itertools, numpy as np
os.environ['FR_REAL'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from alg4 import build_algebras
from herm4 import Herm4, FAMS, moments_multi, ALL_CONS
from validate import F1_matrix, zfull_from_moments
from rep4 import unsvec
import fr2_herm4 as FR
from tm1_analyze import qform

SY = np.array([[0, -1j], [1j, 0]])

def embed(uq, rots):
    """(C^2)^(x)4 -> (C^3)^(x)4 with site isometries rots[i] (3 x 2)"""
    B = rots[0]
    for R in rots[1:]: B = np.kron(B, R)
    return B @ uq

def spinflip_pair(rng, real):
    Y = SY
    for _ in range(3): Y = np.kron(Y, SY)
    if real:
        Y = Y.real
        w, V = np.linalg.eigh(Y); Ep = V[:, w > 0]; Em = V[:, w < 0]
        a = Ep @ rng.standard_normal(Ep.shape[1]); b = Em @ rng.standard_normal(Em.shape[1])
        a /= np.linalg.norm(a); b /= np.linalg.norm(b); uq = (a + b) / np.sqrt(2)
        vq = Y @ uq
    else:
        a = rng.standard_normal(16) + 1j * rng.standard_normal(16); b = rng.standard_normal(16) + 1j * rng.standard_normal(16)
        A2, A1, A0 = b @ Y @ b, 2 * (a @ Y @ b), a @ Y @ a           # (a + l b)^T Y (a + l b) = 0
        lam = np.roots([A2, A1, A0])[0]; uq = a + lam * b; uq /= np.linalg.norm(uq)
        vq = Y @ uq.conj()
    rots = []
    for i in range(4):
        if real: Q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
        else: Q, _ = np.linalg.qr(rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)))
        rots.append(Q[:, :2])
    return embed(uq, rots), embed(vq, rots)

def inner_pair(rng, real):
    """u = a (x) phi, v = a' (x) chi, (phi, chi) a rotated orthonormal basis of span(bcd, b'c'd')"""
    def onb():
        Q, _ = np.linalg.qr(rng.standard_normal((3, 2)) + (0 if real else 1j * rng.standard_normal((3, 2)))); return Q
    p = np.ones(1); q = np.ones(1)
    for s in range(3):
        Q = onb(); p = np.kron(p, Q[:, 0]); q = np.kron(q, Q[:, 1])
    th = rng.uniform(0, 2 * np.pi)
    phi = np.cos(th) * p + np.sin(th) * q; chi = -np.sin(th) * p + np.cos(th) * q
    Qa = onb()
    return np.kron(Qa[:, 0], phi), np.kron(Qa[:, 1], chi)

def x_of(rel, algs, u, v, n=4):
    zf = {c: zfull_from_moments(moments_multi([u.reshape(-1, 1) if ch == 'u' else v.reshape(-1, 1) for ch in c], n), n, algs[0])
          for c in FAMS}
    x, res = rel.x_from_zfulls(zf)
    return x, res

if __name__ == "__main__":
    nsr = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    nsc = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    filters = sys.argv[3].split(';') if len(sys.argv) > 3 else ['all', 'COV', 'PT', 'COV+PT1']
    n = 4; t = -1.0
    algs = build_algebras(); F1 = F1_matrix(algs[0]); t0 = time.time()
    rel = Herm4(n, algs, ALL_CONS, verbose=False)
    print(f"n={n} t={t:+.0f}: build {time.time() - t0:.0f}s ({rel.N} parameters, {len(rel.blocks)} PSD blocks)", flush=True)
    rng = np.random.default_rng(2026)
    xs_struct = []
    for k in (2, 4):
        for j in range(2):
            zf = FR.zero_point_k(n, k, rng, F1); x, res = rel.x_from_zfulls(zf)
            mins, mom, eqr, obj = rel.evaluate(x, 1.0, t)
            print(f"   k={k} #{j}: reconstruction {res:.1e}, q {obj:.1e}, min eig {min(mins.values()):.1e}", flush=True)
            xs_struct.append(x)
    for j in range(4):
        u, v = inner_pair(rng, True); qd = qform(np.outer(u, u.conj()) - np.outer(v, v.conj()), n)
        x, res = x_of(rel, algs, u, v); mins, mom, eqr, obj = rel.evaluate(x, 1.0, t)
        print(f"   inner #{j}: q_4 direct {qd:.1e}, reconstruction {res:.1e}, q {obj:.1e}, min eig {min(mins.values()):.1e}", flush=True)
        xs_struct.append(x)
    xs_sr = []
    for j in range(nsr):
        u, v = spinflip_pair(rng, True); qd = qform(np.outer(u, u.conj()) - np.outer(v, v.conj()), n)
        x, res = x_of(rel, algs, u, v); mins, mom, eqr, obj = rel.evaluate(x, 1.0, t)
        print(f"   spinflip real #{j}: <u,v> {abs(np.vdot(u, v)):.1e}, q_4 direct {qd:.1e}, reconstruction {res:.1e}, q {obj:.1e}, "
              f"min eig {min(mins.values()):.1e}", flush=True)
        xs_sr.append(x)
    xs_sc = []
    for j in range(nsc):
        u, v = spinflip_pair(rng, False); qd = qform(np.outer(u, u.conj()) - np.outer(v, v.conj()), n)
        x, res = x_of(rel, algs, u, v); mins, mom, eqr, obj = rel.evaluate(x, 1.0, t)
        print(f"   spinflip complex #{j}: <u,v> {abs(np.vdot(u, v)):.1e}, q_4 direct {qd:.1e}, reconstruction {res:.1e}, "
              f"q {obj:.1e}, min eig {min(mins.values()):.1e}", flush=True)
        xs_sc.append(x)
    print(f"   zeros built; {time.time() - t0:.0f}s", flush=True)
    # face of the structured + real spin-flip zeros; how far do the complex ones stick out?
    kers_r = FR.kernels(rel, sum(xs_struct) + sum(xs_sr))
    worst = 0.0
    for x in xs_sc:
        for (name, Mf, Nk) in kers_r:
            if Nk.shape[1] == 0: continue
            m = Nk.shape[0]; W = unsvec(Mf @ x, m); worst = max(worst, np.abs(Nk.T @ W @ Nk).max())
    print(f"   complex spin-flip zeros on the real face: max |N^T W(x) N| = {worst:.1e}", flush=True)
    for f in filters:
        FR.solve_reduced(rel, t, FR.kernels(rel, sum(xs_struct)), FR.FILTERS[f], label=f'structured only, {f}')
        FR.solve_reduced(rel, t, kers_r, FR.FILTERS[f], label=f'structured + {nsr} real spin-flip, {f}')
    if nsc:
        kers_c = FR.kernels(rel, sum(xs_struct) + sum(xs_sr) + sum(xs_sc))
        for f in filters:
            FR.solve_reduced(rel, t, kers_c, FR.FILTERS[f], label=f'structured + real + {nsc} complex spin-flip, {f}')
    print(f"done; {time.time() - t0:.0f}s", flush=True)
