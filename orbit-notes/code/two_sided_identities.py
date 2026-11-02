"""Check the second variations at a critical point of R over ALL rank-two operators (t = 1/2).

Let C be critical for R on the rank-two variety, r = R(C), L = Phi_r(C) = Phi(C) - r C (then L vanishes on ran C*
and has range orthogonal to ran C).  For a 1-local Hermitian K = sum_s K_s:
  squeeze    e^{eK} C e^{eK} :  q_r = eps^2 [ 2 q_r(KC) + 2 q_r(CK) + 4 b(K) ] + O(eps^3)
  similarity e^{eK} C e^{-eK}:  q_r = eps^2 [ -4 b(K) ] + O(eps^3),       b(K) = Re <L, K C K>
and for K at a single copy s, in the eigenbasis of K_s (needs only criticality for squeezes at copy s):
  b(K_s) = -(1/2) sum_ab (k_a - k_b)^2 q'_r(X_ab),   X_ab = <a|_s C |b>_s.
Both curves stay in rank two, so at a global minimum of R both second variations are >= 0.
usage: python two_sided_identities.py   (run from orbit-notes/)"""
import sys, itertools
import numpy as np
import torch
sys.path.insert(0, 'code')
n = 2
sys.argv = ['x', str(n), '0.5', '0.1', '0', '0']
import squeeze_hessian as SH

torch.set_default_dtype(torch.float64)
C = np.load('data/variety_critical_n2.npy'); C = C / np.linalg.norm(C)
r = SH.q_np(C, n)
e9 = np.eye(3).reshape(9) / np.sqrt(3)


def Phi(Y):
    X = Y.reshape(3, 3, 3, 3).transpose(0, 2, 1, 3).reshape(9, 9)
    for s in range(n):
        tr = np.tensordot(e9.conj(), X, axes=([0], [s]))
        X = X - 1.5 * np.moveaxis(np.tensordot(e9, tr, axes=0), 0, s)
    return X.reshape(3, 3, 3, 3).transpose(0, 2, 1, 3).reshape(9, 9)


def qr_(Y): return SH.q_np(Y, n) - r * np.vdot(Y, Y).real


L = Phi(C) - r * C
U, S, Vh = np.linalg.svd(C); PA = U[:, :2] @ U[:, :2].conj().T; PB = Vh[:2].conj().T @ Vh[:2]
print(f"R = {r:+.6f}; criticality on the rank-two variety: |Pi_A L| + |L Pi_B| = {np.linalg.norm(PA @ L) + np.linalg.norm(L @ PB):.1e}")
rng = np.random.default_rng(2)


def coef2(curve, h=1e-3):
    return (qr_(curve(h)) - 2 * qr_(curve(0.0)) + qr_(curve(-h))) / (2 * h * h)


def expmh(M, e):
    w, V = np.linalg.eigh(M); return (V * np.exp(e * w)) @ V.conj().T


for trial in range(3):
    Ks = []
    for s in range(n):
        M = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)); Ks.append((M + M.conj().T) / 2)
    K = np.kron(Ks[0], np.eye(3)) + np.kron(np.eye(3), Ks[1])
    b = np.vdot(L, K @ C @ K).real
    sq = coef2(lambda e: expmh(K, e) @ C @ expmh(K, e)); si = coef2(lambda e: expmh(K, e) @ C @ expmh(K, -e))
    print(f"   direction {trial}: squeeze  finite difference {sq:+.6f}, formula {2 * qr_(K @ C) + 2 * qr_(C @ K) + 4 * b:+.6f};"
          f"  similarity  finite difference {si:+.6f}, formula {-4 * b:+.6f}")
w, V = np.linalg.eigh(Ks[0]); K0 = np.kron(Ks[0], np.eye(3))
T4 = C.reshape(3, 3, 3, 3); Xr = np.einsum('ia,ibjc,jd->adbc', V.conj(), T4, V)   # X_ab in the eigenbasis of K_0
qp = lambda Y: np.vdot(Y, Y).real - 0.5 * abs(np.trace(Y)) ** 2 - r * np.vdot(Y, Y).real
val = -0.5 * sum((w[a] - w[b]) ** 2 * qp(Xr[a, b]) for a in range(3) for b in range(3))
print(f"   single copy: b(K_0) = {np.vdot(L, K0 @ C @ K0).real:+.10f}, -(1/2) sum (k_a - k_b)^2 q'_r(X_ab) = {val:+.10f}")
