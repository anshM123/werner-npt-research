"""Checks of three exact facts used in the README.
 (1) Scale invariance: for A in GL_3 at copy s, Q_{A C A^*}(H) = Q_C(A^* H A) and N_{A C A^*}(H) = N_C(A^* H A).
 (2) Thin-point tilt identity: if the range of C is thin at copy s with support V and f perp V, then for h in V and
     H = |h><f| + |f><h|:  N_s(H) = 0 and Q_s(H) = -1/2 q_{n-1}(<h|_s C |f>_s).
 (3) Phase-averaged tilts: for H_theta = e^{i theta}|f><w| + h.c. (f, w orthonormal),
     E_theta Q = 2Re q'(C_ff, C_ww) - 1/2[q'(C_fw) + q'(C_wf)],  E_theta N = 2Re<C_ff, C_ww>."""
import numpy as np
from common import q_direct, squeeze_forms, Bt, tangent_range, d
rng = np.random.default_rng(3)


def herm_coords(H):
    return np.array([np.trace(Bk.conj().T @ H).real for Bk in Bt])


def block(X, n, s, x, y):
    T = X.reshape([d] * (2 * n)); T = np.tensordot(x.conj(), T, axes=([0], [s])); T = np.tensordot(T, y, axes=([n - 1 + s], [0]))
    return T.reshape(d ** (n - 1), d ** (n - 1))


def qpol(A, B, m):
    return 0.25 * (q_direct(A + B, m) - q_direct(A - B, m) - 1j * q_direct(A + 1j * B, m) + 1j * q_direct(A - 1j * B, m))


def act(A, X, n, s):
    T = X.reshape([d] * (2 * n)); T = np.moveaxis(np.tensordot(A, T, axes=([1], [s])), 0, s)
    T = np.moveaxis(np.tensordot(T, A.conj().T, axes=([n + s], [0])), -1, n + s); return T.reshape(d ** n, d ** n)


worst = 0.0
for n in (2, 3):
    D = d ** n
    for t in range(3):
        U = [rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(n)]; a = rng.normal(size=2) + 1j * rng.normal(size=2)
        u = tangent_range(U, a, n); v = rng.normal(size=(2, D)) + 1j * rng.normal(size=(2, D))
        X = np.outer(u[0], v[0].conj()) + np.outer(u[1], v[1].conj()); s = t % n
        # (1)
        A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); H = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); H = H + H.conj().T
        K1, N1 = squeeze_forms(act(A, X, n, s), n, s); K0, N0 = squeeze_forms(X, n, s)
        h1 = herm_coords(H); h0 = herm_coords(A.conj().T @ H @ A)
        e1 = abs(h1 @ K1 @ h1 - h0 @ K0 @ h0) / max(1, abs(h1 @ K1 @ h1)); e2 = abs(h1 @ N1 @ h1 - h0 @ N0 @ h0) / max(1, abs(h1 @ N1 @ h1))
        # (2) thin at s: make U_s rank one (then V = span{|0>, y})
        Ut = [x.copy() for x in U]; yv = rng.normal(size=2) + 1j * rng.normal(size=2); Ut[s] = np.outer(yv, rng.normal(size=2) + 1j * rng.normal(size=2))
        ut = tangent_range(Ut, a, n); Xt = np.outer(ut[0], v[0].conj()) + np.outer(ut[1], v[1].conj())
        y = np.zeros(3, complex); y[1:] = yv / np.linalg.norm(yv); f = np.zeros(3, complex); f[1:] = [-np.conj(y[2]), np.conj(y[1])]
        hv = 0.6 * np.array([1, 0, 0]) + (0.3 + 0.5j) * y
        Hh = np.outer(hv, f.conj()) + np.outer(f, hv.conj()); Kt, Nt = squeeze_forms(Xt, n, s); hh = herm_coords(Hh)
        e3 = abs(hh @ Nt @ hh) / max(1, np.linalg.norm(Nt)); e4 = abs(hh @ Kt @ hh + 0.5 * q_direct(block(Xt, n, s, hv, f), n - 1)) / max(1, np.linalg.norm(Kt))
        # (3) phase average
        fv = np.linalg.qr(rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2)))[0]; f3, w3 = fv[:, 0], fv[:, 1]
        thetas = np.linspace(0, 2 * np.pi, 8, endpoint=False)
        EQ = np.mean([herm_coords(np.exp(1j * th) * np.outer(f3, w3.conj()) + np.exp(-1j * th) * np.outer(w3, f3.conj())) @ K0 @
                      herm_coords(np.exp(1j * th) * np.outer(f3, w3.conj()) + np.exp(-1j * th) * np.outer(w3, f3.conj())) for th in thetas])
        EN = np.mean([herm_coords(np.exp(1j * th) * np.outer(f3, w3.conj()) + np.exp(-1j * th) * np.outer(w3, f3.conj())) @ N0 @
                      herm_coords(np.exp(1j * th) * np.outer(f3, w3.conj()) + np.exp(-1j * th) * np.outer(w3, f3.conj())) for th in thetas])
        Cff, Cww, Cfw, Cwf = block(X, n, s, f3, f3), block(X, n, s, w3, w3), block(X, n, s, f3, w3), block(X, n, s, w3, f3)
        EQf = 2 * qpol(Cff, Cww, n - 1).real - 0.5 * (q_direct(Cfw, n - 1) + q_direct(Cwf, n - 1)); ENf = 2 * np.vdot(Cff, Cww).real
        e5 = abs(EQ - EQf) / max(1, abs(EQ)); e6 = abs(EN - ENf) / max(1, abs(EN))
        errs = (e1, e2, e3, e4, e5, e6); worst = max(worst, max(errs))
        print(f"n={n} copy {s}: scale-invariance {e1:.1e} {e2:.1e} | thin tilt N {e3:.1e} Q {e4:.1e} | phase average Q {e5:.1e} N {e6:.1e}")
print(f"RESULT: maximal relative error {worst:.1e}")
assert worst < 1e-9
print("PASS")
