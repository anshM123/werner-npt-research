"""Strong adversary on the tangent class, with analytic gradients through the eigenvectors (checked against finite
differences to about 1e-6).  C = sum_r u_r v_r^* with u_r = a_r e + sum_t (U_t e_r)@t and free co-range.  Adam steps on
f = softmin over copies of lambda_min(Q_s - lam_s N_s)/||C||^2, minus a penalty that keeps the range and co-range
local-support ratios at every copy above tau.  Modes:
  common  one lam in [-1/2, -0.01] shared by all copies (Condition SC'); f > 0 refutes it;
  exists  each copy its own best lam in [-8, 8] (Condition SC'-inf); f > 0 refutes it;
  min0    copy 0 only, lam = R(C): f = 0 means C is a minimum of R on the squeeze orbit of copy 0;
  minall  every copy, lam = R(C): f = 0 means C is a squeeze minimum at every copy (refutes SC-cw).
usage: python search_strong.py n mode seed iters tau [start-parameter-file]"""
import numpy as np, sys, itertools, time

from common import squeeze_forms, tangent_range, q, Bt, Lmap
from scipy.optimize import minimize_scalar
n, mode, seed, iters, tau = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
d = 3; D = d ** n; rng = np.random.default_rng(seed); T_soft = 0.003
st = list(itertools.product(range(d), repeat=n)); idx = {s: i for i, s in enumerate(st)}
def embed_copy(Hm, s, k):
    ops = [np.eye(3)] * k; ops[s] = Hm; out = ops[0]
    for o in ops[1:]: out = np.kron(out, o)
    return out
def ptr(X, J, k):
    """partial trace over the copies in J of a 3^k x 3^k matrix; returns (matrix, remaining copies)"""
    Tt = X.reshape([3] * (2 * k)); cur = list(range(k))
    for t in sorted(J, reverse=True):
        i = cur.index(t); m = len(cur); Tt = np.trace(Tt, axis1=i, axis2=m + i); cur.pop(i)
    m = len(cur); return Tt.reshape(3 ** m, 3 ** m), cur
def lift(Y, cur, J):
    """adjoint of the partial trace: Y on copies cur -> Y (x) I_J on all n copies (copy order restored)"""
    k = len(cur); full = np.kron(Y, np.eye(3 ** len(J))).reshape([3] * (2 * n)); order = cur + sorted(J)
    perm = np.argsort(order); return full.transpose(list(perm) + [p + n for p in perm]).reshape(D, D)
def phi_grad(X, s, Hm, lam):
    """value and dPhi/dXbar of Phi(X) = (Q_s - lam N_s)(Hm) for fixed Hermitian Hm at copy s"""
    others = [t for t in range(n) if t != s]; val = 0; G = np.zeros((D, D), complex); Hs = embed_copy(Hm, s, n)
    XH = X @ Hs
    for k in range(n):
        for Jp in itertools.combinations(others, k):
            c = (-0.5) ** k
            Y, cur = ptr(X, list(Jp), n); Hl = embed_copy(Hm, cur.index(s), len(cur)); HYH = Hl @ Y @ Hl
            val += c * np.vdot(Y, HYH).real; G += c * lift(HYH, cur, list(Jp))
            Z, cur2 = ptr(XH, list(Jp) + [s], n); val += -0.5 * c * np.vdot(Z, Z).real
            G += -0.5 * c * lift(Z, cur2, list(Jp) + [s]) @ Hs
    HXH = Hs @ X @ Hs; val -= lam * np.vdot(X, HXH).real; G -= lam * HXH
    return val, G
def local_ops(X):
    out = []
    for s in range(n):
        for side in (0, 1):
            Y = X @ X.conj().T if side == 0 else X.conj().T @ X
            Tt = np.moveaxis(Y.reshape([3] * (2 * n)), [s, n + s], [0, n]).reshape(3, D // 3, 3, D // 3)
            ev, V = np.linalg.eigh(np.einsum('aibi->ab', Tt)); out.append((s, side, ev, V))
    return out
def unpack(x):
    z = x[:len(x) // 2] + 1j * x[len(x) // 2:]
    U = [z[4 * t:4 * t + 4].reshape(2, 2) for t in range(n)]; a = z[4 * n:4 * n + 2]; V = z[4 * n + 2:].reshape(2, D)
    return U, a, V
def objective(x, lam_c, want_grad=True):
    U, a, V = unpack(x); u = tangent_range(U, a, n); X = np.outer(u[0], V[0].conj()) + np.outer(u[1], V[1].conj())
    nrm = np.vdot(X, X).real; vals = []; grads = []; lams = []
    LX = Lmap(X, n); Rc = np.vdot(X, LX).real / nrm; dR = (LX - Rc * X) / nrm
    for s in ([0] if mode == "min0" else range(n)):
        K, Kn = squeeze_forms(X, n, s)
        if mode in ("min0", "minall"): lam = Rc
        elif mode == "common": lam = lam_c
        else:
            r = minimize_scalar(lambda l: -np.linalg.eigvalsh(K - l * Kn)[0], bounds=(-8, 8), method='bounded', options={'xatol': 1e-9}); lam = r.x
        w, Vv = np.linalg.eigh(K - lam * Kn); h = Vv[:, 0]; Hm = np.einsum('k,kab->ab', h, Bt)
        if want_grad:
            val, G = phi_grad(X, s, Hm, lam)
            if mode in ("min0", "minall"): G = G - (h @ Kn @ h) * dR
            grads.append((G - val / nrm * X) / nrm)
        vals.append(w[0] / nrm); lams.append((lam, -(h @ Kn @ h) / nrm))
    vals = np.array(vals); wts = np.exp(-(vals - vals.min()) / T_soft); wts /= wts.sum()
    f = vals.min() - T_soft * np.log(np.sum(np.exp(-(vals - vals.min()) / T_soft)))
    # thickness penalty (soft-min over the 2n local-support ratios)
    lo = local_ops(X); ratios = np.array([ev[0] / ev[-1] for _, _, ev, _ in lo])
    pen_w = 200.0; deficit = np.maximum(0, tau - ratios); f_pen = pen_w * np.sum(deficit ** 2); f_tot = f - f_pen
    if not want_grad: return f_tot, f, ratios.min(), X
    GX = sum(wt * g for wt, g in zip(wts, grads))
    for (s, side, ev, Vv), dfc in zip(lo, deficit):
        if dfc <= 0: continue
        Pmin = embed_copy(np.outer(Vv[:, 0], Vv[:, 0].conj()), s, n); Pmax = embed_copy(np.outer(Vv[:, -1], Vv[:, -1].conj()), s, n)
        dmin = Pmin @ X if side == 0 else X @ Pmin; dmax = Pmax @ X if side == 0 else X @ Pmax
        dratio = (dmin * ev[-1] - dmax * ev[0]) / ev[-1] ** 2
        GX += pen_w * 2 * dfc * dratio            # d(-pen)/dXbar = +2 pen_w deficit * dratio
    g_lam = sum(wt * lm[1] for wt, lm in zip(wts, lams)) if mode == "common" else 0.0
    gu = [GX @ V[r] for r in range(2)]; gv = [GX.conj().T @ u[r] for r in range(2)]
    gz = np.zeros(len(x) // 2, complex)
    for t in range(n):
        for z_ in (1, 2):
            tt = [0] * n; tt[t] = z_; i = idx[tuple(tt)]
            for r in range(2): gz[4 * t + 2 * (z_ - 1) + r] = gu[r][i]
    gz[4 * n:4 * n + 2] = [gu[0][0], gu[1][0]]; gz[4 * n + 2:] = np.concatenate(gv)
    return f_tot, f, ratios.min(), 2 * np.concatenate([gz.real, gz.imag]), g_lam
P = 2 * (4 * n + 2 + 2 * D); x = rng.normal(size=P); lam_c = -0.25
if len(sys.argv) > 6: x = np.load(sys.argv[6]) + 0.02 * rng.normal(size=P)
m1 = np.zeros(P); m2 = np.zeros(P); lr = 0.02; t0 = time.time(); best = (-np.inf, None)
for it in range(1, iters + 1):
    f_tot, f, th, g, g_lam = objective(x, lam_c)
    if th >= tau - 1e-3 and f > best[0]: best = (f, x.copy(), lam_c, th)
    m1 = 0.9 * m1 + 0.1 * g; m2 = 0.999 * m2 + 0.001 * g ** 2
    x = x + lr * (m1 / (1 - 0.9 ** it)) / (np.sqrt(m2 / (1 - 0.999 ** it)) + 1e-8)
    if mode == "common": lam_c = float(np.clip(lam_c + 0.01 * np.sign(g_lam), -0.5, -0.01))
    x /= np.linalg.norm(x) / np.sqrt(P)
    if it % 100 == 0: print(f"it {it}: f = {f:+.4e}, thickness {th:.3f}, lam {lam_c:+.3f}, {time.time() - t0:.0f}s", flush=True)
f, xb, lb, thb = best if best[1] is not None else (f, x, lam_c, th)
print(f"n={n} {mode} seed {seed} tau {tau}: best f = {f:+.4e} at thickness {thb:.3f} (lam {lb:+.3f})", flush=True)
U, a, V = unpack(xb); u = tangent_range(U, a, n); X = np.outer(u[0], V[0].conj()) + np.outer(u[1], V[1].conj()); X /= np.linalg.norm(X)
