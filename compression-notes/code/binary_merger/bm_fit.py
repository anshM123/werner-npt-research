"""Fit of a 3 x 3 operator X by ONE atom lambda G(E, F) plus products (the SF form).

Stage A (parametric): joint L-BFGS over (E, F) and P product vectors, several random starts.
Stage B (variable projection): for fixed (E, F) the best weights solve a nonnegative least-squares problem over
[G(E, F)/Tr] + a product dictionary (Frobenius-isometric vectorisation).  The atom is moved by L-BFGS-B on
f(E, F) = min_{c >= 0} ||X - c_0 G^ - sum c_k P_k||^2 with the envelope gradient -2 c_0 Re Tr(R dG^), and products most
aligned with the residual R are added (column generation).  Residuals are ||R||_F / ||X||_F with Tr X = 1.
"""
import numpy as np
import torch
from scipy.optimize import minimize
from bm_core import (RHO_T, D, N, vecr, vecr_batch, rand_c, product_cols, best_products, SepDict, nnls_fit, G)


def norm_t(X):
    X = (X + X.conj().T) / 2
    return X / np.trace(X).real


# ------------------------------------------------------------------ stage A
def param_fit(X, P=80, restarts=4, seed=0, iters=4000, init=None, scale=0.4):
    """joint fit X ~ G(E, F) + sum_k |x_k y_k><x_k y_k|; returns (residual, E, F, xs, ys)"""
    torch.manual_seed(seed)
    Xt = torch.tensor(X, dtype=torch.complex128)
    nX = float(torch.linalg.norm(Xt))
    best = (np.inf, None)
    for r in range(restarts):
        if init is not None and r == 0:
            p = [torch.tensor(q, dtype=torch.complex128).clone().requires_grad_() for q in init]
        else:
            p = [(scale * torch.randn(*s, dtype=torch.complex128)).requires_grad_() for s in ((D, D), (D, D), (P, D), (P, D))]

        def model():
            V = torch.einsum('ka,kb->kab', p[2], p[3]).reshape(-1, N)
            K = torch.kron(p[0], p[1])
            return V.T @ V.conj() + K @ RHO_T @ K.conj().T
        o = torch.optim.LBFGS(p, lr=1, max_iter=iters, history_size=50, line_search_fn='strong_wolfe',
                              tolerance_grad=1e-15, tolerance_change=1e-18)

        def cl():
            o.zero_grad(); loss = torch.linalg.norm(model() - Xt) ** 2; loss.backward(); return loss
        o.step(cl)
        with torch.no_grad():
            rr = float(torch.linalg.norm(model() - Xt)) / nX
        if rr < best[0]:
            best = (rr, [q.detach().numpy().copy() for q in p])
    rr, pp = best
    return rr, pp[0], pp[1], pp[2], pp[3]


# ------------------------------------------------------------------ stage B
def atom_vec(th):
    E = (th[:9] + 1j * th[9:18]).reshape(D, D); Fm = (th[18:27] + 1j * th[27:36]).reshape(D, D)
    Y = G(E, Fm); return vecr(Y / np.trace(Y).real), E, Fm


def th_from(E, Fm):
    return np.concatenate([E.real.ravel(), E.imag.ravel(), Fm.real.ravel(), Fm.imag.ravel()])


def grad_atom(th, Rm, c0):
    """gradient of -2 c0 Re Tr(R G^(th)), G^ = G/Tr G, w.r.t. the 36 real parameters"""
    t = torch.tensor(th, requires_grad=True)
    E = (t[:9] + 1j * t[9:18]).reshape(D, D); Fm = (t[18:27] + 1j * t[27:36]).reshape(D, D)
    K = torch.kron(E, Fm); Y = K @ RHO_T @ K.conj().T; Y = Y / torch.trace(Y).real
    val = -2 * c0 * torch.trace(torch.tensor(Rm, dtype=torch.complex128) @ Y).real
    val.backward()
    return t.grad.numpy().copy()


class OneAtomFit:
    def __init__(self, X, seed=0, nrand=800):
        self.X = norm_t(X)
        self.x = vecr(self.X)
        self.nX = np.linalg.norm(self.x)
        self.g = np.random.default_rng(seed)
        self.dict = SepDict(self.g, nrand)
        self.th = None
        self.coef = None

    def solve(self, th):
        a, _, _ = atom_vec(th)
        cols = np.concatenate([a[:, None], self.dict.V], axis=1)
        coef, rn = nnls_fit(self.x, cols)
        r = self.x - cols @ coef
        return coef, r

    def fobj(self, th):
        coef, r = self.solve(th)
        f = float(r @ r)
        Rm = unvec(r)
        gr = grad_atom(th, Rm, coef[0])
        return f, gr

    def residual(self):
        coef, r = self.solve(self.th)
        self.coef = coef
        return np.linalg.norm(r) / self.nX, coef, r

    def colgen(self, nadd=12, nstart=24):
        _, coef, r = self.residual()
        Rm = unvec(r)
        new = best_products(Rm, self.g, nstart=nstart, keep=nadd)
        xs = [x for v, x, y in new if v > 0]; ys = [y for v, x, y in new if v > 0]
        if xs:
            self.dict.add_many(xs, ys)
        return max(v for v, _, _ in new)

    def atom_step(self, maxiter=60):
        res = minimize(self.fobj, self.th, jac=True, method='L-BFGS-B', options={'maxiter': maxiter, 'ftol': 1e-30, 'gtol': 1e-18})
        th = res.x
        # remove the scale gauge
        E = (th[:9] + 1j * th[9:18]).reshape(D, D); Fm = (th[18:27] + 1j * th[27:36]).reshape(D, D)
        s = np.sqrt(np.linalg.norm(E) / np.linalg.norm(Fm))
        self.th = th_from(E / s, Fm * s)

    def run(self, E0, F0, prods=None, rounds=30, verbose=False, tol=1e-14):
        self.th = th_from(E0, F0)
        if prods is not None:
            self.dict.add_many(prods[0], prods[1])
        hist = []
        for rd in range(rounds):
            self.atom_step()
            gain = self.colgen()
            rr, coef, _ = self.residual()
            hist.append(rr)
            if verbose:
                print(f"   round {rd}: residual {rr:.3e}, atom weight {coef[0]:.3f}, active {int(np.sum(coef > 0))}, "
                      f"dict {self.dict.V.shape[1]}, best new alignment {gain:+.2e}", flush=True)
            if rr < tol:
                break
            if self.dict.V.shape[1] > 3000:
                self.dict.prune(coef[1:], keep_extra=600)
        return hist[-1], hist


def unvec(v):
    from bm_core import unvecr
    return unvecr(v)


def fit_one(X, P=80, restarts=4, seed=0, rounds=30, verbose=False, nrand=800):
    """full pipeline; returns dict with residuals of both stages and the fitted atom"""
    X = norm_t(X)
    rA, E, Fm, xs, ys = param_fit(X, P=P, restarts=restarts, seed=seed)
    fo = OneAtomFit(X, seed=seed, nrand=nrand)
    rB, hist = fo.run(E, Fm, prods=(xs, ys), rounds=rounds, verbose=verbose)
    E = (fo.th[:9] + 1j * fo.th[9:18]).reshape(D, D); Fm = (fo.th[18:27] + 1j * fo.th[27:36]).reshape(D, D)
    return dict(param=rA, vp=rB, hist=hist, E=E, F=Fm, weight=fo.coef[0], fit=fo)


# ------------------------------------------------------------------ structured atom candidates
def vdm_frame(X, iters=300):
    """local filtering (a (x) b) X (a (x) b)^* with both marginals proportional to I (full-rank X)"""
    a = np.eye(D, dtype=complex); b = np.eye(D, dtype=complex)
    Y = X.copy()
    for _ in range(iters):
        Y4 = Y.reshape(D, D, D, D)
        rA = np.einsum('ijkj->ik', Y4); rB = np.einsum('ijil->jl', Y4)
        wA, UA = np.linalg.eigh(rA); wB, UB = np.linalg.eigh(rB)
        if min(wA.min(), wB.min()) <= 0:
            break
        sa = UA @ np.diag(wA ** -0.5) @ UA.conj().T; sb = UB @ np.diag(wB ** -0.5) @ UB.conj().T
        K = np.kron(sa, sb); Y = K @ Y @ K.conj().T; Y /= np.trace(Y).real
        a = sa @ a; b = sb @ b
        if np.abs(wA / wA.mean() - 1).max() < 1e-13 and np.abs(wB / wB.mean() - 1).max() < 1e-13:
            break
    return a, b, Y


def polar_u(M):
    U, s, Vh = np.linalg.svd(M); return U @ Vh


def candidates(X):
    """atom candidates (E, F) from the PT-negative direction, in the original and in the normal-form frame"""
    from bm_core import pt
    out = []
    w, U = np.linalg.eigh(pt(X)); M = U[:, 0].reshape(D, D)       # v = vec(E F^*)
    out += [(M, np.eye(D)), (np.eye(D), M.conj().T)]
    try:
        a, b, Y = vdm_frame(X)
        w, U = np.linalg.eigh(pt(Y)); Mn = U[:, 0].reshape(D, D)  # ~ vec(W^*) for (I (x) W) rho (I (x) W)^*
        W = polar_u(Mn).conj().T
        ai, bi = np.linalg.inv(a), np.linalg.inv(b)
        out += [(ai, bi @ W), (ai @ polar_u(Mn), bi)]
    except np.linalg.LinAlgError:
        pass
    return out


def fit_robust(X, P=80, restarts=4, seed=0, rounds=20, nbest=3, extra=(), verbose=False, nrand=800):
    """multi-start variable projection: parametric restarts + structured candidates + extra (E, F) pairs;
    the nbest starts (by residual after 3 VP rounds) are run to the end."""
    X = norm_t(X)
    starts = []
    rA, E, Fm, xs, ys = param_fit(X, P=P, restarts=restarts, seed=seed)
    starts.append(('param', E, Fm, (xs, ys)))
    for k, (E, Fm) in enumerate(list(candidates(X)) + list(extra)):
        starts.append((f'cand{k}', E, Fm, None))
    trial = []
    for k, (tag, E, Fm, prods) in enumerate(starts):
        fo = OneAtomFit(X, seed=seed + 17 * k, nrand=nrand)
        if prods is None:
            prods = (xs, ys)
        r3, _ = fo.run(E, Fm, prods=prods, rounds=3)
        trial.append((r3, tag, fo))
        if verbose:
            print(f"   start {tag}: residual after 3 rounds {r3:.3e}", flush=True)
    trial.sort(key=lambda z: z[0])
    best = None
    for r3, tag, fo in trial[:nbest]:
        E = (fo.th[:9] + 1j * fo.th[9:18]).reshape(D, D); Fm = (fo.th[18:27] + 1j * fo.th[27:36]).reshape(D, D)
        rB, hist = fo.run(E, Fm, rounds=rounds)
        if verbose:
            print(f"   start {tag}: final {rB:.3e}", flush=True)
        if best is None or rB < best[0]:
            best = (rB, tag, fo)
    rB, tag, fo = best
    E = (fo.th[:9] + 1j * fo.th[9:18]).reshape(D, D); Fm = (fo.th[18:27] + 1j * fo.th[27:36]).reshape(D, D)
    return dict(param=rA, vp=rB, start=tag, E=E, F=Fm, weight=fo.coef[0], fit=fo)
