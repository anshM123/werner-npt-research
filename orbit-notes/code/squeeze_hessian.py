"""Squeeze-critical points and their second variation.

Finds a thick operator C of rank two that is critical for R_t on its squeeze orbit, aiming at a chosen value of R
(minimise |squeeze gradient|^2 + bias (R - target)^2 with a thickness floor, then polish with the gradient alone), and
reports:
  * the squeeze Hessian (log coordinates H_s, G_s = expm(H_s)) and, along its most negative direction, the split into
    the single-copy blocks and the cross-copy part;
  * the rank-one test  min_{s,e} F_s(ee*) = (1-t) q'(X_ee) - r ||X_ee||^2  (negative => a rank-one squeeze descends);
  * the plane test  min_{s,K} Tr_{su(K)} F_s  from the closed formula in NOTES.md;
  * checks of two identities: Tr_HS F_s = q'_r(Tr_s C) - t sum_ab q'(X_ab), and F_s(I - E) = F_s(E).

usage: python squeeze_hessian.py n t target starts seed [save_dir]
"""
import sys, itertools
import numpy as np
import torch
from scipy.optimize import minimize

torch.set_default_dtype(torch.float64)
n, t, TARGET, STARTS, SEED = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
SAVE = sys.argv[6] if len(sys.argv) > 6 else None
D = 3 ** n
TAU = 0.02        # thickness floor used only in the biased phase


def q_np(Y, m):
    """q_t(Y) = sum_J (-t)^|J| ||Tr_J Y||^2 for an operator Y on (C^3)^m."""
    T = Y.reshape((3,) * (2 * m)); tot = 0.0
    for J in itertools.product([0, 1], repeat=m):
        Z = T
        for s in reversed(range(m)):
            if J[s]:
                Z = np.trace(Z, axis1=s, axis2=Z.ndim // 2 + s)
        tot += (-t) ** sum(J) * np.vdot(Z, Z).real
    return tot


def b_np(Y, Z, m):
    """Sesquilinear form of q (conjugate-linear in Y)."""
    return 0.25 * (q_np(Y + Z, m) - q_np(Y - Z, m)) + 0.25j * (q_np(Y - 1j * Z, m) - q_np(Y + 1j * Z, m))


e9 = torch.eye(3, dtype=torch.complex128).reshape(9) / np.sqrt(3)


def R_t(C):
    X = C.reshape((3,) * (2 * n)).permute([i for s in range(n) for i in (s, n + s)]).reshape((9,) * n); Y = X
    for s in range(n):
        tr = torch.tensordot(e9.conj(), Y, dims=([0], [s]))
        Y = Y - 3 * t * torch.moveaxis(torch.tensordot(e9, tr, dims=0), 0, s)
    return torch.vdot(X.reshape(-1), Y.reshape(-1)).real / torch.vdot(X.reshape(-1), X.reshape(-1)).real


def herm8(x):
    M = torch.zeros(3, 3, dtype=torch.complex128) + torch.diag(torch.stack([x[0], x[1], -x[0] - x[1]]).to(torch.complex128))
    iu = torch.triu_indices(3, 3, 1)
    M[iu[0], iu[1]] = torch.complex(x[2:5], x[5:8]); M[iu[1], iu[0]] = torch.complex(x[2:5], -x[5:8])
    return M


def squeeze_R(C, h):
    G = torch.ones(1, 1, dtype=torch.complex128)
    for s in range(n):
        G = torch.kron(G, torch.linalg.matrix_exp(herm8(h[8 * s:8 * s + 8])))
    return R_t(G @ C @ G)


def C_of(p):
    A = torch.complex(p[:2 * D], p[2 * D:4 * D]).reshape(D, 2); B = torch.complex(p[4 * D:6 * D], p[6 * D:]).reshape(D, 2)
    C = A @ B.conj().T
    return C / torch.linalg.norm(C), A, B


def thickness(A, B):
    """Smallest eigenvalue of the one-copy reductions of the projectors onto range and co-range (min over copies)."""
    vals = []
    for M in (A, B):
        Q, _ = torch.linalg.qr(M); P = (Q @ Q.conj().T).reshape((3,) * (2 * n))
        for s in range(n):
            letters = 'abcdefghijklmnop'; inp = [letters[i] for i in range(2 * n)]
            for a in range(n):
                if a != s:
                    inp[n + a] = inp[a]
            vals.append(torch.linalg.eigvalsh(torch.einsum(''.join(inp) + '->' + letters[s] + inp[n + s], P))[0])
    return torch.stack(vals).min()


def blocks(C, s):
    """X[a, b] = <a|_s C |b>_s, operators on the other n-1 copies."""
    T = C.reshape((3,) * (2 * n))
    return np.moveaxis(T, [s, n + s], [0, 1]).reshape(3, 3, 3 ** (n - 1), 3 ** (n - 1))


def F_form(C, s, r, P):
    """F_s(P) = N(P) - r D(P), P Hermitian at copy s (other copies at I), via the block formula."""
    X = blocks(C, s); m = n - 1
    first = sum(P[c, a] * P[b, d] * (b_np(X[c, d], X[a, b], m) - r * np.vdot(X[c, d], X[a, b]))
                for a in range(3) for b in range(3) for c in range(3) for d in range(3))
    Y = sum(P[b, a] * X[a, b] for a in range(3) for b in range(3))
    return (first - t * q_np(Y, m)).real


def rank_one_test(C, r, rng):
    best = np.inf
    for s in range(n):
        X = blocks(C, s)

        def f(x):
            u = x[:3] + 1j * x[3:]; u = u / np.linalg.norm(u); Xu = np.einsum('i,ijxy,j->xy', u.conj(), X, u)
            return (1 - t) * q_np(Xu, n - 1) - r * np.vdot(Xu, Xu).real
        best = min(best, min(minimize(f, rng.standard_normal(6), method='Nelder-Mead', options={'maxiter': 3000}).fun
                             for _ in range(8)))
    return best


def plane_test(C, r, rng):
    m = n - 1; qr_ = lambda Y: q_np(Y, m) - r * np.vdot(Y, Y).real; best = np.inf
    for s in range(n):
        X = blocks(C, s)

        def f(x):
            e = x[:3] + 1j * x[3:]; e = e / np.linalg.norm(e)
            Q, _ = np.linalg.qr(np.column_stack([e, np.eye(3)[:, :2]]).astype(complex)); u1, u2 = Q[:, 1], Q[:, 2]
            Xc = lambda a, b: np.einsum('i,ijxy,j->xy', a.conj(), X, b)
            X11, X22, X12, X21 = Xc(u1, u1), Xc(u2, u2), Xc(u1, u2), Xc(u2, u1)
            Sg, Dl = X11 + X22, X11 - X22
            return (1.5 * qr_(Sg) - 0.5 * qr_(Dl) - t * q_np(Dl, m) - (1 + 2 * t) * (q_np(X12, m) + q_np(X21, m))
                    + r * (np.vdot(X12, X12).real + np.vdot(X21, X21).real))
        best = min(best, min(minimize(f, rng.standard_normal(6), method='Nelder-Mead', options={'maxiter': 3000}).fun
                             for _ in range(8)))
    return best


def to8(Hm):
    """Coordinates of a traceless Hermitian 3 x 3 matrix in the parametrisation herm8."""
    return np.array([Hm[0, 0].real, Hm[1, 1].real, Hm[0, 1].real, Hm[0, 2].real, Hm[1, 2].real,
                     Hm[0, 1].imag, Hm[0, 2].imag, Hm[1, 2].imag])


def restricted_min(Hess, V):
    """Smallest eigenvalue of Hess restricted to the column span of V (Euclidean metric of the coordinates; only the
    sign is coordinate-independent)."""
    A = V.T @ Hess @ V; G = V.T @ V; Li = np.linalg.inv(np.linalg.cholesky(G))
    return np.linalg.eigvalsh(Li @ A @ Li.T)[0]


def top_two_test(C, Hess):
    """Squeeze Hessian restricted to (i) the n directions H_s = |b_s><b_s| - |a_s><a_s| (a_s, b_s: largest and second
    largest eigenvectors of the local support Tr_other(CC* + C*C) at copy s; weights optimised jointly), (ii) the same
    one copy at a time, (iii) all traceless squeezes diagonal in the local-support eigenbases (2 per copy, jointly)."""
    V, W = [], []
    for s in range(n):
        X = blocks(C, s)
        M = np.einsum('abxy,cbxy->ac', X, X.conj()) + np.einsum('baxy,bcxy->ac', X.conj(), X)
        w, U = np.linalg.eigh(M); P = [np.outer(U[:, k], U[:, k].conj()) for k in range(3)]
        z = np.zeros(8 * n); z[8 * s:8 * s + 8] = to8((P[1] - P[2]) / np.sqrt(2)); V.append(z)
        for Hm in ((P[0] - P[1]) / np.sqrt(2), (P[0] + P[1] - 2 * P[2]) / np.sqrt(6)):
            z = np.zeros(8 * n); z[8 * s:8 * s + 8] = to8(Hm); W.append(z)
    V = np.array(V).T; W = np.array(W).T
    single = [restricted_min(Hess, V[:, [s]]) for s in range(n)]
    return restricted_min(Hess, V), single, restricted_min(Hess, W)


def identity_checks(C, r, rng):
    """max deviation of Tr F_s = q'_r(Tr_s C) - t sum q'(X_ab) and of F_s(I - E) = F_s(E) over copies / random e."""
    herm_basis = []
    for a in range(3):
        E = np.zeros((3, 3), complex); E[a, a] = 1; herm_basis.append(E)
    for a, b in itertools.combinations(range(3), 2):
        E = np.zeros((3, 3), complex); E[a, b] = E[b, a] = 1 / np.sqrt(2); herm_basis.append(E)
        E = np.zeros((3, 3), complex); E[a, b] = -1j / np.sqrt(2); E[b, a] = 1j / np.sqrt(2); herm_basis.append(E)
    dev_tr, dev_e = 0.0, 0.0
    for s in range(n):
        X = blocks(C, s); Ts = sum(X[a, a] for a in range(3))
        trF = sum(F_form(C, s, r, E) for E in herm_basis)
        rhs = q_np(Ts, n - 1) - r * np.vdot(Ts, Ts).real - t * sum(q_np(X[a, b], n - 1) for a in range(3) for b in range(3))
        dev_tr = max(dev_tr, abs(trF - rhs))
        e = rng.standard_normal(3) + 1j * rng.standard_normal(3); e /= np.linalg.norm(e); E = np.outer(e, e.conj())
        dev_e = max(dev_e, abs(F_form(C, s, r, np.eye(3) - E) - F_form(C, s, r, E)))
    return dev_tr, dev_e


if __name__ == "__main__":
    torch.manual_seed(SEED); rng = np.random.default_rng(SEED)
    for st in range(STARTS):
        p = torch.randn(8 * D, requires_grad=True)
        for bias in (30.0, 0.0):
            opt = torch.optim.LBFGS([p], lr=1, max_iter=2500, tolerance_grad=1e-15, tolerance_change=1e-18,
                                    history_size=60, line_search_fn='strong_wolfe')

            def closure():
                opt.zero_grad(); C, A, B = C_of(p); h = torch.zeros(8 * n, requires_grad=True)
                g = torch.autograd.grad(squeeze_R(C, h), h, create_graph=True)[0]
                val = (g * g).sum()
                if bias > 0:
                    val = val + bias * (R_t(C) - TARGET) ** 2 + 1e2 * torch.relu(TAU - thickness(A, B)) ** 2
                val.backward(); return val
            opt.step(closure)
        C, A, B = C_of(p.detach()); h = torch.zeros(8 * n, requires_grad=True)
        g = torch.autograd.grad(squeeze_R(C, h), h)[0]; gn = torch.linalg.norm(g).item()
        th = thickness(A, B).item(); r = R_t(C).item()
        head = f"start {st}: |grad| {gn:.1e}  R {r:+.5f}  thickness {th:.2e}"
        if not (gn < 1e-8 and th > 1e-3):
            print(head + "  -> not a thick critical point (slid to the boundary or not converged)", flush=True); continue
        H = torch.autograd.functional.hessian(lambda hh: squeeze_R(C, hh), torch.zeros(8 * n)); H = ((H + H.T) / 2).numpy()
        w, V = np.linalg.eigh(H); v = V[:, 0]
        bd = sum(v[8 * s:8 * s + 8] @ H[8 * s:8 * s + 8, 8 * s:8 * s + 8] @ v[8 * s:8 * s + 8] for s in range(n))
        Cn = C.numpy()
        if SAVE:
            np.save(f'{SAVE}/crit_n{n}_t{t}_R{TARGET}_s{SEED}_{st}.npy', Cn)
        d1, d2 = identity_checks(Cn, r, rng)
        per = [np.linalg.eigvalsh(H[8 * s:8 * s + 8, 8 * s:8 * s + 8])[0] for s in range(n)]
        jt, sg, md = top_two_test(Cn, H)
        print(head + f"\n    Hessian min {w[0]:+.4f}: single-copy part {bd:+.4f}, cross-copy part {w[0] - bd:+.4f}"
              f"\n    smallest eigenvalue of each single-copy block {np.round(per, 5)}"
              f"\n    top-two rebalance: joint {jt:+.4f}, each copy alone {np.round(sg, 4)};"
              f" all local-support-diagonal squeezes jointly {md:+.4f}"
              f"\n    rank-one test {rank_one_test(Cn, r, rng):+.5f}   plane test {plane_test(Cn, r, rng):+.5f}"
              f"\n    identity checks: trace formula {d1:.1e}, F(I-E) = F(E) {d2:.1e}", flush=True)
