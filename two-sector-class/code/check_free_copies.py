"""f free copies (copies 0..f-1, full qutrits) + sector copies (A_s perp B_s, dims (2,1),(1,2),(1,1)).
Checks (i) identity 2^m q_n(C) = sum_mu q_F(X_mu) + sum_nu q_F(Y_nu) + 2^m (qt_F(P) + qt_F(Q)) + 2(-1)^m Re B_F(X_0, Y_0);
(ii) Lemma R_F: for X = x y^*, sum_mu q_F(X_mu) >= 2 sum_{ab} q_F(x_a y_b^*) - q_F(X_0);
(iii) reduction: 2^m q_n(C) >= 2 q_{F+big}(xy^* (+) e^{i phi} x'y'^*) + 2^m(qt_F P + qt_F Q) (phi worst);
(iv) Theorem A_f: adversarial min of q_n/|C|^2 with rank X, Y <= 1, rank P, Q <= 2;
(v) Corollary B_f: range spanned by uA in H_F (x) A and uB in H_F (x) B, co-range arbitrary.
usage: python check_free_copies.py f m trials seed  (f = number of free copies, m = number of sector copies)"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64); cplx = torch.complex128
f, m, TR, SEED = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); n = f + m; D = 3 ** n; DF = 3 ** f
torch.manual_seed(SEED); rng = np.random.default_rng(SEED)
w9 = torch.eye(3, dtype=cplx).reshape(9) / np.sqrt(3)
def make_Phi(k):
    Dk = 3 ** k
    def Phi(C):
        if k == 0: return C
        X = C.reshape((3,) * (2 * k)).permute([i for s in range(k) for i in (s, k + s)]).reshape((9,) * k)
        for s in range(k):
            tr = torch.tensordot(w9.conj(), X, dims=([0], [s])); X = X - 1.5 * torch.moveaxis(torch.tensordot(w9, tr, dims=0), 0, s)
        return X.reshape([3, 3] * k).permute([2 * s for s in range(k)] + [2 * s + 1 for s in range(k)]).reshape(Dk, Dk)
    return Phi
Phi_n, Phi_F = make_Phi(n), make_Phi(f)
q = lambda C: torch.vdot(C.reshape(-1), Phi_n(C).reshape(-1)).real
qF = lambda Z: torch.vdot(Z.reshape(-1), Phi_F(Z).reshape(-1)).real
BF = lambda Z, W: torch.vdot(Z.reshape(-1), Phi_F(W).reshape(-1))
def kron_all(ms):
    out = torch.ones(1, 1, dtype=cplx)
    for M in ms: out = torch.kron(out, M)
    return out
PA = [torch.tensor(x, dtype=cplx) for x in ([[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]])]
def strings(dims):
    two = [k for k, d in enumerate(dims) if d == 2]; out = []
    for mus in itertools.product(range(3), repeat=len(two)):
        ms, it = [], iter(mus)
        for d in dims: ms.append(PA[next(it)] if d == 2 else torch.eye(1, dtype=cplx))
        out.append(kron_all(ms))
    return out
def hat(X, Sg, dA):
    X4 = X.reshape(DF, dA, DF, dA); return torch.einsum('ag,igja->ij', Sg, X4)
def qtF(P, dA, dB):  # sum over sector blocks of q_F
    P4 = P.reshape(DF, dA, DF, dB); return sum(qF(P4[:, a, :, b]) for a in range(dA) for b in range(dB))
def q_big(x, y, xp, yp, phi, dA, dB):  # q_{F + big} of x y^* (+) e^{i phi} x' y'^*, big copy = A (+) B with Werner form
    X4 = (x @ y.conj().T).reshape(DF, dA, DF, dA); Y4 = (xp @ yp.conj().T).reshape(DF, dB, DF, dB)
    tot = sum(qF(X4[:, a, :, b]) for a in range(dA) for b in range(dA)) + sum(qF(Y4[:, a, :, b]) for a in range(dB) for b in range(dB))
    T = torch.einsum('iaja->ij', X4) + np.exp(1j * phi) * torch.einsum('iaja->ij', Y4)
    return tot - 0.5 * qF(T)
idmax, Rmin, red_min, thA, corB = 0.0, 1e9, 1e9, 1e9, 1e9
for trial in range(TR):
    dims = [[(2, 1), (1, 2), (1, 1)][rng.integers(3)] for _ in range(m)]
    Us = [torch.tensor(np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))[0]) for _ in range(m)]
    EA = kron_all([torch.eye(DF, dtype=cplx)] + [Us[s][:, :dims[s][0]] for s in range(m)])
    EB = kron_all([torch.eye(DF, dtype=cplx)] + [Us[s][:, dims[s][0]:sum(dims[s])] for s in range(m)])
    dA, dB = EA.shape[1] // DF, EB.shape[1] // DF
    SA, SB = strings([d[0] for d in dims]), strings([d[1] for d in dims])
    IA, IB = torch.eye(dA, dtype=cplx), torch.eye(dB, dtype=cplx)
    X = torch.randn(DF * dA, DF * dA, dtype=cplx); Y = torch.randn(DF * dB, DF * dB, dtype=cplx)
    P = torch.randn(DF * dA, DF * dB, dtype=cplx); Qm = torch.randn(DF * dB, DF * dA, dtype=cplx)
    C = EA @ X @ EA.conj().T + EB @ Y @ EB.conj().T + EA @ P @ EB.conj().T + EB @ Qm @ EA.conj().T
    def qtF_any(P, d1, d2):
        P4 = P.reshape(DF, d1, DF, d2); return sum(qF(P4[:, a, :, b]) for a in range(d1) for b in range(d2))
    rhs = sum(qF(hat(X, S, dA)) for S in SA) + sum(qF(hat(Y, S, dB)) for S in SB) + 2 ** m * (qtF_any(P, dA, dB) + qtF_any(Qm, dB, dA)) \
        + 2 * (-1) ** m * BF(hat(X, IA, dA), hat(Y, IB, dB)).real
    idmax = max(idmax, (abs(2 ** m * q(C) - rhs) / torch.linalg.norm(C) ** 2).item())
    x = torch.randn(DF * dA, 1, dtype=cplx); y = torch.randn(DF * dA, 1, dtype=cplx)
    lhsR = sum(qF(hat(x @ y.conj().T, S, dA)) for S in SA)
    X4 = (x @ y.conj().T).reshape(DF, dA, DF, dA)
    rhsR = 2 * sum(qF(X4[:, a, :, b]) for a in range(dA) for b in range(dA)) - qF(torch.einsum('iaja->ij', X4))
    Rmin = min(Rmin, ((lhsR - rhsR) / (torch.linalg.norm(x) * torch.linalg.norm(y)) ** 2).item())
    def memberA(p):
        k = 0
        def take(r, c):
            nonlocal k
            v = torch.complex(p[k:k + r * c], p[k + r * c:k + 2 * r * c]).reshape(r, c); k += 2 * r * c; return v
        x, y, xp, yp = take(DF * dA, 1), take(DF * dA, 1), take(DF * dB, 1), take(DF * dB, 1)
        U1, V1, U2, V2 = take(DF * dA, 2), take(DF * dB, 2), take(DF * dB, 2), take(DF * dA, 2)
        C = EA @ (x @ y.conj().T) @ EA.conj().T + EB @ (xp @ yp.conj().T) @ EB.conj().T \
            + EA @ (U1 @ V1.conj().T) @ EB.conj().T + EB @ (U2 @ V2.conj().T) @ EA.conj().T
        nrm = torch.linalg.norm(C); return C / nrm, (x, y, xp, yp, U1 @ V1.conj().T, U2 @ V2.conj().T, nrm)
    npA = 2 * (2 * DF * dA + 2 * DF * dB + 2 * DF * dA * 2 + 2 * DF * dB * 2)
    for st in range(2):
        p = torch.randn(npA, requires_grad=True)
        opt = torch.optim.LBFGS([p], lr=1, max_iter=1000, tolerance_grad=1e-14, tolerance_change=1e-18, history_size=50, line_search_fn='strong_wolfe')
        def closure():
            opt.zero_grad(); v = q(memberA(p)[0]); v.backward(); return v
        opt.step(closure)
        with torch.no_grad():
            Cc, (x, y, xp, yp, Pp, Qp, nrm) = memberA(p); thA = min(thA, q(Cc).item())
            # reduction inequality at the minimiser (worst phase over a grid)
            worst = min(q_big(x, y, xp, yp, ph, dA, dB).item() for ph in np.linspace(0, 2 * np.pi, 73))
            red = (2 ** m * q(Cc * nrm) - 2 * worst - 2 ** m * (qtF_any(Pp, dA, dB) + qtF_any(Qp, dB, dA))).item() / nrm.item() ** 2
            red_min = min(red_min, red)
    def memberC(p):
        uA = EA @ torch.complex(p[:DF * dA], p[DF * dA:2 * DF * dA]); o = 2 * DF * dA
        uB = EB @ torch.complex(p[o:o + DF * dB], p[o + DF * dB:o + 2 * DF * dB]); o += 2 * DF * dB
        fv = torch.complex(p[o:o + 2 * D], p[o + 2 * D:o + 4 * D]).reshape(D, 2)
        C = torch.stack([uA, uB], 1) @ fv.conj().T; return C / torch.linalg.norm(C)
    npC = 2 * DF * dA + 2 * DF * dB + 4 * D
    for st in range(2):
        p = torch.randn(npC, requires_grad=True)
        opt = torch.optim.LBFGS([p], lr=1, max_iter=1000, tolerance_grad=1e-14, tolerance_change=1e-18, history_size=50, line_search_fn='strong_wolfe')
        def closure():
            opt.zero_grad(); v = q(memberC(p)); v.backward(); return v
        opt.step(closure)
        with torch.no_grad(): corB = min(corB, q(memberC(p)).item())
    print(f"trial {trial} dims {dims}: running: id err {idmax:.1e}, LemmaR gap min {Rmin:+.2e}, reduction slack min {red_min:+.2e}, ThmA min {thA:+.2e}, CorB min {corB:+.2e}", flush=True)
