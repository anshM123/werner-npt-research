"""Relative anchors (README, Lemma 11 and Conjecture 12).

For p in H_1 (single excitations over |0...0>): D_c = compression to the column space H_c of D_p = V_p^* M V_p,
R_c = M V_p D_c^{-1} V_p^* M, and Y_c(p) = compression of (M - R_c)^Gamma to H_1 (x) H_c, i.e.
Y_c(p)[(a, y'), (b, y)] = W[(b, y'), (a, y)] with W = M - R_c, rows a, b in H_1 and columns y, y' in H_c.
M = M_1^{(x)n}, M_1[(a,y),(b,y')] = d_ab d_yy' - 1/2 d_ay d_by'.
Run as a script:
(a) Lemma 11's Schur step: for X = |p><v| + |u><w| with u in H_1 and v, w in H_c, min over v of q_n(X) equals
    (conj u (x) conj w)^* Y_c(p) (conj u (x) conj w)   (identity, random data);
(b) evidence for Conjecture 12: Y_c(sum_t w_t |1_t>) >= 0 at random weights, n = 3, 4, 5, every c, and the kernel
    contains conj(p) (x) H_c.  This is evidence, not a proof.
(c) charge conservation: the full R_p agrees with R_c on H_1 (x) H_c (n = 3, every c), so Y_c(p) is a compression
    of (M - R_p)^Gamma; with Proposition 10 this proves Conjecture 12 when at most two weights are nonzero.
adversarial_relative.py minimises lambda_min(Y_c) over the weights (evidence)."""
import numpy as np, itertools

def strings(n, k):
    out = []
    for S in itertools.combinations(range(n), k):
        for lets in itertools.product([1, 2], repeat=k):
            s = [0] * n
            for t, z in zip(S, lets): s[t] = z
            out.append(tuple(s))
    return out

def Mel(a, y, b, yp):
    v = 1.0
    for t in range(len(a)):
        v *= (1.0 if (a[t] == b[t] and y[t] == yp[t]) else 0.0) - (0.5 if (a[t] == y[t] and b[t] == yp[t]) else 0.0)
        if v == 0: return 0.0
    return v

class Relative:
    """Weight-independent tensors for given n, c; Y(pvec) for p given by its 2n coordinates on H_1."""
    def __init__(self, n, c):
        self.n, self.c = n, c
        self.H1 = strings(n, 1); self.Hc = strings(n, c)
        self.idx = [(a, y) for a in self.H1 for y in self.Hc]; pos = {k: i for i, k in enumerate(self.idx)}
        H1, Hc = self.H1, self.Hc
        self.Dt = np.array([[[[Mel(a, x, b, xp) for b in H1] for a in H1] for xp in Hc] for x in Hc])
        self.Gt = np.array([[[Mel(a, y, b, x) for b in H1] for x in Hc] for (a, y) in self.idx])
        self.Mc = np.array([[Mel(a, y, b, yp) for (b, yp) in self.idx] for (a, y) in self.idx])
        self.perm = np.array([[pos[(b, yp)] for (b, y) in self.idx] for (a, yp) in self.idx])
        self.perm2 = np.array([[pos[(a, y)] for (b, y) in self.idx] for (a, yp) in self.idx])
    def parts(self, p):
        D = np.einsum('ijab,a,b->ij', self.Dt, p.conj(), p); G = self.Gt @ p
        W = self.Mc - G @ np.linalg.solve(D, G.conj().T)
        return D, G, W
    def Y(self, p):
        D, G, W = self.parts(p); Y = W[self.perm, self.perm2]; return (Y + Y.conj().T) / 2, D
    def pvec(self, w):
        p = np.zeros(len(self.H1), complex)
        for i, a in enumerate(self.H1):
            t = [s for s in range(self.n) if a[s]][0]
            if a[t] == 1: p[i] = w[t]
        return p

if __name__ == '__main__':
    rng = np.random.default_rng(5); ok = True
    cz = lambda m: rng.normal(size=m) + 1j * rng.normal(size=m)
    # (a) the Schur step
    for n, c in [(3, 2), (4, 2), (4, 3)]:
        R = Relative(n, c); m1, mc = len(R.H1), len(R.Hc)
        p, u, w = cz(m1), cz(m1), cz(mc)
        D, G, W = R.parts(p)
        # q(X) as a quadratic form in conj(v): vec X = p (x) conj v + u (x) conj w
        Pm = np.kron(p.reshape(-1, 1), np.eye(mc)); b = np.kron(u, w.conj())
        H = Pm.conj().T @ R.Mc @ Pm; k = Pm.conj().T @ R.Mc @ b
        qmin = (b.conj() @ R.Mc @ b - k.conj() @ np.linalg.solve(H, k)).real
        Yv, _ = R.Y(p); z = np.kron(u.conj(), w.conj())
        err = abs(qmin - (z.conj() @ Yv @ z).real) / abs(qmin)
        print(f"(a) n={n} c={c}: min_v q_n(X) vs (u (x) w)^* Y (u (x) w): rel. error {err:.1e}"); ok &= err < 1e-10
    # (b) evidence for Conjecture 12
    for n in [3, 4, 5]:
        for c in range(1, n + 1):
            R = Relative(n, c); worst = 9
            for k in range(5):
                wts = rng.random(n) + 0.02
                if k == 3: wts = np.array([1.0] + [0.03] * (n - 1))
                if k == 4: wts = np.array([1.0, 1.0] + [0.05] * (n - 2))
                p = R.pvec(wts); Yv, D = R.Y(p)
                ev = np.linalg.eigvalsh(Yv); worst = min(worst, ev[0])
                ker = np.kron(p.conj().reshape(-1, 1), np.eye(len(R.Hc)))       # conj(p) (x) H_c
                ok &= np.linalg.norm(Yv @ ker) < 1e-9 * max(1, np.linalg.norm(Yv))
            print(f"(b) n={n} c={c}: min eigenvalue of Y_c over 5 weight vectors {worst:+.2e}"); ok &= worst > -1e-12
    # (c) charge conservation: the full R_p agrees with R_c on H_1 (x) H_c, so Y_c(p) is a compression of (M - R_p)^Gamma
    n = 3; Sfull = list(itertools.product(range(3), repeat=n)); N = len(Sfull)
    idxf = [(a, y) for a in Sfull for y in Sfull]; ipf = {k: i for i, k in enumerate(idxf)}
    Mf = np.array([[Mel(a, y, b, yp) for (b, yp) in idxf] for (a, y) in idxf])
    for c in range(1, n + 1):
        R = Relative(n, c); p = cz(len(R.H1))
        Vp = np.zeros((N * N, N), complex)
        for i, a in enumerate(R.H1):
            for j, x in enumerate(Sfull): Vp[ipf[(a, x)], j] += p[i]
        Dp = Vp.conj().T @ Mf @ Vp; Rp = Mf @ Vp @ np.linalg.solve(Dp, Vp.conj().T @ Mf)
        D, G, W = R.parts(p); Rc = R.Mc - W
        sel = [ipf[k] for k in R.idx]
        err = np.abs(Rp[np.ix_(sel, sel)] - Rc).max()
        print(f"(c) n=3 c={c}: max |R_p - R_c| on H_1 (x) H_c = {err:.1e}"); ok &= err < 1e-10
    print("PASS" if ok else "FAIL")
