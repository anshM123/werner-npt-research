"""Three-level certificate of Proposition 7 (single-excitation-range/README.md).

Hard-core operators on the co-range sector (F-convention, gauge sum_t A_t = I): beta_t F = sum_r <w_{t,r}|_t G^r,
a_r^dag G = sum_t |w_{t,r}>_t <0|_t G, w_{t,r} = U_t e_r.  For every real mu
  2Kt = sum_S F_S^*(2I + (2-2mu) Abar_H)F_S - sum_rho tau^*[(2-2mu)I - (1-2mu)J]tau + 2 mu sum_{S'} Pi_{S'},
Pi_{S'} = part of <Theta, S Theta> = sum_{r,r'} <a_{r'}^dag G^r, a_r^dag G^{r'}> landing on excited set S'.
certificate(M, c) maximises a common slack over theta_Q >= 0, letter-local metrics m_{S,t}, child metrics k_{S,u} and mu,
and re-checks every matrix inequality of the solution with numpy eigenvalues.  Riesz vector of tau in this convention:
vec(U_t) (index (letter, frame)).
"""
import numpy as np, cvxpy as cp, itertools, warnings
warnings.filterwarnings('ignore')

class Model:
    """Restricted heavy family: copy 0 heavy (2 letters, U_0 = diag(sqrt(1-d), 1)); copies 1..N light with one letter,
    U_t = sqrt(d/N)|1><e_1| (the invariant subspace in which every excited light copy carries letter 1)."""
    def __init__(self, N, d):
        self.N, self.d, self.n = N, d, N + 1
        n = self.n
        self.U = [np.array([[np.sqrt(1 - d), 0], [0, 1.0]], complex)] + [np.array([[np.sqrt(d / N), 0]], complex) for _ in range(N)]
        self.L = [2] + [1] * N                      # letters per copy
        self.A = [u.conj().T @ u for u in self.U]
    def strings(self, k):
        out = []
        for exc in itertools.combinations(range(self.n), k):
            for lets in itertools.product(*[range(1, self.L[t] + 1) for t in exc]):
                s = [0] * self.n
                for t, z in zip(exc, lets): s[t] = z
                out.append(tuple(s))
        return out
    def level(self, k):
        S = self.strings(k); pos = {s: i for i, s in enumerate(S)}; return S, pos
    def ann(self, c):
        """beta_t: V (x) H_c -> H_{c-1} as matrices (index r*m + j)."""
        S, pS = self.level(c); R, pR = self.level(c - 1); m = len(S); B = []
        for t in range(self.n):
            M = np.zeros((len(R), 2 * m), complex)
            for j, s in enumerate(S):
                if s[t] == 0: continue
                rho = list(s); rho[t] = 0; i = pR[tuple(rho)]
                for r in range(2): M[i, r * m + j] += np.conj(self.U[t][s[t] - 1, r])
            B.append(M)
        return B
    def cre(self, c):
        """a_r^dag: H_c -> H_{c+1} (frame component r of F mapped by creation with w_{t,r}), list over r of matrices m_{c+1} x m_c."""
        S, pS = self.level(c); P, pP = self.level(c + 1); out = []
        for r in range(2):
            M = np.zeros((len(P), len(S)), complex)
            for j, s in enumerate(S):
                for t in range(self.n):
                    if s[t] != 0: continue
                    for z in range(1, self.L[t] + 1):
                        sp = list(s); sp[t] = z; M[pP[tuple(sp)], j] += self.U[t][z - 1, r]
            out.append(M)
        return out
    def ops(self, c):
        S, pS = self.level(c); m = len(S); n = self.n; I = np.eye(2 * m)
        B = self.ann(c); T = sum(b.conj().T @ b for b in B); al = sum(B); AA = al.conj().T @ al
        Vac = np.zeros((2 * m, 2 * m), complex)
        for j, s in enumerate(S):
            Ab = sum((self.A[t].conj() for t in range(n) if s[t] == 0), np.zeros((2, 2)))
            for r in range(2):
                for rp in range(2): Vac[r * m + j, rp * m + j] += Ab[r, rp]
        K = I + Vac - T + 0.5 * AA
        # Theta pairing: <Theta, S Theta> = sum_{r r'} <a_{r'}^dag G^r, a_r^dag G^{r'}>;  G^r = F[r*m:(r+1)*m]
        Cr = self.cre(c); P, pP = self.level(c + 1)
        def E(r):   # selector G^r = E(r) F
            M = np.zeros((m, 2 * m)); M[:, r * m:(r + 1) * m] = np.eye(m); return M
        ThS = sum((Cr[rp] @ E(r)).conj().T @ (Cr[r] @ E(rp)) for r in range(2) for rp in range(2))
        # split by child pattern S' (excited set of level c+1 strings)
        supp = lambda s: frozenset(t for t in range(n) if s[t])
        children = sorted(set(supp(s) for s in P), key=lambda x: sorted(x))
        Pi = {}
        for Sp in children:
            proj = np.diag([1.0 if supp(s) == Sp else 0.0 for s in P])
            Pi[Sp] = sum((proj @ Cr[rp] @ E(r)).conj().T @ (proj @ Cr[r] @ E(rp)) for r in range(2) for rp in range(2))
        return dict(S=S, m=m, B=B, T=T, AA=AA, Vac=Vac, K=(K + K.conj().T) / 2, ThS=ThS, Pi=Pi, I=I)



class GeneralModel(Model):
    """Full single-excitation model for an arbitrary list of range blocks U_t (2 x 2, gauge sum U_t^* U_t = I)."""
    def __init__(self, U):
        self.n = len(U); self.N = self.n - 1; self.d = None
        self.U = [np.asarray(u, complex) for u in U]; self.L = [u.shape[0] for u in self.U]; self.A = [u.conj().T @ u for u in self.U]

def gauge(X):
    S = sum(x.conj().T @ x for x in X); w, V = np.linalg.eigh(S); return [x @ (V @ np.diag(w ** -0.5) @ V.conj().T) for x in X]

def certificate(M, c, three_level=True, solver="CLARABEL", fix_mu=None, no_theta=False, real=False, solver_opts=None):
    n = M.n; o = M.ops(c); S_list = o['S']; m = o['m']
    if real:
        o = {k_: (np.real(v) if isinstance(v, np.ndarray) else ({kk: np.real(vv) for kk, vv in v.items()} if isinstance(v, dict) else v)) for k_, v in o.items()}
    supp = lambda s: frozenset(t for t in range(n) if s[t])
    pats = sorted(set(supp(s) for s in S_list), key=lambda x: sorted(x))
    idx = {Sx: [r * m + j for r in range(2) for j, s in enumerate(S_list) if supp(s) == Sx] for Sx in pats}   # (r, string) order
    strs = {Sx: [s for s in S_list if supp(s) == Sx] for Sx in pats}
    ys = [np.array([M.U[t][z, r] for z in range(M.L[t]) for r in range(2)]) for t in range(n)]
    if real: ys = [np.real(y) for y in ys]          # F-convention Riesz vector vec(U_t), index (z, r)
    def embed_local(Sx, t):
        """basis matrices E[(a,b)] mapping a (letter_t, r) Hermitian block into the pattern block of Sx (order (r, string))."""
        st = strs[Sx]; k = len(st); D = 2 * k; dl = 2 * M.L[t]; E = {}
        for a in range(dl):
            for b in range(dl):
                za, ra = divmod(a, 2); zb, rb = divmod(b, 2); Mx = np.zeros((D, D))
                for i, si in enumerate(st):
                    for jx, sj in enumerate(st):
                        if si[t] - 1 != za or sj[t] - 1 != zb: continue
                        if any(si[u] != sj[u] for u in range(n) if u != t): continue
                        Mx[ra * k + i, rb * k + jx] = 1.0
                E[a, b] = Mx
        return E
    Qs = list(itertools.combinations(range(n), c + 1))
    th = {Q: (0.0 if no_theta else cp.Variable(nonneg=True)) for Q in Qs}; mu = (fix_mu if fix_mu is not None else cp.Variable()) if three_level else 0.0; lam = cp.Variable()
    mv = {}; kv = {}; cons = []; exprs = []
    for Sx in pats:
        for t in sorted(Sx): mv[Sx, t] = cp.Variable((2 * M.L[t], 2 * M.L[t]), symmetric=True) if real else cp.Variable((2 * M.L[t], 2 * M.L[t]), hermitian=True)
        if three_level:
            for u in range(n):
                if u not in Sx: kv[Sx, u] = cp.Variable((len(idx[Sx]), len(idx[Sx])), symmetric=True) if real else cp.Variable((len(idx[Sx]), len(idx[Sx])), hermitian=True)
    # (P)
    for Sx in pats:
        D = len(idx[Sx]); k = D // 2; H = [u for u in range(n) if u not in Sx]
        AbH = sum((M.A[u].conj() for u in H), np.zeros((2, 2)))
        def frame_op(X):   # X on V -> pattern block (order (r, string))
            return np.kron(X, np.eye(k))
        budget = 2 * np.eye(D) + (2 - 2 * mu) * frame_op(AbH) if three_level else 2 * np.eye(D) + 2 * frame_op(AbH)
        tot = 0
        for t in sorted(Sx):
            E = embed_local(Sx, t); tot = tot + sum(mv[Sx, t][a, b] * E[a, b] for (a, b) in E)
        for u in H:
            Q = tuple(sorted(tuple(Sx) + (u,))); tot = tot + 2 * th[Q] * frame_op(sum(M.A[x].conj() for x in Q))
            if three_level: tot = tot + kv[Sx, u]
        cons.append(budget - tot - lam * np.eye(D) >> 0); exprs.append(budget - tot)
    # (L)
    R_sets = sorted(set(supp(s) for s in M.strings(c - 1)), key=lambda x: sorted(x))
    for R in R_sets:
        free = [t for t in range(n) if t not in R]; kf = len(free)
        coefI = (2 - 2 * mu) if three_level else 2.0; coefJ = (1 - 2 * mu) if three_level else 1.0
        blocks = [[None] * kf for _ in range(kf)]
        for i in range(kf):
            for j in range(kf):
                th_ij = th[tuple(sorted(tuple(R) + (free[i], free[j])))] if i != j else None
                if i == j: Nij = coefI - coefJ - sum(th[tuple(sorted(tuple(R) + (free[i], free[l])))] for l in range(kf) if l != i)
                else: Nij = -coefJ + th_ij
                blk = Nij * np.outer(ys[free[i]], ys[free[j]].conj())
                blocks[i][j] = (mv[frozenset(R) | {free[i]}, free[i]] if i == j else 0) - blk
        D = sum(2 * M.L[t] for t in free)
        Lb = cp.bmat(blocks); cons.append(Lb - lam * np.eye(D) >> 0); exprs.append(Lb)
    # (C)
    if three_level:
        for Sp, Pi in o['Pi'].items():
            members = [u for u in sorted(Sp) if (Sp - {u}) in idx]
            ind = sum((idx[Sp - {u}] for u in members), []); Pis = Pi[np.ix_(ind, ind)]
            sizes = [len(idx[Sp - {u}]) for u in members]
            blocks = [[None] * len(members) for _ in members]
            for i, u in enumerate(members):
                for j, w in enumerate(members):
                    oi, oj = sum(sizes[:i]), sum(sizes[:j])
                    blk = 2 * mu * Pis[oi:oi + sizes[i], oj:oj + sizes[j]]
                    blocks[i][j] = blk + (kv[Sp - {u}, u] if i == j else 0)
            Cb = cp.bmat(blocks); cons.append(Cb - lam * np.eye(len(ind)) >> 0); exprs.append(Cb)
    prob = cp.Problem(cp.Maximize(lam), cons); prob.solve(solver=solver, **(solver_opts or {}))
    chk = min(np.linalg.eigvalsh((lambda X: (X + X.conj().T) / 2)(np.asarray(E.value if hasattr(E, 'value') else E)))[0] for E in exprs)
    muv = (mu.value if hasattr(mu, 'value') else mu) if three_level else 0.0
    thmin = min([0.0] + [th[Q].value for Q in Qs if hasattr(th[Q], 'value')]); thmax = max([0.0] + [th[Q].value for Q in Qs if hasattr(th[Q], 'value')])
    return lam.value, muv, prob.status, chk, thmin, thmax

