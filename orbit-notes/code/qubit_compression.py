"""Product-qubit compressions (orbit-notes/NOTES.md, section "Product-qubit compressions").

For a rank-two C on (C^3)^{(x)n} and leak vectors l_1..l_n in C^3, put Q_s = l_s^perp, P = P_{Q_1} (x) ... (x) P_{Q_n}
and C_Q = P C P.  On operators supported on Q_1 (x) ... (x) Q_n the endpoint form is positive semidefinite, so
R(C_Q) >= 0.  Conjecture CC_n: some choice of Q gives R(C_Q) <= R(C) (this would imply the endpoint at n copies).
R(X) = q_n(1/2, X)/||X||_F^2 throughout.

Subcommands (all floating point; inner problem = minimise R(C_Q) over the leak vectors by L-BFGS from random starts):
  points FILE...                 min_Q R(C_Q) against R(C) at saved operators
  adversary n starts seed steps  outer gradient ascent on G(C) = min_Q R(C_Q) - R(C) (a positive G refutes CC_n);
                                 the inner minimiser is refreshed every few steps and re-checked at the end
  nearzero n seed                C = (U0 + e E1)(V0 + e E2)^* around known zeros, e = 0.1, 0.03, 0.01
  onecopy n starts seed VARIANT  the same ascent for one-copy compressions only:
                                 fixed  = two-sided compression of copy 0,
                                 rows   = rows of copy 0 only,
                                 any    = two-sided compression of the best single copy,
                                 stage2 = copy 0 already a qubit, compress copy 1 (n = 2)
  verify FILE restarts           min_Q R(C_Q) with all eigenvector starts (marginals of CC*+C*C, CC*, C*C) plus
                                 random restarts; used for data/compression_counterexample_n3.npy
  independent restarts FILE...   rows and columns compressed to independent planes, R(P_A C P_B) (still q >= 0)
  orbit n trials seed            minimise R(G C G^*) over the squeeze orbit from G = I (where does descent go?)
See also qubit_compression_check.py (independent numpy/scipy check of the counterexample).
"""
import sys, numpy as np, torch
torch.set_default_dtype(torch.float64)
CPX = torch.complex128
T = 0.5


def phi(C, n):
    d = 3 ** n; X = C.reshape((3,) * (2 * n)); eye = torch.eye(3, dtype=C.dtype)
    for s in range(n):
        tr = torch.diagonal(X, dim1=s, dim2=n + s).sum(-1)
        rest = [a for a in range(2 * n) if a not in (s, n + s)]; dest = rest + [s, n + s]
        X = X - T * (tr.unsqueeze(-1).unsqueeze(-1) * eye).permute(*[dest.index(a) for a in range(2 * n)])
    return X.reshape(d, d)


def R(C, n):
    return torch.vdot(C.reshape(-1), phi(C, n).reshape(-1)).real / torch.vdot(C.reshape(-1), C.reshape(-1)).real


def kron_list(ms):
    out = torch.ones((1, 1), dtype=CPX)
    for m in ms: out = torch.kron(out, m)
    return out


def qproj(l):
    l = l / torch.linalg.norm(l); return torch.eye(3, dtype=CPX) - torch.outer(l, l.conj())


class Compressor:
    """sites: copies that are compressed; mode 'two' (P C P) or 'rows' (P C)."""
    def __init__(self, n, sites, mode, gen):
        self.n, self.sites, self.mode, self.gen = n, list(sites), mode, gen

    def apply(self, C, L):           # L: (len(sites), 3) complex
        ms = [torch.eye(3, dtype=CPX) for _ in range(self.n)]
        for k, s in enumerate(self.sites): ms[s] = qproj(L[k])
        P = kron_list(ms)
        return P @ C @ P if self.mode == 'two' else P @ C

    def inner(self, C, starts, init=None):
        best = (np.inf, None); Cd = C.detach(); m = len(self.sites)
        inits = ([init] if init is not None else []) + [
            torch.complex(torch.randn(m, 3, generator=self.gen), torch.randn(m, 3, generator=self.gen)) for _ in range(starts)]
        for L0 in inits:
            p = torch.cat([L0.real.reshape(-1), L0.imag.reshape(-1)]).clone().requires_grad_(True)
            opt = torch.optim.LBFGS([p], max_iter=300, tolerance_grad=1e-12, tolerance_change=1e-15, line_search_fn='strong_wolfe')

            def cl():
                opt.zero_grad(); CQ = self.apply(Cd, torch.complex(p[:3 * m].reshape(m, 3), p[3 * m:].reshape(m, 3)))
                nq = torch.vdot(CQ.reshape(-1), CQ.reshape(-1)).real
                v = R(CQ, self.n) + 1e3 * torch.relu(1e-4 - nq); v.backward(); return v
            opt.step(cl)
            with torch.no_grad():
                L = torch.complex(p[:3 * m].reshape(m, 3), p[3 * m:].reshape(m, 3)); CQ = self.apply(Cd, L)
                if torch.vdot(CQ.reshape(-1), CQ.reshape(-1)).real.item() > 1e-5:
                    v = R(CQ, self.n).item()
                    if v < best[0]: best = (v, L.clone())
        return best


def rank_two(p, d):
    U = torch.complex(p[:2 * d].reshape(d, 2), p[2 * d:4 * d].reshape(d, 2))
    V = torch.complex(p[4 * d:6 * d].reshape(d, 2), p[6 * d:].reshape(d, 2))
    C = U @ V.conj().T; return C / torch.linalg.norm(C)


def ascent(n, comps, starts, gen, steps=150, lr=0.02, every=10, pre=None, label=''):
    """maximise min over compressors c in comps of [min_Q R(c(C)) ] - R(C)."""
    d = 3 ** n; worst = -np.inf
    for st in range(starts):
        p = torch.randn(8 * d, generator=gen).requires_grad_(True)
        build = (lambda p: rank_two(p, d)) if pre is None else (lambda p: pre(rank_two(p, d)))
        Ls = [None] * len(comps); opt = torch.optim.Adam([p], lr=lr)
        for it in range(steps):
            C = build(p)
            if it % every == 0:
                for k, c in enumerate(comps):
                    Ls[k] = c.inner(C, starts=(10 if it % (5 * every) == 0 else 3), init=Ls[k])[1]
            opt.zero_grad()
            vals = torch.stack([R(c.apply(C, Ls[k].detach()), n) for k, c in enumerate(comps)])
            G = -0.002 * torch.logsumexp(-vals / 0.002, 0) - R(C, n); (-G).backward(); opt.step()
        C = build(p).detach()
        v = min(c.inner(C, starts=40, init=Ls[k])[0] for k, c in enumerate(comps)); g = v - R(C, n).item()
        worst = max(worst, g)
        print(f"{label} start {st}: R(C) = {R(C, n).item():+.5f}  min R(compressed) = {v:+.5f}  gap = {g:+.5f}", flush=True)
    print(f"{label} n = {n}: largest gap found = {worst:+.6f}  (positive refutes the statement)", flush=True)


def main():
    cmd = sys.argv[1]
    if cmd == 'points':
        gen = torch.Generator().manual_seed(3)
        for f in sys.argv[2:]:
            C0 = np.load(f); n = round(np.log(C0.shape[0]) / np.log(3)); C = torch.tensor(C0 / np.linalg.norm(C0), dtype=CPX)
            v, _ = Compressor(n, range(n), 'two', gen).inner(C, starts=30)
            print(f"{f.replace(chr(92), '/').split('/')[-1]:30s} n={n} R(C) = {R(C, n).item():+.5f}  min_Q R(C_Q) = {v:+.5f}"
                  f"  gap = {v - R(C, n).item():+.5f}", flush=True)
    elif cmd == 'adversary':
        n, starts, seed, steps = (int(x) for x in sys.argv[2:6])
        lr, every = (0.02, 10) if steps <= 200 else (0.006, 6)
        gen = torch.Generator().manual_seed(seed)
        ascent(n, [Compressor(n, range(n), 'two', gen)], starts, gen, steps, lr, every, label='CC')
    elif cmd == 'nearzero':
        n, seed = int(sys.argv[2]), int(sys.argv[3]); gen = torch.Generator().manual_seed(seed)
        rng = np.random.default_rng(seed); d = 3 ** n
        e = lambda *idx: np.eye(d)[np.ravel_multi_index(idx, (3,) * n)]
        zeros = {'P_01 (x) |0><1|': (np.stack([e(0, 0), e(1, 0)], 1), np.stack([e(0, 1), e(1, 1)], 1)),
                 '|00><00| - |11><11|': (np.stack([e(0, 0), e(1, 1)], 1), np.stack([e(0, 0), -e(1, 1)], 1)),
                 '|10><10| - |01><01|': (np.stack([e(1, 0), e(0, 1)], 1), np.stack([e(1, 0), -e(0, 1)], 1))}
        comp = Compressor(n, range(n), 'two', gen)
        for name, (U0, V0) in zeros.items():
            print(f"zero {name}: R = {R(torch.tensor(U0 @ V0.T, dtype=CPX), n).item():+.1e}")
            for trial in range(3):
                E1 = rng.standard_normal((d, 2)) + 1j * rng.standard_normal((d, 2)); E2 = rng.standard_normal((d, 2)) + 1j * rng.standard_normal((d, 2))
                for eps in (0.1, 0.03, 0.01):
                    C = torch.tensor((U0 + eps * E1) @ (V0 + eps * E2).conj().T); C = C / torch.linalg.norm(C)
                    v, _ = comp.inner(C, starts=40); r = R(C, n).item()
                    print(f"   trial {trial} eps {eps:5.2f}: R(C) = {r:+.3e}  min_Q R(C_Q) = {v:+.3e}  ratio = {v / r:.3f}", flush=True)
    elif cmd == 'onecopy':
        n, starts, seed = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); var = sys.argv[5]
        gen = torch.Generator().manual_seed(seed); pre = None
        if var == 'fixed': comps = [Compressor(n, [0], 'two', gen)]
        elif var == 'rows': comps = [Compressor(n, [0], 'rows', gen)]
        elif var == 'any': comps = [Compressor(n, [s], 'two', gen) for s in range(n)]
        elif var == 'stage2':
            P0 = torch.kron(torch.diag(torch.tensor([1., 1., 0.], dtype=CPX)), torch.eye(3 ** (n - 1), dtype=CPX))
            pre = lambda C: (P0 @ C @ P0) / torch.linalg.norm(P0 @ C @ P0)
            comps = [Compressor(n, [s], 'two', gen) for s in range(1, n)]
        ascent(n, comps, starts, gen, 600, 0.006, 6, pre=pre, label='one-copy ' + var)
    elif cmd == 'verify':
        f, NR = sys.argv[2], int(sys.argv[3]); gen = torch.Generator().manual_seed(1234)
        C0 = np.load(f); n = round(np.log(C0.shape[0]) / np.log(3)); C0 = C0 / np.linalg.norm(C0); C = torch.tensor(C0, dtype=CPX)
        comp = Compressor(n, range(n), 'two', gen); inits = []
        for M in (C0 @ C0.conj().T + C0.conj().T @ C0, C0 @ C0.conj().T, C0.conj().T @ C0):
            Tn = M.reshape((3,) * (2 * n)); ev = []
            for s in range(n):
                Lt = 'abcdefgh'; inp = [Lt[i] for i in range(2 * n)]
                for a in range(n):
                    if a != s: inp[n + a] = inp[a]
                ev.append(np.linalg.eigh(np.einsum(''.join(inp) + '->' + Lt[s] + inp[n + s], Tn))[1])
            for ch in np.ndindex(*(3,) * n):
                inits.append(torch.tensor(np.stack([ev[s][:, ch[s]] for s in range(n)]), dtype=CPX))
        vals = [comp.inner(C, starts=0, init=I0)[0] for I0 in inits] + [comp.inner(C, starts=1)[0] for _ in range(NR)]
        vals = np.sort(np.array(vals))
        print(f"{f.replace(chr(92), '/').split('/')[-1]}: n = {n}, R(C) = {R(C, n).item():.10f}")
        print(f"min of R(P C P) over {len(vals)} inner runs ({len(inits)} eigenvector starts, {NR} random) = {vals[0]:.10f};"
              f" gap = {vals[0] - R(C, n).item():+.10f}; five smallest {np.round(vals[:5], 9).tolist()}")
    elif cmd == 'independent':
        NR = int(sys.argv[2]); gen = torch.Generator().manual_seed(5); I3 = torch.eye(3, dtype=CPX)
        for f in sys.argv[3:]:
            C0 = np.load(f); n = round(np.log(C0.shape[0]) / np.log(3)); C = torch.tensor(C0 / np.linalg.norm(C0), dtype=CPX)
            best = np.inf
            for k in range(NR):
                x = torch.randn(12 * n, generator=gen).requires_grad_(True)

                def build(x):
                    A = [qproj(torch.complex(x[12 * s:12 * s + 3], x[12 * s + 3:12 * s + 6])) for s in range(n)]
                    B = [qproj(torch.complex(x[12 * s + 6:12 * s + 9], x[12 * s + 9:12 * s + 12])) for s in range(n)]
                    return kron_list(A) @ C @ kron_list(B)
                opt = torch.optim.LBFGS([x], max_iter=300, tolerance_grad=1e-12, tolerance_change=1e-15, line_search_fn='strong_wolfe')

                def cl():
                    opt.zero_grad(); X = build(x); nq = torch.vdot(X.reshape(-1), X.reshape(-1)).real
                    v = R(X, n) + 1e3 * torch.relu(1e-4 - nq); v.backward(); return v
                opt.step(cl)
                with torch.no_grad():
                    X = build(x)
                    if torch.vdot(X.reshape(-1), X.reshape(-1)).real.item() > 1e-5: best = min(best, R(X, n).item())
            print(f"{f.replace(chr(92), '/').split('/')[-1]:36s} n={n} R(C) = {R(C, n).item():+.6f}  min over independent row/column"
                  f" planes of R(P_A C P_B) = {best:+.6f}  ({NR} starts)", flush=True)
    elif cmd == 'orbit':
        n, trials, seed = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); gen = torch.Generator().manual_seed(seed); d = 3 ** n
        basis = []
        for i in range(3):
            for j in range(3):
                E = torch.zeros(3, 3, dtype=CPX)
                if i == j and i < 2: E[i, i] = 1; E[2, 2] = -1; basis.append(E)
                elif i < j: E[i, j] = 1; E[j, i] = 1; basis.append(E); F = torch.zeros(3, 3, dtype=CPX); F[i, j] = 1j; F[j, i] = -1j; basis.append(F)
        Bs = torch.stack(basis)
        for tr in range(trials):
            C = rank_two(torch.randn(8 * d, generator=gen), d)
            p = torch.zeros(8 * n, requires_grad=True)
            G = lambda p: kron_list([torch.linalg.matrix_exp(torch.einsum('k,kij->ij', p[8 * s:8 * s + 8].to(CPX), Bs)) for s in range(n)])
            opt = torch.optim.LBFGS([p], max_iter=3000, tolerance_grad=1e-11, tolerance_change=1e-15, line_search_fn='strong_wolfe')

            def cl():
                opt.zero_grad(); g = G(p); v = R(g @ C @ g.conj().T, n); v.backward(); return v
            opt.step(cl)
            with torch.no_grad():
                g = G(p); Cf = g @ C @ g.conj().T; Cf = Cf / torch.linalg.norm(Cf)
                M = (Cf @ Cf.conj().T + Cf.conj().T @ Cf).reshape((3,) * (2 * n)); th = []
                for s in range(n):
                    Lt = 'abcdefgh'; inp = [Lt[i] for i in range(2 * n)]
                    for a in range(n):
                        if a != s: inp[n + a] = inp[a]
                    Ms = torch.einsum(''.join(inp) + '->' + Lt[s] + inp[n + s], M).numpy()
                    if not np.all(np.isfinite(Ms)):
                        th.append(float('nan')); continue
                    w = np.linalg.eigvalsh((Ms + Ms.conj().T) / 2); th.append(w[0] / w.sum())
            print(f"trial {tr}: R(C) = {R(C, n).item():+.4f} -> orbit minimum {R(Cf, n).item():+.2e}, |p| = {p.norm().item():.1f},"
                  f" smallest local eigenvalue per copy {' '.join(f'{x:.1e}' for x in th)}", flush=True)


if __name__ == '__main__':
    main()
