"""Per-site algebra for the four-replica relaxation (local dimension d = 3, four replicas -> 81-dim per site).

Conventions
  replicas are 0,1,2,3 (the text calls them 1,2,3,4);  perm_op(s) sends the content of replica i to replica s[i];
  perm_op(a) @ perm_op(b) = perm_op(a o b) with (a o b)(i) = a(b(i)).
  An algebra is given by a real spanning set of 81 x 81 matrices closed under transpose.  Alg computes an orthogonal Q
  with Q^T M Q = (+)_nu rho_nu(M) (x) I_{m_nu}  (columns ordered (irrep nu, irrep-basis index e, copy c)).
  Fourier coordinates of an algebra element X:  X_nu[k, l] = Tr(B_k^T X B_l) / m_nu  (B_k = columns of copy-set k).
  "hat" coordinates (used for the variables):  Xhat_nu = m_nu * X_nu, so that Tr[X Y] = sum_nu Tr[Xhat_nu Y_nu].
"""
import itertools, numpy as np

d = 3
NR = 4
DIM = d ** NR

def mul(a, b): return tuple(a[b[i]] for i in range(len(a)))
def inv(a):
    r = [0] * len(a)
    for i, ai in enumerate(a): r[ai] = i
    return tuple(r)
def cycles(s):
    seen = set(); c = 0
    for i in range(len(s)):
        if i not in seen:
            c += 1; j = i
            while j not in seen: seen.add(j); j = s[j]
    return c

S4 = list(itertools.permutations(range(4)))
S3 = [s for s in S4 if s[3] == 3]                       # fixes replica 3
S2S2 = [s for s in S4 if set(s[:2]) == {0, 1}]          # preserves {0,1} and {2,3}

def perm_op(s, nr=NR):
    D = d ** nr
    P = np.zeros((D, D))
    for x in itertools.product(range(d), repeat=nr):
        y = [0] * nr
        for i in range(nr): y[s[i]] = x[i]
        P[np.ravel_multi_index(y, (d,) * nr), np.ravel_multi_index(x, (d,) * nr)] = 1.0
    return P

def pt(M, reps, nr=NR):
    """partial transpose on the replicas in reps"""
    T = M.reshape((d,) * (2 * nr)); perm = list(range(2 * nr))
    for r in reps: perm[r], perm[nr + r] = perm[nr + r], perm[r]
    return T.transpose(perm).reshape(d ** nr, d ** nr)

def ptrace_last(M, nr=NR):
    D1 = d ** (nr - 1)
    return np.einsum('aibi->ab', M.reshape(D1, d, D1, d))

def embed_last(M):
    return np.kron(M, np.eye(d))

PERM = {s: perm_op(s) for s in S4}


class Alg:
    """Fourier decomposition of a real semisimple *-algebra of 'real type' given by a spanning set."""
    def __init__(self, Ms, name, seed=2024, tol=1e-7, order_key=None):
        self.name = name
        D = Ms[0].shape[0]; self.D = D
        rng = np.random.default_rng(seed)
        A = sum(rng.standard_normal() * M for M in Ms); H = (A + A.T) / 2
        ev, V = np.linalg.eigh(H)
        sc = max(1.0, np.abs(ev).max())
        spaces = []; st = 0
        for k in range(1, D + 1):
            if k == D or ev[k] - ev[k - 1] > tol * sc:
                spaces.append(V[:, st:k]); st = k
        Rr = sum(rng.standard_normal() * M for M in Ms)
        ns = len(spaces); par = list(range(ns))
        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]; x = par[x]
            return x
        for i in range(ns):
            for j in range(i + 1, ns):
                if np.linalg.norm(spaces[i].T @ Rr @ spaces[j]) > 1e-6 or np.linalg.norm(spaces[j].T @ Rr @ spaces[i]) > 1e-6:
                    par[find(i)] = find(j)
        groups = {}
        for i in range(ns): groups.setdefault(find(i), []).append(i)
        blocks = []
        for g in groups.values():
            m = spaces[g[0]].shape[1]; f = len(g)
            assert all(spaces[i].shape[1] == m for i in g), "eigenspace dims differ inside a block"
            B1 = spaces[g[0]]; Bs = [B1]
            for i in g[1:]:
                T = spaces[i].T @ Rr @ B1
                s2 = np.trace(T.T @ T) / m
                assert np.allclose(T.T @ T, s2 * np.eye(m), atol=1e-8 * max(1.0, s2)), "alignment failed"
                Bs.append(spaces[i] @ T / np.sqrt(s2))
            blocks.append(dict(f=f, m=m, B=Bs))
        if order_key is not None:
            blocks.sort(key=lambda b: order_key(self, b))
        else:
            blocks.sort(key=lambda b: (-b['m'], b['f']))
        self.blocks = blocks
        self.nb = len(blocks)
        self.f = [b['f'] for b in blocks]
        self.m = [b['m'] for b in blocks]
        self.Q = [np.hstack(b['B']) for b in blocks]          # D x (f m), columns (e, c)
        self.dim = sum(f * f for f in self.f)
        # verification: every spanning element is block diagonal of the form rho (x) I
        err = 0.0
        Qall = np.hstack(self.Q)
        assert Qall.shape == (D, D) and np.allclose(Qall.T @ Qall, np.eye(D), atol=1e-9)
        for M in Ms:
            Y = Qall.T @ M @ Qall
            rec = np.zeros_like(Y); off = 0
            for b, rho in zip(blocks, self.coeffs(M)):
                f, m = b['f'], b['m']
                rec[off:off + f * m, off:off + f * m] = np.kron(rho, np.eye(m)); off += f * m
            err = max(err, np.abs(Y - rec).max())
        self.verify_err = err
        assert err < 1e-9, f"{name}: Fourier decomposition check failed ({err})"

    def coeffs(self, X):
        """Fourier coordinates X_nu (f x f) of an algebra element X"""
        out = []
        for b, Q in zip(self.blocks, self.Q):
            f, m = b['f'], b['m']
            Y = (Q.T @ X @ Q).reshape(f, m, f, m)
            out.append(np.einsum('kclc->kl', Y) / m)
        return out

    def unit(self, nu, i, j):
        b = self.blocks[nu]
        return b['B'][i] @ b['B'][j].T

    def rep(self, M):
        """rho_nu(M) for an algebra element M (same as coeffs)"""
        return self.coeffs(M)

    def from_coeffs(self, C):
        X = np.zeros((self.D, self.D))
        for nu, Cn in enumerate(C):
            f = self.f[nu]
            for i in range(f):
                for j in range(f):
                    if Cn[i, j] != 0: X += Cn[i, j] * self.unit(nu, i, j)
        return X


def site_map(phi, src, tgt):
    """hat-coordinate tensors of a linear map phi (dense 81x81 -> dense) from algebra src into algebra tgt:
    T[(mu, lam)][k, l, i, j] = Tr(B'_k^T phi(E^lam_ij) B'_l) / m_lam,  so that  (phi Z)hat_mu = sum_lam T . Zhat_lam."""
    T = {}
    resid = 0.0
    for lam in range(src.nb):
        f = src.f[lam]
        for i in range(f):
            for j in range(f):
                Y = phi(src.unit(lam, i, j))
                C = tgt.coeffs(Y)
                # check that phi(E) lies in the target algebra
                Yr = tgt.from_coeffs(C)
                resid = max(resid, np.abs(Y - Yr).max())
                for mu in range(tgt.nb):
                    if (mu, lam) not in T:
                        T[(mu, lam)] = np.zeros((tgt.f[mu], tgt.f[mu], f, f))
                    T[(mu, lam)][:, :, i, j] = C[mu] * tgt.m[mu] / src.m[lam]
    assert resid < 1e-9, f"map image not in target algebra ({resid})"
    return T


def site_functional(X, alg):
    """per-site functional: value on Zhat block lam is sum_ij Zhat[i,j] * c[lam][i,j]  with c = X_lam^T"""
    return [C.T.copy() for C in alg.coeffs(X)]


def build_algebras():
    A = Alg([PERM[s] for s in S4], 'A(S4)')
    B1 = Alg([pt(PERM[s], [3]) for s in S4], 'B31')
    B2 = Alg([pt(PERM[s], [2, 3]) for s in S4], 'B22',
             order_key=lambda alg, b: (-b['m'], b['f'], float(np.trace(b['B'][0].T @ PERM[(1, 0, 2, 3)] @ b['B'][0]))))
    S3_3 = list(itertools.permutations(range(3)))
    A3 = Alg([perm_op(s, 3) for s in S3_3], 'A3(S3)')
    return A, B1, B2, A3


if __name__ == "__main__":
    A, B1, B2, A3 = build_algebras()
    for alg in (A, B1, B2, A3):
        print(f"{alg.name}: irreps (f, m) = {list(zip(alg.f, alg.m))}, dim = {alg.dim}, check err = {alg.verify_err:.1e}")
