"""Matrix-free universal-anchor test (NOTES.md, Section 5; used for the n = 5 numerics).
Y_p = (M^{(x)n} - R_p)^Gamma vanishes on K = pbar (x) C^{3^n}; p is a universal anchor iff Y_p >= 0, i.e. iff the
smallest eigenvalue of Y_p on K^perp is >= 0 (K^perp is invariant since Y_p is Hermitian).
Vectors of row (x) col are arrays Z[a, y] (a = row string, y = column string), shape (3,)*n + (3,)*n.
One copy: M_1[a, y, b, y'] = d_ab d_yy' - 1/2 d_ay d_by';  A_1 = (M_1)^Gamma: A_1[a, y', b, y] = M_1[b, y', a, y].
G_x = M(p (x) e_x);  D[x, z] = sum_a conj(p_a) G_z[a, x];  H_x = sum_z Dinv[x, z] conj(G_z);
(R^Gamma Z)[a, y'] = sum_x (H_x Z^T G_x)[a, y'].
Usage: python anchor_mf.py n kind [seed]   kind in {rand, randc, W, GHZ, dicke2}."""
import numpy as np, sys, time
from scipy.sparse.linalg import LinearOperator, eigsh

def M1():
    m = np.zeros((3, 3, 3, 3))
    for a in range(3):
        for y in range(3):
            for b in range(3):
                for yp in range(3):
                    m[a, y, b, yp] = float(a == b and y == yp) - 0.5 * float(a == y and b == yp)
    return m

def apply_local(Z, m, n):
    T = Z.reshape((3,) * (2 * n))
    for t in range(n):
        T = np.tensordot(m, T, axes=([2, 3], [t, n + t]))
        T = np.moveaxis(T, [0, 1], [t, n + t])
    return T.reshape(3 ** n, 3 ** n)

def build(p, n):
    N = 3 ** n; m = M1()
    dt = complex if np.iscomplexobj(p) else float
    G = np.empty((N, N, N), dt)
    for x in range(N):
        Z = np.zeros((N, N), dt); Z[:, x] = p
        G[x] = apply_local(Z, m, n)
    D = np.einsum('a,zax->xz', p.conj(), G)
    Dinv = np.linalg.inv(D)
    H = np.einsum('xz,zab->xab', Dinv, G.conj())
    mGam = m.transpose(2, 1, 0, 3)
    return dict(n=n, N=N, G=G, H=H, mGam=mGam, phat=p.conj() / np.linalg.norm(p), dt=dt, D=D)

def Y_matvec(Zv, B):
    N = B['N']; Z = Zv.reshape(N, N)
    AZ = apply_local(Z, B['mGam'], B['n'])
    W = np.einsum('xay,by->xab', B['H'], Z, optimize=True)
    RZ = np.einsum('xab,xbc->ac', W, B['G'], optimize=True)
    return (AZ - RZ).reshape(-1)

def min_eig_complement(B, shift=50.0, tol=1e-5, ncv=40):
    N = B['N']; ph = B['phat']
    def proj(Zv):
        Z = Zv.reshape(N, N); return (Z - np.outer(ph, ph.conj() @ Z)).reshape(-1)
    def mv(v):
        v = np.asarray(v).reshape(-1)
        pv = proj(v)
        return proj(Y_matvec(pv, B)) + shift * (v - pv)
    op = LinearOperator((N * N, N * N), matvec=mv, dtype=B['dt'])
    vals = eigsh(op, k=1, which='SA', tol=tol, ncv=ncv, maxiter=20000, return_eigenvectors=False)
    return vals[0]

def vector(kind, n, rng):
    N = 3 ** n
    if kind == 'rand': p = rng.normal(size=N)
    elif kind == 'randc': p = rng.normal(size=N) + 1j * rng.normal(size=N)
    else:
        p = np.zeros(N)
        for idx in range(N):
            s = np.base_repr(idx, 3).zfill(n)
            ones = s.count('1'); zeros = s.count('0'); twos = s.count('2')
            if kind == 'W' and ones == 1 and zeros == n - 1: p[idx] = 1
            if kind == 'GHZ' and (zeros == n or ones == n or twos == n): p[idx] = 1
            if kind == 'dicke2' and ones == 2 and zeros == n - 2: p[idx] = 1
    return p / np.linalg.norm(p)

if __name__ == '__main__':
    n = int(sys.argv[1]); kind = sys.argv[2]; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    rng = np.random.default_rng(seed)
    t0 = time.time(); p = vector(kind, n, rng); B = build(p, n)
    lm = min_eig_complement(B)
    print(f"n={n} {kind} seed={seed}: min eig of Y_p on complement of kernel = {lm:+.5f}  [{time.time()-t0:.0f}s]", flush=True)
