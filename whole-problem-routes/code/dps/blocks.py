"""Float block machinery for n copies of the per-copy algebra B_3(3) (strands A, B, B').

Block coordinates of one copy: concatenation over spins j = 0..3 of vec_C(m_j x m_j block), 15 numbers.
For n copies the coordinate tensor is 15^n (C order), and block beta = (j_1..j_n) is read off by idx[beta].
"""
import itertools
import numpy as np
import scipy.sparse as sps
import brauer as br

SPINS = [0, 1, 2, 3]


def per_copy_float():
    Bex = br.spin_m0_bases()
    U = {}
    for j in SPINS:
        Bj = np.array(Bex[j].tolist(), dtype=float)
        q, r = np.linalg.qr(Bj)
        U[j] = q
    m = {j: U[j].shape[1] for j in SPINS}
    off = {}
    o = 0
    for j in SPINS:
        off[j] = o
        o += m[j] ** 2
    assert o == br.NDIAG
    Dm = [M.astype(float) for M in br.DMAT]
    # sanity: each diagram is block diagonal on the m=0 space, cross-spin blocks vanish
    Uall = np.hstack([U[j] for j in SPINS])
    for M in Dm:
        X = Uall.T @ M @ Uall
        o1 = 0
        for j in SPINS:
            o2 = 0
            for k in SPINS:
                if j != k:
                    assert np.abs(X[o1:o1 + m[j], o2:o2 + m[k]]).max() < 1e-12
                o2 += m[k]
            o1 += m[j]
    Phi = np.zeros((br.NDIAG, br.NDIAG))
    for a, M in enumerate(Dm):
        col = []
        for j in SPINS:
            col.append((U[j].T @ M @ U[j]).reshape(-1))
        Phi[:, a] = np.concatenate(col)
    assert np.linalg.matrix_rank(Phi) == br.NDIAG
    Phinv = np.linalg.inv(Phi)

    def ptmat(strands):
        pi = br.pt_perm(strands)
        T = np.zeros((br.NDIAG, br.NDIAG))
        for a in range(br.NDIAG):
            T[pi[a], a] = 1.0
        return Phi @ T @ Phinv

    nm = br.named()
    S = Dm[nm['S_BB']]
    sigma = {j: U[j].T @ S @ U[j] for j in SPINS}
    return dict(U=U, m=m, off=off, Phi=Phi, Phinv=Phinv, ptmat=ptmat, sigma=sigma, Dm=Dm, named=nm)


def global_blocks(pc, n):
    """list of blocks: (beta, b, idx (b x b flat global indices), Pi (orthonormal basis of Sym range), bplus)"""
    m, off = pc['m'], pc['off']
    out = []
    for beta in itertools.product(SPINS, repeat=n):
        ms = [m[j] for j in beta]
        b = int(np.prod(ms))
        rows = list(itertools.product(*[range(x) for x in ms]))
        idx = np.zeros((b, b), dtype=np.int64)
        for R, r in enumerate(rows):
            for C_, c in enumerate(rows):
                flat = 0
                for s in range(n):
                    k = off[beta[s]] + r[s] * ms[s] + c[s]
                    flat = flat * br.NDIAG + k
                idx[R, C_] = flat
        sig = np.array([[1.0]])
        for j in beta:
            sig = np.kron(sig, pc['sigma'][j])
        w, v = np.linalg.eigh(sig)
        assert np.allclose(np.abs(w), 1)
        Pi = v[:, w > 0]
        out.append(dict(beta=beta, b=b, idx=idx, Pi=Pi, bp=Pi.shape[1]))
    return out


def kron_chain_mats(M, n):
    """sparse factors K_s = I (x) .. M (at position s) .. (x) I, s = 0..n-1"""
    I = sps.identity(br.NDIAG, format='csr')
    out = []
    for s in range(n):
        f = None
        for t in range(n):
            x = sps.csr_matrix(M) if t == s else I
            f = x if f is None else sps.kron(f, x, format='csr')
        out.append(f)
    return out


def coeff_power(v, n):
    out = np.array([1.0])
    for _ in range(n):
        out = np.kron(out, v)
    return out
