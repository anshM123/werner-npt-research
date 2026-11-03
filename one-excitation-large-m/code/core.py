"""Independent check of the one-excitation recursion and bridge (written from the definitions in two-sector-class/README.md, not from the code of the working notes).

Conventions
-----------
* Operators are torch complex128 matrices.
* qform(Z, dims) = sum_J (-1/2)^{|J|} ||Tr_J Z||_F^2  (the endpoint form q_n on factors `dims`).
* Free copy is always the first tensor factor (dimension 3).
* One-free-copy one-excitation class at level m: V = C^3 (x) (A (+) C e), A = (C^2)^{(x) m}.
  A vector of V is stored with index i*(N+1)+a, i in {0,1,2} (free), a in {0..N-1} for A, a = N for e,
  where N = 2^m.  C = R L^*, R = [R_A; R_B], L = [L_A; L_B], two columns.
  Blocks: X = R_A L_A^* (on C^3 (x) A), P = R_A L_B^*, Q = R_B L_A^*, Y = R_B L_B^* (3x3).
"""
import itertools, math
import numpy as np
import torch

torch.set_default_dtype(torch.float64)
CD = torch.complex128


def fro2(T):
    return (T.real ** 2 + T.imag ** 2).sum()


def qform(Z, dims, w=-0.5):
    """sum_J w^{|J|} ||Tr_J Z||^2 over all subsets J of the factors."""
    n = len(dims)
    T = Z.reshape(*dims, *dims)
    out = [0.0]

    def rec(T, k, kept, coef):
        if k == n:
            out[0] = out[0] + coef * fro2(T)
            return
        nrem = kept + (n - k)
        rec(T, k + 1, kept + 1, coef)
        Tt = torch.diagonal(T, dim1=kept, dim2=kept + nrem).sum(-1)
        rec(Tt, k + 1, kept, coef * w)

    rec(T, 0, 0, 1.0)
    return out[0]


def qbil(Z1, Z2, dims, w=-0.5):
    """Sesquilinear polarisation (conjugate-linear in Z1) of qform."""
    n = len(dims)
    T1 = Z1.reshape(*dims, *dims)
    T2 = Z2.reshape(*dims, *dims)
    out = [0.0]

    def rec(A, B, k, kept, coef):
        if k == n:
            out[0] = out[0] + coef * (A.conj() * B).sum()
            return
        nrem = kept + (n - k)
        rec(A, B, k + 1, kept + 1, coef)
        At = torch.diagonal(A, dim1=kept, dim2=kept + nrem).sum(-1)
        Bt = torch.diagonal(B, dim1=kept, dim2=kept + nrem).sum(-1)
        rec(At, Bt, k + 1, kept, coef * w)

    rec(T1, T2, 0, 0, 1.0)
    return out[0]


def q1(Z):
    return fro2(Z) - 0.5 * (torch.trace(Z).abs() ** 2)


# ---------------------------------------------------------------- README-level quantities at level m

def blocks_from_factors(R, L, N):
    """R, L: (3*(N+1), 2).  Returns R_A (3N,2), R_B (3,2), L_A, L_B."""
    k = R.shape[1]
    Rr = R.reshape(3, N + 1, k)
    Lr = L.reshape(3, N + 1, k)
    RA = Rr[:, :N, :].reshape(3 * N, k)
    RB = Rr[:, N, :]
    LA = Lr[:, :N, :].reshape(3 * N, k)
    LB = Lr[:, N, :]
    return RA, RB, LA, LB


def blocks(R, L, N):
    RA, RB, LA, LB = blocks_from_factors(R, L, N)
    X = RA @ LA.conj().T
    P = RA @ LB.conj().T
    Q = RB @ LA.conj().T
    Y = RB @ LB.conj().T
    return X, P, Q, Y


def qtilde_P(P, N):
    """P: (3N, 3) maps C^3 (x) e -> C^3 (x) A.  q~(P) = ||P||^2 - 1/2 ||Tr_0 P||^2."""
    Pr = P.reshape(3, N, 3)
    tr0 = torch.einsum('iai->a', Pr)
    return fro2(P) - 0.5 * fro2(tr0)


def qtilde_Q(Q, N):
    """Q: (3, 3N) maps C^3 (x) A -> C^3 (x) e."""
    Qr = Q.reshape(3, 3, N)
    tr0 = torch.einsum('iib->b', Qr)
    return fro2(Q) - 0.5 * fro2(tr0)


def Sigma_full(X, m):
    """2^m q_{m+1}(X) on factors (3, 2, ..., 2)."""
    return (2 ** m) * qform(X, (3,) + (2,) * m)


def Gamma_X(X, m):
    """2 q_{0,K}(X): free copy and all of A merged into one copy."""
    return 2 * qform(X, (3, 2 ** m))


PAULI = [torch.eye(2, dtype=CD),
         torch.tensor([[0, 1], [1, 0]], dtype=CD),
         torch.tensor([[0, -1j], [1j, 0]], dtype=CD),
         torch.tensor([[1, 0], [0, -1]], dtype=CD)]


def Sigma_full_literal(X, m):
    """README literal definition: sum over mu in {1,2,3}^m of q_1(X^_mu), X^_mu = Tr_A[(1 (x) sigma_mu) X]."""
    N = 2 ** m
    Xr = X.reshape(3, N, 3, N)
    tot = 0.0
    for mu in itertools.product((1, 2, 3), repeat=m):
        S = PAULI[mu[0]]
        for k in mu[1:]:
            S = torch.kron(S, PAULI[k])
        # Tr_A[(1 (x) S) X]_{ij} = sum_{a,b} S_{ba} X_{(i,a),(j,b)}
        Xh = torch.einsum('ba,iajb->ij', S, Xr)
        tot = tot + q1(Xh)
    return tot


def Gamma_X_literal(X, m):
    """README literal: 2 sum_{alpha,beta} q_1(X_{alpha beta}) - q_1(X^_0)."""
    N = 2 ** m
    Xr = X.reshape(3, N, 3, N)
    tot = 0.0
    for a in range(N):
        for b in range(N):
            tot = tot + 2 * q1(Xr[:, a, :, b])
    X0 = torch.einsum('iaja->ij', Xr)
    return tot - q1(X0)


def star_parts(R, L, m):
    N = 2 ** m
    X, P, Q, Y = blocks(R, L, N)
    return dict(Sig=Sigma_full(X, m), Gam=Gamma_X(X, m), y=q1(Y), p=qtilde_P(P, N), q=qtilde_Q(Q, N))


def star_ratio_from_parts(d):
    """sqrt(q1 Y)(sqrt Gamma - sqrt Sigma)/sqrt(p q), computed stably."""
    Sig = torch.clamp(d['Sig'], min=0.0)
    Gam = torch.clamp(d['Gam'], min=0.0)
    num = torch.sqrt(torch.clamp(d['y'], min=0.0)) * (Gam - Sig) / (torch.sqrt(Gam) + torch.sqrt(Sig) + 1e-300)
    den = torch.sqrt(torch.clamp(d['p'] * d['q'], min=1e-300))
    return num / den


# ---------------------------------------------------------------- embedding into (C^3)^{(x)(m+1)} for convention checks

def embed_isometry(m):
    """Isometry J : C^3 (x) (A (+) C e) -> (C^3)^{(x)(m+1)}, K_s = span{0,1}, e_s = |2>."""
    N = 2 ** m
    dim_big = 3 ** (m + 1)
    J = torch.zeros(dim_big, 3 * (N + 1), dtype=CD)
    for i in range(3):
        for a in range(N + 1):
            if a < N:
                bits = [(a >> (m - 1 - s)) & 1 for s in range(m)]
            else:
                bits = [2] * m
            idx = i
            for b in bits:
                idx = idx * 3 + b
            J[idx, i * (N + 1) + a] = 1.0
    return J


def qn_big(C, m):
    J = embed_isometry(m)
    Cb = J @ C @ J.conj().T
    return qform(Cb, (3,) * (m + 1))


# ---------------------------------------------------------------- merged forms on C^3 (x) C^2_new (x) C^D

def Gam_m(X, D):
    """Gamma = 2 q_2 with new (x) H merged."""
    return 2 * qform(X, (3, 2 * D))


def Gbar_m(X, D):
    """Gamma_bar = 4 q_3 on (free, new, H)."""
    return 4 * qform(X, (3, 2, D))


def H_m(X, D):
    return Gbar_m(X, D) - Gam_m(X, D)


def Hbil_m(Z1, Z2, D):
    return 4 * qbil(Z1, Z2, (3, 2, D)) - 2 * qbil(Z1, Z2, (3, 2 * D))


def internal_B(R, L, D, R2=None, L2=None):
    """B_{(i,l),(k,j)} = <r_i (x) l_l | K | r2_k (x) l2_j> = Hbil(r_i l2_j^*, r2_k l_l^*).
    With R2=R, L2=L this is the internal two-qubit Gram (R^* (x) L^*) K (R (x) L)."""
    if R2 is None:
        R2, L2 = R, L
    B = torch.zeros(4, 4, dtype=CD)
    for i in range(2):
        for l in range(2):
            for k in range(2):
                for j in range(2):
                    Z1 = torch.outer(R[:, i], L2[:, j].conj())
                    Z2 = torch.outer(R2[:, k], L[:, l].conj())
                    B[2 * i + l, 2 * k + j] = Hbil_m(Z1, Z2, D)
    return B


def gram_sum(pairs, D):
    """Gram of the purification sum_p sqrt(K)(R_p (x) L_p):  sum_{p,p'} (R_p (x) L_p)^* K (R_p' (x) L_p')."""
    B = torch.zeros(4, 4, dtype=CD)
    for (R1, L1) in pairs:
        for (R2, L2) in pairs:
            B = B + internal_B(R1, L1, D, R2, L2)
    return B


SY = np.array([[0, -1j], [1j, 0]])
SYY = np.kron(SY, SY)


def concurrence(B):
    """Homogeneous Wootters concurrence of an (unnormalised) PSD 4x4 matrix (pure cost 2|det x|)."""
    B = np.asarray(B.detach().numpy() if torch.is_tensor(B) else B, dtype=complex)
    B = 0.5 * (B + B.conj().T)
    w, V = np.linalg.eigh(B)
    w = np.clip(w, 0, None)
    sB = (V * np.sqrt(w)) @ V.conj().T
    M = sB @ SYY @ sB.conj()
    lam = np.sort(np.linalg.svd(M, compute_uv=False))[::-1]
    return max(0.0, lam[0] - lam[1] - lam[2] - lam[3])


def witness(A):
    """W_A = (|vec A><vec A|)^Gamma, W[(k,j),(i,l)] = A_kl conj(A_ij)."""
    A = np.asarray(A, dtype=complex)
    W = np.zeros((4, 4), dtype=complex)
    for k in range(2):
        for j in range(2):
            for i in range(2):
                for l in range(2):
                    W[2 * k + j, 2 * i + l] = A[k, l] * np.conj(A[i, j])
    return W


# ---------------------------------------------------------------- product frame for 2I - F on the new qubit

def octahedral_frame():
    """Six atoms (u, v), weight 1: antipodal pairs of the x, y, z Bloch axes, both orders."""
    atoms = []
    for P in PAULI[1:]:
        w, V = torch.linalg.eigh(P)
        up, dn = V[:, 1], V[:, 0]
        atoms.append((up, dn))
        atoms.append((dn, up))
    return [(1.0, u, v) for (u, v) in atoms]


def contract_new(RA, u, D):
    """R_{A,u} = (1 (x) <u| (x) 1) R_A : rows (free, new, H) -> (free, H)."""
    Rr = RA.reshape(3, 2, D, -1)
    return torch.einsum('a,iahk->ihk', u.conj(), Rr).reshape(3 * D, -1)


def rand_c(*shape, gen=None):
    return torch.randn(*shape, generator=gen) + 1j * torch.randn(*shape, generator=gen)


def rand_isometry(n, k, gen=None):
    Qm, _ = torch.linalg.qr(rand_c(n, k, gen=gen).to(CD))
    return Qm


# ---------------------------------------------------------------- explicit replica kernel K (fast path)
_KCACHE = {}


def swap_perm(dims, s):
    """Permutation matrix swapping factor s between the two replicas of (x)dims."""
    n = len(dims)
    full = list(dims) + list(dims)
    N = int(np.prod(full))
    idx = np.arange(N).reshape(full)
    perm_axes = list(range(2 * n))
    perm_axes[s], perm_axes[n + s] = perm_axes[n + s], perm_axes[s]
    src = idx.transpose(perm_axes).reshape(-1)
    P = torch.zeros(N, N, dtype=CD)
    P[torch.arange(N), torch.tensor(src)] = 1.0
    return P


def K_matrix(D):
    """K = 2 (I - F_free/2)(I - F_new)(I - F_H) on (C^3 (x) C^2 (x) C^D)^{(x)2}."""
    if D not in _KCACHE:
        dims = (3, 2, D)
        I = torch.eye((6 * D) ** 2, dtype=CD)
        F0, F1, F2 = swap_perm(dims, 0), swap_perm(dims, 1), swap_perm(dims, 2)
        _KCACHE[D] = 2 * (I - 0.5 * F0) @ (I - F1) @ (I - F2)
    return _KCACHE[D]


def kron2(A, B):
    return torch.einsum('ij,kl->ikjl', A, B).reshape(A.shape[0] * B.shape[0], A.shape[1] * B.shape[1])


def internal_B_fast(R, L, D, R2=None, L2=None):
    if R2 is None:
        R2, L2 = R, L
    K = K_matrix(D)
    return kron2(R, L).conj().T @ K @ kron2(R2, L2)


def gram_sum_fast(pairs, D):
    T = sum(kron2(R, L) for (R, L) in pairs)
    return T.conj().T @ K_matrix(D) @ T
