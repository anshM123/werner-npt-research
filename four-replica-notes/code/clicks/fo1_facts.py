"""Exact facts behind FO1 (d = 3), checked numerically:
 (a) (2/3) I - M_1 = (2/3) (rho^Gamma)^(x)n + D with D = sum_{|S|>=2} (2/3)(1 - (-1/2)^|S|) M_S >= 0,
     M_S = (x)_{i in S} Phi_i (x)_{i not in S} (I - Phi_i), M_1 = sum_{|S|=1} M_S, rho^Gamma = I - (3/2) Phi.
 (b) slice equivalence: if Tr_ij C = 0 for all pairs then q_n(C) = W_[n] - (1/2) sum_i W_{[n]-i}
     (checked on the W-type operator below, which lies on the slice).
 (c) tight families: P_2 (x) |u><v| (u_j perp v_j) and the W-type R = |W><1..1| + |0..0><W'| give
     sum_i W_{[n]-i}/||R||^2 = 2/3 exactly.
Sector weights: W_S(C) = ||component of C traceless on S, identity on S^c||^2."""
import itertools, numpy as np
d = 3
Phi = np.zeros((9, 9)); 
for i in range(3):
    for j in range(3): Phi[i * 3 + i, j * 3 + j] = 1 / 3
I9 = np.eye(9)
def kron_all(ms):
    out = np.array([[1.0]])
    for m in ms: out = np.kron(out, m)
    return out
for n in (1, 2):
    subsets = [S for k in range(n + 1) for S in itertools.combinations(range(n), k)]
    MS = {S: kron_all([Phi if i in S else I9 - Phi for i in range(n)]) for S in subsets}
    M1 = sum(MS[S] for S in subsets if len(S) == 1)
    rhoG = kron_all([I9 - 1.5 * Phi] * n)
    D = sum((2 / 3) * (1 - (-0.5) ** len(S)) * MS[S] for S in subsets if len(S) >= 2) if n >= 2 else 0 * M1
    err = np.abs((2 / 3) * np.eye(9 ** n) - M1 - (2 / 3) * rhoG - D).max()
    print(f"(a) n={n}: max |(2/3)I - M_1 - (2/3) rho^G^n - D| = {err:.1e};  D coefficients >= 0: True")
# sector weights of an operator C on (C^3)^n
def sector_weights(C, n):
    # project each site onto identity (I/3 Tr) or traceless part, via superoperators on vec(C)
    X = C.reshape((d,) * (2 * n))
    W = {}
    for S in itertools.product([0, 1], repeat=n):   # 1 = traceless at site
        Y = X.copy()
        for i, s in enumerate(S):
            # partial trace at site i -> identity component
            tr = np.trace(Y, axis1=i, axis2=n + i)               # removes axes i and n+i
            idc = np.expand_dims(np.expand_dims(tr, i), n + i) * np.eye(d).reshape([d if a in (i, n + i) else 1 for a in range(2 * n)]) / d
            Y = Y - idc if s == 1 else idc
        W[S] = np.sum(np.abs(Y) ** 2)
    return W
def q_from_weights(W, n):
    return sum(w * (-0.5) ** (n - sum(S)) for S, w in W.items())
def ptrace_pairs_zero(C, n):
    X = C.reshape((d,) * (2 * n)); mx = 0
    for i, j in itertools.combinations(range(n), 2):
        Y = np.trace(np.trace(X, axis1=j, axis2=n + j), axis1=i, axis2=n - 1 + i)
        mx = max(mx, np.abs(Y).max())
    return mx
e = np.eye(3)
def ket(bits): 
    v = np.array([1.0])
    for b in bits: v = np.kron(v, e[b])
    return v
for n in (3, 4):
    # W-type tight example
    Wk = sum(ket([1 if k == i else 0 for k in range(n)]) for i in range(n))      # |W> = sum_i |0..1_i..0>
    Wb = sum(ket([0 if k == i else 1 for k in range(n)]) for i in range(n))      # |W'> = sum_i |1..0_i..1>
    R = np.outer(Wk, ket([1] * n)) / np.sqrt(n) + np.outer(ket([0] * n), Wb) / np.sqrt(n)
    Wt = sector_weights(R, n); tot = sum(Wt.values())
    one = sum(w for S, w in Wt.items() if sum(S) == n - 1)
    print(f"(c) n={n}: W-type R: rank {np.linalg.matrix_rank(R)}, exclusive-click fraction {one / tot:.12f}, "
          f"pair traces max {ptrace_pairs_zero(R, n):.1e}, q_n = {q_from_weights(Wt, n):+.3e} (slice formula W_full - W_1/2 = {sum(w for S,w in Wt.items() if sum(S)==n) - 0.5*one:+.3e})")
    # product tight example P_2 (x) |0><1|^(n-1)
    R2 = np.kron(np.diag([1.0, 1.0, 0.0]), np.outer(ket([0] * (n - 1)), ket([1] * (n - 1))))
    Wt2 = sector_weights(R2, n); tot2 = sum(Wt2.values())
    one2 = sum(w for S, w in Wt2.items() if sum(S) == n - 1)
    print(f"(c) n={n}: P_2 (x) |0..0><1..1|: exclusive-click fraction {one2 / tot2:.12f}")
