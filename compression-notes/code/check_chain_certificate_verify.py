"""Section 9: rigorous check of the certificate written by check_chain_certificate.py.
Claim: X = (A (x) B) rho^{(x)M} (A (x) B)^* lies in K_2(6) = cone{separable, (a (x) b) rho^{(x)2} (a (x) b)^*}.
All parameters (A, B, the atom parameters) are float64 numbers, i.e. exact dyadic rationals, and the weights w_i used
below are float64 numbers >= 0.  The identity  X = sum_i w_i G_i + Z  holds by definition of Z.  Every matrix is
recomputed in float64 with an entrywise error bound |fl(s) - s| <= c K u (sum of |products|), u = 2^-53, K the number of
products per entry, c a safety factor.  Then Z (exact, Hermitian) satisfies
    Tr Z >= t_lo > 0  and  ||Z||_2 <= z_hi,
and if z_hi^2 <= t_lo^2 / (36 - 1), then Z / Tr Z lies in the Gurvits-Barnum ball (Tr sigma^2 <= 1/(N - 1), N = 36), so
Z is separable [Gurvits and Barnum, Phys. Rev. A 66, 062311 (2002)].  Hence X is in K_2(6).
usage: python check_chain_certificate_verify.py M   (NumPy only; about a second)"""
import sys, itertools, numpy as np
M = int(sys.argv[1])
u = 2.0 ** -53
SAFE = 4.0
D = 6; N = D * D
def filt_with_bound(A, B, m):
    """float value and entrywise error bound of (A (x) B) rho^{(x)m} (A (x) B)^*, via sum_S (-1/2)^{|S|} (A(x)B) F_S (A(x)B)^*."""
    subsets = [S for r in range(m + 1) for S in itertools.combinations(range(m), r)]
    At, Bt = A.reshape([D] + [3] * m), B.reshape([D] + [3] * m)
    aA, aB = np.abs(At), np.abs(Bt)
    Y = np.zeros((N, N), complex); Yabs = np.zeros((N, N))
    for S in subsets:
        Sc = [s for s in range(m) if s not in S]
        pa = [0] + [1 + s for s in S] + [1 + s for s in Sc]
        def parts(T):
            return T.transpose(pa).reshape(D, 3 ** len(S), 3 ** len(Sc))
        Ap, Bp, aAp, aBp = parts(At), parts(Bt), parts(aA), parts(aB)
        M1 = np.einsum('axc,eyc->aexy', Ap, Ap.conj()); M2 = np.einsum('byc,fxc->bfyx', Bp, Bp.conj())
        YS = np.einsum('aexy,bfyx->abef', M1, M2).reshape(N, N)
        N1 = np.einsum('axc,eyc->aexy', aAp, aAp); N2 = np.einsum('byc,fxc->bfyx', aBp, aBp)
        YSabs = np.einsum('aexy,bfyx->abef', N1, N2).reshape(N, N)
        Y += (-0.5) ** len(S) * YS; Yabs += 0.5 ** len(S) * YSabs
    K = 9 ** m
    return Y, SAFE * (K + 4 * len(subsets) + 20) * 2 * u * Yabs * (1 + 1e-6)
def prod_with_bound(x, y):
    v = np.kron(x, y); av = np.abs(v)
    return np.outer(v, v.conj()), SAFE * 40 * u * np.outer(av, av)
T = np.load(f"chain_cert_target_M{M}.npz"); A, B = T["A"], T["B"]
X, EX = filt_with_bound(A, B, M)
trX = np.trace(X).real
print(f"M = {M}: Tr X = {trX:.6e}; max entrywise error bound on X {EX.max():.2e}", flush=True)
At = np.load(f"chain_cert_atoms_M{M}.npz", allow_pickle=True)
kinds, coef = At["kinds"], At["coef"]
n_at = len(kinds)
assert np.all(coef >= 0)
Zc = X.copy(); EZ = EX.copy(); absum = np.abs(X).copy()
for k in range(n_at):
    if kinds[k] == 'f':
        G, E = filt_with_bound(At["P"][k], At["Q"][k], 2)
    else:
        G, E = prod_with_bound(At["x"][k], At["y"][k])
    w = coef[k] * trX / np.trace(G).real          # any float w >= 0 is allowed; this one makes Z small
    assert w >= 0
    Zc = Zc - w * G
    EZ += w * E * (1 + 1e-9)
    absum += w * np.abs(G)
# rounding of the n_at products w*G and of the running sum (recursive summation), on top of the atom errors
EZ += SAFE * (n_at + 4) * u * absum * (1 + 1e-9)
Zh = (Zc + Zc.conj().T) / 2                        # exact Z is Hermitian; |Z - Zh| <= max(EZ, EZ^T) entrywise
EZh = np.maximum(EZ, EZ.T) + 2 * u * np.abs(Zh)       # + rounding of the Hermitian part itself
t_lo = np.trace(Zh).real - np.trace(EZh) * (1 + 1e-9)
t_lo -= SAFE * N * u * np.abs(np.trace(Zh).real)   # rounding of the trace sum
z_hi = np.linalg.norm(Zh) * (1 + SAFE * N * N * u) + np.linalg.norm(EZh) * (1 + 1e-9)
print(f"atoms {n_at} (filterings {int(np.sum(kinds == 'f'))}, products {int(np.sum(kinds == 'p'))}); eps = {float(At['eps']):.1e}", flush=True)
print(f"remainder Z: Tr Z >= {t_lo:.6e} (eps Tr X / 36 = {float(At['eps']) * trX / 36 * 36:.3e} expected for Tr Z), "
      f"||Z||_2 <= {z_hi:.6e}, max error bound {EZh.max():.2e}", flush=True)
ratio = z_hi ** 2 / t_lo ** 2 if t_lo > 0 else np.inf
print(f"purity bound ||Z||_2^2 / (Tr Z)^2 <= {ratio:.10f};  Gurvits-Barnum threshold 1/(N-1) = {1 / (N - 1):.10f}", flush=True)
if t_lo > 0 and ratio <= 1 / (N - 1):
    print(f"CERTIFIED: X = sum of {n_at} nonnegative multiples of filterings of rho^2 and products + a separable remainder, "
          f"so X lies in K_2(6).", flush=True)
else:
    print("NOT CERTIFIED", flush=True)
