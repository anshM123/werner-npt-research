"""Rigorous check of a one-atom (SF-form) certificate for a sum of two atoms.
Claim: X = G(A1, B1) + G(A2, B2) equals lam G(E, F) + (products) + Z with Z separable, i.e. X lies in
Sigma = {lam G(E, F) + separable}.  This is the verifier of compression-notes/code/check_cone_certificate_verify.py with
the target replaced by a sum of two filterings of one copy (each recomputed with its own entrywise rounding bound) and
an assertion that the certificate uses exactly ONE filtering atom.  All parameters are float64 numbers (exact dyadic
rationals); the identity X = sum_i w_i G_i + Z holds by definition of Z.  Every matrix is recomputed in float64 with an
entrywise error bound |fl(s) - s| <= c K u (sum of |products|), u = 2^-53.  If Tr Z >= t_lo > 0, ||Z||_2 <= z_hi and
z_hi^2 <= t_lo^2/(N - 1), N = 9, then Z/Tr Z is in the Gurvits-Barnum ball, so Z is separable.
usage: python verify_bm_cert.py target.npz atoms.npz      (NumPy only)
  target.npz: A1, B1, A2, B2;  atoms.npz: kinds ('f', 'p'), P, Q, x, y, coef, eps"""
import sys, itertools, numpy as np
TFILE, AFILE = sys.argv[1], sys.argv[2]
u = 2.0 ** -53
SAFE = 4.0
D = 3; N = 9


def filt_with_bound(A, B, m=1):
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


T = np.load(TFILE)
X1, E1 = filt_with_bound(T["A1"], T["B1"]); X2, E2 = filt_with_bound(T["A2"], T["B2"])
X = X1 + X2
EX = (E1 + E2) * (1 + 1e-9) + SAFE * 2 * u * (np.abs(X1) + np.abs(X2))      # + rounding of the sum
trX = np.trace(X).real
print(f"target G(A1,B1) + G(A2,B2): Tr X = {trX:.6e}; max entrywise error bound on X {EX.max():.2e}", flush=True)
At = np.load(AFILE, allow_pickle=True)
kinds, coef = At["kinds"], At["coef"]
n_at = len(kinds)
assert np.all(coef >= 0)
n_f = int(np.sum(kinds == 'f'))
assert n_f <= 1, "an SF certificate uses at most one filtering atom"
Zc = X.copy(); EZ = EX.copy(); absum = np.abs(X).copy()
for k in range(n_at):
    if kinds[k] == 'f':
        Gk, Ek = filt_with_bound(At["P"][k], At["Q"][k])
    else:
        Gk, Ek = prod_with_bound(At["x"][k], At["y"][k])
    w = coef[k] * trX / np.trace(Gk).real          # any float w >= 0 is allowed; this one makes Z small
    assert w >= 0
    Zc = Zc - w * Gk
    EZ += w * Ek * (1 + 1e-9)
    absum += w * np.abs(Gk)
EZ += SAFE * (n_at + 4) * u * absum * (1 + 1e-9)
Zh = (Zc + Zc.conj().T) / 2
EZh = np.maximum(EZ, EZ.T) + 2 * u * np.abs(Zh)
t_lo = np.trace(Zh).real - np.trace(EZh) * (1 + 1e-9)
t_lo -= SAFE * N * u * np.abs(np.trace(Zh).real)
z_hi = np.linalg.norm(Zh) * (1 + SAFE * N * N * u) + np.linalg.norm(EZh) * (1 + 1e-9)
print(f"atoms {n_at} (filterings {n_f}, products {int(np.sum(kinds == 'p'))}); eps = {float(At['eps']):.1e}", flush=True)
print(f"remainder Z: Tr Z >= {t_lo:.6e} (eps Tr X = {float(At['eps']) * trX:.3e} expected), ||Z||_2 <= {z_hi:.6e}, "
      f"max error bound {EZh.max():.2e}", flush=True)
ratio = z_hi ** 2 / t_lo ** 2 if t_lo > 0 else np.inf
print(f"purity bound ||Z||_2^2 / (Tr Z)^2 <= {ratio:.10f};  Gurvits-Barnum threshold 1/(N-1) = {1 / (N - 1):.10f}", flush=True)
if t_lo > 0 and ratio <= 1 / (N - 1):
    print(f"CERTIFIED: X = (one filtering of rho) + {n_at - n_f} products + a separable remainder, so X is in "
          f"{{lam G(E, F) + Sep}} (binary merger holds for this pair).", flush=True)
else:
    print("NOT CERTIFIED", flush=True)
