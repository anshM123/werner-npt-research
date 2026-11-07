# The thickness bonus at three copies (Theorem 7.29)

This folder proves Theorem 7.29 of [`whole-problem-routes/README.md`](../../README.md): for every real unit vector y on
three copies,

    λ_min(𝒫_{3/4}^{⊗3}(yyᵀ)) ≥ (1/64)(1 + (98/27)(τ₁ + τ₂ + τ₃)) ≥ ∏_s (¼ + ½τ_s),

where τ_s is the least eigenvalue of the one-copy reduced state of y at copy s. This is Conjecture TC at n = 3. The
proof is computer-assisted. The certificate is a list of rational numbers, and its validity is decided in exact
arithmetic.

## What the certificate is

A unit least eigenvector w_s of each reduced state enters as an extra variable. With t = 1/64, κ = 98/27 and
W = |w₀|²|w₁|²|w₂|², the polynomial

    P(x, y, w) = f(x, y)|y|²W − t|x|²|y|⁴W − tκ|x|²|y|² Σ_s ⟨w_s|ρ^y_s|w_s⟩ W/|w_s|²

is written exactly as an average over permutations of the copies of

    ⟨Ψ|Z₀ ⊗ I_w + Z₁|Ψ⟩ + Tr(Λ(x, y) G(y, w₀)) |w₁|²|w₂|²,    Ψ = x ⊗ y ⊗ y ⊗ w,

with G = |w₀|²ρ^y_0 − ⟨w₀|ρ^y_0|w₀⟩I (positive semidefinite exactly when w₀ is a least eigenvector) and
Λ(x, y)_{ij} = ⟨e_i ⊗ x ⊗ y|L|e_j ⊗ x ⊗ y⟩. At each copy the strands are x, y, y′, w.
- Z₀ acts on (x, y, y′) at every copy.
- Z₁ acts on (x, y, y′, w) at copy 0 and on (x, y, y′) at copies 1 and 2.
- L acts on (e_i, x, y) at copy 0 and on (x, y) at copies 1 and 2.

All three are positive semidefinite rational combinations of Brauer diagrams.

`cert/tc3_certificate.json.xz` stores, for each operator and each spin sector, an integer matrix R and a rational
symmetric matrix Y. The Gram block of the operator on the m = 0 vectors of that sector is R Y Rᵀ. The file also stores
the integer spin bases and a rational right inverse Psi of the map from diagram coefficients to Gram blocks, so the
operators can be decoded to diagram coefficients. It also stores the 14 relations among the 105 four-strand diagrams
(at d = 3 they span only 91 dimensions) and the rational combination of relations used in the identity.

## Verification

- `verify_tc3.py`: the exact verifier, written separately from the scripts that produced the certificate. It
  rebuilds everything except the operators from scratch: the diagrams, the target polynomial, the multiplier terms, the
  group of relabelings and its orbits, the spin bases, and the forced kernels (from exact rational samples of the product
  zeros). It then checks three things. (A) Every relation used is the zero operator. (B) All 190 orbit sums of target
  minus certificate vanish exactly; the quadratic form of a product of diagrams on real Ψ depends only on its orbit
  under exchanging the two x legs, permuting the four y legs, exchanging the two legs of each w_s, and permuting the
  copies (the identity is averaged over copy permutations). (C) Z₀, Z₁ and L are positive semidefinite. Each Gram block Γ
  must annihilate its forced kernel K exactly, and Γ + KKᵀ must be positive definite (exact LDLᵀ). Log:
  `logs/verify_tc3.log`.
- `check_float.py`: decodes the operators in floating point and compares both sides of the identity at random real
  points, with f computed from the definition of 𝒫_{3/4} (relative difference 3.4·10⁻¹⁵, `logs/check_float.log`).

## How it was made

- `exact_alg.py`: exact per-copy data for 2, 3 and 4 strands (diagrams, integer spin bases, the Gram-block map Phi, a
  right inverse Psi, relations), using python-flint.
- `tc3_setup.py`: the target, the orbits, the forced kernels and the complement bases R, then a floating-point SDP
  (CVXPY with Clarabel) at κ = 98/27 that maximises the least eigenvalue of all blocks Y (`logs/tc3_setup.log`).
- `tc3_round.py`: rounds Y, solves for an exact correction in 189 variables (diagonal entries of Y and relation
  coefficients; the 190 orbit equations have rank 189 on the reduced variables, and the remaining one holds exactly),
  and checks the identity and the positivity of every Y exactly (`logs/tc3_round.log`).
- `tc3_export.py`: writes the certificate file.

## Numerics around the theorem (floating point)

- `kappa_scan.py`: the best constant κ that this certificate scheme reaches: 5.0000 at two copies (the exact best
  constant there) and 4.1253 at three copies (`logs/kappa_scan_n2.log`, `logs/kappa_scan_n3.log`).
- `tc3_numerics.py`: the best linear constant is 5 at two copies (maximally entangled vector) and 5 to search accuracy at
  three copies. It also checks the three-copy input that is thin at one copy and has τ = ¼ at the others, every one of
  whose conditional vectors has Schmidt coefficients (1/√2, 1/√2, 0) (`logs/tc3_numerics.log`).

`run_checks.sh` reruns the verifiers and the numerics; `run_checks.sh all` also regenerates the certificate (about
three minutes). Python with NumPy, SciPy and python-flint; CVXPY with Clarabel for `tc3_setup.py` and `kappa_scan.py`;
PyTorch for `tc3_numerics.py`.
