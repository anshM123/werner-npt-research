# One copy carrying the letter 2, every weight vector: the energy theorem and Theorem W

Ansh Mishra and Aryan Senthilkumar, 6 November 2026. Nothing here is refereed.

Notation as in `relative-anchors/README.md` (Lemmas B, N; Theorems T, C; Proposition U1) and
`relative-anchors/code/typeI_G_m1_relax/NOTES.md` (Lemma R1, the relaxed form Ψ, Conjecture R). Throughout m = |T| = 1,
T = {t}, σ = s_t, U the other u copies, w = w_U.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **The energy theorem.** Normalise S_U = 1 and let μ be any probability measure on [0, 1] (in the physical case, mass x_h = s_h at the point x_h). For every σ > 0 the energy E = A − B − 8R_K of Section 2 is strictly positive on every nonzero amplitude with a ⊥ w | **proved**, computer-assisted: an analytic reduction (Section 3) and 45 exact polynomial certificates (`one_letter_certificate.js`) |
| 2 | **The source bound.** ‖R/2‖_F² ≤ (29/72)(B + 8R_K) for physical weights | **proved** by hand (Section 4): an exact decomposition of R/2 and a rational variance identity |
| 3 | **Corollary R.** Ψ ≥ A/30 + (29/30)E > A/30 > 0; so Conjecture R of `typeI_G_m1_relax/` holds for every u ≥ 2 and every positive weight vector (also with every λ replaced by 5/12) | **proved** (rows 1 and 2, Section 4) |
| 4 | **Theorem W.** For every k ≥ 2, every positive weight vector and every copy t, the type I block of (R∞)_{k−1} with T = {t} is ⪰ 0 with kernel exactly span{p ⊗ f_h : h ∈ U} | **proved**: u = 1 is Proposition U1; for u ≥ 2, Corollary R with Lemma R1 gives G ⪰ G_rel ≻ 0, and Theorem C gives the block |
| 5 | Checks written separately from the certificate, from our own objects (`psi_exact.py` of `typeI_G_m1_relax/`) | **exact** (Section 5) |
| 6 | The M-matrix structure (Z), (D) of Conjecture M (`typeI_G_u3/`); type I blocks of (R∞)_{k−1} with \|T\| ≥ 2 and \|U\| ≥ 3 at general weights; the middle conditions; Conjecture SF; the endpoint itself | **open** |

Theorem W replaces Theorems EU and TL (one or two U-weight levels) by every weight vector. Nothing here says that
Conjecture 12 of `single-excitation-range/` holds in general, or anything about the NPT problem beyond these blocks.

## 1. The relaxed form in measure coordinates

By Lemma R1 of `typeI_G_m1_relax/`, after the exact elimination of Q and the splitting X = X′ + a wᵀ + w cᵀ of
Theorem C (a ⊥ w), the relaxed form is

    Ψ(a, c, ρ, ζ) = A(a, c, ρ, ζ) − ¼ Σ_{r,h} R_rh² / λ_rh,   λ_rh = 5/12 (r ≠ h),  λ_hh = 3/2 − s_h/(σ + S/2 + s_h) > 5/12,

where A is the q-eliminated form at X = Ω = a wᵀ + w cᵀ (X′ = 0) and R = Π(σ𝒢 + ℋ)Π is the projected source of the
hard part. Normalise S = Σ_h s_h = 1, put x_h = s_h and E[q] = Σ_h x_h q_h, and use the amplitudes A_h = a_h/w_h,
C_h = c_h/w_h, Z_h = ζ_h/w_h (so E A = 0) and ρ. All expressions below are then averages over the probability measure
μ = Σ_h x_h δ_{x_h}. The same formulas make sense for any probability measure on [0, 1] ("formal" measures, with
locations and masses unrelated); the energy theorem holds on that larger domain.

## 2. The decomposition of the source and the energy

Put K = σ + ½, F = E[1/(K + x)], G = 2 + F, p = E C, k = K − σ/2 = (σ + 1)/2, b = σρ + √σ E Z and
L = 2 + σ/(k + σ) + E[1/(k + x)]. Define the field

    V = σ[(A + C)/(K + x) − p/(2K(K + x))] + √σ Z/(2(k + x)) − b/(4k(k + x)),

and the three rank-one terms P_j ⊗ Q_j (centered: P_j = PP_j − E PP_j, Q_j = QQ_j − E QQ_j) with

    PP_0 = σ/(G(K + x)),                             QQ_0 = x A/(K + x),
    PP_1 = (σ/G)(x C/(K + x) + p/(2(K + x))),        QQ_1 = 1/(K + x),
    PP_2 = [√σ x Z/(2(k + x)) + b/(4(k + x))]/L,      QQ_2 = 1/(k + x).

In each term one factor is fixed and the other is linear in the amplitudes. Then, as u × u matrices,

    R/2 = Π diag(x ∘ V) Π − K_rank,    K_rank = Σ_j (w ∘ P_j)(w ∘ Q_j)ᵀ,

with w_h = √x_h (checked exactly against `psi_exact.py`, Section 5). Put

    B = E[(x + ¼)V²] + (E V)²/4,    R_K = ‖K_rank‖_F² = Σ_{j,l} E[P_j P_l] E[Q_j Q_l],    E = A − B − 8R_K.

## 3. The energy theorem

**Theorem.** For every probability measure μ on [0, 1], every σ > 0 and every nonzero amplitude with E A = 0,
E = A − B − 8R_K > 0 (real or complex amplitudes).

*Proof (computer-assisted in the last step).* Write L = √σ(A + C₀)/2, D = √σ(A − C₀)/2 with C = C₀ + p, Z = Z₀ + z,
E A = E C₀ = E Z₀ = 0, and put q = ρ − √σ z, h = (√σ ρ + z)/(σ + 1) + √σ p/(2K), δ = (3σ + 1)/2, f = 1/(K + x),
g = 1/(k + x), α = 1/(2 + E f), γ = 1/(2 + σ/δ + E g).
1. **A coherent square.** The shift Z₀ + z ↦ Z₀ + z + η, ρ ↦ ρ + √σ η leaves V and the centered rank part unchanged
   and adds |η|² to A with cross term 2 Re(η̄ h). Hence E = |h|² + E_red(L, D, Z₀; p, q) exactly.
2. **The rank part with the means kept.** By the orthogonality of the symmetric and antisymmetric parts of the rank
   terms, Cauchy–Schwarz for covariances and the variance bounds Var f ≤ 1/(4K²(K + 1)²),
   Var g ≤ 1/(4k²(k + 1)²), one gets 8R_K ≤ E[U_R(x)] with an explicit quadratic U_R in which the means p, q stay
   inside the squares.
3. **The collective terms.** The zero-mean collective terms are bounded by weighted Cauchy–Schwarz with the same
   resolvent variances and a few elementary polynomial inequalities (each a sum of squares plus a positive constant).
   Together with the exact local part this gives a pointwise lower bound
   E_red ≥ ∫ (v, y)* P(x; α, γ, F, B) (v, y) dμ with v = (L, D, Z₀), y = (p, q) and an explicit 5 × 5 matrix P
   (a constant shift of the mean/fluctuation source is allowed because E v = 0).
4. **Nine vertices.** For fixed (σ, x), P is jointly matrix-concave in (α, γ, F, B): everything is affine except
   negative semidefinite quadratic terms. The physical values lie on the convex curves α = 1/(2 + F),
   F ∈ [1/(K + 1), 1/K], and γ = 1/(2 + σ/δ + B), B ∈ [1/(k + 1), 1/k]. Each curve lies in the triangle formed by
   its chord and its two endpoint tangents. So the physical point is a convex combination of the nine products of
   triangle vertices, and P is bounded below by the same combination of the nine vertex matrices.
5. **Certificates.** After the positive congruence diag(d, d, e, d/√σ, e/√σ) (d = K + x, e = k + x) each vertex matrix,
   with its denominator cleared, is an integer polynomial matrix in σ and x. In the order (D, Z, q, L, p) its five
   leading principal minors have only nonnegative coefficients in σ and x, and each has a nonzero coefficient
   polynomial free of x. So all nine families are positive definite for σ > 0, x ∈ [0, 1] (45 minors, 20 406 positive
   coefficients, none negative; `one_letter_certificate.js`, `logs/one_letter_certificate.log`).

Hence ∫ (v, y)* P (v, y) dμ > 0 for nonzero (v, y), and adding |h|² proves the theorem. ∎

The program regenerates the nine families and the 45 minors in exact integer arithmetic, compares 2 025 matrix entries
and 405 determinants with a separate rational implementation, checks the reduction of steps 1–3 against the energy at
fixed measures, and checks the triangle and concavity step. Doubling the source-centering correction gives 280
negative coefficients in one family, and an explicit point where the altered matrix is indefinite.

## 4. The source bound and Corollary R

Let N_D = ‖Π diag(x ∘ V) Π‖_F² = E[x(1 − 2x)V²] + (E[xV])². With b = 13/72 and d = 11/2592,

    (5/18) B − N_D = 2 Var((x − b)V) + (E[(x − 2b)V])² + d (E V² + (E V)²) ≥ 0,

an exact identity. So ‖R/2‖_F ≤ √N_D + √R_K, and Cauchy–Schwarz with weight 9/20 gives
‖R/2‖_F² ≤ (29/20)(5/18)B + (29/9)R_K = (29/72)(B + 8R_K). Since λ_rh ≥ 5/12,

    Ψ ≥ A − (12/5)‖R/2‖_F² ≥ A − (29/30)(B + 8R_K) = A/30 + (29/30)E > A/30 > 0

on every nonzero amplitude (A = E + B + 8R_K > 0). The same holds with every λ replaced by 5/12. This is Corollary R,
and with Lemma R1 and Theorem C it proves Theorem W.

## 5. Checks written separately

All exact; the Python checks use `psi_exact.py` of `typeI_G_m1_relax/`, written before the energy theorem.
- `crossbase.py` (`logs/crossbase_s5.log`): the energy model of the programs (`relax_comparison.js`, the same function
  as in the certificate) equals our base A and our ‖R/2‖_F² at 16 random rational weight vectors and amplitudes,
  u = 2…5.
- `bridge.py` (`logs/bridge_s1.log`): at 30 random rational weight vectors with Σ s_h = 1 (u = 2…6, σ from 1/40 000
  to 40 000): the decomposition R/2 = Π diag(x ∘ V) Π − K_rank holds exactly, the variance identity holds exactly, the
  source bound holds, E is positive definite as a quadratic form on {a ⊥ w}, and Ψ − A/30 − (29/30)E ⪰ 0.
- `test_A30.py` (`logs/test_A30_s1.log`): Ψ − A/30 is positive definite at 72 rational weight vectors (u = 2, …, 6 and 8; small
  and large σ; one dominant free weight), also with every λ = 5/12. Control (`logs/test_A30_control_half.log`): the
  stronger Ψ − A/2 is not positive definite at several of them.
- `relax_comparison.js` (`logs/relax_comparison.log`): the energy model against the relaxed form at 12 rational
  configurations, entry by entry.
- `stress_profiles.py` (`logs/stress_profiles.log`): the reduction of steps 1–4 at 40 further random measures
  (physical and formal), as full matrix inequalities.

## Files

| File | Content |
|---|---|
| `one_letter_certificate.js` | the certificate of the energy theorem (Section 3), with its own checks and the negative control |
| `relax_comparison.js` | the energy model, and its comparison with the relaxed form of Lemma R1 |
| `run_js.py` | runs either program with Node (no packages); prints the verdict as JSON |
| `stress_profiles.py` | the reduction checked at 40 further random measures |
| `crossbase.py`, `bridge.py`, `test_A30.py` | the checks of Section 5 |
| `run_all.sh` | every step, with logs in `logs/` |
