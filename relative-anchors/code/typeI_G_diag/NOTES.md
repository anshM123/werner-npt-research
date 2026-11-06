# The diagonal blocks of G at every weight vector (Theorem DG)

Ansh Mishra and Aryan Senthilkumar, 6 November 2026. Nothing here is refereed.

Notation as in `relative-anchors/README.md` (Theorem C and the Schur complement G of a type I block of (R∞)_{k−1}) and
`../typeI_G_m1_all/NOTES.md` (Theorem W, its energy theorem and source bound). T is the set of copies carrying the
letter 2 (m = |T| ≥ 1), U the free copies (u = |U| ≥ 1), σ = Σ_{t∈T} s_t, S = Σ_{h∈U} s_h, θ_t = s_t/σ.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Theorem DG.** For every k, every positive weight vector and every t ∈ T, the diagonal block of G at t satisfies Ĝ_tt − (3/4)(1 − θ_t) I ≻ 0. Numerically 3/4 is the best constant (Section 4) | **proved**, computer-assisted (Section 2; `diag_certificate.js`) |
| 2 | **Lemma GG (weighted Gershgorin).** If also ‖Ĝ_ij‖_op ≤ (3/4)√(θ_iθ_j) for all i ≠ j in T, then G ≻ 0, and with Theorem C the type I block is ⪰ 0 with the kernel of Theorem T | **proved** (Section 3) |
| 3 | The off-diagonal bound of row 2 | **open**; numerically the ratio ‖Ĝ_ij‖_op/√(θ_iθ_j) is at most 0.258 (Section 4) |
| 4 | Type I blocks of (R∞)_{k−1} with \|T\| ≥ 2 at general weights; the middle conditions; Conjecture 12 in general; the endpoint | **open** |

For m = 1 the bound of row 1 reads Ĝ_tt ≻ 0, which is Theorem W. A first version of this note (commit 64a3707) had the
constant 2/3; the same certificate goes through with 3/4. Nothing here says that Conjecture 12 of
`single-excitation-range/` holds in general, or anything about the NPT problem beyond these blocks.

## 1. The frame

Ĝ = PᵀGP, where P acts on the letter-2 states of each copy t by ρ_t ↦ ρ_t and ζ_t = c z_t w_U/|w_U| + ζ′_t, with
ζ′_t ⊥ w_U and c = √(1 + S/σ) (`../typeI_G_u3/numG.py`, functions `frame` and `Ghat`). In these coordinates the letter
norm of the energy formulation is ρ_t² + E|Z_t|² − |EZ_t|²/(t + 1) (normalising S = 1, t = σ), so the identity block of
the statement is this norm.

## 2. Proof of Theorem DG

The argument extends the energy theorem of `../typeI_G_m1_all/` to one active letter-2 copy among spectators.
1. **Formal statement.** Normalise S = 1 and let μ be any probability measure on [0, 1]. Fix a = s_t ≤ σ. Only copy t
   carries letter-2 data (ρ, Z); every pair form with a spectator vanishes, and the gauge projection costs exactly
   −(1 − θ)E|ℓ|²/4. The claim is that, on every nonzero amplitude with E A = 0,

       A − B − 8R − H(1 − θ)(ρ² + E|Z|² − |EZ|²/(t + 1)) > 0   with H = 3/4,

   for the soft form A, source B and centred rank norm R of the energy theorem, restricted to these data.
2. **Point matrix.** After the coherent square, the rank bound with the means kept and the collective bounds of the energy
   theorem, the left side is at least E[(q, m)* P₆(s) (q, m)], with q = (L, D, Z₀), m = (p, r, h) and an explicit
   symmetric 6 × 6 matrix P₆(s). Centring the ℓ-source as well is needed for this.
3. **Nine vertices.** For fixed (t, a, s), P₆ is jointly matrix-concave in the four resolvent statistics
   (α, γ, F, B). The curves α = 1/(2 + F) and γ = 1/(2 + a/δ + B) lie in their chord-and-tangent triangles, so P₆ is
   bounded below by a convex combination of nine vertex matrices.
4. **Certificate.** After the congruences diag(1, 1, √(t/a), 1/√t, √(t/a), √(t/a)) and diag(d, d, e, d, e, e)
   (d = K + s, e = k + s), every vertex matrix is a rational matrix in (a, b = t − a, s) with positive denominators. In
   the order (h, D, Z₀, r, L, p) its six leading principal minors have numerators with only positive coefficients, and
   each has a positive monomial free of b and s: 54 minors, 533 525 positive coefficients, none negative. So every
   vertex matrix is positive definite for a > 0, b ≥ 0, s ≥ 0.
5. **Physical weights.** With the source bound ‖R/2‖_F² ≤ (29/72)(B + 8R_K) and λ ≥ 5/12 of `../typeI_G_m1_all/`, the
   form left after eliminating the hard residual is at least the formal energy of step 1 plus
   H(1 − θ)(ρ² + E|Z|² − |EZ|²/(t + 1)) plus a positive multiple of B + 8R. Minimising over the other variables gives
   row 1. ∎

The program `diag_certificate.js` implements the field and source matrices of step 2 and the bounds used there
(resolvent variance bounds with V_F = 1/(4K²(K + 1)²) and V_G = 1/(4k²(k + 1)²), the rank bound with one active
Z-rank, the row-gauge bound), and regenerates the nine families and the 54 minors in exact integer arithmetic. Its own
guards (all PASS, `logs/diag_certificate.log`): 1 296 matrix entries and 216 determinants against a separate rational
implementation; 6 original-energy lower bounds (9 to 12 dimensions); 6 full decompositions (45 coefficient identities
each); 18 concavity gaps; 4 860 identities reducing a = t to the energy theorem of Theorem W.

**Controls.** Doubling the centring constant gives a point matrix with inertia (5, 1) at t = a = 10 000, s = 0 (in the
program). With H = 4/5 the program fails (negative coefficients in the third minor of the first family), and with the
earlier constant H = 2/3 it passes (`control_H.py`, `logs/control_H.log`).

## 3. Lemma GG

Let x = (x_t)_{t∈T}. If Ĝ_tt ⪰ (3/4)(1 − θ_t)I (strictly) and ‖Ĝ_ij‖ ≤ (3/4)√(θ_iθ_j) for i ≠ j, then

    x*Ĝx ≥ (3/4)[Σ_t (1 − θ_t)|x_t|² − Σ_{i≠j} √(θ_iθ_j)|x_i||x_j|] = (3/4)[Σ_t |x_t|² − (Σ_t √θ_t |x_t|)²] ≥ 0

by Cauchy–Schwarz and Σ_t θ_t = 1, with strict inequality for x ≠ 0 from the strict diagonal bound. ∎

## 4. Checks written separately

- `diag_check.py` (`logs/diag_check.log`): our own G (`../typeI_G_u3/numG.py`, floating point, built from the closed
  form of `cf.py`), 729 diagonal blocks at 240 random weight vectors with m = 2, …, 4 and u = 3, …, 7 (spread, wide,
  one dominant letter-2 weight, all letter-2 weights small): min λ_min(Ĝ_tt)/(1 − θ_t) = 0.750012, approached as
  θ_t → 0. Runs at u = 1 and u = 2 gave 0.7512 and 0.7505. So 3/4 appears to be the best constant (we have not
  proved that the infimum equals 3/4).
- `offdiag_ratio.py` (`logs/offdiag_ratio.log`): at the same kind of weight vectors, 821 off-diagonal blocks have
  ‖Ĝ_ij‖_op/√(θ_iθ_j) ≤ 0.258. So the route of Lemma GG has a lot of room (0.258 against 3/4), but the off-diagonal
  bound is not proved.

## Files

| File | Content |
|---|---|
| `diag_certificate.js` | the certificate of Theorem DG (nine families, 54 minors) with its guards and the centring control |
| `control_H.py` | the same program with H = 4/5 (fails) and H = 2/3 (passes) |
| `run_js.py` | runs it with Node (no packages); prints the verdict as JSON |
| `diag_check.py`, `offdiag_ratio.py` | the checks of Section 4 |
| `run_all.sh` | every step, with logs in `logs/` |
