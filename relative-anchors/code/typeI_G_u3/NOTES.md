# G ≻ 0 for the type I blocks of (R∞)_{k−1} with |U| ≥ 3: the case k = 4, a positivity pattern, and the obstruction to uniformity

Ansh Mishra and Aryan Senthilkumar, 4–5 November 2026. Nothing here is refereed. Notation as in
`relative-anchors/README.md` (Lemmas B, N; Theorem T; Propositions T2, U0, U1; Theorems C, U2),
`relative-anchors/code/typeI_u2/NOTES.md` (Parts A, B) and `relative-anchors/code/typeI_G_u2/NOTES.md`.
U = copies without a letter 2, u = |U|; T = copies carrying one letter 2, m = |T|; s_x = w_x²; σ = s_T, S = s_U.
All forms are in units of 2Y.

## Status

| Statement | Status |
|---|---|
| **Theorem U3₁.** For k = 4, \|T\| = 1, \|U\| = 3 and all positive weights, G ≻ 0. Hence every type I block of (R∞)₃ with \|U\| = 3 is ⪰ 0 with kernel exactly span{p ⊗ f_h : h ∈ U} | **proved**, computer-assisted, exact: closed form over Q(s); every leading principal minor of G_z is a ratio of polynomials with only positive coefficients; Sylvester's criterion. A second exact certificate: G_z is a symmetric Z-matrix and G_z v > 0 for an explicit v > 0 (Section 4b). No floating point and no Bernstein subdivision is used |
| **Corollary.** For k = 4 and all positive weights the whole condition (R∞)₃ = (R∞)_{k−1} holds (every block ⪰ 0, kernels as in the README) | **proved** (Theorem U3₁ with Theorem T, Proposition T2, Propositions U0, U1, Theorem U2 and Theorem C) |
| Closed form of G for (u, m) = (3, 1) equals direct exact elimination of the whole block | **exact** check: 12 rational weight vectors × 2 gauges |
| The same certificate for (u, m) = (1, 1), (2, 1), (1, 2), (2, 2) | **exact**; an independent re-proof of Proposition U1 and Theorem U2 for m ≤ 2 |
| The same certificate for (u, m) = (1, 3) | **exact** (G ≻ 0 there is also Proposition U1) |
| **Conjecture P.** For every (u, m), all leading principal minors of G_z (order ρ_0, z_{0,·}, ρ_1, z_{1,·}, …) are ratios of polynomials in the squared weights with only positive coefficients | **open**; true in all six cases with u + m ≤ 4; it would give G ≻ 0 for every k |
| For m = 1: G_z is a symmetric Z-matrix and G_z v > 0 for v = (4s_t(s_t + S), 3(2s_t + S), …, 3(2s_t + S)), every entry of den·G_z v having only positive coefficients; so D_v G_z D_v is strictly diagonally dominant | **exact** for u ≤ 3; **numerical** for u ≤ 8. As Conjecture M (Section 4b) it would prove the case m = 1 for every u |
| (u, m) = (3, 2) and (4, 1) (the two open cases at k = 5) by the same pipeline | **not completed** in the time (Section 5) |
| u = 3 for every m by the U2 pipeline ((A₃), (B₃)) | **open**; the five-parameter certificate (A₃) is beyond the U2 code in time and memory (Section 5) |
| Negative controls for (u, m) = (3, 1) | **exact** witnesses: G_z − W/10 (W the identity of the rescaled frame) and G_z − e_ρe_ρᵀ/1000 are indefinite at explicit rational points |
| The brief's merged comparison as an exact identity, or as a one-sided bound | **false**: the (ρ, w_U)-block differs from the merged u = 1 block by an indefinite matrix of size up to 0.37 (numerical, Section 4) |
| G ≻ 0 and the weighted Gershgorin constants of U1/U2 for u ≤ 7, m ≤ 4 | **numerical** (Section 6): sup A = 0.30, inf B = 3/4, dominance ratio ≤ 0.293 |
| G ≻ 0 uniformly in u (every u ≥ 3, every m), i.e. the rest of (R∞)_{k−1} for k ≥ 5 | **open** |

Nothing here says anything about Conjecture 12 or the endpoint problem beyond these blocks. The conditions
(R∞)_j with 1 ≤ j ≤ k − 2 remain open for general weights, also at k = 4.

## 1. Setting

Fix a type I block of (R∞)_{k−1}: T ≠ ∅ carries one letter 2, U = [k] ∖ T ≠ ∅. Its states are P_{r,h}, R_t, Z_{t,h},
Q (README, Remark after Proposition U0); the closed form is `cf.py` (copied from `typeI_G_u2/`; `check_cf.py` compares
it exactly with the Lemma B builder and in floating point with the direct construction). By Theorem C and Section A3 of
`typeI_u2/NOTES.md`: if the coordinates P_{g,h} (h ∈ U) are deleted for one fixed g ∈ U (a complement of the kernel
span{p ⊗ f_h}, whose vectors have P_{g,h}-coordinate w_g c_h), the principal block on all remaining P and on Q is positive
definite, and the type I block is ⪰ 0 with kernel exactly span{p ⊗ f_h} iff the Schur complement G on the letter-2
states (ρ_t, ζ_{t,h}) is positive definite. G does not depend on g (checked again below).

**Scaled coordinates.** Put X_{r,h} = w_r w_h x_{r,h} and ζ_{t,h} = w_t w_h z_{t,h}, with ρ and q unscaled. Modulo the
part common to all states, the states carry the parity vectors e_r + e_h (P_{r,h}), e_t + e_h (Z_{t,h}) and 0 (R_t, Q)
in (Z/2)^k, so in these coordinates every entry of the block is even in every w_i, i.e. a rational function of the
squared weights s_i. The scripts check this exactly (an odd exponent raises an error). Let G_z be the Schur complement in
these coordinates. G_z = Dg G Dg with Dg = diag(1, w_t w_h, …) invertible, so G ≻ 0 ⟺ G_z ≻ 0. The entries of G_z are
homogeneous of degree 0 (ρρ), 1 (ρz) and 2 (zz), so G_z(λs) = D_λ G_z(s) D_λ with D_λ = diag(1, λ, …): positivity is
invariant under scaling all weights.

## 2. Theorem U3₁ (k = 4: u = 3, m = 1)

**Theorem U3₁.** Let k = 4, T = {t} one copy, U the other three copies, and all weights positive. Then G ≻ 0.

*Proof (computer-assisted, exact; `closed_m1.py`, `build_m1.py`, `prove_m1.py`).* Number the copies t = 0 and
U = {1, 2, 3}.
1. **Closed form.** `cf.TypeI` is built over Q(w_0, …, w_3) (python-flint), rescaled as in Section 1 and converted to
   Q(s_0, …, s_3). The gauge coordinates P_{1,h} are deleted and the other 10 variables (P_{0,h}, the P_{r,h} with
   r ∈ {2, 3}, and Q) are eliminated by exact Gaussian elimination over Q(s_0, …, s_3) (under a minute). Every pivot is
   positive on the open orthant, since the eliminated block is positive definite there (Theorem C), so the rational
   function obtained is the Schur complement at every positive weight vector. G_z is a symmetric 4 × 4 matrix in the
   coordinates (ρ, z_1, z_2, z_3) with one common denominator
   den = 8 · f₇ · f₁₇ · (s_0 + S)² · e₄²,
   where f₇, f₁₇, e₄ are irreducible of degrees 7, 17 and 4 (120, 1140 and 35 terms). Each of these polynomials has only
   positive coefficients. We did not identify f₇, f₁₇, e₄ further.
2. **Leading minors.** With Nt = den·G_z (integer coefficients after one positive scaling), the leading principal minors
   m_r = det Nt[:r, :r] are computed by Laplace expansion. The irreducible factors of den are pulled out by exact
   division at each step, and every expansion is checked against a direct determinant at a rational point. Each m_r is a
   product of powers of the factors of den times a core c_r:

   | r | exponents of (f₇, f₁₇, s_0 + S, e₄) | degree of c_r | terms of c_r | negative coefficients |
   |---|---|---|---|---|
   | 1 | (0, 0, 0, 0) | 34 | 7 769 | 0 |
   | 2 | (1, 1, 2, 0) | 44 | 14 188 | 0 |
   | 3 | (2, 2, 4, 0) | 54 | 23 374 | 0 |
   | 4 | (3, 3, 6, 2) | 56 | 22 099 | 0 |

3. **Signs.** minor_r(G_z) = m_r/den^r = c_r · ∏_i f_i^{x_ir − r E_i} / 8^r (E_i the exponents in den). Every polynomial
   in this expression has only positive coefficients, so it is positive whenever all s_i > 0. By Sylvester's criterion
   G_z ≻ 0, hence G ≻ 0, for all positive weights. ∎

**Second proof (`check_zpattern.py`, `check_mvec.py`).** Every off-diagonal numerator of G_z has only negative
coefficients and the common denominator only positive ones, so G_z is a symmetric Z-matrix on the open orthant. With
v = (4s_0(s_0 + S), 3(2s_0 + S), 3(2s_0 + S), 3(2s_0 + S)) > 0, every entry of den·G_z v is a polynomial with only
positive coefficients (8 435, 7 770, 7 770, 7 770 terms), so G_z v > 0 entrywise. Then D_v G_z D_v (D_v = diag(v)) has
nonpositive off-diagonal entries and positive row sums v_i(G_z v)_i, so it is strictly diagonally dominant with positive
diagonal, hence positive definite (Gershgorin). So G_z ≻ 0. ∎

**Corollary.** For k = 4 and all positive weights, (R∞)₃ holds: the letter-1 block by Theorem T, the type II blocks by
Proposition T2, and the type I blocks with |U| = 0, 1, 2, 3 by Proposition U0, Proposition U1, Theorem U2 and Theorem
U3₁ with Theorem C. (At k = 4 a type I block has |T| = 4 − |U| ≥ 1, so these are all of them.)

**Checks and controls.**
- `check_closed_m1.py`: at 12 random rational weight vectors the closed form equals the Schur complement obtained by
  exact elimination of the whole block over Q (Fractions), for the gauges P_{1,h} and P_{2,h}: 24 exact equalities. So
  the closed form is right, and G does not depend on the gauge.
- `build_m1.py`: every Laplace expansion equals a direct exact determinant at a rational point.
- `controls_m1.py`: (P) the exact LDLᵀ pivots of G_z at s = (1/2, 1/4, 1/6, 1/12) are positive. (C1) Let
  W = diag(1, s_0(diag(s_U) − s_U s_Uᵀ/(s_0 + S))), the identity of the rescaled frame of `typeI_G_u2/` written in
  z-coordinates. G_z − W/10 is indefinite at s = (1, 1/300, 1/300, 1/300): there is no uniform spectral gap in that frame;
  G degenerates as S → 0, as in Proposition U1. (C2) G_z − e_ρe_ρᵀ/1000 is indefinite at s = (1, 10⁻⁴, 10⁻⁴, 10⁻⁴). So
  the certificate is exact where it has to be: near S → 0 there is no room.

## 3. The same certificate for small (u, m), and Conjecture P

`closed_gen.py u m` computes G_z for T = {0, …, m − 1}, U = {m, …, m + u − 1} by the same method (for m ≥ 2 the
entries have different denominators; the common denominator is their lcm). `check_closed_gen.py u m` compares it
exactly with direct elimination (6 rational points, two gauges). `minors_pos.py u m` computes every leading principal
minor of G_z with the denominator factors pulled out (all minors of the top r rows, by expansion along row r) and checks
the coefficient signs. It prints PASS only if the expansions check, the constant of the denominator is positive, and
every denominator factor and every leading-minor core has only positive coefficients. Then G ≻ 0 for that (u, m) and
all positive weights, by Sylvester's criterion.

| (u, m) | k | size of G_z | variables | cores (degree, terms) of the leading minors | all coefficients positive | time |
|---|---|---|---|---|---|---|
| (1, 1) | 2 | 2 | 2 | (8, 8), (10, 7) | yes | < 1 s |
| (2, 1) | 3 | 3 | 3 | (15, 135), (23, 251), (25, 230) | yes | < 1 s |
| (1, 2) | 3 | 4 | 3 | (8, 44), (14, 89), (26, 322), (32, 430) | yes | < 1 s |
| (2, 2) | 4 | 6 | 4 | (16, 951) … (64, 32 173) | yes | about 1 min |
| (1, 3) | 4 | 6 | 4 | (8, 164) … (56, 23 413) | yes | about 40 s |
| (3, 1) | 4 | 4 | 4 | (34, 7 769) … (56, 22 099) | yes | about 1.5 min |

The cases u = 1 and u = 2 are Proposition U1 and Theorem U2 (for every m); for u + m ≤ 4 the table re-proves them by a
different certificate. The case (3, 1) reproduces Theorem U3₁ through the generic code.

**Conjecture P.** For every u, m ≥ 1, every leading principal minor of G_z in the order (ρ_0, z_{0,·}, ρ_1, z_{1,·}, …)
is a ratio of polynomials in the squared weights with only positive coefficients.

Conjecture P implies G ≻ 0 for every type I block, hence, with the README results, the whole condition (R∞)_{k−1} for
every k and all weights. We have no proof; it holds in all six cases with u + m ≤ 4. The pattern suggests a
subtraction-free (combinatorial) expansion of these minors, as for the principal minors of a weighted Laplacian plus a
positive diagonal. For m = 1 Section 4b makes this precise for u ≤ 3; for m ≥ 2 we have no such structure.

## 4. The suggested route (comparison with the merged u = 1 system) and why it does not give a uniform proof

The brief proposed: split ζ_t into its w_U-component z_t and the part ζ′_t ⊥ w_U; compare the (ρ, z) part with the
u = 1 system in which the U-copies are merged into one copy of weight √S; bound the ⊥ part; close with the u = 1
margins. Findings (`ex1.py`, `ex_equal.py`, `sweep_u.py`; floating point, frame (ρ_t, c z_t, ζ′_t), c = √(1 + S/σ)):

1. **The ⊥ part is stiff and weakly coupled.** Over the sweep (u ≤ 7, m ≤ 4) the ζ′-block has eigenvalues in
   [0.719, 1.000] and its coupling to (ρ, z) is at most 1.8·10⁻² in absolute value. At equal U-weights the coupling
   vanishes exactly (the ζ′_t span copies of the standard representation of the permutations of U) and the ζ′-block is
   about 0.72 to 0.79. We did not derive these values; a uniform bound for this block would still need the (C_u)
   inverse, which enters every entry of G.
2. **The (ρ, z) part is not the merged u = 1 block, and not on one side of it.** The difference has entries up to 0.37
   in size (largest for u = 7 and spread weights) and eigenvalues between −0.37 and +0.07. It is present at equal
   U-weights, so it is not an effect of the spread of the U-weights alone. So neither an exact identity nor a one-sided
   comparison holds, and the margins of Proposition U1 (dominance ratio 0.29, diagonal constant 3/4) cannot simply be
   transferred.
3. **Its size depends on u.** At equal U-weights (`ex_equal.py`) the smallest eigenvalue of the (ρ, z) block decreases
   with u in every example. Example: one T-copy with s_t = 0.004 and S = 1 gives 0.579 (u = 1), 0.530, 0.505, 0.489,
   0.478, 0.470, 0.464, 0.460 (u = 8) and 0.448 (u = 12). The decrements (0.049, 0.026, 0.016, 0.011, 0.008, 0.006,
   0.005) shrink a little faster than 1/u, so a positive limit as u → ∞ is plausible, but we have not established it. A
   uniform proof has to control this u-dependence; the u = 1 constants alone do not.
4. **What causes the difference.** Two kinds of U-directions have no counterpart in the merged system:
   - the trace direction Π = I − w_U w_Uᵀ/S of the hard part X′ of Theorem C's splitting. It is invariant under the
     permutations of equal-weight U-copies (so it survives at equal weights), and it couples to ρ at order one through
     ½(tr X + Σ_t ρ_t)²;
   - the soft modes a w_Uᵀ and w_U cᵀ (a, c ⊥ w_U), with eigenvalue of order σ/S (Part A of `typeI_u2/NOTES.md`). In
     the rescaled frame an error of order σ in their treatment becomes an error of order one in the z-entries, so they
     have to be treated exactly to first order. The quantitative form of Theorem C (Corollary Q, constant 9·10⁻⁵ against
     about 0.42 numerically) is far too weak to be used perturbatively.
5. **Consequence.** A uniform argument has to control C_tᵀC^{−1}C_t′ (Section 1 of `typeI_G_u2/NOTES.md`), the inverse of
   the whole (C_u) block applied to the couplings, uniformly in u and exactly on the soft directions. We have no such
   estimate. The relaxation of Section B5 of `typeI_u2/NOTES.md` (hard part relaxed, trace direction kept) is numerically
   positive, but with a gap of only about 5·10⁻⁴ σ/|s| and with couplings indexed by U-coordinates and T-copies at once,
   which we did not reduce to finitely many parameters.

### 4b. m = 1: an M-matrix structure, and a route to every u

For m = 1 (T = {t}), and u = 1, 2, 3, exactly (`check_zpattern.py`, `check_mvec.py`):
- **(Z)** every off-diagonal entry of G_z (ρ–z_h and z_h–z_h′) is minus a ratio of polynomials with only positive
  coefficients, and the diagonal entries are positive such ratios. So G_z is a symmetric Z-matrix on the open orthant.
- **(D)** with v = (4s_t(s_t + S), 3(2s_t + S), …, 3(2s_t + S)), every entry of den·G_z v has only positive
  coefficients. So G_z v > 0 entrywise.

(Z) and (D) give G ≻ 0 at once: D_v G_z D_v has nonpositive off-diagonal entries and positive row sums v_i(G_z v)_i, so it
is strictly diagonally dominant with positive diagonal. This is the second proof of Theorem U3₁ in Section 2, and it
re-proves the cases (1, 1), (2, 1). It also explains Conjecture P for these cases: D_v G_z D_v is a weighted graph
Laplacian with subtraction-free edge weights plus the subtraction-free diagonal diag(v_i(G_z v)_i), so by the
all-minors matrix-forest theorem every principal minor is a sum of products of edge weights and root weights.

The vector v is (ρ, ζ) = (1, β w_t w_U) up to a positive factor, β = 3(2s_t + S)/(4s_t(s_t + S)): the letter-2 row plus
a multiple of the merged direction ζ ∝ w_U. We found β numerically: the admissible window for β is
max_h |G_z[h, ρ]|/r_h < β < G_z[ρ, ρ]/|Σ_h G_z[ρ, h]| (r_h the row sums of the z-block); in units of 1/s_t its ends go
from (1, 2) as S/s_t → 0 to about (0.5–0.55, 1) as S/s_t → ∞ (`ex_beta2.py`, u ≤ 5), and β = 1/s_t approaches both
ends in the two limits (`ex_beta.py`), so a weight-dependent β is needed.

**Conjecture M.** For m = 1, every u ≥ 1 and all positive weights, (Z) and (D) hold (as inequalities on the open
orthant). Numerically they hold for u ≤ 8 (`ex_mvec.py`: 150 weight vectors per u ≤ 5 and 60 for u = 6, 7, 8, spreads
up to e^{±7}, one U-weight blown up or all shrunk): every off-diagonal entry is negative and every ratio
(G_z v)_i/(G_ii v_i) is positive. In the most degenerate samples (a tiny U-weight) these ratios go down to about 10⁻⁹,
close to the resolution of the floating-point computation, so those samples are weak evidence. The window of β has ratio
of ends at most 0.56 for u ≤ 7 (`ex_mnum.py`).

Conjecture M would prove G ≻ 0 for m = 1 and every u, i.e. every type I block of (R∞)_{k−1} with exactly one copy
carrying the letter 2, for every k. What a proof needs: (D) is the action of G on one explicit vector, which only
requires the (C_u) block applied to the source generated by ρ and the merged direction; that source lies in the algebra
generated by diag(s_U) and w_U, so (D) should reduce to per-coordinate inequalities with moment sums, as in the proof
of Theorem C. (Z) needs the signs of all entries of G_z, i.e. the (C_u) block applied to the source of a single
coordinate z_h; this is the harder half. We did not attempt either in the time.

For m ≥ 2 the Z-pattern fails (entries of G_tt′ have coefficients of both signs), although Conjecture P holds in the
computed cases.

## 5. What else was tried, and what failed

- **The U2 pipeline for u = 3 and every m** (`closed_u3.py`). The closed forms of G_tt and G_tt′ need the configuration
  m = 3, U = {3, 4, 5}. The exact elimination over Q(s_0, …, s_5) did 10 of its 16 pivots in 3.5 minutes, with the cost
  per pivot growing (81 s for the tenth), and was stopped (`logs/abandoned/closed_u3.log`). Even with the closed forms, (A₃) is a
  polynomial inequality in five parameters (θ, φ, μ and two for the shape of w_U). Extrapolating the degrees from u = 2
  ((24, 16, 38, 28) after the triangle map) gives tens of millions of tensor Bernstein coefficients, beyond the 3 GB
  budget of the U2 code. (B₃) has four parameters and might be feasible, but alone it proves nothing for m ≥ 2.
- **(u, m) = (4, 1)** with the certificate of Section 3 (`closed_gen.py 4 1`, also with s_0 = 1). The exact elimination
  (17 pivots in five, resp. four, variables) slowed down sharply in the U-block (more than two minutes on the sixth
  pivot) and was stopped (`logs/abandoned/closed_u4_m1*.log`). This is a matter of computing time, not of method; nothing is
  claimed for (4, 1).
- **(u, m) = (3, 2)** (`closed_gen.py 3 2 dehom`, s_0 = 1, four variables). The elimination did 9 of its 13 pivots
  (140 s, the ninth alone 85 s) and was stopped during the tenth (`logs/abandoned/closed_u3_m2_dehom.log`). Even if it finished,
  the Laplace scheme for an 8 × 8 matrix needs 255 minors with much larger polynomials than in the (2, 2) case (which
  took about a minute for 63 minors); we estimate hours. Nothing is claimed for (3, 2).
- **Bareiss elimination and gcd-reduced LDLᵀ** for the minors of the 4 × 4 case (`pivots_m1.py`): correct but much
  slower than the Laplace scheme with early division by the denominator factors; replaced.
- **The merged comparison** (Section 4): not an identity and not one-sided.

## 6. Numerical sweep (`sweep_u.py`; floating point; no claim rests on it)

Samples: log-uniform weights with spreads up to e^{±5}, plus families with S or σ shrunk by up to e^{−10}, one U-weight
blown up by up to e^{8}, or equal U-weights; 60 samples per (u, m) (25 when u + m ≥ 8).

| quantity (worst over u ≤ 7, m ≤ 4, u + m ≤ 9) | value |
|---|---|
| sup ‖Ĝ_tt′‖_F σ/√(s_t s_t′) (A, Frobenius) | 0.3035 (u = 5, m = 2); u = 1: 0.283 |
| inf λ_min(Ĝ_tt)/(1 − s_t/σ), m ≥ 2 (B) | 0.7500 |
| weighted dominance ratio max_t Σ_{t′} ‖Ĝ_tt′‖₂ √(s_t′/s_t)/λ_min(Ĝ_tt) | ≤ 0.2927 |
| λ_min(Ĝ) | > 0 in every sample (smallest 7.4·10⁻⁸, u = 2, m = 1, S → 0) |
| ζ′-block eigenvalues; ζ′–(ρ, z) coupling | [0.719, 1.000]; ≤ 1.8·10⁻² |
| (ρ, z) block minus merged u = 1 block | entries up to 0.374; eigenvalues in [−0.374, 0.069] |

So the weighted block Gershgorin route of Propositions U1 and U2 (H = 2/5) is numerically sound for every u tested,
with the u = 1 constants; what is missing for each fixed u ≥ 3 is the closed forms and certificates in u + 2 parameters,
and for all u at once an argument that does not need closed forms.

## 7. Files

| File | Content |
|---|---|
| `run_all.sh` | reruns every step below, one Python process at a time, prints PASS/FAIL per script and ALL PASS; logs in `logs/` (about 7 minutes; `logs/abandoned/` keeps the logs of the stopped runs of Section 5) |
| `cf.py`, `lemmab.py`, `fastdirect.py`, `check_cf.py` | copied from `typeI_G_u2/`: closed form of the type I block, exact check against the Lemma B builder, float check against the direct construction |
| `flintfield.py`, `bern.py`, `test_bern.py`, `cu_form.py` | copied: rational function field on python-flint; Bernstein module (only a fallback in `prove_m1.py`, never reached) and its tests; the (C_u) form (not used) |
| `closed_m1.py 3` → `closed_m1_u3.pkl` | exact G_z for u = 3, m = 1 |
| `check_closed_m1.py 3` | closed form = direct exact elimination, 12 points × 2 gauges |
| `build_m1.py` → `build_m1.pkl` | leading minors, factor stripping, expansion checks |
| `prove_m1.py` | the sign certificate of Theorem U3₁ |
| `controls_m1.py` | point check and negative controls |
| `closed_gen.py u m [dehom]`, `check_closed_gen.py u m [n] [dehom]`, `minors_pos.py u m [dehom]` | the same pipeline for general (u, m) (Section 3) |
| `check_zpattern.py` | exact Z-matrix pattern (Z) for m = 1, u ≤ 3; prints the mixed pattern for m = 2 (Section 4b) |
| `check_mvec.py` | exact (D) for m = 1, u ≤ 3 with the explicit v (second proof of Theorem U3₁; Section 4b) |
| `numG.py`, `sweep_u.py`, `ex_mnum.py`, `ex_mvec.py`, `ex_beta.py`, `ex_beta2.py`, `ex1.py`, `ex_equal.py` | floating-point exploration (Sections 4, 4b, 6); `sweep_u.py`, `ex_mnum.py` and `ex_mvec.py` are rerun by `run_all.sh` |
| `closed_u3.py`, `pivots_m1.py`, `inspect_m1.py`, `ex_signs.py`, `ex_mmatrix.py` | abandoned or superseded exploration (Section 5); not used |
