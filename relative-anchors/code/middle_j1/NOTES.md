# (R∞)₁ for every number of copies: the blocks with |T| = 2 (and |T| = 3 in (R∞)₂), and reductions of the rest

Ansh Mishra and Aryan Senthilkumar, 5 November 2026. Nothing here is refereed.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Theorem D₁.** For every k ≥ 3 and every positive weight vector, every type I block of (R∞)₁ with \|T\| = 2 is positive definite (Section 2) | **proved**: closed form, moment bounds and reduction by hand; the last step is an exact coefficient-sign certificate for eight explicit rational functions of (α, β, e, λ) (`cert_t2.py`, under a second) |
| 2 | Checks of Theorem D₁: closed form against the Lemma B block (30 blocks, k = 3, …, 8), moment bounds and reduction at 300 rational configurations, positive definiteness at 60 rational weight vectors (k ≤ 9), negative control (B → (1 + 1/1000)B indefinite at an explicit point); the certified formulas equal the true f at k = 3 and bound it at k = 4, …, 8 | **exact check** (`check_t2.py`, `check_cert_formulas.py`) |
| 3 | The builder used here (Lemma B, `j1lib.py`) equals `ylib.Builder` on the letter-1 and type I blocks of (R∞)₁ for k = 4, 5, 6; block sizes 35, 14, 3 (k = 5) and 51, 17, 3 (k = 6); the kernels p ⊗ (columns) and positive definiteness after deleting one coordinate per kernel vector, at the weight vectors used | **exact check** (`check_ylib.py`) |
| 4 | \|T\| = 1: structure (Lemma T1.1), closed form of every entry and an inverse-free congruent form X on 2k coordinates (Lemma T1.3) | **proved** (by hand); **exact check** k = 3, …, 8 (`check_t1_closed.py`, `check_t1.py`) |
| 5 | \|T\| = 1: if the R–Q block is positive definite and the Schur complement S_P onto the states (1@r, 2@t) is a Laplacian with negative off-diagonal entries, the block is ⪰ 0 with kernel exactly p ⊗ \|2@t⟩ (Proposition T1.2) | **proved** |
| 6 | **Conjecture L** (the hypothesis of row 5, every k, all weights) | **proved for k = 3, 4** (exact symbolic certificates, `sym_t1.py`, `bareiss_t1.py 4`); k = 5 **not completed** (two symbolic runs stopped); **exact check** at 120 rational weight vectors, k = 3, …, 10 (`check_t1.py`); **open** in general |
| 7 | Letter-1 block: vacuum rows eliminated exactly (Lemma G), and an inverse-free congruent form (Lemma L1.1) | **proved** (by hand); **exact check** k = 4, …, 7 (`check_l1_closed.py`) |
| 8 | Letter-1 block: Conjecture L1 (sequential column-group Laplacian property) implies the block is ⪰ 0 with kernel exactly p ⊗ (letter-1 part of H₁) | **proved** (Section 4) |
| 9 | **Conjecture L1** | **exact check** at 32 rational weight vectors, k = 4, …, 7, two orders each (`check_l1.py`); **numerical**, 1 108 stages, k = 4, …, 7 (`explore_l1b.py`); **open** |
| 10 | The gap of the letter-1 block off its kernel is at least 1.74 min_h x_h (units 2^{k−1}) in all samples | **numerical** (`explore_l1.py`) |
| 11 | **Theorem D₂.** For every k ≥ 4 and every positive weight vector, every type I block of (R∞)₂ with \|T\| = 3 (size 4) is positive definite (Section 5) | **proved**: closed form, moment bounds and reduction by hand (joint concavity of the matrix-fractional function); the last step is an exact certificate of 8 + 64 vertex inequalities in (a₁, a₂, a₃, e, λ) (`cert_t3j2.py`, about 52 minutes); **exact checks** (`check_t3j2.py`, `check_cert_formulas.py`) |
| 12 | (R∞)₁ for every k ≥ 5 | **open**: the letter-1 blocks and the \|T\| = 1 blocks remain (reduced to Conjectures L1 and L); the type II blocks are Proposition T2′ of the README and the \|T\| = 2 blocks are Theorem D₁; k = 4 is Theorem A4 of the README |

Nothing here says anything about Conjecture 12 of `single-excitation-range/` in general, or about the NPT (Werner
distillability) problem. The new uniform-in-k statements are Theorem D₁ (the blocks of (R∞)₁ with |T| = 2) and
Theorem D₂ (the blocks of (R∞)₂ with |T| = 3), for every k and every positive weight vector; no other block of a
middle condition is proved here for k ≥ 5.

Notation as in `relative-anchors/README.md` (Sections 1 and 2: Lemmas 1, B, N, G). k copies, weights w_t > 0,
s_t = w_t², |s| = Σ_t s_t; we normalise |s| = 1 where convenient and write x_t = s_t/|s|. (R∞)₁ is Y compressed to
G1₁ ⊕ G0₂. Its torus blocks are the letter-1 block (size (3k² − k)/2), the type I blocks with |T| = 1 (size 3k − 1)
and |T| = 2 (size 3), and the type II blocks (positive definite for every k by Proposition T2′ of the README).

## 1. Files and how they relate

- `j1lib.py`: a direct implementation of Lemma B (entries of 2^Z Y, Z = k − j) for the blocks of (R∞)₁ (its `entry`
  is Lemma B for every j; Section 5 uses it with j = 2), exact (Fractions) or floating point; `bscale` multiplies the
  B-part (every term containing an E⁻¹) for negative controls.
- `ylib.py`, `flintfield.py`: copies of `relative-anchors/code/ylib.py` (the exact builder in the rescaled basis) and of
  the small rational-function field on top of python-flint.
- `lazyb.py`: a lazy copy of `ylib.Builder` (D⁻¹ only on the column sectors that occur), for the symbolic runs.
- `check_ylib.py`: row 3 of the table.
- `cert_t2.py`, `check_t2.py`: Theorem D₁ (Section 2). The file names refer to |T|: `t2` for |T| = 2 in (R∞)₁,
  `t3j2` for |T| = 3 in (R∞)₂, `t1` for |T| = 1, `l1` for the letter-1 block.
- `check_t1_closed.py`, `check_t1.py`, `sym_t1.py`, `bareiss_t1.py`: Section 3 (the blocks with |T| = 1);
  `bigstack.py` runs a script in a thread with a 128 MB stack (used for `bareiss_t1.py 5`, see Section 3).
- `check_l1_closed.py`, `check_l1.py`: Section 4 (the letter-1 block).
- `cert_t3j2.py`, `check_t3j2.py`: Section 5 (Theorem D₂, the blocks of (R∞)₂ with |T| = 3).
- `check_cert_formulas.py`: evaluates the very formulas certified by `cert_t2.py` and `cert_t3j2.py` against the
  blocks computed from Lemma B (Sections 2 and 5).
- `explore_*.py`: explorations quoted in the text, in floating point (and one SymPy expansion, `explore_t2d.py`); no
  claim rests on them; their outputs are in `logs/explore_*.log`.
- `guard.sh`: a copy of `relative-anchors/code/typeI_G_m1/guard.sh` (kills its child above 2.8 GB or after a time
  limit).
- `run_all.sh`: reruns every check, one Python process at a time under `guard.sh`, prints PASS or FAIL per step and
  ALL PASS at the end. `logs/run_all.out` and `logs/run_*.log` are from a complete rerun in a fresh copy of this
  folder (ALL PASS, about 55 minutes). Two earlier runs were stopped by hand during `cert_t3j2.py` and are kept as
  `logs/abandoned/run_all_aborted.out` and `logs/abandoned/run_all_aborted2.out`: the first to make the
  end-point condition of the Bernstein certificates explicit (Lemma D₁.5), the second to add
  `check_cert_formulas.py`; every step that had finished in them passed.

## 2. The type I blocks with |T| = 2: Theorem D₁

Let T = {t, t′}, U = [k] ∖ T (so |U| = k − 2), and normalise |s| = 1. Put α = x_t, β = x_{t′}, e = Σ_{u∈U} x_u =
1 − α − β, and

    F_i = Σ_{u∈U} x_u/(c_i − 2x_u),   c_0 = 2, c_1 = 2 − β, c_2 = 2 − α, c_3 = 2 − α − β.

The block has three states for every k: X₁ = (2@t, 2@t′), X₂ = (2@t′, 2@t) and X₃ = (vac, 2@t 2@t′).

**Lemma D₁.1 (closed form).** In units 2^{k−1} and in the basis X₁, X₂, X₃ the block is

    M = [[2, c, −A], [c, 2, −B], [−A, −B, 2/Γ₃]],

with g(x) = x/(2(1 − x)) and

- c = κ/(2Γ), κ = (1 + g(α))(1 + g(β)), Γ = 1 + ½(g(α) + g(β) + F_0);
- A = a′/Γ′, a′ = 1 + α/(2 − 2α − β), Γ′ = ½(1 + a′) + ½F_1;
- B = b′/Γ″, b′ = 1 + β/(2 − 2β − α), Γ″ = ½(1 + b′) + ½F_2;
- Γ₃ = 1 + ½F_3.

*Proof.* Write Ψ = y₁ e_a ⊗ |2@t′⟩ + y₂ e_b ⊗ |2@t⟩ + q vac ⊗ |2@t 2@t′⟩ with a = (2 at t), b = (2 at t′), and
evaluate the quadratic form after Lemma B term by term.
- Swap term: only (a, b) and (b, a) contribute, with σ_b⁻ψ_a = y₁|vac⟩ and σ_a⁻ψ_b = y₂|vac⟩, so it is y₁y₂.
- Diagonal term: the column of ψ_a has no letter at t and that of ψ_b none at t′, so it is 2y₁² + 2y₂².
- B-term: 𝒯_aψ_a = 𝒯_bψ_b = 0, while 𝒯_bψ_a = y₁v_{t′} and 𝒯_aψ_b = y₂v_t with v_t = ½w − w_t e_t in the letter-1
  sector of H₁, where E₁ = 2 diag(1 − x_u) + ½wwᵀ. So the term is −2y₁y₂⟨v_{t′}, E₁⁻¹v_t⟩.
- Vacuum coupling: σ_a⁻ maps |2@t 2@t′⟩ to |2@t′⟩, and 𝒯_a|2@t 2@t′⟩ = ½w′ − w_t e_t on the strings |1@u 2@t′⟩
  (u ≠ t′), where w′ is w without the entry t′ and E₂ = diag(2 − 2x_u − β) + ½w′w′ᵀ; J⁻ maps |1@u 2@t′⟩ to
  w_u|2@t′⟩. So the y₁q coefficient is −2(1 − ⟨w′, E₂⁻¹(½w′ − w_te_t)⟩), and symmetrically for y₂q.
- Vacuum diagonal: J⁺|2@t 2@t′⟩ = w_U on the strings |1@u 2@t 2@t′⟩ (u ∈ U), where E₃ = diag(2 − 2x_u − α − β)
  + ½w_Uw_Uᵀ, so the q² coefficient is 2(1 − ½⟨w_U, E₃⁻¹w_U⟩).

Each E is a diagonal matrix D plus ½wwᵀ, and Sherman–Morrison gives ⟨w, E⁻¹y⟩ = ⟨w, D⁻¹y⟩/(1 + ½⟨w, D⁻¹w⟩).
- With γ = Σ_u g(x_u) over all copies, Γ = 1 + γ/2, a = 1 + g(α) and b = 1 + g(β): ⟨v_{t′}, D⁻¹v_t⟩ = γ/4 −
  (g(α) + g(β))/2 and ⟨w, D⁻¹v_t⟩ = γ/2 − g(α), which give ⟨v_{t′}, E₁⁻¹v_t⟩ = ½(1 − ab/Γ). So
  M₁₂ = ½ − ½(1 − ab/Γ) = κ/(2Γ). Here γ = g(α) + g(β) + F_0, since g(x_u) = x_u/(c_0 − 2x_u).
- With γ′ = Σ_{u≠t′} x_u/(2 − 2x_u − β) = (a′ − 1) + F_1: 1 − (½γ′ − (a′ − 1))/(1 + ½γ′) = a′/(1 + ½γ′) = a′/Γ′.
  Symmetrically M₂₃ = −b′/Γ″.
- 2(1 − ½F_3/(1 + ½F_3)) = 2/Γ₃. ∎

`check_t2.py` (a) compares this closed form with the block computed from Lemma B, exactly, at rational weights for
k = 3, …, 8 (30 blocks), and `check_ylib.py` compares that builder with ylib.

**Lemma D₁.2 (Sylvester).** c < 2, and det M = (2/Γ₃) f with

    f = 4 − c² − Γ₃(A² + B² − cAB).

So M ≻ 0 if and only if f > 0.

*Proof.* Since α + β ≤ 1 we have αβ ≤ (1 − α)(1 − β), so g(α)g(β) ≤ ¼ and κ ≤ 5/4 + g(α) + g(β) < 4Γ. So c < 2, the
leading minors 2 and 4 − c² are positive, and expanding the determinant gives the formula. ∎

**Lemma D₁.3 (moment bounds).** Put Z = F_3/e. Then Z₀ ≤ Z ≤ Z₁ with Z₀ = 1/(2 − α − β) and Z₁ = 1/(α + β), and
for i = 0, 1, 2, with d_0 = α + β, d_1 = α, d_2 = β and h_i(z) = z/(1 + d_i z),

    e·ℓ_i(Z) ≤ F_i ≤ e·h_i(Z),

where ℓ_i is the chord of h_i between Z₀ and Z₁ (the affine function that agrees with h_i at Z₀ and Z₁).

*Proof.* For u ∈ U put z_u = 1/(c_3 − 2x_u). Since 0 < x_u ≤ e and c_3 − 2e = α + β, z_u ∈ (Z₀, Z₁]. Since
c_i − c_3 = d_i, 1/(c_i − 2x_u) = h_i(z_u). So F_3 = Σ_u x_u z_u and F_i = Σ_u x_u h_i(z_u), with Σ_u x_u = e. Each
h_i is concave and increasing on [Z₀, Z₁]. Jensen's inequality with the weights x_u/e gives F_i ≤ e·h_i(Z), and
h_i ≥ ℓ_i on [Z₀, Z₁] gives F_i ≥ Σ_u x_u ℓ_i(z_u) = e·ℓ_i(Z), because ℓ_i is affine. ∎

Both bounds are attained at the two extreme configurations: all of e on one copy (Z = Z₁) and e spread over many
tiny copies (Z → Z₀).

**Lemma D₁.4 (reduction to eight inequalities).** Fix the configuration and put λ = (Z − Z₀)/(Z₁ − Z₀) ∈ [0, 1]. For
each choice ε = (ε_0, ε_1, ε_2) ∈ {lo, hi}³ let g_ε be f evaluated with Γ₃ = 1 + eZ/2 and with F_i replaced by
e·ℓ_i(Z) (lo) or e·h_i(Z) (hi). Then f ≥ min_ε g_ε.

*Proof.* With Γ₃ fixed, f is a concave function of each of c, A, B separately (its second derivatives are −2, −2Γ₃,
−2Γ₃). c is a decreasing function of F_0, A of F_1 and B of F_2, so by Lemma D₁.3 the point (c, A, B) lies in a box
whose vertices are the eight choices ε. A function that is concave in each variable separately attains its minimum
over a box at a vertex. ∎

**Lemma D₁.5 (certificate, `cert_t2.py`).** For every ε, g_ε > 0 for all α, β, e > 0 (with α + β + e = 1) and all
λ ∈ [0, 1].

*Proof (exact, computer-assisted, under a second).* Replace every constant 1 by α + β + e. Then g_ε is a rational
function of (α, β, e, λ), homogeneous of degree 0 in (α, β, e), so the constraint α + β + e = 1 can be dropped. In
exact arithmetic (python-flint) g_ε = N_ε/D_ε, where
- D_ε factors into polynomials with only positive coefficients (degrees 1, 4, 4, 5);
- N_ε (894 to 1 333 terms, total degree 27 or 28, degree 6 or 7 in λ), written in the Bernstein basis of [0, 1] in λ,
  has Bernstein coefficients that are nonzero polynomials in (α, β, e) with only nonnegative coefficients (all of
  them nonzero; the script requires at least the first and the last to be nonzero, which is what positivity on the
  closed interval needs).

So N_ε and D_ε are positive for α, β, e > 0 and λ ∈ [0, 1]. ∎

A certificate only certifies the formulas it contains. `check_cert_formulas.py` executes the part of `cert_t2.py`
that builds the eight functions g_ε (so the same symbolic objects are used) and evaluates them exactly: at k = 3 the
configuration has λ = 1, all eight vertices coincide, and every g_ε equals the f of the block computed from Lemma B
by `j1lib.py` (30 rational weight vectors); for k = 4, …, 8 every g_ε is positive and the smallest is at most the true
f (30 more).

**Theorem D₁.** For every k ≥ 3 and every positive weight vector, every type I block of (R∞)₁ with |T| = 2 is positive
definite.

*Proof.* By Lemmas D₁.1 and D₁.2 it suffices to show f > 0. Since every weight is positive, e > 0. By Lemma D₁.4,
f ≥ min_ε g_ε, and by Lemma D₁.5 each g_ε is positive at the λ of the configuration. ∎

Remarks.
- The reduction loses almost nothing: numerically (`explore_t2h.py`) min_ε g_ε / min(1 − α, 1 − β) ≥ 2.43 over the
  triangle, the same constant as for f itself over random configurations (`explore_t2b.py`, `explore_t2c.py`). The
  block degenerates only at the two corners where one copy of T carries all the weight; there f → 0 linearly in
  1 − max(α, β), and the kernel direction is (0, 1, 1) (or (1, 0, 1)).
- Simpler reductions fail, which is why F_3 is kept exact: with independent bounds on F_0, F_1, F_2 and F_3 bounded
  through F_1 alone, the inequality fails when α, β → 0 (`explore_t2g.py`); the cone relaxation (F_0 ≤ F_1, F_2 ≤ F_3
  and a nonnegative mixed second difference d₃) fails, since in the cone variables (F_0, d₁, d₂, d₃) the
  coefficients of d₃, d₁³, d₂³ and others are negative for all (α, β) (`explore_t2d.py`, 58 monomials not of one
  positive sign); and
  Γ₃ ≤ (4 − c²)/(A² + B² − cAB) is not a concave constraint in (F_0, F_1, F_2) (`explore_t2e.py`), so convexity of
  the admissible set cannot be used directly.
- Negative control (`check_t2.py` (d)): with B replaced by (1 + 1/1000)B the block is indefinite at k = 3,
  w = (1/98, 1, 1/70), next to the corner β → 1; this scaled block agrees with ylib's A − (1001/1000)B there.
- The proof uses k only through e > 0: the copies of U enter only through the four sums F_i, and Lemma D₁.3 controls
  them for any number of copies. This replaces the moment reductions of Theorem B0: the four sums are values of one
  Stieltjes-type function of the measure Σ_u x_u δ_{x_u}, and keeping one of them exact confines the other three to
  explicit intervals.

## 3. The type I blocks with |T| = 1

Let T = {t} and σ = s_t. The block has 3k − 1 states:
- R_c = (2@t, 1@c) for c ∈ [k] (the row group of the row 2@t; c = t is allowed);
- P_r = (1@r, 2@t) for r ∈ [k] (the column group of the column 2@t);
- Q_c = (vac, 2@t 1@c) for c ≠ t.

Write Ψ = Σ_c ρ_c R_c + Σ_r π_r P_r + Σ_c q_c Q_c, and N_r = 2 for r ≠ t, N_t = 1.

**Lemma T1.1 (structure).** In units 2^{k−1}:
1. the block on the R-states is diagonal, with entries 2 (c ≠ t) and 1 (c = t);
2. the block on the P-states is Π = diag(N) − (N∘w)(N∘w)ᵀ/(2|s| − σ);
3. Σ_r w_r P_r = p ⊗ |2@t⟩ is in the kernel of the block.

In the rescaled basis of ylib (R_c and Q_c divided by w_c, P_r by w_r) the first two read diag(2s_c, s_t) and
Π̃ = diag(Ns) − (Ns)(Ns)ᵀ/(2|s| − σ), and the kernel vector is the all-ones vector on the P-states. Since
Σ_r N_rs_r = 2|s| − σ, Π̃ is the weighted Laplacian π̃ ↦ Σ_{r<r′} (N_rN_{r′}s_rs_{r′}/(2|s| − σ))(π̃_r − π̃_{r′})².

*Proof.* (1) The part of the form of Lemma B that is quadratic in ρ involves ψ_{(2,t)} = Σ_c ρ_c|1@c⟩ only. Its
swap term and its B-term vanish, because σ_{(2,t)}⁻ and both parts of 𝒯_{(2,t)} need a letter 2 at t in the column.
What is left is the diagonal term Σ_c (2 − [c = t])ρ_c². (2) Here ψ_{(1,r)} = π_r|2@t⟩. The swap terms vanish (no
letter 1 in the column), the diagonal term is Σ_r (2 − [r = t])π_r², and 𝒯_{(1,r)}|2@t⟩ = w_rN_r|2@t⟩, where E₁ on
the one-dimensional sector spanned by |2@t⟩ is 2(|s| − σ) + σ. (3) is the identity Σ_r w_r𝒯_{(r,1)} = E_c of
Lemma N. ∎

`check_t1.py` (a) checks (1)–(3) exactly in the rescaled basis at every weight vector it uses.

**Lemma T1.3 (closed form, and a reduction without inverses).** Normalise |s| = 1, so |w| = 1, and put
g_c = x_c/(2(1 − x_c)), γ = Σ_c g_c, τ = (1 + g_t)/(1 + γ/2), γ′ = Σ_{u∈U} x_u/(2 − σ − 2x_u), τ′ = 1/(1 + γ′/2).
In units 2^{k−1}, with the states ordered R, P, Q, the block is

    [[D_R, Y_RP, Y_RQ], [Y_RPᵀ, Π, Y_PQ], [Y_RQᵀ, Y_PQᵀ, Y_QQ]],

- D_R = diag(N), Π as in Lemma T1.1;
- Y_RP[c, r] = (τ_c/(2(1 − x_c)))(δ_rc − w_rw_c), with τ_c = τ for c ≠ t and τ_t = τ − 2, that is
  Y_RP = ½ diag(τ_c/(1 − x_c)) (I − wwᵀ), the analogue of the cross block −½τPΛ of Theorem B0;
- Y_PQ[r, c] = (τ′/(2 − σ − 2x_c))(N_rw_rw_c − (2 − σ)δ_rc) for c ∈ U;
- Y_QQ = 2(I + ½N′)⁻¹ (Lemma G), N′ = J⁻C⁻¹J⁺ on the strings |2@t 1@c⟩, with C = 2 − σ − 2x_c − 2x_u on the strings
  |2@t 1@c 1@u⟩: N′[c, c] = Σ_{u∈U∖c} x_u/(2 − σ − 2x_c − 2x_u), N′[c, c′] = w_cw_{c′}/(2 − σ − 2x_c − 2x_{c′});
- Y_RQ = −(I + ½N₁)⁻¹H, where N₁ = J⁻C₂⁻¹J⁺ on the letter-1 strings of H₁, with C₂ = 2(1 − x_c − x_u) on the pair
  strings |1@c 1@u⟩ (the matrix of Theorem B0: N₁[c, c] = Σ_{u≠c} x_u/(2(1 − x_c − x_u)),
  N₁[c, c′] = w_cw_{c′}/(2(1 − x_c − x_{c′}))), and
  H[c, c′] = δ_cc′(1 + σ/(2(1 − σ − x_c))) for c ∈ U, H[t, c′] = w_tw_{c′}/(2(1 − σ − x_{c′})).

The congruence ρ = (I + ½N₁)ρ′ and the elimination of q with Y_QQ⁻¹ = ½(I + ½N′) give the explicit form

    X = [[X₁₁, X₁₂], [X₁₂ᵀ, X₂₂]],   X₁₁ = (I + ½N₁)D_R(I + ½N₁) − ½H(I + ½N′)Hᵀ,
    X₁₂ = (I + ½N₁)Y_RP + ½H(I + ½N′)Y_PQᵀ,   X₂₂ = Π − ½Y_PQ(I + ½N′)Y_PQᵀ,

on (ρ′, π), with no inverse left. The block is ⪰ 0 with kernel exactly p ⊗ |2@t⟩ if and only if X ⪰ 0 with kernel
exactly span(0, w).

*Proof.* The entries follow from Lemma B as in Lemma D₁.1: E₁⁻¹ on the letter-1 strings and E₂⁻¹ on the strings
|2@t 1@u⟩ by Sherman–Morrison (the latter gives τ′ and the vector τ′w_u/(2 − σ − 2x_u)), and for the R–Q coupling the
identity 𝒯_{(2,t)}φ = ½J⁺q̃ − w_tΣ_c q_c|1@t 1@c⟩ (q̃ = Σ_c q_c|1@c⟩) together with J⁻E₂⁻¹J⁺ = 2(I − (I + ½N₁)⁻¹)
and J⁻E₂⁻¹ = (I + ½N₁)⁻¹J⁻C₂⁻¹ on the pair strings (Woodbury, as in Lemma G). The reduction is a congruence followed
by a Schur complement with a positive definite eliminated block. ∎

`check_t1_closed.py` checks every block of this closed form against the Lemma B block, the identity for X, and
X(0, w) = 0, exactly at 18 rational weight vectors with |w| = 1, k = 3, …, 8.

So the |T| = 1 problem is an explicit 2k × 2k form whose entries contain pair sums (through N₁ and N′) and, in
X₁₁ and X₁₂, sums over a third copy (products of N₁, N′ with diagonal matrices). We did not find an SOS identity for
it; Conjecture L below is the structure we observe.

**Proposition T1.2 (reduction).** Let C be the block on the R- and Q-states and S_P the Schur complement onto the
P-states, in the rescaled basis. If C ≻ 0 and S_P[r, r′] < 0 for all r ≠ r′, then the block is ⪰ 0 and its kernel is
exactly p ⊗ |2@t⟩.

*Proof.* The kernel vector (1, 0, 0) (ones on P, zero on R and Q) shows that the coupling of P to R ∪ Q annihilates
1, so S_P1 = 0. A symmetric matrix with zero row sums and negative off-diagonal entries is the Laplacian of a complete
graph with positive edge weights: vᵀS_Pv = Σ_{r<r′} |S_P[r, r′]|(v_r − v_{r′})², which is ⪰ 0 with kernel span(1).
With C ≻ 0 the block is ⪰ 0, and its kernel is the image of ker S_P, which is span(1, 0, 0). ∎

**Conjecture L.** For every k ≥ 3, every copy t and every positive weight vector: C ≻ 0, and S_P is a Laplacian with
negative off-diagonal entries.

Conjecture L would give every type I block of (R∞)₁ with |T| = 1, for every k, with the kernel p ⊗ |2@t⟩.

What we have:
- **Proved for k = 3 and k = 4** (`sym_t1.py 3`, `sym_t1.py 4`, exact, 10 s): over Q(s) with s_1 = 1 (the entries
  of the rescaled block are homogeneous of degree 1), all pivots of C are ratios of polynomials with only positive
  coefficients, S_P1 = 0 identically, and the numerator of every S_P[r, r′] (r ≠ r′) has only negative coefficients
  (about 4 900 terms, degree 30 or 31, at k = 4) over a denominator with only positive coefficients. For k = 4 this is
  a second, independent proof of the |T| = 1 rows of Theorem A4 of the README. `bareiss_t1.py 4` gives a third one,
  by the leading-minor certificate of `middle_k4/` (all ten leading-minor cores of the gauge-fixed block have only
  positive coefficients, 13 s).
- **k = 5: not completed.** Two exact symbolic attempts over Q(s₀, s₂, s₃, s₄) (four free parameters) were stopped:
  the rational-function elimination of `sym_t1.py 5` spent more than five minutes in its first vacuum pivot; the
  leading-minor certificate `bareiss_t1.py 5` (gauge-fixed size 13) found cores with only positive coefficients for
  D_2, …, D_7 (D_7: degree 27, 12 622 terms; not cross-checked, since the run never reached its final point checks)
  and was stopped while computing D_8, after three minutes at 1.3 GB. Two earlier runs ended with a segmentation fault
  inside the polynomial library at the same step; the run through `bigstack.py` did not. Neither is part of
  `run_all.sh`; the logs are `logs/sym_t1_k5.log` and `logs/bareiss_t1_k5.log`.
- **Exact checks** (`check_t1.py 10 15`): C ≻ 0 and the Laplacian property of S_P hold exactly at 120 rational weight
  vectors, 15 for each k = 3, …, 10, in six families (random, spread over powers of 3, one tiny weight, two tiny
  weights, one dominant weight, all weights equal but one), with t chosen at random. The smallest ratio
  |S_P[r, r′]|/S_P[r, r] seen is 1.0·10⁻⁸ (at extreme weights), so the property is tight in some limits and a proof
  must be exact there.
- **Not a consequence of the Laplacian Π̃.** S_P = Π̃ − K with K = BᵀC⁻¹B ⪰ 0 and K1 = 0. The off-diagonal entries of K
  take both signs, and S_P[r, r′]/Π̃[r, r′] ranges over [6·10⁻⁵, 130] in floating-point samples (`explore_t1d.py`).
  The entries of S_P depend on how the other weights are distributed, not only on their total (`explore_t1c.py`);
  structurally, the R–Q coupling goes through E₂⁻¹ on the pair sector, that is through (I + ½J⁻C₂⁻¹J⁺)⁻¹ on all k
  copies, and the Q-block through the pair kernel J⁻C⁻¹J⁺ on the strings |2@t 1@c⟩ (Lemma G), so pair functions of
  the weights enter and no reduction to a few sums over single copies is apparent.
- **Partial eliminations** (floating point, `explore_t1b.py`, 400 samples, k = 3, …, 9): eliminating only the
  Q-states (this is X₂₂ of Lemma T1.3) or only the R-states also leaves a matrix on P with nonpositive off-diagonal
  entries and zero row sums. So Conjecture L could be approached in two steps: X₂₂ is a Laplacian (an explicit
  statement about pair sums), and the correction from ρ′ does not destroy the sign pattern.
- **Negative control** (`check_t1.py` (c)): with B replaced by (1 + 1/1000)B the gauge-fixed block (one P-state
  deleted) is not positive definite at k = 4, t = 1, w = (1/39, 1, 1/26, 1/39), where the true one is.

## 4. The letter-1 block

The block has the states P_{r,c} = (1@r, 1@c) (r, c ∈ [k]) and Q_{cc′} = (vac, 1@{c, c′}), (3k² − k)/2 in all, and
the kernel vectors p ⊗ |1@c⟩: in the rescaled basis (P_{r,c} divided by w_rw_c, Q_{cc′} by w_cw_{c′}) these are the
indicators 1_c of the column groups {P_{r,c} : r ∈ [k]}. For k = 3 we have j = 1 = k − 2, and the block has the
extra kernel vector of Proposition K2; everything below is for k ≥ 4.

**The vacuum rows (exact, by Lemma G).** In units 2^{k−1} the vacuum block is 2(I + ½N₂)⁻¹ with N₂ = J⁻C₃⁻¹J⁺ on
the pair sector of H₂, C₃ = 2 diag(1 − x_c − x_{c′} − x_u) on the triples. So N₂ has diagonal entries
Σ_{u∉{c,c′}} x_u/(2(1 − x_c − x_{c′} − x_u)) and, between pairs that share one copy c, entries
w_{c′}w_{c″}/(2(1 − x_c − x_{c′} − x_{c″})). Eliminating the vacuum rows therefore costs no inverse:

    min_φ 2^{k−1}⟨Ψ, YΨ⟩ = F_PP(X) − ½⟨b(X), (I + ½N₂) b(X)⟩,   b(X) = Σ_r (σ_{(1,r)}⁺ − 𝒯_{(1,r)}ᵀE₂⁻¹J⁺) ψ_r,

where ψ_r = Σ_c X[r, c]|1@c⟩ and F_PP is the G1–G1 part of Lemma B (with E₁⁻¹ by Sherman–Morrison). The triple
functions 1/(|s| − s_c − s_{c′} − s_u) enter through N₂, and E₂⁻¹J⁺ = C₂⁻¹J⁺(I + ½N₁)⁻¹, with N₁ = J⁻C₂⁻¹J⁺ the k × k
matrix of Theorem B0, enters through b.

**Lemma L1.1 (an inverse-free form).** For every row r and pair P,

    (I + ½N₁) Y[(1@r, ·), (vac, P)] = −(I + ½N₁) e(r, P) + J⁻C₂⁻¹𝒯_{(1,r)}|P⟩,

where e(r, P) = e_{c*} if P = {r, c*} and 0 otherwise, and 𝒯_{(1,r)}|P⟩ = [r ∈ P] ½ Σ_{u≠c*} w_u|{c*, u}⟩ +
2w_r[r ∉ P]|P⟩. So the congruence X = X′(I + ½N₁) on the letter-1 rows (ψ_r = (I + ½N₁)ψ′_r for every r), followed by
the elimination of the vacuum states with the explicit inverse ½(I + ½N₂), turns the letter-1 block into an explicit
quadratic form in X′ whose only inverse is the Sherman–Morrison inverse of E₁ inside F_PP. The kernel becomes
{w c′ᵀ : c′ ∈ R^k} again, since (I + ½N₁) is invertible.

*Proof.* Lemma B gives the coupling as −⟨1@c|σ_{(1,r)}⁻ − J⁻E₂⁻¹𝒯_{(1,r)}|P⟩, and J⁻E₂⁻¹ = (I + ½N₁)⁻¹J⁻C₂⁻¹ on the
pair strings (Woodbury, E₂ = C₂ + ½J⁺J⁻). ∎

`check_l1_closed.py` checks the vacuum block 2(I + ½N₂)⁻¹ and Lemma L1.1 exactly against the Lemma B block
(rational weights with |w| = 1, k = 4, …, 7). The resulting form contains pair sums (N₁, N₂ and C₂⁻¹), sums over a
third copy (the products N₁ · diag · N₁ and the triple denominators of N₂), and the crossed terms X[r, c]X[c, r] of
the partial transpose. We did not find an SOS identity for it in the time available.

**A structural pattern (Conjecture L1).** Eliminate the vacuum states, then the column groups one at a time in any
order, each in the gauge in which one of its coordinates is dropped. Conjecture L1: at every stage, the diagonal block
of the current Schur complement on the next column group is a Laplacian with negative off-diagonal entries (its row
sums vanish automatically, because 1_c is in the kernel).

Conjecture L1 implies that the letter-1 block of (R∞)₁ is ⪰ 0 with kernel exactly p ⊗ (letter-1 part of H₁): the
pivot blocks of the gauge-fixed block in this order are the vacuum block (positive definite by Lemma G) and the
group blocks with one vertex deleted, which are positive definite for a connected Laplacian. This is the same
mechanism as Conjecture L, where there is a single column group.

Evidence (not proofs):
- **exact** (`check_l1.py 7 8`): every stage passes at 32 rational weight vectors (8 for each k = 4, …, 7; random,
  spread over powers of 2, one tiny weight, one dominant weight), each in the natural order and in a random order of
  the column groups (64 runs; `logs/check_l1_full.log`, 3 minutes). `run_all.sh` reruns the smaller set
  `check_l1.py 6 4` (24 runs);
- **numerical** (`explore_l1b.py`): 1 108 stages (k = 4, …, 7, random orders, 200 weight vectors including weights
  down to 10⁻⁷): the largest off-diagonal entry over the diagonal is −2.7·10⁻¹², the row sums vanish to 3·10⁻¹⁴;
- **numerical** (`explore_l1.py`): after eliminating everything except one column group, that group's block is a
  Laplacian; after eliminating only the vacuum, the entries between two states of the same column are negative, while
  the other entries have both signs; the smallest eigenvalue of the block off its kernel is at least 1.74 min_h x_h
  (units 2^{k−1}) in 150 samples (an observed lower bound only; we did not determine how the gap behaves as the
  smallest weight tends to 0, in contrast with Corollary S of the README for (R∞)_{k−2});
- **negative control** (`check_l1.py`): with B replaced by (1 + 1/1000)B the gauge-fixed block is not positive
  definite at k = 4, w = (1, 1, 1/30, 1/30).

## 5. The same method one level up: the blocks of (R∞)₂ with |T| = 3

The proof of Theorem D₁ uses only that the block has no letter 1, so that every entry is a Sherman–Morrison
expression in Stieltjes-type sums of the measure Σ_{u∉T} x_u δ_{x_u}. The same holds for the type I blocks of (R∞)_j
with |T| = j + 1. For j = 2 the method goes through.

Let T = {t₁, t₂, t₃}, a_i = x_{t_i}, x_T = a₁ + a₂ + a₃, U = [k] ∖ T, e = 1 − x_T. The block has four states for every
k: X_i = (2@t_i, 2 on T ∖ t_i), i = 1, 2, 3, and X₀ = (vac, 2 on T).

**Lemma D₂.1 (closed form).** In units 2^{k−2} the block is

    M = [[C, −a], [−aᵀ, 2/Γ_T]],   C = [[2, c₁₂, c₁₃], [c₁₂, 2, c₂₃], [c₁₃, c₂₃, 2]],   a = (A₁, A₂, A₃),

where, for {i, j, l} = {1, 2, 3}, with g_l(y) = y/(2 − 2y − a_l):
- c_ij = (1 + g_l(a_i))(1 + g_l(a_j)) / (2Γ^{(l)}), Γ^{(l)} = 1 + ½(g_l(a_i) + g_l(a_j) + F^{(l)}),
  F^{(l)} = Σ_{u∈U} x_u/(2 − 2x_u − a_l);
- A_i = (1 + g′_i)/(1 + ½(g′_i + F′_i)), g′_i = a_i/(2 − a_i − x_T), F′_i = Σ_{u∈U} x_u/(2 − 2x_u − (x_T − a_i));
- Γ_T = 1 + ½F_T, F_T = Σ_{u∈U} x_u/(2 − 2x_u − x_T).

*Proof.* As for Lemma D₁.1, with the letters 2 not involved in a given term acting as spectators: the B-term between
X_i and X_j lives on the strings |1@u 2@t_l⟩ (u ≠ t_l), where E = diag(2 − 2x_u − a_l) + ½wwᵀ, and gives
c_ij = ½ − ½(1 − (1 + g_l(a_i))(1 + g_l(a_j))/Γ^{(l)}); the vacuum coupling of X_i lives on the strings
|1@u 2 on T ∖ t_i⟩ (u ∉ T ∖ t_i), and the vacuum diagonal on |1@u 2 on T⟩ (u ∈ U). ∎

**Theorem D₂.** For every k ≥ 4 and every positive weight vector, every type I block of (R∞)₂ with |T| = 3 is
positive definite.

*Proof.*
1. **Sylvester.** M ≻ 0 if and only if C ≻ 0 and f := 2 − Γ_T aᵀC⁻¹a > 0.
2. **Moment bounds.** With z_u = 1/(2 − 2x_u − x_T) ∈ (Z₀, Z₁], Z₀ = 1/(2 − x_T), Z₁ = 1/x_T (since
   2 − 2e − x_T = x_T), every sum above is Σ_u x_u h_d(z_u) with h_d(z) = z/(1 + dz), where d = x_T − a_l for F^{(l)},
   d = a_i for F′_i, and d = 0 for F_T. So with Z = F_T/e and λ = (Z − Z₀)/(Z₁ − Z₀) ∈ [0, 1], Lemma D₁.3 gives
   e·ℓ_d(Z) ≤ F ≤ e·h_d(Z) for each of the six sums.
3. **Concavity.** For fixed Γ_T, the map (C, a) ↦ aᵀC⁻¹a is jointly convex on {C ≻ 0} × R³ (the matrix-fractional
   function), so f is jointly concave in (c₁₂, c₁₃, c₂₃, A₁, A₂, A₃) wherever C ≻ 0. Each c_ij and each A_i is a
   decreasing function of its sum, so given Z the point (c, a) lies in a box with 64 vertices (each sum at e·ℓ_d(Z)
   or e·h_d(Z)). If C ≻ 0 at the 8 vertices of the c-part of the box, then C ≻ 0 on the whole box (C is affine in c
   and the positive definite matrices form a convex set), f is concave on the box, and f ≥ the minimum of f over the
   64 vertices.
4. **Certificate** (`cert_t3j2.py`, exact, about 52 minutes). In homogeneous coordinates (a₁, a₂, a₃, e) and λ, at each
   of the 8 c-vertices the minors 4 − c₁₂² and det C, and at each of the 64 vertices P = 2 det C − Γ_T aᵀadj(C)a
   (so f = P/det C), are ratios N/D in which D is a product of polynomials with only positive coefficients and N, in
   the Bernstein basis of λ, has coefficients that are nonzero polynomials in (a₁, a₂, a₃, e) with only nonnegative
   coefficients (the largest N has 134 068 terms, degree 56, degree 13 in λ). So all of them are positive
   for a₁, a₂, a₃, e > 0 and λ ∈ [0, 1]. Since every weight is positive, e > 0 when k ≥ 4. ∎

Checks (`check_t3j2.py`): the closed form against the Lemma B block (j1lib's `entry`, which is Lemma B for every
j; exact, k = 4, …, 8, 15 blocks) and, for k ≤ 6, against ylib; the moment bounds, C ≻ 0 and f ≥ min of the vertex
values at 200 rational configurations; positive definiteness at 40 rational weight vectors (k ≤ 9); a negative
control (B → (1 + 1/1000)B indefinite at k = 4, w = (1, 1/160, 1/161, 1/162), next to the corner where one copy of T
carries all the weight). Numerically the reduction keeps the margin: min f/(1 − max_i a_i) is
about 1.07 over the vertices and 1.12 for the true blocks (`explore_t3j2.py`); f → 0 only when one copy of T carries
all the weight.

`check_cert_formulas.py` ties the certified formulas to the true blocks as for Theorem D₁: at k = 4 (λ = 1) all 64
vertex values 2 − Γ_T aᵀC⁻¹a, computed with exact linear algebra from the objects built in `cert_t3j2.py`, equal the
value of the block computed from Lemma B, and P/det C agrees with them (20 rational weight vectors); for k = 5, …, 8
the smallest vertex value is positive and at most the true value (20 more).

For k = 4 this is the row j = 2, |T| = 3 of Theorem B4 of the README (there by a certificate in s₁, s₂, s₃); Theorem D₂
is uniform in k. For j ≥ 3 the same scheme (the blocks of (R∞)_j with |T| = j + 1, of size j + 2) needs
2^{(j+1)(j+2)/2} vertex inequalities in j + 2 homogeneous variables and λ; we have not tried it.

## 6. What did not work, and what is left

- For |T| = 2: the cone relaxation, the box relaxation, a bound of F_3 through F_1 alone, and a convexity argument in
  the moment space all fail or do not apply (Section 2, remarks). Keeping F_3 exact and using concavity in the
  entries is what works; the same device proves Theorem D₂.
- For |T| = 1 and the letter-1 block we found no reduction to finitely many sums over single copies: the entries of
  S_P change when the other weights are redistributed with the same total (`explore_t1c.py`), and the closed forms
  contain pair functions of the weights. A proof of Conjecture L (or L1) has to
  control BᵀC⁻¹B, which contains the inverse of I + ½N₁ over all k copies, or work directly with the inverse-free
  forms of Lemmas T1.3 and L1.1, whose entries contain pair and triple sums. A naive M-matrix argument is ruled out:
  the correlations b_rᵀC⁻¹b_{r′} are negative for most pairs (`explore_t1d.py`), so S_P is not a Laplacian for the
  trivial reason C⁻¹ ≥ 0, B ≥ 0.
- No candidate SOS identity was written down or tested for the reduced forms of Lemmas T1.3 and L1.1; the time went
  into the Laplacian structure (Conjectures L and L1) instead.
- The direct symbolic certificates for |T| = 1 at k = 5 (four free parameters) did not finish (Section 3).

Open: Conjecture L (all k), Conjecture L1 (all k), and with them (R∞)₁ for every k ≥ 5.
At k = 4 everything is Theorem A4 of the README.

Possible next steps, in the order we would try them:
- a first step towards Conjecture L: X₂₂ = Π − ½Y_PQ(I + ½N′)Y_PQᵀ is, in the rescaled basis, a Laplacian (observed
  in `explore_t1b.py`); its off-diagonal entries are explicit sums over one copy of pair functions (Lemma T1.3), a
  statement of the size of (C0) in Theorem B0;
- the method of Theorems D₁ and D₂ for the blocks of (R∞)_j with |T| = j + 1 for j ≥ 3 (the vertex count grows as
  2^{(j+1)(j+2)/2}, so a symmetry reduction of the vertex set would be needed);
- for the letter-1 block, the first stage of Conjecture L1 (after the vacuum elimination of Lemma L1.1, the entries
  between two states of the same column are negative), which is again an explicit statement.
