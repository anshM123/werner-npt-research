# The condition G ≻ 0 of Theorem C: every number of letter-2 copies at equal U-weights, and what blocks the general case

Ansh Mishra and Aryan Senthilkumar, 5 November 2026. Nothing here is refereed.

Notation as in `relative-anchors/README.md` (Sections 1–3: Lemmas B, N; Theorem C; Proposition U1; Theorems U2, U3₁,
EU, TL) and the notes in `relative-anchors/code/typeI_u2/`, `typeI_G_u2/`, `typeI_G_m1/`, `typeI_G_u3/`. A type I
block of (R∞)_{k−1} is given by the set T ≠ ∅ of copies carrying one letter 2 (m = |T|) and U = [k] ∖ T (u = |U|);
s_x = w_x², σ = s_T, S = s_U, K = σ + S/2. Forms are in units of 2Y. G is the Schur complement on the letter-2 states
(ρ_t, ζ_{t,h}) after eliminating the P-states and Q modulo the kernel span{p ⊗ f_h : h ∈ U} (Theorem C).

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Theorem EM.** For every u ≥ 3, every m ≥ 2, every positive weight vector on T and every common weight s > 0 of the copies of U: G ≻ 0 | **proved**, computer-assisted: closed forms symbolic in u (Sections 2, 3), weighted block Gershgorin with H = 1/2 on the invariant and the standard part separately (Section 4), five exact coefficient-sign certificates (Section 5, about 2 s) |
| 2 | **Corollary EM.** For every k, every type I block of (R∞)_{k−1} in which the copies of U all have the same weight (arbitrary weights on T, any u ≥ 1, m ≥ 1) is ⪰ 0 with kernel exactly span{p ⊗ f_h : h ∈ U} | **proved** (Theorem EM for u ≥ 3, m ≥ 2; Theorem EU for m = 1; Theorem U2 for u = 2; Proposition U1 for u = 1; then Theorem C) |
| 3 | Closed forms of G_tt and G_tt′ at equal U-weights: invariant part (2 × 2 in (ρ_t, z_t)) and standard part (1 × 1), as rational functions of (u, σ, s, s_t, s_t′) | **exact** symbolic computation (`equal_m.py`, 57 s); **exact** agreement with direct elimination of the whole block at 44 rational configurations, 268 blocks (u = 3…7, m = 1…4), and gauge independence at 10 of them (`check_equal_m.py`) |
| 4 | The five Gershgorin conditions (A_inv), (A_std), (B_inv) (two inequalities), (B_std) with H = 1/2 | **proved**: every numerator coefficient positive, every denominator factor with positive coefficients, in orthant variables (Section 5); the same holds for H = 2/5, 29/100 and 7/10 |
| 5 | End-to-end: at 21 rational configurations (u = 3…7, m = 2…5) the Lemma B builder (`fastj1.py`) equals the closed form `cf.py`; Y(p ⊗ f_h) = 0; all pivots of the gauge-reduced block are positive; G ≻ 0; the five Gershgorin quantities computed from that G are positive and **equal the certified rational functions** at the same point | **exact check** (`check_e2e.py`) |
| 6 | Negative controls: the certificate fails for H = 1/4 (A_inv) and H = 4/5 (B_inv, B_std); explicit rational configurations where these inequalities are false (from the exact G of the whole block); with B → (1 + 1/1000)B the gauge-reduced block has negative pivots at three explicit configurations where the true block has all pivots positive | **exact** (`controls_em.py`) |
| 7 | Step 1 for general U-weights (closed forms through a few symmetric quantities): G_tt is not a function of (s_t, σ, Σs_h, Σs_h²): exact witness w_U = (1, 9, 10) vs (5, 6, 11) | **exact** (`witness_moments.py`); consistent with the pair-structure obstruction of `typeI_G_m1/NOTES.md`, Section 4 |
| 8 | Where a uniform bound fails (general U-weights): the Frobenius-norm constant of (A) grows with u (largest found 0.293 at u = 3, 0.309 at u = 6, 0.329 at u = 10, 0.354 at u = 16, 0.388 at u = 24), while the spectral-norm constant stays at 0.2777–0.2778 (≈ 5/18, the u = 1 value) and inf B = 3/4 | **numerical** (`ex_Agrowth.py`, `sweep_general.py`; Section 8) |
| 9 | G ≻ 0 for general U-weights with u ≥ 3 (all of (R∞)_{k−1} for k ≥ 5) | **open** |

So the type I blocks of (R∞)_{k−1} are now settled, for every k, when |U| ≤ 2 (any weights), when |T| = 1 and the
U-weights take at most two values (Theorems EU, TL), when k = 4 (Theorem U3₁), and, new here, when all copies of U
have the same weight, for every |T|. For k ≥ 5 they remain open for general weights on U with |U| ≥ 3. Nothing here
says anything about Conjecture 12 or the NPT problem beyond these blocks.

## 1. The formula for G

Work in the gauge Σ_{t∈T} w_t ψ_t = 0 (ψ_t = row t of X), a complement of the kernel (`typeI_G_u2/NOTES.md`,
Section 1). The eliminated variables are X_U (all u² entries), q and the rows ψ_t; the block on (X_U, q) is the (C_u)
block C of Theorem C (positive definite), and the ψ_t enter only through Σ_t (|ψ_t|² + ψ_t·ℓ_t) with

    ℓ_t = −ζ_t − 2 𝒯_{(t,2)}ᵀ E_{T∖t}⁻¹ (𝒯_{(t,1)} − w_t) y_t,     y_t = ρ_t e_t + ζ_t,

whose minimum under the gauge is −¼Σ_t|ℓ_t|² + |Σ_t w_tℓ_t|²/(4σ). Hence

    G_tt  = diag(2, I) − ¼(1 − s_t/σ) L_tᵀL_t − C_tᵀ C⁻¹ C_t,
    G_tt′ = ½ e_ρe_ρᵀ − P_tt′ + (w_t w_t′/(4σ)) L_tᵀL_t′ − C_tᵀ C⁻¹ C_t′        (t ≠ t′),

with L_t: y_t ↦ ℓ_t, C_t the coupling of y_t to (X_U, q) (the terms ρ_t tr X, −2q(2Kρ_t + w_t w_U·ζ_t)/(2K − s_t) and
−2Σ_{r∈U}⟨𝒯_{(t,2)}ψ_r, E_{T∖t}⁻¹𝒯_{(r,1)}y_t⟩), and P_tt′ = 𝒯_{(t′,2)}ᵀE_{T∖tt′}⁻¹𝒯_{(t,2)} the pair term. Here
E_{T∖t} = diag(d) + ½w̃w̃ᵀ on U ∪ {t} with d_t = K + s_t/2, d_h = s_h + K − s_t/2, and E_{T∖tt′} = diag(s_x + K −
(s_t + s_t′)/2) + ½w̃w̃ᵀ on U ∪ {t, t′}. G_tt depends only on (s_t, σ, w_U), G_tt′ only on (s_t, s_t′, σ, w_U).

## 2. Equal U-weights: invariant and standard parts

Let s_h = s = ω² for all h ∈ U, so S = us. The whole block, the gauge and C are invariant under the permutations of
U, and C ≻ 0, so the minimiser over the eliminated variables is unique and equivariant. The letter-2 data split into
- the **invariant part**: ρ_t and ζ_t = w_t ω z_t 1_U (two coordinates (ρ_t, z_t) per copy t); its minimiser has
  X_U = s(a I + b(J − I)) and some q;
- the **standard part** (u ≥ 3): ρ_t = 0, ζ_t = w_t ω z_t f with Σ_h f_h = 0; its minimiser has X_U = s(α diag f +
  β f1ᵀ + γ 1fᵀ) and q = 0 (the standard representation occurs exactly three times in R^{U×U} when u ≥ 3, and not in
  q); the value is |f|² times a form in the z_t.

The two parts are G-orthogonal (Schur's lemma), and G is the direct sum of the 2m × 2m matrix G^inv on the invariant
coordinates and of G^std ⊗ I_{u−1}, with G^std an m × m matrix. So **G ≻ 0 ⟺ G^inv ≻ 0 and G^std ≻ 0.** The blocks
(t, t′) of G^inv and G^std are those of Section 1 restricted to the two parts.

## 3. The closed forms (`equal_m.py`, `check_equal_m.py`)

`equal_m.py` evaluates every term of Section 1 on the two families, symbolically in u, with T = {t₁, t₂} carrying the
weights s₁, s₂ and further T-copies of zero data (total weight σ ≥ s₁ + s₂). It reuses the validated machinery of
`typeI_G_m1/equal_u.py` for the (C_u) form and for the sums over r ∈ U (split into r = r′ and r ≠ r′, with
Σ_r f_r = 0, Σ_r f_r² = |f|², Σ_{r≠r′} f_r f_r′ = −|f|², valid for every u). Vectors on U ∪ {t, t′} stay in
span{e_t, e_t′, 1_U} (invariant) or span{f} (standard), and E⁻¹ is applied by Sherman–Morrison. The square roots
w_t, ω are carried as auxiliary generators; every exponent of them is even at the end (checked: an odd exponent
raises an error), so all entries are rational in (u, σ, s, s₁, s₂). The elimination of (a, b, q), resp. (α, β, γ), is
exact over the rational function field (python-flint). Result: G^inv as a 4 × 4 matrix on (ρ₁, z₁, ρ₂, z₂) and G^std
as a 2 × 2 matrix on (z₁, z₂); the diagonal blocks give G_tt (they do not involve the other weight, checked) and the
off-diagonal blocks give G_tt′. Numerators have up to 2 576 terms and degree ≤ 14 in each variable; the common
denominators are products of the factors of E_{T∖t}, E_{T∖tt′} (2K − s_t, the Sherman–Morrison factors 2 + g) and one
factor of degree 5 in (u, σ, s) from C (the factor of Theorem EU).

Validation (`check_equal_m.py`): at 44 rational configurations (u = 3, …, 7; m = 1, …, 4; T-amplitudes and the common
U-amplitude drawn from {1/20, …, 10}), the whole block is built by `cf.TypeI` over Q, the coordinates P_{r0,h} are
deleted, everything but the letter-2 states is eliminated exactly, and the projections of G_z (z-coordinates
ζ_{t,h} = w_t w_h z_{t,h}) on the two parts are compared with the closed forms: 268 blocks (t, t′), all equal. The
exact G_z is invariant under the transposition of two U-copies, and at 10 configurations a second gauge (P_{r1,h}
deleted) gives the same G_z.

## 4. The frame and the weighted block Gershgorin criterion

Frame (as in Proposition U1 and Theorem U2): in the invariant part, ζ_t = c ẑ_t w_U/|w_U| with c = √(1 + S/σ), i.e.
ẑ_t = κ_t z_t with κ_t² = s_t S σ/(σ + S); in the standard part ẑ_t = w_t ω |f| z_t (orthonormal). Let Ĝ^inv, Ĝ^std be
G^inv, G^std in these coordinates (congruences, so positivity is unchanged), θ_t = s_t/σ, and H a constant. Consider

- **(A_inv)** ‖Ĝ^inv_tt′‖_F < H √(s_t s_t′)/σ and **(A_std)** |Ĝ^std_tt′| < H √(s_t s_t′)/σ, for t ≠ t′;
- **(B_inv)** Ĝ^inv_tt − H(1 − θ_t) I ≻ 0 and **(B_std)** Ĝ^std_tt > H(1 − θ_t), for every t.

**Lemma W.** If m ≥ 2 and these hold for all t ≠ t′, then G ≻ 0.

*Proof.* For Ĝ^inv and a vector with 2-blocks x_t: xᵀĜx ≥ Σ_t λ_min(Ĝ_tt)|x_t|² − Σ_{t≠t′}‖Ĝ_tt′‖₂|x_t||x_t′|, and
2|x_t||x_t′| ≤ √(s_t′/s_t)|x_t|² + √(s_t/s_t′)|x_t′|² gives xᵀĜx ≥ Σ_t |x_t|²(λ_min(Ĝ_tt) − Σ_{t′≠t}‖Ĝ_tt′‖₂√(s_t′/s_t)).
By (A_inv) the sum is < (H/σ)Σ_{t′≠t}s_t′ = H(1 − θ_t), and by (B_inv) this is ≤ λ_min(Ĝ_tt). So Ĝ^inv ≻ 0. The same
argument with scalars gives Ĝ^std ≻ 0 from (A_std), (B_std). Section 2 gives G ≻ 0. ∎

Applying the criterion to the two parts separately matters: the full off-diagonal block Ĝ_tt′ also contains u − 1
copies of the standard entry, so its Frobenius norm grows like √u; the spectral norm does not (Section 8).

## 5. The certificates (`prove_em.py`)

Put s_t = a, s_t′ = b, σ = a + b + r (r ≥ 0 the weight of the other T-copies) for (A), and s_t = a, σ = a + r (r > 0,
since m ≥ 2) for (B); the U-weight is s and u = 3 + x. Each condition, multiplied out without square roots, is a
rational function of (a, b, r, s, x):

    A_inv = H² ab/σ² − [g₀₀² + g₀₁²/κ_t′² + g₁₀²/κ_t² + g₁₁²/(κ_t²κ_t′²)],    A_std = H² ab/σ² − g²/(ab s²),
    B_inv,00 = g₀₀ − H r/σ,    B_inv,det = (g₀₀ − H r/σ)(g₁₁ − H (r/σ) κ_t²) − g₀₁g₁₀,    B_std = g/(a s) − H r/σ,

with g the corresponding entries of G^inv_tt′, G^std_tt′, G^inv_tt, G^std_tt in z-coordinates (B_inv,det is κ_t² det of
the 2 × 2 matrix of (B_inv), and B_inv,00 > 0 with B_inv,det > 0 is (B_inv)). For H = 1/2, after reduction:

| condition | numerator terms | degrees in (a, b, r, s, x) | negative coefficients |
|---|---|---|---|
| A_inv | 58 905 | (33, 33, 32, 32, 32) | 0 |
| A_std | 14 950 | (23, 23, 22, 22, 22) | 0 |
| B_inv,00 | 664 | (13, 0, 14, 13, 13) | 0 |
| B_inv,det | 1 687 | (19, 0, 19, 18, 18) | 0 |
| B_std | 352 | (11, 0, 11, 10, 10) | 0 |

Every irreducible factor of every reduced denominator has only positive coefficients. Each numerator has a monomial
free of x and r (for (A)), resp. free of x (for (B)), and so does every denominator factor. Hence all five are positive
for a, b, s > 0, r ≥ 0 (resp. r > 0) and x ≥ 0, that is, for every configuration with m ≥ 2 and every real u ≥ 3.
The same certificate goes through for H = 2/5, 29/100 and 7/10 (`run_all.sh`), so the window of admissible constants
is wide: numerically (`ex_frame.py`) sup ‖Ĝ^inv_tt′‖_F σ/√(s_ts_t′) = 0.2832 (θ = φ = 1/2, S/σ → 0, all u),
sup |Ĝ^std_tt′| σ/√(s_ts_t′) = 0.0555, inf λ_min(Ĝ^inv_tt)/(1 − θ) = 0.7508 and inf Ĝ^std_tt/(1 − θ) = 0.7507 (θ → 0).

## 6. Theorem EM and the corollary

**Theorem EM.** Let u ≥ 3 and m ≥ 2, let all copies of U carry the same weight s > 0, and let the weights of T be
arbitrary positive numbers. Then G ≻ 0.

*Proof.* By Section 2, G ≻ 0 ⟺ G^inv ≻ 0 and G^std ≻ 0. The blocks of G^inv and G^std are the closed forms of
Section 3 (an exact symbolic computation, valid for every real u ≥ 3; every pivot of the eliminated block is a
leading principal minor ratio of the positive definite (C_u) block restricted to an invariant subspace, so it is
positive and the closed forms are the Schur complements at every positive weight vector). By Section 5 the five
conditions of Section 4 hold with H = 1/2 for every pair t ≠ t′ and every t. By Lemma W, G ≻ 0. ∎

**Corollary EM.** For every k and every type I block of (R∞)_{k−1} whose copies without a letter 2 all have the same
weight (arbitrary weights on T), the block is ⪰ 0 with kernel exactly span{p ⊗ f_h : h ∈ U}.

*Proof.* G ≻ 0 by Theorem EM (u ≥ 3, m ≥ 2), Theorem EU (m = 1, every u), Theorem U2 (u = 2) and Proposition U1
(u = 1). Theorem C and the reduction of `typeI_u2/NOTES.md`, Section A3, give the block. ∎

This settles every type I block whose set U has equal weights, whatever the weights on T. It does not give
(R∞)_{k−1} for general weights: U ranges over all subsets, and a weight vector for which every possible U with
|U| ≥ 3 has equal weights is constant when k ≥ 4, a case already covered by Section 4 of the README.

## 7. Checks and controls

- `check_e2e.py` (exact): at 21 configurations (u = 3…7, m = 2…5) the whole block from the Lemma B builder
  `fastj1.FastJ1` (sector-by-sector E⁻¹, not the closed form) equals `cf.TypeI` entry by entry; Y(p ⊗ f_h) = 0; with
  P_{r0,h} deleted all LDLᵀ pivots are positive; G ≻ 0; the five Gershgorin quantities computed from this G are
  positive and equal the values of the certified rational functions of `prove_em.py` at the same (a, b, r, s, x).
- `controls_em.py` (exact):
  - the certificate fails for H = 1/4 (A_inv has 139 negative coefficients) and for H = 4/5 (B_inv,det, B_std);
  - the inequalities are false there: (A_inv) with H = 1/4 at w_T = (1, 1), ω = 1/10, u = 3 (value −0.0035);
    (B_std) with H = 4/5 at w_T = (1/10, 1), ω = 1, u = 3 (−0.042); (B_inv,det) with H = 4/5 at w_T = (1/10, 1),
    ω = 1/10, u = 3 (−5.6·10⁻⁶); all computed from the exact G of the whole block; with H = 1/2 all values there are
    positive;
  - B → (1 + 1/1000)B: at w_T = (1/100, 1/100), ω = 1, u = 3; w_T = (1, 1/100), ω = 1/100, u = 3; and
    w_T = (1/30, 1/20, 1/10), ω = 1, u = 4, the gauge-reduced block has 4, 4 and 1 negative pivots, while the true
    block has all pivots positive (smallest 9.4·10⁻⁵, 3.9·10⁻⁴, 6.4·10⁻³). The scaled builder differs from the
    closed form, so the comparison of `check_e2e.py` is not vacuous.

- `emcheck.py` (exact, written separately; `run_emcheck.sh`): at weight vectors equal on a random set U₀ of u ≥ 3
  copies and arbitrary on the other m = k − u ≥ 2 copies, every torus group of (R∞)_{k−1} is built from Lemma B
  (`j1lib.J1` for k ≤ 6, `fastj1.FastJ1` for k = 7, 8), its known kernel vectors (`crosscheck.known_kernel`) are
  annihilated exactly, the group without a set of pivot coordinates of them is positive definite by the exact
  test `crosscheck.pd_cert`, and the kernel dimensions add up to |H_{k−1}| + k − 1. Ten points at k = 5, eight at
  k = 6, four at k = 7 and two at k = 8 (u = 3, …, 6, m = 2, …, 5): all pass (`logs/emcheck_k*.log`). This uses
  none of the closed forms of Sections 1–3.

## 8. General U-weights: the closed-form question and where a uniform bound fails

**Closed forms through a few symmetric quantities (Step 1).** For general U-weights this is not possible in the
naive sense. Exact witness (`witness_moments.py`): the U-amplitudes (1, 9, 10) and (5, 6, 11) give squared weights
with equal first and second power sums (182 and 16 562; both are a, b, a + b with a² + ab + b² = 91), and the ρ_t ρ_t
entry of G differs, for m = 1 (w_T = (2)) and m = 2 (w_T = (2, 3)). So G_tt is not a function of (s_t, σ, S, Σs_h²).
At u = 3 the U-weights are determined by three power sums, so this is as far as a finite witness at u = 3 goes; the
structural reason for every u is the pair structure of the crossed term of the (C_u) form (`typeI_G_m1/NOTES.md`,
Section 4). At equal U-weights everything collapses to rational functions of (u, σ, s, s_t, s_t′), which is what
Sections 2–5 use; with finitely many weight levels the same holds with the multiplicities as parameters (as in
Theorem TL), but the number of unknowns grows with the square of the number of levels.

**The constants for general U-weights (Step 3, numerical, no claim rests on it).** In the frame of Proposition U1
(numG.Ghat: ρ_t, c z_t, ζ′_t ⊥ w_U), with A_F = ‖Ĝ_tt′‖_F σ/√(s_ts_t′), A₂ the same with the spectral norm, and
B = λ_min(Ĝ_tt)/(1 − θ):
- `sweep_general.py` (u = 3, 4, 5, 6, 8, 10; m = 2, 3; seven families: log-uniform, S shrunk, σ shrunk, one U-weight
  dominant, one tiny, two values, equal): sup A_F = 0.3165 (u = 10, two values), inf B = 0.7500, weighted dominance
  ratio ≤ 0.2916, λ_min(Ĝ) > 0 in every sample;
- `ex_Agrowth.py` (m = 2, random search over two-level, spread and one-dominant U-weights): the largest A_F found is
  0.293, 0.309, 0.329, 0.354, 0.388 for u = 3, 6, 10, 16, 24, always with S ≪ σ, while the largest A₂ is
  0.2777–0.2778 for every u (≈ 5/18, the u = 1 value) and the smallest B is ≥ 0.79 in these samples (0.75 in the sweep).

So the **Frobenius-norm version of (A_u) is not uniform in u**: the u − 1 directions ⊥ w_U carry many small
off-diagonal entries (at equal U-weights these are the u − 1 copies of the standard entry, of size ≤ 0.0555 each),
and their Frobenius norm grows with u. Whether it eventually exceeds the diagonal constant 3/4 we do not know, but a
uniform proof cannot rely on it. The spectral-norm version, A₂ ≤ 5/18 + o(1) together with B ≥ 3/4, is numerically
uniform. At equal U-weights the decomposition of Section 2 is exactly what replaces the Frobenius norm of the full
block by the norms of the two parts; for general weights the ζ′-part is not decoupled (coupling up to 1.8·10⁻² in
`typeI_G_u3/NOTES.md`, Section 6), and a uniform argument would need (i) a bound of the spectral norm of Ĝ_tt′ that
does not pass through the Frobenius norm of the ζ′-part, and (ii) control of C_tᵀC⁻¹C_t′ exactly on the soft modes of
C (`typeI_G_u3/NOTES.md`, Section 4). Neither was attempted beyond these observations.

## 9. Files

| File | Content |
|---|---|
| `run_all.sh`, `guard.sh` | reruns every step, one Python process at a time under `guard.sh` (killed above 2 500 000 KB or after its time limit); PASS or FAIL per step, ALL PASS at the end; `logs/run_all.out`, one log `logs/run_<step>.log` per step (about 3 minutes) |
| `equal_m.py` → `equal_m.pkl` | Sections 2–3: the closed forms of G^inv, G^std at equal U-weights, symbolic in u (forms in SymPy, elimination in python-flint) |
| `check_equal_m.py` | Section 3: exact validation against direct elimination (44 configurations, 268 blocks; gauge independence) |
| `prove_em.py` | Section 5: the five coefficient-sign certificates (`python prove_em.py [Hnum Hden]`, default H = 1/2) |
| `check_e2e.py` | Section 7: end-to-end exact checks with the Lemma B builder, and equality with the certified functions |
| `controls_em.py` | Section 7: negative controls (H = 1/4, 4/5, explicit witnesses, B → (1 + 1/1000)B) |
| `witness_moments.py` | Section 8: G_tt is not a function of (s_t, σ, S, Σs_h²) |
| `sweep_general.py`, `ex_Agrowth.py`, `ex_frame.py` | floating point / exploration: the constants for general and for equal U-weights (Section 8, Section 5); no claim rests on them |
| `emlib.py` | shared helpers: loading the closed forms, exact G of a whole block from `cf.py` or from the Lemma B builder, the Gershgorin quantities computed directly |
| `ex_factors.py`, `ex_witness.py` | exploration: factor structure of the closed forms; locating the control witnesses |
| `cf.py`, `equal_u.py`, `flintfield.py`, `numG.py`, `fastdirect.py` | copied from `relative-anchors/code/typeI_G_m1/` (cf.py, equal_u.py, flintfield.py), `typeI_G_u3/` (numG.py) and `typeI_G_u2/` (fastdirect.py, for the state list) |
| `fastj1.py`, `j1lib.py` | the fast exact Lemma B builder (FastJ1) and the reference builder it wraps (copied) |
| `emcheck.py`, `crosscheck.py`, `run_emcheck.sh` | the second check of Section 7 (whole torus groups from Lemma B; `crosscheck.py` supplies the states, the known kernel vectors and the exact test pd_cert) |
| `logs/` | output of a complete run of `run_all.sh` |
