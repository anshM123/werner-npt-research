# The condition G ≻ 0 of Theorem C at two U-weight levels: the marked copy for every k, and two-valued weights

Ansh Mishra and Aryan Senthilkumar, 6 November 2026. Nothing here is refereed.

Notation as in `relative-anchors/README.md` (Sections 1–3: Lemmas B, N; Theorem C; Proposition U1; Theorems U2, U3₁,
EU, TL, EM) and `relative-anchors/code/typeI_G_em/NOTES.md` (Section 1: the formula for G; Section 4: the frame and
Lemma W). A type I block of (R∞)_{k−1} is given by the set T ≠ ∅ of copies carrying one letter 2 (m = |T|) and
U = [k] ∖ T (u = |U|); s_x = w_x², σ = s_T, S = s_U, K = σ + S/2. Forms are in units of 2Y. G is the Schur complement
on the letter-2 states (ρ_t, ζ_{t,h}) after eliminating the P-states and Q modulo the kernel span{p ⊗ f_h : h ∈ U}
(Theorem C). Here U = U₁ ⊔ U₂, with u_i = |U_i| copies of squared weight s_i and S_i = u_i s_i.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Lemma S (two levels).** G = G^inv ⊕ (G^{std,1} ⊗ I_{u₁−1}) ⊕ (G^{std,2} ⊗ I_{u₂−1}). G^inv has 3 × 3 blocks (t, t′) in (ρ_t, z_t¹, z_t²) and G^{std,i} is m × m. The minimisers lie in explicit families with 7 unknowns (invariant part) and 5 unknowns (standard part) | **proved** by hand (Section 2); equivariance **checked exactly** (Section 3) |
| 2 | Closed forms of the blocks (t, t′), symbolic in the level sizes: (a) u₁ = 1, u₂ ≥ 2 and **arbitrary** T-weights (s_t, s_t′, σ free); (b) T-weights in {s₁, s₂}: u₁ = 1; u₁, u₂ ≥ 2 both symbolic; u₁ = 2, …, 16 fixed with u₂ symbolic; (c) all standard parts | **exact** symbolic computation (`tl_m.py`). (a)–(c) agree **exactly** with direct elimination of whole blocks: 685 blocks (t, t′) at 152 rational configurations for (a), for (b) with u₁ = 1 or both sizes symbolic, and for (c); the fixed-u₁ forms at u₁ = 2, 3, 5 directly (24 configurations) and at every u₁ ≤ 16 against the symbolic forms (Section 3) |
| 3 | **Theorem M1.** U = {h₀} ⊔ U₂ with u₂ ≥ 2, m ≥ 2, all copies of U₂ and of T of one squared weight and h₀ of any squared weight: G ≻ 0 | **proved**, computer-assisted: exact S_m-splitting into D − O and D + (m − 1)O, Sylvester, 10 coefficient-sign certificates (Section 4, under 1 s) |
| 4 | **Corollary M1.** For every k and every weight vector whose squared weights are all equal but one, (R∞)_{k−1} holds completely: every torus block is ⪰ 0 with exactly its known kernel | **proved** (row 3, or row 5, with the theorems listed in Section 6). For this condition it extends Theorem I1 of the marked-copy notes (k ≤ 29) to every k |
| 5 | **Theorem TV.** u = u₁ + u₂ ≥ 3 with u₁, u₂ ≥ 1, m ≥ 2, all squared weights of U and of T in {s₁, s₂}: G ≻ 0, for **all** level sizes | **proved**, computer-assisted: weighted block Gershgorin (Lemma W) with H = 1/2 in a two-level frame; 129 exact coefficient-sign certificates (19 for u₁ = 1, 20 with both sizes symbolic, 90 for the fixed sizes u₁ = 2, …, 16; Section 5), and the third pivot of (B_inv) with both sizes symbolic by determinant and exact interpolation in u₁ (Proposition B3, Section 10). Before Section 10 it was proved for min(u₁, u₂) ≤ 16 |
| 6 | **Corollary TV.** For **every k**, (R∞)_{k−1} holds completely at every weight vector whose squared weights take at most two values | **proved** (row 5 with the theorems listed in Section 6). Before Section 10 it was proved for k ≤ 35, and for every k when one of the two values is carried by at most 16 copies |
| 7 | **Proposition B3.** The third pivot of (B_inv) with both sizes symbolic, for all real u₁, u₂ ≥ 2, s_t = s₁ and s_t = s₂. This was the one condition of Lemma W left open for min(u₁, u₂) ≥ 17 (the off-diagonal conditions (A), all standard parts and the first two pivots of (B_inv) were already certified for all sizes) | **proved**, computer-assisted, exact (Section 10): det A = Ñ/D_c for the 3 × 3 matrix A of (B_inv). Ñ is determined by exact interpolation in u₁ from 77 exact slices u₁ = 2, …, 78 for each type (11 910 terms each) under a rigorous degree bound (B = 76, deg_{u₁} R = 40); every coefficient of Ñ is ≥ 0 (139 005 terms for each type), and D_c > 0. The direct elimination had needed more than 2.5 GB (Section 8) |
| 8 | End-to-end: at 18 two-valued configurations (k ≤ 8) the Lemma B builder equals the closed form `cf.py`, Y(p ⊗ f_h) = 0, all pivots of the gauge-reduced block are positive, G ≻ 0, and all Gershgorin quantities computed from this G are positive; for u₁ = 1 they **equal the certified rational functions** at the same point | **exact check** (`check_e2e2.py`) |
| 9 | Negative controls: the certificate fails for H = 1/10 and H = 4/5, and explicit exact witnesses show these inequalities false; B → (1 + 1/1000)B gives negative pivots at four configurations where the true block is positive definite; Theorem M1 with I/1000 subtracted fails, with an explicit indefinite point | **exact** (`controls2.py`, Section 7) |
| 10 | G ≻ 0 at two U-levels with **arbitrary** T-weights (the statement aimed at) | **open**: closed forms for u₁ = 1 done and validated (row 2a); the Gershgorin conditions hold numerically with room (sup A = 0.2885, inf B = 0.7503); the off-diagonal certificate is out of reach in size (Section 8) |
| 11 | Floating-point exploration of the Gershgorin constants (arbitrary and two-valued T-weights); factor and sign structure of the u₁ = 1 forms | **numerical / exploration** (`ex_numframe.py`, `ex_twoval.py`, `ex_factor_one.py`, `ex_signs_one.py`); no claim rests on them |
| 12 | Checks and controls of Proposition B3: Ñ/D_c equals det A computed directly from the closed forms at 5 points far outside the grid, for both types; for a larger H there are explicit rational points with a negative third pivot, and the slice through such a point fails the sign test; one slice changed by the factor 1 + 1/1000 is detected by the vanishing test of the interpolation; every exponent of D_c is needed | **exact** (`em3points.py`, `em3neg.py`; Section 10.4) |

So for every k, (R∞)_{k−1} now holds completely at the weights (σ, 1, …, 1) and at every weight vector whose squared
weights take at most two values (rows 5–7). For general weights the type I blocks with |U| ≥ 3 and k ≥ 5 remain open,
and so does the two-level case with arbitrary weights on T. Nothing here says anything about Conjecture 12 or the NPT
problem beyond these blocks.

## 1. The formula for G

As in `typeI_G_em/NOTES.md`, Section 1 (valid for every weight vector): in the gauge Σ_{t∈T} w_t ψ_t = 0,

    G_tt  = diag(2, I) − ¼(1 − s_t/σ) L_tᵀL_t − C_tᵀ C⁻¹ C_t,
    G_tt′ = ½ e_ρe_ρᵀ − P_tt′ + (w_t w_t′/(4σ)) L_tᵀL_t′ − C_tᵀ C⁻¹ C_t′        (t ≠ t′).

Here C is the (C_u) block of Theorem C on (X_U, q), which is positive definite, and C_t is the coupling of
y_t = ρ_t e_t + ζ_t to (X_U, q). L_t is the map y_t ↦ ℓ_t = −ζ_t − 2𝒯_{(t,2)}ᵀE_{T∖t}⁻¹(𝒯_{(t,1)} − w_t)y_t, and P_tt′
is the pair term ⟨𝒯_{(t′,2)}y_t, E_{T∖tt′}⁻¹𝒯_{(t,2)}y_t′⟩. G_tt depends only on (s_t, σ, w_U) and G_tt′ only on
(s_t, s_t′, σ, w_U). So the blocks can be computed with two copies A, B of T carrying data and further copies of zero
data of total weight σ − s_A − s_B. Since 2^Z Y does not change when all weights are multiplied by one positive
factor, we put s₂ = 1.

## 2. Two levels: the splitting (Lemma S)

The type I block, the gauge and C are invariant under the group S_{U₁} × S_{U₂} of the permutations of U that keep
the weights. Since C ≻ 0, the minimiser over the eliminated variables is unique for given letter-2 data and depends
equivariantly on it; hence G commutes with S_{U₁} × S_{U₂}. Per copy t the letter-2 data space is R ⊕ R^{U₁} ⊕ R^{U₂}.
It consists of three copies of the trivial representation (ρ_t and the level sums z_t¹, z_t²) and of the standard
representation std_i of S_{U_i} for each level with u_i ≥ 2. By Schur's lemma

    G = G^inv ⊕ (G^{std,1} ⊗ I_{u₁−1}) ⊕ (G^{std,2} ⊗ I_{u₂−1}),

so **G ≻ 0 ⟺ G^inv ≻ 0 and G^{std,i} ≻ 0 for each level with u_i ≥ 2.**

Coordinates are ζ_{t,h} = w_t w_h z_{t,h}. Invariant data have z_{t,h} = z_t^i for h ∈ U_i, and
G^inv[z^i, z^j] = Σ_{h∈U_i, h′∈U_j} G_z[z_{t,h}, z_{t′,h′}]. Standard data of level i have z_{t,h} = z_t f_h with f
supported on U_i and Σf = 0, and G^{std,i}_tt′ = G_z[z_{t,h}, z_{t′,h}] − G_z[z_{t,h}, z_{t′,h′}] (h ≠ h′ in U_i).

The minimisers, with X_rh = w_r w_h x_rh:
- **Invariant data.** x_hh = p_i (h ∈ U_i), x_rh = o_ij (r ∈ U_i, h ∈ U_j, r ≠ h; o₁₁ is absent when u₁ = 1) and q:
  seven unknowns, six when u₁ = 1. These are the invariant elements of R^{U×U} ⊕ R_q, one per orbit.
- **Standard data of level 1** (f on U₁). The minimiser is x = α diag f + β f1ᵀ + γ 1fᵀ on U₁ × U₁, x_rh = δ f_r on
  U₁ × U₂, x_hr = ε f_r on U₂ × U₁, zero on U₂ × U₂, and q = 0. When u₁ = 2 the three matrices are dependent and γ is
  dropped. The std₁-isotypic part of R^{U×U} has multiplicity 3 + 1 + 1 (2 + 1 + 1 if u₁ = 2), and q is trivial.
  Level 2 is the same with the levels exchanged. ∎

## 3. The closed forms (`tl_m.py`) and their validation (`check_tl.py`)

`tl_m.py` evaluates every term of Section 1 on these families, symbolically in the level sizes. It generalises the
two-level reduction of `typeI_G_m1/two_level_fast.py` (SymPy sparse polynomial rings over a rational function field)
to two letter-2 copies A, B, including the gauge terms ℓ_t and the pair term.
- Vectors on U ∪ {A, B} are combinations of e_A, e_B, 1_{U₁}, 1_{U₂}, e_r, e_r′ and f.
- Each sum over r, r′ ∈ U is split by the levels of r and r′ and by r = r′ or r ≠ r′ (u_i and u_i(u_j − δ_ij) terms).
- For the standard part the identities Σ_r f_r = 0, Σ_r f_r² = |f|² and Σ_{r≠r′} f_r f_r′ = −|f|² inside the level
  of f are used; they hold for every size.
- E⁻¹ is applied by Sherman–Morrison.
- The square roots w_A, w_B, ω₁ are carried as generators. Every exponent of them is even at the end (checked: an odd
  exponent raises an error).
- The elimination of the 7 (resp. 5) unknowns is exact over the rational function field (python-flint).
Output: the blocks (A, A) and (A, B) of G^inv (3 × 3) and of G^{std,i} (1 × 1).

Forms computed (s₂ = 1; variables in brackets):
- `tl_inv_one.pkl`: u₁ = 1, u₂ ≥ 2, arbitrary T-weights [u₂, σ, s₁, a = s_A, b = s_B]. The entries of G_tt have up
  to 12 007 terms and those of G_tt′ up to 90 907 terms, over a common denominator of 101 315 terms. That denominator
  factors into σ, the linear factors 2K − s_t and 2K − s_t + 2 and a cubic in s_t from E_{T∖t} (for t = A, B), a
  quartic from E_{T∖AB}, and one factor of degree 10 in (u₂, σ, s₁) (`ex_factor_one.py`).
- Modes T11, T12 and T22, where (s_A, s_B) = (s₁, s₁), (s₁, 1), (1, 1) and σ is free:
  - `tl_inv_one_T*`: u₁ = 1 [u₂, σ, s₁], at most 2 478 terms;
  - `tl_inv_gen_T*`: u₁, u₂ ≥ 2 [u₁, u₂, σ, s₁], at most 34 768 terms, 2–3 minutes each;
  - `tl_inv_fN_T11`, `tl_inv_fN_T22`: u₁ = N = 2, …, 16 fixed [u₂, σ, s₁];
  - standard parts: `tl_std2_one_T*` (level 2; u₁ = 1, u₂ ≥ 3), `tl_std2_onetwo_T*` (u₁ = 1, u₂ = 2),
    `tl_std1_gen_T*` (level 1; u₁ ≥ 3, u₂ ≥ 1) and `tl_std1_two_T*` (u₁ = 2, u₂ ≥ 1).

**Validation** (`check_tl.py`, exact). Each form is checked at 8 rational configurations: level sizes in the range of
the form, k ≤ 10, T-amplitudes and ω₁ from {1/20, …, 10}, the T-weights required by the mode, and for T12 extra
T-copies of other weights. At each configuration:
- the whole block is built by `cf.TypeI` over Q, the coordinates P_{r0,h} are deleted, and everything but the
  letter-2 states is eliminated exactly;
- the projections of G_z on the parts are compared with the closed forms: 685 blocks (t, t′) at 152 configurations
  (19 forms), all equal;
- the exact G_z is invariant under a transposition inside each level, and the diagonal closed forms do not involve
  s_B;
- at a third of the points, a second gauge (P_{r1,h} deleted) gives the same G_z.
The fixed-u₁ forms are the same computation with u₁ a number. They agree exactly with direct elimination for
u₁ = 2, 3, 5 (`check_tl.py`, 4 configurations for each of the six forms). For every N = 2, …, 16 they also agree exactly with the symbolic forms at
u₁ = N, at four random rational points per form (`check_fN.py`).

## 4. Theorem M1: one marked copy in U, all other weights equal

**Theorem M1.** Let U = {h₀} ⊔ U₂ with u₂ ≥ 2 and m ≥ 2. Let every copy of U₂ and of T have squared weight 1, and h₀
squared weight s₁ > 0. Then G ≻ 0.

*Proof (computer-assisted, exact; `prove_eqT.py`).* Here σ = m, and the block also commutes with the permutations of
T. Let D = G_tt and O = G_tt′ (t ≠ t′; the same for all pairs; O is symmetric, checked) on the invariant part, and d,
o the same on the standard part of U₂. Write x_t = x̄ + x̃_t with Σ_t x̃_t = 0. Then the form of G^inv is
m x̄ᵀ(D + (m − 1)O)x̄ + Σ_t x̃_tᵀ(D − O)x̃_t, and likewise for the standard part. So G ≻ 0 if and only if
D − O ≻ 0, D + (σ − 1)O ≻ 0, d − o > 0 and d + (σ − 1)o > 0.

We prove these for every real σ ≥ 2, u₂ ≥ 2 and s₁ > 0. Put u₂ = 2 + x and σ = 2 + y; for the standard part put
u₂ = 3 + x, with a separate form for u₂ = 2. The three LDLᵀ pivots of D − O and of D + (σ − 1)O (order ρ, z¹, z²) and
the four scalars are rational functions of (x, y, s₁):
- every numerator has only positive coefficients (120, 285, 340; 1 771, 4 959, 5 394; 20, 1 330; 10, 136 terms) and a
  monomial free of x and y;
- every denominator factors into polynomials with positive coefficients, each with a monomial free of x and y.
So all are positive for x, y ≥ 0 and s₁ > 0, and Sylvester's criterion gives the four conditions. ∎

## 5. Theorem TV: two-valued weight vectors

**Theorem TV.** Let u₁, u₂ ≥ 1 with u₁ + u₂ ≥ 3, and let m ≥ 2. Suppose the squared weights of all copies of U and
of T lie in {s₁, s₂}: U₁ has weight s₁, U₂ weight s₂, and the two values are distributed in any way on T. Then
G ≻ 0. (Until Section 10 this was proved under the extra condition min(u₁, u₂) ≤ 16.)

**The frame.** For each copy t:
- the coordinate ρ_t is kept;
- in the invariant part, y¹ = S₁z¹ + S₂z² and y² = z¹ − z², with ẑ¹ = κ₁y¹, ẑ² = κ₂y², κ₁² = s_tσ/(S(σ + S)) and
  κ₂² = s_tS₁S₂/S. Then ζ_t = c ẑ¹ w_U/|w_U| + ẑ² n_⊥, where c² = 1 + S/σ and n_⊥ is the unit invariant vector
  orthogonal to w_U. This is the frame of Proposition U1 and Theorem EM, with the second invariant direction left
  unscaled;
- the standard parts use ẑ = w_t ω_i |f| z (orthonormal).
Write θ_t = s_t/σ, and let H be the Gershgorin constant.

**Lemma W** (`typeI_G_em/NOTES.md`, Section 4, applied to each part separately). Suppose that for all t ≠ t′
‖Ĝ^inv_tt′‖_F < H√(s_ts_t′)/σ and |Ĝ^{std,i}_tt′| < H√(s_ts_t′)/σ, and that for all t
Ĝ^inv_tt − H(1 − θ_t)I ≻ 0 and Ĝ^{std,i}_tt > H(1 − θ_t). Then G ≻ 0. (This is the weighted block Gershgorin
criterion with weights √s_t: Σ_{t′≠t}‖Ĝ_tt′‖√(s_t′/s_t) < (H/σ)Σ_{t′≠t}s_t′ = H(1 − θ_t).)

*Proof of Theorem TV (computer-assisted, exact; `prove_gersh.py`).* By homogeneity s₂ = 1 and s₁ = s > 0.
- **The pairs.** Every pair (s_t, s_t′) is (s, s), (s, 1), (1, s) or (1, 1). The norm for (1, s) is that for (s, 1),
  since Ĝ_t′t = Ĝ_tt′ᵀ.
- **The orthant variables.** For (A) put σ = s_t + s_t′ + r, where r ≥ 0 is the weight of the other copies of T; for
  (B) put σ = s_t + r with r > 0 (since m ≥ 2). The level sizes are u_i = u_i⁰ + x_i.
- **The conditions**, as rational functions of the closed forms of Section 3:

      A_inv = H²s_ts_t′/σ² − Σ_{i,j} G_y,tt′[i,j]²/(κ_{t,i}²κ_{t′,j}²)   (κ_{t,0} = 1, G_y = NᵀG_zN, z = Ny),
      B_inv: the three LDLᵀ pivots of G_y,tt − H(r/σ) diag(1, κ_{t,1}², κ_{t,2}²),
      A_std = H²s_ts_t′/σ² − (G^std_tt′)²/(s_ts_t′s_i²),   B_std = G^std_tt/(s_ts_i) − H r/σ.

  They are computed exactly in the variables of the closed forms, and only the final rational function is composed
  into the orthant variables (x₁, x₂, r, s).
- **The certificate.** The numerator has only positive coefficients and a monomial free of the shift variables and r.
  Every irreducible factor of the denominator has coefficients of one sign and a monomial free of them, and the total
  sign is +. So each condition holds for all x_i ≥ 0, r ≥ 0 and s > 0.

The certified cases (H = 1/2 throughout):
- **u₁ = 1, u₂ ≥ 2.** Invariant part with u₂ = 2 + x₂; standard part of U₂ with u₂ = 3 + x₂ and with u₂ = 2. All 19
  conditions are certified (largest numerator 23 422 terms; `logs/run_prove_gersh_one.log`). The same holds with
  H = 2/5.
- **u₁, u₂ ≥ 2, both symbolic** (invariant part with u_i = 2 + x_i; standard part of U₁ with u₁ = 3 + x₁ or u₁ = 2,
  and u₂ = 1 + x₂). Certified: the off-diagonal conditions (A) for all three pairs; all ten standard-part conditions;
  the first two pivots of (B_inv) for s_t = s and s_t = 1. The numerators have up to 635 135 terms (the z-block of
  the pair (s, 1), degrees 58, 58, 60, 61 in x₁, x₂, r, s); the (B_inv) pivots have 27 027 and 78 539. To stay below the
  memory limit, each Frobenius condition is certified in two parts: the entries in the row or column of ρ, and the
  2 × 2 z-block. Each part is held to the budget H²/2, so together they give the condition with H. Each group runs
  in its own process. The third pivot of (B_inv), for s_t = s and s_t = 1, is Proposition B3 (Section 10). It is
  proved by determinant and interpolation, since its direct elimination exceeds the memory limit (Section 8).
- **u₁ = N fixed, N = 2, …, 16, u₂ = 2 + x₂ symbolic.** All three pivots of (B_inv) for s_t = s and s_t = 1 are
  certified, 90 conditions. The third pivot has 11 910 terms for every N.
- **u₁ ≥ 2, u₂ = 1.** Exchange the levels and divide all weights by s₁. The weights stay two-valued in the same two
  values, and the conditions are invariant (degree 0 in the weights, frame-independent norms). This is the case
  u₁ = 1 with s → 1/s.
- **The standard part of U₂ when u₁, u₂ ≥ 2.** The same exchange turns it into the standard part of level 1 of the
  configuration (u₂, u₁) with s → 1/s. This is covered, since the level-1 certificates hold for every size u₂ ≥ 1 of
  the other level.
- **(B_inv) with 2 ≤ u₂ ≤ 16 and any u₁ ≥ 2** follows in the same way from the fixed-N certificates. Both types,
  s_t = s and s_t = 1, are certified, so the exchange of the two values on T is covered. Proposition B3 covers all
  u₁, u₂ ≥ 2 directly, so for min(u₁, u₂) ≤ 16 the fixed-N certificates are now a second, independent proof of the
  third pivot.
- **u₁ = 0 or u₂ = 0** (one level) is Theorem EM.
Hence every condition of Lemma W holds for all level sizes, and G ≻ 0. ∎

## 6. The corollaries

**Corollary TV.** Let the squared weights of a weight vector on k copies take at most two values. Then, for every
k, every torus block of (R∞)_{k−1} is ⪰ 0 with exactly its known kernel.

*Proof.* The letter-1 block is Theorem T and the type II blocks are Proposition T2. For the type I blocks:
- U = ∅ is Proposition U0, |U| = 1 is Proposition U1 and |U| = 2 is Theorem U2;
- for |U| ≥ 3, Theorem C reduces the block to G ≻ 0;
- if |T| = 1, G ≻ 0 holds by Theorems EU and TL;
- if |T| ≥ 2 and U carries one value, by Theorem EM;
- if |T| ≥ 2 and U carries both values, T carries values from the same two, and Theorem TV applies; with
  Proposition B3 (Section 10) it holds for all level sizes. ∎

Before Section 10 the last case was covered only for min(u₁, u₂) ≤ 16. That gave the corollary for k ≤ 35 (since
min(u₁, u₂) ≥ 17 and |T| ≥ 2 force k ≥ 36), and for every k when one of the two values is carried by at most 16
copies.

**Corollary M1.** For every k and every weight vector whose squared weights are all equal but one, (R∞)_{k−1} holds
completely. This is a special case of Corollary TV. For the type I blocks with the odd copy in U, |U| ≥ 3 and |T| ≥ 2,
Theorem M1 gives a second, independent proof of G ≻ 0 (exact splitting, no Gershgorin step). For this condition it
extends Theorem I1 of the marked-copy notes (`../marked_copy/typeI/NOTES.md`, k ≤ 29) to every k. The other conditions (R∞)_j
with j ≤ k − 2 are not touched here.

## 7. Checks and controls

- `check_tl.py` (exact): Section 3.
- `check_fN.py` (exact): (1) the fixed-u₁ forms equal the symbolic forms at u₁ = N, entry by entry, at 4 random
  rational points per form, N = 2, …, 16 (1 800 entry values); (2) the three (B_inv) pivots at H = 1/2, evaluated from the symbolic forms,
  are positive at 288 rational points with (u₁, u₂) ∈ {(17, 17), (17, 40), (25, 19), (50, 50), (100, 30), (200, 200)},
  s ∈ [10⁻³, 10³], r ∈ [10⁻³, 50], both types of s_t. The smallest value of (third pivot)/(H(1 − θ)κ₂²) found is
  0.5. This was evidence for row 7 before Section 10; it is consistent with Proposition B3.
- Proposition B3 (Section 10.4, exact): `em3points.py` (Ñ/D_c against det A computed directly far outside the grid)
  and `em3neg.py` (negative controls for the larger H, the vanishing test and the candidate denominator).
- `check_e2e2.py` (exact) checks 18 two-valued configurations: u₁ = 1, 2, 3, u₂ = 2, …, 5, (m₁, m₂) ∈ {(0, 2), (1, 1),
  (2, 0), (1, 2), (2, 1), (0, 3)}, k ≤ 8. At each one:
  - the whole block from the Lemma B builder `fastj1.FastJ1` (sector-by-sector E⁻¹, not the closed form) equals
    `cf.TypeI` entry by entry;
  - Y(p ⊗ f_h) = 0, all LDLᵀ pivots of the block with P_{r0,h} deleted are positive, and G ≻ 0;
  - all Gershgorin quantities (A_inv for every ordered pair, the three B_inv pivots, both standard parts) computed
    from this G are positive at H = 1/2;
  - for u₁ = 1 the A_inv values and B_inv pivots equal the certified rational functions of `prove_gersh.py` at the
    same point (x₂, r, s).
- `controls2.py` (exact):
  - with H = 1/10 three (A) conditions and with H = 4/5 eight (B) conditions are not certified;
  - the inequalities are false there, computed from the exact G of the whole block: (A_inv) with H = 1/10 at
    w_T = (20, 20), w_U = (20, 1, 1) (value −0.0024), and (B) with H = 4/5 at w_T = (1/10, 1), w_U = (1/10, 1, 1)
    (−0.042); with H = 1/2 all values there are positive;
  - B → (1 + 1/1000)B: at w_T = (1/100, 1/100), w_U = (1/100, 1, 1); w_T = (1, 1), w_U = (1/100, 1, 1);
    w_T = (1, 1, 1), w_U = (1/30, 1, 1, 1); and w_T = (1/20, 1), w_U = (1/20, 1/20, 1, 1), the gauge-reduced block has
    6, 3, 4 and 3 negative pivots. The true block has all pivots positive there (smallest 2.0·10⁻⁸, 4.0·10⁻⁵,
    4.8·10⁻⁴, 8.3·10⁻⁴). The scaled builder differs from the closed form, so the comparison of `check_e2e2.py` is not
    vacuous;
  - Theorem M1 with I/1000 subtracted from D − O and D + (σ − 1)O: the certificate fails, and D − O − I/1000 is
    indefinite at u₂ = 2, σ = 2, s₁ = 10⁻⁴, where D − O ≻ 0.
- Floating point (no claim rests on it), two-valued family (`ex_twoval.py`; u₁ ≤ 3, u₂ ≤ 8, m ≤ 3, s from e⁻⁶ to e⁶):
  sup A_inv = 0.140, sup A_std = 0.031, inf B = 0.751. So H = 1/2 has room on both sides, and the split budget H²/2
  (that is, 0.354 per part) is far from tight.

## 8. Where it stops

**The third pivot of (B_inv) with both sizes symbolic (row 7; now proved in Section 10).** For u₁, u₂ ≥ 2 symbolic,
the first two pivots of
(B_inv) are certified (27 027 and 78 539 terms in (x₁, x₂, r, s)). For the third pivot the exact elimination exceeded
the 2.5 GB limit in three variants (`logs/abandoned/`):
- in the (ρ, y¹, y²) coordinates with late composition (peak above 2.52 GB, for s_t = s and for s_t = 1);
- in the z-coordinates (a congruent form);
- with the entries composed into orthant variables first (peak 3.1 GB).
At fixed u₁ the third pivot has only 11 910 terms, for every u₁ = 2, …, 16. So the obstacle is the size of the
intermediate gcds, not of the result. Two routes looked feasible:
- reconstruct the pivot by interpolation in u₁ from fixed-u₁ results, with a rigorous degree bound;
- a fraction-free elimination that tracks the known denominator factors.
The first one, applied to det A instead of the pivot, is carried out in Section 10 and completes Theorem TV and
Corollary TV for every k. The second was not needed.

**Arbitrary weights on T (row 10).** The original target was two U-levels with arbitrary T-weights. For u₁ = 1 the closed
forms exist (`tl_inv_one.pkl`, validated), and the Gershgorin route is numerically sound. Over 540 samples
(`ex_numframe.py`: u₁ = 1, …, 4, u₂ = 2, …, 8, m = 2, 3, log-uniform and extreme families):
- sup ‖Ĝ^inv_tt′‖_F σ/√(s_ts_t′) = 0.2885 (at S ≪ σ);
- sup |Ĝ^std_tt′| σ/√(s_ts_t′) = 0.0555;
- inf λ_min(Ĝ_tt)/(1 − θ_t) = 0.7503 for both parts.
So H = 1/2 has room, as at one level (0.2832 and 0.7508 in Theorem EM). It stops at the size of the off-diagonal
certificate:
- with s_t, s_t′ and σ free, the entries of G_tt′ have numerators of 70 657 to 90 907 terms over a common denominator
  of 101 315 terms in five variables, so the squared Frobenius condition would have on the order of 10⁸ terms;
- before the orthant substitution the numerators have coefficients of both signs (`ex_signs_one.py`), so no linear,
  sign-based bound is available;
- the blocks cannot be bounded term by term either: G_z[ρ_t, ρ_t′] has the factor s_ts_t′ (`ex_signs_one.py`), while
  the pair term alone is of order s_t + s_t′, so the bound needs the exact cancellation between the pair term and
  C_tᵀC⁻¹C_t′.

**Other attempts.**
- A first run of the u₁, u₂ ≥ 2 certificate in a single process exceeded 2.5 GB at the (s, 1) condition. This led to
  late composition, the split Frobenius budget and one process per condition group.
- A direct exact route for mixed T-values, the 6 × 6 T-invariant part of the S_{m₁} × S_{m₂} splitting
  (`abandoned/prove_two.py`), did not finish its second pivot in five minutes, and Lemma W replaced it. Its other
  conditions, D_i − O_ii ≻ 0 and D_i + (m_i − 1)O_ii ≻ 0 for u₁ = 1, did certify.

## 9. Files

| File | Content |
|---|---|
| `run_all.sh`, `guard.sh` | reruns every step, one Python process at a time under `guard.sh` (killed above 2 500 000 KB or after its time limit); PASS or FAIL per step, ALL PASS at the end; one log `logs/run_<step>.log` per step. Sections 1–6 of the script: 86 steps, 44 minutes on a loaded machine (complete run: `logs/run_all_frozen.out`). Section 7 (Section 10 of these notes): 9 steps, 63 minutes; `NEW_ONLY=1 bash run_all.sh` runs only section 7, and `logs/run_all.out` is such a run |
| `em3lib.py`, `em3det.py` | Section 10: the matrix A, the factors, D*, the degree bound, the slices (`em3lib.py`); the steps bounds, slices, interp and final (`python em3det.py bounds\|slices\|interp s\|1`, `python em3det.py final`); data in `em3_out/` |
| `em3points.py`, `em3neg.py` | Section 10.4: Ñ/D_c against det A computed directly far outside the grid; the negative controls |
| `em3probe.py`, `em3probe2.py`, `em3prof.py` | exploration for Section 10 (sizes, factor structure, the cost of a slice in the two sets of variables, which agree; no claim rests on them); logs in `logs/em3/` |
| `tl_m.py` → `tl_*.pkl` | Sections 2–3: the closed forms (`python tl_m.py inv\|std1\|std2 CASE [T11\|T12\|T22]`; CASE one, gen, fN, two, onetwo) |
| `check_tl.py` | Section 3: exact validation of the closed forms against direct elimination |
| `prove_eqT.py` | Section 4: Theorem M1 (`python prove_eqT.py [eps]`) |
| `prove_gersh.py` | Section 5: Theorem TV (`python prove_gersh.py one\|gen\|fN [Hnum Hden] [only=TAG] [bz]`) |
| `check_fN.py`, `check_e2e2.py`, `controls2.py` | Section 7 (the fixed-u₁ forms against the symbolic ones, and (B_inv) at large sizes; end-to-end checks; negative controls) |
| `ex_numframe.py`, `ex_twoval.py`, `ex_factor_one.py`, `ex_signs_one.py` | exploration: the Gershgorin constants for arbitrary and two-valued T-weights; factor and sign structure of the u₁ = 1 forms |
| `cf.py`, `flintfield.py`, `fastj1.py`, `j1lib.py`, `fastdirect.py`, `emlib.py`, `numG.py` | copied from `relative-anchors/code/typeI_G_em/` (closed form of the block, rational functions, Lemma B builder, helpers) |
| `abandoned/`, `logs/abandoned/` | the attempts of Section 8 (not run by `run_all.sh`) |
| `logs/` | output of a complete run of `run_all.sh` (all sections, from a clean copy; `logs/run_all.out`, ALL PASS) |
| `run_fN_check.sh`, `logs/fN_check/` | a check of Proposition B3 by the other method: the fixed-u₁ certificate of Section 5 (`prove_gersh.py fN`, all three pivots of (B_inv) by elimination with coefficient signs, u₂ ≥ 2 symbolic) at u₁ = 17, 18, 20, 24, 30, 40, 50, 77, 100, 200 and 500, all certified (`tl_m.py` has these N in its case list) |

## 10. Proposition B3: the third pivot of (B_inv) by determinant and interpolation

**Proposition B3.** Let s_t ∈ {s, 1} with s > 0, let r > 0 and σ = s_t + r, and let u₁, u₂ ≥ 2 be real. Let A be
the 3 × 3 matrix of (B_inv) in the frame of Section 5,

    A = G_y,tt − H (r/σ) diag(1, κ₁², κ₂²),   H = 1/2,

with G_y,tt = NᵀG_z,tt N from the closed form `tl_inv_gen_T11` (s_t = s) or `tl_inv_gen_T22` (s_t = 1). Then
det A > 0.

The first two LDLᵀ pivots p₁, p₂ of A are positive (Section 5: `prove_gersh.py gen only=B_s12` and `only=B_112`), and
the leading 2 × 2 minor is p₁p₂. So the third pivot det A/(p₁p₂) is positive and A ≻ 0. This is the condition
(B_inv) of Lemma W in congruent form: A = K^{1/2}(Ĝ^inv_tt − H(1 − θ_t)I)K^{1/2} with K = diag(1, κ₁², κ₂²), and
1 − θ_t = r/σ. No base point is needed: det A > 0 alone gives the third pivot.

### 10.1 The argument

Orthant variables: u_i = 2 + x_i, σ = s_t + r, s₁ = s. The region is 𝓡 = {x₁, x₂ ≥ 0, r > 0, s > 0}. The steps are
checked by `em3det.py` (modes bounds, slices, interp and final).

1. *The matrix.* `em3lib.build_A` forms A in the variables (u₁, u₂, σ, s₁) of the closed form, with the same lines as
   `prove_gersh.run_inv`, and composes each entry into (x₁, x₂, r, s). The change of variables is affine and
   invertible, so this gives the same rational functions as the late composition of `prove_gersh.py`. The first two
   pivots of this A are equal to the functions certified there (checked exactly in `em3det.py final`, and certified
   again there).
2. *The entry denominators.* They factor over Q into nine irreducible polynomials f₀, …, f₈. Each has positive
   coefficients and a monomial free of x₁, x₂, r. So every entry of A is continuous on 𝓡, and no entry denominator
   vanishes on a slice x₁ = c ≥ 0. For s_t = s their degrees in x₁ are 0, 1, 14, 1, 1, 1, 2, 1, 1; the largest, f₂,
   has 3 060 terms and degree 14 in each variable. For s_t = 1 there are again nine irreducible factors, with the same
   two properties and the same exponents e_f (`logs/run_em3_bounds_1.log`).
3. *A multiple of the denominator and the degree bound.* For each f let
   e_f = max over permutations p of Σ_a mult_f(den A_{a,p(a)}) (exact dynamic programming over subsets,
   `em3lib.maxperm`), and D* = Π f^{e_f}. Every Leibniz term D* Π_a A_{a,p(a)} equals
   (D*/Π_a den A_{a,p(a)}) Π_a num A_{a,p(a)} and is a polynomial. So T = D* det A is a polynomial, and

       deg_{x₁} T ≤ B = deg_{x₁} D* + max_p Σ_a (deg_{x₁} num A_{a,p(a)} − deg_{x₁} den A_{a,p(a)}).

4. *The candidate denominator.* D_c = Π f^{d_f}, where d_f is the largest exponent of f in the denominator of det A
   on the slices x₁ = 0, 7 and 31. Then d_f ≤ e_f, so R = D*/D_c is a polynomial. D_c is a product of the f, so
   D_c > 0 on 𝓡.
5. *Slices.* For c = 0, 1, …, B, det A at x₁ = c is computed exactly. The entries are specialised at u₁ = 2 + c in the
   variables of the closed form, the determinant is taken over Q(u₂, σ, s₁) (as M₀₀ times the determinant of the
   Schur complement of M₀₀), and the reduced result is composed into (x₂, r, s) (`em3lib.OrigSlicer`). This equals the
   determinant computed in the orthant variables (checked exactly at x₁ = 0 and 40, `em3probe2.py`). Then
   P_c = D_c(c, ·) det A(c, ·) must be a polynomial in (x₂, r, s); this is an exact division of D_c(c, ·) by the
   denominator of the slice.
6. *Interpolation.* For every monomial x₂^a r^b s^e that occurs in some P_c, its B + 1 coefficients are interpolated in
   x₁ at the nodes 0, …, B, coefficient by coefficient, with the inverse of the Vandermonde matrix over Q. Every
   coefficient of x₁^l with l > L = B − deg_{x₁} R must vanish. The result is a polynomial Ñ(x₁, x₂, r, s) of degree at
   most L in x₁ with Ñ(c, ·) = P_c for c = 0, …, B.
7. *The identity.* T − ÑR has degree at most B in x₁, since deg_{x₁}(ÑR) ≤ L + deg_{x₁} R = B. It vanishes at every
   node: specialising each Leibniz term at x₁ = c is allowed (no entry denominator vanishes there), so
   T(c, ·) = D*(c, ·) det A(c, ·) = R(c, ·) P_c = R(c, ·) Ñ(c, ·). A polynomial of degree at most B in x₁ that
   vanishes at B + 1 values of x₁, as a polynomial in (x₂, r, s), is zero. So T = ÑR and det A = Ñ/D_c identically.
8. *Signs.* Every coefficient of Ñ is ≥ 0, and some coefficient of a monomial free of x₁, x₂ and r is > 0. So Ñ > 0 on
   𝓡 (even for r ≥ 0), and with D_c > 0, det A > 0 on 𝓡. ∎

### 10.2 The numbers

| | s_t = s (T11) | s_t = 1 (T22) |
|---|---|---|
| numerator terms of the entries of A | 22 605 to 34 712 | the same |
| irreducible factors of the entry denominators | 9 | 9 |
| exponents e_f of D* | 3, 3, 3, 4, 6, 6, 6, 4, 1 | 3, 3, 3, 4, 6, 6, 6, 4, 1 |
| exponents d_f of D_c (seen on x₁ = 0, 7, 31) | 3, 1, 1, 2, 2, 4, 6, 2, 1 | 3, 1, 1, 2, 4, 2, 6, 2, 1 |
| deg_{x₁} of D*, D_c, R | 78, 38, 40 | 78, 38, 40 |
| B, the bound for deg_{x₁} T (the bounds in x₂, r, s) | 76 (76, 77, 82) | 76 (76, 77, 76) |
| L = B − deg_{x₁} R | 36 | 36 |
| grid | x₁ = 0, …, 76 (u₁ = 2, …, 78): 77 slices | the same |
| det A on each slice: numerator, denominator | 11 910 and 12 996 terms | the same |
| P_c = D_c det A on each slice | a polynomial with 11 910 terms, all positive | the same |
| time of one slice; of the step | 10 to 25 s; 1 553 s | 10 to 23 s; 1 441 s |
| monomials x₂^a r^b s^e interpolated | 11 910 | 11 910 |
| nonzero coefficients above x₁^L | none | none |
| Ñ: terms; negative coefficients; degrees in (x₁, x₂, r, s) | 139 005; none; (36, 36, 39, 42) | 139 005; none; (36, 36, 39, 36) |
| a positive coefficient free of x₁, x₂, r | yes | yes |
| time of the interpolation step | 7 s | 6 s |

Times are wall-clock times on a machine shared with other heavy jobs, one Python process at a time, each under 2.5 GB
(`logs/run_em3_*.log`). The whole of section 7 of `run_all.sh` (nine steps) took 63 minutes; the final assembly takes
612 s, mostly for forming the first two pivots in four variables twice (from this A, and through `prove_gersh.py` for
the comparison). Every slice polynomial P_c, for both types, has
only positive coefficients. The slice determinant is computed in the variables of the closed form because it is 1.2
to 1.4 times faster there than in the orthant variables (`em3probe2.py`: 9.0 s against 12.5 s at x₁ = 0, 13.8 s
against 16.9 s at x₁ = 40); taking the common denominator factor of the entries out first is slower (37.5 s against
14.9 s at x₁ = 40, `em3prof.py`).

### 10.3 Consequences

Proposition B3 is the condition of row 7. With Section 5, every condition of Lemma W holds for all u₁, u₂ ≥ 2 and both
types of s_t. So Theorem TV holds for all level sizes, and Corollary TV for every k (Section 6). Proposition B3 is
proved for all real u₁, u₂ ≥ 2; for min(u₁, u₂) ≤ 16 it repeats, independently, the fixed-N certificates of Section 5.

This says nothing about T-weights outside the two values of U (row 10), about three or more weight values, or about
Conjecture 12 and the NPT problem.

### 10.4 Checks and controls

- `em3det.py final` (exact): for both types, the setup recomputed from scratch equals the saved one; D_c divides D*;
  D_c and every entry denominator are positive on 𝓡; B, deg_{x₁} R and L are as saved; the first two pivots of A
  equal the functions of `prove_gersh.py gen only=B_s12` and `only=B_112`, and are certified again; every slice file
  is present and the interpolation passed. Both types pass.
- `em3points.py` (exact): at five points with x₁ outside {0, …, 76}, (u₁, u₂, s, r) = (300, 17, 2/7, 5/3),
  (1000, 1000, 1000, 1/1000), (79, 5000, 1/1000, 50), (2003/7, 41/3, 7/5, 3) and (100 000, 23, 3/11, 1/7), the closed
  form is evaluated over Q and A, det A and the three pivots are formed directly, with code independent of
  `em3lib.py` (as in `check_fN.py`, part 2). For both types det A = Ñ/D_c exactly, and all three pivots are positive.
  The values of det A there range from 5.8·10⁻¹⁴ (s_t = s at (79, 5000, 1/1000, 50)) to 1.1.
- `em3neg.py` (exact):
  1. *A larger H.* With H = 1 the third pivot fails at explicit rational points with u₁ = u₂ = 17. For s_t = s at
     (s, r) = (1/1000, 2) the first two pivots are 0.500 and 4.5·10⁻⁶ and det A = −9.5·10⁻¹²; for s_t = 1 at
     (s, r) = (1, 50) they are 0.510 and 5.2·10⁻⁴ and det A = −5.2·10⁻⁴. With H = 1/2 all three pivots are positive at
     both points. On the slice x₁ = 15 through them, det A at H = 1 has a numerator of 11 910 terms with 7 727
     (s_t = s) and 4 946 (s_t = 1) negative coefficients, and its value at the point equals the direct determinant.
  2. *The vanishing test.* With the slice x₁ = 38 multiplied by 1 + 1/1000, the interpolant of the first 1 500
     monomials has 60 000 nonzero coefficients above x₁^36, that is every coefficient of x₁^37, …, x₁^76 of every one
     of them; unchanged it has none. Both types.
  3. *The candidate denominator.* Lowering any one of the nine exponents d_f by 1, D_c det A is no longer a polynomial
     on the slice x₁ = 0. Both types.
- `em3probe2.py` (exploration, exact): the slice computed in the variables of the closed form equals the slice computed
  in the orthant variables, at x₁ = 0 and 40 (s_t = s).
