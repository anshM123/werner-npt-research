# (R∞)₂ at two weight values for every k: the type I blocks with one and two letters 2

5 November 2026. Nothing here is refereed.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | The torus blocks of (R∞)₂ treated here, their sizes and the kernel vectors p ⊗ y: family 1 (type I, \|T\| = 2, one letter 1), 4k − 4 states, kernel p ⊗ \|2@t 2@τ⟩; family 2 (type I, \|T\| = 1, two letters 1), k(k − 1)/2 + k(k − 1) + (k − 1)(k − 2)/2 states, kernel p ⊗ \|2@t 1@c⟩ (c ∈ U); family 3 (the letter-1 block), k·k(k − 1)/2 + k(k − 1)(k − 2)/6 states, kernel p ⊗ \|1@P⟩ (P a pair). At k = 5, 6, 7 (two exact rational weight vectors each, one for family 3 at k = 7; two U-weight levels, arbitrary T-weights in family 1) the nullity of each block equals the number of these vectors, so nothing else appears, and the gauge-fixed block is positive definite | **exact check** at points (`explore_kernel.py`, `run_all.sh` steps `kernels`, `kernels_k7`); not a proof |
| 2 | **Lemma P2** (family 1, P-states first). For every k ≥ 3, every T = {t, τ} and every positive weight vector, the family-1 block is ⪰ 0 with kernel exactly p ⊗ \|2@t 2@τ⟩ if and only if the explicit matrix Φ = C − Bᵀ diag(1/N) B of size 3k − 4 is positive definite | **proved** (by hand, Section 2); identities **checked exactly** (`check_t2.py`) |
| 3 | Closed forms of all ten sub-blocks of family 1 (Lemma F1), all weights | **proved** (by hand from Lemma B, Section 2); **exact check** entry by entry against the Lemma B builder `j1lib.py` at 6 random rational weight vectors (k = 4, 5, 6; every copy of U treated as its own orbit, so no symmetry is used), with the B-part scaled by 1 and by 1001/1000; `j1lib.py` = `ylib.py` on the blocks (`check_ylib_r2.py`) |
| 4 | **Theorem F1a.** Family 1 for every k ≥ 4 when the weights of U = [k] ∖ {t, τ} are all equal, with **arbitrary** weights on t and τ: ⪰ 0 with kernel exactly p ⊗ \|2@t 2@τ⟩ | **proved**, computer-assisted, exact, symbolic in u = k − 2 ≥ 2 and in the three weights (Section 2.4) |
| 5 | **Theorem F1b.** Family 1 for every k ≥ 5 whenever all k weights take at most two values **and** either t and τ have equal weights, or U lies in one level, or one of the two U-levels has at most M_1 = 40 copies (M_1 = 12 in `run_all.sh`, extended to 40 by `run_fext.sh`, Section 5b). In particular every two-value weight vector with 5 ≤ k ≤ 83 | **proved**, computer-assisted, exact (Section 2.4: regimes, certificates, coverage) |
| 5b | The remaining two-value configurations of family 1: t and τ with different weights and both U-levels with more than M_1 = 40 copies (possible only for k ≥ 84). The certificate of the regime tw = (0, 1), (2+, 2+) passed 7 of the 8 leading minors of the invariant form (minor 7: numerator degree 158 with 91 162 terms) and was stopped by hand during minor 8 after about 23 minutes in all (memory up to 1.8 GB), to leave time for the other checks (`logs/attempts/stopped_regimes.txt`) | **open** beyond M_1 (Section 2.4) |
| 6 | Family 1 with arbitrary weights on t and τ and U-weights on two levels (the full target for this family) | **open**. Attempted: 8 × 8 invariant form in five variables; the standard part passed its first two minors (53 124 terms) and the run was stopped during the third (Section 2.6) |
| 7 | Closed forms of all six sub-blocks of family 2 (Lemma F2), all weights | **proved** (by hand from Lemma B and Woodbury, Section 3); **exact check** entry by entry against `j1lib.py` at 6 two-level weight vectors (k = 4, 5, 6), B-part scaled by 1 and 1001/1000 (`check_f2.py`) |
| 8 | Symmetry reduction of family 2 (Lemma F2.2): parts, test vectors, pattern sums | **proved**; **exact checks**: reduced blocks against explicit projections of the `j1lib.py` block at 13 size pairs (both scalings, irreps on either level), kernels, and dimension counts for every size pair with 5 ≤ k ≤ 60 (`check_f2.py`, `dims.py`) |
| 9 | **Theorem F2a.** Family 2 for every k ≥ 4 when the U-weights are all equal, with **arbitrary** weight on t: ⪰ 0 with kernel exactly span{p ⊗ \|2@t 1@c⟩ : c ∈ U} | **proved**, computer-assisted, exact (Section 3.3) |
| 10 | **Theorem F2b.** Family 2 for every k ≥ 5 whenever all k weights take at most two values **and** U lies in one level or one of the two U-levels has at most M_2 = 40 copies (M_2 = 12 in `run_all.sh`, extended to 40 by `run_fext.sh`, Section 5b). In particular every two-value weight vector with 5 ≤ k ≤ 82 | **proved**, computer-assisted, exact (Section 3.3: the two-level regimes, certificates, coverage) |
| 10b | The remaining two-value configurations of family 2: both U-levels with more than M_2 = 40 copies (possible only for k ≥ 83). Only the invariant part is missing there: its certificate in the regime (3+, 3+) passed 7 of its 9 leading minors (the last one with 109 975 terms) and was stopped by the 2.5 GB limit during minor 8 (`logs/attempts/stopped_regimes.txt`). The fixed-size regimes (m, 3+) and (3+, m) of `run_all.sh` and `run_fext.sh` give M_2 | **open** beyond M_2 (Section 3.5) |
| 11 | Family 2 with arbitrary weight on t and U-weights on two levels (the full target for this family) | **open**. Attempted (four variables): the R-block and the vacuum block of the invariant part are certified, the run was stopped while eliminating the vacuum states (Section 3.5) |
| 12 | Family 3 (the letter-1 block of (R∞)₂) at two weight values | **open**; only the exact checks of row 1 (k = 5, 6, 7). Not attempted symbolically |
| 13 | Negative controls: with B → (1 + 1/1000)B the gauge-fixed block is not positive definite at explicit rational points (families 1 and 2, including points where all weights take two values), and the certificates of the one-level cases fail | **exact** (`neg_t2.py`, `neg_f2.py`, `run_all.sh` steps `*_neg_*`) |
| 14 | Consequence: for 5 ≤ k ≤ 26 and every weight vector with at most two distinct weights, every torus block of (R∞)₂ **except the letter-1 block** is ⪰ 0 with kernel exactly the part of p ⊗ H₂ it contains; for larger k the same holds except also for the blocks of rows 5b and 10b | **proved**, given rows 5 and 10 and Theorem DJ and Proposition T2′ of the README (themselves computer-assisted, not refereed); (R∞)₂ itself is **not** proved for any k ≥ 5, since the letter-1 block is open |
| 15 | Final run of `run_all.sh` (one Python process at a time; the two regimes of rows 5b and 10b are not part of it, see `attempts.sh`) | see the last line of `logs/run_all.out` and Section 6 |
| 16 | A second check written separately (Section 7): for family 1 the Schur complement of Lemma P2 taken directly from the block built from the definitions, projected, equals the certified symbolic forms exactly at 7 points (k = 7, 8, weight ratios from 1/900 to 40); for family 2 the spectra of the certified parts, with the irrep dimensions as multiplicities, add up to the spectrum of that block at 5 points (k = 7, 8); kernels and gauge-fixed positivity exactly at these 12 blocks | **exact check** at points (`crosscheck.py`, `run_crosscheck.sh`); the spectral comparison in floating point |

Nothing here says anything about Conjecture 12 of `single-excitation-range/` in general, or about the NPT (Werner
distillability) problem. Computer-assisted means: exact rational-function arithmetic in python-flint and sign checks of
coefficients; no floating point enters any check or proof.

## 1. Setting

Notation as in `relative-anchors/README.md`, Sections 1 and 2 (repository `werner-npt-research`). The group has k
copies, s_c = w_c², L = |s|. (R∞)₂ is Y compressed to G1₂ ⊕ G0₃; it preserves the number of letters 1 and, at each
copy, the number of letters 2. Every state has three letters in row and column together. Its torus blocks are:
the letter-1 block (family 3); the type I blocks with |T| = 1 (family 2), |T| = 2 (family 1) and |T| = 3 (Theorem DJ
of the README, all weights); the type II blocks (Proposition T2′, all weights). We work with the entries of Lemma B
in units 2^{k−2} (the `entry` of `j1lib.py` with Z = k − 2 on G1₂ and the same unit on G0₃) and in the rescaled frame
of `ylib.py`: each state is multiplied by the product of the weights over its letters 1, so that every entry is
rational in the squared weights and homogeneous of degree 1 in them. Positive diagonal rescaling keeps inertia, and
the kernel vectors p ⊗ y become indicators of column groups.

Two weight values: the copies outside T form levels U₀, U₁ (sizes u₀, u₁) with squared weights 1 and q (q > 0; the
homogeneity makes the normalisation 1 harmless). In family 1 the copies t, τ of T carry squared weights a, b; in
family 2 the copy t carries a. "Two weight values in all" means a, b (or a) ∈ {1, q}.

Two symmetries reduce the cases (used in `coverage.py`):
- (S1) dividing all squared weights by q and exchanging the names of the levels maps two-value configurations to
  two-value configurations (q ↦ 1/q, sizes exchanged, a T-copy of level i goes to level 1 − i, a part on level i goes
  to the same part on level 1 − i);
- (S2) in family 1, exchanging t and τ permutes the states (R ↔ R′) and exchanges a and b.

## 2. Family 1: T = {t, τ}

### 2.1 States and closed forms

U = [k] ∖ {t, τ}, u = k − 2. States (4k − 4):
- R_c = (2@t, 2@τ 1@c), c ∈ U ∪ {t};
- R′_c = (2@τ, 2@t 1@c), c ∈ U ∪ {τ};
- P_r = (1@r, 2@t 2@τ), r ∈ [k];
- Q_c = (vac, 2@t 2@τ 1@c), c ∈ U.

Write J^±_{(x)} for J^± with w_x replaced by −w_x, ε_x = −1 at the flipped copy and +1 elsewhere, and, on any
single-letter sector with pair energies C(x, z), N_{XY} = J⁻_X C⁻¹ J⁺_Y and A = I + ½N_{∅∅}.

**Lemma F1 (closed forms; every k ≥ 3, all positive weights).** Unrescaled, in units 2^{k−2}:
1. Y_RR = diag(2 − [c = t]) and Y_R′R′ = diag(2 − [c = τ]).
2. Y_PP = Π = diag(N) − (N∘w)(N∘w)ᵀ/(2L − a − b), with N_r = 2 on U and 1 on T.
3. Y_RR′[c, c′] = ½[c = c′ ∈ U] − ¼(N_{(τ)(t)} − ½N_{(τ)∅}A⁻¹N_{∅(t)})[c, c′], on the singles of [k] with
   C(x, z) = 2(L − s_x − s_z).
4. Y_RP[c, r] = (½ − [r = t])[c = r] − ½[c = r]·ρ − ½ν_{r,c} w_r w_c (ε_c − ρ)/Δ_c, with Δ_x = 2L − b − 2s_x
   (x ≠ τ), γ = Σ_{x≠τ} s_x/Δ_x, γ̃ = Σ_{x≠τ} ε_x s_x/Δ_x (flip at t), ρ = γ̃/(2 + γ), ν_{r,c} = 2 if r ∉ {τ, c},
   1 if r = τ, 0 if r = c. Y_R′P: the same with t ↔ τ, a ↔ b.
5. Y_QQ = 2(I + ½N_Q)⁻¹ with N_Q = J⁻C₄⁻¹J⁺ on {|2@t 2@τ 1@c⟩}, C₄ = 2L − a − b − 2s_x − 2s_y on
   |2@t 2@τ 1@{x, y}⟩.
6. Y_RQ = −E + ½A_τ⁻¹N^τ_{∅(t)}, on the singles {|2@τ 1@x⟩ : x ≠ τ} with C′(x, z) = 2L − b − 2s_x − 2s_z,
   A_τ = I + ½N^τ_{∅∅}, E the embedding of U. Y_R′Q: t ↔ τ, a ↔ b.
7. Y_PQ[r, c] = −[r = c] + (2/σ_Q)(½[r = c]γ_Q + ν′_r w_r w_c/f_c), f_x = 2L − a − b − 2s_x,
   γ_Q = Σ_{x∈U} s_x/f_x, σ_Q = 2 + γ_Q, ν′_r = 2 on U ∖ c, 1 on T, 0 at r = c.

*Proof.* Lemma B with j = 2.
- 𝒯_{(2,t)} needs a letter 2 at t, and the columns of the R-states have none there; so every B-term with an R-state
  and its own row vanishes, and the swap term between two R-states vanishes for the same reason. This gives item 1.
- 𝒯_{(2,τ)}|2@τ 1@c⟩ = ½J⁺|1@c⟩ − w_τ|1@{τ, c}⟩ = ½J⁺_{(τ)}|1@c⟩, and likewise for t. The B-term between R_c and
  R′_{c′} is ¼⟨J⁺_{(τ)}e_c, E₂⁻¹J⁺_{(t)}e_{c′}⟩ on the letter-1 pairs; with E₂ = C + ½J⁺J⁻ and Woodbury,
  E₂⁻¹ = C⁻¹ − ½C⁻¹J⁺A⁻¹J⁻C⁻¹, which gives item 3. The swap term is ½[c = c′ ∈ U].
- On |2@t 2@τ⟩ the operator 𝒯_{(1,r)} is w_r N_r and E = 2L − a − b: item 2.
- Between R_c and P_r: the swap term is (½ − [r = t])[c = r]; the B-term lives in the sector of |2@τ 1@x⟩, where
  E = diag(Δ) + ½wwᵀ (Sherman–Morrison), 𝒯_{(2,t)}|2@t 2@τ⟩ = ½w̃ (w with w_t ↦ −w_t) and
  𝒯_{(1,r)}|2@τ 1@c⟩ = ½[r = c]w + ν_{r,c}w_r e_c. This gives item 4.
- Item 5 is Lemma G. Item 6: σ⁻_{(2,t)} gives −[c = c′]; 𝒯_{(2,t)}|2@t 2@τ 1@c′⟩ = ½J⁺_{(t)}|2@τ 1@c′⟩ and
  J⁻E₃⁻¹ = A_τ⁻¹J⁻C′⁻¹ (from (I + ½J⁻C′⁻¹J⁺)J⁻ = J⁻C′⁻¹E₃). Item 7: σ⁻_{(1,r)} gives −[r = c], and the B-term lives
  in the sector of |2@t 2@τ 1@x⟩ (Sherman–Morrison). ∎

`check_t2.py` (a) evaluates the closed forms with every copy of U in its own level (then the reduced matrices of the
algebra of Section 2.3 are the full blocks) and compares all ten sub-blocks entry by entry with the rescaled
`j1lib.py` block at six random rational weight vectors with k = 4, 5, 6, with the B-part scaled by 1 and 1001/1000:
all equal.

### 2.2 Lemma P2

**Lemma P2.** For every k ≥ 3, every T = {t, τ} and every positive weight vector, order the states P first and write
the block as [[Π, B], [Bᵀ, C]] (C on the R-, R′- and Q-states). Then the block is ⪰ 0 with kernel exactly
p ⊗ |2@t 2@τ⟩ if and only if Φ = C − Bᵀ diag(1/N) B (size 3k − 4) is positive definite.

*Proof.* As Lemma P of `middle_j1_t1/NOTES.md`. By item 2, Πw = 0, and (w, 0) = p ⊗ |2@t 2@τ⟩ is in the kernel
(Lemma N), so Bᵀw = 0. Since Σ_r N_r s_r = 2L − a − b, every b ⊥ w has Π diag(1/N) b = b − (N∘w)(wᵀb)/(2L − a − b)
= b, so Π diag(1/N) B = B and the congruence K = [[I, −diag(1/N)B], [0, I]] gives KᵀMK = diag(Π, Φ). Π =
N^{1/2}(I − vvᵀ)N^{1/2} with v = N^{1/2}w/|N^{1/2}w| is ⪰ 0 with kernel span(w). ∎

`check_t2.py` (b) checks Πw = 0, Bᵀw = 0 and Π diag(1/N) B = B exactly at the six points of (a).

### 2.3 The symmetry reduction

The group G = S_{u₀} × S_{u₁} permutes the copies inside each level; t and τ are fixed. The coordinates of Φ are
indexed by single copies: R by U ∪ {t}, R′ by U ∪ {τ}, Q by U. Each such space is the sum of one trivial
representation per orbit ({t}, {τ}, U₀, U₁) and one standard representation per level of size ≥ 2.

**Lemma F1.2.** Φ ≻ 0 if and only if its invariant part Φ_inv (on the orbit indicators: (R_t, R_{U₀}, R_{U₁},
R′_τ, R′_{U₀}, R′_{U₁}, Q_{U₀}, Q_{U₁}), 8 × 8 for two levels, 5 × 5 for one) and, for each level of size ≥ 2, its
standard part (3 × 3 on (R, R′, Q)) are positive definite.

*Proof.* Schur's lemma for the decomposition above; the rescaling D = diag(w) is constant on orbits. ∎

Computation (`t2lev.py`): an invariant operator between two single-index spaces is stored by its reduced matrix
M[o_y][o_x] = Σ_{x ∈ o_x} O[y, x] and its scalar O[h, h] − O[h, h′] on each standard part. Products are products of
reduced matrices, transposes use the Gram matrices diag(orbit sizes), inverses are inverses of the reduced matrices
and of the scalars; in the rescaled frame (O₁O₂)~ = Õ₁D⁻²Õ₂ and (O⁻¹)~ = D²Õ⁻¹D². The closed forms of Lemma F1 are
entered as functions of the orbit labels of the two indices and of their equality, with every sum over a third copy
written with the orbit sizes. Level sizes are fixed (1) or symbolic (u = 2 + x, x ≥ 0 real). All divisions are by
quantities that are positive at admissible sizes (energies of existing strings, reduced matrices and scalars of
positive definite operators, orbit sizes), and the off-diagonal value of a level of fixed size 1 is never used, so
the symbolic forms specialise correctly at every admissible integer size.

`check_t2.py` (c): at 12 size pairs ((1,1), (2,1), (1,2), (2,2), (3,1), (3,2), (2,3), (4,1), (3,3), (2,0), (3,0),
(4,0)) with random rational weights and random placement of the levels, the invariant and standard forms of Φ and of
the gauge-fixed block (P_t deleted; B scaled by 1 and 1001/1000) equal the projections of the `j1lib.py` block on the
explicit test vectors (orbit indicators, and e_h − e_{h′} on a level), exactly.

### 2.4 Certificates

`cert_t2.py SPECS [tw=…]` builds Φ_inv and the standard forms over Q(variables) and checks every leading principal
minor (one pass of symmetric elimination, D_m = D_{m−1}·pivot): numerator and denominator must be products of
irreducible polynomials with one-signed coefficients, each with a monomial free of the x_i, of total sign +1. Then the
minor is positive for all a, b, q > 0 and x_i ≥ 0. In every certified regime all numerators and denominators were
one-signed as they stand (no factorisation needed). By Sylvester's criterion the forms are positive definite on the
whole regime, and Lemmas P2 and F1.2 give the block statement.

**Theorem F1a.** For every k ≥ 4, every T = {t, τ} and all weights with w_c equal for all c ∈ U (w_t, w_τ arbitrary),
the family-1 block is ⪰ 0 with kernel exactly p ⊗ |2@t 2@τ⟩.

*Proof.* Regime u = 2 + x, variables (a, b, x) (`cert_t2.py 2+`): the 5 invariant and 3 standard minors pass; the
largest has numerator degree 40 with 10 425 terms. ∎

**Theorem F1b.** Let k ≥ 5, T = {t, τ}, and let the weights take at most two distinct values. Suppose that t and τ
have equal weights, or that U lies in one level, or that one of the two U-levels has at most M_1 copies (M_1 as in
row 5 of the status table). Then the family-1 block is ⪰ 0 with kernel exactly p ⊗ |2@t 2@τ⟩. Both U-levels have more
than M_1 copies only if k − 2 ≥ 2M_1 + 2, so every two-value weight vector with k ≤ 2M_1 + 3 is covered.

*Proof.* If all U-weights are equal this is Theorem F1a. Otherwise U meets both levels, and by (S1), (S2) the T-weights
can be taken as (a, b) = (1, 1) or (1, q). Certified regimes (variables q and the x_i; all in `run_all.sh`): tw = (0, 0) at
(u₀, u₁) = (1, 2+), (2+, 1), (2+, 2+), and tw = (0, 1) at (1, 2+) and at (m, 2+) for
2 ≤ m ≤ M_1 (for tw = (0, 1), (S1) followed by (S2) exchanges u₀ and u₁, so these cover min(u₀, u₁) ≤ M_1).
`coverage.py` checks that every configuration of the statement with 5 ≤ k ≤ 60 lies in the orbit of a certified
regime (the regimes are uniform in the sizes, so this is bookkeeping). Largest minors: the determinant of the 8 × 8
invariant form for tw = (0, 0), (2+, 2+), numerator degree 119 with 40 932 terms over a denominator of degree 110 with
30 856 terms (`logs/run_f1_cert_A_2_2.log`, about two minutes). All numerators and denominators are one-signed as they
stand. ∎

*Not covered.* The regime tw = (0, 1), (2+, 2+) (t and τ with different weights, both U-levels with at least two
copies): in the recorded attempt its first 7 leading minors pass (minor 7: numerator degree 158 with 91 162 terms over a
denominator of degree 150 with 76 076 terms, after about 12 minutes); minor 8, the determinant, had not finished after
about 11 more minutes (memory up to 1.8 GB) and the run was stopped by hand (`attempts.sh`; the lines it displayed are in `logs/attempts/stopped_regimes.txt`).
So the configurations with t and τ of different weights and both U-levels larger than M_1 remain open.

### 2.5 Controls

- `neg_t2.py`: with B → (1 + 1/1000)B the reduced forms of the gauge-fixed block (P_t deleted) fail at 74 of 126
  candidate points while the true ones pass; at 9 of them (one level, two levels, and one point with all weights
  taking two values, w = (1/30, 1, 1, 1, 1/30)) the full scaled gauge-fixed block built by `j1lib.py` is not positive
  definite while the true one is (exact LDL).
- `run_all.sh` step `f1_neg_cert_one`: the certificate of the gauge-fixed block at one level fails with the scaled
  B-part (a minor with a factor of mixed signs); step `f1_gauge_cert_one`: the same certificate passes for the true
  block (a second proof of Theorem F1a without Lemma P2).
- For the block itself B → (1 + 1/1000)B always fails, trivially, on the kernel vector; the controls are on the
  gauge-fixed block.

### 2.6 Not covered: arbitrary T-weights with two U-levels

With a, b free and two U-levels there are five variables (a, b, q, x₀, x₁). The standard part of level 0 passed its
first two minors (numerator degree 40, 53 124 terms) and the third had not finished after several minutes; the
invariant 8 × 8 form would be far larger. Stopped (`logs/dev/dev_cert_t2_2lev_std.log`). So for family 1 the stated
target (arbitrary T-weights, U-weights on two levels) is **not** reached; what is proved is Theorem F1a (arbitrary
T-weights, one U-level) and Theorem F1b (two values in all).

## 3. Family 2: T = {t}

### 3.1 States and closed forms

U = [k] ∖ {t}. States:
- R_P = (2@t, 1@P), P a pair of [k];
- P_{r,c} = (1@r, 2@t 1@c), r ∈ [k], c ∈ U;
- Q_P = (vac, 2@t 1@P), P a pair of U.

Notation: C(P) = 2(L − s_P) on letter-1 pairs, C₃(T) = 2(L − s_T) on letter-1 triples, C′(P) = 2L − a − 2s_P on the
pairs of the t-sector, C₄(T) = 2L − a − 2s_T on its triples; Δ_x = 2L − a − 2s_x, γ = Σ_{x∈U} s_x/Δ_x, σ = 2 + γ;
ν_{r,c} = 2 if r ∉ {t, c}, 1 if r = t, 0 if r = c; ε_x = −1 at x = t. On the singles of [k]: N = J⁻C⁻¹J⁺, A = I + ½N,
Ǎ = D A⁻¹ D, and M̃ = (A⁻¹N_{∅(t)})~; on the singles of U in the t-sector: A_t = I + ½J⁻C′⁻¹J⁺ and Ǎ_t = D A_t⁻¹ D;
on the letter-1 pairs of [k]: N₂ = J⁻C₃⁻¹J⁺, Â₂ = D⁻¹(I + ½N₂)D⁻¹, Ň₂ = D⁻¹N₂_{∅(t)}D⁻¹; on the pairs of U in the
t-sector: Ŵ₄ = D⁻¹(I + ½J⁻C₄⁻¹J⁺)D⁻¹.

**Lemma F2 (closed forms; every k ≥ 4, all positive weights).** Rescaled, in units 2^{k−2}:
1. RR: diagonal, s_P(2 − [t ∈ P]).
2. PP: (½ − [r = r′])[c′ = r][c = r′]s_r s_{r′} + [r = r′][c = c′](2 − [r ∈ {t, c}])s_r s_c − β, with
   β = [r′ = c][r = c′] s_c s_{c′} γ/(2σ) + [r′ = c] s_r s_c s_{c′} ν_{r,c′}/(Δ_{c′}σ) + [r = c′] s_{c′} s_{r′} s_c ν_{r′,c}/(Δ_c σ)
   + s_r s_{r′} ν_{r′,c} ν_{r,c′}([c = c′] s_c/Δ_c − s_c s_{c′}/(Δ_c Δ_{c′} σ)).
3. RP: (½ − [r = t])[P = {c, r}] s_r s_c − ¼[r ∈ P] s_r M̃[P ∖ r, c]
   − [r ∉ P](s_r/C(P))([c ∈ P] ε_{P∖c} s_{P∖c} s_c − ½(s_z M̃[x, c] + s_x M̃[z, c])), P = {x, z}.
4. PQ: −[r ∈ P] s_r Ǎ_t[c, P ∖ r] + [r ∉ P] ν_{r,P} s_r (s_z Ǎ_t[c, x] + s_x Ǎ_t[c, z])/C′(P), ν = 1 at r = t, 2 on U ∖ P.
5. RQ: −[P = P′] s_P + ½(Â₂⁻¹ Ň₂ D²)[P, P′].
6. QQ: 2Ŵ₄⁻¹.

*Proof.* Lemma B with j = 2, as in Section 2.1. The R-states have no letter 2 at t in the column, so RR is diagonal
with the values 2 − n_t. PP: the B-term lives in the sector of |2@t 1@x⟩, E = diag(Δ) + ½w_Uw_Uᵀ (Sherman–Morrison),
with 𝒯_{(1,r)}|2@t 1@c⟩ = ½[r = c]w_U + ν_{r,c}w_r e_c. RP: 𝒯_{(2,t)}|2@t 1@c⟩ = ½J⁺_{(t)}|1@c⟩ and
𝒯_{(1,r)}|1@P⟩ = ½[r ∈ P]J⁺|1@P∖r⟩ + 2[r ∉ P]w_r|1@P⟩; with E₂⁻¹ = C⁻¹ − ½C⁻¹J⁺A⁻¹J⁻C⁻¹ one gets
J⁻E₂⁻¹J⁺_{(t)} = A⁻¹N_{∅(t)} and E₂⁻¹J⁺_{(t)} = C⁻¹(J⁺_{(t)} − ½J⁺A⁻¹N_{∅(t)}). PQ: J⁻E₃⁻¹ = A_t⁻¹J⁻C′⁻¹ in the
t-sector and A_t⁻¹N′ = 2(I − A_t⁻¹), which cancels the A-part on r ∈ P (the same computation as the P–Q entry of
Theorem 2V). RQ: J⁻E₃⁻¹ = (I + ½N₂)⁻¹J⁻C₃⁻¹ on the letter-1 triples (k ≥ 4, so C₃ > 0) and
𝒯_{(2,t)}|2@t 1@P′⟩ = ½J⁺_{(t)}|1@P′⟩. QQ: Lemma G. With the B-part scaled by b, the A-parts are the swap and diagonal
terms of PP, the first term of RP, −[r ∈ P][c = P ∖ r]s_r s_c in PQ, −[P = P′]s_P in RQ and 2D² in QQ. ∎

`check_f2.py` (a): at (u₀, u₁) = (2,1), (1,2), (2,2), (3,1), (2,3), (3,0) with random rational weights, every entry of
RR, PP, RP, PQ (from the label formulas) and of RQ, QQ (built as explicit matrices on all pairs) equals the rescaled
`j1lib.py` block, with B scaled by 1 and by 1001/1000.

### 3.2 The symmetry reduction

G = S_{u₀} × S_{u₁}, t fixed. Parts (irrep on level 0; the parts on level 1 are the same with the levels exchanged)
and their test vectors (P: functions of (r, c); R, Q: functions of the pair; f = e_a − e_b, g = e_c − e_d with named
copies; h = e_{ac} + e_{bd} − e_{ad} − e_{bc} on four named copies of level 0):

| part | irrep | R-vectors | P-vectors | Q-vectors |
|---|---|---|---|---|
| inv | trivial | {t,U₀}, {t,U₁}, {U₀,U₀}, {U₀,U₁}, {U₁,U₁} | (t,U₀), (t,U₁), (c,c) on U₀, (U₀,U₀) off, (U₀,U₁), (U₁,U₀), (c,c) on U₁, (U₁,U₁) off | {U₀,U₀}, {U₀,U₁}, {U₁,U₁} |
| std | (std₀, triv) | f on {t,U₀}; f_x + f_z on {U₀,U₀}; f on the U₀-member of {U₀,U₁} | f_c at r = t; f_r at r = c; f_r and f_c on (U₀,U₀) off; f_r on (U₀,U₁); f_c on (U₁,U₀) | f_x + f_z on {U₀,U₀}; f on {U₀,U₁} |
| two | ((u₀−2,2), triv) | h | h(r, c) | h |
| anti | ((u₀−2,1,1), triv) | – | f₁(r)f₂(c) − f₂(r)f₁(c) | – |
| mix | (std₀, std₁) | f_x g_z + f_z g_x | f_r g_c, g_r f_c | f_x g_z + f_z g_x |

At small sizes the vectors that vanish or coincide are dropped (u₀ = 2: f_c = −f_r on (U₀,U₀) off and the
{U₀,U₀}-vectors vanish; `keep_for2`). The kernel vectors are the column indicators: in inv, (t,U₀) + (c,c)_{U₀} +
(U₀,U₀)off + (U₁,U₀) and the same for level 1; in std, f_c[r = t] + f_r[r = c] + f_c(U₀,U₀)off + f_c(U₁,U₀). Gauge:
delete the vectors on the diagonal states P_{c,c} (each kernel vector has coefficient 1 on its own gauge coordinate
and 0 on the other).

**Lemma F2.2.** The test vectors of each part are images of a seed under independent equivariant maps (Gram
determinants positive on every regime, checked symbolically by `cert_f2.py`), and the dimension counts add up to
k(k − 1)/2, k(k − 1) and (k − 1)(k − 2)/2 states (`dims.py`, every size pair with 5 ≤ k ≤ 60, small sizes
included). Hence the block is ⪰ 0 iff every reduced form is, and the kernel dimensions add. Each reduced block is a
pattern sum over the joint equality patterns of the indices of two states (named copies, anonymous copies of each
level, t named), weighted by falling factorials in the level sizes, as in Lemma 3 of `middle_j1_l1/NOTES.md`
(`f2sym.py`). Products and inverses of the pair-space operators (Â₂, Ŵ₄) are taken in the reduced bases:
R_{XY} = R_X G⁻¹ R_Y and R_{X⁻¹} = G R_X⁻¹ G. ∎

`check_f2.py` (b): at 13 size pairs ((2,1), (1,2), (2,2), (3,1), (1,3), (3,2), (2,3), (4,1), (4,2), (2,4), (3,3),
(4,0), (5,0)), both scalings, every part (irreps on either level, named copies placed on concrete copies): the
reduced blocks RR, RP, RQ, PP, PQ, QQ and the Gram matrices equal the projections of the rescaled `j1lib.py` block,
the kernel vectors are annihilated, and the dimension counts add up.

### 3.3 Certificates

`cert_f2.py PART U0 U1 [tw=…]`: builds the reduced blocks over Q(a or q, x_i), checks the Gram determinants and the
kernel identities, deletes the gauge coordinates, eliminates the R-states (the R-block is diagonal with positive
entries) and checks every leading principal minor of the Schur complement in the order (P-vectors, then Q-vectors),
with the same sign test as in Section 2.4. Offsets: u = 3 + x for the parts that involve triples inside a level (a
factor such as L − 3 = a + (u₀ − 3) + u₁q is one-signed only for u₀ ≥ 3), u = 4 + x for the part two, fixed sizes
1 and 2 otherwise.

**Theorem F2a.** For every k ≥ 4, every t and all weights with w_c equal on U = [k] ∖ {t} (w_t arbitrary), the
family-2 block is ⪰ 0 with kernel exactly span{p ⊗ |2@t 1@c⟩ : c ∈ U}.

*Proof.* Regimes u = 3 + x (inv, std, anti) and u = 4 + x (two), variables (a, x): all pass, in about a second
each (`run_all.sh` steps `f2_one_*`). The kernel has dimension 1 (inv) + (u − 1) (std) = u, the number of the
vectors p ⊗ |2@t 1@c⟩. ∎

**Theorem F2b.** Let k ≥ 5, t ∈ [k], and let the weights take at most two distinct values, with U = [k] ∖ {t} in one
level or with one of the two U-levels of at most M_2 copies (M_2 as in row 10 of the status table). Then the family-2
block is ⪰ 0 with kernel exactly span{p ⊗ |2@t 1@c⟩ : c ∈ U}. Both U-levels have more than M_2 copies only if
k − 1 ≥ 2M_2 + 2, so this covers every weight vector with at most two distinct weights and k ≤ 2M_2 + 2 (= 26).

*Proof.* If U lies in one level this is Theorem F2a. Otherwise, by (S1), the parts can be taken on level 0 with the
weight of t equal to that of level 0 (tw = 0) or of level 1 (tw = 1), and for inv and mix (which are mapped to
themselves by (S1)) tw = 0 suffices. Certified regimes (all in `run_all.sh`): inv at (1, 3+), (3+, 1), (2, 2), (2, 3+),
(3+, 2), and at (m, 3+) and (3+, m) for 3 ≤ m ≤ M_2; mix at (2, 2), (2, 3+), (3+, 2), (3+, 3+); and for tw = 0 and 1: std at (2, 2), (2, 3+), (3+, 1), (3+, 2),
(3+, 3+), two at (4+, 1), (4+, 2), (4+, 3+), anti at (3+, 1), (3+, 2), (3+, 3+). `coverage.py M_2` checks that every configuration of the statement with 5 ≤ k ≤ 60 is in
the orbit of a certified regime, for every part. Every numerator and denominator is one-signed as it stands. The
largest certificates are the regimes (3+, 3+) of std: for tw = 1 the last minor (5 P-vectors and 2 Q-vectors) has
numerator degree 140 with 63 948 terms over a denominator of degree 132 with 52 384 terms (about ten minutes). ∎

### 3.4 Controls

- `neg_f2.py`: with B → (1 + 1/1000)B the reduced forms of the gauge-fixed block fail at 64 of 90 candidate points
  while the true ones pass; at 6 of them (two levels, including points where all weights take two values, for example
  w = (1/30, 1, 1, 1, 1/30)) the full scaled gauge-fixed block from `j1lib.py` is not positive definite while the
  true one is. The failing parts are inv and std.
- `run_all.sh` step `f2_neg_cert_one`: the one-level inv certificate fails with the scaled B-part.

### 3.5 Not covered

*The inv part with both U-levels large, at two weight values.* In the regime (3+, 3+) (variables q, x₀, x₁) the
gauge-fixed invariant form has 6 P-vectors and 3 Q-vectors after the elimination of the R-states. Its first seven
leading minors pass (minor 7 has numerator degree 166 with 109 975 terms over a denominator of degree 152 with 79 079
terms, after about 19 minutes); during minor 8 the process passed 2.5 GB and was stopped by `guard.sh`
(`attempts.sh`; the lines it displayed are in `logs/attempts/stopped_regimes.txt`). The other order (vacuum states first) was slower in
a development run. So the inv part is certified only with one level of fixed size at most M_2; the configurations
with both U-levels larger than M_2 remain open. A smaller certificate (for example the congruence that removes the
inverse of Â₂ from the R–Q block, or evaluation and interpolation of the last two minors) is the obvious next step.

*Arbitrary t-weight with two U-levels.* With a free (variables a, q, x₀, x₁) the inv part at (3+, 3+) passed the Gram,
kernel, R-block and vacuum-block checks (vacuum minors up to 16 868 terms) and was stopped while eliminating the
vacuum states (`logs/dev/dev_cert_f2_inv_33_freea_stopped.log`). So for family 2 the stated target (arbitrary weight
on t, U-weights on two levels) is **not** reached; Theorems F2a and F2b are what is proved.

## 4. Family 3: the letter-1 block

States (1@r, 1@P) with P a pair, and (vac, 1@T) with T a triple. At k = 5, 6, 7 and two-level rational weights the
nullity is exactly k(k − 1)/2 (the vectors p ⊗ |1@P⟩; no vector of the type of Proposition K2 appears for k ≥ 5) and
the block with one row deleted per column group is positive definite (`explore_kernel.py`, exact). The symbolic
reduction was not attempted: its states carry three indices, so the parts include irreps with up to three rows on
each level and mixed products of them, and the vacuum block needs the reduced inverse of an operator on triples.
This family is **open**.

## 5. What this gives for (R∞)₂

For a weight vector with at most two distinct weights and k ≥ 5: the type I blocks with |T| = 2 (within the range of
Theorem F1b) and with |T| = 1 (within the range of Theorem F2b) are ⪰ 0 with kernel exactly the part of p ⊗ H₂ they
contain; the type I blocks with |T| = 3 are positive definite (Theorem DJ of the README, all weights); the type II
blocks are positive definite (Proposition T2′, all weights). So at two weight values every torus block of (R∞)₂ except
the letter-1 block holds for 5 ≤ k ≤ 82 (k ≤ 26 with `run_all.sh` alone); for larger k every block except the
letter-1 block, the |T| = 2 blocks of row 5b and the |T| = 1 blocks of row 10b. The letter-1 block is not treated
here (see `code/middle_j2_l1/` and `code/middle_j2_l1_sym/`). Beyond two values, Theorems
F1a and F2a allow arbitrary weights on the copies of T when the remaining weights are equal.

## 5b. The fixed-size regimes up to m = 40 (`run_fext.sh`)

`run_fext.sh` runs the same certificates as `run_all.sh` at the fixed-size regimes with 13 ≤ m ≤ 40: for family 2
the invariant part at (m, 3+) and (3+, m), for family 1 the regime tw = (0, 1) at (m, 2+), 84 certificates in all
(each 6 to 30 seconds; the largest sizes do not grow with m, only the coefficients). Then `coverage.py 40 M1=40
KMAX=100` checks that every configuration with 5 ≤ k ≤ 100 inside the statements of Theorems F1b and F2b with
M_1 = M_2 = 40 lies in a certified regime (`coverage.py` takes the largest k as an option; the default 60 is the
one used by `run_all.sh`). All steps pass (`logs/run_fext.out`).

## 6. Files and the final run

| File | Content |
|---|---|
| `r2lib.py` | the states of the three families for `j1lib.py`, the kernel vectors p ⊗ y, exact rank and LDL |
| `explore_kernel.py` | row 1: sizes, kernels and gauge-fixed positivity at k = 5, 6 and k = 7, exact, at points |
| `check_ylib_r2.py` | `j1lib.py` (Lemma B, j = 2) equals `ylib.py` (from the definitions) on the three families, k = 5, 6 |
| `t2lev.py` | family 1: the closed forms of Lemma F1 in the invariant-operator algebra, Φ (Lemma P2) and the gauge-fixed block, leading minors |
| `check_t2.py` | family 1: closed forms entry by entry, Lemma P2 identities, reduced forms against projections (Sections 2.1–2.3) |
| `cert_t2.py` | family 1 certificates (Theorems F1a, F1b) |
| `neg_t2.py` | family 1 negative controls at points |
| `f2sym.py` | family 2: the closed forms of Lemma F2, the parts and test vectors, pattern sums, reduced forms |
| `check_f2.py` | family 2: closed forms entry by entry, reduced forms against projections, kernels, dimensions |
| `cert_f2.py` | family 2 certificates (Theorems F2a, F2b) |
| `neg_f2.py` | family 2 negative controls at points |
| `coverage.py` | coverage of the statements of Theorems F1b (with M_1) and F2b (with M_2) by the certified regimes of `run_all.sh`, 5 ≤ k ≤ 60 |
| `dims.py` | family 2: the test vectors times the irrep dimensions add up to the numbers of states, 5 ≤ k ≤ 60 |
| `check_ylib_reduced.py` | the reduced forms of families 1 and 2 against projections of the block built by `ylib.py` directly (two-level points, k = 5, 6) |
| `j1lib.py`, `ylib.py`, `flintfield.py`, `l1sym.py`, `cert.py`, `guard.sh` | copied unchanged from `relative-anchors/code/` (`middle_j1/`, the code folder, `middle_j1_l1/`, `no_letter1/guard.sh`); `l1sym.py` supplies small helpers, `cert.py` the sign test |
| `run_all.sh` | reruns every check and every certificate of rows 1–13, one Python process at a time under `guard.sh`; PASS/FAIL per step, ALL PASS at the end |
| `attempts.sh` | the two regimes of rows 5b and 10b (not finished; not part of `run_all.sh`) |
| `crosscheck.py`, `run_crosscheck.sh` | the second check of Section 7; output `logs/run_crosscheck.out`, `logs/run_crosscheck.log` |
| `logs/run_all.out`, `logs/run_*.log` | the final run of `run_all.sh` |
| `logs/attempts/` | `stopped_regimes.txt`: the lines displayed by the two unfinished regimes (rows 5b, 10b); `first_run_all.out`, `first_run_extra.out`: the first complete run, in which these regimes were still steps of the runner and the fixed-size regimes were run separately |
| `logs/dev/` | development runs, including the stopped attempts of Sections 2.6 and 3.5 |

The status of rows 5, 10 and 14 rests on the final run: see the last line of `logs/run_all.out`.

## 7. A second check

`crosscheck.py` was written separately from `t2lev.py` and `f2sym.py`; it uses them only to build the symbolic forms
that the certificates use, in the certified regime, and evaluates those at the point. The blocks are built by
`ylib.py` (from the definitions, not from Lemma B), with squared weights 1 and q on the two levels.

- *Family 1 (exact).* With P_t deleted, the Schur complement of the P-block is taken directly (exact rational
  elimination). By the argument of Lemma P2 this is Φ, but it is computed without the formula Φ = C − Bᵀ diag(1/N) B.
  Its projections on the orbit indicators and on e_h − e_{h′} equal the certified invariant and standard forms
  exactly at 7 configurations: one level with free weights on T (k = 7); tw = (0, 0) at (2, 4), (1, 5), (5, 1);
  tw = (0, 1) at (3, 3), (2, 4), (1, 4) (k = 7, 8), with q = 1/900, 1/40, 2/7, 5/3, 30 and 40, and the copies of the two
  levels interleaved.
- *Family 2 (spectral).* The union of the generalized eigenvalues of the certified reduced forms of all parts (each
  with its irrep dimension as multiplicity; a part on level 1 through the level exchange (S1), with the eigenvalues
  multiplied by q²) equals the spectrum of the block, to 10⁻⁹ relative to its largest eigenvalue, at 5
  configurations: one level with free weight on t (k = 7), and (3, 4), (2, 5), (6, 1), (3, 3) with both t-weights,
  q = 1/900, 1/30, 17/5 and 25. This checks at once that the parts are complete and that every reduced form is right.
- *End to end (exact).* At the same 12 blocks the column-group indicators are in the kernel, the nullity equals their
  number, and the block with one P-state per column group deleted is positive definite (exact LDL).
