# One marked copy: the type I blocks with a letter 1

6 November 2026. Nothing here is refereed.

Squared weights s_0 = σ > 0 (copy 0 is *marked*), s_1 = … = s_{k−1} = 1. In (R∞)_j a *type I block* is labelled by a
set T ≠ ∅ of copies: every copy of T carries the letter 2 exactly once (in the row or in the column), no other copy
carries it. It contains a letter 1 iff |T| ≤ j. Its known kernel is the part of p ⊗ H_j in it,
P_T := span{p ⊗ y : y ∈ H_j with the letter 2 on every copy of T and j − |T| letters 1 off T}.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Lemma R′ (reduction).** For each T, the block commutes with S_{T′} × S_{U′} (T′, U′ the copies of T, U other than copy 0). It is ⪰ 0 with kernel exactly P_T iff, for each T′-mode (sym: S_{T′}-invariant; anti: the irrep (τ′ − 1, 1)) and each irrep index i, the reduced form on the twisted vectors (at most 13 + 3 vectors when copy 0 ∈ T, 22 + 6 when copy 0 ∈ U) is ⪰ 0 with kernel exactly the twisted part of P_T | **proved** by hand (Section 3); **exact checks** against the Lemma B block for k ≤ 9 and at k = 10 (row 5) |
| 2 | **Lemma A′ (closed form)** of the block in the rescaled frame with letters 2, and **Lemma E′** (the energy operator on strings with letters 2: one Johnson operator if copy 0 carries a 2, the 2 × 2 form of Lemma E of `marked_copy/` with the constant κ₂ added if not) | **proved** by hand (Section 2); **exact checks** (row 5) |
| 3 | **Theorem I1 (fixed k).** Let 4 ≤ k ≤ 29 and σ > 0. For every j with 3 ≤ j ≤ k − 1 and every T with 1 ≤ \|T\| ≤ j, the type I block of (R∞)_j labelled by T, at squared weights (σ, 1, …, 1), is ⪰ 0 with kernel exactly P_T; hence the same for every weight vector on k copies whose squared weights are all equal but one. This includes (R∞)_{k−1} with the copy of the other weight in U, \|T\| ≥ 2 and \|U\| ≥ 3 | **proved**, computer-assisted, exact: 76 736 reduced forms over Q(σ) for 8 ≤ k ≤ 29 (and 268 for 4 ≤ k ≤ 7), each certified for every σ > 0 by its pivots, all ratios of polynomials in σ with only nonnegative coefficients (Sections 4 and 8) |
| 4 | **Corollary I2.** For 8 ≤ k ≤ 29 and every weight vector on k copies whose squared weights are all equal but one, all of (R∞) holds (every torus block of every (R∞)_j, j = −1, …, k, is ⪰ 0 with exactly its known kernel). Hence Conjecture 12 of `single-excitation-range/` holds for every p ∈ H_1 supported on at most 29 copies whose nonzero weights are all equal but one | **proved**, given row 3 and the earlier theorems listed in Section 5 (Theorem MK for the letter-1 blocks) |
| 5 | Checks: (a) every entry of every reduced form equals the projection of the rescaled Lemma B block (`fastj1.py`, exact) onto the explicit twisted vectors, for k = 3, …, 9 and every j ≥ 1 and at k = 10 for j = 3, 8, 9, every T-size and case, B scaled by 1 and by 1001/1000; the twisted spaces split into exactly the two expected U′-irreps and the T′-irrep of the mode, and the dimension counts add up to the number of states of every block; (b) symbolic forms against numeric ones and the Newton window against all Johnson nodes, every form with k ≤ 10 and every 20th form with 18 ≤ k ≤ 24; every sector denominator has nonnegative coefficients; (c) whole blocks end to end, without the reduction: at k = 8 (σ = 4, 1/9) and k = 9 (σ = 9/4) for every j = 3, …, k − 1, at k = 10 (σ = 1/4; j = 3, 8, 9) and k = 12 (σ = 4; j = 3, 10, 11), every T-size and case: ⪰ 0 with kernel exactly P_T by an exact positive-definiteness certificate | **exact checks**, not proofs (Sections 6 and 8) |
| 6 | Negative controls, B → (1 + 1/1000)B: the certificates fail for many forms at every k = 8, …, 12, with an explicit σ; in every end-to-end block ⟨κ, Ỹκ⟩ < 0 for a known kernel vector κ | **exact** (Section 6) |
| 7 | Coverage: for each k = 4, …, 29 every (j, \|T\|, case, mode, i) with a nonzero twisted space has a certificate with the expected kernel; for k ≤ 8 every torus group of every (R∞)_j is one of the block types used, and the type I blocks with a letter 1 have exactly the states used by the checks | **exact bookkeeping** (`coverage_t.py`, `blocks_t.py`) |
| 8 | Families symbolic in the sizes (every k): begun in `families/` (Theorems AT, S1, SD, S2 there: the anti mode and the sym-mode families with few large sizes, for every k); the rest of k ≥ 30; two or more marked copies; general weights | **open** (Section 7; `families/NOTES.md`) |

Nothing here says that Conjecture 12 of `single-excitation-range/` holds in general, or anything about the NPT (Werner
distillability) problem. Computer-assisted means exact rational arithmetic (python-flint `fmpq`, `fmpq_poly`) and sign
checks of coefficients. No floating point enters a proof (it only proposes the integer matrices of the end-to-end
positive-definiteness test, which is then exact).

Notation is that of `relative-anchors/README.md` (Sections 1, 2, 4, 4b) in the repository `werner-npt-research`.

## 1. Setting

**The two cases.** Let τ = |T| and u = k − τ.
- *Case MT* (copy 0 ∈ T): T′ = T − {0} has τ′ = τ − 1 copies, U′ = U has u′ = u copies.
- *Case MU* (copy 0 ∈ U): T′ = T has τ′ = τ copies, U′ = U − {0} has u′ = u − 1 copies.

All copies of T′ ∪ U′ have squared weight 1, so the group G = S_{T′} × S_{U′} permutes the states of the block and
commutes with it.

**Reduction to two configurations.** 2^Z Y does not change when all weights are multiplied by one positive factor
(in Lemma B, 𝒯_a and J^± are linear and E is quadratic in the weights, and every B-term has the form 𝒯ᵀE⁻¹𝒯,
J⁻E⁻¹𝒯 or J⁻E⁻¹J⁺). A permutation of the copies maps the block labelled by T at weights s to the block labelled by
πT at the permuted weights. So for every weight vector whose squared weights are all equal but one, every type I
block is the block at (σ, 1, …, 1) with T = {0, …, τ − 1} (case MT) or T = {1, …, τ} (case MU). (At σ = 1 both
apply.)

**Known kernel.** By Lemma N, Σ_r w_r 𝒯_{(r,1)} = E_c, which is the kernel p ⊗ H_j of (R∞)_j; the vectors p ⊗ y in the
block are those with y carrying the letter 2 on all of T. In the rescaled frame below, p ⊗ y is (up to a positive
factor) the indicator of its column group {(1@r, y) : r ∈ [k]}.

## 2. The closed form and the energy operator

**Rescaled frame.** Each state (a, y) is multiplied by μ(a, y) = Π_t w_t^{N_t}, N_t the number of letters 1 at copy t
in the row and the column (as in `ylib.py`; letters 2 are not counted), and Ỹ = D (2^Z Y) D with D = diag(μ). With
ψ_a the column vector of the row a and φ that of the vacuum row, write ψ_a = w_a^{[ℓ_a = 1]} D_w ψ̃_a and
φ = D_w φ̃, where D_w|y⟩ = Π_t w_t^{[y_t = 1]}|y⟩.

**Lemma A′.** With Π = D_w² = diag(σ^{[y_0 = 1]}), s_a^{ℓ} := σ if a = (1@0) and 1 otherwise, j̄^± the unweighted
raising and lowering operators, t_{(r,1)} = ½σ_{(r,1)}⁻ j̄⁺ + N_r (N_r = 3/2, 1/2, 1 on the letters 0, 1, 2 at r),
t_{(t,2)} = ½σ_{(t,2)}⁻ j̄⁺ − ½|1⟩⟨2|_t and Ê = D_w⁻¹ E D_w⁻¹:

    2^Z Q(Ψ) = Σ_{a,b} (½ − [r_a = r_b]) s_a^ℓ s_b^ℓ ⟨σ_b⁻ψ̃_a, Π σ_a⁻ψ̃_b⟩ + Σ_a s_a^ℓ ⟨ψ̃_a, Π (2 − n_{r_a}) ψ̃_a⟩
             − Σ_{a,b} s_a^ℓ s_b^ℓ ⟨t_b ψ̃_a, Ê_j⁻¹ t_a ψ̃_b⟩ − 2 Σ_a s_a^ℓ (⟨ψ̃_a, Π σ_a⁻ φ̃⟩ − ⟨j̄⁺ψ̃_a, Ê_{j+1}⁻¹ t_a φ̃⟩)
             + 2 ⟨φ̃, Π φ̃⟩ − ⟨j̄⁺φ̃, Ê_{j+2}⁻¹ j̄⁺ φ̃⟩.

*Proof.* Substitute into the quadratic form of Lemma B (README, Section 2) and use: J⁺D_w = D_w j̄⁺;
σ_a⁻ D_w = w_a^{[ℓ_a = 1]} D_w σ_a⁻ (lowering a letter 2 does not change D_w); and 𝒯_a = w_a^{[ℓ_a = 1]} D_w t_a D_w⁻¹.
The last one is Lemma N: for a = (r, 1), D_w σ⁻ j̄⁺ D_w⁻¹ = Σ_t (w_t/w_r) σ_r⁻|1⟩⟨0|_t, so w_r D_w t_a D_w⁻¹ =
½σ⁻J⁺ + w_r N_r; for a = (t, 2), D_w⁻¹|1⟩⟨2|_t D_w = w_t⁻¹|1⟩⟨2|_t, so D_w t_a D_w⁻¹ = ½σ⁻J⁺ − ½w_t|1⟩⟨2|_t. Then
⟨D_w x, E⁻¹ D_w y⟩ = ⟨x, Ê⁻¹ y⟩ and ⟨D_w x, D_w y⟩ = ⟨x, Π y⟩. ∎

Under w_0 → −w_0 the block is conjugated by (−1)^{N_0}, and D(−w_0) = (−1)^{N_0} D(w_0), so Ỹ is even in w_0: all
entries are rational functions of σ. Lemma A′ is checked exactly against the block of `fastj1.py` (Section 6).

**Lemma E′ (the energy operator with letters 2).** E preserves the positions of the letters 2 and the number of
letters 1. On the strings with a given pattern of letters 2 let n be the number of unmarked copies without a letter 2,
κ₂ the number of unmarked copies with a letter 2, and c′ the number of letters 1 on unmarked copies.
- If copy 0 carries a 2: Ê = E = κ₂ + σ + (3/2)(n − c′) + ½c′ + ½ DU, DU = j̄_U⁻ j̄_U⁺ the Johnson operator of the n
  Boolean unmarked copies; on the irrep (n − i′, i′) of the c′-subsets, DU = (n − c′ − i′)(c′ + 1 − i′).
- If copy 0 is Boolean, with c letters 1 in all, split into A (letter 0 at copy 0) and B (letter 1 there):

      Ê = [[κ₂ + 2σ + 2(n − c) + ½J_U⁺J_U⁻,  ½J_U⁺σ_0⁻], [½σ_0⁺J_U⁻,  (κ₂ + 2(n − c + 1) + ½σ + ½J_U⁺J_U⁻)/σ]].

  On the irrep i′ put μ = (n − c + 1 − i′)(c − i′), ν = (n − c + 2 − i′)(c − 1 − i′), a = κ₂ + 2σ + 2(n − c) + μ/2,
  e = κ₂ + 2(n − c + 1) + σ/2 + ν/2 and Δ = ae − σμ/4. Then Ê⁻¹ = Δ⁻¹[[e, −½σJ_U⁺σ_0⁻], [−½σσ_0⁺J_U⁻, σa]] there;
  if the irrep occurs on one side only, Ê⁻¹ = 1/a on A or σ/e on B.

*Proof.* The letters 2 are inert for J^± and contribute s_t each to 𝒩 = Σ_t s_t N_t and to
C = Σ_t s_t(2|0⟩⟨0|_t + |2⟩⟨2|_t). First case: E = 𝒩 + ½J⁻J⁺ (Lemma N), all letters 1 sit on copies of weight 1, so
D_w = I on these strings and J⁻J⁺ = DU. Second case: E = C + ½J⁺J⁻ (Lemma 1); this is Lemma E of
`marked_copy/NOTES.md` with u replaced by n and the constant κ₂ added to both diagonal blocks (multiply out). The
determinant is
Δ = σ² + σ((5/2)κ₂ + 2B₀ + n − c) + (κ₂ + A₀)(κ₂ + B₀), A₀ = 2(n − c) + μ/2, B₀ = 2(n − c + 1) + ν/2, with all
coefficients ≥ 0 when the irrep occurs on both sides (then n − c ≥ i′ ≥ 0); so Δ, a, e and the scalar of the first
case are > 0 for σ > 0. ∎

## 3. The reduction (Lemma R′)

**Which irreps occur.** A state distinguishes at most one copy of T′ (the copy of a row letter there; all other copies
of T′ carry the letter 2 in the column). So as an S_{T′}-module the state space is a sum of trivial modules and copies
of C^{T′} = (τ′) ⊕ (τ′ − 1, 1). As an S_{U′}-module it is a sum of permutation modules of subsets of U′ (two-row
irreps) and of C^{U′} ⊗ (subsets of U′) (Young modules M^{(a, b, 1)}: at most three rows, λ₃ ≤ 1), as in Section 4
of the README.

**Twisted vectors.** For 0 ≤ i ≤ u′/2 let H_i = (S_2 ≀ S_i) × S_m act on U′ (pairs (U′_1, U′_2), …, rest R,
m = u′ − 2i) with the sign χ_i on each pair swap. By Pieri, among the irreps above, Ind_{H_i}^{S_{u′}} χ_i contains
exactly (u′ − i, i) and (u′ − i − 1, i, 1), each once. On T′: the *sym* mode takes the S_{T′}-invariants (the irrep
(τ′), once); the *anti* mode takes the vectors with the sign under the swap of two copies Ta, Tb of T′ and invariant
under the permutations of the others (Ind of sgn ⊗ 1 from S_2 × S_{τ′−2} contains (τ′ − 1, 1) once and otherwise only
(τ′ − 2, 1, 1), which does not occur). So every irrep of G that occurs has a one-dimensional space of twisted vectors
for exactly one pair (mode, i), and the reduced form Q_{mode,i} on the twisted vectors V_{mode,i} is the orthogonal sum
of the forms of the two irreps on their multiplicity spaces (up to positive factors). Since P_T is G-invariant, the
block is ⪰ 0 with kernel exactly P_T iff every Q_{mode,i} is ⪰ 0 with kernel exactly P_T ∩ V_{mode,i}. (As in
Section 2 of `marked_copy/NOTES.md`.)

V_{mode,i} is spanned by the twisted orbit sums that are not zero. With q0 = j − τ − i the Dicke count of the base
column on R, singlets |10⟩ − |01⟩ on the pairs and Dicke vectors on R, they are (x, y, z, y0 ∈ {0, 1}; a row letter 1 at
a copy of R or of a pair is followed by the column letters x there and y at the partner, with the sign −1 when the copy
is the second of its pair):
- *Case MT, sym* (13): A_x (row 1@R), B_xy (row 1@pair), C (row 1@T′, column 2 there), D_z (row 2@T′, column z there),
  CM (row 1@0, column 2 at 0), DM_z (row 2@0, column z at 0), G (vacuum row). Every other copy of T carries the letter
  2; the Dicke counts are q0 − x, q0 + 1 − x − y, q0, q0 + 1 − z, q0, q0 + 1 − z, q0 + 1.
- *Case MT, anti* (3, needs τ′ ≥ 2): Ca and Da_z: the rows of C and D_z at Ta minus those at Tb.
- *Case MU, sym* (22): the same rows with the column letter y0 at copy 0 added: O_{y0} (row 1@0), A_{x y0}, B_{xy y0},
  C_{y0} (row 1@T′), D_{z y0} (row 2@T′), G_{y0}; Dicke counts lowered by y0.
- *Case MU, anti* (6, needs τ′ ≥ 2): Ca_{y0}, Da_{z y0}.

Distinct names have disjoint supports, so the nonzero ones form an orthogonal basis of V_{mode,i}. The known kernel
there: in case MT, sym, κ = A_0 + A_1 + B_10 − B_01 + C + CM; in case MU, sym, κ_{y0} = O_{y0} + A_{0y0} + A_{1y0} +
B_{10y0} − B_{01y0} + C_{y0} (y0 = 0, 1; those that are nonzero); nothing in the anti modes (P_T is S_{T′}-invariant).

**The orbit calculus** (`tlib.py`) is that of `marked_copy/mlib.py` and of the E3 calculus `orbit2.py`, joined:
a frame holds copy 0, at most two copies of U′ with their partners and at most two copies of T′; a frame vector is a
dict (letters on the frame, Dicke offset) → coefficient, with singlets on the other pairs, the letter 2 on the other
copies of T′ and a Dicke vector on the rest of R. The entries of Q_{mode,i} are sums over the H-orbits of row sites
and of ordered pairs of row sites (sizes such as m, 2i, τ′, 2im, τ′(τ′ − 1)), each computed in one frame by Lemma A′.
Ê⁻¹ is applied sector by sector (pattern of letters 2) by Lemma E′, each function of the irrep being interpolated
(Newton) in DU over the genuine Johnson nodes i′ ∈ [i − 2, i + 2]. Every application asserts that the remainder after
these nodes vanishes, so the interpolation is exact on the vector (a frame with p distinguished pairs, r distinguished
copies of R and b Boolean copies of T′ has components only on i′ ∈ [i − p, i + p + r + b], and p + r + b ≤ 2 in
every frame used). In the symbolic mode σ is the variable of Q[σ]; each sector returns its result times D, the product
of the denominators (Δ, a, e or the scalar of Lemma E′) at the nodes used, and the entries are assembled as reduced
rational functions.

## 4. Certificates for fixed k (`cert_t.py`)

For each k and each (j, τ, case, mode, i) with 3 ≤ j ≤ k − 1, 1 ≤ τ ≤ j and V_{mode,i} ≠ 0:
1. Q = Q_{mode,i} over Q(σ); Q is symmetric and every entry denominator has no root in (0, ∞).
2. Q κ = 0 identically in σ for every known kernel vector κ.
3. Gauge: delete one coordinate per κ (CM in case MT; O_0, O_1, or C_0, C_1, … in case MU), with a determinant of the
   kernel vectors on these coordinates that is a nonzero constant (it is ±1). Then for every σ > 0 the coordinate
   subspace W on the other coordinates is a complement of span(κ), and Q ⪰ 0 with kernel exactly span(κ) iff Q|_W is
   positive definite.
4. Symmetric elimination of Q|_W in the order of the names; every pivot is a reduced rational function of σ that is
   > 0 on (0, ∞). In every form every pivot has numerator and denominator with only nonnegative coefficients (the
   witness search and the Sturm fallback of `signtest.py` were never needed).

*Why this is a proof for every σ > 0.* The entries of Q are polynomials in σ divided by products of the sector
denominators D, which have nonnegative coefficients (Lemma E′; also checked for every sector met for k ≤ 10 and in
the sampled forms with 18 ≤ k ≤ 24), so the
rational functions are the true reduced forms at every σ > 0. The r-th leading principal minor of Q|_W equals, in
Q(σ), the product of the first r pivots, each positive on (0, ∞), so Sylvester's criterion applies at every σ > 0.
With Lemma R′ and the reduction of Section 1 this proves Theorem I1.

### 4a. Theorem I1 and Corollary I2

**Theorem I1.** For 4 ≤ k ≤ 29, every σ > 0, every j with 3 ≤ j ≤ k − 1 and every T with 1 ≤ |T| ≤ j, the type I
block of (R∞)_j labelled by T, at squared weights (σ, 1, …, 1), is ⪰ 0 with kernel exactly P_T. By Section 1 the same
holds for every weight vector on k copies whose squared weights are all equal but one. (For j = k − 1 this includes
the blocks with the copy of the other weight in U, |T| ≥ 2 and |U| ≥ 3, which were open for k ≥ 8.)

*Proof.* Section 1 reduces T to the two configurations, Lemma R′ to the reduced forms, and the certificates of Section 4
(steps `cert_*` of `run_all.sh` for 8 ≤ k ≤ 20, of `run_extra.sh` for 4 ≤ k ≤ 7 and 21 ≤ k ≤ 28, of
`run_extra2.sh` for k = 29) cover every one of them; `coverage_t.py` checks that no (j, τ, case, mode, i) with a
nonzero twisted space is missing for any k = 4, …, 29. ∎

Sizes: 76 736 reduced forms for 8 ≤ k ≤ 29 (k = 8: 216 forms in 50 blocks (j, |T|, case); k = 20: 3 320 forms in
374 blocks; k = 28: 8 736 forms in 750 blocks; k = 29: 9 665 forms), at most 20 vectors after the gauge; the largest
pivot degree in σ is 80 at k = 8 and stays at 99 for every k from 14 on, because only five Johnson nodes ever enter.
The certificates take about 12 minutes for 8 ≤ k ≤ 20, then from 3 minutes (k = 21) to 10 minutes (k = 29) per k.

**Corollary I2.** For 8 ≤ k ≤ 29 and every weight vector on k copies whose squared weights are all equal but one, all
of (R∞) holds: every torus block of every (R∞)_j, j = −1, …, k, is ⪰ 0 with exactly its known kernel. Hence
Conjecture 12 of `single-excitation-range/` holds for every p ∈ H_1 supported on at most 29 copies whose nonzero weights
are all equal but one.

*Proof.* By Section 1 it suffices to take the squared weights (σ, 1, …, 1).
- (R∞)_{−1}, (R∞)_0 and (R∞)_k: Proposition M, Theorems B0 and K (every weight vector).
- (R∞)_1: the Corollary of Section 3h of the README (Theorems 2V, EQ, 2L, D₁, Proposition T2′: all of (R∞)_1 at
  weights with at most two values).
- (R∞)_2: Corollary M2 (k ≥ 5, weights all equal but one).
- (R∞)_j, 3 ≤ j ≤ k − 1: the letter-1 block by Theorem MK (k ≤ 40); the blocks without a letter 1 by Theorem DJ; the
  type II blocks by Proposition T2′; the type I blocks with a letter 1 by Theorem I1.
The kernels are exactly the known ones block by block, and they add up to those of (R∞)_j. For exactly k nonzero
weights, Conjecture 12 follows with Lemma 13 of `single-excitation-range/` and the remark that (R2) always holds, as
in Corollaries 4 and 7₂; for at most seven nonzero weights it is Corollary 7₂ (two weight values include all equal but
one). ∎

## 5. Sources used by Corollary I2

All in `relative-anchors/README.md`: Proposition M, Theorem B0, Theorem K (every weight vector); Theorems 2V, EQ, 2L,
D₁ and Proposition T2′ (all of (R∞)_1 at two weight values, Section 3h); Corollary M2 (Section 4b); Theorem MK
(Section 4b, letter-1 blocks for k ≤ 40); Theorem DJ (Section 3h); Proposition T2′; Lemma 13 and the remark on (R2)
of `single-excitation-range/`; Corollary 7₂ (Section 3f). The type I blocks of (R∞)_{k−1} covered earlier (Corollary
EM, Theorems EU and TL, U2, Propositions U0 and U1) are certified again here as part of Theorem I1.

## 6. Checks and controls

All exact unless said otherwise; `run_all.sh` reruns them (`logs/run_all.out`).
- `check_explicit.py 3 8 1 2,1/3` (`logs/run_explicit_3_8.log`): for k = 3, …, 8, every j = 1, …, k − 1, every τ and
  both cases, at σ = 4 and 1/9 (w_0 = 2, 1/3), with B scaled by 1 and by 1001/1000: every entry of every reduced form
  of `tlib.py` equals v_xᵀỸv_y, with Ỹ the rescaled block of `fastj1.FastJ1` (Lemma B, every E⁻¹ inverted over Q
  sector by sector) and v_x the explicit twisted vectors of `explicit_t.py` (built from the states, not from the
  frames); the norms agree and the sets of nonzero twisted vectors agree (166 blocks, 608 forms, 150 100 entries).
  Completeness: the class sum of the transpositions of U′ maps every V_{mode,i} into itself with eigenvalues only the
  content sums of (u′ − i, i) and (u′ − i − 1, i, 1); the class sum of T′ acts by the content sum of (τ′) (sym) or
  (τ′ − 1, 1) (anti); and Σ_{mode,i} Σ_λ dim(T′-irrep) · dim λ · (multiplicity, by exact ranks) equals the number of
  states of every block. The known vectors p ⊗ y are annihilated by the block.
- `check_symnum.py 3 10` (`logs/run_symnum_3_10.log`): for every form with k ≤ 10 and j ≥ 1, the symbolic form at
  σ ∈ {3/7, 5/2} equals the numeric one exactly, the numeric form with the window equals the one with all genuine
  Johnson nodes, and every sector denominator D met has nonnegative coefficients (1 350 forms, 200 204 entries, 6 206
  denominators).
- `e2e_t.py` (`logs/run_e2e_8.log`, `logs/run_e2e_9.log`): no reduction. For k = 8 at σ = 4 and 1/9 and k = 9 at
  σ = 9/4, every j = 3, …, k − 1, every τ = 1, …, j and both cases, the whole rescaled block from `fastj1.py`: the
  vectors p ⊗ y are annihilated exactly, and the block without a set of pivot coordinates of them (exact RREF) is
  positive definite by the exact test of `crosscheck.py` (integer congruence Z A Zᵀ, strict diagonal dominance in exact
  integer arithmetic). So the block is ⪰ 0 with kernel exactly P_T there (100 block checks at k = 8, 66 at k = 9;
  largest block 385 states at k = 8). `run_extra.sh` adds k = 10 at σ = 1/4 for j = 3, 8, 9 and k = 12 at σ = 4 for
  j = 3, 10, 11 (Section 8).
- Negative controls, B → (1 + 1/1000)B. `cert_t.py 8 12 1001/1000` (`logs/run_neg_cert_8_12.log`, expected FAIL): the
  kernel identities are not required, the same gauge is used, and at every k = 8, …, 12 many forms fail, with an
  explicit witness σ (10⁻⁶ in the forms listed) where a pivot of the gauge-fixed form is ≤ 0. In `e2e_t.py` the scaled
  block has ⟨κ, Ỹκ⟩ < 0 exactly for a known kernel vector κ of every block, so it is not positive semidefinite (as it
  must: Yκ = 0 gives ⟨κ, (A − 1.001B)κ⟩ = −0.001⟨κ, Bκ⟩).
- `coverage_t.py 8 20` (`logs/run_coverage_8_20.log`): the bookkeeping of Theorem I1 for 8 ≤ k ≤ 20; the later
  scripts redo it for 4 ≤ k ≤ 7, 8 ≤ k ≤ 28 and 8 ≤ k ≤ 29, together with `blocks_t.py` (Section 8).

## 7. What is not done

- **Families symbolic in the sizes** (step 2 of the plan): begun in `families/` (see `families/NOTES.md`). With copy 0 in T the U-part is the
  equal-weight calculus of Theorem E3 (`orbit2.SymF2`) and σ enters only through copy 0, so the families
  [q0, t, i, τ′ fixed or symbolic] with σ free are the natural next step; the generic family has four symbolic sizes
  and σ, and at equal weights already needed Proposition C (a determinant certificate). Case MU adds the letter y0 at
  copy 0 and the 2 × 2 energy of Lemma E′. For the case named in the plan, (R∞)_{k−1} with copy 0 in U, the sizes are
  rigid: t = m − q0 = −i, so only i ∈ {0, 1, 2} occur, q0 = u′ − i, and τ is free; an every-k proof there needs, for each
  i ≤ 2 and both modes, families in (q0, τ, σ) (with q0 or τ also fixed at small values), three variables. In case MT
  at j = k − 1, t = 1 − i and i ≤ 3 (that case is Corollary EM for every k). Here these blocks are proved for k ≤ 29
  only (Theorem I1).
- **k beyond the extension of Section 8.** The cost per k grows roughly like k³ (the number of forms) times a slowly
  growing cost per form; the letter-1 blocks used by Corollary I2 are proved for k ≤ 40 (Theorem MK) and k ≤ 55
  (Theorem M″ of `../middle/NOTES_MI.md`), so the fixed-k route cannot give more than k ≤ 55 without the
  letter-1 families.
- Two or more marked copies, two weight values with both levels large, general weights.

## 8. Larger k

Three more scripts were run after `run_all.sh`, in this order and one Python process at a time:
- `run_extra.sh` (`logs/run_extra.out`): the certificates of Section 4 for each k = 21, …, 28 (`logs/run_cert_<k>.log`,
  `logs/cert_<k>.pkl`) and for 4 ≤ k ≤ 7, with their coverage (`coverage_t.py 4 7`, `coverage_t.py 8 28`); the block
  bookkeeping `blocks_t.py 8 8 28` (below); `check_symnum.py 18 24 20` (every 20th form for 18 ≤ k ≤ 24: symbolic
  equals numeric at two values of σ, the window equals all nodes, every sector denominator has nonnegative
  coefficients; 1 379 forms, 311 520 entries, 6 615 denominators); and whole blocks end to end at k = 10 (σ = 1/4,
  j = 3, 8, 9; 40 blocks) and k = 12 (σ = 4, j = 3, 10, 11; 48 blocks, the largest with 1 045 states), each ⪰ 0 with
  kernel exactly P_T by the exact test, and the scaled control negative;
- `run_check2.sh` (`logs/run_check2.out`): `check_explicit.py` at k = 9 (every j ≥ 1, σ = 9/4, B scaled by 1 and by
  1001/1000) and at k = 10 (j = 3, 8, 9, σ = 1/4): the reduced forms against the Lemma B block and the dimension
  counts, now with every orbit type present at once (m ≥ 2, i ≥ 2 and τ′ ≥ 2 together need k ≥ 9); k = 9: 72 blocks,
  380 forms, 51 208 entries; k = 10: 40 blocks, 221 forms, 21 388 entries;
- `run_extra2.sh` (`logs/run_extra2.out`): the certificates for k = 29 and the coverage and block bookkeeping for
  8 ≤ k ≤ 29.

`blocks_t.py` is the bookkeeping of every block type: (a) for k ≤ 8 and every j, the torus groups of all of (R∞)_j
(every row and column, any letters, grouped by the number of letters 2 at each copy) are exactly the letter-1 block,
the type II blocks and the type I blocks, and every type I block with a letter 1 has exactly the states used by
`explicit_t.py` and `e2e_t.py`; (b) for every k in the range, every (j, |T|, case) with 3 ≤ j ≤ k − 1 and every form of
it is certified, and the result used for every other class of block is listed.

Results: every step of the three scripts passed (each ends with ALL PASS). The certificates for 21 ≤ k ≤ 29 have the
same shape as before (every pivot of every form with nonnegative coefficients, largest pivot degree 99): 3 821, 4 367,
4 963, 5 608, 6 307, 7 059, 7 869, 8 736 and 9 665 forms for k = 21, …, 29. With `run_all.sh` this gives Theorem I1
and Corollary I2 for every k ≤ 29 as stated in Section 4a.

Each further k costs about 1.15 times the previous one (k = 29: about 10 minutes), and the letter-1 blocks
that Corollary I2 needs are proved only for k ≤ 40 (Theorem MK; k ≤ 55 with Theorem M″ of
`../middle/NOTES_MI.md`), so this route can reach at most k ≤ 55 and not every k.

## Files

| File | Content |
|---|---|
| `tlib.py` | the orbit calculus (Sections 2, 3): twisted vectors of both cases and modes, orbits, Lemma E′ with Newton interpolation (window asserted), numeric and symbolic in σ, the known kernel vectors, the list of forms |
| `explicit_t.py` | explicit objects, independent of the frames: the states of a type I block, the rescaled Lemma B block from `fastj1.py`, the explicit twisted vectors, the class sums, irrep dimensions and content sums, the known kernel vectors |
| `cert_t.py` | the fixed-k certificates (Section 4); with B scaled, the negative control |
| `check_explicit.py`, `check_symnum.py`, `e2e_t.py`, `coverage_t.py`, `blocks_t.py` | the checks and the bookkeeping of Sections 6 and 8 |
| `fastj1.py`, `j1lib.py`, `crosscheck.py`, `ufield.py`, `signtest.py`, `guard.sh` | copied unchanged (FastJ1 builder, Lemma B operators, the exact PD test, rational functions in σ, sign test, the memory and time guard) |
| `run_all.sh` | the checks of Section 6 and the certificates for 8 ≤ k ≤ 20 with their coverage, one Python process at a time under `guard.sh`: `bash run_all.sh > logs/run_all.out` (12 steps, 16 minutes, ALL PASS) |
| `run_extra.sh`, `run_check2.sh`, `run_extra2.sh` | run after `run_all.sh`, in this order (Section 8): the extension to 21 ≤ k ≤ 28 and 4 ≤ k ≤ 7 with more checks; the explicit check at k = 9, 10; k = 29 and the final coverage. Outputs `logs/run_extra.out`, `logs/run_check2.out`, `logs/run_extra2.out` |
| `logs/` | the logs of the final runs; `logs/cert_*.pkl` hold the certified gauge-fixed forms (names, kernel vectors, gauge, pivot signs) |
