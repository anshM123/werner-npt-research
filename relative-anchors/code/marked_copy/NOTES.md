# The letter-1 block of every (R∞)_j with one marked copy

5 November 2026. Nothing here is refereed.

In `relative-anchors/README.md` (Section 4b) Theorems M, B and I of these notes are called Theorems MK, MB and MI.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Theorem M (fixed k).** Let 2 ≤ k ≤ 32 and σ > 0, and let the squared weights be s_0 = σ, s_1 = … = s_{k−1} = 1. For every j = −1, …, k the letter-1 block of (R∞)_j is ⪰ 0, and its kernel is exactly K_j: p ⊗ (letter-1 part of H_j), plus span(vac ⊗ J⁻\|1…1⟩) when j = k − 2, plus {u ⊗ J⁻\|1…1⟩ : u ∈ span(e_r), u ⊥ p} when j = k − 1 | **proved**, computer-assisted, exact: 3 909 reduced forms, each certified for every σ > 0 by its leading principal minors, all ratios of polynomials in σ with only nonnegative coefficients (Sections 2–4) |
| 2 | **Corollary M.** The same for every weight vector on k ≤ 32 copies in which all squared weights but at most one are equal | **proved**, given row 1 (permutations and homogeneity, Section 1) |
| 3 | Lemma R (reduction): the S_{k−1}-twisted invariants (at most 16 per i), the closed form of the block in the rescaled frame, and the 2 × 2 inverse of the energy operator on each isotypic part | **proved** by hand (Section 2); **exact checks** against the Lemma B builder (row 6) |
| 4 | **Theorem B (every k, both ends).** For every k ≥ 2, every σ > 0 and every j with j ≤ 8 or j ≥ k − 8, the letter-1 block of (R∞)_j with squared weights (σ, 1, …, 1) is ⪰ 0 with kernel exactly K_j (and, by Corollary M's argument, the same for every weight vector with all squared weights but one equal) | **proved**, computer-assisted, exact: families symbolic in k and σ (110 band families with i fixed and one symbolic size: 109 proved, 1 empty), plus Theorem M for the members with k ≤ 32 (Section 5) |
| 5 | **Theorem F.** Every family with one symbolic size (q0 = j − i ≥ 4, t = k − 1 − j − i ≥ 4 or i ≥ 2, the other two fixed small): the reduced form is ⪰ 0 with kernel exactly the twisted part of K_j, for every member and every σ > 0 | **proved**, computer-assisted, exact: 45 nonempty families of 52; pivots ratios of polynomials in (x, y, z, σ) with coefficients of one sign after the shifts (Section 5) |
| 5a | **Theorem I.** For every k ≥ 2, every j = −1, …, k and every σ > 0, the form on the vectors invariant under the permutations of the unmarked copies (i = 0) is ⪰ 0 with kernel exactly the invariant part of K_j | **proved**, computer-assisted, exact: the family [q0 ≥ 4, t ≥ 4, i = 0], the families of Theorem F and Theorem M (Section 5) |
| 6 | The reduced forms equal the projections of the block built by `j1lib.py` (Lemma B) onto the explicit twisted vectors | **exact check**: k = 2, …, 6, every j and i, 3 to 6 values of σ, B scaled by 1 and by 1001/1000 (all entries); the saved certified forms at k ≤ 6 against the same projections; the symbolic forms against the numeric ones for k ≤ 16; the Newton window against all Johnson nodes for k ≤ 16 |
| 7 | Completeness of the reduction: the class sum of the transpositions splits each twisted space into the two irreps, and Σ_i Σ_λ dim(λ)·mult = number of states | **exact check**, k ≤ 6, two values of σ, every j; the spectrum of the block equals the union of the spectra of the parts (floating point, not used in a proof) |
| 8 | End to end at points: the known kernel vectors are annihilated by the Lemma B block, its exact nullity equals their number, and it is ⪰ 0 | **exact check**, k ≤ 6, every j, three values of σ; not a proof |
| 9 | Negative controls, B → (1 + 1/1000)B: for every k = 2, …, 16 and every j the certificate of the reduced form at i = 0 fails, with an explicit σ (10⁻⁶; for j = −1 and k ≥ 3, 10⁴ or 10⁶) where a pivot is ≤ 0; at these points (k ≤ 6) the restriction of the scaled Lemma B block to the kept twisted vectors is not positive definite while the true one is; the family certificate of Theorem F no longer passes for 42 of its 45 nonempty families | **exact** |
| 10 | k ≥ 33 and 9 ≤ j ≤ k − 9, on the parts i ≥ 1 (needs the families [q0 ≥ 4, t ≥ 4, i = b] for every b ≥ 1, i.e. [q0 ≥ 4, t ≥ 4, i = 1] and the generic [q0 ≥ 4, t ≥ 4, i ≥ 2]; development runs stopped, Section 6); two or more marked copies; general weights | **open** |
| 11 | A second check written separately (Section 7). (a) A separate builder of the block from Lemma B (a Gram form of E⁻¹ on the vectors T_r\|y⟩ and J⁺\|y⟩, flint matrices) agrees with `j1lib.py` entry by entry for k = 4, 5, 6; with it the block is ⪰ 0 with kernel exactly K_j (known vectors annihilated exactly, their number C(k, j) + [j = k − 2] + (k − 1)[j = k − 1], the block without pivot coordinates of them positive definite by an exact test) for k = 7, 8, 9, every j, at 5, 5 and 3 values of σ. (b) The reduced forms of `mlib.py` equal the exact projections of this block onto the explicit twisted vectors for k = 7, 8, 9, every j and i (20 392 entries). (c) The forms of all 155 nonempty families used by Theorems B, F and I equal N0 times the numeric forms at members with symbolic sizes up to 40 (k up to about 66; 1 192 points) | **exact check** (`markcheck.py`, `projcheck.py`, `famcheck.py`, `run_check.sh`); not a proof |
| 12 | Continued in `middle/`: the irrep split of each reduced form, the families [q0 ≥ 4, t ≥ 4, i = b] for b = 1, …, 4 (Theorems N1 to N4 there), hence every k and every j on the parts i ≤ 4 (Corollary P4 there; Theorem MS of the README), and every j for 33 ≤ k ≤ 40 (Theorem M′ there). Row 10 is then open only for k ≥ 41, 9 ≤ j ≤ k − 9, i ≥ 5 | **proved**, computer-assisted (see `middle/NOTES.md`); the second check (Section 7, `run_check.sh`) also compares those families with `mlib.py` at far members and checks 40 of their members with 41 ≤ k ≤ 60 exactly |
| 13 | Continued in `middle/NOTES_MI.md`: the three-row parts of the generic family with i symbolic (Theorem G3 there), the eleven edge families with i symbolic (Theorem E), the families i = 5, 6 (Theorems N5, N6), hence every part except the two-row parts with i ≥ 7, q0 ≥ 4, t ≥ 4 for every k (Theorem MF′ there; Theorem MP of the README), and every j for 41 ≤ k ≤ 55 (Theorem M″). Row 10 is then open only for k ≥ 56 on those two-row parts | **proved**, computer-assisted (see `middle/NOTES_MI.md`); the second check `run_check_mi.sh` (`mfcheck.py`) checks 240 members with 56 ≤ k ≤ 100 exactly |
| 14 | Continued in `typeI/`: the type I blocks with a letter 1 of (R∞)_j, 3 ≤ j ≤ k − 1, for 4 ≤ k ≤ 29 (Theorem I1 there), hence all of (R∞) and Conjecture 12 for supports of at most 29 copies with all nonzero weights but one equal (Corollary I2) | **proved**, computer-assisted (see `typeI/NOTES.md`); second check `run_markall.sh` (`markall.py`: every torus group from Lemma B at k = 8, 9, 10) |
| 14b | (R∞)_{k−1} completely for every k at weights all equal but one (Corollary M1 of `../typeI_G_tv/`: Theorem M1 for the type I blocks with the copy of the other weight free of the letter 2, Corollary EM when it carries one); for j = k − 1 this extends Theorem I1 of `typeI/` from k ≤ 29 to every k | **proved**, computer-assisted (see `../typeI_G_tv/NOTES.md`); second check `run_markall2.sh` |
| 15 | Continued in `middle/NOTES_G2.md`: the two-row parts of the generic family [q0 ≥ 4, t ≥ 4, i ≥ 7] (Theorem G2 there, by determinant and inertia with exact interpolation), hence every letter-1 block of every (R∞)_j for every k at weights all equal but one (Theorem MF) | **proved**, computer-assisted (see `middle/NOTES_G2.md`) |
| 16 | Continued in `typeI/families/`: the type I blocks with a letter 1 with the sizes symbolic, for every k: the anti mode (Theorem AT), the three-row isotypic parts (Theorem Q3) and the sym-mode families with few large sizes (Theorems S1, SD, S2, S3, SD2, MD); about 58 % of the reduced forms with 30 ≤ k ≤ 120 are covered | **proved**, computer-assisted, for those families; the rest **open** |

Nothing here says that Conjecture 12 of `single-excitation-range/` holds in general, or anything about the NPT (Werner
distillability) problem. Computer-assisted means exact rational arithmetic (python-flint `fmpq`, `fmpq_poly`,
`fmpq_mpoly`) and sign checks of coefficients. No floating point enters a proof.

## 1. Setting and statement

Notation is that of `relative-anchors/README.md`, Sections 1 and 2, in the repository `werner-npt-research`. There are
k ≥ 2 copies; copy 0 is *marked*, with squared weight s_0 = σ > 0, and the u = k − 1 copies of U = {1, …, u} have
squared weight 1. The letter-1 block of (R∞)_j has the states (1@r, 1@Y) with r ∈ [k] and |Y| = j, and (vac, 1@Z) with
|Z| = j + 1. We use Lemma B of the README in units 2^Z (Z = k − j) and the rescaled frame of `ylib.py`: each state is
divided by Π_t w_t^{N_t}, N_t the number of letters 1 at copy t in its row and column. Positive diagonal rescaling keeps
the inertia and maps kernels to kernels.

K_j denotes the space of the known kernel vectors (Lemma N, Proposition K2 and Theorem T of the README):
- p ⊗ |1@Y⟩ for every j-subset Y (in the rescaled frame: the indicator of the column group {(1@r, Y) : r ∈ [k]}, up to a
  positive factor);
- at j = k − 2 also vac ⊗ J⁻|1…1⟩ (rescaled: Σ_h s_h (vac, 1…1 − e_h));
- at j = k − 1 also u ⊗ J⁻|1…1⟩ for u ∈ span(e_r) with u ⊥ p (rescaled: Σ_{r,h} d_r s_h (1@r, 1…1 − e_h) with
  Σ_r d_r s_r = 0).

**Theorem M.** For 2 ≤ k ≤ 32, every j = −1, …, k and every σ > 0, the letter-1 block of (R∞)_j is ⪰ 0 with kernel
exactly K_j.

**Corollary M.** For k ≤ 32 the same holds for every weight vector whose squared weights are all equal except possibly
one. *Proof.* A permutation of the copies maps the block onto the block at the permuted weights, and in the rescaled
frame the block is homogeneous of degree j + 1 in the squared weights (Lemma R of `code/middle_k56/NOTES.md`, with
T = ∅). Dividing all squared weights by the common value of the unmarked ones and moving the marked copy to position 0
gives the weights (σ, 1, …, 1). ∎

At σ = 1 this is Theorem E1 of the README for k ≤ 32 (proved there for every k); for j ∈ {−1, 0, k − 1, k} it also
follows from Proposition M, Theorem B0 (k ≥ 3), Theorem T and Theorem K (all weights), for j = 1 from Theorem 2V, for
j = 2 and 5 ≤ k ≤ 10 from Theorem L (two weight values), and for 4 ≤ k ≤ 6 from Theorems 4, 5₂ and 6₂. The proof here
is independent of these and uniform in j.

## 2. The reduction (Lemma R)

**Closed form.** With π_y = σ^{y_0} (y_0 the column letter at the marked copy), s_a = σ if the row letter sits at the
marked copy and 1 otherwise, j̄^± the unweighted raising and lowering operators, t_a = ½σ_a⁻j̄⁺ + N_a
(N_a = (3/2)|0⟩⟨0|_a + ½|1⟩⟨1|_a) and Ê_c = D_w⁻¹E_cD_w⁻¹ (D_w|y⟩ = w_0^{y_0}|y⟩), the entries of the rescaled block are
- Ỹ[(a, y′), (b, y)] = (½ − [a = b]) s_a s_b π_{y−a} [y′ = y − a + b] + [a = b][y′ = y] s_a π_y (2 − n_a)
  − s_a s_b ⟨y′| t_bᵀ Ê_j⁻¹ t_a |y⟩;
- Ỹ[(a, y′), (vac, y)] = −[y′ = y − a] π_y + s_a ⟨y′| j̄⁻ Ê_{j+1}⁻¹ t_a |y⟩;
- Ỹ[(vac, y′), (vac, y)] = 2[y′ = y] π_y − ⟨y′| j̄⁻ Ê_{j+2}⁻¹ j̄⁺ |y⟩.

*Proof.* The proof of Lemma A of `code/middle_j2_l1/NOTES.md` (written there for j = 2, from Lemmas B and N, for every
positive weight vector) applies verbatim to every j: D_w⁻¹J^±D_w = j̄^±, 𝒯_a = w_a D_w t_a D_w⁻¹ (t_a in the normal order
of Lemma N) and E = D_w Ê D_w; the products of weights are μ(a, y′)μ(b, y) = s_b π_y = s_a s_b π_{y−a} when
y′ = y − a + b. The formulas are also checked exactly against `j1lib.py` (Section 4). ∎

**The energy operator (Lemma E).** On the strings with c letters 1 write A for those with letter 0 at the marked copy
(c-subsets of U) and B for those with letter 1 there ((c − 1)-subsets of U). Then

    Ê_c = [[2σ + 2(u − c) + ½J_U⁺J_U⁻,  ½J_U⁺σ_0⁻], [½σ_0⁺J_U⁻,  (2(u − c + 1) + ½σ + ½J_U⁺J_U⁻)/σ]].

On the isotypic part of the S_u-irrep (u − i′, i′) put μ = (u − c + 1 − i′)(c − i′) (the eigenvalue of J_U⁺J_U⁻ on the
A-part and of J_U⁻J_U⁺ on the B-part), ν = (u − c + 2 − i′)(c − 1 − i′) (that of J_U⁺J_U⁻ on the B-part),
a = 2σ + 2(u − c) + μ/2, e = 2(u − c + 1) + σ/2 + ν/2 and Δ = ae − σμ/4. Then

    Ê_c⁻¹ = Δ⁻¹ [[e,  −½σJ_U⁺σ_0⁻], [−½σσ_0⁺J_U⁻,  σa]]   on that part.

If the irrep occurs on only one side, μ = 0 there and the inverse is 1/a on A, σ/e on B.

*Proof.* C = Σ_t s_t·2|0⟩⟨0|_t and ½J⁺J⁻ = ½(J_U⁺ + w_0σ_0⁺)(J_U⁻ + w_0σ_0⁻) give the blocks of E_c; conjugation by D_w
divides the off-diagonal blocks by w_0 and the B-block by σ. On the irrep, the matrix [[a, ½X], [½X*, e/σ]] with
X = J_U⁺σ_0⁻, XX* = μ on A and X*X = μ on B has the stated inverse (multiply out). With A₀ = 2(u − c) + μ/2 and
B₀ = 2(u − c + 1) + ν/2, Δ = σ² + (2B₀ + u − c)σ + A₀B₀; when the irrep occurs on both sides, c ≤ u and μ, ν ≥ 0, so
every coefficient is ≥ 0 and Δ > 0 for σ > 0. ∎

**Twisted invariants.** Let S_u act on U, let H_i = (S₂ ≀ S_i) × S_m (pairs (1, 2), (3, 4), …, and the rest R,
m = u − 2i), and χ_i the sign on each pair swap. Put q0 = j − i and t = u − j − i, so m = q0 + t.
- Every S_u-irrep in the state space is (u − i, i) or (u − i − 1, i, 1) for some i: the states with the row letter in U
  span C^U ⊗ (subsets of U), whose orbits are Young modules M^{(a, b, 1)} (dominance gives at most three rows and
  λ₃ ≤ 1); the others are permutation modules of subsets (two rows).
- Pieri (as in Section 4 of the README, with k replaced by u): among these irreps, Ind_{H_i}^{S_u} χ_i contains exactly
  (u − i, i) and (u − i − 1, i, 1), each once. So each irrep that occurs has a one-dimensional space of
  (H_i, χ_i)-twisted vectors for exactly one i, and Ỹ ⪰ 0 with kernel exactly K_j if and only if, for every i, the form
  Q_i on the twisted invariants V_i is ⪰ 0 with kernel exactly K_j ∩ V_i (Ỹ commutes with S_u; distinct isotypic
  parts are Ỹ-orthogonal; on V_i the form is the orthogonal sum of the forms of the two irreps up to positive factors).
- V_i is spanned by the twisted orbit sums that are not zero: with singlets |10⟩ − |01⟩ on the pairs and Dicke vectors on
  R,
  - O_{y0}: row 1@0, column y0 at copy 0, singlets, Dicke(R, q0 − y0);
  - A_{x,y0}: row 1@a, a ∈ R, column x at a, y0 at copy 0, singlets, Dicke(R − a, q0 − y0 − x);
  - B_{xy,y0}: row 1@a, a in a pair (a, a′), column x at a, y at a′, y0 at copy 0, singlets on the other pairs,
    Dicke(R, q0 + 1 − x − y − y0), with sign −1 when a is the second site of its pair;
  - G_{z0}: vacuum row, column z0 at copy 0, singlets, Dicke(R, q0 + 1 − z0).

  These are 16 types; distinct types have disjoint supports, so the nonzero ones form an orthogonal basis of V_i.
- K_j ∩ V_i is spanned by κ_{y0} = O_{y0} + A_{0,y0} + A_{1,y0} + B_{10,y0} − B_{01,y0} (y0 = 0, 1, when nonzero: in the
  rescaled frame p ⊗ H_j is 1_r ⊗ H_j, and the twisted part of the column space is spanned by the two products
  |y0⟩ ⊗ singlets ⊗ Dicke); at j = k − 2, i = 0 also K2 = σG_0 + G_1; at j = k − 1 also T0 = uσO_0 + uO_1
  − σ²A_{10} − σA_{01} − σA_{11} (i = 0, the trivial part of {d : Σ d_r s_r = 0}) and T1 = σB_{110} + B_{011} + B_{101}
  + B_{111} (i = 1, the standard part).

**The orbit calculus** (`mlib.py`). The entries of Q_i are sums over H_i-orbits of sites and of ordered pairs of sites of
{0} ∪ U (12 pair orbits, three site orbits, with explicit sizes). A vector in a frame (the marked copy and at most two
sites of U with their partners) is stored as letters on the frame and a Dicke offset; J_U^± act on Dicke vectors by
explicit factors, and singlets are killed by J_U^±. Ê_c⁻¹ is applied by Lemma E, with each of e/Δ, −½σ/Δ, σa/Δ
interpolated (Newton) in J_U⁻J_U⁺ over the nodes i′ = i − 2, …, i + 2 that exist. A vector of a frame with p
distinguished pairs and r distinguished sites of R has components only in the irreps i′ ∈ [i − p, i + p + r]
(Littlewood–Richardson) and p + r ≤ 2, so the window is exact; every numeric run also checks this, since the
interpolation asserts that the remainder after the window nodes vanishes. In exact numeric mode (σ rational) the cost is
independent of k. In symbolic mode σ is the variable of Q[σ]: every Ê_c⁻¹ carries the denominator D_c (the product of the
Δ, a or e at the nodes used), and the entries of Q_i come out as reduced rational functions of σ.

## 3. Certificates for fixed k (`cert_fixed.py`)

For each k = 2, …, 32, each j = −1, …, k and each i = 0, …, ⌊(k − 1)/2⌋ with V_i ≠ 0:
1. Q_i over Q(σ) (symbolic mode). Q_i is symmetric and every denominator of an entry is > 0 on (0, ∞).
2. Q_i κ = 0 identically in σ for every listed kernel vector κ of K_j ∩ V_i.
3. Gauge: a set S of coordinates, one per kernel vector, such that the determinant of the kernel vectors restricted to S
   is a nonzero polynomial with no root in (0, ∞) (sign test). Then for every σ > 0 the coordinate subspace W on the
   other coordinates is a complement of span(κ), and Q_i ⪰ 0 with kernel exactly span(κ) iff Q_i restricted to W is
   positive definite (as Lemma R (3) of `code/middle_k56/NOTES.md`).
4. Symmetric elimination of Q_i on W over Q(σ) in the order of the names; each pivot is a reduced rational function of
   σ, checked to be > 0 for every σ > 0. In all 3 909 forms every pivot has numerator and denominator with only
   nonnegative coefficients (the witness search and Sturm fallback of `signtest.py` were never needed).

*Why this is a proof for every σ > 0.* The entries of Q_i are polynomials in σ divided by products of the positive
polynomials D_c, so the rational functions give the true reduced form at every σ > 0. The r-th leading principal minor
is a polynomial in the entries, equal in Q(σ) to the product of the first r pivots, and each pivot is positive on
(0, ∞); so every leading principal minor is positive at every σ > 0, and Sylvester's criterion applies. With Lemma R
this proves Theorem M.

Sizes: at most 16 twisted vectors per i, at most 3 known kernel vectors; the degrees in σ of the pivots grow up to 46
and then stay there (k ≥ 12), because only five Johnson nodes ever enter. The whole run takes about 5 minutes.

## 4. Checks and controls (fixed k)

All exact unless said otherwise; `run_all.sh` reruns them.
- `check_builder.py` (a): for k = 2, …, 6, every j = −1, …, k, every i, σ ∈ {4, 1/9, 9/4} (and also 1/100, 49/16, 25
  for k ≤ 5), with the B-part scaled by 1 and by 1001/1000: every entry of the reduced form of `mlib.py` equals
  v_aᵀỸv_b, Ỹ the rescaled block built by `j1lib.py` (Lemma B, every E⁻¹ inverted over Q sector by sector) and v_a the
  explicit twisted vectors of `explicit.py` (built from the states, not from the frames); the norms agree (23 628
  entries).
- `check_builder.py` (b): for k = 2, …, 6, every j and σ ∈ {4, 1/9, 1/100}: the known kernel vectors are annihilated by
  the Lemma B block, its exact nullity equals their number, and it is ⪰ 0 (exact elimination; 90 blocks).
- `check_builder.py` (c), (d): for k ≤ 16, every j and i and two values of σ, the symbolic forms evaluated at σ equal the
  numeric ones, and the numeric forms with all Johnson nodes equal those with the window (140 428 entries each).
- `check_spectrum.py`: for k ≤ 6, every j, σ ∈ {1/4, 9}: the class sum of the transpositions of S_u preserves V_i and
  acts there with the content sums of (u − i, i) and (u − i − 1, i, 1) only; with m_λ(i) the dimensions of its
  eigenspaces, Σ_i Σ_λ dim(λ) m_λ(i) equals the number of states (exact ranks), so no irrep is missed. In floating point,
  the spectrum of Ỹ equals the union of the spectra of the parts, each with multiplicity dim λ (60 blocks).
- `coverage.py`: every (k, j, i) with k ≤ 32 and V_i ≠ 0 (3 909; the other 2 480 have V_i = 0) has a certificate with
  status "PD on complement" and exactly the expected kernel vectors; the gauge-fixed forms saved by the certificates, evaluated at σ ∈ {4, 1/9}, equal the
  restrictions of the Lemma B block to the kept explicit twisted vectors (k ≤ 6).
- Negative controls, B → (1 + 1/1000)B (`cert_fixed.py 2 16 1001/1000`, expected FAIL; `neg_points.py`). For every
  k = 2, …, 16 and every j the certificate fails at i = 0, naming σ = 10⁻⁶ (for j = −1 and k ≥ 3: σ = 10⁴ at k = 3, 10⁶
  for k ≥ 4) where a pivot of the gauge-fixed form is ≤ 0. For k = 3, …, 6 and every j, at σ = 10⁻⁶ (10⁶ for j = −1)
  the restriction of the scaled Lemma B block to the kept twisted vectors is not positive definite (exact LDL), while
  that of the true block is. The failures sit at σ → 0 or ∞, where the gap of the block tends to 0 (for j = k − 2
  this is Corollary S of the README).

## 5. Families symbolic in the sizes (`famlib.py`, `prove_fam.py`, `prove_band.py`)

**Parameters.** A reduced form is labelled by (k, j, i), or equivalently by q0 = j − i, t = k − 1 − j − i and i
(m = q0 + t unpaired copies of U). If V_i ≠ 0 then q0 ≥ −1 and t ≥ −2: every twisted vector needs a Dicke count in
[0, m] (or [0, m − 1]); the smallest lower bound comes from B_{000} (count q0 + 1 ≥ 0) and the weakest upper bound from
B_{111} (count q0 − 2 ≤ m, so t ≥ −2; t = −2 occurs only for B_{111}). A *family* fixes some of q0, t, i to integers and
leaves the others symbolic with q0 ≥ 4, t ≥ 4, i ≥ 2; σ > 0 is always symbolic. Which twisted vectors are nonzero is
decided by thresholds q0 ≥ c, t ≥ c′, i ≥ c″ with c ≤ 2, c′ ≤ 1, c″ ≤ 1 (and m ≥ 1, t + i ≥ 0), so it is constant on
a family.

**The symbolic calculus** (`famlib.py`) is that of Section 2 over Q(sizes, σ) (python-flint `fmpq_mpoly`):
- *Formal continuation*, as in `orbit3.Sym3` and `symf.SymF` of `relative-anchors/code`: the Dicke weights are the
  rational functions 2^{−p} C(m − r, q0 + d)/C(m, q0 + sh) (sh = 1 if q0 = −1, sh = t if t < 0, else 0), and the raising
  and lowering factors (q0 + d + 1) and (t − r − d + 1) are polynomials. A key outside the range of a symbolic size has
  weight 0, and the factor that would bring it back into the range vanishes; so the part of a vector inside the range
  evolves under the true operators, and such keys never contribute. At an integer point the form equals the numeric one
  of `mlib.py` divided by N0 = 2^i C(m, q0 + sh) > 0.
- *Nodes.* Ê⁻¹ uses the window nodes i′ = i + e (|e| ≤ 2) that are genuine for the part (A: e ≥ −i, e ≤ q0 + δ,
  e ≤ t − δ; B: e ≥ −i, e ≤ q0 + δ − 1, e ≤ t − δ + 1, at level j + δ). For q0 ≥ 4, t ≥ 4, i ≥ 2 the conditions on a
  symbolic size always hold, so on a family the set of genuine nodes, and whether a node lies on both sides or on one
  side only, is decided by the fixed sizes alone and is the same at every member. Genuine nodes have distinct values.
- *No division by zero at members.* The denominators that occur are q0 + sh + s and t − sh + s (s ≥ 1), m − s (s < r;
  frames with r unpaired sites are used only when m ≥ r), differences of distinct genuine node values, and the Δ, a, e
  of Lemma E at genuine nodes (> 0 for σ > 0). So at every member and every σ > 0 the symbolic computation specialises
  to the exact one.

**Certificate** (`prove_fam.prove`), per family:
1. the form over Q(sizes, σ);
2. sanity: at the 3^s sample points (each symbolic size at its minimum, minimum + 1 and minimum + 3) and σ ∈ {2/7, 5},
   N0 times the form equals the numeric form of `mlib.py` exactly (which asserts its Newton window);
3. the known kernel vectors (κ_{y0}; K2 when t = 1, i = 0; T0 when t = 0, i = 0; T1 when t = −1, i = 1) are annihilated
   identically;
4. gauge: delete O_{y0} for κ_{y0}, G_1 for K2, A_{10} for T0, B_{111} for T1; the determinant of the kernel vectors on
   these coordinates is 1, or −σ² − uσ when T0 is present, so it is sign-definite;
5. symmetric elimination of the rest, the next pivot always being the smallest diagonal entry that passes the test
   below (the pivots of any order are the ratios of the leading principal minors of a symmetric permutation of the
   matrix, so Sylvester's criterion applies). The test: after q0 = 4 + x, t = 4 + y, i = 2 + z for the symbolic sizes,
   numerator and denominator are polynomials in (x, y, z, σ) whose coefficients all have one common sign and whose part
   free of x, y, z is not zero. Then the pivot is > 0 for all x, y, z ≥ 0 and σ > 0. Every entry denominator passes the
   same test.

**Theorem F (one symbolic size).** Of the 52 families with exactly one symbolic size (the others fixed with
q0 ∈ {−1, 0, 1, 2, 3}, t ∈ {−2, …, 3}, i ∈ {0, 1}), 7 have V_i = 0, and for each of the other 45, for every member
(k, j, i) and every σ > 0, Q_i ⪰ 0 with kernel exactly K_j ∩ V_i.

**Theorem B (the bands at both ends, every k).** For every k ≥ 2, every σ > 0 and every j with j ≤ 8 or j ≥ k − 8,
the letter-1 block of (R∞)_j with squared weights (σ, 1, …, 1) is ⪰ 0 with kernel exactly K_j. With Theorem M, every j
is covered for k ≤ 32; so the statement is open only for k ≥ 33 and 9 ≤ j ≤ k − 9.

*Proof.* By Lemma R it suffices to treat each i with V_i ≠ 0.
- *Bottom band, j ≤ 8.* Then q0 = j − i ≥ −1 gives i ≤ 9. If t ≥ 4, (k, j, i) lies in the family [q0 = j − i, t ≥ 4,
  i fixed], and the families with q0 + i ≤ 8 are all proved (`prove_band.py bottom 8`: 55 families, all nonempty).
  If t ≤ 3, then
  k = j + i + t + 1 ≤ 8 + 9 + 3 + 1 = 21, and Theorem M applies.
- *Top band, k − j ≤ 8.* Then t + i = k − 1 − j ≤ 7 and t ≥ −2 give i ≤ 9. If q0 ≥ 4, (k, j, i) lies in the family
  [q0 ≥ 4, t fixed, i fixed] with t + i ≤ 7, all proved (`prove_band.py top 8`: 55 families, one of them,
  [q0 ≥ 4, t = −2, i = 0], empty). If q0 ≤ 3, then
  k = q0 + (t + i) + i + 1 ≤ 3 + 7 + 9 + 1 = 20, and Theorem M applies. ∎

`coverage_fam.py` checks this bookkeeping for every k ≤ 60: every nonempty (k, j, i) of the bands lies in a proved
family (3 052 forms with 33 ≤ k ≤ 60) or has k ≤ 32 (2 381 forms), no nonempty form lies in a family reported empty,
and always q0 ≥ −1, t ≥ −2. The band families take about 9 minutes; their largest pivots have about 2 000 terms (as
polynomials in the symbolic size and σ).

**Theorem I (the invariant part, every k and every j).** For every k ≥ 2, every j = −1, …, k and every σ > 0, the
form of the letter-1 block of (R∞)_j on the vectors invariant under the permutations of the unmarked copies (the
twisted space V_0) is ⪰ 0 with kernel exactly K_j ∩ V_0.

*Proof.* V_0 is the multiplicity space of the trivial irrep of S_u (H_0 = S_u, χ_0 = 1, and (u − 1, 0, 1) is not a
partition). A member (k, j, 0) with V_0 ≠ 0 has q0 ≥ −1 and t ≥ −2. If q0 ≥ 4 and t ≥ 4 it lies in the family
[q0 ≥ 4, t ≥ 4, i = 0] (`prove_lowirr.py 0`); if exactly one of q0, t is ≥ 4 it lies in a family of Theorem F; if
both are ≤ 3, then k = q0 + t + 1 ≤ 7 and Theorem M applies. `coverage_irr.py` checks this for every k ≤ 60. ∎

The family [q0 ≥ 4, t ≥ 4, i = 0] has six pivots (largest about 3 200 terms in (x, y, σ)) and takes about a second.

Negative control: `prove_fam.py 1 1 1001/1000` (B → (1 + 1/1000)B in the symbolic calculus) fails, as it must: 42 of
the 45 nonempty families no longer pass the certificate (3 still do; there the gap exceeds what the change of B
removes). The explicit points where the scaled block is not positive semidefinite are those of Section 4.

## 6. What is not covered

- **Every k for the middle values of j.** For k ≥ 33 and 9 ≤ j ≤ k − 9 the statement is open. Its part i = 0 is
  Theorem I; the parts with i ≥ 1 there have q0 ≥ 4 and t ≥ 4, so they need the families [q0 ≥ 4, t ≥ 4, i = 1] and
  [q0 ≥ 4, t ≥ 4, i ≥ 2] (three and four variables with σ). In development runs (`logs/dev/`) the families
  [q0 = −1, t ≥ 4, i ≥ 2] and [q0 = 0, t ≥ 4, i ≥ 2] passed, but for [q0 = 1, t ≥ 4, i ≥ 2] the eliminated entries
  grew to about 10⁴ terms in three variables, and the run was stopped after a few minutes in both pivot orders tried.
  A probe of [q0 ≥ 4, t ≥ 4, i = 1] (greedy order) reached 11 of its 14 pivots in about 6 minutes, all of them passing
  the sign test, with pivots of about 4·10⁴ terms, and was stopped by its 12-minute limit (`logs/dev/probe_ss1.log`);
  so the parts i = 1 are not claimed. The generic family [q0 ≥ 4, t ≥ 4, i ≥ 2], with four variables, was not
  attempted. Splitting V_i into its two irreps
  (at most 12 + 4 vectors instead of 16; on the states with the row letter in U the common kernel of the four
  equivariant maps (1@r, Y) ↦ Y, [r ∉ Y](Y ∪ r), [r ∈ Y](Y − r), [r ∈ Y]Y contains the (u − i − 1, i, 1) part, and
  equals it when the multiplicities are generic), or a determinant certificate in the style of Proposition C of the
  README, would be the next things to try.
- More than one marked copy, two weight values with both levels large, and general weights: not attempted here.
- Nothing here concerns the type I or type II blocks with a marked copy, or Conjecture 12 as such.

## 7. A second check

`markcheck.py`, `projcheck.py` and `famcheck.py` were written separately from the certificates. `markcheck.py` uses
none of the reduction: only the operators of `j1lib.py`, python-flint, numpy and the exact positive-definiteness test of
`code/middle_k56/crosscheck.py` (copied here).
- *A second builder.* The B-part of Lemma B is a Gram form of E⁻¹ on the vectors T_r|y⟩ and J⁺|y⟩ (E is symmetric and
  J⁻ = (J⁺)ᵀ), arranged with the index pattern of `J1.entry`; it is computed with flint matrices, one energy sector at a
  time. It agrees with `J1.matrix` entry by entry for k = 4, 5, 6, every j and two values of σ.
- *End to end, exact.* For k = 7 and 8 at σ ∈ {1/100, 1/9, 49/25, 4, 400}, for k = 9 at σ ∈ {1/100, 4, 400}, and every
  j = −1, …, k: the known kernel vectors (unrescaled frame) are annihilated exactly, their number is
  C(k, j) + [j = k − 2] + (k − 1)[j = k − 1], and the block without a set of pivot coordinates of them is positive
  definite by the exact test (A_int = D·A for the lcm D of the denominators; an integer rounding Z of the inverse of a
  floating-point Cholesky factor; Z A_int Zᵀ, which is congruent to A_int, strictly diagonally dominant in exact
  integer arithmetic). So the block is ⪰ 0 with kernel exactly K_j at these points. Negative control: with
  B → (1 + 1/1000)B the test refuses the block for some j at σ = 10⁻⁶ (k = 5 and 7), and at k = 5 the scaled builder
  agrees with `J1` scaled in the same way.
- *The reduced forms.* For k = 7 and 8 at σ ∈ {4, 1/9, 1/100}, for k = 9 at σ ∈ {4, 1/100}, every j and every i: every
  entry of the numeric reduced form of `mlib.py` equals v_aᵀỸv_b, with Ỹ the rescaled block of the second builder and
  v_a the explicit twisted vectors of `explicit.py`, and the norms agree (20 392 entries). With B scaled by
  1 + 1/1000 in `mlib.py` alone, 240 of the 256 entries of the form at k = 7, j = 3, i = 1 differ.
- *The families far from the sample points.* The certificates of Section 5 compare each family form with the numeric
  one at the sizes min, min + 1 and min + 3. `famcheck.py` recomputes the forms of all 155 nonempty families used by
  Theorems B, F and I (45 with one symbolic size, 55 of the bottom band, 54 of the top band, and [q0 ≥ 4, t ≥ 4, i = 0])
  and compares them, times N0, with the numeric forms entry by entry at symbolic sizes 6, 11, 23 and 40 (i: 2, 6, 13;
  for the family with two symbolic sizes the pairs in {5, 23} × {6, 31}) and σ ∈ {3/11, 7/2}: 1 192 points, with k up
  to about 66.

`run_check.sh` runs all of this, one Python process at a time; the output is `logs/run_check.out`.

## Files

| File | Content |
|---|---|
| `mlib.py` | the orbit calculus with one marked copy (Section 2), numeric (σ rational) and symbolic in σ: the 16 twisted vectors, the orbits, Lemma E with Newton interpolation (window asserted), the known kernel vectors |
| `explicit.py` | explicit objects at small k, independent of the frames: the rescaled Lemma B block from `j1lib.py`, the explicit twisted vectors, the known kernel vectors of the whole block |
| `cert_fixed.py` | the fixed-k certificates (Section 3); with B scaled, the negative control |
| `check_builder.py` | the reduced forms against the Lemma B builder, end-to-end kernels, symbolic against numeric, window against all nodes (Section 4) |
| `check_spectrum.py` | completeness of the reduction: class-sum split, exact dimension counts, spectra (Section 4) |
| `coverage.py` | coverage of every (k, j, i) with k ≤ 32; the saved certified forms against the Lemma B builder |
| `neg_points.py` | the negative controls at explicit points, from the Lemma B builder |
| `famlib.py` | the orbit calculus with the sizes symbolic (Section 5) |
| `prove_fam.py` | the family certificate; `prove_fam.py 1 1` runs every family with one symbolic size |
| `prove_band.py` | the families of the two bands (Theorem B) |
| `coverage_fam.py` | the bookkeeping of Theorem B for k ≤ 60 |
| `prove_lowirr.py` | the families [q0 ≥ 4, t ≥ 4, i = b] for given b (b = 0 for Theorem I) |
| `coverage_irr.py` | the bookkeeping of Theorem I for k ≤ 60 |
| `j1lib.py`, `flintfield.py`, `ufield.py`, `signtest.py`, `guard.sh` | copied unchanged from `relative-anchors/code/` (`middle_j1/`, the code folder itself, `middle_j2_l1/`, `middle_k56/`, `no_letter1/`) |
| `markcheck.py`, `projcheck.py`, `famcheck.py`, `crosscheck.py`, `run_check.sh` | the second check of Section 7 (`bash run_check.sh > logs/run_check.out`); `crosscheck.py` is copied from `code/middle_k56/` for its exact test |
| `run_all.sh` | reruns everything, one Python process at a time under `guard.sh`: `bash run_all.sh > logs/run_all.out` (14 steps, about 19 minutes), ending with ALL PASS |
| `logs/run_all.out`, `logs/run_*.log` | the final run (its pickles `logs/*.pkl`, with the certified fixed-k forms and the family results, are regenerated by `run_all.sh` and are not in the repository); `logs/run_check.out` and `logs/run_check_*.log` are the second check |
| `logs/dev/`, `dev/` | development runs and scripts, not used by any proof: the first full run (fixed k only, ALL PASS), a full run whose band steps did not start (an import error in `prove_fam.py`, fixed before the final run), the stopped two-size runs, the probe of [q0 ≥ 4, t ≥ 4, i = 1], timing probes, first smoke tests |
