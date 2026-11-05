# All of (R∞) at five and six copies when the weights take at most two values

5 November 2026. Nothing here is refereed.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Theorem A.** k = 5, weights with at most two distinct values. (a) The letter-1 block of (R∞)₃ is ⪰ 0 with kernel exactly span{p ⊗ \|1@Y⟩ : \|Y\| = 3} ⊕ span{vac ⊗ J⁻\|11111⟩}. (b) For \|T\| = 1, 2, 3, the type I block of (R∞)₃ labelled by T is ⪰ 0 with kernel exactly span{p ⊗ y : y a column of its letter-1 rows}. (c) For \|T\| = 2, the type I block of (R∞)₄ labelled by T is ⪰ 0 with kernel exactly span{p ⊗ f_h : h ∈ U} | **proved**, computer-assisted, exact: one certificate per configuration (52 in all), each valid for every weight ratio q > 0 (Sections 2–4) |
| 2 | **Theorem 5₂.** For k = 5 and every weight vector with at most two distinct weights, all of (R∞) holds: (R∞)_j ⪰ 0 for j = −1, …, 5; the kernel of (R∞)₃ is exactly p ⊗ H₃ ⊕ span{vac ⊗ J⁻\|11111⟩} | **proved**, given row 1 and the results of `relative-anchors/README.md` and `code/middle_j2_l1/` listed in Section 4 (all computer-assisted in part, none refereed) |
| 3 | **Corollary 5₂.** Conjecture 12 of `single-excitation-range/` holds for every p ∈ H₁ supported on at most five copies whose nonzero weights take at most two distinct values | **proved**, given row 2, Corollary 4 of `relative-anchors/` and Lemma 13 with the remark that (R2) always holds (`single-excitation-range/`) |
| 4 | Lemma R (reductions): permutations of the copies, homogeneity of degree j + 1 − \|T\| of the rescaled block, and the gauge | **proved** by hand (Section 2); the reduction of every block and every two-value weight vector to a certified configuration is **checked exactly** for all 1 152 pairs (block, level pattern) at one weight ratio each (`coverage.py`) |
| 5 | The certified blocks (built over Q(q) by `ylib.py`, from the definitions) equal the blocks built by `j1lib.py` (Lemma B) | **exact check** at 156 (configuration, q) pairs, every configuration, q ∈ {4/9, 25/4, 1/10⁴}; equal to the `ylib.py` point builder at 52 points; the two builders agree at general weights (two weight vectors per family, B-part scaled by 1 and by 1001/1000); homogeneity checked (`check_builders.py`) |
| 6 | The certified leading minors equal the determinants of the leading principal submatrices of the gauge-fixed Lemma B block | **exact check**: all 1 400 certified minors, each at two values of q (2 800 comparisons; `check_minors.py`) |
| 7 | End to end at points: the known vectors are annihilated, the nullity equals their number, the gauge-fixed block is positive definite in two different gauges | **exact check** at 260 (configuration, q) pairs, q from 10⁻⁶ to 10⁶ (`check_end2end.py`); not a proof |
| 8 | Negative controls, B → (1 + 1/1000)B: every certificate of a two-level configuration fails (42 of 42), and at an explicit rational point per configuration the scaled gauge-fixed block is not positive definite while the true one is | **exact** (`run_all.sh` steps `neg_cert_*`, `neg_points.py`) |
| 9 | Observation: in all 52 certificates at k = 5 and all 76 at k = 6 every leading principal minor is a ratio of two polynomials in q with only nonnegative coefficients (the Sturm fallback was never used) | **exact** observation; consistent with Conjecture SF of Section 3g of the README; not used beyond the sign test |
| 10 | **Theorem B.** k = 6, weights with at most two distinct values. Each of the following blocks is ⪰ 0 with kernel exactly the span of its known kernel vectors (the vectors p ⊗ y, y a column of its letter-1 rows, and for (c) also vac ⊗ J⁻\|111111⟩): (a) the letter-1 block of (R∞)₃; (b) the type I blocks of (R∞)₃ with \|T\| = 1, 2, 3; (c) the letter-1 block of (R∞)₄; (d) the type I blocks of (R∞)₄ with \|T\| = 1, 2, 3, 4; (e) the type I blocks of (R∞)₅ with \|T\| = 2, 3 | **proved**, computer-assisted, exact: one certificate per configuration up to the level exchange, 76 in all (56 direct, 20 after an isotypic reduction), each valid for every q > 0 (Section 7) |
| 11 | **Theorem 6₂.** For k = 6 and every weight vector with at most two distinct weights, all of (R∞) holds; the kernel of (R∞)₃ is exactly p ⊗ H₃ and that of (R∞)₄ is exactly p ⊗ H₄ ⊕ span{vac ⊗ J⁻\|111111⟩} | **proved**, given row 10 and the cited results (Section 7) |
| 12 | **Corollary 6₂.** Conjecture 12 of `single-excitation-range/` holds for every p ∈ H₁ supported on at most six copies whose nonzero weights take at most two distinct values | **proved**, given rows 3 and 11 and Lemma 13 there |
| 13 | k = 6 checks: certified blocks against `j1lib.py` and `ylib.py`, certified minors and reduced parts against the Lemma B block, end to end at points, coverage, negative controls (every two-level certificate fails with B → (1 + 1/1000)B) | **exact checks** (Section 7); not proofs |
| 14 | A second, independent certificate for every configuration that is certified directly (52 at k = 5, 56 at k = 6): the isotypic reduction of Section 7, with the complement W in place of the deletion gauge; and denser samples of the exact coverage identity at k = 6 (1 088 more pairs) | **exact**, all pass (`run_extra.sh`, `logs/run_extra.out`, ALL PASS); not needed for any proof |
| 15 | k ≥ 7 at two weight values; three or more weight values; general weights | **open**; not attempted beyond one timing probe at k = 7 (Section 8) |
| 16 | A second check written separately (Section 9), using only `j1lib.py`: at 12 random two-value points at k = 5 and 24 at k = 6, every torus block of every (R∞)_j, j = −1, …, k (grouped by the letters 2 at each copy; entries between groups vanish on a random sample), has its known kernel vectors annihilated and is positive definite without a set of pivot coordinates of them (exact certificate), so it is ⪰ 0 with exactly that kernel, and the kernel dimensions add up to \|H_j\| (+ 1 at j = k − 2, + k − 1 at j = k − 1); the blocks saved by all 128 certificates are ⪰ 0 with exactly their known kernels at 31 values of q from 10⁻⁸ to 10⁸ | **exact check** at points (`crosscheck.py`, `run_crosscheck.sh`); not a proof |

Nothing here says that Conjecture 12 holds in general, or anything about the NPT (Werner distillability) problem.
Computer-assisted means: exact rational-function arithmetic in one variable (python-flint `fmpq_poly`), exact rational
arithmetic, and sign checks of coefficients. No floating point enters any proof (floating point is used only to
propose candidate points for the negative controls, whose signs are then evaluated exactly).

## 1. Setting

Notation is that of `relative-anchors/README.md` (Sections 1 and 2) in the repository `werner-npt-research`. In
Sections 1 to 6 there are k = 5 copies (Section 7 treats k = 6 with the same tools), s_t = w_t² and
p = Σ_t w_t |1_t⟩. (R∞)_j is Y = A − B compressed to G1_j ⊕ G0_{j+1}. Its torus blocks
are the letter-1 block, the type I blocks (labelled by the set T of copies that carry the letter 2 exactly once, in
row or column) and the type II blocks. U = [k] ∖ T.

We work in the rescaled frame of `ylib.py`: Ỹ = DYD with D = diag(μ), μ(X) = Π_t w_t^{N_t(X)}, where N_t(X) is the
number of letters 1 at copy t in the row and the column of the state X. In this frame every entry is a rational function
of the squared weights (Section 3f of the README; it is also what `ylib.py` computes, from M, V_p, D_p, R_p and the
partial transpose). A vector v is in the kernel of Y if and only if D⁻¹v is in the kernel of Ỹ, and positive diagonal
rescaling keeps the inertia.

The five block families of Theorem A, with their sizes and known kernels:

| family | block | states | known kernel (rescaled frame) | gauge-fixed size |
|---|---|---|---|---|
| L3 | letter-1 block of (R∞)₃: (1@r, 1@Y), \|Y\| = 3; (vac, 1@Z), \|Z\| = 4 | 50 + 5 = 55 | the indicators of the 10 column groups {(1@r, 1@Y) : r}; the vector (s_h)_h on the states (vac, 11111 − e_h) (Proposition K2) | 44 |
| T31 | type I, (R∞)₃, T = {t}: (2@t, 1@Y), \|Y\| = 3; (1@r, 2@t 1@X), X ∌ t, \|X\| = 2; (vac, 2@t 1@Z), \|Z\| = 3 | 10 + 30 + 4 = 44 | 6 column-group indicators | 38 |
| T32 | type I, (R∞)₃, \|T\| = 2 | 12 + 15 + 3 = 30 | 3 column-group indicators | 27 |
| T33 | type I, (R∞)₃, \|T\| = 3 | 9 + 5 + 2 = 16 | 1 column-group indicator | 15 |
| T42 | type I, (R∞)₄, \|T\| = 2 (states P_{r,h}, R_t, Z_{t,h}, Q of the Remark in Section 3 of the README) | 15 + 2 + 6 + 1 = 24 | 3 column-group indicators (p ⊗ f_h) | 21 |

The column-group indicators are the vectors p ⊗ y = Σ_r w_r (1@r, y): since μ(1@r, y) = w_r Π_{t : y_t = 1} w_t, the
rescaled vector D⁻¹(p ⊗ y) is a constant times the indicator of the column group {(1@r, y) : r ∈ [k]}. That p ⊗ y is in
the kernel is Lemma N of the README (Σ_r w_r 𝒯_{(r,1)} = E_c). For the vacuum vector, μ(vac, 11111 − e_h) =
Π_{t≠h} w_t, so D⁻¹(vac ⊗ J⁻|11111⟩) = (Π_t w_t)⁻¹ Σ_h s_h (vac, 11111 − e_h).

Before any certificate we checked exactly, at two-value rational points, that each block has no kernel beyond these
vectors (`explore_kernel.py`, three weight ratios per representative configuration, and `check_end2end.py`, five
ratios per configuration): in every case the nullity equals the number of known vectors.

## 2. Reductions (Lemma R)

**Lemma R.** Fix one of the families, with j and m = |T|, and put d = j + 1 − m.
1. *Permutations.* A permutation π of the copies, applied to rows and columns at once and with the weights permuted
   along, maps the block labelled by T at the weights s onto the block labelled by π(T) at the weights s ∘ π⁻¹,
   entry by entry, in the rescaled frame.
2. *Homogeneity.* Multiplying every squared weight by c > 0 multiplies the rescaled block by c^d.
3. *Gauge.* Let v_1, …, v_n be the known kernel vectors of the block and delete coordinates δ_1, …, δ_n such that
   v_i(δ_i) ≠ 0 and v_l(δ_i) = 0 for l ≠ i. Then the block is ⪰ 0 with kernel exactly span{v_i} if and only if the
   principal submatrix on the remaining coordinates is positive definite.

*Proof.* (1) M^{⊗k}, the map p ↦ R_p and the partial transpose commute with permutations of the copies, and so do the
torus labels and μ. (2) In Lemma B, E_c is homogeneous of degree 2 in w and 𝒯_a of degree 1, so 2^Z Y is homogeneous of
degree 0; every state of the block carries the same number d of letters 1 (row and column together: the block has
j + 1 letters, m of them equal to 2), so μ(X)μ(X′) has degree 2d in w, that is d in s. (3) The coordinate subspace W
on the kept coordinates is a complement of K = span{v_i}: a combination Σ c_i v_i that vanishes on every δ_l has
c_l v_l(δ_l) = 0 for every l. For x = w + v with w ∈ W and v ∈ K ⊆ ker Ỹ, ⟨x, Ỹx⟩ = ⟨w, Ỹw⟩. So Ỹ ⪰ 0 if and only if
Ỹ|_W ⪰ 0, and ker Ỹ = K if and only if Ỹ|_W is positive definite (a w ∈ W with ⟨w, Ỹw⟩ = 0 lies in ker Ỹ, and
ker Ỹ ∩ W = 0 when ker Ỹ = K). ∎

The gauge used: for every column group delete the state (1@4, y) (the last copy, which lies in U in every family), and
for L3 also the state (vac, 11111 − e_0). Each indicator is supported on its own column group and the vacuum vector on
the vacuum states, so the condition of (3) holds (asserted in `k5lib.known_kernel`). `check_end2end.py` also uses the
gauge with (1@0, y) deleted.

**Configurations.** Put T₀ = {0, …, m − 1} and U₀ = {m, …, 4}. Configuration (a, b), 0 ≤ a ≤ m, 0 ≤ b ≤ 5 − m, has
squared weight 1 on the copies 0, …, a − 1 and m, …, m + b − 1, and squared weight q on the others. Now take any block of
the family, labelled by T, and any weight vector with at most two distinct squared weights: α on a set A of copies and β
on the rest (A = [5] or A = ∅ if all weights are equal). Put a = |T ∩ A| and b = |U ∩ A|, and let π map T ∩ A, T ∖ A,
U ∩ A, U ∖ A in increasing order onto 0, …, a − 1, then a, …, m − 1, then m, …, m + b − 1, then the rest. By (1) and
(2) the block is, up to the relabelling of the states by π, α^d times the block of configuration (a, b) at q = β/α. So
it suffices to prove, for every family and every configuration (a, b), that the gauge-fixed block is positive definite
for every q > 0. All 52 configurations are certified (6 for L3, 10 for T31, 12 each for T32, T33, T42); the level
exchange, which would halve this, is not needed. `coverage.py` runs over every block of every family (every T with
|T| = m) and every level pattern A ⊆ [5] (1 152 pairs), checks that (a, b) is certified, and checks the identity above
exactly, entry by entry, at α = 4/9, β = 25/16, with the left side built by `j1lib.py` (Lemma B) and the right side
from the certified block.

## 3. The certificate

`cert5.py FAMILY a b` works over Q(q) with exact univariate rational functions (`ufield.py`, python-flint `fmpq_poly`,
reduced fractions with monic denominators):
1. the rescaled block over Q(q), built by `ylib.py` with the squared weights 1 and q (the inverse of each block of
   D_p by exact Gauss–Jordan over Q(q)); every denominator of an entry is checked to be > 0 on (0, ∞);
2. the known kernel vectors are annihilated, identically in q;
3. the gauge of Section 2;
4. symmetric Gaussian elimination of the gauge-fixed block in the natural order of `k5lib.block_states` (letter-2 rows,
   then letter-1 rows, then vacuum states); the r-th leading principal minor is the product of the first r pivots,
   as a reduced rational function, and it is checked to be > 0 for every q > 0 (`signtest.py`): numerator and
   denominator have only nonnegative coefficients (and are not zero); if not, an exact evaluation at candidate points
   looks for a q > 0 where the minor is ≤ 0 (a witness, which ends the run with FAIL), and if there is none, an exact
   Sturm count of the roots in (0, ∞) decides. In the true certificates the first test always succeeded; the other
   two are used only by the negative controls (a Sturm count at degree 150 takes minutes, so the witness is looked
   for first).

*Why the rational functions give the true block at every q > 0.* The blocks of D_p are polynomial in q and positive
definite for every q > 0 (Lemma 1 of the README), so the entries of their inverses over Q(q) have denominators dividing
their determinants, which do not vanish on (0, ∞); every entry of Ỹ is a polynomial in these entries and in polynomial
data, so its value at q₀ > 0 is the value of the rational function at q₀ (step 1 also checks the denominators
directly). A leading minor is a polynomial in the entries, so its value at q₀ is the value at q₀ of the rational
function computed in step 4 (the identity "minor = product of pivots" holds in Q(q), and each such product is
defined at q₀ because its denominator is positive there). So every leading principal minor of the gauge-fixed block is
positive at every q > 0, and by Sylvester's criterion the gauge-fixed block is positive definite for every q > 0.

## 4. Proofs

*Proof of Theorem A.* By Lemma R it suffices that, for every family and every configuration (a, b), the gauge-fixed
block is positive definite for every q > 0, which is what the 52 certificates of Section 3 show (`run_all.sh`,
`logs/run_cert_k5_*.log`). Lemma R (3) then gives the kernels. ∎

*Proof of Theorem 5₂.* The weights take at most two values, so the weights of every subset of the copies do too.
Every torus block of every (R∞)_j at k = 5 is covered:
- (R∞)_{−1}: Proposition M. (R∞)_0: Theorem B0 (k ≥ 3, all weights).
- (R∞)_1: the Corollary of Section 3h of the README (two weight values, every k ≥ 4: Theorems 2V, EQ, 2L, D₁ and
  Proposition T2′).
- (R∞)_2: the Corollary of Theorem L in `code/middle_j2_l1/` (5 ≤ k ≤ 10, two weight values: Theorem L for the
  letter-1 block, Theorems F1b and F2b of `code/middle_j2/` for |T| = 2 and 1, Theorem DJ for |T| = 3, Proposition
  T2′ for type II).
- (R∞)_3: Theorem A (a) and (b) for the letter-1 block and |T| = 1, 2, 3; Theorem DJ for |T| = 4 (no letter 1); type
  II by Proposition T2′.
- (R∞)_4 = (R∞)_{k−1}: Theorem T (letter-1 block); |T| = 1 by Theorems EU and TL (the four U-weights take at most two
  values; the weight of t is arbitrary); |T| = 2 by Theorem A (c); |T| = 3, 4, 5 by Theorem U2, Proposition U1 and
  Proposition U0; type II by Proposition T2.
- (R∞)_5 = (R∞)_k: Theorem K.

Each block is ⪰ 0, so each (R∞)_j is. In (R∞)₃ the type II blocks and the block |T| = 4 are positive definite, and
the columns y ∈ H₃ of the letter-1 rows are distributed over the letter-1 block and the type I blocks with |T| ≤ 3
(the block of (1@r, y) is labelled by the positions of the letters 2 of y), so its kernel is exactly
p ⊗ H₃ ⊕ span{vac ⊗ J⁻|11111⟩}. ∎

*Proof of Corollary 5₂.* Local unitaries that fix |0⟩ bring p to Σ_t w_t|1_t⟩ with w_t ≥ 0, and the weights w_t are
the norms of the single-copy parts of p, so their nonzero values are unchanged. If at most four weights are nonzero,
this is Corollary 4 of `relative-anchors/` (all weights). If exactly five are nonzero, let G be these five copies and
p = p_G ⊗ |0⟩^{⊗m}. By Theorem 5₂, (R∞) holds for p_G; (R2) always holds (the remark before Theorem 16 of
`single-excitation-range/`); so by Lemma 13 there, Y_c(p) ⪰ 0 for every m ≥ 0 and every c. ∎

This rests on the cited theorems as they are stated in the README and in `code/middle_j2_l1/`; several of them are
computer-assisted, and none is refereed.

## 5. Checks and controls

All exact. `run_all.sh` reruns everything, one Python process at a time under `guard.sh` (2.5 GB), and ends with
ALL PASS (`logs/run_all.out`: 247 steps, 1 917 s for k = 5 and k = 6 together; the k = 5 part takes about 7
minutes).
- *Builders* (`check_builders.py`). (a) For every family and every configuration, the block saved by the certificate,
  evaluated at q ∈ {4/9, 25/4, 1/10⁴}, equals the block built by `j1lib.py` (Lemma B, every E⁻¹ inverted over Q
  sector by sector) at the weights w = √s, rescaled by μ and divided by 2^Z. (b) It equals the block built by
  `ylib.py` over Q at q = 9/16. (c) At two general weight vectors (five distinct squared weights) the two builders agree
  for every family, with the B-part scaled by 1 and by 1001/1000 (so the negative controls use the right block).
  (d) Homogeneity of degree j + 1 − |T|.
- *Minors* (`check_minors.py`). Every saved leading minor, at q ∈ {4/9, 49/16}, equals the determinant (flint) of the
  corresponding leading principal submatrix of the gauge-fixed Lemma B block: 1 400 minors, 2 800 comparisons.
- *End to end* (`check_end2end.py`). For every configuration and q ∈ {1/4, 9/4, 1/10⁶, 10⁶, 121/100}, from the Lemma B
  block: the known vectors are annihilated, the nullity (exact rank over Q) equals their number, and the gauge-fixed
  block is positive definite in the gauge of the certificates and in the gauge with (1@0, y) deleted.
- *Coverage* (`coverage.py`), Section 2.
- *Negative controls* (B → (1 + 1/1000)B). `cert5.py FAMILY a b 1001/1000` prints FAIL for each of the 42 configurations
  with two levels, naming a leading minor and a rational q at which it is ≤ 0 (`run_all.sh` steps `neg_cert_*`, expected
  FAIL). `neg_points.py` gives, for each of them, a rational square q (10⁻⁶ for 38 configurations, 10⁴ for 4) at
  which the gauge-fixed block built by `j1lib.py` with `bscale` = 1001/1000 is not positive definite while the true
  one is. At one level (equal weights) the scaled certificates still pass: there the gap of the block exceeds what a
  change of 1/1000 in B can remove. The failures occur as q → 0 or q → ∞, where the gap of the true block tends to
  zero (for the letter-1 block this is Corollary S of Section 3e of the README).

## 6. Sizes

From the final run (`logs/run_cert_k5_*.log`): the number of configurations, of leading minors checked, the largest
degrees in q of a numerator and of a denominator of a leading minor, and the total time.

| family | states → gauge-fixed | configurations | minors | largest num. degree | largest den. degree | time |
|---|---|---|---|---|---|---|
| L3 | 55 → 44 | 6 | 264 | 274 | 160 | 15 s |
| T31 | 44 → 38 | 10 | 380 | 248 | 172 | 28 s |
| T32 | 30 → 27 | 12 | 324 | 178 | 145 | 8 s |
| T33 | 16 → 15 | 12 | 180 | 86 | 78 | 2 s |
| T42 | 24 → 21 | 12 | 252 | 123 | 84 | 2 s |

With one weight variable the certificates are small: the four-parameter certificates at k = 5 that did not finish
(Sections 3, 3f and 3h of the README) are not needed for two weight values.

## 7. Six copies

At k = 6 and two weight values, the blocks not covered by earlier results were (Section 3h of the README and the
results cited in Section 4): in (R∞)₃ = (R∞)_{k−3} and (R∞)₄ = (R∞)_{k−2} the letter-1 block and every type I block
that contains a letter 1 (|T| ≤ 3, resp. |T| ≤ 4), and in (R∞)₅ = (R∞)_{k−1} the type I blocks with |T| = 2, 3
(|U| = 4, 3). We name the families jJmM: the type I block of (R∞)_J with |T| = M (M = 0: the letter-1 block). The same
scripts are used with `k=6`.

**Known kernels.** They are the column-group indicators (p ⊗ y), and for j4m0 (j = k − 2) also the vacuum vector of
Proposition K2. For j3m0 (j = k − 3) there is no vacuum kernel vector. Before the certificates, `explore_kernel.py`
checked exactly, at 228 two-value points (three ratios for each representative configuration of each family), that
every block has nullity exactly the number of known vectors and that the gauge-fixed block is positive definite.

**Level exchange.** At k = 6 only one configuration of each pair (a, b), (m − a, 6 − m − b) is certified
(`k5lib.rep_configs`). Configuration (a, b) at ratio q has squared weight 1 on a set A₀ and q on its complement;
dividing all squared weights by q and sorting the complement first (Lemma R (1), (2)) gives configuration
(m − a, 6 − m − b) at ratio 1/q. So one of the two suffices. `coverage.py` uses this for every pattern whose own
configuration is not certified.

**The isotypic reduction (cert_sym.py).** For the four largest families (j4m1, j4m0, j3m1, j3m0, with 80 to 135
states) the direct certificate becomes slow: a direct run takes about 1 s per configuration at 20 states and up to
about a minute at 60 states, growing roughly like the fifth power of the size. There we use the symmetry of the
configuration.
Let A be the set of copies of squared weight 1 in configuration (a, b) and
G = Sym(T ∩ A) × Sym(T ∖ A) × Sym(U ∩ A) × Sym(U ∖ A). G permutes the states of the block, and Ỹ is
G-invariant (Lemma R (1); checked identically in q on the generating transpositions in every run).
- *A G-stable complement.* Let K be the span of the known kernel vectors and W the space of vectors x with
  ⟨1_y, x⟩ = 0 for every column group y (and, for j4m0, ⟨1_vac, x⟩ = 0, 1_vac the indicator of the vacuum states).
  G permutes the column groups and the vacuum states, so W is G-stable. The pairing between the vectors 1_y, 1_vac and
  the kernel vectors is diagonal with positive entries (k on 1_y, Σ_h s_h on the vacuum vector), so W ∩ K = 0 and
  dim W + dim K = n. As in Lemma R (3), the block is ⪰ 0 with kernel exactly K if and only if Ỹ is positive definite
  on W.
- *Young symmetrizers.* For an irrep λ = (λ₁, …, λ₄) of G (one partition per factor), fill each Young diagram row by
  row with the copies of its factor in increasing order, and let e_λ be the product of the row symmetrizer and the
  column antisymmetrizer. On a permutation module, e_λ x e_λ is a multiple of e_λ for every x in the group algebra, so
  e_λ is zero on the other isotypic parts and has rank one on each copy of λ. Hence e_λW has dimension mult(λ, W), and
  its vectors are s₀ ⊗ c for one fixed vector s₀ of λ. Irreps of symmetric groups are absolutely irreducible and
  self-dual, so Ỹ acts on the λ-isotypic part of W as ⟨·,·⟩ ⊗ Q_λ, and for any basis B_λ of e_λW,
  B_λỸB_λᵀ = |s₀|² C Q_λ Cᵀ with C invertible. So Ỹ|_W ≻ 0 if and only if B_λỸB_λᵀ ≻ 0 for every λ. Every run checks
  Σ_λ dim(λ) · dim(e_λW) = dim W, so no irrep is missed. (This is the argument of Lemma Y of `code/middle_j2_l1/`,
  with a different complement.)
- *Certificate.* `cert_sym.py` builds Ỹ over Q(q) with `ylib.py`, checks the kernel identities and the G-invariance,
  multiplies Ỹ by den, the lcm of the entry denominators (a polynomial that is > 0 on (0, ∞)), builds a basis of each
  e_λW by applying e_λ to a basis of W and row-reducing over Q (`symred.py`), forms F_λ = B_λ(den Ỹ)B_λᵀ over Q[q],
  and checks every leading principal minor of every F_λ with the sign test of Section 3. Since den > 0 on (0, ∞), each
  B_λỸB_λᵀ is then positive definite for every q > 0.

Sizes (final run, `logs/run_cert_k6_*.log`; certified size: the gauge-fixed block, or dim W for the reduced families;
for these, the parts with nonzero e_λW summed over the configurations, and the largest part; minors and time summed
over the configurations):

| family | block | states | known kernel | certified size | configurations | method | minors | largest num. degree | time |
|---|---|---|---|---|---|---|---|---|---|
| j4m4 | (R∞)₄, \|T\| = 4 | 20 | 1 | 19 | 8 | direct | 152 | 106 | 7 s |
| j3m3 | (R∞)₃, \|T\| = 3 | 21 | 1 | 20 | 8 | direct | 160 | 115 | 7 s |
| j5m3 | (R∞)₅, \|T\| = 3 | 31 | 3 | 28 | 8 | direct | 224 | 157 | 9 s |
| j5m2 | (R∞)₅, \|T\| = 2 | 35 | 4 | 31 | 8 | direct | 248 | 205 | 11 s |
| j4m3 | (R∞)₄, \|T\| = 3 | 39 | 3 | 36 | 8 | direct | 288 | 224 | 23 s |
| j3m2 | (R∞)₃, \|T\| = 2 | 50 | 4 | 46 | 8 | direct | 368 | 324 | 92 s |
| j4m2 | (R∞)₄, \|T\| = 2 | 60 | 6 | 54 | 8 | direct | 432 | 393 | 225 s |
| j4m1 | (R∞)₄, \|T\| = 1 | 80 | 10 | 70 | 6 | reduced: 30 parts, largest 16 | 202 | 394 | 8 s |
| j4m0 | (R∞)₄, letter-1 | 96 | 15 + 1 | 80 | 4 | reduced: 27 parts, largest 10 | 108 | 162 | 4 s |
| j3m1 | (R∞)₃, \|T\| = 1 | 90 | 10 | 80 | 6 | reduced: 30 parts, largest 18 | 226 | 600 | 10 s |
| j3m0 | (R∞)₃, letter-1 | 135 | 20 | 115 | 4 | reduced: 27 parts, largest 13 | 150 | 293 | 8 s |

In all 76 certificates every leading minor (of the gauge-fixed block, or of a part F_λ) has numerator and denominator
with only nonnegative coefficients; no witness search and no Sturm count was needed.

**Theorem B** (row 10 of the status table). *Proof.* By Lemma R and the level exchange, every block of the eleven
families at every two-value weight vector reduces to a certified configuration, and each certificate shows that the
gauge-fixed block (direct families) or Ỹ restricted to W (reduced families) is positive definite for every q > 0. ∎

**Theorem 6₂.** For k = 6 and every weight vector with at most two distinct weights, (R∞)_j ⪰ 0 for j = −1, …, 6.
*Proof.* Every torus block is covered:
- (R∞)_{−1}, (R∞)_0, (R∞)_1, (R∞)_2 and (R∞)_6 as in the proof of Theorem 5₂ (Proposition M, Theorem B0, the
  Corollary of Section 3h of the README, the Corollary of Theorem L in `code/middle_j2_l1/` for 5 ≤ k ≤ 10, Theorem K);
- (R∞)₃: Theorem B (a), (b); |T| = 4 by Theorem DJ; type II by Proposition T2′;
- (R∞)₄: Theorem B (c), (d); |T| = 5 by Theorem DJ; type II by Proposition T2′;
- (R∞)₅ = (R∞)_{k−1}: Theorem T (letter-1 block); |T| = 1 by Theorems EU and TL; |T| = 2, 3 by Theorem B (e);
  |T| = 4, 5, 6 by Theorem U2, Proposition U1 and Proposition U0; type II by Proposition T2.
The kernels of (R∞)₃ and (R∞)₄ follow as in the proof of Theorem 5₂. ∎

**Corollary 6₂.** Conjecture 12 of `single-excitation-range/` holds for every p ∈ H₁ supported on at most six copies
whose nonzero weights take at most two distinct values. *Proof.* With exactly six nonzero weights, Theorem 6₂, the
remark that (R2) always holds and Lemma 13 there; with at most five, Corollary 5₂. ∎

**Checks at k = 6** (final run; all exact):
- `check_builders.py`: every certified block (both methods save the block over Q(q)) equals the Lemma B block at
  three rational squares q per configuration (228 comparisons) and the `ylib.py` point block at q = 9/16 (76); the
  builders agree at two general weight vectors for every family, with B scaled by 1 and by 1001/1000; homogeneity.
- `check_minors.py` (direct families): every certified minor (1 872) equals the direct minor of the gauge-fixed Lemma B
  block at two values of q.
- `check_sym.py` (reduced families): at two values of q, every certified part F_λ (114 parts) equals den · B_λỸB_λᵀ
  with Ỹ the Lemma B block, every certified minor (686) equals the direct minor of that matrix, the basis vectors lie
  in W, and every part is positive definite.
- `check_end2end.py`: at five values of q per representative configuration (10⁻⁶ to 10⁶; 380 points), the known
  vectors are annihilated, the nullity equals their number, and the gauge-fixed block is positive definite in two
  gauges.
- `coverage.py`: every one of the 8 576 pairs (block, level pattern) of the eleven families reduces to a certified
  configuration (directly or through the level exchange); the exact identity of Section 2 is checked for every 40th
  pair (a sample: these Lemma B builds are slow at k = 6).
- Negative controls: with B → (1 + 1/1000)B every certificate of a two-level representative configuration fails
  (65 of 65, `neg_cert_k6_*`, expected FAIL), and `neg_points.py` gives an explicit rational point for each (q = 10⁻⁶
  for 58 of them, 10⁴ for 7) where the scaled gauge-fixed block built by `j1lib.py` is not positive definite while the
  true one is.

## 8. What is not covered

- Three or more distinct weights, and general weights, at k = 5 and 6: the blocks of Theorems A and B, the letter-1
  and type I blocks of (R∞)₂ with a letter 1, the |T| = 1 and letter-1 blocks of (R∞)₁ (Conjectures L and L1), and the
  type I blocks of (R∞)_{k−1} with |T| = 1, 2 (k = 5) or |T| = 1, 2, 3 (k = 6). With more weight values the
  certificates would have more variables; this was not attempted.
- k ≥ 7 at two weight values. At k = 7 the open blocks would be those of (R∞)₃, (R∞)₄ and (R∞)₅ with a letter 1 and
  the type I blocks of (R∞)₆ with |T| = 2, 3, 4 (eighteen families, about 150 configurations up to the level
  exchange; (R∞)₂ is covered for k ≤ 10 by `code/middle_j2_l1/`). The largest is the letter-1 block of (R∞)₃ with
  7·35 + 35 = 280 states. One timing probe was run (`logs/probe_k7/`): for that block in configuration (0, 3)
  (three copies of squared weight 1, four of weight q), `cert_sym.py j3m0 0 3 k=7` passes in 15 s (W of dimension 245,
  11 parts, every leading minor with only nonnegative coefficients). This is one configuration of one family, without
  the checks of Sections 5 and 7, and it is not claimed as a result. By this probe the certificates for all of k = 7
  (about 150 configurations up to the level exchange) would take well under an hour; the Lemma B checks at k = 7
  would take longer.
- Nothing here bears on Conjectures L, L1, P or SF of the README as such: the certificates prove the block statements
  at two weight values, not the structural properties behind those conjectures.

## Files

| File | Content |
|---|---|
| `k5lib.py` | block states, known kernel vectors and gauge, configurations, the builders (`ylib.py` over any ring, `j1lib.py` at rational weights, rescaled), exact rank and LDL over Q |
| `cert5.py` | the certificate for one family and one configuration (Section 3); with `bscale` ≠ 1 the negative control |
| `signtest.py` | the sign test of the certificates (coefficient signs; exact witness search; Sturm counts from `ufield.py`) |
| `explore_kernel.py` | exact kernel and positivity check at a few points (run first, before the certificates) |
| `check_builders.py`, `check_minors.py`, `check_end2end.py`, `coverage.py`, `neg_points.py` | the checks of Sections 5 and 7 (`k=6`, `reps=1`: the representative configurations at k = 6) |
| `saved.py` | loads the blocks and minors saved by the certificates |
| `symred.py`, `cert_sym.py`, `check_sym.py` | isotypic reduction by Young symmetrizers, the certificate after it, and its checks (Section 7) |
| `dev_k6.sh`, `dev_k6_sym.sh` | development drivers for the k = 6 certificates (superseded by `run_all.sh`) |
| `j1lib.py`, `ylib.py`, `flintfield.py` | copied unchanged from `relative-anchors/code/middle_j2/` |
| `ufield.py` | copied unchanged from `code/middle_j2_l1/` (exact univariate rational functions, Sturm counts, sign test) |
| `guard.sh` | copied unchanged from `relative-anchors/code/no_letter1/` (one Python process, killed above 2.5 GB) |
| `run_all.sh` | reruns everything; `bash run_all.sh > logs/run_all.out` |
| `crosscheck.py`, `run_crosscheck.sh` | the second check of Section 9 (`bash run_crosscheck.sh > logs/run_crosscheck.out`, after `run_all.sh`, whose saved blocks it reads) |
| `run_extra.sh` | a second, independent certificate (isotypic reduction) for every directly certified configuration, and denser coverage samples at k = 6; not needed for any proof; `bash run_extra.sh > logs/run_extra.out` |
| `logs/run_all.out`, `logs/run_*.log`, `logs/cert_k5_*.pkl`, `logs/cert_k6_*.pkl` | the final run (the pickles hold the certified blocks, minors and reduced parts; about 89 MB, regenerated by `run_all.sh`; not in the repository) |
| `logs/run_crosscheck.out`, `logs/run_crosscheck_*.log` | the second check |
| `logs/extra/`, `logs/run_extra.out` | the run of `run_extra.sh` |
| `logs/probe_k7/` | the k = 7 timing probe of Section 8 |
| `logs/dev/` | development runs; `logs/dev/stopped_run/` is an earlier full run, stopped during the k = 6 negative controls (one of them spent 200 s in a Sturm count before `signtest.py` looked for witnesses first) |

## 9. A second check

`crosscheck.py` was written separately and uses none of `k5lib.py`, `cert5.py`, `cert_sym.py`, `symred.py`,
`signtest.py`, `saved.py` or the check scripts: only the Lemma B builder `j1lib.py`, python-flint and numpy. The torus
blocks, their known kernel vectors and the gauge are set up there from the definitions.

- *An exact test of positive definiteness.* For a rational symmetric matrix A let A_int be A times the lcm of its
  denominators. A floating-point Cholesky factor L of A gives Z, a rounding of L⁻¹ to integers row by row (one power of
  two per row); Z is lower triangular with nonzero diagonal, so B = Z A_int Zᵀ is congruent to A_int. If B, rescaled by
  these powers of two, is strictly diagonally dominant with positive diagonal (exact integer arithmetic), then A is
  positive definite. Floating point only chooses Z. Negative controls: the test refuses the singular letter-1 block of
  (R∞)₃ at k = 5, the gauge-fixed block minus 1.01 λ_min vvᵀ/|v|² (vᵀAv < 0 exactly), the block with one kernel
  vector left out of the gauge, and the gauge-fixed block with B → (1 + 1/1000)B at two-value points where it is not
  positive definite.
- *Every block of every condition, from Lemma B.* At 12 random two-value points at k = 5 and 24 at k = 6 (levels at
  random positions, weights from a list of rationals, q from 1/81 to 225/4), all states of G1_j ⊕ G0_{j+1} are
  enumerated for every j = −1, …, k and grouped by the number of letters 2 at each copy (row and column together),
  which the torus symmetry conserves; 300 random entries between different groups are checked to vanish. In every
  group the known kernel vectors (p ⊗ y for the columns y of its letter-1 rows; vac ⊗ J⁻|1…1⟩ in the letter-1 group of
  (R∞)_{k−2}; u ⊗ J⁻|1…1⟩ with u ⊥ p in the letter-1 group of (R∞)_{k−1}) are annihilated exactly, and the group
  without a set of pivot coordinates of these vectors (an exact row reduction picks them) is positive definite by the
  test above. So every group is ⪰ 0 with exactly that kernel, and the kernel dimensions add up to |H_j|, plus 1 at
  j = k − 2 and k − 1 at j = k − 1. This checks Theorems 5₂ and 6₂ at these points directly, including the blocks
  covered by the results cited in Sections 4 and 7.
- *The saved blocks.* For all 128 certificates, the rescaled block saved over Q(q) is evaluated at 31 values of q from
  10⁻⁸ to 10⁸; the column-group indicators (and, for the letter-1 block of (R∞)_{k−2}, the vector (s_h) on the states
  (vac, 1…1 − e_h)) are annihilated, and the block without a set of pivot coordinates of them is positive definite
  (exact elimination over Q).

`run_crosscheck.sh` runs all of this, one Python process at a time; the output is `logs/run_crosscheck.out`.
