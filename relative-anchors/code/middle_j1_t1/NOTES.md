# (R∞)₁, the type I blocks with one copy carrying the letter 2: equal U-weights and two U-weight levels for every k

Ansh Mishra and Aryan Senthilkumar, 5 November 2026. Nothing here is refereed.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Lemma P (P-first reduction).** For every k ≥ 2, every copy t and every positive weight vector, the type I block of (R∞)₁ with T = {t} (size 3k − 1) is ⪰ 0 with kernel exactly p ⊗ \|2@t⟩ if and only if the explicit matrix Φ = C − Bᵀ diag(1/N) B of size 2k − 1 is positive definite. Φ is written out in closed form (Section 1) | **proved** (by hand, from Lemma T1.1 of `middle_j1`); **exact check** of the reduction and of the closed form of Φ at 18 rational weight vectors, k = 3, …, 8 (`check_phi_closed.py`), and of the reduction at the 32 weight vectors of `check_lev.py` |
| 2 | **Theorem EQ (equal U-weights).** For every k ≥ 3, every t and every weight vector whose weights on U = [k] ∖ {t} are all equal (the weight σ of t arbitrary), the block is ⪰ 0 with kernel exactly p ⊗ \|2@t⟩ | **proved**, computer-assisted, exact, symbolic in u = k − 1 and in the two weights (Section 2); two independent certificates (through Φ, and through the gauge-fixed block) |
| 3 | **Theorem 2L (two U-weight levels).** The same for every k ≥ 3 and every weight vector whose U-weights take at most two distinct values (u₁ copies of one value, u₂ of the other, any u₁, u₂, any σ) | **proved**, computer-assisted, exact, symbolic in u₁, u₂ and the three weights (Section 3) |
| 4 | k = 3 and k = 4, all positive weights: a third proof of these cases (already proved in `middle_j1`), through Lemma P and the inverse-free congruent form Φ′ of Φ | **proved**, computer-assisted, exact (`cert_lev.py F 2 1 1 fixed=1,1`, `cert_pivots.py 3`, `cert_pivots.py 4`, `bareiss_phic.py 4`; Section 4) |
| 5 | k = 5, all positive weights (four free parameters), through Φ′ (size 9 instead of 13) | **not completed**: the leading-minor cores D₁, D₂ have only positive coefficients; the run passed 2.5 GB at the next step and was stopped (Section 4) |
| 5b | Two U-copies of arbitrary weight and u₃ ≥ 2 equal ones (every k ≥ 5) | **not completed**: the first three leading minors of the 7 × 7 invariant form are certified; the elimination then crashed in the polynomial library (three runs) or was stopped (Section 3b) |
| 6 | The symbolic level reductions equal the exact projections of the Lemma B block (`j1lib.py`): 64 comparisons at 32 rational weight vectors, each with B scaled by 1 and by 1 + 1/1000, for Φ and for the gauge-fixed block (one level, u = 2, 3, 4, 6; two levels, (u₁, u₂) = (1, 2), (2, 2), (2, 3), (3, 1), (1, 4), and with fixed sizes (1, 2), (1, 3), (1, 1); three levels (1, 1, u₃), u₃ = 2, 3, 5; four single copies, k = 5); `j1lib.py` equals `ylib.Builder` on six \|T\| = 1 blocks, k = 3, …, 6; the closed forms of Lemma T1.3 equal the Lemma B block at 18 points; the two implementations of the equal-U reduction agree | **exact check** (`check_lev.py`, `check_j1_ylib.py`, `check_t1_closed.py`, `cmp_impl.py`) |
| 7 | Negative controls: with B → (1 + 1/1000)B the certificate of Theorem EQ fails, and at explicit rational weight vectors (one level, u = 2; two levels, (2, 2) and (1, 3); three levels, (1, 1, 2)) the scaled gauge-fixed block is not positive definite while the true one is (both from the symbolic forms and from `j1lib.py`) | **exact** (`cert_lev.py G 1 2 1001/1000`, `neg_lev.py`) |
| 8 | For general U-weights the diagonally normalised Φ has its smallest eigenvalue far below the equal-U value (for example 0.07 against 0.61 at k = 5, σ = 0.01), at U-weights concentrated on one or two copies, and down to about 10⁻⁶ at extreme weights; so equal U-weights is not the worst case, a reduction of general weights to Theorem EQ by a monotonicity argument cannot work, and lossy bounds cannot work near the degenerate limits | **numerical**, 40-digit arithmetic (`explore_gap.py`, `explore_phi.py`); no claim rests on it |
| 9 | The block for k ≥ 5 and general weights; three or more U-weight levels | **open** |

Nothing here says anything about Conjecture 12 of `single-excitation-range/` in general, or about the NPT (Werner
distillability) problem. What is new: Lemma P (every k, all weights, a reduction), and the type I blocks of (R∞)₁ with
|T| = 1 for every k when the U-weights take at most two values (Theorems EQ and 2L), with arbitrary weight on the
copy t. For general weights and k ≥ 5 the block remains open.

Notation as in `relative-anchors/README.md` (Sections 1, 2) and `middle_j1/NOTES.md` (Section 3). Forms are in units
2^{k−1}. Computer-assisted means: exact rational-function arithmetic in python-flint, and sign checks of coefficients;
no floating point enters any proof.

## 1. Setting and the P-first reduction

One copy t carries the letter 2, σ = s_t, U = [k] ∖ {t}, u = |U| = k − 1. The block has 3k − 1 states:
- R_c = (2@t, 1@c), c ∈ [k];
- P_r = (1@r, 2@t), r ∈ [k];
- Q_c = (vac, 2@t 1@c), c ∈ U.

N_r = 2 for r ∈ U and N_t = 1. Order the states P first, then R and Q, and write the block as
M = [[Π, B], [Bᵀ, C]], with C the block on the R- and Q-states and B the coupling of the P-states to them.

**Lemma P.** For every k ≥ 2, every t and every positive weight vector, put

    Φ = C − Bᵀ diag(1/N) B      (size 2k − 1).

Then M ⪰ 0 with kernel exactly p ⊗ |2@t⟩ if and only if Φ ≻ 0.

*Proof.* By Lemma T1.1 of `middle_j1/NOTES.md`, Π = diag(N) − (N∘w)(N∘w)ᵀ/(2|s| − σ), and (w, 0), which is
p ⊗ |2@t⟩, lies in the kernel of M; so Πw = 0 and Bᵀw = 0. Since Σ_r N_r s_r = 2|s| − σ, every b ⊥ w satisfies

    Π diag(1/N) b = b − (N∘w)(wᵀb)/(2|s| − σ) = b,

so Π diag(1/N) B = B. With K = [[I, −diag(1/N)B], [0, I]] this gives KᵀMK = diag(Π, Φ). Further
Π = N^{1/2}(I − vvᵀ)N^{1/2} with v = N^{1/2}w/|N^{1/2}w|, so Π ⪰ 0 with kernel span(w). Hence M ⪰ 0 iff Φ ⪰ 0, and
ker M = K(span(w) ⊕ ker Φ), which is span((w, 0)) iff Φ ≻ 0. ∎

Remarks.
- Φ is also the Schur complement, onto the R- and Q-states, of the gauge-fixed block in which the state P_t is
  deleted: shifting π along w to π_t = 0 does not change the form.
- With the closed forms of Lemma T1.3 (|s| = 1) the products are explicit:

      Φ_RR = D_R − ¼ Λ (I − wwᵀ) N⁻¹ (I − wwᵀ) Λ,
      Φ_RQ = −(I + ½N₁)⁻¹H + ¼ τ′(2 − σ) Λ (E − w w_Uᵀ) D′,
      Φ_QQ = 2(I + ½N′)⁻¹ − τ′²(2 − σ) D′ (½(2 − σ) I − w_U w_Uᵀ) D′,

  with Λ = diag(τ_c/(1 − x_c)) (τ_c = τ on U and τ − 2 at t), D′ = diag_{c∈U} 1/(2 − σ − 2x_c), E the embedding
  R^U → R^[k], and τ, τ′, N₁, N′, H as in Lemma T1.3. They follow from Y_RP = ½Λ(I − wwᵀ),
  Y_PQ = τ′(N w w_Uᵀ − (2 − σ)E)D′, wᵀNw = 2 − σ and EᵀN⁻¹E = ½I.
- The congruence by diag(I + ½N₁, I + ½N′) (both positive definite) removes the two inverses:

      Φ′_RR = A₁D_RA₁ − ¼A₁Λ(I − wwᵀ)N⁻¹(I − wwᵀ)ΛA₁,   Φ′_RQ = −HA₂ − A₁Y_RPN⁻¹Y_PQA₂,
      Φ′_QQ = 2A₂ − A₂Y_PQᵀN⁻¹Y_PQA₂,   A₁ = I + ½N₁, A₂ = I + ½N′,

  and Φ ≻ 0 iff Φ′ ≻ 0. Its entries contain pair sums (through N₁, N′) and sums over a third copy.
- Compared with Conjecture L of `middle_j1` (eliminate R and Q, then a Laplacian on P), Lemma P eliminates P with an
  explicit inverse and leaves a problem without kernel. By itself it proves nothing for general weights; row 8 of the
  table shows that Φ comes close to singular in several limits.

`check_phi_closed.py` checks Π diag(1/N) B = B, wᵀB = 0, the closed form of Φ and Φ ≻ 0 exactly against the Lemma B
builder `j1lib.py` at 18 rational weight vectors with |w| = 1, k = 3, …, 8 (three per k, t at different positions).

## 2. Equal U-weights: Theorem EQ

**Theorem EQ.** Let k ≥ 3, t ∈ [k], and let w_h = ω for all h ∈ U, with ω > 0 and w_t > 0 arbitrary. Then the type I
block of (R∞)₁ with T = {t} is ⪰ 0 with kernel exactly p ⊗ |2@t⟩.

*Proof (computer-assisted, exact, symbolic in u; `lev.py`, `cert_lev.py`).*
1. **Coordinates.** Write p = σ and q = ω² for the unnormalised squared weights, L = p + uq; the closed forms of
   Lemma T1.3 are used with x_t = p/L and x_h = q/L (the unrescaled block is homogeneous of degree 0 in the
   weights). Work in the rescaled frame D = diag(μ), μ = w_c on R_c, P_c and Q_c, as in `ylib`: for an operator O
   between the coordinate spaces R^[k] and R^U put Õ = DOD. Then (O₁O₂)~ = Õ₁D⁻²Õ₂, (O⁻¹)~ = D²Õ⁻¹D² and
   (Oᵀ)~ = Õᵀ, and every Õ that occurs has entries rational in p, q.
2. **Symmetry.** The permutations of U preserve the weights, so they act on the block (by permuting states) and on
   Φ, and D commutes with them. As a representation, the coordinates (ρ, q) of Φ form R^[k] ⊕ R^U = 3 trivial + 2
   standard (u ≥ 2). By Schur's lemma Φ ≻ 0 iff its invariant part F_inv (3 × 3, on the vectors e_t, 1_U in R^[k]
   and 1_U in R^U) and its standard part F_std (2 × 2, on (f, f) with Σ_h f_h = 0) are positive definite.
3. **The reduced forms, symbolically in u.** An invariant operator is stored by its matrix in the bases (e_t, 1_U)
   and (1_U) and by its scalar on the standard part: a k × k invariant O with entries α = O[t,t], β = O[t,h],
   β′ = O[h,t], δ = O[h,h], ε = O[h,h′] (h ≠ h′ in U) is [[α, uβ], [β′, δ + (u − 1)ε]] and δ − ε. Products,
   transposes (Gram matrices diag(1, u) and (u)) and inverses are then 2 × 2 computations. `lev.py` evaluates the
   closed forms of Lemma T1.3 and Lemma P in this algebra, with u = 2 + x, over Q(p, q, x).
4. **Validation.** At 8 rational points (u = 2, 3, 4, 6) F_inv and F_std equal, exactly, the invariant and standard
   projections of DΦD computed from the Lemma B block of `j1lib.py` (`check_lev.py`; the same for the gauge-fixed
   block of step 6, with B scaled by 1 and by 1 + 1/1000).
5. **Certificate 1** (`cert_lev.py F 1 2`). Every leading principal minor of F_inv and of F_std is a constant times a
   product of irreducible polynomials in (p, q, x) with only positive coefficients, each with a monomial free of x.
   So every leading minor is positive for p, q > 0 and every real u ≥ 2, and Φ ≻ 0 (Sylvester). The determinant of
   F_inv has numerator degree 27 (core factor: degree 21, 76 terms); that of F_std degree 24 (core: degree 21,
   71 terms). With Lemma P this proves the theorem.
6. **Certificate 2, independent of Lemma P** (`cert_lev.py G 1 2`). Delete the state P_t. Since p ⊗ |2@t⟩ is in the
   kernel and has P_t-coordinate w_t ≠ 0, the block is ⪰ 0 with kernel exactly p ⊗ |2@t⟩ iff this gauge-fixed block
   is positive definite. Its invariant part (4 × 4, on (ρ_t, ρ_U, π_U, q_U)) and standard part (3 × 3) have all
   leading minors of the same form as in step 5 (the same cores of degree 21 appear in the last minors). ∎

Remarks.
- The README proves the blocks at equal weights for all copies (w_t = ω, Section 4 there). Theorem EQ allows any
  weight on t. At equal U-weights the smallest eigenvalue of the diagonally normalised Φ is about 0.6 for small
  σ/|s| at k = 5, 8, 12, is of order σ/|s| at k = 3, and tends to 0 as σ → |s| (`explore_phi.py`, floating point).
- An earlier implementation of the same reduction written separately for this case (`equ_sym.py`, with its own
  certificate `cert_equ.py`) gives the same minors (`cmp_impl.py` compares the two exactly).

**Negative control.** In every entry, multiply the B-part (the terms with an E⁻¹) by 1001/1000. The A-parts read off
from Lemma B are: all of D_R; diag(N) in Π; δ_cr(½ − [r = t]) in Y_RP; −δ in Y_PQ and in Y_RQ; 2I in Y_QQ. The
scaled closed forms equal the `j1lib.py` block with `bscale = 1001/1000` (`check_lev.py`). With them:
- the certificate of step 6 fails (`cert_lev.py G 1 2 1001/1000`: the third invariant minor has a factor with mixed
  signs);
- at u = 2, w_t = 1, ω = 1/100 a leading minor of the scaled forms is negative while all true ones are positive, and
  the scaled gauge-fixed block from `j1lib.py` is not positive definite while the true one is (`neg_lev.py`).

(For the block itself B → (1 + 1/1000)B always fails trivially, on the kernel vector p ⊗ |2@t⟩; the controls above
are on the gauge-fixed block, where the kernel is removed.)

## 3. Two U-weight levels: Theorem 2L

**Theorem 2L.** Let k ≥ 3, t ∈ [k], U = U₁ ⊔ U₂ with |U_i| = u_i ≥ 0, and let the squared weights be q_i on U_i
(q₁, q₂ > 0 arbitrary) and σ on t (arbitrary). Then the type I block of (R∞)₁ with T = {t} is ⪰ 0 with kernel exactly
p ⊗ |2@t⟩.

*Proof (computer-assisted, exact, symbolic in u₁, u₂; `lev.py`, `cert_lev.py`).* If u₁ = 0 or u₂ = 0 this is Theorem
EQ. Otherwise the symmetry is S_{u₁} × S_{u₂}; the invariant vectors are e_t, 1_{U₁}, 1_{U₂} in R^[k] and 1_{U₁},
1_{U₂} in R^U, and each level with u_i ≥ 2 contributes a standard part. An invariant k × k operator is stored by its
3 × 3 matrix in the basis (e_t, 1_{U₁}, 1_{U₂}) and two scalars; entries between different levels count with the
size of the other level (Gram matrix diag(1, u₁, u₂)). By Lemma P it suffices that the invariant part F_inv of DΦD
(5 × 5, on (ρ_t, ρ₁, ρ₂, q₁, q₂)) and the standard parts F_std,i (2 × 2) for the levels with u_i ≥ 2 are positive
definite. The reduced forms equal the exact projections from `j1lib.py` at 16 two-level weight vectors (10 with
symbolic sizes, 6 with u₁ = 1 or u₁ = u₂ = 1 fixed; `check_lev.py`). Certificates (each: every leading principal
minor is a ratio of polynomials, or of products of irreducible polynomials, with one-signed coefficients, positive,
each with a monomial free of the x_i):
- **u₁, u₂ ≥ 2** (u_i = 2 + x_i): F_inv (`cert_lev.py F 2 2 2 only=inv`; the numerator of the determinant has degree
  75 and 96 946 terms, all positive; factored, its core has degree 66 and 66 038 terms; about 6 minutes); F_std,1 for
  u₁ ≥ 2, u₂ ≥ 1 (`cert_lev.py F 2 2 1 only=std`) and F_std,2 for u₁ ≥ 1, u₂ ≥ 2 (`cert_lev.py F 2 1 2 only=std`).
- **u₁ = 1, u₂ ≥ 2**: F_inv with u₁ = 1 exactly (`cert_lev.py F 2 1 2 fixed=1,0 only=inv`; last minor numerator of
  degree 67, 6 427 terms, all positive) and F_std,2 as above. With u₁ treated as a real variable in [1, 2) the fourth
  minor has a factor with mixed signs, so the integer case u₁ = 1 is done separately.
- **u₁ ≥ 2, u₂ = 1**: this is the previous case with the names of the two levels exchanged (the set of weight vectors
  covered is the same).
- **u₁ = u₂ = 1** (k = 3, all weights): F_inv is the whole Φ (`cert_lev.py F 2 1 1 fixed=1,1`; last minor degree 32,
  373 terms, all positive). ∎

The case "one copy of U with an arbitrary weight, all other U-weights equal" is contained in Theorem 2L (u₁ = 1); this
includes one small U-weight, the regime where the exact checks of `middle_j1` find the smallest Laplacian margins.

## 3b. Two U-copies of arbitrary weight, the others equal: attempted, not completed

The next family would be U = {a, b} ⊔ U₃ with squared weights q₁ at a, q₂ at b and q₃ on all of U₃, |U₃| = u₃ ≥ 2
(k ≥ 5), σ arbitrary. It contains every concentrated limit of Section 5 with at most two dominant U-copies, and with
k = 3, 4 and Theorems EQ, 2L it would give the block for every k whenever all but at most two U-weights are equal.
By Lemma P and the symmetry S_{u₃} it reduces to a 7 × 7 invariant form (ρ_t, ρ_a, ρ_b, ρ₃, q_a, q_b, q₃) and a
2 × 2 standard form, symbolic in (q₁, q₂, q₃, x₃) with p = 1 and u₃ = 2 + x₃. The reduced forms are validated against
`j1lib.py` at u₃ = 2, 3, 5 (`check_lev.py`). What we have:
- Fraction-free elimination of den·F_inv (`bareiss_phic.py levels 3 1,1,2 1,1,0 inv F`): den has five irreducible
  factors, all with positive coefficients; the cores of D₁ (degree 16, 493 terms), D₂ (degree 33, 4 843 terms) and
  D₃ (degree 50, 20 473 terms) have only positive coefficients. Three runs (one directly, two in a thread with a
  128 MB stack through `bigstack.py`, copied from `middle_j1`) ended with a segmentation fault inside the polynomial
  library while computing D₄, at about 1.4 GB (`logs/dev/bareiss_l3_112_inv_F_segfault.log`). The same elimination
  on Φ′ instead of Φ has larger entries (D₁: 2 988 terms) and was stopped by hand.
- Pivots over the rational-function field (`pivots_lev.py 3 1,1,2 1,1,0 inv`): pivots 1 and 2 are ratios of
  polynomials with only positive coefficients; stopped by hand during pivot 3 to leave time for the final run
  (`logs/dev/pivots_l3_112_inv_stopped.log`).

So the first three leading minors of the invariant part are certified, which is not a proof of anything about the
block. `run_all.sh` starts both attempts again with LONG=1.

## 4. Fixed k, all weights: k = 3, 4 again, and k = 5

Through Lemma P, the statement at a fixed k for all weights is: Φ′ (size 2k − 1) ≻ 0 on the open orthant of squared
weights. In `lev.py` this is the case of k − 1 levels of size one; the rescaled entries are homogeneous of degree 0 in
(p, q_1, …, q_{k−1}), so p = 1.
- k = 3: `cert_lev.py F 2 1 1 fixed=1,1` (Φ) and `cert_pivots.py 3` (Φ′): all pivots are ratios of polynomials with
  only positive coefficients.
- k = 4: `cert_pivots.py 4` (pivots of Φ′ over Q(q₁, q₂, q₃), about 4 minutes; the last pivot has numerator degree 57,
  32 316 terms) and `bareiss_phic.py 4` (fraction-free elimination with the denominator factors pulled out, about
  40 seconds; the core of D₇ has degree 62, 32 316 terms). All cores have only positive coefficients. Both scripts
  also compare the symbolic Φ′ with Φ′ computed from `j1lib.py` and the closed A₁, A₂ (exactly, at two points), and
  the minors with direct determinants at a point.
- k = 5: `bareiss_phic.py 5` (Φ′ of size 9 over Q(q₁, …, q₄)); not completed, see below.

These are further independent proofs of k = 3 and k = 4 (besides Conjecture L at k = 3, 4 in `middle_j1` and Theorem
A4 of the README). The matrix is smaller (2k − 1 instead of 3k − 2), but the certificate is not: at k = 4 the largest
core has degree 62 and 32 316 terms, against degree 62 and 27 536 terms for the gauge-fixed block in Section 3f of the
README.

**k = 5: not completed.** `bareiss_phic.py 5` builds Φ′ (size 9, over Q(q₁, …, q₄)) in 4 seconds; the common
denominator has only factors with positive coefficients (linear forms, and two factors of degrees 3 and 5); the cores
of D₁ (degree 18, 7 101 terms) and D₂ (degree 35, 73 614 terms) have only positive coefficients. During the next
Bareiss step the process passed 2.5 GB and was stopped by `guard.sh` (`logs/bareiss_phic_k5_stopped.log`, after
about 3 minutes). From k = 3 to k = 4 the largest core grew from 401 to 32 316 terms; at k = 5 the cores of D₇, D₈, D₉
would plausibly have 10⁶ terms or more, so this direct certificate needs either much more memory or a smarter
elimination (for example evaluation and interpolation of the cores, or a block elimination that keeps Φ_RR, which has
no pair sums, separate). In the three-level family (Section 3b) the entries of Φ itself were several times smaller than
those of Φ′ (first core 493 against 2 988 terms); the same elimination on Φ at k = 5 was not tried for lack of time.
It is not part of `run_all.sh` (LONG=1 starts it).

## 5. Observations for general weights (floating point; no claim rests on them)

`explore_gap.py` (40-digit arithmetic) computes the smallest eigenvalue of Dg^{−1/2}ΦDg^{−1/2}, Dg = diag(Φ), at
fixed σ = s_t/|s| over random U-weights:

| k | σ | equal U-weights | smallest found | where |
|---|---|---|---|---|
| 4 | 0.01 | 0.619 | 0.059 | two U-copies with half each, one tiny |
| 5 | 0.01 | 0.606 | 0.072 | two U-copies dominant, two tiny |
| 7 | 0.01 | 0.599 | 0.252 | two or three U-copies dominant |
| 4, 5 | 0.3 | 0.50 | 0.33 | one U-copy dominant |
| 4, 5, 7 | 0.9 | 0.13 to 0.15 | 0.035 to 0.040 | one U-copy dominant |

So the worst U-weights are concentrated, not equal. `explore_phi.py` (also 40 digits) shows the two hard limits:
- σ → 0 with essentially two U-copies: at k = 3 and equal U-weights the normalised gap is 6.2·10⁻⁴ at σ = 10⁻⁴ and
  5.8·10⁻² at σ = 10⁻²; among random samples, 4.3·10⁻⁶ at k = 4 with σ ≈ 6·10⁻⁷, two U-weights near ½ and one tiny;
- σ → |s| (the copy t carries almost all the weight): 7.6·10⁻⁷ at k = 3, about 10⁻⁵ at k = 5, 6; at equal U-weights
  the gap is about 2·10⁻³ at σ = 0.999 for every k tested (3, 5, 8, 12).

At equal U-weights and k ≥ 5 the gap stays near 0.6 as σ → 0 (k = 5, 8, 12). A proof for general weights must be exact
in these limits; lossy Gershgorin or Schur bounds of the kind used in Theorem C of the README cannot see a gap of order
σ or 1 − σ. The concentrated limits are covered by Theorem 2L only when the remaining U-weights are equal.

(An earlier double-precision version of `explore_phi.py` reported a negative eigenvalue, −10⁻⁴, at k = 3 near
w = (1, 0.0123, 0.0001); this is a cancellation artefact of the closed form: Φ ≻ 0 for all weights at k = 3
(Section 4).)

## 6. What did not work, and what is left

**One pattern behind the certificates (not a mechanism).** Every certificate above has the same shape: in the rescaled
frame, the leading principal minors of Φ (or Φ′, or the gauge-fixed block) are ratios of polynomials with only
positive coefficients in the squared weights, and, in the level families, also in the shifted level sizes
x_i = u_i − 2. This is Conjecture SF of the README (Section 3g) seen through Lemma P: Φ is a Schur complement of the gauge-fixed
block (Section 1, first remark), and Schur complements preserve the property (Lemma 1 of Section 3g there), so
Conjecture SF for the gauge-fixed |T| = 1 block of (R∞)₁ would imply that every principal minor of DΦD is
subtraction-free, hence Φ ≻ 0 for every k. We did not find a reason for the pattern; the level computations are
consistent with it but do not explain it.

- With u₁ ∈ [1, 2) as a real parameter the two-level minors are not one-signed (Section 3); the integer case u₁ = 1 is
  separate. The same happens with the coordinate order (ρ_t, ρ₁, q₁, ρ₂, q₂) and offsets (1, 1).
- The sign pattern of Φ is not that of a Z-matrix up to signs: in the rescaled frame Φ_RR[t, c] < 0 and
  Φ_RR[c, c′] > 0 for c ≠ c′ in U, which no signature makes nonpositive. So the Z-matrix certificate of Theorems EU
  and TL of the README does not transfer to Φ; leading minors were used instead.
- The direct certificates grow fast: k = 5 for all weights (Section 4) and the three-level family (Section 3b) both
  ran into the memory limit or into crashes of the polynomial library. A smarter elimination (evaluation and
  interpolation of the cores, or eliminating the R-states first through the explicit form of Φ_RR, which is a diagonal
  matrix plus a matrix of rank two and contains no pair sums) is the obvious next step for both.
- Open: the block for general weights at every k ≥ 5; three or more U-weight levels, in particular the family of
  Section 3b, which contains every concentrated limit of Section 5 with at most two dominant U-copies.

## Files

- `j1lib.py`, `ylib.py`, `flintfield.py`, `check_t1_closed.py`: copies from `middle_j1` (unchanged), the Lemma B
  builder, the exact builder of `relative-anchors/code`, the rational-function field on python-flint, and the check of
  Lemma T1.3.
- `tlib.py`: the closed forms of Lemma T1.3 for any number type, and Φ of Lemma P.
- `lev.py`: the level reductions (symbolic): F_inv, F_std (Φ), G_inv, G_std (gauge-fixed block, B-part scaled by any
  factor), C_inv (the congruent form Φ′); `det` and the sign check used by the certificates.
- `check_j1_ylib.py`, `check_phi_closed.py`, `check_lev.py`: exact checks (rows 1 and 6).
- `cert_lev.py`: the certificates of Theorems EQ and 2L; `neg_lev.py`: negative controls.
- `cert_pivots.py`, `bareiss_phic.py`: fixed k, all weights (Section 4); `bareiss_phic.py levels …` and
  `pivots_lev.py`: the level families by elimination (Section 3b); `bigstack.py` (copied from `middle_j1`): runs a
  script in a thread with a 128 MB stack.
- `equ_sym.py`, `cert_equ.py`, `cmp_impl.py`: the separate first implementation for equal U-weights and its comparison
  with `lev.py`.
- `explore_phi.py`, `explore_gap.py`: floating-point explorations (Section 5).
- `guard.sh`: runs one Python process and kills it above 2.5 GB or after a time limit; `run_all.sh`: reruns every
  check, one Python process at a time, PASS or FAIL per step and ALL PASS at the end; its output is
  `logs/run_all.out`, the per-step logs `logs/run_*.log` (both from a complete rerun in a fresh copy of this folder,
  ALL PASS). `logs/explore_*.log`: the explorations;
  `logs/bareiss_phic_k5_stopped.log`: the stopped k = 5 run; `logs/dev/`: logs of the development runs quoted in the
  text (the same scripts before the final run, the crashed and stopped three-level runs).
