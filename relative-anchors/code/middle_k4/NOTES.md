# The middle conditions (R∞)₁ and (R∞)₂ at k = 4: exact certificates for the letter-1 and type I blocks

Ansh Mishra and Aryan Senthilkumar, 5 November 2026. Nothing here is refereed.

Notation as in `relative-anchors/README.md` (Sections 1–3e): k = 4 copies, weights w_t > 0, s_t = w_t²,
Y = A − B, (R∞)_j = Y compressed to G1_j ⊕ G0_{j+1}. A torus block is labelled by j and by the set T of copies that
carry the letter 2 exactly once (T = ∅ is the letter-1 block). Type II blocks (some copy carries the letter 2 in row
and column) are positive definite by Proposition T2′ and are not treated here.

## Status

Every torus block of (R∞)₁ and (R∞)₂ at k = 4 that is not of type II, one line per orbit under permutations of the
copies. "Known kernel" lists the vectors that the certificate shows to be the whole kernel of the block.

**Summary.** (R∞)₁ at k = 4 is proved for all positive weights (Theorem A). (R∞)₂ at k = 4 is proved on every
block except the orbit of the four type I blocks with |T| = 1 (Theorem B); that orbit is open. In particular the
letter-1 block of (R∞)₂ (the block of Section 3e of the README, with the soft modes of Corollary S) is proved at
k = 4. So Conjecture 12 for four copies is reduced to one orbit of blocks of size 21 (gauge-fixed 18).

| j | T | blocks in orbit | size | known kernel | gauge-fixed size | certificate | status |
|---|---|---|---|---|---|---|---|
| 1 | ∅ (letter-1) | 1 | 22 | p ⊗ \|1@c⟩, c = 0..3 (4) | 18 | (a) all 18 leading principal minors (natural order) have only positive coefficients after the denominator factors are pulled out; cores up to degree 93, 45 335 terms | **proved** (exact, computer-assisted; 10.4 min) |
| 1 | {t} | 4 | 11 | p ⊗ \|2@t⟩ (1) | 10 | (a), natural order; cores up to degree 62, 27 536 terms | **proved** (exact; 17 s) |
| 1 | {t, t′} | 6 | 3 | none | 3 | (a), natural order | **proved** (exact; < 1 s) |
| 2 | ∅ (letter-1) | 1 | 28 | p ⊗ \|1@P⟩, P a pair (6), and vac ⊗ J⁻\|1111⟩ (Proposition K2) | 21 | (a), natural order: all 21 leading principal minors have only positive coefficients after the denominator factors are pulled out; cores up to degree 121, 47 905 terms | **proved** (exact, computer-assisted; 10.7 min) |
| 2 | {t} | 4 | 21 | p ⊗ \|2@t, 1@c⟩, c ≠ t (3) | 18 | none: (a) not run (time budget, Section 5) | **open**; exact point checks only (gauge-fixed block ≻ 0 at 33 rational weight vectors, exact LDL) |
| 2 | {t, t′} | 6 | 12 | p ⊗ \|2@t, 2@t′⟩ (1) | 11 | (a), natural order; cores up to degree 79, 59 548 terms | **proved** (exact; 2.4 min) |
| 2 | {t, t′, t″} | 4 | 4 | none | 4 | (a), natural order | **proved** (exact; < 1 s) |

| Other statements | Status |
|---|---|
| Torus structure at k = 4: (R∞)₁ has 88 states in 15 torus blocks (1 letter-1, 10 type I, 4 type II), (R∞)₂ has 224 states in 31 blocks (1, 14, 16) | **exact** check at a rational point (`num_k4.py`); the block structure itself is the torus symmetry of the README, Section 1 |
| The blocks used equal the dense construction from M, V_p, D_p (definitions) | **numerical** (3 weight vectors, deviation ≤ 9·10⁻¹⁶, `check_direct.py`) |
| The symbolic blocks (s₀ = 1) equal the exact point builder, also at points with s₀ ≠ 1 after scaling | **exact** (6 points × 7 blocks, `check_direct.py`) |
| The certified leading minors equal the leading minors of the directly built block at random rational points | **exact** (4 points per certificate, `check_minors_points.py`) |
| Known kernel vectors annihilated identically; gauge-fixed blocks ≻ 0 at 33 rational weight vectors (spreads up to e^{±5} in s, near-degenerate points) | **exact** (`build_sym.py`, `num_k4.py`) |
| Negative controls: with B replaced by (1 + 1/1000)B every block fails, with an explicit rational witness vector | **exact** (`controls.py`) |
| Conjecture 12 for every p ∈ H₁ supported on at most four copies (all of (R∞) at k = 4) | **open**: it follows from what is proved here and in the README once the four blocks j = 2, T = {t} are certified (Section 5) |

Nothing here says anything about Conjecture 12 in general or about the endpoint (NPT distillability) problem beyond
the blocks listed.

## 1. The blocks at k = 4

A state is (row, column) with row = a single letter ℓ@r or the vacuum. The blocks that are not of type II are:

- **(R∞)₁, letter-1** (22 states): (1@r, 1@c) for r, c ∈ [4], and (vac, 1@{c, c′}).
- **(R∞)₁, T = {t}** (11): (2@t, 1@c) for every c (c = t allowed), (1@r, 2@t) for every r, (vac, 2@t 1@c), c ≠ t.
- **(R∞)₁, T = {t, t′}** (3): (2@t, 2@t′), (2@t′, 2@t), (vac, 2@t 2@t′).
- **(R∞)₂, letter-1** (28): (1@r, 1@P) for r ∈ [4] and P a pair, and (vac, 1@[4]∖h) for every hole h.
- **(R∞)₂, T = {t}** (21): (2@t, 1@P) for every pair P, (1@r, 2@t 1@c) with c ≠ t, (vac, 2@t 1@{c, c′}), c, c′ ≠ t.
- **(R∞)₂, T = {t, t′}** (12): (2@t, 2@t′ 1@c) and (2@t′, 2@t 1@c) with c ≠ t′ (resp. t), (1@r, 2@t 2@t′), and
  (vac, 2@t 2@t′ 1@c) with c ∉ {t, t′}.
- **(R∞)₂, T = {t, t′, t″}** (4): (2@t, 2 on the other two), three of them, and (vac, 2@T).

`num_k4.py` checks exactly, at a rational weight vector, that all entries of (R∞)_j between different torus labels
vanish and that these blocks, their images under permutations of the copies and the type II blocks exhaust
(R∞)₁ and (R∞)₂. (The block structure itself holds identically: Y commutes with the torus action, README Section 1.)

## 2. Reductions

**Permutations of the copies.** Let π permute the four copies, in row and column at once. It commutes with M^{⊗4} and
with the partial transpose, maps p_G(w) to p_G(w ∘ π⁻¹), preserves G1_j and G0_{j+1}, and permutes the torus labels.
So the block (j, T) at the weights w is, up to a permutation of its states, the block (j, π(T)) at the weights
w ∘ π⁻¹. All certificates below are symbolic in all four weights, so the representatives T = ∅, {0}, {0, 1},
{0, 1, 2} cover every block of the orbit.

**Rescaling and homogeneity.** As in `relative-anchors/code/ylib.py`, put x = Ψ/μ with μ_X = ∏_t w_t^{N_t(X)}, N_t(X)
the number of letters 1 at copy t in row and column. Then μYμ has entries in Q(s₀, …, s₃) (ylib asserts that every
exponent is even). On a block every state has the same total number n₁ = j + 1 − |T| of letters 1, so under w → λw
the rescaled block is multiplied by λ^{2n₁}: its entries are homogeneous of degree n₁ in s. Hence it suffices to take
s₀ = 1, and every computation below is over Q(s₁, s₂, s₃).

**Known kernel and gauge.** In rescaled coordinates:
- p ⊗ y (y a column string of the block) is the all-ones vector on the four states (1@r, y), r = 0..3 (its
  coefficient on (1@r, y) is w_r, and μ = w_r · ∏_{y_t = 1} w_t);
- vac ⊗ J⁻|1111⟩ (j = 2, letter-1 block) is x_(vac, hole h) = s_h / ∏_t w_t, zero elsewhere.

`build_sym.py` checks that these vectors are annihilated, as identities of rational functions. Delete the coordinate
(1@1, y) for every y, and for the letter-1 block of (R∞)₂ also (vac, hole 0). Each known kernel vector is the only one
that is nonzero on its deleted coordinate, so every vector x can be written as x = κ + x′ with κ in the span of the
known kernel and x′ supported on the kept coordinates, and ⟨x, Yx⟩ = ⟨x′, Yx′⟩. Hence: **the block is ⪰ 0 with kernel
exactly the span of the known vectors iff the gauge-fixed principal submatrix is ≻ 0.**

## 3. The certificate (method (a) of the brief)

`bareiss_cert.py j T order`. Let M_k be the gauge-fixed block in the given order and den the lcm of the denominators
of its entries, den = cst · ∏ f_i^{E_i} with f_i irreducible. Every f_i that occurs has only positive coefficients:
they are |s| (degree 1), the determinants of the compressions of E_c to the relevant column sectors (degrees 2, 3, 4, 6)
and similar factors. Fraction-free (Bareiss) elimination of Nt = den·M_k produces the leading minors D_r(Nt) exactly
(Sylvester's identity; the divisions by the previous pivot are exact). Every intermediate entry is stored as
∏ f_i^{x_i} · core with the factors f_i divided out exactly after each step. Then

    D_r(M_k) = core_r · ∏_i f_i^{x_ir − r E_i} / cst^r.

If cst > 0, every f_i and every core_r have only positive coefficients, then every leading principal minor of M_k is
positive at every point of the open orthant (s₁, s₂, s₃ > 0, and by homogeneity every s ∈ (0, ∞)⁴), and M_k ≻ 0 by
Sylvester's criterion. The script checks every D_r(Nt) against a direct exact determinant at a random integer point,
and `check_minors_points.py` compares the certified formula with leading minors of the block built directly over Q
(independently of the symbolic pipeline) at four rational points per block. No floating point enters the proof.

An independent route, `elim_minors.py` (symmetric elimination over the rational function field, D_r = D_{r−1}·pivot_r
as reduced fractions), gives the same cores; it is about six times slower, and is rerun only for the two smallest blocks.

## 4. Results

All in the natural order of the states (rows grouped by the row letter, then the vacuum states last; the exact order is
printed at the top of every log). In every certified block every leading-minor core has only positive coefficients, so no other
order was needed. (The run for the letter-1 block of (R∞)₁ used the same script with division with remainder in place
of exact division, kept as `bareiss_cert_divrem_version.py.txt`; the arithmetic and the output are the same.)

| block | size | den | cores of D_1 … D_n (degree, terms) | time |
|---|---|---|---|---|
| j = 1, T = {0, 1} | 3 | 16 · \|s\| · [2] · [3]³ | (0, 1), (8, 165), (22, 2 296) | < 1 s |
| j = 1, T = {0} | 10 | 16 · [1]² · [3]³ · [6] (two distinct linear factors) | (2, 3), (3, 2), (3, 1), (12, 165), (19, 680), (26, 1 771), (26, 1 768), (48, 13 234), (55, 19 534), (62, 27 536) | 17 s |
| j = 1, letter-1 | 18 | 8 · \|s\| · [3] · [4] · [6] | (6, 31), (8, 26), (10, 20), (11, 12), (19, 231), (24, 449), (26, 404), (27, 354), (34, 1 100), (41, 2 380), (45, 3 096), (44, 2 452), (63, 12 077), (77, 25 788), (90, 45 335), (91, 45 101), (92, 42 987), (93, 39 341) | 10.4 min |
| j = 2, T = {0, 1, 2} | 4 | 16 · \|s\| · [2]³ · [3]³ | (0, 1), (6, 84), (18, 1 330), (31, 5 980) | < 1 s |
| j = 2, T = {0, 1} | 11 | 16 · (two linear) · [2] · [3]⁴ · [6] | (2, 3), (3, 2), (10, 112), (17, 546), (18, 544), (24, 1 328), (43, 9 137), (50, 14 188), (56, 20 817), (74, 50 100), (79, 59 548) | 2.4 min |
| j = 2, letter-1 | 21 | 8 · \|s\| · [4] · [6] | (9, 81), (12, 69), (14, 48), (16, 34), (18, 20), (20, 10), (35, 791), (44, 1 936), (46, 1 703), (48, 1 633), (50, 1 408), (52, 1 160), (66, 4 788), (80, 12 216), (88, 17 156), (89, 16 121), (90, 15 006), (91, 12 341), (102, 22 099), (112, 34 162), (121, 47 905) | 10.7 min |

([d] denotes an irreducible factor of degree d with only positive coefficients; the exact factors and exponents are
in the logs and in the `cert_*.pkl` files.)

**Theorem A ((R∞)₁ at k = 4).** For k = 4 and every positive weight vector, the condition (R∞)₁ holds: every torus
block is ⪰ 0, the kernel is exactly p ⊗ H₁ (the letter-1 block contributes p ⊗ |1@c⟩, the blocks T = {t} contribute
p ⊗ |2@t⟩, the blocks |T| = 2 and the type II blocks are positive definite).

*Proof.* The type II blocks by Proposition T2′. The other blocks, up to permutations of the copies (Section 2), are
the three representatives j = 1, T = ∅, {0}, {0, 1}; for each, the gauge reduction of Section 2 and the certificate of
Section 3 (rows 1–3 of the table). ∎

**Theorem B (part of (R∞)₂ at k = 4).** For k = 4 and every positive weight vector:
- the blocks of (R∞)₂ with |T| = 3 are positive definite, and those with |T| = 2 are ⪰ 0 with kernel exactly
  p ⊗ |2@t 2@t′⟩;
- the letter-1 block of (R∞)₂ is ⪰ 0 with kernel exactly p ⊗ (letter-1 part of H₂) ⊕ span(vac ⊗ J⁻|1111⟩). This is the goal (Proposition K2⁺) of `relative-anchors/code/middle_k2/NOTES.md` at k = 4, for all positive weights. The certificate is exact on the whole open orthant, as Corollary S requires: the leading minors are positive there and are allowed to tend to 0 at its boundary.

*Proof.* As for Theorem A, with the rows of the table for j = 2. ∎

What is missing for all of (R∞)₂ at k = 4 is the orbit of the four blocks with |T| = 1 (Section 5). If they are
⪰ 0 with kernel p ⊗ |2@t 1@c⟩, then (R∞)₂ holds at k = 4 with kernel exactly p ⊗ H₂ ⊕ span(vac ⊗ J⁻|1111⟩), and with
Proposition M, Theorem B0, Theorem A, the README's results for (R∞)₃ (Theorems T, C, U2, U3₁, Propositions T2, U0,
U1) and Theorem K, all of (R∞) holds at k = 4. By Lemma 13 and the remark that (R2) always holds, together with
Theorem 16 for smaller supports, that would give Conjecture 12 for every p ∈ H₁ supported on at most four copies. We
do not claim this.

## 5. What was not completed, and why

- **The block j = 2, T = {t} (size 21, gauge-fixed 18) is not certified.** It was not run within the time budget;
  nothing is claimed for it beyond the exact point checks (positive definite after the gauge at 33
  rational weight vectors, `num_k4.py`). The obstacle is computing time, not a sign failure that we observed. The
  cost estimate comes from the other blocks: at equal size, the leading-minor cores of the type I blocks are about
  ten times larger than those of the letter-1 blocks (e.g. D₁₀: 27 536 terms for j = 1, T = {0}, against 2 380 for the
  letter-1 block of (R∞)₁), and the letter-1 block of (R∞)₁ (size 18) already took 10.4 minutes, dominated by the last
  five steps. For T = {0} the denominator has five irreducible factors (one linear, two of degree 3, one of degree 4,
  one of degree 6) and the entries have numerator degree up to 8. We expect a run of the order of hours with the present
  code. `run_all.sh` attempts it with `FULL=1`.
- Suggested routes for that block: (i) run `bareiss_cert.py 2 0 nat` without a time limit (memory stayed below 0.6 GB
  in all our runs, so memory is not the issue); (ii) first eliminate the (2@t, ·) states, which play the role of the
  letter-2 states of Theorem C, and certify the Schur complement on the remaining states, as in Theorem U3₁; (iii) the
  Z-matrix route (c) was not tried.
- **Speed.** Elimination over the rational function field (`elim_minors.py`, reduced fractions with gcds at every
  step) is correct but about six times slower than the Bareiss scheme: 112 s against 17 s for j = 1, T = {0}; for
  j = 2, T = {0, 1} it needed 216 s for the first nine minors and was stopped (`logs/abandoned/`), while the Bareiss
  scheme did all eleven in 144 s. Replacing division with remainder by flint's exact division did not change the time.

## 6. Numerical facts (floating point; no claim rests on them)

From `num_k4.py` (33 weight vectors with s₀ = 1 and s₁, s₂, s₃ log-uniform with spreads up to e^{±5}, plus equal
weights and two near-degenerate points): the smallest eigenvalue of the gauge-fixed blocks goes down to 1.7·10⁻¹⁰
(j = 1, letter-1), 9.3·10⁻⁷ (j = 1, T = {0}), 1.1·10⁻³ (j = 1, |T| = 2), 1.5·10⁻¹³ (j = 2, letter-1), 1.1·10⁻⁹
(j = 2, T = {0}), 1.5·10⁻⁶ (j = 2, |T| = 2), 1.9·10⁻³ (j = 2, |T| = 3), always at the near-degenerate points. In the
rescaled coordinates these numbers mix the gap with the scaling of the coordinates; they only say that the blocks
degenerate at the boundary of the orthant, consistent with Corollary S. This is why only certificates that are exact
on the whole open orthant (positive coefficients) are used.

## 7. Other orders (an observation)

The certificates use the natural order, but the positivity does not seem to depend on it. Running `bareiss_cert.py`
on other orders gives the same result:
- block j = 1, T = {0}: every leading-minor core has only positive coefficients in eight orders. These are the
  reversed order, vacuum states first, ordered by column, and five random orders.
- block j = 2, T = {0, 1}: the same holds in two random orders.

Logs are in `logs/other_orders/`; the first line of each log gives the order. So in these blocks probably every
principal minor, not only the leading ones, is a ratio of polynomials with only positive coefficients. We have no
explanation, and nothing above depends on this.

## 8. Files

| File | Content |
|---|---|
| `run_all.sh` | reruns every step, one Python process at a time, PASS/FAIL per step and ALL PASS; `LONG=0` skips the large certificates |
| `ylib.py`, `flintfield.py` | copied from `relative-anchors/code/` (exact builder of the (R∞)_j blocks in the rescaled basis; rational function field) |
| `mid.py` | copied from `relative-anchors/code/middle_k2/` (not used by the certificates) |
| `k4lib.py` | block states, known kernel vectors, gauge, symbolic and point builders |
| `num_k4.py` | torus structure (exact), known kernels and gauge-fixed positivity at 33 rational points (exact), floating-point spectra |
| `build_sym.py` → `blocks_sym.pkl` | exact blocks over Q(s₁, s₂, s₃), kernels annihilated identically, denominators |
| `check_direct.py` | blocks against the dense construction from the definitions (float) and symbolic blocks against the exact point builder |
| `bareiss_cert.py` → `cert_*.pkl` | the certificate of Section 3 |
| `elim_minors.py` | independent elimination over the rational function field (cross-check, small blocks) |
| `check_minors_points.py` | certified minors against direct exact minors at rational points |
| `controls.py` | negative controls with exact witnesses |
| `logs/` | all logs; `logs/abandoned/` keeps the logs of stopped runs; `logs/other_orders/` the runs of Section 7 |
