# Seven copies at two weight values

5 November 2026. Nothing here is refereed.

This folder extends `code/middle_k56/` from k = 5, 6 to k = 7. The setting, the reductions (Lemma R), the level
exchange, the isotypic reduction, the certificate and the checks are those of `code/middle_k56/NOTES.md` (Sections 1–3,
5, 7 and 9), run with `k=7`. The scripts are copies of those in `code/middle_k56/`; the only change is a fix to
`check_builders.py` (Section 4). Notation and names are as there: the family jJmM is the type I block of (R∞)_J with
|T| = M copies carrying the letter 2 once (M = 0: the letter-1 block), and a configuration (a, b) puts squared weight 1
on a copies of T and on b copies outside T, and q elsewhere.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Theorem A₇.** k = 7, weights with at most two distinct values. Each of the following blocks is ⪰ 0 with kernel exactly the span of its known kernel vectors (p ⊗ y for the columns y of its letter-1 rows; for the letter-1 block of (R∞)₅ also vac ⊗ J⁻\|1111111⟩): (a) the letter-1 block and the type I blocks with \|T\| = 1, 2, 3 of (R∞)₃; (b) the letter-1 block and the type I blocks with \|T\| = 1, …, 4 of (R∞)₄; (c) the letter-1 block and the type I blocks with \|T\| = 1, …, 5 of (R∞)₅; (d) the type I blocks of (R∞)₆ with \|T\| = 2, 3, 4 | **proved**, computer-assisted, exact: 148 certificates (18 families, every configuration up to the level exchange), all by the isotypic reduction of `code/middle_k56/` Section 7 (`cert_sym.py`); every leading minor of every part (4 955 minors in 886 parts) has numerator and denominator with only nonnegative coefficients |
| 2 | **Theorem 7₂.** For k = 7 and every weight vector with at most two distinct weights, all of (R∞) holds: (R∞)_j ⪰ 0 for j = −1, …, 7; the kernels of (R∞)₃ and (R∞)₄ are exactly p ⊗ H₃ and p ⊗ H₄, and that of (R∞)₅ is exactly p ⊗ H₅ ⊕ span{vac ⊗ J⁻\|1111111⟩} | **proved**, given row 1 and the cited results (Section 2) |
| 3 | **Corollary 7₂.** Conjecture 12 of `single-excitation-range/` holds for every p ∈ H₁ supported on at most seven copies whose nonzero weights take at most two distinct values | **proved**, given row 2, Corollary 6₂ of `code/middle_k56/` and Lemma 13 with the remark that (R2) always holds (`single-excitation-range/`) |
| 4 | The checks of `code/middle_k56/` at k = 7: coverage of every (block, level pattern) pair; every certified part and minor against the Lemma B block; the saved blocks against both builders; homogeneity; negative controls; end to end at points (with the faster exact routines of `fast/`, Section 3) | **exact checks** (Section 3); not proofs |
| 5 | A second check written separately (`crosscheck.py`, as in `code/middle_k56/` Section 9): at 24 random two-value points at k = 7, every torus block of every (R∞)_j, j = −1, …, 7, is ⪰ 0 with exactly its known kernel; the 148 saved blocks at 31 values of q from 10⁻⁸ to 10⁸ | **exact check** at points (Section 5); not a proof |
| 6 | k ≥ 8 at two weight values; three or more weight values; general weights | **open** |

## 1. The families and the certificates

At k = 7 and two weight values the blocks not covered by earlier results are those of (R∞)₃, (R∞)₄ and (R∞)₅ that
contain a letter 1, and the type I blocks of (R∞)₆ = (R∞)_{k−1} with |T| = 2, 3, 4 (|U| = 5, 4, 3); (R∞)₂ is covered
for k ≤ 10 by `code/middle_j2_l1/`. These are the eighteen families below. Every configuration is certified up to the
level exchange (`k5lib.rep_configs`), all by `cert_sym.py`: the block is built over Q(q) with `ylib.py`, the known
kernel vectors are checked identically, and the form on each e_λW (Young symmetrizers of the level-preserving
permutations on the invariant complement W of the kernel) is certified by its leading principal minors.

Sizes (`logs/run_cert_k7_*.log`; parts and minors summed over the configurations; degree: the largest numerator
degree of a minor; time summed over the configurations):

| family | block | states | known kernel | dim W | configurations | parts (largest) | minors | degree | time |
|---|---|---|---|---|---|---|---|---|---|
| j3m0 | (R∞)₃, letter-1 | 280 | 35 | 245 | 4 | 34 (16) | 198 | 564 | 58 s |
| j3m1 | (R∞)₃, \|T\| = 1 | 160 | 15 | 145 | 7 | 46 (21) | 342 | 1 041 | 75 s |
| j3m2 | (R∞)₃, \|T\| = 2 | 75 | 5 | 70 | 9 | 59 (21) | 327 | 1 196 | 73 s |
| j3m3 | (R∞)₃, \|T\| = 3 | 26 | 1 | 25 | 10 | 44 (11) | 166 | 390 | 59 s |
| j4m0 | (R∞)₄, letter-1 | 266 | 35 | 231 | 4 | 34 (15) | 188 | 431 | 51 s |
| j4m1 | (R∞)₄, \|T\| = 1 | 190 | 20 | 170 | 7 | 46 (22) | 386 | 970 | 83 s |
| j4m2 | (R∞)₄, \|T\| = 2 | 120 | 10 | 110 | 9 | 65 (27) | 474 | 1 584 | 95 s |
| j4m3 | (R∞)₄, \|T\| = 3 | 64 | 4 | 60 | 10 | 64 (21) | 356 | 1 007 | 81 s |
| j4m4 | (R∞)₄, \|T\| = 4 | 26 | 1 | 25 | 10 | 44 (11) | 166 | 335 | 61 s |
| j5m0 | (R∞)₅, letter-1 | 154 | 21 + 1 | 132 | 4 | 31 (10) | 118 | 172 | 27 s |
| j5m1 | (R∞)₅, \|T\| = 1 | 132 | 15 | 117 | 7 | 46 (17) | 286 | 486 | 53 s |
| j5m2 | (R∞)₅, \|T\| = 2 | 105 | 10 | 95 | 9 | 65 (24) | 420 | 1 066 | 77 s |
| j5m3 | (R∞)₅, \|T\| = 3 | 76 | 6 | 70 | 10 | 64 (23) | 410 | 953 | 83 s |
| j5m4 | (R∞)₅, \|T\| = 4 | 48 | 3 | 45 | 10 | 49 (17) | 285 | 564 | 72 s |
| j5m5 | (R∞)₅, \|T\| = 5 | 24 | 1 | 23 | 9 | 35 (11) | 135 | 269 | 53 s |
| j6m2 | (R∞)₆, \|T\| = 2 | 48 | 5 | 43 | 9 | 53 (15) | 219 | 314 | 49 s |
| j6m3 | (R∞)₆, \|T\| = 3 | 44 | 4 | 40 | 10 | 58 (15) | 251 | 314 | 60 s |
| j6m4 | (R∞)₆, \|T\| = 4 | 38 | 3 | 35 | 10 | 49 (14) | 228 | 307 | 60 s |

In all 148 certificates every leading minor has numerator and denominator with only nonnegative coefficients; no
witness search and no Sturm count was needed (as at k = 5 and 6; compare Conjecture SF of Section 3g of the README).

## 2. Proofs

**Theorem A₇.** *Proof.* As for Theorem B of `code/middle_k56/` (Section 7 there): by Lemma R and the level exchange
every block of the eighteen families at every two-value weight vector reduces to a certified configuration
(`coverage.py` checks this for every pair (block, level pattern)), and each certificate shows that Ỹ restricted to W is
positive definite for every q > 0. ∎

**Theorem 7₂.** *Proof.* Every torus block is covered:
- (R∞)_{−1}, (R∞)_0, (R∞)_1, (R∞)_2 and (R∞)_7 as in the proof of Theorem 5₂ of `code/middle_k56/` (Proposition M,
  Theorem B0, the Corollary of Section 3h of the README for (R∞)_1 at two weight values, Corollary L for 5 ≤ k ≤ 10,
  Theorem K);
- (R∞)₃: Theorem A₇ (a); |T| = 4 by Theorem DJ (no letter 1); type II blocks by Proposition T2′;
- (R∞)₄: Theorem A₇ (b); |T| = 5 by Theorem DJ; type II by Proposition T2′;
- (R∞)₅ = (R∞)_{k−2}: Theorem A₇ (c); |T| = 6 by Theorem DJ; type II by Proposition T2′;
- (R∞)₆ = (R∞)_{k−1}: Theorem T (letter-1 block); |T| = 1 by Theorems EU and TL; |T| = 2, 3, 4 by Theorem A₇ (d);
  |T| = 5, 6, 7 by Theorem U2, Proposition U1 and Proposition U0; type II by Proposition T2.
The kernels follow as in the proof of Theorem 5₂: each block has exactly the known kernel, and these add up to
p ⊗ H_j, plus the vector of Proposition K2 in (R∞)₅. ∎

**Corollary 7₂.** *Proof.* With exactly seven nonzero weights, Theorem 7₂, the remark that (R2) always holds and
Lemma 13 of `single-excitation-range/`; with at most six, Corollary 6₂. ∎

## 3. Checks at k = 7

All exact. The scripts are those of `code/middle_k56/`, run with `k=7`, one Python process at a time under `guard.sh`.
- `coverage.py` (`logs/run_coverage_k7.log`): every one of the 47 872 pairs (block, level pattern) of the eighteen
  families (23 936 of them through the level exchange) reduces to a certified configuration; the exact identity of
  Section 2 of `code/middle_k56/` is checked for every 800th pair (71 pairs).
- `check_sym.py` (`logs/run_sym_k7.log`): at two values of q, every certified part F_λ (886 parts, 1 772 (part, q)
  pairs) equals den · B_λỸB_λᵀ with Ỹ the Lemma B block, every certified minor equals the direct minor of that matrix
  (9 910 comparisons), the basis vectors lie in W, and every part is positive definite.
- `check_builders.py` (`logs/run_k7_builders_builders_k7.log`, the rerun of Section 4): every certified block equals the
  Lemma B block at three rational squares q per configuration (444 comparisons) and the `ylib.py` point block at
  q = 9/16 (148); the two builders agree at two general weight vectors for every family, with B scaled by 1 and by
  1001/1000 (72 comparisons); homogeneity of degree j + 1 − |T| (36 checks).
- Negative controls: with B → (1 + 1/1000)B every certificate of a two-level representative configuration fails
  (130 of 130, `logs/run_neg_cert_k7_*.log`), and `neg_points.py` gives an explicit rational point for each of them
  (q = 10⁻⁶ for 114, q = 10⁴ for 16) where the scaled gauge-fixed block built by `j1lib.py` is not positive definite
  while the true one is.
- End to end (`fast/`, `logs/run_end2end_k7_fast.log`): at five values of q per representative configuration
  (1/4, 9/4, 10⁻⁶, 10⁶, 121/100; 740 (configuration, q) points), the known vectors are annihilated, the exact nullity
  equals their number, and the gauge-fixed block is positive definite in two gauges. The step `end2end_k7` of the first
  run did the same with the pure-Python elimination of `k5lib.py`; it was too slow at the extreme ratios (about 90 s
  per elimination of a 245-dimensional block at q = 10⁶ even with the fast routines) and was stopped after 5 546 s,
  before finishing its first family, so `logs/run_k7.out` records it as FAIL. The check was run instead with the
  same `check_end2end.py` on top of two faster exact routines, kept in `fast/`: `fastj1.py`, a Lemma B builder that
  groups the E⁻¹ terms (the same entries as `j1lib.J1.matrix`; checked entry by entry against it on every torus group
  of every (R∞)_j for k = 3, …, 6, with B scaled by 1 and by 1001/1000, by `validate.py`), and `fastpd.py`, the exact
  symmetric elimination by recursive Schur complements in python-flint (the same result, including the index of the
  first nonpositive pivot, as the sequential elimination; checked on random rational matrices by `validate_pd.py`).
  `k5lib.py` in `fast/` uses them for its Lemma B blocks, eliminations and products; nothing else differs.

## 4. The change to check_builders.py

In the first run (`logs/run_k7.out`) the step `builders_k7` stopped with an error after its first family: part (c)
of `check_builders.py` compares the two builders at two general weight vectors, which were written as lists of six
weights cut to length k, so at k = 7 the Lemma B builder was called with six weights for seven copies (a KeyError in
`j1lib.py`). Parts (a) and (b) had passed for that family. A seventh weight was added to both lists (11/5 and 6/5),
with an assertion on the length; nothing else changed. The step was then rerun on its own with the corrected script
(`run_k7_builders.sh`, `logs/run_k7_builders.out`). The other steps of the first run did not use this script.

## 5. A second check

`crosscheck.py` (the second check of `code/middle_k56/`, Section 9 there), run with `run_crosscheck_k7.sh`:
- negative controls of the exact positive-definiteness test;
- every torus block of every (R∞)_j, j = −1, …, 7, at 24 random two-value points (levels at random positions, weights
  from a list of rationals): the states of G1_j ⊕ G0_{j+1} grouped by the number of letters 2 at each copy (4 480
  states in 253 groups for (R∞)₃, 8 512 in 414 for (R∞)₄, 9 856 in 526 for (R∞)₅), entries between groups vanishing on a
  random sample, every group ⪰ 0 with exactly its known kernel (annihilated exactly; positive definite without a set
  of pivot coordinates of them, by the exact test), and the kernel dimensions adding up to |H_j| (+ 1 at j = 5,
  + 6 at j = 6): 216 conditions, all pass (`logs/run_crosscheck_all_blocks_k7.log`, 3 651 s);
- the 148 saved blocks (`logs/cert_k7_*.pkl`, the rescaled blocks over Q(q)) at all 31 values of q from 10⁻⁸ to
  10⁸: the column-group indicators (and, for j5m0, the vector (s_h) on the states (vac, 1…1 − e_h)) are annihilated,
  and the block without a set of pivot coordinates of them is positive definite by an exact test: the integer
  congruence certificate pd_cert (4 583 of the 4 588 cases), and the exact elimination of `fast/fastpd.py` in the 5
  cases where pd_cert gave no certificate (`logs/run_crosscheck_saved_blocks_k7_fast.log`, 673 s; the entries are
  evaluated with flint polynomials, and at one value of q per block 400 entries are compared with the evaluation in
  Fraction arithmetic). In `run_crosscheck_k7.sh` this step used the pure-Python elimination; it took about 23
  minutes per block of 280 states and was stopped after 2 768 s, so `logs/run_crosscheck_k7.out` records it as
  FAIL. It was rerun with `fast/run_saved_fast.sh`.

## Files

| File | Content |
|---|---|
| `k5lib.py`, `ylib.py`, `j1lib.py`, `cert5.py`, `cert_sym.py`, `symred.py`, `signtest.py`, `saved.py`, `flintfield.py`, `ufield.py`, `guard.sh` | copied unchanged from `code/middle_k56/` |
| `coverage.py`, `check_sym.py`, `check_minors.py`, `check_end2end.py`, `neg_points.py`, `explore_kernel.py` | copied unchanged from `code/middle_k56/` |
| `check_builders.py` | from `code/middle_k56/`, with the fix of Section 4 |
| `crosscheck.py` | copied unchanged from `code/middle_k56/` (the second check) |
| `run_k7.sh`, `run_k7_builders.sh`, `run_crosscheck_k7.sh` | the run at k = 7, the rerun of the builder step, the second check |
| `fast/` | `k5lib.py` using `fastj1.py` (Lemma B builder) and `fastpd.py` (exact elimination), their checks `validate.py` and `validate_pd.py`, `crosscheck.py` using `fastpd.py` for blocks of size above 24 (and, for the saved blocks, flint evaluation and pd_cert first), copies of the other scripts that `check_end2end.py` imports, and the runners `run_end2end_fast.sh` and `run_saved_fast.sh` |
| `logs/run_k7.out`, `logs/run_*.log` | the first run (148 certificates and the checks; the step `builders_k7` failed as described in Section 4, and the step `end2end_k7` was stopped as described in Section 3) |
| `logs/run_end2end_k7_fast.log`, `logs/validate_*.log` | the end-to-end check with the routines of `fast/`, and the checks of those routines |
| `logs/run_k7_builders.out`, `logs/run_k7_builders_*.log` | the rerun of the builder step |
| `logs/run_crosscheck_k7.out`, `logs/run_crosscheck_*.log` | the second check (its saved-block step stopped and rerun with `fast/run_saved_fast.sh`: `logs/run_crosscheck_saved_blocks_k7_fast.log`) |
| `logs/cert_k7_*.pkl` | the saved blocks, minors and parts (regenerated by `run_k7.sh`; not in the repository) |
