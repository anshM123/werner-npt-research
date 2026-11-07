# Eight copies at two weight values

6 November 2026. Nothing here is refereed.

This folder extends `code/middle_k7/` from k = 7 to k = 8. The setting, the reductions (Lemma R), the level exchange,
the isotypic reduction, the certificate and the checks are those of `code/middle_k56/NOTES.md` (Sections 1–3, 5, 7
and 9), run with `k=8`. The scripts are those of `code/middle_k7/` and its folder `fast/`: the Lemma B builder of
`fast/fastj1.py` and the exact elimination of `fast/fastpd.py` are used throughout (`k5lib.py` here is the one of
`fast/`), and the second check uses the flint evaluation and the integer congruence certificate of Section 5 of
`code/middle_k7/NOTES.md`. Notation as there: the family jJmM is the type I block of (R∞)_J with |T| = M copies
carrying the letter 2 once (M = 0: the letter-1 block), and a configuration (a, b) puts squared weight 1 on a copies of
T and on b copies outside T, and q elsewhere.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Theorem A₈.** k = 8, weights with at most two distinct values. Each of the following blocks is ⪰ 0 with kernel exactly the span of its known kernel vectors (p ⊗ y for the columns y of its letter-1 rows; for the letter-1 block of (R∞)₆ also vac ⊗ J⁻\|11111111⟩): (a) the letter-1 block and the type I blocks with \|T\| = 1, 2, 3 of (R∞)₃; (b) the same with \|T\| = 1, …, 4 of (R∞)₄; (c) with \|T\| = 1, …, 5 of (R∞)₅; (d) with \|T\| = 1, …, 6 of (R∞)₆; (e) the type I blocks of (R∞)₇ with \|T\| = 2, 3, 4, 5 | **proved**, computer-assisted, exact: 266 certificates (26 families, every configuration up to the level exchange), all by the isotypic reduction of `code/middle_k56/` Section 7 (`cert_sym.py`); every leading minor of every part (11 307 minors in 1 912 parts) has numerator and denominator with only nonnegative coefficients |
| 2 | **Theorem 8₂.** For k = 8 and every weight vector with at most two distinct weights, all of (R∞) holds: (R∞)_j ⪰ 0 for j = −1, …, 8; the kernels of (R∞)₃, (R∞)₄, (R∞)₅ are exactly p ⊗ H_j, and that of (R∞)₆ is exactly p ⊗ H₆ ⊕ span{vac ⊗ J⁻\|11111111⟩} | **proved**, given row 1 and the cited results (Section 2) |
| 3 | **Corollary 8₂.** Conjecture 12 of `single-excitation-range/` holds for every p ∈ H₁ supported on at most eight copies whose nonzero weights take at most two distinct values | **proved**, given row 2, Corollary 7₂ of `code/middle_k7/` and Lemma 13 with the remark that (R2) always holds (`single-excitation-range/`) |
| 4 | The checks of `code/middle_k56/` at k = 8 (Section 3) | **exact checks**, not proofs |
| 5 | A second check written separately (`crosscheck.py`, Section 5) | **exact check** at points, not a proof |
| 6 | k ≥ 9 at two weight values for 3 ≤ j ≤ k − 2 (the other conditions hold at two weight values for every k: README, Sections 3 and 3h); three or more weight values; general weights | **open** |

## 1. The families and the certificates

At k = 8 and two weight values the blocks not covered by results for every k are those of (R∞)₃, (R∞)₄, (R∞)₅ and
(R∞)₆ that contain a letter 1, and the type I blocks of (R∞)₇ = (R∞)_{k−1} with |T| = 2, 3, 4, 5 (|U| = 6, 5, 4, 3).
(R∞)₂ holds at two weight values for every k (Corollary L∞), and (R∞)₇ also by Corollary TV, which was proved after
this run; Theorem A₈ (e) is an independent certificate of those blocks at k = 8. These are the twenty-six families
below. Every configuration is certified up to the level exchange (`k5lib.rep_configs`), all by `cert_sym.py`: the block
is built over Q(q), the known kernel vectors are checked identically, and the form on each e_λW (Young symmetrizers
of the level-preserving permutations on the invariant complement W of the kernel) is certified by its leading
principal minors.

Sizes (`logs/run_cert_k8_*.log`; parts and minors summed over the configurations; degree: the largest numerator
degree of a minor; time summed over the configurations):

| family | block | states | known kernel | dim W | configurations | parts (largest) | minors | degree | time |
|---|---|---|---|---|---|---|---|---|---|
| j3m0 | (R∞)₃, letter-1 | 518 | 56 | 462 | 5 | 55 (18) | 333 | 783 | 806 s |
| j3m1 | (R∞)₃, \|T\| = 1 | 259 | 21 | 238 | 8 | 62 (21) | 458 | 1 293 | 870 s |
| j3m2 | (R∞)₃, \|T\| = 2 | 105 | 6 | 99 | 11 | 80 (21) | 451 | 1 364 | 959 s |
| j3m3 | (R∞)₃, \|T\| = 3 | 31 | 1 | 30 | 12 | 56 (11) | 210 | 379 | 1 162 s |
| j4m0 | (R∞)₄, letter-1 | 616 | 70 | 546 | 5 | 55 (19) | 373 | 780 | 901 s |
| j4m1 | (R∞)₄, \|T\| = 1 | 385 | 35 | 350 | 8 | 68 (25) | 606 | 1 636 | 1 126 s |
| j4m2 | (R∞)₄, \|T\| = 2 | 210 | 15 | 195 | 11 | 103 (31) | 757 | 2 587 | 1 444 s |
| j4m3 | (R∞)₄, \|T\| = 3 | 95 | 5 | 90 | 12 | 92 (21) | 502 | 1 280 | 1 460 s |
| j4m4 | (R∞)₄, \|T\| = 4 | 32 | 1 | 31 | 13 | 65 (11) | 241 | 445 | 1 456 s |
| j5m0 | (R∞)₅, letter-1 | 476 | 56 | 420 | 5 | 55 (16) | 306 | 459 | 749 s |
| j5m1 | (R∞)₅, \|T\| = 1 | 357 | 35 | 322 | 8 | 68 (23) | 566 | 1 174 | 1 142 s |
| j5m2 | (R∞)₅, \|T\| = 2 | 245 | 20 | 225 | 11 | 103 (33) | 844 | 2 440 | 1 508 s |
| j5m3 | (R∞)₅, \|T\| = 3 | 150 | 10 | 140 | 12 | 100 (27) | 722 | 1 692 | 1 532 s |
| j5m4 | (R∞)₅, \|T\| = 4 | 78 | 4 | 74 | 13 | 93 (21) | 518 | 1 070 | 1 184 s |
| j5m5 | (R∞)₅, \|T\| = 5 | 31 | 1 | 30 | 12 | 56 (11) | 210 | 335 | 955 s |
| j6m0 | (R∞)₆, letter-1 | 232 | 29 | 203 | 5 | 45 (10) | 174 | 182 | 322 s |
| j6m1 | (R∞)₆, \|T\| = 1 | 203 | 21 | 182 | 8 | 62 (17) | 370 | 503 | 620 s |
| j6m2 | (R∞)₆, \|T\| = 2 | 168 | 15 | 153 | 11 | 99 (25) | 622 | 1 186 | 936 s |
| j6m3 | (R∞)₆, \|T\| = 3 | 130 | 10 | 120 | 12 | 100 (24) | 638 | 1 163 | 1 086 s |
| j6m4 | (R∞)₆, \|T\| = 4 | 92 | 6 | 86 | 13 | 93 (23) | 596 | 976 | 1 156 s |
| j6m5 | (R∞)₆, \|T\| = 5 | 57 | 3 | 54 | 12 | 62 (17) | 360 | 581 | 1 012 s |
| j6m6 | (R∞)₆, \|T\| = 6 | 28 | 1 | 27 | 11 | 44 (11) | 173 | 247 | 852 s |
| j7m2 | (R∞)₇, \|T\| = 2 | 63 | 6 | 57 | 11 | 71 (15) | 292 | 329 | 817 s |
| j7m3 | (R∞)₇, \|T\| = 3 | 59 | 5 | 54 | 12 | 80 (15) | 334 | 359 | 1 007 s |
| j7m4 | (R∞)₇, \|T\| = 4 | 53 | 4 | 49 | 13 | 83 (15) | 363 | 404 | 1 178 s |
| j7m5 | (R∞)₇, \|T\| = 5 | 45 | 3 | 42 | 12 | 62 (14) | 288 | 307 | 1 025 s |

In all 266 certificates every leading minor has numerator and denominator with only nonnegative coefficients (every
denominator of a minor is a constant); no witness search and no Sturm count was needed (as at k = 5, 6 and 7;
compare Conjecture SF of Section 3g of the README).

One certificate, `cert_k8_j4m0_0_0` (the letter-1 block of (R∞)₄ at equal weights, configuration (0, 0)), was
stopped by the 2.5 GB limit of `guard.sh` in the first run (`logs/first_attempt_run_cert_k8_j4m0_0_0.log`) and rerun
alone with a 5.2 GB limit (`fix_big.sh`, `logs/fix_big.out`), where it passed in 214 s; the log
`logs/run_cert_k8_j4m0_0_0.log` is that of the rerun. (At equal weights this block is also Theorem E1.)

## 2. Proofs

**Theorem A₈.** *Proof.* As for Theorem B of `code/middle_k56/` (Section 7 there): by Lemma R and the level exchange
every block of the twenty-six families at every two-value weight vector reduces to a certified configuration
(`coverage.py` checks this for every pair (block, level pattern)), and each certificate shows that Ỹ restricted to W
is positive definite for every q > 0. ∎

**Theorem 8₂.** *Proof.* Every torus block is covered:
- (R∞)_{−1}, (R∞)_0, (R∞)_1, (R∞)_2 and (R∞)_8 as in the proof of Theorem 5₂ of `code/middle_k56/` (Proposition M,
  Theorem B0, the Corollary of Section 3h of the README for (R∞)_1 at two weight values, Corollary L for 5 ≤ k ≤ 10,
  Theorem K);
- (R∞)₃: Theorem A₈ (a); |T| = 4 by Theorem DJ (no letter 1); type II blocks by Proposition T2′;
- (R∞)₄: Theorem A₈ (b); |T| = 5 by Theorem DJ; type II by Proposition T2′;
- (R∞)₅: Theorem A₈ (c); |T| = 6 by Theorem DJ; type II by Proposition T2′;
- (R∞)₆ = (R∞)_{k−2}: Theorem A₈ (d); |T| = 7 by Theorem DJ; type II by Proposition T2′;
- (R∞)₇ = (R∞)_{k−1}: Theorem T (letter-1 block); |T| = 1 by Theorems EU and TL; |T| = 2, …, 5 by Theorem A₈ (e)
  (and by Corollary TV); |T| = 6, 7, 8 by Theorem U2, Proposition U1 and Proposition U0; type II by Proposition T2.
The kernels follow as in the proof of Theorem 5₂: each block has exactly the known kernel, and these add up to
p ⊗ H_j, plus the vector of Proposition K2 in (R∞)₆. ∎

**Corollary 8₂.** *Proof.* With exactly eight nonzero weights, Theorem 8₂, the remark that (R2) always holds and
Lemma 13 of `single-excitation-range/`; with at most seven, Corollary 7₂. ∎

## 3. Checks at k = 8

All exact, one Python process at a time under `guard.sh`.
- The fast routines (`logs/run_fastpd.log`, `logs/run_fastj1_k*.log`): `fastpd.py` gives the same result and the same
  first nonpositive pivot as the sequential elimination on 30 random rational matrices (20 positive definite);
  `fastj1.py` equals `j1lib.J1` entry by entry on every torus group of every (R∞)_j for k = 4, 5, 6 (576, 1 536 and
  1 984 groups; 34 620, 229 012 and 723 670 entries), with B scaled by 1 and by 1001/1000.
- `coverage.py` (`logs/run_coverage_k8.log`): every one of the 238 592 pairs (block, level pattern) of the twenty-six families (110 756 of
  them through the level exchange) reduces to a certified configuration; the exact identity of Section 2 of
  `code/middle_k56/` is checked for 1 204 of the pairs. (The first run of this step, `logs/first_run_coverage_k8.log`,
  read the log of `cert_k8_j4m0_0_0` while that certificate was being rerun, so it saw no PASS line for that
  configuration; the rerun by `fix_k8.sh` is the log above.)
- `check_sym.py` (`logs/run_sym_k8.log`): at two values of q, every certified part F_λ (1 912 parts, 3 824 (part, q)
  pairs) equals den · B_λỸB_λᵀ with Ỹ the Lemma B block, every certified minor equals the direct minor of that matrix
  (22 614 comparisons), the basis vectors lie in W, and every part is positive definite.
- `check_builders.py` (`logs/run_builders_k8.log`): every certified block equals the Lemma B block at three rational squares q per
  configuration (798 comparisons) and the `ylib.py` point block at q = 9/16 (266); the two builders agree at two
  general weight vectors for every family, with B scaled by 1 and by 1001/1000 (104 comparisons); homogeneity of
  degree j + 1 − |T| (52 checks). (The first run of this step, `logs/first_run_builders_k8.log`, was stopped two
  seconds after its start, before writing anything; the log above is that of the rerun by `fix_k8.sh`.)
- Negative controls (`neg_points.py`, `logs/run_neg_points_first_k8.log`, 4 760 s): on the first two-level
  representative configuration of each of the twenty-six families (the run over all 240 configurations would take
  about 18 hours at k = 8), the gauge-fixed block with the B-part multiplied by 1001/1000 is not positive definite at
  q = 10⁻⁶ while the true block is positive definite there (exact symmetric elimination): 26 of 26 caught.
- End to end (`logs/run_end2end_first_k8_h1.log`, `logs/run_end2end_first_k8_h2.log`, `logs/run_end2end_first_k8_h1rest.log`): on the first two-level representative configuration of a family, at five (configuration, q) points, the known kernel
vectors are annihilated exactly and the gauge-fixed block is positive definite in two gauges. All twenty-six families
pass. The second half (thirteen families) is in `logs/run_end2end_first_k8_h2.log`. The first four families of the
first half (j3m0–j3m3) are in `logs/run_end2end_first_k8_h1.log`; that run was then stopped by the 2.5 GB limit of
`guard.sh` at j4m0. The other nine families (j4m0–j5m3, the largest blocks) were rerun under `guard_big.sh` by
`finish_k8_v2.sh` (`logs/run_end2end_first_k8_h1rest.log`, PASS, 19 061 s in all; the summary of that script is
`logs/run_k8_v2.out`). This check is a second look at points; the certificates themselves are checked by the steps
above. The run over all 240 configurations would take about thirty hours at k = 8.

## 4. A second check

`crosscheck.py` (the second check of `code/middle_k56/`, Section 9 there, in the form of `code/middle_k7/fast/`):
- negative controls (`logs/run_crosscheck_neg.log`): a singular block, a perturbed gauge-fixed block and a gauge
  missing one kernel vector are refused by pd_cert and by pd_cert2, and blocks with B scaled by 1001/1000 are refused
  at two of four weight ratios by each;
- every torus block of every (R∞)_j, j = −1, …, 8, at twelve random two-value points in three runs of four (seeds
  8208, 8209, 8210; `logs/run_crosscheck_blocks_{a,b,c}.log`): every group is positive semidefinite with kernel exactly
  the span of its known vectors;
- the blocks saved by the certificates, at 31 values of q in [10⁻⁸, 10⁸] (`logs/run_crosscheck_saved.log`): the known
  kernel vectors are annihilated and the rest is positive definite, for every saved block. The exact test is pd_cert.
  Where its floating step fails, pd_cert2 is used: the same exact test, with the rows and columns first scaled by
  powers of 2 and a shifted floating Cholesky step. This happens at q = 10⁸ on the largest blocks, which are positive
  definite there but too ill-conditioned for one double-precision Cholesky step. The exact LDL of `fastpd.py` is the
  last resort. The log gives the three counts for every block: all 266 blocks pass, pd_cert2 was needed at one
  value of q in seven blocks, and the exact LDL was never needed (1.5 hours in all).
The last step took three attempts. The first was stopped by the 2.5 GB limit of `guard.sh` at the 616-state blocks of
j4m0, after every earlier block had passed (`logs/run_k8_finish.out`). The second, under `guard_big.sh`
(`finish_k8_rest.sh`), was stopped by hand after 2.5 hours in the exact LDL of cert_k8_j4m0_0_4 at q = 10⁸, after its
first 40 blocks had passed (`logs/stopped_run_crosscheck_saved.log`). That is why pd_cert2 was added. The third run,
`finish_k8_v2.sh`, repeated the negative controls and every block. (The log of the first attempt was overwritten when
the second was started twice by mistake; both copies were stopped and it was started again.)

## Files

| File | Content |
|---|---|
| `k5lib.py` (using `fastj1.py` and `fastpd.py`), `ylib.py`, `j1lib.py`, `cert5.py`, `cert_sym.py`, `symred.py`, `signtest.py`, `saved.py`, `flintfield.py`, `ufield.py`, `guard.sh` | as in `code/middle_k7/` and its `fast/` |
| `coverage.py`, `check_sym.py`, `check_minors.py`, `check_end2end.py`, `check_builders.py`, `neg_points.py`, `explore_kernel.py`, `crosscheck.py`, `validate.py`, `validate_pd.py` | as in `code/middle_k7/` and its `fast/` |
| `run_k8.sh` | the run at k = 8 (`logs/run_k8.out`) |
| `fix_big.sh`, `guard_big.sh` | the rerun of the certificate stopped by the memory limit (Section 1; `logs/fix_big.out`) |
| `fix_k8.sh` | the rerun of the two steps of `run_k8.sh` that did not count (Section 3; `logs/run_k8_fix.out`) |
| `finish_k8.sh` | what was left after `run_k8.sh` was stopped in its step `neg_points_k8`: the second check in full, and `neg_points.py` and `check_end2end.py` on the first representative configuration of every family (the end-to-end check in two halves) (`logs/run_k8_finish.out`) |
| `finish_k8_rest.sh` | second attempt at the saved-blocks step and the remaining end-to-end families (stopped by hand, see above) |
| `finish_k8_v2.sh` | negative controls, saved-blocks step and remaining end-to-end families, with pd_cert2 (`logs/run_k8_v2.out`) |
| `k8table.py` | the table of Section 1 from the certificate logs |
| `logs/` | all logs of these runs |
| `logs/cert_k8_*.pkl` | the saved blocks, minors and parts (regenerated by `run_k8.sh`; not in the repository) |
