# The two-row part of the generic family: the determinant route

6 November 2026. Nothing here is refereed.

## Status

Notation as in `NOTES_MI.md`: squared weights s_0 = σ > 0 (copy 0 marked), s_1 = … = s_{k−1} = 1, u = k − 1, the
letter-1 block of (R∞)_j with known kernel K_j, irrep index i, q0 = j − i, t = k − 1 − j − i. The *two-row part* of
Q_i is the form Q₂ on the twelve vectors O_{y0}, G_{z0}, U_{x,y0}, W_{x,y0} of Lemma S (`NOTES.md`, Section 2); its
known kernel is span(κ_0, κ_1); Q denotes the gauge-fixed form, Q₂ without U_{1,0} and U_{1,1} (10 × 10).

| # | Statement | Status |
|---|---|---|
| 1 | **Theorem G2.** For every member of the generic family [q0 ≥ 4, t ≥ 4, i ≥ 7] and every σ > 0, Q₂ ⪰ 0 with kernel exactly span(κ_0, κ_1) | **proved**, computer-assisted, exact: det Q = Ñ/D_c with Ñ a polynomial in (x, y, z, σ) with only nonnegative coefficients, determined by exact interpolation from 10 197 exact slices, plus Q ≻ 0 at one point (Sections 2–4) |
| 2 | **Theorem MF.** For every k ≥ 2, every j and every σ > 0, the letter-1 block of (R∞)_j at squared weights (σ, 1, …, 1) is ⪰ 0 with kernel exactly K_j; hence the same for every weight vector whose squared weights are all equal but one | **proved**, given row 1, Theorem G3 and Theorem MF′ of `NOTES_MI.md` (Section 5) |
| 3 | Checks: Ñ against the determinant computed directly at four points far outside the grid; the reduced forms of the family at 240 members with 56 ≤ k ≤ 100 (`mfcheck.py`, region g7, run before this proof) | **exact checks** (Section 6) |
| 4 | Negative control, B → (1 + 1/1000)B: the kernel vectors are no longer annihilated, and the determinant on the slice (t, i) = (4, 7) has coefficients of both signs and is negative at (q0, σ) = (4, 1000) | **exact** (Section 6) |
| 5 | The type I blocks with a letter 1 for k ≥ 30 at these weights (stage 2 of `../typeI/`), two or more marked copies, general weights | **open** |

Nothing here says that Conjecture 12 of `single-excitation-range/` holds in general, or anything about the NPT (Werner
distillability) problem. Computer-assisted means exact rational arithmetic (python-flint `fmpq`, `fmpq_mpoly`,
`fmpq_mat`) and sign checks of coefficients. No floating point enters a proof.

## 1. Why another route

Section 9 of `NOTES_MI.md` explains where the symmetric elimination of the generic two-row part stopped: in the four
variables (x, y, z, σ) (q0 = 4 + x, t = 4 + y, i = 7 + z) the remaining entries reach 10⁵ to 4·10⁵ terms after four
pivots, and each later step takes minutes to hours. The determinant of Q is a single rational function, and it is
cheap to compute exactly once two of the four variables are fixed: at fixed (t, i) the 10 × 10 determinant over
Q(q0, σ) takes 6 to 24 seconds, its numerator has degree 52 in q0 and in σ, and after q0 = 4 + x every one of its
1 492 coefficients is positive. The proof below computes det Q on enough slices (t, i) to determine it, and uses
determinant and inertia (as in Proposition C of the README) instead of the pivots.

## 2. The argument

Let 𝓡 = {(x, y, z, σ) : x, y, z ≥ 0 real, σ > 0}. The family form is the formal continuation of `famlib.py` (as in
`NOTES_MI.md`, Section 2), so at every member and every σ > 0 it is the true reduced form.

1. *Preliminaries* (as in `prove_gen.py`): N0 times the form equals the exact numeric form of `mlib.py` at the 27
   sanity points and two values of σ; every entry denominator is one-signed after the shifts, so Q is continuous on
   𝓡; the split Tᵀ Q_i T = Q₂ ⊕ Q₃ holds identically; κ_0 and κ_1 are annihilated identically and the gauge
   determinant is 1. So Q₂ ⪰ 0 with kernel exactly span(κ_0, κ_1) wherever Q ≻ 0.
2. *A multiple of the denominator.* The entry denominators of Q have twelve irreducible factors f (`q0 + 1`, `t + 1`,
   `t + 2` and nine factors of degree 2 in each of q0, t, i, σ, one for each node of Lemma E). For each f let
   e_f = max over permutations p of Σ_a mult_f(den Q_{a,p(a)}), and D* = Π f^{e_f}. Every term of the Leibniz
   expansion of det Q has its denominator dividing D*, so T := D* det Q is a polynomial, and for each variable v,
   deg_v T ≤ B_v := deg_v D* + max_p Σ_a (deg_v num Q_{a,p(a)} − deg_v den Q_{a,p(a)}). Both maxima over permutations
   are computed exactly by dynamic programming over subsets. Here B_t = 98 and B_i = 102.
3. *The candidate denominator.* D_c = Π f^{d_f} with the exponents seen on slices (d_f ≤ e_f for every f), so
   R := D*/D_c is a polynomial, of degree 42 in t and 42 in i.
4. *Slices.* For every (y, z) with 0 ≤ y ≤ 98 and 0 ≤ z ≤ 102 (t = 4 + y, i = 7 + z): the entries of Q are specialised
   exactly, det Q is computed by exact symmetric elimination over Q(q0, σ), multiplied by D_c, checked to be a
   polynomial, and written with q0 = 4 + x (`g2slices.py`). That is 10 197 slices.
5. *Interpolation.* For each column y, every coefficient of x^a σ^b is interpolated in z over z = 0..102 (inverse
   Vandermonde over Q); the coefficients of z^l with l > 60 vanish. Then every coefficient of x^a z^l σ^b with l ≤ 60
   is interpolated in y over y = 0..98; the coefficients of y^k with k > 56 vanish (`g2interp.py`). The result is a
   polynomial Ñ(x, y, z, σ) of degree at most 56 in y and 60 in z with Ñ = det Q · D_c on every slice.
6. *Identity.* T − Ñ R is a polynomial of degree at most 98 in y and 102 in z (deg Ñ R ≤ 56 + 42 and 60 + 42), and it
   vanishes on every slice, because there T = det Q · D* = (det Q · D_c) R = Ñ R. A polynomial of degree at most 98 in
   y and 102 in z that vanishes for y = 0..98 and z = 0..102 (as a polynomial in the other two variables) is zero. So
   T = Ñ R and det Q = Ñ/D_c identically.
7. *Signs.* Every coefficient of Ñ is nonnegative and one coefficient free of x, y, z is positive, so Ñ > 0 on 𝓡;
   every factor of D_c is one-signed positive after the shifts, so D_c > 0 on 𝓡. Hence det Q > 0 on 𝓡.
8. *Inertia.* Q is positive definite at (q0, t, i, σ) = (4, 4, 7, 1) (exact symmetric elimination). 𝓡 is connected,
   Q is continuous on 𝓡 and never singular there, so no eigenvalue changes sign: Q ≻ 0 on all of 𝓡, in particular at
   every member and every σ > 0. With step 1 this is Theorem G2. ∎

Steps 1 to 3 and 7 (for D_c) and 8 are checked by `g2final.py`, which also checks that every slice of the grid was
saved, that every column passed step 5 and that the interpolation in y passed; steps 4 and 5 are `run_g2.sh`
(`g2slices.py`, `g2interp.py`).

## 3. The numbers

From `logs/g2/g2_final.log`, `logs/g2/g2_stageB.log` and the worker logs:
- The family form has the 16 twisted vectors of `famlib.py`; the two-row part has 12, and Q (gauge-fixed) is 10 × 10.
  The preliminaries of step 1 hold (27 sanity points at σ = 2/7 and 5; every entry denominator one-signed; the split
  block-diagonal; κ_0, κ_1 annihilated; gauge determinant 1).
- The entry denominators of Q have 12 irreducible factors. The exponents of D* are 8, 4, 8, 8, 4, 3, 2, 4, 4, 2, 6, 1
  and those of D_c are 6, 4, 4, 4, 4, 2, 1, 2, 2, 1, 1, 1, in the order of `g2core.Slicer.factors`: (t + 1) is the
  fifth, (q0 + 1) the sixth and (t + 2) the last, and the others are the nine of degree 2 in each variable.
- The degree bounds of T are 95, 98, 102, 95 in (q0, t, i, σ); R has degrees 43, 42, 42, 42. So the grid is
  y = 0..98, z = 0..102 and the low degrees are 56 in y and 60 in z.
- Every one of the 10 197 slices gives a polynomial det Q · D_c with 1 492 terms after q0 = 4 + x (6 to 24 seconds a
  slice). Eight workers ran for about 5.5 hours each, five taking the columns upwards from y = 0 and three
  downwards from y = 98, a slice already saved being skipped; 181 slices were computed twice where they met.
- One upward worker stopped at the slice (64, 87): another worker was writing the same file at that moment, and the
  rename failed with a Windows file-lock error (`logs/g2/g2_worker_4.log`; the slice itself was saved by the other
  worker). Its remaining columns were done by the downward workers. Every slice file of the grid is present, and
  every column passed stage A, reading all 103 of its slices (`g2final.py`, step 6).
- Stage A gives no coefficient above z^60 in any of the 99 columns (1 492 monomials each). Stage B gives no
  coefficient above y^56. Ñ has 1 117 575 nonzero coefficients, all positive, of degrees at most 52, 54, 56, 52 in
  (x, y, z, σ), and a positive coefficient free of x, y and z. Q is positive definite at (4, 4, 7, 1).

## 4. Theorem MF

By Theorem MF′ of `NOTES_MI.md` the letter-1 block at squared weights (σ, 1, …, 1) is ⪰ 0 with kernel exactly K_j on
every isotypic part of S_{k−1}, for every k, j and σ > 0, except possibly the two-row parts (k − 1 − i, i) with
i ≥ 7, q0 = j − i ≥ 4 and t = k − 1 − j − i ≥ 4. Those are exactly the members of the generic family
[q0 ≥ 4, t ≥ 4, i ≥ 7], and on them the two-row part is Theorem G2 (by Lemma S of `NOTES.md` the vectors
O, G, U, W span the twisted vectors of the two-row part, and K_j ∩ V_i = span(κ_0, κ_1) there since i ≥ 2). So the
statement holds on every isotypic part, which is Theorem MF. The reduction from (σ, 1, …, 1) to any weight vector
whose squared weights are all equal but one is a permutation of the copies and a common rescaling, as in
`marked_copy/NOTES.md`. ∎

## 5. Checks and controls

- *Test points* (`g2points.py`, `logs/g2/g2_points.log`). Stage B also evaluates Ñ at four points far outside the
  grid, (x, y, z, σ) = (3, 131, 157, 2/7), (40, 250, 9, 5), (0, 0, 400, 1/1000) and (17, 1000, 1000, 1000) (so t and i
  up to 1 004 and 1 007). At each, det Q · D_c computed directly from the form (all four variables specialised, a
  10 × 10 rational determinant) equals Ñ exactly (values of order 10¹⁵⁹ to 10²⁷³).
- *Members of the family* (`mfcheck.py`, region g7, in `../run_check_mi.sh`; run before this proof, as a check of
  the open family). At sampled members with 56 ≤ k ≤ 100 and five values of σ the reduced form Q_i is ⪰ 0 with
  exactly the expected nullity (exact rank and elimination).
- *Negative control* (`g2neg.py`, `logs/g2/g2_neg.log`, which ends with FAIL as it should). With B replaced by
  (1 + 1/1000)B the kernel vectors κ_0, κ_1 are no longer annihilated, so the gauge argument of step 1 is not
  available; and on the slice (t, i) = (4, 7) the determinant of the gauge-fixed form, after q0 = 4 + x, has 214
  negative coefficients out of 1 683, and is negative at (q0, σ) = (4, 1000).
- *The vanishing test* (`g2perturb.py`, `logs/g2/g2_perturb.log`, which ends with FAIL as it should). Every one of
  the 99 columns, interpolated in z over 103 values, has no coefficient above z^60. With the single slice (y, z) =
  (0, 51) multiplied by 1 + 1/1000, all 62 664 coefficients above z^60 of that column (1 492 monomials, 42 powers)
  are nonzero.

## Files

| File | What it does |
|---|---|
| `g2lib.py` | the gauge-fixed family form with the preliminaries of step 1; the maxima over permutations (dynamic programming over subsets) |
| `g2core.py` | slices: the four-variable entries specialised in two variables, the exact determinant over the other two |
| `g2slices.py` | step 4 for a set of columns (one worker); then step 5 in z for each finished column (its mode `check` is not used) |
| `g2interp.py` | step 5: stage A (one column, in z) and stage B (all columns, in y; the signs of Ñ; Ñ at four test points) |
| `g2final.py` | steps 1–3, 7 and 8 and the bookkeeping of the run |
| `g2points.py` | the check of Section 5: det Q · D_c computed directly at the four test points against Ñ |
| `g2neg.py` | the negative control of Section 5 |
| `g2perturb.py` | the control of the vanishing test of Section 5 (one slice perturbed) |
| `run_g2.sh`, `run_g2_rev.sh` | the slice workers (five taking the columns from y = 0, three from y = 98; a slice already saved is skipped) and stage B |
| `logs/g2/` | the logs of the run (`g2_worker_*.log`, `g2_rworker_*.log`, `g2_stageB.log`, `run_g2*.out`), of `g2final.py`, `g2points.py`, `g2neg.py` and `g2perturb.py` |
