# One marked copy: the type I blocks with a letter 1, families symbolic in the sizes

6 November 2026. Nothing here is refereed. Authors: Ansh Mishra and Aryan Senthilkumar.

Stage 4 (this folder): Sections 1 to 7 are stage 3, unchanged except for the status rows 8 to 18, a note at the head of
Section 6, and the files table;
Sections 8 to 11 add the irrep split (Lemma S′), the MT families with two and three symbolic sizes, the MU families with
q0 ≤ 3 and t ≤ 3, the three-row parts of every family, the coverage with every result, the controls, and the measured
costs of what is still open.

This continues `relative-anchors/code/marked_copy/typeI/NOTES.md` of the repository `werner-npt-research` (cited as
*typeI*), whose notation is used throughout: squared weights s_0 = σ > 0 (copy 0 is marked), s_1 = … = s_{k−1} = 1;
T the copies that carry one letter 2; case MT (copy 0 ∈ T) and case MU (copy 0 ∈ U); T′ and U′ the unmarked copies
of T and U, τ′ = |T′|; the T′-modes sym and anti; the irrep index i; the reduced forms Q_{mode,i} on the twisted
vectors (Lemma R′ of typeI) and their known kernel vectors κ (sym mode only). By Lemma R′ the type I block labelled by
T is ⪰ 0 with kernel exactly P_T iff every Q_{mode,i} is ⪰ 0 with kernel exactly the span of its κ (anti: positive
definite).

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Lemma Z (sizes).** A form has four sizes q0 = j − τ − i, t = m − q0, i, τ′ with k = τ′ + q0 + t + 2i + 1; nonzero twisted spaces need q0 ≥ −1, t ≥ −1 (MT) or ≥ −2 (MU); 3 ≤ j ≤ k − 2 is t + i ≥ 2 (MT) or ≥ 1 (MU) together with lower bounds on j | **proved** by hand (Section 1) |
| 2 | **Lemma N (structure).** With ν = τ′ + 2i, for i ≥ 2 every entry of a sym form is a polynomial of degree ≤ 2 in i over Q(q0, t, ν, σ) with denominator free of i, and every entry of an anti form depends on (i, τ′) only through ν | **proved** by hand (Section 1); **exact check** on the four generic forms and six edge forms, and (with the calculus denominators) on all 62 anti families with i, τ′ symbolic |
| 3 | **The formal continuation** `tfamlib.py` and **Lemma V** (validity when every recorded denominator is one-signed after the shifts) | **proved** by hand (Section 2); **exact checks** against `tlib.py`: the sample points of every certificate, and 1 476 far (member, σ) points with 221 121 entries, k up to 67, including all four generic forms |
| 4 | **Theorem AT (anti mode, every k).** For every k, every 3 ≤ j ≤ k − 2, every T with τ′ ≥ 2 and every σ > 0, every Q_{anti,i} is positive definite (the block is positive definite on the (τ′ − 1, 1)-isotypic part) | **proved**, computer-assisted, exact: all 150 anti families (68 MT, 82 MU), including the generic ones (Section 3) |
| 5 | **Theorem S1 (every k).** Every sym form of the task with at least three of q0 ≤ 3, t ≤ 3, i ≤ 1, τ′ ≤ 1 is ⪰ 0 with kernel exactly the span of its κ | **proved**, computer-assisted, exact: 89 MT and 67 MU families with one symbolic size (Section 3) |
| 6 | **Theorem SD (every k).** Every MT sym form of the task with q0 ≤ 3 and t ≤ 3 is ⪰ 0 with kernel exactly span(κ) | **proved**, computer-assisted, exact: the 22 families [MT sym; q0 = a, t = c, i ≥ 2, τ′ ≥ 2] by determinant and interpolation with Lemma N (Proposition D, Section 4; 142 s), with Theorems S1 and I1 |
| 7 | **Theorem S2 (case MT, every k).** Every MT sym form of the task with exactly two large sizes, not i and τ′, and with i = 0 or τ′ = 0, is ⪰ 0 with kernel exactly span(κ) | **proved**, computer-assisted, exact: 19 families by determinant and interpolation in one size (Proposition D1, Section 4b; 355 s) |
| 8 | Coverage for k ≤ 120 (stage 3) | **exact bookkeeping** (`coverage_st3.py`, Section 5): of the 18 722 836 forms of the task with k ≤ 120, 73 467 by Theorem I1, 9 138 903 in the 347 certified families, 9 510 466 open (all sym; smallest open k = 30); every anti form covered |
| 9 | Negative controls B → (1 + 1/1000)B (stage 3) | **exact**: the certificate fails for 83 of the 89 MT sym families with one symbolic size (74 with ⟨κ, Qκ⟩ < 0 at a member, so not ⪰ 0) and for both generic anti families (Section 7) |
| 10 | **Lemma S′ (irrep split).** On a sym form with i ≥ 1 the basis change A, B0, B1 → U = A + B0, W = B0 + B1, S = i·A − (n + 1)·B0 + (m − n)·B1 on each split block gives Q = Q₂ ⊕ Q₃ (Q₃ on the S vectors: the (u′ − i − 1, i, 1) part), with the known kernel in the first part | **proved** by hand (Lemma S of `middle/` with spectators, Section 8.1); in every family where it is used, block-diagonality, the change-of-basis determinant m + 1 + i and the kernel vectors are also **verified identically** |
| 11 | **Theorem S3 (case MT, every k).** Every MT sym form with exactly two large sizes, not i and τ′, whose fixed sizes among i, τ′ equal 1, is ⪰ 0 with kernel exactly span(κ); with S1, SD, S2: **every MT sym form with at most two large sizes** | **proved**, computer-assisted, exact: the 19 families of `lists/MT_sym_2e_large.txt` by Proposition D1 (Section 8.2) |
| 12 | **Theorem MD (case MU, every k).** Every MU sym form with q0 ≤ 3 and t ≤ 3 is ⪰ 0 with kernel exactly span(κ) | **proved**, computer-assisted, exact: the 24 families [MU sym; q0 = a, t = c, i ≥ 2, τ′ ≥ 2] by Proposition D on the split form (Q₃ by elimination; Section 8.3; 16 694 s of worker time, the largest family 48 minutes), with Theorems S1 and I1 |
| 12b | **MU, two large sizes, not i and τ′.** For 11 of the 29 MU sym families with exactly two symbolic sizes, not i and τ′ (`lists/MU_sym_2e_cheap.txt` without [q0 = 3, t ≥ 4, i = 0, τ′ ≥ 2], [q0 ≥ 4, t = 2, i = 0, τ′ ≥ 2], [q0 ≥ 4, t = 3, i = 0, τ′ ≥ 2]), every member is ⪰ 0 with kernel exactly span(κ) | **proved**, computer-assisted, exact: Proposition D1 on the split form (Section 8.5b; 5 938 s); the other 18 are **open** (three stopped by the memory cap in the interpolation, 15 not run) |
| 13 | **Theorem SD2 (case MT, every k).** Every MT sym form with i ≥ 2, τ′ ≥ 2 and exactly one of q0 ≤ 3, t ≤ 3 is ⪰ 0 with kernel exactly span(κ) | **proved**, computer-assisted, exact: the 10 families of `lists/MT_sym_3d.txt` by Proposition D2 (i and one of q0, t sliced, with the split; Section 8.4; 9 712 slices, about 5 hours of worker time), with Theorem I1 |
| 14 | **The other MT forms with three large sizes** (q0, t large and exactly one of i, τ′ large; four families) | [q0 ≥ 4, t ≥ 4, i ≥ 2, τ′ = 0] and [q0 ≥ 4, t ≥ 4, i = 0, τ′ ≥ 2]: **proved**, computer-assisted, exact, by Proposition D3 (two sliced sizes, with the split; Section 8.5); [q0 ≥ 4, t ≥ 4, i ≥ 2, τ′ = 1] and [q0 ≥ 4, t ≥ 4, i = 1, τ′ ≥ 2]: **open** (prepared, about 9 400 slices of about 2.5 s each; deferred, resumable) |
| 15 | **Theorem Q3 (three-row parts, every k).** For every form of the task (every k, 3 ≤ j ≤ k − 2, every T) and every σ > 0 the type I block is positive definite on the isotypic part of every S_{T′} × S_{U′}-irrep ((τ′), (u′ − i − 1, i, 1)), in both cases (in particular in the generic families, where Q₃ has 2 (MT) and 4 (MU) vectors) | **proved**, computer-assisted, exact: the certified families, and Q₃ by elimination in the 83 sym families of `lists/q3_list.txt` (59 with a split block, pivots one-signed, at most 16 752 terms; 221 s), with Lemma S′ verified identically (Section 8.6) |
| 16 | Coverage for k ≤ 120 with every result | **exact bookkeeping** (`coverage_st3.py`, Section 9): for 30 ≤ k ≤ 120, 10 798 344 forms in certified families (57.9 %) and 7 851 025 open (all sym, all with Q₃ certified); smallest open k = 30; every anti form covered |
| 17 | Negative controls B → (1 + 1/1000)B for the determinant routes | **exact** (`negdet.py`, Section 10): on 15 families of the new routes the scaled form has ⟨κ, Qκ⟩ < 0 at the minimal member (all 11 with kernel vectors), and the determinant of its gauge-fixed first part has negative coefficients in 13 and is negative at σ = 10⁻⁶ or 10⁶ in 11; 13 of 15 caught, the other two are one-vector families without kernel vector |
| 18 | What is left | **open** for k ≥ 30 (Theorem I1 covers k ≤ 29): the two-row parts Q₂ of the sym forms in [MT sym; q0 ≥ 4, t ≥ 4, i ≥ 2, τ′ = 1], [MT sym; q0 ≥ 4, t ≥ 4, i = 1, τ′ ≥ 2] and the generic MT family, and in 18 MU families with two symbolic sizes (the 15 of `lists/MU_sym_2e_large.txt` and three of Section 8.5b), the 14 MU families with three and the generic MU family. Measured costs (Section 11): about 6.5 h each for the two MT families, about 74 h and an Ñ of about 5 · 10⁷ terms for the generic MT family; 1 to 4 h each for the 15 MU families with two sizes, about 50 h each with three, about 4 · 10⁴ h for the generic MU family |

Nothing here says that Conjecture 12 of `single-excitation-range/` holds in general, or anything about the NPT (Werner
distillability) problem. Computer-assisted means exact rational arithmetic (python-flint `fmpq`, `fmpq_mpoly`,
`fmpq_mat`) and sign checks of coefficients. No floating point enters a proof.

## 1. The sizes and the structure lemma

**Lemma Z (sizes).** Let (j, τ, case, mode, i) be a form at k copies, u′ = k − τ (MT) or k − τ − 1 (MU),
m = u′ − 2i, q0 = j − τ − i (the Dicke count of the base column) and t = m − q0. Then in both cases

    k = τ′ + q0 + t + 2i + 1,   j = q0 + τ′ + i + [case MT],   k − j = t + i + [case MU].

The twisted space is nonzero only if q0 ≥ −1, t ≥ −1 (MT) or t ≥ −2 (MU), i ≥ 0, τ′ ≥ 0 (MT) or τ′ ≥ 1 (MU), and
τ′ ≥ 2 in the anti mode. The forms of the task, 3 ≤ j ≤ k − 2 and τ ≤ j, are those with
- case MT: q0 + τ′ + i ≥ 2, t + i ≥ 2, q0 + i ≥ 0;
- case MU: q0 + τ′ + i ≥ 3, t + i ≥ 1, q0 + i ≥ 0 (and τ′ ≥ 1).
At j = k − 1 (Corollary M1, not treated here) t + i = 1 in case MT and t + i = 0 in case MU.

*Proof.* Substitution (tlib.Calc: τ′ = τ − [MT], u′ = k − τ − [MU], m = u′ − 2i, q0 = j − τ − i). The ranges follow
from the norms of the twisted vectors (tlib.norm_of): every one is a product of an orbit factor (m, 2i, τ′ or 1) and a
binomial C(m − r, q0 + c) with r ∈ {0, 1} and c ∈ {−2, …, 1}, nonzero iff q0 ≥ −c and t ≥ c + r; in case MT the
smallest t allowed by some vector is −1 (B_{11}), in case MU it is −2 (B_{111}). ∎

So each form has four sizes (q0, t, i, τ′) besides σ, and they are tied to k only through k = τ′ + q0 + t + 2i + 1. A
*family* fixes each size at a value below its symbolic minimum (q0 ∈ {−1, …, 3}, t ∈ {−2, …, 3}, i ∈ {0, 1},
τ′ ∈ {0, 1}) or leaves it symbolic with q0 ≥ 4, t ≥ 4, i ≥ 2, τ′ ≥ 2 (shifted: q0 = 4 + x, t = 4 + y, i = 2 + z,
τ′ = 2 + w). Every form lies in exactly one family, its *canonical* family. A family is written
[case mode; q0 …, t …, i …, τ′ …].

**Lemma N (structure).** Let ν = τ′ + 2i. For a family with i and τ′ symbolic, every entry of the sym form is a
polynomial of degree at most 2 in i over Q(q0, t, ν, σ) with a denominator free of i, and every entry of the anti form
depends on (i, τ′) only through ν. The same holds for every denominator recorded in the calculus (Section 2).

*Proof.* In every frame of the orbit calculus (typeI, Section 3) the Dicke weights normalised by 2^i C(m, q0) are
2^{−pD} C(m − r, q0 + d)/C(m, q0) and the raising and lowering factors are q0 + d + 1 and t − r − d + 1: no i, no τ′.
The energy sector of a key has n = m + 2i + bT Boolean unmarked copies, κ₂ = τ′ − bT unmarked letters 2 and
c = q0 + i + δ letters 1 (bT, δ fixed by the frame), so every quantity of Lemma E′ at the node i′ = i + e is a
polynomial in (q0, t, σ) and κ₂ + 2(n − c) = ν + 2t + bT − 2δ, κ₂ + (3/2)(n − c) + c/2 = ν + (3/2)t + q0/2 + bT/2 − δ:
only ν enters, and the DU eigenvalue (n − c − i′)(c + 1 − i′) = (t + bT − δ − e)(q0 + δ + 1 − e) has no i. For i ≥ 2
all five window nodes satisfy i + e ≥ 0, so the node sets do not depend on i. The orbit sizes are 1, m, 2i, τ′,
m(m − 1), 2im, 2i(2i − 2), mτ′, 2iτ′ and τ′(τ′ − 1) in the sym mode, polynomials of degree ≤ 2 in i at fixed ν, and the
constant 2 in the anti mode. ∎

This is the analogue of the structure lemma of Proposition C of `relative-anchors/README.md`. It is checked exactly by
`structure.py` (step `structure`): on the generic forms of all four (case, mode) and on six forms with q0 or t fixed,
the rewritten entries have denominators free of i and degree in i exactly 2 (sym) and 0 (anti).

Consequence used below (*the anti reduction*): in the anti mode the member (q0, t, i, τ′) with i ≥ 2 has the form of the
member (q0, t, 2, τ′ + 2i − 4), and τ′ + 2i − 4 ≥ τ′ ≥ 2. So a family with i symbolic is certified at i = 2 with τ′
symbolic, after Lemma N has been checked on its form G with i and τ′ symbolic (`tprove.py` does both). In detail: let
H be the form of the family with i = 2 (the same computation as G with i = 2: every condition on i is then true in both).
Every entry of G, and every denominator recorded in its calculus (Section 2), rewritten with τ′ = ν − 2i, is free of i
(`structure.py --dens`, step `structure_anti` of `run_extra.sh`, for all 62 anti families with i and τ′ symbolic). So
at a member x = (q0, t, i, τ′) the recorded denominators of G take the values of those of H at
x′ = (q0, t, 2, τ′ + 2i − 4), which are nonzero (the certificate of H checks Lemma V), hence G is the true form at x
(Lemma V) and G(x) = H(x′), which is positive definite by the certificate of H.

## 2. The formal continuation (`tfamlib.py`)

`tfamlib.py` is the calculus of `tlib.py` with each of q0, t, i, τ′ fixed or symbolic and σ symbolic, over the
rational functions of `flintfield.py` (fmpq_mpoly), as `famlib.py` is for `mlib.py`. The frames, the twisted vectors
(the same code as `tlib.Frame.psi`, `phi`), the operators t_a, Π, j̄^±, σ⁻ and the assembly of the form are those of
`tlib.py`. What changes:
- the normalisation: the numeric form of `tlib.py` at a member equals N0 times the family form there,
  N0 = 2^i C(m, q0 + sh) > 0 (sh = −q0 for a fixed negative q0, sh = t for a fixed negative t, else 0);
- the Dicke weights and the factors q0 + d + 1, t − r − d + 1 are rational functions or polynomials in the sizes; the
  range checks 0 ≤ q0 + d ≤ m − r are made only for fixed sizes;
- the energy sector of a key is (kind, bT, δ) (Lemma N) and Lemma E′ is evaluated with the sizes symbolic;
- the Newton nodes of the window e = −2, …, 2: every condition of genuineness (i + e ≥ 0; e ≤ q0 + δ; e ≤ t + bT − δ;
  and the one-sided versions e ≤ q0 + δ − 1, e ≤ t + bT − δ + 1) is decided exactly when it involves fixed sizes only,
  and taken as true when it involves a symbolic size (it is true for large members); a node is two-sided when both
  versions hold;
- every denominator created on the way (the node values of Lemma E′ used as divisors, the differences of node values in
  the divided differences, the factors of the Dicke weights) is recorded.
The presence of a twisted vector is decided at the minimal member and must agree with the member where every
symbolic size is increased by 100 (every condition of tlib.norm_of is a lower bound on a nonnegative combination of
the sizes, so presence is then the same at every member); the orbits must exist at every member.

**Lemma V (validity).** Suppose every recorded denominator and every entry denominator of the family form is
one-signed after the shifts (all coefficients of one sign and a nonzero part free of the size variables). Then at every
member and every σ > 0 the family form, evaluated there, is the reduced form of `tlib.py` divided by N0.

*Proof.* Fix a member and σ > 0. The one-signed denominators do not vanish there, so the family computation can be
evaluated step by step at the member, and each step is the corresponding step of `tlib.py` except in three places.
(a) *Phantoms.* A key out of the range of a symbolic size has Dicke weight 0 at the member (the weight
2^{−pD} C(m − r, q0 + d)/C(m, q0 + sh), written as a product of linear factors as in `famlib.py`, contains a factor that
vanishes there), and the only maps from such a key back into the range are the Dicke raising and lowering steps, whose
factors q0 + d + 1 and t − r − d + 1 vanish at the boundary. So the part of every vector inside the range evolves as
under the true operators, and phantoms never contribute to a dot product. (b) *Extra nodes.* A node of the window that
is not genuine at this member carries no irrep there, and the part of the vector inside the range lies in the sum of
the eigenspaces of the genuine window nodes (the window argument of typeI, Section 3: a frame vector lives on
i′ ∈ [i − p, i + p + r + b] with p + r + b ≤ 2; `tlib.py` asserts it). So the Newton polynomial through the larger
set of nodes, which agrees with the true function at the genuine nodes, gives the same result; the node values are
distinct because their differences are recorded denominators. (c) *One-sided nodes.* A node that is two-sided for large members but one-sided
at this member is i′ = c on the A side or i′ = n − c + 1 on the B side, where μ = 0, so the two-sided formulas e/Δ and
σa/Δ reduce to 1/a and σ/e; the cross term (−σ/2)/Δ multiplies a vector without a component on that irrep (j̄_U^∓
maps the isotypic part of an irrep that does not exist at the other level to zero). Hence every entry evaluated at the
member equals the entry of `tlib.py` divided by N0. ∎

The certificates check the hypothesis of Lemma V for every family, and compare N0 times the family form with the
exact numeric form of `tlib.py` (which asserts its Newton window) at every sample point: each symbolic size at its
minimum, minimum + 1 and minimum + 3, at σ = 2/7 and 5. `tfamcheck.py` (step `famcheck`) compares far members: for
the four generic families (all four sizes symbolic) and twenty families with one, two or three symbolic sizes (both
cases and modes), every member
whose symbolic sizes lie in q0 ∈ {6, 11, 17}, t ∈ {5, 13, 19}, i ∈ {3, 7, 9}, τ′ ∈ {3, 8, 12} (k ≤ 70), at
σ ∈ {3/11, 7/2, 40}: 1 476 (member, σ) points and 221 121 entries, k up to 67, all equal; the forms are symmetric and
every entry and intermediate denominator is one-signed. (One of the twenty, [MU anti; t = −2], has no twisted vector at
all, in agreement with Lemma Z, so its comparison is empty.)

## 3. The certificates by elimination (`tprove.py`)

For a family: (1) the form over Q(sizes, σ), symmetric; (2) the sanity check against `tlib.py`; (3) the hypothesis of
Lemma V; (4) the known kernel vectors annihilated identically; (5) the gauge: one deleted coordinate per kernel vector
(CM in case MT, O_0 and O_1 in case MU), with a constant nonzero determinant of the kernel vectors there; (6) greedy
symmetric elimination of the gauge-fixed form, each pivot a ratio of two polynomials one-signed after the shifts, of
the same sign. Then at every member and every σ > 0 each leading principal minor of the gauge-fixed form, in the order
taken, is a product of positive pivots, so the form is positive definite there (Sylvester), and Q_{mode,i} ⪰ 0 with
kernel exactly the span of the κ. In the anti mode a family with i symbolic is certified at i = 2 (Section 1).

The families are those of `famlist.py`: every pattern of fixed values and symbolic sizes whose member far out
(symbolic sizes at minimum + 100) is a form of the task with a nonzero twisted space. Counts:

| | 1 symbolic size | 2 | 3 | 4 (generic) |
|---|---|---|---|---|
| MT sym | 89 | 60 (22 with i, τ′ symbolic) | 14 | 1 |
| MT anti | 23 | 33 | 11 | 1 |
| MU sym | 67 | 53 (24 with i, τ′ symbolic) | 14 | 1 |
| MU anti | 31 | 38 | 12 | 1 |

**Lemma L (listing).** For every k ≥ 30, every form of the task (3 ≤ j ≤ k − 2, τ ≤ j) with a nonzero twisted space
lies in a family listed by `famlist.py`.

*Proof.* Let x be such a form and F its canonical family. Each condition of Lemma Z (j ≥ 3, k − j ≥ 2, τ ≤ j, and
τ′ ≥ 1 in case MU) is a lower bound c ≤ 3 on a nonnegative integer combination of the sizes. If the combination has a
positive coefficient on a symbolic size of F, it holds at the far member of F (that size is ≥ 102, the fixed sizes are
≥ −2); otherwise it is the same at x and at the far member. The same holds for the presence of each twisted vector
(Section 2). So the far member of F is a form of the task with a nonzero twisted space, and F is listed. ∎

**Theorem AT (the anti mode, every k).** For every k ≥ 4, every j with 3 ≤ j ≤ k − 2, every T with |T| ≤ j and
τ′ ≥ 2, and every σ > 0, the type I block of (R∞)_j labelled by T, at squared weights (σ, 1, …, 1), is positive
definite on the isotypic part of the S_{T′}-irrep (τ′ − 1, 1); that is, every Q_{anti,i} is positive definite. By
the reduction of typeI, Section 1, the same holds for every weight vector whose squared weights are all equal but one.

*Proof.* For k ≤ 29 this is part of Theorem I1. For k ≥ 30, by Lemma L every such form lies in a listed anti family,
and all 150 listed anti families (68 of case MT, 82 of case MU) are certified (steps `anti_MT`, `anti_MU`). ∎

The anti forms have 3 vectors (case MT: Ca, Da_0, Da_1) or 6 (case MU: Ca_{y0}, Da_{z y0}) and no kernel. In case MT
the 68 families take 1 s in all; the generic one [MT anti; q0 ≥ 4, t ≥ 4, i ≥ 2, τ′ ≥ 2] (certified at i = 2) has
entries of at most 111 terms and pivots of at most 923 terms. In case MU the 82 families take 109 s; the generic one
has entries of at most 741 terms and its largest pivot has 9 604 terms (33 s). Every pivot has numerator and
denominator with positive coefficients only after the shifts and a nonzero part free of the size variables
(`logs/run_anti_MT.log`, `logs/run_anti_MU.log`).

**Theorem S1 (sym forms with at most one large size, every k).** For every k ≥ 4, every j with 3 ≤ j ≤ k − 2, every
σ > 0 and every sym form whose sizes satisfy at least three of q0 ≤ 3, t ≤ 3, i ≤ 1, τ′ ≤ 1, Q_{sym,i} ⪰ 0 with kernel
exactly the span of its known kernel vectors.

*Proof.* If all four hold, k ≤ 1 + 3 + 3 + 2 + 1 = 10 and Theorem I1 applies. Otherwise the canonical family has one
symbolic size; for k ≤ 29 Theorem I1 applies, and for k ≥ 30 Lemma L and the certificates of the 89 + 67 listed
families with one symbolic size (steps `sym1_MT`, `sym1_MU`). ∎

In case MT the 89 families have at most 13 vectors (12 after the gauge CM), entries of at most 62 terms and pivots of
at most 1 945 terms, 73 s in all. In case MU the 67 families have at most 22 vectors (20 after the gauge O_0, O_1);
entries of at most 143 terms, pivots of at most 5 237 terms, 3 522 s in all; the slowest are the families with i
symbolic and all 22 vectors present, [MU sym; q0 = 2 or 3, t = 1, 2, 3, i ≥ 2, τ′ = 1], from 243 s to 508 s each.
Every pivot has positive coefficients only after the shifts (`logs/run_sym1_MT.log`, `logs/run_sym1_MU.log`).


## 4. The determinant route with the structure lemma (`tdet.py`)

With i and τ′ both symbolic the elimination of Section 3 is too slow already in case MT (on
[MT sym; q0 = 1, t = 1, i ≥ 2, τ′ ≥ 2] the seventh of twelve pivots has 3 625 terms after 21 s and the next row
updates did not finish in 100 s). Lemma N gives a cheaper route, as Proposition C of the README does at equal weights.

**Proposition D.** Let F = [case sym; q0 = a, t = c, i ≥ 2, τ′ ≥ 2] with a, c fixed, and Qg its gauge-fixed form
(n × n), for which the preliminaries of Section 3 hold (sanity, Lemma V, kernel vectors, gauge). Suppose:
1. with τ′ = ν − 2i every entry of Qg is p_{ab}/d_{ab} with d_{ab} ∈ Q[ν, σ] free of i and deg_i p_{ab} ≤ 2 (Lemma N);
2. B is the maximum over permutations π of Σ_a deg_i p_{a,π(a)}; for each irreducible factor f of the d_{ab},
   e_f is the maximum over π of Σ_a mult_f(d_{a,π(a)}) (both by dynamic programming over subsets); D_c = Π f^{d_f}
   with d_f ≤ e_f; and for z = 0, …, B the determinant of Qg at i = 2 + z, computed exactly over Q(ν, σ), times D_c,
   is a polynomial P_z;
3. Ñ(z, ν, σ), the interpolant of degree ≤ B in z with Ñ(z, ·) = P_z (inverse Vandermonde over Q), has after
   ν = τ′ + 2i, i = 2 + z, τ′ = 2 + w only nonnegative coefficients and a positive one free of z and w; and every
   factor f with d_f > 0 is, after the same substitution, a polynomial with positive coefficients and a monomial free
   of z and w;
4. Qg is positive definite at (i, τ′, σ) = (2, 2, 1).
Then Qg is positive definite at every member of F and every σ > 0, so Q_{sym,i} ⪰ 0 with kernel exactly span(κ) there.

*Proof.* Put D* = Π f^{e_f}. For every permutation π, D* Π_a Qg_{a,π(a)} = (D*/Π_a d_{a,π(a)}) Π_a p_{a,π(a)} is a
polynomial (the exponent of f in the product of the d's is at most e_f) of degree at most B in i, because D* is free
of i. So T := D* det Qg is a polynomial with deg_i T ≤ B, and R := D*/D_c is a polynomial free of i. Specialising
i = 2 + z does not touch any denominator, so T(2 + z, ·) = R · D_c det Qg(2 + z, ·) = R P_z = R Ñ(z, ·) for
z = 0, …, B. T − RÑ has degree at most B in i and vanishes at B + 1 values of i, so T = RÑ and det Qg = Ñ/D_c
identically. By 3, det Qg > 0 on the real region 𝓡 = {i ≥ 2, τ′ ≥ 2, σ > 0}. Qg is continuous on 𝓡 (its entry
denominators are one-signed), 𝓡 is connected and det Qg never vanishes there, so the inertia of Qg is constant on 𝓡;
by 4 it is positive definite. ∎

`tdet.py prove` checks 1–4 (and the preliminaries, recomputed from scratch). It also computes the slice z = B + 1 and
checks that the interpolant through the B + 2 nodes has no term in z^{B+1}, a consistency check of the bound B that
the proof does not need.

**Theorem SD (case MT with q0 ≤ 3 and t ≤ 3, every k).** For each of the 22 pairs (a, c) of
`lists/MT_sym_2d.txt` Proposition D applies to [MT sym; q0 = a, t = c, i ≥ 2, τ′ ≥ 2] (step `det_MT` of
`run_extra.sh`, 142 s in all). Hence for every k ≥ 4, every j with 3 ≤ j ≤ k − 2 and every σ > 0, every MT sym form
with q0 ≤ 3 and t ≤ 3 is ⪰ 0 with kernel exactly span(κ): by Theorem I1 if i ≤ 1 and τ′ ≤ 1 (then k ≤ 10), by Theorem
S1 if exactly one of i ≥ 2, τ′ ≥ 2 holds, and by Proposition D if both hold (Lemma L puts the form in a listed family).

The runs (`logs/run_det_MT.log`): B, the number of slices computed (B + 1 for the proof and one for the check), the
terms of Ñ(z, ν, σ) and its degrees in (z, ν, σ), the terms after the substitution (all coefficients positive in every
family, with a positive one free of z and w), and the time. In all 22 the coefficient of z^{B+1} vanishes, every factor
of D_c is positive on the region and Qg is positive definite at (2, 2, 1). The gauge-fixed forms have 12 vectors when
all are present.

| (q0, t) | B | slices | Ñ terms | degrees (z, ν, σ) | terms after substitution | time (s) |
|---|---|---|---|---|---|---|
| −1, 1 | 4 | 6 | 829 | (4, 19, 17) | 1 725 | 0.1 |
| −1, 2 | 4 | 6 | 941 | (4, 20, 18) | 1 976 | 0.1 |
| −1, 3 | 4 | 6 | 1 367 | (4, 24, 22) | 3 220 | 0.1 |
| 0, 0 | 8 | 10 | 953 | (6, 17, 17) | 1 424 | 0.1 |
| 0, 1 | 12 | 14 | 9 163 | (10, 43, 41) | 17 860 | 2.4 |
| 0, 2 | 12 | 14 | 10 626 | (10, 46, 44) | 21 500 | 2.7 |
| 0, 3 | 12 | 14 | 12 645 | (10, 50, 48) | 27 072 | 6.3 |
| 1, −1 | 2 | 4 | 8 | (2, 2, 2) | 16 | 0.0 |
| 1, 0 | 10 | 12 | 1 648 | (8, 20, 20) | 2 394 | 0.2 |
| 1, 1 | 14 | 16 | 15 316 | (12, 51, 49) | 29 968 | 7.9 |
| 1, 2 | 14 | 16 | 21 521 | (12, 60, 58) | 46 873 | 11.3 |
| 1, 3 | 14 | 16 | 21 520 | (12, 60, 58) | 46 873 | 18.8 |
| 2, −1 | 2 | 4 | 12 | (2, 2, 2) | 16 | 0.0 |
| 2, 0 | 10 | 12 | 1 834 | (8, 21, 21) | 2 710 | 0.2 |
| 2, 1 | 14 | 16 | 15 952 | (12, 52, 50) | 31 605 | 8.9 |
| 2, 2 | 14 | 16 | 23 040 | (12, 62, 60) | 51 330 | 13.9 |
| 2, 3 | 14 | 16 | 22 274 | (12, 61, 59) | 49 068 | 21.5 |
| 3, −1 | 2 | 4 | 13 | (2, 2, 2) | 16 | 0.0 |
| 3, 0 | 10 | 12 | 1 852 | (8, 21, 21) | 2 710 | 0.2 |
| 3, 1 | 14 | 16 | 15 952 | (12, 52, 50) | 31 605 | 8.7 |
| 3, 2 | 14 | 16 | 23 040 | (12, 62, 60) | 51 330 | 15.2 |
| 3, 3 | 14 | 16 | 22 274 | (12, 61, 59) | 49 068 | 21.7 |

The degree of Ñ in z is B − 2 in every family with B ≥ 8, so the bound B is not sharp but the slices are cheap (0.1 to
1.8 s each). In case MU the same route works but is slower: on [MU sym; q0 = 1, t = 1, i ≥ 2, τ′ ≥ 2] (18 vectors
after the gauge) B = 26 and one slice takes 31 s (development run, Section 6), about 14 minutes per family.

### 4b. Slices in one size (`tdet1.py`)

For the families with two symbolic sizes u, v that are not both i and τ′ there is no polynomial structure, and the
route is that of `relative-anchors/code/middle_j2/f2det.py`.

**Proposition D1.** Let F be a sym family with exactly two symbolic sizes u, v, not {i, τ′}, Qg its gauge-fixed form
(n × n) with the preliminaries of Section 3, f the irreducible factors of its entry denominators,
e_f = max_π Σ_a mult_f(den Qg_{a,π(a)}) and D* = Π f^{e_f}. Let B = deg_v D* + max_π Σ_a (deg_v num − deg_v den) of
Qg_{a,π(a)}. Suppose that for y = 0, …, B the determinant of Qg at v = v_min + y, computed exactly over Q(u, σ), times
D*(v = v_min + y) is a polynomial P_y; that the interpolant Ñ(y, u, σ) of degree ≤ B through these values has, after
u = u_min + x, only nonnegative coefficients and a positive one free of x and y; that every factor f is one-signed
positive after the shifts; and that Qg is positive definite at the minimal member and σ = 1. Then Qg is positive
definite at every member of F and every σ > 0.

*Proof.* T := D* det Qg is a polynomial (each Leibniz term is) with deg_v T ≤ B. Specialising v at the nodes does not
touch the entry denominators (they are one-signed, so nonzero there), so T(v_min + y, ·) = P_y = Ñ(y, ·) for
y = 0, …, B, and T − Ñ (a polynomial of degree ≤ B in v) vanishes identically: det Qg = Ñ/D*, which is positive on the
real region; continuity, connectedness and the base point finish the proof as for Proposition D. ∎

(Taking the full D* as the denominator, instead of a candidate from the slices as in `f2det.py`, avoids two pitfalls
met in development: a factor that depends on v alone is constant on every slice and invisible there, and two factors
can share a factor after the specialisation. The price is a larger Ñ.) `tdet1.py` also computes one more slice and
checks that the interpolant through the B + 2 nodes has no term in y^{B+1}.

**Theorem S2 (case MT, two large sizes with i = 0 or τ′ = 0, every k).** For each of the 19 families of `lists/MT_sym_2e_small.txt` (case MT, two symbolic sizes, not both i and τ′, and i = 0 or τ′ = 0 fixed; at most 10 vectors) Proposition D1 applies (step `det1_MT` of `run_extra.sh`). Hence for every k ≥ 4, every j with 3 ≤ j ≤ k − 2 and every σ > 0, every MT sym form with exactly two large sizes (≥ 4 for q0, t; ≥ 2 for i, τ′), not i and τ′, and with i = 0 or τ′ = 0, is ⪰ 0 with kernel exactly span(κ) (Lemma L puts it in one of these families when k ≥ 30, Theorem I1 covers k ≤ 29). Together with Theorems S1 and SD: every MT sym form with at most two large sizes is covered except those with exactly two large sizes, not i and τ′, whose sizes i and τ′ that are not large all equal 1 (the 19 families of `lists/MT_sym_2e_large.txt`, 13 vectors each, left open in Section 6).

The runs (`logs/run_det1_MT.log`, 355 s in all): the family, its number of vectors, the sliced size v with the
range y = 0..B of the proof (one more slice is computed for the check), the terms of Ñ after the shifts and its degrees
in (y, other size, σ), and the time. In all 19 every coefficient of Ñ is nonnegative with a positive one free of the
sizes, the extra node shows no term in y^{B+1}, every factor of D* is positive after the shifts, and Qg is positive
definite at the minimal member and σ = 1.

| family | vectors | sliced size: 0..B | Ñ terms | degrees | time (s) |
|---|---|---|---|---|---|
| [q0 ≥ 4, t ≥ 4, i = 0, τ′ = 0] | 6 | q0: 0..30 | 12356 | (29, 33, 31) | 0.5 |
| [q0 ≥ 4, t ≥ 4, i = 0, τ′ = 1] | 9 | q0: 0..64 | 104169 | (62, 68, 65) | 25.8 |
| [q0 ≥ 4, t ≥ 4, i = 1, τ′ = 0] | 10 | q0: 0..55 | 73072 | (54, 60, 57) | 22.6 |
| [q0 = 0, t ≥ 4, i = 0, τ′ ≥ 2] | 8 | t: 0..54 | 33300 | (54, 53, 51) | 8.3 |
| [q0 = 1, t ≥ 4, i = 0, τ′ ≥ 2] | 9 | t: 0..72 | 74868 | (72, 71, 69) | 37.9 |
| [q0 = 2, t ≥ 4, i = 0, τ′ ≥ 2] | 9 | t: 0..74 | 80920 | (74, 73, 71) | 43.9 |
| [q0 = 3, t ≥ 4, i = 0, τ′ ≥ 2] | 9 | t: 0..74 | 80920 | (74, 73, 71) | 44.9 |
| [q0 ≥ 4, t = 2, i = 0, τ′ ≥ 2] | 9 | q0: 0..63 | 54236 | (61, 66, 64) | 23.8 |
| [q0 ≥ 4, t = 3, i = 0, τ′ ≥ 2] | 9 | q0: 0..62 | 51897 | (60, 65, 63) | 17.9 |
| [q0 = −1, t ≥ 4, i ≥ 2, τ′ = 0] | 3 | t: 0..13 | 663 | (13, 14, 12) | 0.0 |
| [q0 = 0, t ≥ 4, i ≥ 2, τ′ = 0] | 8 | t: 0..39 | 13510 | (39, 38, 36) | 1.8 |
| [q0 = 1, t ≥ 4, i ≥ 2, τ′ = 0] | 10 | t: 0..60 | 46620 | (60, 60, 57) | 20.9 |
| [q0 = 2, t ≥ 4, i ≥ 2, τ′ = 0] | 10 | t: 0..66 | 60791 | (66, 66, 63) | 32.6 |
| [q0 = 3, t ≥ 4, i ≥ 2, τ′ = 0] | 10 | t: 0..66 | 60791 | (66, 66, 63) | 33.0 |
| [q0 ≥ 4, t = −1, i ≥ 2, τ′ = 0] | 1 | q0: 0..2 | 16 | (2, 3, 2) | 0.0 |
| [q0 ≥ 4, t = 0, i ≥ 2, τ′ = 0] | 6 | q0: 0..17 | 1832 | (17, 19, 18) | 0.1 |
| [q0 ≥ 4, t = 1, i ≥ 2, τ′ = 0] | 10 | q0: 0..44 | 20161 | (41, 46, 43) | 8.0 |
| [q0 ≥ 4, t = 2, i ≥ 2, τ′ = 0] | 10 | q0: 0..52 | 36987 | (51, 57, 54) | 12.6 |
| [q0 ≥ 4, t = 3, i ≥ 2, τ′ = 0] | 10 | q0: 0..57 | 47007 | (56, 62, 59) | 19.7 |

## 5. Coverage for k ≤ 120 (`coverage_st3.py`)

`coverage_st3.py 120 29 30` enumerates every (j, τ, case, mode, i) with 3 ≤ j ≤ k − 1 and 1 ≤ τ ≤ j for
k = 4, …, 120 through its sizes, keeps those with a nonzero twisted space (a closed form of `tlib.norm_of`, compared
with `tlib.norm_of` itself on every form with k ≤ 30), sets aside j = k − 1 (Corollary M1), and assigns each form with
3 ≤ j ≤ k − 2 to Theorem I1 (k ≤ 29), to a certified family containing it, or to the open set. It also checks that every
nonempty form lies in the size ranges of Lemma Z, that its canonical family is listed by `famlist.py` when k ≥ 30
(Lemma L), and that no family reported empty contains a nonempty form.

Results with every certified family (step `coverage_all` of `run_extra.sh`, `logs/run_coverage_all.log`):

| | k ≤ 29 (Theorem I1) | 30 ≤ k ≤ 120, in certified families | 30 ≤ k ≤ 120, open |
|---|---|---|---|
| MT sym | (in the 73 467) | 157 703 | 4 683 764 (canonical family with 2, 3, 4 symbolic sizes: 83 310, 1 768 790, 2 831 664) |
| MT anti | (in the 73 467) | 4 418 101 | 0 |
| MU sym | (in the 73 467) | 5 006 | 4 826 702 (2, 3, 4 symbolic sizes: 214 327, 1 780 711, 2 831 664) |
| MU anti | (in the 73 467) | 4 558 093 | 0 |
| all | 73 467 | 9 138 903 | 9 510 466 |

So for 4 ≤ k ≤ 120 the 18 722 836 forms of the task split into 73 467 covered by Theorem I1, 9 138 903 covered by the
347 certified families (49.0 % of those with k ≥ 30) and 9 510 466 open, all in the sym mode; 69 148 forms with
j = k − 1 were set aside. Per k: at k = 30, 5 343 covered and 5 041 open; at k = 60, 38 653 and 40 136; at k = 120,
296 973 and 308 776. The smallest k with an open form is 30 in both cases. Every check of the script passed.

The step `coverage` of `run_all.sh` (with the anti families and the sym families with one symbolic size only) gives:
for k = 4, …, 120, 73 467 forms of the task with k ≤ 29 (Theorem I1), 8 987 297 with 30 ≤ k ≤ 120 in certified
families, 9 662 072 open, and 69 148 forms with j = k − 1 set aside. Every anti form with 30 ≤ k ≤ 120 is covered
(4 418 101 in case MT, 4 558 093 in case MU); in the sym mode 6 097 (MT) and 5 006 (MU) are covered, and the open ones
have canonical families with 2, 3 and 4 symbolic sizes (MT: 234 916, 1 768 790, 2 831 664; MU: 214 327, 1 780 711,
2 831 664). The closed form of the presence agrees with `tlib.norm_of` on all 141 120 forms with k ≤ 30.

## 6. What is open, with measured costs

*(Stage 3. Several of the families listed here are proved in Section 8; Section 11 gives what is still open after stage 4.)*

For every k, the type I blocks with a letter 1 at one marked copy are now proved on the anti part (Theorem AT) and on
every sym form whose canonical family is certified (Theorems S1, SD and S2). What is not proved for every k is
exactly the sym forms whose canonical family is one of the following (all with at least two symbolic sizes):

- case MT: the 19 families with two symbolic sizes of `lists/MT_sym_2e_large.txt` (not both i and τ′; the sizes i, τ′ that are fixed equal 1; 13 vectors), the 14 families with three symbolic sizes and the generic family [MT sym; q0 ≥ 4, t ≥ 4, i ≥ 2, τ′ ≥ 2];
- case MU: all 53 families with two symbolic sizes (24 with i and τ′ symbolic, 29 others), the 14 with three and the
  generic family [MU sym; q0 ≥ 4, t ≥ 4, i ≥ 2, τ′ ≥ 2].
For k ≤ 29 all of these forms are covered by Theorem I1; the smallest k with an open form is 30 (Section 5).

The measurements (`run_extra.sh`, steps `measure_*`; times of one process on the shared machine):

*Case MT, i and τ′ symbolic (Proposition D extended by one or two more interpolations).* With q0 or t symbolic as
well, the denominators of the entries depend on that size, so the route is that of `../../middle/NOTES_G2.md` (Section 2): the
degree bound B_v of T = D* det Qg in each sliced size v from the factored denominators, interpolation over B_v + 1
nodes in each sliced size (with the full D* as in Section 4b, or a candidate D_c from the slices as in `../../middle/NOTES_G2.md` with the
two pitfalls of Section 4b handled), then the signs. The slices are determinants over Q(ν, σ) at fixed (q0, t, i):

| family | gauge-fixed size | degree bounds of T | slices (B + 1 per sliced size) | one slice | slice time in all (one process) |
|---|---|---|---|---|---|
| [MT sym; q0 = 1, t ≥ 4, i ≥ 2, τ′ ≥ 2] | 12 | t: 106, i: 14 (ν: 105, σ: 103) | 107 · 15 = 1 605 | 1.5 s | 40 min |
| [MT sym; q0 ≥ 4, t = 1, i ≥ 2, τ′ ≥ 2] | 12 | q0: 69, i: 14 (ν: 71, σ: 69) | 70 · 15 = 1 050 | 0.5 s | 9 min |
| [MT sym; q0 ≥ 4, t ≥ 4, i ≥ 2, τ′ ≥ 2] (generic) | 12 | q0: 106, t: 112, i: 14 (ν: 111, σ: 109) | 107 · 113 · 15 = 181 365 | 1.9 s | 96 h |

So the ten MT families with three symbolic sizes that include i and τ′ (q0 or t fixed at −1, …, 3) need about
10 · (9 to 40) min of slices, roughly 5 hours of one process (slices with smaller fixed sizes are cheaper, as in the
table of Section 4), and the generic MT family about 96 hours of one process (12 hours on eight), not counting the
interpolations. These bounds are not sharp (in Theorem SD the degree of Ñ in z is B − 2), so the grids are larger than
the polynomials need, but the proof needs B. The other four MT families with three symbolic sizes ([MT sym; q0 ≥ 4,
t ≥ 4, i ∈ {0, 1}, τ′ ≥ 2] and [MT sym; q0 ≥ 4, t ≥ 4, i ≥ 2, τ′ ∈ {0, 1}]) have no polynomial structure in their third
size; by analogy with the bounds above (about 110 in q0 and in t) they need grids of about 110 · 110 slices over
Q(third size, σ), several hours each (not measured).

*Two symbolic sizes, not both i and τ′ (cases MT and MU).* Proposition D1 (Section 4b): slices in one size, each a
determinant over Q(other size, σ). Measured (`tdet.py measure_plain`; B + 2 slices of the smaller bound):

| family | gauge-fixed size | degree bounds of T | one slice | slices in all |
|---|---|---|---|---|
| [MT sym; q0 ≥ 4, t ≥ 4, i = 1, τ′ = 0] | 9 | q0: 55, t: 60 (σ: 58) | 0.3 s (t fixed) | 57 slices in q0, 22.6 s |
| [MT sym; q0 = 1, t ≥ 4, i ≥ 2, τ′ = 0] | 9 | t: 60, i: 62 (σ: 58) | 0.2 s (i fixed) | 62 slices in t, 20.9 s |
| [MU sym; q0 ≥ 4, t ≥ 4, i = 0, τ′ = 1] | 12 | q0: 118, t: 117 (σ: 120) | 5.4 s (t fixed) | 119 slices in t, about 11 min |
| [MT sym; q0 ≥ 4, t ≥ 4, i = 1, τ′ = 1] | 12 | q0: 94, t: 100 | 1.6 to 2.7 s (q0 fixed) | 96 slices in q0, about 4 min |

The first two MT rows are families of Theorem S2: their last column gives the slices and the time of the certificate in `logs/run_det1_MT.log` (which slices in the size with the smaller bound; `measure_plain` times one slice with the other size fixed). The last row is the
largest of the 19 open MT families of `lists/MT_sym_2e_large.txt` (13 vectors); its numbers are from a development run of
`tdet1.py` (`logs/dev/run_det1_MT_killed.log`) that was stopped after more than 80 of its 96 slices for lack of time, not
because anything failed. At a few minutes each (slices plus the interpolation of larger Ñ) these 19 families need
about one to two hours of one process. The 29 MU families of this kind
need minutes to an hour each by Proposition D1 (119 slices of 5 s for the one measured), about 6 hours in all on one
process. On these families the elimination of Section 3 is slower: in a development run with 25 s of elimination per
family, two of the first three MT families of this kind (9 and 10 vectors) were stopped by the limit after about 50 s
each, where Proposition D1 needs about 20 s.

*Case MU, i and τ′ symbolic (Proposition D).* In development, on [MU sym; q0 = 1, t = 1, i ≥ 2, τ′ ≥ 2] (18 × 18
after the gauge, 14 irreducible denominator factors, exponents of D* up to 12): B = 26, so 27 slices (28 with the
check), the first one 31 s: about 14 minutes per family, about 6 hours for the 24 such families on one process (the
slices are independent). Not run here for lack of time.

*Case MU, generic family.* `tdet.py measure MU n n`: 22 vectors (20 after the gauge), entries of the form up to 12 931
terms (in (q0, t, i, ν, σ) after Lemma N: numerators up to 3 391 terms, denominators up to 2 268), 18 irreducible
denominator factors with exponents of D* up to 14, degree bounds of T: q0 259, t 258, i 28 (ν 262, σ 262); one slice
at (q0, t, i) = (4, 4, 2), a 20 × 20 determinant over Q(ν, σ), takes 152 s (7 315 / 6 903 terms). The grid would have
260 · 259 · 29 ≈ 1.95 · 10⁶ slices, about 8 · 10⁴ hours: out of reach by this route. The 14 MU families with three
symbolic sizes sit between this and the MT numbers (by analogy, not measured: [MU sym; q0 = 1, t ≥ 4, i ≥ 2, τ′ ≥ 2] would need
about 260 · 29 slices of 20 to 30 s, a few days of one process). A different idea is needed for case MU
with three or four large sizes (for instance the irrep split of `middle/splitlib.py`, which in the letter-1 block
separated 16 vectors into 12 + 4).


## 7. Checks and controls

All exact. `run_all.sh` (`logs/run_all.out`, 10 steps, 4 036 s, ALL PASS) and `run_extra.sh` (`logs/run_extra.out`, 12
steps, 748 s, ALL PASS) rerun them.
- *Structure lemma* (`structure.py`, step `structure`): Lemma N on the generic forms of all four (case, mode) and on
  [MT anti; 0, 1, ·, ·], [MU anti; 1, −1, ·, ·], [MT sym; 0, 1, ·, ·], [MU sym; 1, −1, ·, ·], [MT sym; 2, ·, ·, ·] and
  [MU sym; ·, 2, ·, ·] (a dot for a symbolic size): denominators free of i and degree in i exactly 2 (sym) or 0 (anti).
- *The family forms against `tlib.py`*: in every certificate at the sample points (Section 2), and at far members in
  step `famcheck` (1 476 (member, σ) points, 221 121 entries, k up to 67). `tlib.py` itself is the calculus checked
  in typeI against the Lemma B blocks of `fastj1.py` for k ≤ 10 (typeI, Sections 6 and 8).
- *Presence of the twisted vectors*: the closed form used by the coverage equals `tlib.norm_of` on every form with
  k ≤ 30 (step `coverage`).
- *Negative controls*, B → (1 + 1/1000)B in the calculus (steps `neg_sym1_MT`, `neg_anti`): with the same gauge, the
  certificate fails for 83 of the 89 MT sym families with one symbolic size; for 74 of them the scaled form has
  ⟨κ, Qκ⟩ < 0 at the minimal member and σ = 1 (exact), so it is not ⪰ 0 there; the six that still pass are the
  one-vector families [MT sym; q0 = 1, 2, 3, t = −1, i ≥ 2, τ′ = 0, 1] (the single vector B_{11}, no kernel vector:
  the margin there exceeds the change of B, as for the one-vector families of `middle/NOTES_MI.md`). The certificates of
  both generic anti families fail with B scaled: at some step no remaining diagonal entry passes the sign test (the
  smallest is Da_1, 789 terms, in case MT and Da_{10}, 8 650 terms, in case MU). These are failures of the
  certificate, not witnesses of indefiniteness (the anti forms have no kernel vector).
- *Consistency of the determinant routes*: the extra node z = B + 1 (Proposition D, all 22 families) and y = B + 1
  (Proposition D1, all 19 families) gives no new term. The anti reduction to i = 2: step `structure_anti` (62 forms,
  their entries and calculus denominators free of i after τ′ = ν − 2i).
- *Hungarian method*: the degree bounds of the measurements with n > 18 use the Hungarian method in exact integer
  arithmetic; for every matrix with n ≤ 18 met in the measurements it agrees with the dynamic programming.

## 8. The irrep split and the determinant routes with two sliced sizes (stage 4)

Sections 1 to 7 are the previous stage, unchanged (their logs, from the rerun of stage 3 from a clean copy, are in `logs/`). This section and the next three
add the families left open in Section 6 that could be done with the measured costs, with the split below.

### 8.1 The irrep split (Lemma S′, `tsplit.py`, `tsbuild.py`)

**Lemma S′.** Consider a sym form (case MT or MU) with i ≥ 1. Call a *block* the three vectors A_x, B_{x0}, B_{x1}
(case MT) or A_{x y0}, B_{x 0 y0}, B_{x 1 y0} (case MU), and n = q0 − x (MT) or n = q0 − y0 − x (MU) the Dicke count of
its A; the block is *split* when all three are nonzero. On a split block put

    S = i·A − (n + 1)·B0 + (m − n)·B1,     U = A + B0,     W = B0 + B1.

Then S lies in the isotypic part of the S_{U′}-irrep (u′ − i − 1, i, 1), while U, W, the B vectors of the blocks that
are not split, and every twisted vector whose row letter is not on U′ (case MT: C, D_z, CM, DM_z, G; case MU: O_{y0},
C_{y0}, D_{z y0}, G_{y0}) lie in that of (u′ − i, i). Hence in this basis Q = Q₂ ⊕ Q₃, Q₃ on the S vectors (at most 2 in
case MT, 4 in case MU), and the known kernel vectors lie in the first part: κ = U_0 − W_0 + U_1 + C + CM (MT),
κ_{y0} = O_{y0} + U_{0 y0} − W_{0 y0} + U_{1 y0} + C_{y0} (MU) on split blocks.

*Proof.* This is Lemma S of `middle/NOTES.md` with spectators. S_{U′} permutes U′ only; the letters 2 of T′, the row
letter when it sits on T′ or on copy 0, and the column letter y0 at copy 0 are fixed by it. A twisted vector whose row
letter is not on U′ has its U′-part in the permutation module of the subsets of U′ carrying a letter 1 in the column,
which contains only two-row irreps. The states of a block (row letter 1 at r ∈ U′, letter x there in the column, the
other letters fixed) form the permutation module of the pairs (r, S), and the two equivariant maps (r, S) ↦ S and
(r, S) ↦ S ∪ {r} into permutation modules of subsets act on A, B0, B1 exactly as in Lemma S (the norms of
`tlib.norm_of` are in the ratios (m − n)/i : (m − n)/(n + 1) : 1 there). So S spans their common kernel, which is the
three-row part of the block, and U, W span its orthogonal complement. The reduced form is the restriction of an
S_{U′}-equivariant form, so distinct isotypic parts are orthogonal for it (Lemma R′ of typeI). ∎

As in `middle/`, no certificate relies on this argument: `tsbuild.build_split` checks, for every family it is used on,
that the change of basis has determinant m + 1 + i on each split block (one-signed after the shifts), that the block of
TᵀQT between the two parts vanishes identically as rational functions in the sizes and σ, and that the kernel vectors,
written in the new basis, have no S component and are annihilated by Q₂ identically. Then at every member and every
σ > 0, Q ⪰ 0 with kernel exactly span(κ) iff the gauge-fixed first part Q₂g is positive definite and Q₃ is positive
definite. In every family run here Q₃ is certified by the elimination of Section 3 (pivots one-signed after the
shifts), and Q₂g by a determinant route. The gauge deletes CM (case MT) and O_0, O_1 (case MU) as before.

For the slices the order of the rows matters: the determinant on a slice is computed by elimination in the order of
the rows, and putting U, W first is the fastest order measured (`dev/probe_order.py`: on a slice of
[MT sym; q0 = 1, t ≥ 4, i ≥ 2, τ′ ≥ 2], 2.2 s with U, W first, 12.8 s in the natural order, 3.2 to 5.4 s for three
other orders and for a dynamic choice of the smallest pivot).

### 8.2 Theorem S3: the 19 MT families with two symbolic sizes left open in Section 6

**Theorem S3 (case MT, every k).** Every MT sym form of the task with exactly two large sizes, not i and τ′, and with the
fixed one of i, τ′ (or both) equal to 1, is ⪰ 0 with kernel exactly span(κ).

*Proof.* Proposition D1 (Section 4b, `tdet1.py`, without the split) applies to each of the 19 families of
`lists/MT_sym_2e_large.txt` (`run_det1_large.sh`: three workers, `logs/run_det1_MT_large_w*.log`, joined in
`logs/run_det1_MT_large.log`; results `res/det1_MT_large_w*.pkl`). Lemma L puts every such form with k ≥ 30 in one of
them, and Theorem I1 covers k ≤ 29. ∎

In all 19 every coefficient of Ñ is nonnegative with a positive one free of the sizes, the extra node shows no term in
y^{B+1}, every factor of D* is positive after the shifts and Qg is positive definite at the minimal member and σ = 1:

| family | vectors | sliced size: 0..B | Ñ terms | degrees | time (s) |
|---|---|---|---|---|---|
| [q0 ≥ 4, t ≥ 4, i = 1, τ′ = 1] | 13 | q0: 0..94 | 331 485 | (92, 100, 96) | 251.4 |
| [q0 = 1, t ≥ 4, i = 1, τ′ ≥ 2] | 13 | τ′: 0..101 | 201 954 | (100, 102, 98) | 248.7 |
| [q0 ≥ 4, t = 1, i = 1, τ′ ≥ 2] | 13 | q0: 0..69 | 66 816 | (65, 70, 68) | 70.4 |
| [q0 = −1, t ≥ 4, i ≥ 2, τ′ = 1] | 4 | t: 0..23 | 2 898 | (23, 24, 22) | 0.3 |
| [q0 = 2, t ≥ 4, i ≥ 2, τ′ = 1] | 13 | t: 0..106 | 231 391 | (106, 105, 102) | 359.4 |
| [q0 ≥ 4, t = 0, i ≥ 2, τ′ = 1] | 8 | q0: 0..23 | 3 814 | (23, 25, 24) | 0.4 |
| [q0 ≥ 4, t = 3, i ≥ 2, τ′ = 1] | 13 | q0: 0..92 | 170 312 | (90, 97, 94) | 186.8 |
| [q0 = −1, t ≥ 4, i = 1, τ′ ≥ 2] | 4 | t: 0..23 | 2 898 | (23, 24, 22) | 0.2 |
| [q0 = 2, t ≥ 4, i = 1, τ′ ≥ 2] | 13 | τ′: 0..105 | 225 634 | (104, 106, 102) | 295.1 |
| [q0 ≥ 4, t = 2, i = 1, τ′ ≥ 2] | 13 | q0: 0..76 | 97 941 | (74, 80, 78) | 75.9 |
| [q0 = 0, t ≥ 4, i ≥ 2, τ′ = 1] | 11 | t: 0..74 | 80 546 | (74, 72, 70) | 49.0 |
| [q0 = 3, t ≥ 4, i ≥ 2, τ′ = 1] | 13 | t: 0..106 | 231 391 | (106, 105, 102) | 361.0 |
| [q0 ≥ 4, t = 1, i ≥ 2, τ′ = 1] | 13 | q0: 0..67 | 63 835 | (63, 69, 66) | 66.6 |
| [q0 = 0, t ≥ 4, i = 1, τ′ ≥ 2] | 11 | τ′: 0..76 | 90 204 | (75, 77, 73) | 51.8 |
| [q0 = 3, t ≥ 4, i = 1, τ′ ≥ 2] | 13 | τ′: 0..105 | 225 634 | (104, 106, 102) | 304.4 |
| [q0 ≥ 4, t = 3, i = 1, τ′ ≥ 2] | 13 | q0: 0..92 | 165 677 | (90, 96, 94) | 168.3 |
| [q0 = 1, t ≥ 4, i ≥ 2, τ′ = 1] | 13 | t: 0..100 | 195 890 | (100, 99, 96) | 280.6 |
| [q0 ≥ 4, t = −1, i ≥ 2, τ′ = 1] | 1 | q0: 0..2 | 16 | (2, 3, 2) | 0.0 |
| [q0 ≥ 4, t = 2, i ≥ 2, τ′ = 1] | 13 | q0: 0..82 | 124 322 | (80, 87, 84) | 122.1 |

The 19 runs took 2 892 s of worker time in all (three workers, 11:07 to 11:25). The degrees are those of Ñ in (y, the other size, σ) after the shifts, as in Section 4b.

With Theorems S1, SD and S2 this covers every MT sym form with at most two large sizes.

### 8.3 Theorem MD: case MU with q0 ≤ 3 and t ≤ 3 (Proposition D with the split, `tsdet.py D`)

`tsdet.py` replaces the builder of `tdet.py` by `tsbuild.build_split` followed by the elimination of Q₃, and hands
the gauge-fixed first part Q₂g to Proposition D unchanged (Section 4: the structure lemma on Q₂g, which holds since U
and W have constant coefficients; D*, the bound B of deg_i T, the slices i = 2 + z over Q(ν, σ), the interpolation with
one extra node as a check, the signs after ν = τ′ + 2i, i = 2 + z, τ′ = 2 + w, the base point). With the split the
slices are smaller: on [MU sym; q0 = 1, t = 1, i ≥ 2, τ′ ≥ 2] the first part is 15 × 15 after the gauge (18 × 18
before the split) and a slice takes 18 to 21 s instead of 31 s; B = 26 in both.

**Theorem MD (case MU, every k).** Every MU sym form of the task with q0 ≤ 3 and t ≤ 3 is ⪰ 0 with kernel exactly
span(κ).

*Proof.* Only i and τ′ can be large. If at most one of them is (i ≤ 1 or τ′ = 1), the canonical family has at most one
symbolic size, and Theorem S1 (k ≥ 30) or Theorem I1 (k ≤ 29) applies. If i ≥ 2 and τ′ ≥ 2, then for k ≥ 30 the form
lies in [MU sym; q0 = a, t = c, i ≥ 2, τ′ ≥ 2] for one of the 24 pairs (a, c) of `lists/MU_sym_2d.txt` (Lemma L), and
Proposition D on the split form applies to each of these 24 families (`run_sd.sh lists/MU_sym_2d_fams.txt`, logs
`logs/sd_MU_sym_*_n_n.log`, results `res/sd_MU_sym_*_n_n.pkl`); for k ≤ 29 Theorem I1 applies. ∎

In every family Q₃ (0 to 4 vectors) is positive definite by elimination, the coefficient of z^{B+1} vanishes, every
coefficient of Ñ after the substitution is nonnegative with a positive one free of z and w, every factor of D_c is
positive on the region, and Q₂g is positive definite at (i, τ′, σ) = (2, 2, 1):

| family | vectors | parts | Q₂g | B | Ñ terms | degrees (z, ν, σ) | terms after substitution | negative | time (s) | result |
|---|---|---|---|---|---|---|---|---|---|---|
| [q0 = −1, t = 1, i ≥ 2, τ′ ≥ 2] | 3 | 3 + 0 | 3 | 4 | 933 | (4, 20, 18) | 1 976 | 0 | 0.1 | PD |
| [q0 = −1, t = 2, i ≥ 2, τ′ ≥ 2] | 3 | 3 + 0 | 3 | 4 | 1 040 | (4, 21, 19) | 2 250 | 0 | 0.1 | PD |
| [q0 = −1, t = 3, i ≥ 2, τ′ ≥ 2] | 3 | 3 + 0 | 3 | 4 | 1 039 | (4, 21, 19) | 2 250 | 0 | 0.1 | PD |
| [q0 = 0, t = 0, i ≥ 2, τ′ ≥ 2] | 8 | 8 + 0 | 7 | 12 | 6 281 | (10, 36, 35) | 10 976 | 0 | 0.8 | PD |
| [q0 = 0, t = 1, i ≥ 2, τ′ ≥ 2] | 12 | 11 + 1 | 10 | 16 | 21 178 | (14, 57, 54) | 40 358 | 0 | 19.4 | PD |
| [q0 = 0, t = 2, i ≥ 2, τ′ ≥ 2] | 12 | 11 + 1 | 10 | 16 | 24 489 | (14, 61, 58) | 48 730 | 0 | 38.8 | PD |
| [q0 = 0, t = 3, i ≥ 2, τ′ ≥ 2] | 12 | 11 + 1 | 10 | 16 | 23 641 | (14, 60, 57) | 46 539 | 0 | 42.2 | PD |
| [q0 = 1, t = −1, i ≥ 2, τ′ ≥ 2] | 6 | 6 + 0 | 5 | 10 | 1 657 | (8, 20, 22) | 2 382 | 0 | 0.2 | PD |
| [q0 = 1, t = 0, i ≥ 2, τ′ ≥ 2] | 16 | 14 + 2 | 12 | 22 | 27 336 | (18, 58, 59) | 45 509 | 0 | 57.5 | PD |
| [q0 = 1, t = 1, i ≥ 2, τ′ ≥ 2] | 20 | 17 + 3 | 15 | 26 | 68 509 | (22, 83, 82) | 125 171 | 0 | 575.8 | PD |
| [q0 = 1, t = 2, i ≥ 2, τ′ ≥ 2] | 20 | 17 + 3 | 15 | 26 | 75 797 | (22, 87, 86) | 142 649 | 0 | 796.9 | PD |
| [q0 = 1, t = 3, i ≥ 2, τ′ ≥ 2] | 20 | 17 + 3 | 15 | 26 | 87 398 | (22, 93, 92) | 171 766 | 0 | 1246.6 | PD |
| [q0 = 2, t = −2, i ≥ 2, τ′ ≥ 2] | 1 | 1 + 0 | 1 | 2 | 21 | (2, 3, 4) | 30 | 0 | 0.0 | PD |
| [q0 = 2, t = −1, i ≥ 2, τ′ ≥ 2] | 8 | 7 + 1 | 6 | 12 | 2 615 | (10, 23, 26) | 3 703 | 0 | 0.4 | PD |
| [q0 = 2, t = 0, i ≥ 2, τ′ ≥ 2] | 18 | 15 + 3 | 13 | 24 | 35 802 | (20, 63, 65) | 59 108 | 0 | 174.4 | PD |
| [q0 = 2, t = 1, i ≥ 2, τ′ ≥ 2] | 22 | 18 + 4 | 16 | 28 | 91 998 | (24, 92, 92) | 170 980 | 0 | 1404.4 | PD |
| [q0 = 2, t = 2, i ≥ 2, τ′ ≥ 2] | 22 | 18 + 4 | 16 | 28 | 103 028 | (24, 97, 97) | 198 000 | 0 | 2170.8 | PD |
| [q0 = 2, t = 3, i ≥ 2, τ′ ≥ 2] | 22 | 18 + 4 | 16 | 28 | 112 290 | (24, 101, 101) | 221 536 | 0 | 2918.6 | PD |
| [q0 = 3, t = −2, i ≥ 2, τ′ ≥ 2] | 1 | 1 + 0 | 1 | 2 | 20 | (2, 3, 4) | 30 | 0 | 0.0 | PD |
| [q0 = 3, t = −1, i ≥ 2, τ′ ≥ 2] | 8 | 7 + 1 | 6 | 12 | 2 539 | (10, 23, 26) | 3 703 | 0 | 0.5 | PD |
| [q0 = 3, t = 0, i ≥ 2, τ′ ≥ 2] | 18 | 15 + 3 | 13 | 24 | 35 865 | (20, 63, 65) | 59 108 | 0 | 176.8 | PD |
| [q0 = 3, t = 1, i ≥ 2, τ′ ≥ 2] | 22 | 18 + 4 | 16 | 28 | 94 151 | (24, 93, 93) | 176 176 | 0 | 1806.7 | PD |
| [q0 = 3, t = 2, i ≥ 2, τ′ ≥ 2] | 22 | 18 + 4 | 16 | 28 | 105 306 | (24, 98, 98) | 203 721 | 0 | 2413.8 | PD |
| [q0 = 3, t = 3, i ≥ 2, τ′ ≥ 2] | 22 | 18 + 4 | 16 | 28 | 114 668 | (24, 102, 102) | 227 695 | 0 | 2848.6 | PD |

The 24 runs took 16 694 s of worker time (one to three workers, 11:23 to 13:52); the largest Ñ has 114 668 terms
(227 695 after the substitution), and every family has B ≤ 28 and Q₃ with at most 4 vectors.

### 8.4 Proposition D2: i, τ′ and one of q0, t symbolic (`tdet2.py`)

**Proposition D2.** Let F = [case sym; q0, t, i ≥ 2, τ′ ≥ 2] with exactly one of q0, t symbolic, call it v, and let the
preliminaries of `tsbuild.build_split` hold for F, with Q₃ positive definite by elimination. Write the entries of Q₂g
with τ′ = ν − 2i (Lemma N: p/d, d free of i, deg_i p ≤ 2), let D* = Π f^{e_f} (e_f the largest sum of multiplicities of
f along a permutation), B_i the largest sum of deg_i p along a permutation, B_v = deg_v D* + max_π Σ (deg_v num − deg_v den),
and D_c = Π f^{d_f} with d_f ≤ e_f, R = D*/D_c (free of i) and L_v = B_v − deg_v R. Suppose that for every
y = 0, …, B_v and z = 0, …, B_i the determinant of Q₂g at v = v_min + y, i = 2 + z, computed exactly over Q(ν, σ),
times D_c there, is a polynomial P_{y,z}; that the interpolant Ñ(y, z, ν, σ) of these values, of degree ≤ B_i in z and
≤ B_v in y, has no term y^k with k > L_v; that after ν = τ′ + 2i, i = 2 + z, τ′ = 2 + w every coefficient of Ñ is
≥ 0 and one free of y, z, w is > 0; that every factor f with d_f > 0 is positive after the shifts; and that Q₂g is
positive definite at the minimal member and σ = 1. Then at every member of F and every σ > 0, Q ⪰ 0 with kernel exactly
span(κ).

*Proof.* T := D* det Q₂g is a polynomial (each Leibniz term is, as in Proposition D) with deg_i T ≤ B_i (D* is free of
i) and deg_v T ≤ B_v. On the grid, T = R · D_c det Q₂g = R P_{y,z} = R Ñ (the entry denominators are one-signed, so the
specialisation is exact). T − RÑ has degree ≤ B_i in i and ≤ max(B_v, deg_v R + L_v) = B_v in v and vanishes on the
(B_v + 1) × (B_i + 1) grid, identically in (ν, σ); so it is zero and det Q₂g = Ñ/D_c. By the signs det Q₂g > 0 on the
connected real region {v ≥ v_min, i ≥ 2, τ′ ≥ 2, σ > 0}, where Q₂g is continuous; with the base point Q₂g ≻ 0 there.
With Q₃ ≻ 0 and Lemma S′ (checked identically) this is the statement. ∎

`tdet2.py` takes d_f = e_f for the factors free of ν and σ (constant on every slice) and otherwise the largest
multiplicity of the specialised factor in the denominator of the slice determinant on four probe slices away from the
smallest sizes, capped at e_f (a first version that probed the smallest slices saw multiplicities larger than e_f: two
specialised factors can coincide there, as noted in Section 4b; the proof does not depend on how d_f is chosen, only on
the checks above). The slices are spread over worker processes (`run_d2.sh`): a worker claims a column y by creating a
file, computes and saves its slices (skipping saved ones) and interpolates the column in z (stage A); the final step
recomputes the preliminaries, checks that every slice and column is present, interpolates in y (stage B, in chunks of
monomials) and checks the signs.

**Theorem SD2 (case MT, every k).** Every MT sym form of the task with i ≥ 2, τ′ ≥ 2 and exactly one of q0 ≤ 3, t ≤ 3 is
⪰ 0 with kernel exactly span(κ).

*Proof.* For k ≥ 30 it lies in one of the ten families of `lists/MT_sym_3d.txt` (Lemma L), and Proposition D2 applies
to each (`run_d2.sh lists/MT_sym_3d.txt`; logs `logs/d2_*`, results `res/det2_*.pkl`); Theorem I1 covers k ≤ 29. ∎

The runs (`run_d2.sh lists/MT_sym_3d.txt`: two workers from 11:14 and a third from 11:19; logs `logs/d2_prep_*`,
`logs/d2_w*_*`, `logs/d2_final_*`; results `res/det2_*.pkl`). In every family every slice gave a polynomial D_c det,
no coefficient of Ñ above y^{L_v} is nonzero, every coefficient of Ñ after the substitution is nonnegative with a
positive one free of y, z, w, every factor of D_c is positive after the shifts, Q₃ is positive definite by elimination
and Q₂g at the minimal member and σ = 1. The bounds are those of T = D* det Q₂g; L_v = B_v − deg_v R is the degree
in v that Ñ was found to have at most (every coefficient above it vanishes); the worker time is the sum over the
workers of the time spent on the family's slices. The final steps of the first three families were run as soon as
their slices were complete (`logs/d2_final_*.computed.log`); the final loop of `run_d2.sh` then only reports their
marker `final_ok` (a resumed run of `run_d3.sh` would do the same; copies of its two final logs are kept likewise).

| family | Q₂g + Q₃ | bounds | L | slices | Ñ terms after substitution | negative | degrees (y, z, w, σ) | worker time (s) | result |
|---|---|---|---|---|---|---|---|---|---|
| [q0 = −1, t ≥ 4, i ≥ 2, τ′ ≥ 2] | 4 + 0 | t: 25, i: 4 | t: 25 | 26 · 5 = 130 | 30 650 | 0 | (25, 26, 26, 24) | 2 | PD |
| [q0 = 0, t ≥ 4, i ≥ 2, τ′ ≥ 2] | 9 + 1 | t: 79, i: 12 | t: 59 | 80 · 13 = 1 040 | 728 008 | 0 | (58, 56, 56, 54) | 940 | PD |
| [q0 = 1, t ≥ 4, i ≥ 2, τ′ ≥ 2] | 10 + 2 | t: 102, i: 14 | t: 72 | 103 · 15 = 1 545 | 1 524 120 | 0 | (71, 69, 68, 66) | 3 755 | PD |
| [q0 = 2, t ≥ 4, i ≥ 2, τ′ ≥ 2] | 10 + 2 | t: 104, i: 14 | t: 73 | 105 · 15 = 1 575 | 1 604 790 | 0 | (72, 70, 69, 67) | 4 430 | PD |
| [q0 = 3, t ≥ 4, i ≥ 2, τ′ ≥ 2] | 10 + 2 | t: 104, i: 14 | t: 73 | 105 · 15 = 1 575 | 1 604 790 | 0 | (72, 70, 69, 67) | 4 371 | PD |
| [q0 ≥ 4, t = −1, i ≥ 2, τ′ ≥ 2] | 1 + 0 | q0: 2, i: 2 | q0: 2 | 3 · 3 = 9 | 25 | 0 | (2, 3, 2, 2) | 0 | PD |
| [q0 ≥ 4, t = 0, i ≥ 2, τ′ ≥ 2] | 6 + 1 | q0: 22, i: 10 | q0: 18 | 23 · 11 = 253 | 15 426 | 0 | (18, 20, 19, 19) | 2 | PD |
| [q0 ≥ 4, t = 1, i ≥ 2, τ′ ≥ 2] | 10 + 2 | q0: 66, i: 14 | q0: 49 | 67 · 15 = 1 005 | 376 506 | 0 | (45, 48, 47, 45) | 660 | PD |
| [q0 ≥ 4, t = 2, i ≥ 2, τ′ ≥ 2] | 10 + 2 | q0: 81, i: 14 | q0: 56 | 82 · 15 = 1 230 | 740 622 | 0 | (54, 58, 57, 55) | 1 345 | PD |
| [q0 ≥ 4, t = 3, i ≥ 2, τ′ ≥ 2] | 10 + 2 | q0: 89, i: 14 | q0: 55 | 90 · 15 = 1 350 | 695 642 | 0 | (53, 57, 56, 54) | 2 679 | PD |

In all, 9 712 slices and 18 184 s of worker time (11:14 to 12:59 on three workers); the final steps took 0 to 19 s each.
Q₃ has 0, 1 or 2 vectors (S_0, S_1) and its pivots have at most 321 terms. As in Theorem SD, the families with q0 = 2
and q0 = 3 give Ñ with the same number of terms.

### 8.5 Proposition D3: three symbolic sizes without the structure lemma (`tdet3.py`)

The four MT families [MT sym; q0 ≥ 4, t ≥ 4, i ∈ {0, 1}, τ′ ≥ 2] and [MT sym; q0 ≥ 4, t ≥ 4, i ≥ 2, τ′ ∈ {0, 1}] have
three symbolic sizes and no polynomial structure in any of them. **Proposition D3** is Proposition D2 with two sliced
sizes v1, v2 (the pair with the smallest grid (B_v1 + 1)(B_v2 + 1)), the third size u and σ free on the slices,
L_vk = B_vk − deg_vk R for both sliced sizes (stage A in z = v2 − min checks that the terms above z^{L_v2} vanish, stage B
those above y^{L_v1}), and the signs after u = u_min + x; the proof is the same (T − RÑ has degree ≤ B_v1 in v1 and
≤ B_v2 in v2 and vanishes on the grid). One caution in the final step: the shift u = u_min + x mixes the powers of u, so
the sign check is made on complete polynomials in (u, σ) for each pair of powers of y and z (each column is split once
by the power of z, and stage B treats one power of z at a time). Likewise in Proposition D2 the substitution
ν = 6 + 2z + w mixes the powers of z and ν, and stage B treats complete groups of one power of σ.

The runs (`run_d3.sh lists/MT_sym_3e.txt`: two workers; logs `logs/d3_prep_*`, `logs/d3_w*_*`, `logs/d3_final_*`;
results `res/det3_MT_sym_n_n_n_0.pkl`, `res/det3_MT_sym_n_n_0_n.pkl`). In both proved families every slice gave a
polynomial D_c det, no coefficient above z^{L} (stage A) or y^{L} (stage B) is nonzero, every coefficient of Ñ after the
shift is nonnegative with a positive one free of y, z, x, every factor of D_c is positive after the shifts, Q₃ is
positive definite by elimination and Q₂g at the minimal member and σ = 1. The other two families, with 10 + 2 vectors
and grids of about 9 400 slices of about 2.5 s (about 6.5 hours of one process each), were prepared (`meta.pkl`, Q₃
certified, base point checked) but deferred for lack of time (folder files `STOP`; Section 11).

| family | Q₂g + Q₃ | bounds of T | L | slices | Ñ terms after the shift | negative | degrees (y, z, x, σ) | worker time (s) | result |
|---|---|---|---|---|---|---|---|---|---|
| [q0 ≥ 4, t ≥ 4, i ≥ 2, τ′ = 0] | 7 + 2 | t: 58, q0: 54 (free i: 59) | t: 44, q0: 40 | 59 · 55 = 3 245 | 380 778 | 0 | (43, 39, 42, 39) | 1 012 | PD |
| [q0 ≥ 4, t ≥ 4, i = 0, τ′ ≥ 2] | 8 + 0 | t: 74, q0: 70 (free τ′: 74) | t: 57, q0: 53 | 75 · 71 = 5 325 | 1 089 786 | 0 | (57, 51, 56, 54) | 4 409 | PD |
| [q0 ≥ 4, t ≥ 4, i ≥ 2, τ′ = 1] | 10 + 2 | t: 98, q0: 93 (free i: 98) | t: 70, q0: 65 | 99 · 94 = 9 306 | – | – | – | – | deferred |
| [q0 ≥ 4, t ≥ 4, i = 1, τ′ ≥ 2] | 10 + 2 | τ′: 97, q0: 95 (free t: 100) | τ′: 69, q0: 67 | 98 · 96 = 9 408 | – | – | – | – | deferred |

### 8.5b MU families with two symbolic sizes, not i and τ′ (Proposition D1 with the split, `tsdet.py D1`)

Of the 29 MU sym families with two symbolic sizes that are not i and τ′ (`lists/MU_sym_2e.txt`), the 14 with at most
14 vectors (`lists/MU_sym_2e_cheap.txt`: those with i = 0 fixed, and those whose fixed q0 or t is small) were run by
Proposition D1 on the split form (`run_sd.sh lists/MU_sym_2e_cheap.txt`; logs `logs/sd_MU_sym_*.log`, results
`res/sd_MU_sym_*.pkl`); the other 15 (18 to 22 vectors, `lists/MU_sym_2e_large.txt`) are left open (Section 11).

Eleven of the 14 are proved: in each, every coefficient of Ñ after the shifts is nonnegative with a positive one free
of the sizes, the extra node shows no term in y^{B+1}, every factor of D* is positive after the shifts, Q₂g is
positive definite at the minimal member and σ = 1, and Q₃ (one vector, in the three families with a split block:
[MU sym; q0 = 0, t ≥ 4, i = 1, τ′ ≥ 2], [MU sym; q0 = 0, t ≥ 4, i ≥ 2, τ′ = 1] and
[MU sym; q0 ≥ 4, t = −1, i ≥ 2, τ′ = 1]) is positive definite by elimination. The other three,
[MU sym; q0 = 3, t ≥ 4, i = 0, τ′ ≥ 2], [MU sym; q0 ≥ 4, t = 2, i = 0, τ′ ≥ 2] and [MU sym; q0 ≥ 4, t = 3, i = 0,
τ′ ≥ 2] (bounds 141 to 146), are not proved: in the first two every slice was computed and gave a polynomial D* det,
but the interpolation step of `tdet1.py`, which keeps the B + 2 ≈ 143 slice polynomials (the full D* times the
determinant, dense in the free size and σ) in memory, went above the 1.2 GB cap and the guard stopped it; the third,
with the same bounds, was stopped by hand before that step. They need the slice-saving form of Propositions D2 and D3
(slices on disk, interpolation by groups), or D_c in place of D* (Section 11).

| family | vectors | sliced size: 0..B | Ñ terms | degrees | negative | time (s) | result |
|---|---|---|---|---|---|---|---|
| [q0 ≥ 4, t ≥ 4, i = 0, τ′ = 1] | 14 | t: 0..117 | 554 215 | (117, 115, 119) | 0 | 1464.9 | PD |
| [q0 = 2, t ≥ 4, i = 0, τ′ ≥ 2] | 14 | t: 0..135 | 483 470 | (135, 138, 137) | 0 | 2152.9 | PD |
| [q0 ≥ 4, t = 1, i = 0, τ′ ≥ 2] | 14 | q0: 0..115 | 294 935 | (110, 114, 113) | 0 | 1063.4 | PD |
| [q0 = 1, t ≥ 4, i = 0, τ′ ≥ 2] | 13 | t: 0..110 | 270 620 | (110, 113, 111) | 0 | 697.8 | PD |
| [q0 = 0, t ≥ 4, i = 1, τ′ ≥ 2] | 12 | t: 0..96 | 169 634 | (95, 95, 92) | 0 | 292.7 | PD |
| [q0 = 0, t ≥ 4, i ≥ 2, τ′ = 1] | 12 | t: 0..90 | 137 199 | (89, 89, 86) | 0 | 249.7 | PD |
| [q0 = 0, t ≥ 4, i = 0, τ′ ≥ 2] | 8 | t: 0..58 | 42 651 | (58, 59, 56) | 0 | 15.3 | PD |
| [q0 ≥ 4, t = −1, i ≥ 2, τ′ = 1] | 8 | q0: 0..29 | 6 776 | (29, 32, 33) | 0 | 0.9 | PD |
| [q0 = −1, t ≥ 4, i = 1, τ′ ≥ 2] | 3 | t: 0..19 | 1 748 | (19, 20, 18) | 0 | 0.1 | PD |
| [q0 = −1, t ≥ 4, i ≥ 2, τ′ = 1] | 3 | t: 0..19 | 1 748 | (19, 20, 18) | 0 | 0.1 | PD |
| [q0 ≥ 4, t = −2, i ≥ 2, τ′ = 1] | 1 | q0: 0..3 | 30 | (3, 4, 4) | 0 | 0.0 | PD |

The eleven runs took 5 938 s of worker time (12:59 to 13:46); the degrees are those of Ñ in (y, the other size, σ)
after the shifts, as in Section 4b.

### 8.6 The three-row parts (`q3cert.py`)

By Lemma R′ and Lemma S′, for a sym form with a split block the isotypic part of the S_{T′} × S_{U′}-irrep
((τ′), (u′ − i − 1, i, 1)) of the block is carried by the S vectors: the block is positive definite there iff Q₃ is
(each split block contributes one dimension to the multiplicity of (u′ − i − 1, i, 1), and the other vectors none, by
the proof of Lemma S). The anti forms have no split block (their row letter is on T′). For a family whose whole form is
certified, Q ⪰ 0 with kernel exactly span(κ) and κ in the first part give Q₃ ≻ 0. For the families still open,
`q3cert.py` certifies Q₃ directly: the preliminaries and the split of `tsbuild.build_split`, then the elimination of
Section 3 on Q₃ alone.

**Theorem Q3 (every k).** For every k ≥ 4, every j with 3 ≤ j ≤ k − 2, every T and every σ > 0, the type I block of
(R∞)_j labelled by T at squared weights (σ, 1, …, 1) is positive definite on the isotypic part of every irrep
((τ′), (u′ − i − 1, i, 1)) of S_{T′} × S_{U′}, in both cases MT and MU.

*Proof.* For k ≤ 29 it is part of Theorem I1. For k ≥ 30 the form lies in its canonical family F, which is listed
(Lemma L). An anti form has no split block, so there is nothing to prove. A sym form either has no split block (then
the part is zero) or Q₃ ≠ 0. If F is certified as a whole (Theorems S1, SD, S2, S3, MD, SD2 and the families of
Sections 8.5 and 8.5b that are proved), Q ⪰ 0 with kernel exactly span(κ), and κ lies in the first part (Lemma S′),
so Q₃ ≻ 0. Every other listed sym family is in `lists/q3_list.txt` (all sym families of case MT with three or four
symbolic sizes and of case MU with two, three or four), and `q3cert.py` certifies its Q₃ (`logs/run_q3cert.log`): of
the 83 families, 59 have Q₃ positive definite by elimination (one to four pivots, all one-signed after the shifts; the
largest, in the generic MU family, has 16 752 terms) and 24 have no split block, 221 s in all, ending with PASS. ∎

So for every k the type I blocks with a letter 1 at one marked copy are proved positive definite on every three-row
isotypic part of S_{U′}; what is open lies in the two-row parts Q₂ of the families of Section 11.

## 9. Coverage for k ≤ 120 with every result (`coverage_st3.py`)

`coverage_st3.py 120 29 30 res/*.pkl res3/q3.pkl` (`logs/run_coverage_st4.log`, which ends with PASS), with every result of stages 3 and 4 (413 families certified as a whole; the 83 three-row certificates of Section 8.6 cover no form and are only counted):

| | 30 ≤ k ≤ 120, in certified families | 30 ≤ k ≤ 120, open | open, by symbolic sizes of the canonical family | open forms whose Q₃ is certified |
|---|---|---|---|---|
| MT sym | 1 688 971 | 3 152 496 | 3: 320 832, 4: 2 831 664 | 3 152 496 |
| MT anti | 4 418 101 | 0 | – | – |
| MU sym | 133 179 | 4 698 529 | 2: 86 154, 3: 1 780 711, 4: 2 831 664 | 4 698 529 |
| MU anti | 4 558 093 | 0 | – | – |
| all | 10 798 344 | 7 851 025 | | |

So for 4 ≤ k ≤ 120 the forms of the task split into 73 467 covered by Theorem I1 (k ≤ 29), 10 798 344 covered by certified families (57.9 % of those with k ≥ 30; 49.0 % in stage 3) and 7 851 025 open, all in the sym mode; 69 148 forms with j = k − 1 were set aside. Per k:

| k | covered by families | open |
|---|---|---|
| 30 | 7 374 | 3 010 |
| 40 | 16 009 | 8 060 |
| 60 | 48 529 | 30 260 |
| 80 | 108 049 | 75 260 |
| 100 | 202 569 | 151 060 |
| 120 | 340 089 | 265 660 |

The smallest k with an open form is 30 in case MT and 30 in case MU. Every open form lies in a family whose three-row part is certified (Theorem Q3), so what is open is the two-row parts Q₂ only. The checks of the script (the closed form of the presence against `tlib.norm_of` for k ≤ 30, the size ranges of Lemma Z, Lemma L, no family reported empty containing a nonempty form) all passed.

## 10. Checks and controls of Sections 8 and 9

All exact.
- *The split*: in every family of Sections 8.3 to 8.5 `tsbuild.build_split` verified the block-diagonality of TᵀQT
  identically, the change-of-basis determinant m + 1 + i one-signed, and the kernel vectors in the first part and
  annihilated; and the preliminaries of Section 3 (sanity against `tlib.py` at the sample points, Lemma V).
- *Consistency of Propositions D and D1 with the split*: the extra node z = B + 1 gives no new term in all 24 families
  of Theorem MD, and y = B + 1 none in the 11 families of Section 8.5b and the 19 of Theorem S3. In Propositions D2
  and D3 the vanishing of the terms above y^{L} (and z^{L} in D3) is part of the proof, and every slice checked that
  D_c det is a polynomial.
- *The slice determinants factor trivially* (`logs/dev/probe_factor.log`): on three slices of D2 families, D_c det is
  σ² times the cube of a linear factor in ν times one irreducible factor (Section 11 (d)).
- *Negative controls*, B → (1 + 1/1000)B (`negdet.py`, `logs/run_negdet.log`, which ends with FAIL as it should):
  15 families covering every route of Section 8 (four of Theorem S3, four of Theorem SD2, the two of Section 8.5
  and five of Theorem MD; `lists/neg_list.txt`). With B scaled the split stays block-diagonal (it does not depend on
  the scaling). In all 11 families with kernel vectors ⟨κ, Qκ⟩ < 0 at the minimal member and σ = 1 (−0.005 to −0.023),
  so the scaled form is not ⪰ 0 and no certificate can succeed. With the same gauge, the determinant of the gauge-fixed
  first part at the minimal member, as a polynomial in σ, has negative coefficients in 13 of the 15 families and is
  negative at σ = 10⁻⁶ or at σ = 10⁶ in 11 of them (the two families without kernel vectors and with 4 and 3 vectors
  included), so the sign condition of the determinant routes fails there. The two families not caught are the
  one-vector families [MT sym; q0 ≥ 4, t = −1, i ≥ 2, τ′ ≥ 2] and [MU sym; q0 = 2, t = −2, i ≥ 2, τ′ ≥ 2] (no kernel
  vector; the margin exceeds the change of B, as for the one-vector families of Section 7). 13 of 15 caught, 14 s.
- *Edits during the runs.* `tdet2.py` and `tdet3.py` were changed while their slice workers were running, only in the
  final step (stage B by complete groups of one power of σ, resp. of z; the marker `final_ok`) and by the stop-file
  check of the slice loop (each family's worker is a new process, so later families loaded the new file); the
  computation of a slice was never changed, and every final step ran with the files as they are now. One early edit
  of `tdet2.py`, made in place, left a syntax error for about a minute, during which no worker started (no traceback
  in any log); every later edit was made on a copy (`dev/*_edit*.py`) and installed after a syntax check.
- *The order of elimination on slices* (`logs/dev/probe_order.log`) and the letter-1 comparison of Section 11
  (`logs/dev/probe_l1.log`) are development probes, used by no proof.

- *Rerun from a clean copy* (`run_verify_st4.sh`, `logs/verify_st4.out`, 6 November 2026, two worker processes): every
  step the claims of this stage rest on was rerun from a fresh copy, with the stage-3 results taken from the rerun of
  stage 3. Everything passed. That is: the large det-1 families; the 24 families of the first Section 8.5 list; the 10
  families of Proposition D2 and the 2 proved families of Proposition D3; the q₃ certificate and the coverage count;
  and 10 of the 11 proved MU families of Section 8.5b, with the negative controls ending in FAIL as they should. The
  eleventh family, [MU sym; 2, n, 0, n], was stopped by the 1.2 GB memory guard of the rerun after 2 264 s (that
  attempt is kept as `logs/sd_MU_sym_2_n_0_n.killed_1p2GB.log`). Run alone under a 2.2 GB guard (`guard_big22.sh`, the same guard with a larger cap) it passed
  (`logs/sd_MU_sym_2_n_0_n.log`, `logs/rerun_2n0n.out`). The check written separately (`run_st4check.sh`,
  `logs/st4check_21.log`, `logs/st4check_22.log`, `logs/st4check_wide.log`) passed for all three seeds.

## 11. What is open, with measured costs, and how to resume

After stage 4 the following holds for every k (with Theorem I1 for k ≤ 29): the anti mode (Theorem AT); the
three-row parts of every sym form (Theorem Q3); and every sym form whose canonical family is certified as a whole. In
case MT these are the 89 + 60 families with one or two symbolic sizes (Theorems S1, SD, S2, S3) and 12 of the 14 with
three (Theorem SD2 and Section 8.5); in case MU the 67 families with one symbolic size (Theorem S1), the 24 with i and
τ′ symbolic and q0, t fixed (Theorem MD) and 11 of the other 29 with two (Section 8.5b). What is not proved
for every k is exactly the two-row part Q₂ of the sym forms whose canonical family is one of the following:
- case MT: [q0 ≥ 4, t ≥ 4, i ≥ 2, τ′ = 1] and [q0 ≥ 4, t ≥ 4, i = 1, τ′ ≥ 2] (prepared for Proposition D3: about
  9 400 slices of about 2.5 s each, about 6.5 hours of one process each), and the generic family
  [MT sym; q0 ≥ 4, t ≥ 4, i ≥ 2, τ′ ≥ 2] (with the split: Q₂g 10 × 10, bounds (q0, t, i) = (99, 104, 14), slices of
  about 1.7 s with U, W first, 100 · 105 · 15 ≈ 1.6 · 10⁵ slices, about 74 hours of one process, and an Ñ of about 5 · 10⁷
  terms, item (d) below);
- case MU: the three families with two symbolic sizes of Section 8.5b stopped by the memory cap (25 to 45 minutes of
  slices each, once the interpolation keeps the slices on disk), the 15 families with two symbolic sizes of
  `lists/MU_sym_2e_large.txt` (18 to 22 vectors, 1 to 4 hours each by Proposition D1 with the split, about 30 to 40
  hours in all), the 14 with three symbolic sizes (about
  50 hours each: [MU sym; q0 = 1, t ≥ 4, i ≥ 2, τ′ ≥ 2] measured at 4 914 slices of 38 s) and the generic family (about
  4 · 10⁴ hours).
For k ≤ 29 all of these are covered by Theorem I1; the smallest k with an open form is 30 (Section 9).

*Case MU: the ideas for a cheaper certificate, measured* (`tsmeasure.py`, `logs/dev/split_MUgen.log`,
`logs/dev/split_MU_more.log`; times of one process on the shared machine).

(a) *The irrep split* (Lemma S′). It makes Q₃ cheap everywhere, even in the generic family with five variables, and
shrinks the first part, the bounds and the slices, but by factors of 1.5 to 2, not by orders of magnitude:

| family | without the split | with the split |
|---|---|---|
| [MU sym; q0 = 1, t = 1, i ≥ 2, τ′ ≥ 2] (Proposition D) | 18 × 18 after the gauge; B = 26; slice 31 s | Q₂g 15 × 15, Q₃ 3 × 3 (elimination, < 0.1 s); B = 26; slice 18 to 21 s; certificate 576 s |
| [MU sym; q0 ≥ 4, t ≥ 4, i = 1, τ′ = 1] (Proposition D1) | not measured | Q₂g 16 × 16, Q₃ 4 × 4 (elimination 0.3 s, pivots ≤ 976 terms); bounds q0 182, t 183, σ 181; slice 75 s (t fixed); 184 slices, about 3.8 h |
| [MU sym; q0 = 1, t ≥ 4, i ≥ 2, τ′ ≥ 2] (three sizes) | by analogy 260 · 29 slices of 20 to 30 s | Q₂g 15 × 15, Q₃ 3 × 3 (elimination < 0.1 s, pivots ≤ 697 terms); bounds t 181, i 26 (ν 181, σ 179); slice 38 s; 182 · 27 = 4 914 slices, about 52 h |
| [MU sym; q0 ≥ 4, t ≥ 4, i ≥ 2, τ′ ≥ 2] (generic) | 20 × 20; bounds (q0, t, i) = (259, 258, 28); slice 152 s; 1.95 · 10⁶ slices | Q₂g 16 × 16; Q₃ 4 × 4 positive definite by elimination in 43 s (pivots one-signed, at most 16 752 terms); bounds (222, 223, 28), ν 222, σ 221; slice 95 s; 223 · 224 · 29 ≈ 1.45 · 10⁶ slices, about 3.8 · 10⁴ hours |

So for the generic MU family the three-row part (u′ − i − 1, i, 1) is now proved for every k (it is Q₃), but the
two-row part Q₂ is not.

(b) *A Schur complement onto the T′ rows.* The rows without a T′ letter (O, U, W, G) carry the letter-1 two-row form
of `middle/` with every energy raised by κ₂ = τ′: with B = 0 the entries of `tlib.py` and of `mlib.py` (at k′ = k − τ′,
j′ = j − τ′) agree exactly (`dev/probe_l1b.py`). One would hope X ⪰ L₀ (X these rows of the type I form, L₀ the
letter-1 form, which is proved for every k with kernel span(κ′_{y0})), so that only a Schur complement on the six C, D
rows would be left. This fails: X − L₀ is indefinite, with 7 positive and 5 negative eigenvalues (exact inertia from the
characteristic polynomial) at the members (q0, t, i, τ′) = (3, 3, 2, 2), (4, 5, 2, 3), (5, 4, 3, 1) and σ = 1/3, 5,
1/7 (`logs/dev/probe_l1.log`). The reason is that the B-part Σ_{a,b} s_a s_b ⟨t_b ψ_a, Ê⁻¹ t_a ψ_b⟩ pairs t_b ψ_a with
t_a ψ_b; it is not a Gram form, so it does not decrease when Ê increases. Without the comparison, the Schur complement
needs X⁻¹ (12 × 12 in five variables), a problem at least of the size of the generic two-row family of NOTES_G2
(10 197 slices in four variables) with κ₂ = τ′ entering every energy. Not pursued.

(c) *Sharper degree bounds.* On generic slices the denominators have much smaller multiplicities than D* allows: in the
ten families of Theorem SD2, deg_v R = 20 to 34, so the degree of Ñ in the sliced size is at most L_v = B_v − deg_v R,
0.6 to 0.75 of B_v
(Section 8.4), and similarly in Section 8.5. A rigorous bound with D_c in place of D* would shrink each sliced dimension
by that factor; the proofs above still need B_v + 1 nodes, because D_c is only shown to work through the full grid. No
independent proof that D_c det is a polynomial was found. Even with L_v in place of B_v the generic MU family would
need about 10⁶ slices.

(d) *The size of the determinant itself.* For the generic families the determinant is too large to interpolate at all:
in the D2 family [MT sym; q0 = 1, t ≥ 4, i ≥ 2, τ′ ≥ 2], Ñ(y, z, ν, σ) has 757 135 terms (1.5 · 10⁶ after the
substitution), and in the generic family its degree in q0 would be about 50 to 70 (the values of L_v in the families
with t fixed), so the generic MT determinant would have about 5 · 10⁷ terms, and its 1.6 · 10⁵ slices about 4 · 10⁸
coefficients (tens of gigabytes); the generic MU determinant, with degrees about 150 to 220 in q0, t, ν, σ and 28 in i,
would have more than 10⁹. So any route through the full determinant of a generic family is out of reach in memory as well
as in time, and these families need a further exact decomposition of the first part, beyond Lemma S′. None was found
here, and in case MT there is no sign of one: the slice polynomials D_c det of three D2 families factor as σ² times
the cube of a linear factor in ν times one irreducible factor with about 2 200 terms (`dev/probe_factor.py`,
`logs/dev/probe_factor.log`), so the determinant of the first part does not split further on slices.

(e) *Not tried.* A slice is now a determinant over Q(ν, σ), computed by elimination with reduced fractions in two
variables. Fixing σ as well, at about deg_σ + 1 rational values, and interpolating in σ would replace it by
univariate determinants over Q(ν); this might cut the time of a slice several-fold (not measured). It would bring the
MU families with three symbolic sizes within reach of a long parallel run, but not the generic ones (item (d)).

*How to resume.* Every run script of Section 8 is resumable, and nothing was left running. A family whose log ends with
PASS (`run_sd.sh`) or whose folder holds the marker `final_ok` (`tdet2.py`, `tdet3.py`; `REDO_FINAL=1` recomputes) is
skipped; saved slices are skipped; `run_d2.sh` and `run_d3.sh` remove the claims of a previous run before starting
(`JOIN=1` adds workers to a run in progress instead). A file `STOP` in a family folder, or `STOP_ALL` in this folder,
makes the slice workers stop after their current column (the way the runs here were stopped).
- the two deferred MT families of Section 8.5: `rm -f STOP_ALL d3/*/STOP; NW=4 bash run_d3.sh lists/MT_sym_3e.txt >
  logs/run_d3_MT_resume.out` (prep and finished families are skipped);
- the MU families with two symbolic sizes, not i and τ′ (Proposition D1 with the split): `NW=4 bash run_sd.sh
  lists/MU_sym_2e_cheap.txt` reruns only the three of Section 8.5b that are not proved (the 11 proved ones are
  skipped), but `tdet1.py` must first keep its slices on disk and interpolate by groups, as `tdet2.py` does, or it
  will stop at the memory cap again; then `lists/MU_sym_2e_large.txt` (15 families with 18 to 22 vectors, 1 to 4
  hours each; the largest measured, [MU sym; q0 ≥ 4, t ≥ 4, i = 1, τ′ = 1], 184 slices of 75 s), with the same change;
- the generic MT family, the 14 MU families with three symbolic sizes and the generic MU family: out of reach by the
  determinant routes (items (a) to (d) above), with Q₃ proved for all of them (Theorem Q3).

## Files

| File | Content |
|---|---|
| `tfamlib.py` | the formal continuation of `tlib.py` (Section 2) |
| `tprove.py` | the certificates by elimination (Section 3); with `--bscale` the negative control |
| `tdet.py` | the determinant route with the structure lemma (Proposition D, Section 4) and the measurements of Section 6 |
| `tdet1.py` | the determinant route with slices in one size (Proposition D1, Section 4b) |
| `structure.py` | the exact check of Lemma N |
| `famlist.py`, `mklists.py`, `lists/` | the families of Section 3 and their lists |
| `tfamcheck.py` | the family forms against `tlib.py` at far members |
| `coverage_st3.py` | the bookkeeping of Section 5 |
| `summary.py`, `d1table.sh`, `smoke.sh` | a summary of the result pickles; a table of a `tdet1.py` log; quick tests of the code paths of `run_extra.sh` (not part of any proof) |
| `tlib.py`, `ufield.py`, `signtest.py`, `cert_t.py`, `coverage_t.py`, `flintfield.py` | copied unchanged from the repository (`typeI/`, `middle/`); `guard.sh` copied with the memory cap lowered to 1.2 GB (it was 2.0 GB in stage 3) and the PID written to LOG.pid |
| `run_all.sh`, `run_extra.sh` | the steps of Sections 1 to 7 (stage 3), one Python process at a time under `guard.sh`; their logs (rerun from a clean copy) are in `logs/` |
| `tsplit.py`, `tsbuild.py` | Lemma S′: the new basis, the congruence, the kernel vectors in the new basis; the split family form with all preliminaries (Section 8.1) |
| `tsdet.py` | Propositions D and D1 on the first part of the split form, Q₃ by elimination (Section 8.3) |
| `tdet2.py`, `tdet3.py` | Propositions D2 and D3: two sliced sizes, workers that claim columns, stage A per column, stage B and the signs in the final step (Sections 8.4, 8.5) |
| `q3cert.py` | the three-row parts Q₃ of the open families (Section 8.6) |
| `negdet.py` | the negative controls of the determinant routes (Section 10) |
| `tsmeasure.py`, `tables.py` | the measurements of Section 11 and the tables of Section 8 (not part of any proof) |
| `run_det1_large.sh`, `run_sd.sh`, `run_d2.sh`, `run_d3.sh`, `run_st4.sh` | the parallel runs of Section 8 (each worker one Python process under `guard.sh`; resumable, see Section 11); `run_st4.sh` runs them in order with the controls and the coverage |
| `lists/` | the family lists, with the new `MT_sym_3d.txt`, `MT_sym_3e.txt`, `MU_sym_2d_fams.txt`, `MU_sym_2e_cheap.txt`, `MU_sym_2e_large.txt`, `q3_list.txt`, `neg_list.txt` |
| `d2/`, `d3/` | the slices and columns of Propositions D2 and D3, one folder per family (`meta.pkl`, `s_y_z.pkl`, `col_y.pkl`, `tz/` for D3, claim files, the markers `final_ok` and `STOP`); 1.8 GB and 1.0 GB, not in the repository (the runs write them) |
| `res3/`, `claims/` | the three-row certificates (`res3/q3.pkl`, read by the coverage only for counting); the family claims of `run_sd.sh` |
| `st3check.py`, `run_st3check.sh`, `run_st4check.sh` | a check written separately from the certificates: sampled members of every certified family, the numeric reduced form of `tlib.py` ⪰ 0 with exactly the span of its kernel vectors as kernel (exact rank and elimination); stage 3: 2 082 members with 30 ≤ k ≤ 90; stage 4: 1 056 members of the 66 families certified there, 30 ≤ k ≤ 126; all pass |
| `run_verify_st4.sh`, `lists/MU_sym_2e_cheap_proved.txt`, `lists/MT_sym_3e_proved.txt` | the rerun of the steps of Section 8 from a clean copy (Section 10) |
| `logs/`, `res/` | the logs of the reruns from a clean copy (`run_all.sh`, `run_extra.sh` of stage 3 and `run_verify_st4.sh` of stage 4) and the certified results (pickles read by the coverage; not in the repository, the runs write them); `logs/first_run_st4/` the logs of the first run of stage 4, from which the times in the tables of Section 8 are taken; `logs/dev/` development runs used by no proof; `dev/` development probes (`probe_l1.py`, `probe_l1b.py` with a copy of `middle/mlib.py`, `probe_order.py`, `probe_factor.py`), the helpers that wrote tables into these notes (`fill_notes.py`, `coverage_text.py`, `table_d2.txt`) and removed absolute paths from tracebacks in three logs (a helper not included here), the edited copies of the scripts that were checked before replacing the files in use, and `orchestrate*.sh`, the helpers that handed the process slots freed by one run to the next during the session (none of these is part of a proof) |
