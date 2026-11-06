# One marked copy: the type I blocks with a letter 1, families symbolic in the sizes

6 November 2026. Nothing here is refereed. Authors: Ansh Mishra and Aryan Senthilkumar.

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
| 8 | Coverage for k ≤ 120 | **exact bookkeeping** (`coverage_st3.py`, Section 5): of the 18 722 836 forms of the task with k ≤ 120, 73 467 by Theorem I1, 9 138 903 in the 347 certified families, 9 510 466 open (all sym; smallest open k = 30); every anti form covered |
| 9 | Negative controls B → (1 + 1/1000)B | **exact**: the certificate fails for 83 of the 89 MT sym families with one symbolic size (74 with ⟨κ, Qκ⟩ < 0 at a member, so not ⪰ 0) and for both generic anti families (Section 7) |
| 10 | The sym forms whose canonical family has at least three symbolic sizes (case MT: 14 + 1 families, and 19 of the families with two) or at least two (case MU: 53 + 14 + 1 families); for k ≤ 29 they are part of Theorem I1 | **open**; measured costs in Section 6 (case MT: about 5 hours of one process for the ten three-size families with i, τ′ symbolic, about 96 hours for the generic family; case MU: about 12 hours for the two-size families, out of reach for the generic one by this route) |

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
| `tlib.py`, `ufield.py`, `signtest.py`, `cert_t.py`, `coverage_t.py`, `flintfield.py` | copied unchanged from the repository (`typeI/`, `middle/`); `guard.sh` copied with the memory cap lowered to 2.0 GB and the PID written to LOG.pid |
| `run_all.sh`, `run_extra.sh` | every step, one Python process at a time under `guard.sh`; `bash run_all.sh > logs/run_all.out`, then `bash run_extra.sh > logs/run_extra.out`, each ending with ALL PASS |
| `st3check.py`, `run_st3check.sh` | a check written separately from the certificates: three sampled members of every certified family with 30 ≤ k ≤ 90 (two seeds; 2 082 members), the numeric reduced form of `tlib.py` ⪰ 0 with exactly the span of its kernel vectors as kernel (exact rank and elimination); all pass |
| `logs/`, `res/` | the logs of the final runs and the certified results (pickles read by the coverage); `logs/dev/` holds development runs and scripts used by no proof, `logs/smoke/` the quick tests of `smoke.sh` |
