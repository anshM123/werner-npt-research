# The letter-1 block with one marked copy: the parts of irrep index i ≥ 5

5–6 November 2026. Nothing here is refereed.

## Status

Squared weights s_0 = σ > 0 (copy 0 is marked), s_1 = … = s_{k−1} = 1; u = k − 1; the letter-1 block of (R∞)_j; K_j
its known kernel. For an irrep index i: q0 = j − i, t = k − 1 − j − i, m = q0 + t; V_i the (H_i, χ_i)-twisted
invariants and Q_i the reduced form on V_i. The parts of index i are the isotypic parts of the S_u-irreps (u − i, i)
(the *two-row part*) and (u − i − 1, i, 1) (the *three-row part*). A family [q0 ≥ a, t ≥ c, i ≥ b] is the set of
members (k, j, i) with these sizes; symbolic sizes are shifted to their minimum (q0 = a + x, t = c + y, i = b + z).

| # | Statement | Status |
|---|---|---|
| 1 | **Theorem G3 (the generic family, three-row part).** For every member of [q0 ≥ 4, t ≥ 4, i ≥ 5] and every σ > 0, the form Q₃ of Q_i on the S vectors (the three-row part) is positive definite | **proved**, computer-assisted, exact: 4 pivots, each a ratio of polynomials in (x, y, z, σ) with positive coefficients only (largest 4 533 terms); 11 s (Section 3) |
| 2 | **Theorem E (the edge families, i symbolic).** For a = −1, …, 3 and every member of [q0 = a, t ≥ 4, i ≥ 6], for c = −2, …, 3 and every member of [q0 ≥ 4, t = c, i ≥ 5], and every σ > 0: Q_i ⪰ 0 with kernel exactly K_j ∩ V_i (the span of the κ_{y0} that are nonzero) | **proved**, computer-assisted, exact: 11 families, both parts, three variables each, pivots with positive coefficients only (largest 23 137 terms); 29 minutes (Section 4) |
| 3 | **Theorems N5 and N6.** For b = 5, 6, every member of [q0 ≥ 4, t ≥ 4, i = b] and every σ > 0, Q_b ⪰ 0 with kernel exactly span(κ_0, κ_1) = K_j ∩ V_b | **proved**, computer-assisted, exact: 4 + 10 pivots each, with positive coefficients after the shifts (largest 40 986 terms); 24 and 18 minutes (Section 5; N6 in the extension `run_extra_mi.sh`) |
| 4 | **Corollary T (all three-row parts, every k).** For every k ≥ 2, every j = −1, …, k and every σ > 0, on the isotypic parts of all S_{k−1}-irreps (k − 2 − i, i, 1) the letter-1 block of (R∞)_j is ⪰ 0 with kernel exactly the part of K_j there (for i ≥ 2 that part is 0, so the block is positive definite there) | **proved**, given rows 1 and 2, Lemma S and Theorems MK, MB, MS, F (Section 6); bookkeeping checked for k ≤ 120 |
| 5 | **Corollary P6 (two-row parts with i ≤ 6, q0 ≤ 3 or t ≤ 3, every k).** For every k ≥ 2, every j and every σ > 0, on the isotypic part of (k − 1 − i, i) the block is ⪰ 0 with kernel exactly the part of K_j there whenever i ≤ 6, or q0 = j − i ≤ 3, or t = k − 1 − j − i ≤ 3 | **proved**, given rows 2 and 3 and Theorems MK, MB, MS, F (Section 6); bookkeeping checked for k ≤ 120 |
| 6 | **Theorem M″ (fixed k = 41, …, 55).** For 41 ≤ k ≤ 55, every j and every σ > 0, the letter-1 block of (R∞)_j is ⪰ 0 with kernel exactly K_j (and the same for every weight vector on these k with all squared weights but at most one equal) | **proved**, computer-assisted, exact: the fixed-k certificates of `cert_fixed.certify` for the 1 775 forms with 41 ≤ k ≤ 50 and the 1 282 forms with 51 ≤ k ≤ 55 outside the proved families (11 minutes), plus the families (Section 7) |
| 7 | **Theorem MF′ (what is now proved).** For every k ≥ 2, every j and every σ > 0, the letter-1 block of (R∞)_j is ⪰ 0 with kernel exactly K_j on every isotypic part except possibly the two-row parts (k − 1 − i, i) with i ≥ 7, j − i ≥ 4 and k − 1 − j − i ≥ 4 (these occur only for k ≥ 23); for k ≤ 55 also on those, so for every k ≤ 55 and every j the whole statement holds | **proved**, given rows 1–6 (Section 6) |
| 8 | Checks: the family forms against `mlib.py` at members far from the sample points (k up to 76); exact eliminations of 1 401 sampled members with 41 ≤ k ≤ 60 of the generic family and the edge families (Q_i ⪰ 0 with kernel exactly K_j ∩ V_i at five values of σ, including the two-row part of the generic family, which is not proved); coverage of every (k, j, i) with k ≤ 120 | **exact checks**, not proofs (Section 8) |
| 9 | Negative controls, B → (1 + 1/1000)B: the certificates of Theorems G3, N5 and N6 and of 10 of the 11 edge families fail (the one-vector family [q0 ≥ 4, t = −2, i ≥ 5] still passes); at the sampled members the scaled form is not positive definite at σ = 10⁻⁶ or 10⁶ (all 828 members of the generic family, 512 of the 573 edge members) | **exact** (Section 8) |
| 10 | The two-row parts (k − 1 − i, i) with i ≥ 7, q0 ≥ 4, t ≥ 4 and k ≥ 56 (so "Theorem MF", every letter-1 block of every (R∞)_j for every k at these weights, is **not** proved); the generic family on the two-row part was attempted and stopped at its fifth pivot (Section 9); two or more marked copies; general weights; the type I blocks with a letter 1 at a marked copy (stage 2) | **open** |

Nothing here says that Conjecture 12 of `single-excitation-range/` holds in general, or anything about the NPT (Werner
distillability) problem. Computer-assisted means exact rational arithmetic (python-flint `fmpq`, `fmpq_poly`,
`fmpq_mpoly`) and sign checks of coefficients. No floating point enters a proof.

These notes continue `NOTES.md` of this folder (the middle band, `relative-anchors/code/marked_copy/middle/`); read
it for Lemma S and the family certificate with the split. The scripts listed in Section Files are new; the others and
the saved results (`logs/*.pkl`, `logs/frozen/`) are those of `NOTES.md`, unchanged. Theorems MK, MB, MS and
M′ below are those of the README (Section 4b) and of `NOTES.md`; Theorem F is that of
`relative-anchors/code/marked_copy/NOTES.md`.

## 1. What was open, and what is done here

Proved before: the statement (the letter-1 block of (R∞)_j with squared weights (σ, 1, …, 1) is ⪰ 0 with kernel
exactly K_j) for k ≤ 40 and every j (Theorems MK, M′), for every k when j ≤ 8 or j ≥ k − 8 (Theorem MB), for every k
and j on the parts i ≤ 4 (Theorem MS, from the families [q0 ≥ 4, t ≥ 4, i = b], b ≤ 4) and on the parts with both
q0 ≤ 3 and t ≤ 3 (Theorem F). Open: k ≥ 41, 9 ≤ j ≤ k − 9, i ≥ 5, except the parts with both q0 ≤ 3 and t ≤ 3.
Finishing it needs the generic family [q0 ≥ 4, t ≥ 4, i ≥ 5] (four variables x, y, z, σ) on both parts, and the edge
families [q0 = a, t ≥ 4, i ≥ 6] (a ≤ 3) and [q0 ≥ 4, t = c, i ≥ 5] (c ≤ 3).

Done here:
- the three-row part of the generic family, for every i ≥ 5 (Theorem G3, Section 3);
- all eleven edge families, with i symbolic, on the whole form (Theorem E, Section 4);
- the families [q0 ≥ 4, t ≥ 4, i = 5] and [q0 ≥ 4, t ≥ 4, i = 6] (Theorems N5, N6, Section 5);
- hence every three-row part for every k (Corollary T), and every two-row part with i ≤ 6 or q0 ≤ 3 or t ≤ 3 for
  every k (Corollary P6, Section 6);
- fixed-k certificates for 41 ≤ k ≤ 55 of every form outside the proved families (Theorem M″, Section 7).

`run_all_mi.sh` does all of this except Theorem N6 and the certificates for 51 ≤ k ≤ 55, which are the extension
`run_extra_mi.sh` (run after it, as the `run_extra*.sh` of `NOTES.md`).

Not done: the two-row part of the generic family. What is missing is exactly the two-row parts with i ≥ 7, q0 ≥ 4,
t ≥ 4 and k ≥ 56 (Section 9).

## 2. Family certificates with chosen lower bounds (`prove_gen.py`)

`prove_gen.py` is the certificate of Section 3 of `NOTES.md` (`prove_split.py` with the gauge that deletes
U_{1,y0}), with two changes.

*Lower bounds.* A symbolic size has a chosen minimum: q0 ≥ 4, t ≥ 4, and i ≥ 5 or i ≥ 6. The minimums are written into
the shared dictionary `famlib.SMIN` before anything is computed, so they govern (a) which twisted vectors are present
(decided at the minimal point; presence is monotone in each size), (b) the sanity points (each symbolic size at its
minimum, minimum + 1 and minimum + 3, each at σ ∈ {2/7, 5}; 27 points with three symbolic sizes), and (c) the shifts of
the sign test (q0 = 4 + x, t = 4 + y, i = 5 + z or 6 + z). The formal continuation of `famlib.py` (Section 5 of
`relative-anchors/code/marked_copy/NOTES.md`) is valid whenever the symbolic sizes satisfy q0 ≥ 4, t ≥ 4, i ≥ 2: every
condition that involves a symbolic size (the genuine Johnson nodes of the window and whether a node lies on one side or
on both, the presence of vectors, the orbits that exist) then holds, and the denominators at members (q0 + sh + s,
t − sh + s, m − s, differences of distinct node values, and the Δ, a, e of Lemma E) are nonzero. A larger minimum of i
keeps all of this, so on each family the symbolic form specialises to the true reduced form at every member and every
σ > 0. As before, the certificate checks that every entry denominator is one-signed after the shifts, and compares
N0 times the form with `mlib.py` at the sanity points (exactly; `mlib.py` asserts its Newton window).

*Parts.* With `--parts=3` only the three-row part Q₃ (the S vectors) is eliminated, and the conclusion is only that Q₃
is positive definite on the family. By default Q₃ and then the gauge-fixed Q₂ are eliminated, as before.

The other steps are unchanged: the change of basis T of `splitlib.py` (determinant m + 1 + i per split block,
one-signed), the identity Tᵀ Q T = Q₂ ⊕ Q₃ (the off-diagonal block vanishes identically in the sizes and σ), the known
kernel vectors written in the new basis, lying in the first part and annihilated identically, the gauge deleting
U_{1,y0} (or A_{1,y0} if that block is not split, or else O_{y0}) with a one-signed determinant (it is 1 in every family
here), and greedy symmetric elimination in which every pivot is a ratio of two polynomials whose coefficients, after the
shifts, have one common sign and whose part free of the size variables is nonzero. (In every case here the signs are
positive: denominators are stored with a positive leading coefficient in lex order, which the shifts preserve.) *Why
this is a proof* is as in Section 3 of `NOTES.md`: at every member and every σ > 0 each leading principal minor
of a symmetric permutation of an eliminated part is a product of positive pivots, so the part is positive definite
(Sylvester); with the gauge and Tᵀ Q T = Q₂ ⊕ Q₃, Q_i ⪰ 0 with kernel exactly the span of the known kernel vectors.

## 3. The three-row part of the generic family (`prove_gen.py logs/gen3.pkl n,n,n --imin=5 --parts=3`)

**Theorem G3.** For every k, j, i with q0 = j − i ≥ 4, t = k − 1 − j − i ≥ 4 and i ≥ 5, and every σ > 0, the form Q₃
of Q_i on S_{00}, S_{01}, S_{10}, S_{11} is positive definite.

*Proof (computer-assisted).* The certificate of Section 2 with q0 = 4 + x, t = 4 + y, i = 5 + z and σ symbolic. All
sixteen twisted vectors are present and all four blocks split (parts 12 + 4). The form over Q(x, y, z, σ) has entries
of at most 4 220 terms (after the congruence: 4 849 in the two-row part, 2 349 in the three-row part); the 27 sanity
points (q0, t ∈ {4, 5, 7}, i ∈ {5, 6, 8}) agree with `mlib.py` at σ = 2/7 and 5; every entry denominator is one-signed;
the 12 × 4 off-diagonal block vanishes identically. The pivots of Q₃, in the order taken, with the terms of numerator /
denominator: S_{10} 1 171/809, S_{01} 4 097/3 109, S_{00} 4 135/3 110, S_{11} 4 533/3 595. Each has numerator and
denominator with positive coefficients only and a nonzero part free of x, y, z. 11 seconds (`logs/run_gen3.log`). ∎

By Lemma S of `NOTES.md` the S vectors span the twisted vectors of the three-row part, so Theorem G3 says that
for these members the letter-1 block is positive definite on the isotypic part of (u − i − 1, i, 1).

## 4. The edge families with i symbolic (`prove_gen.py logs/edges_q.pkl ... --imin=6`, `logs/edges_t.pkl ... --imin=5`)

**Theorem E.** Let a ∈ {−1, 0, 1, 2, 3} and c ∈ {−2, −1, 0, 1, 2, 3}. For every member of [q0 = a, t ≥ 4, i ≥ 6], for
every member of [q0 ≥ 4, t = c, i ≥ 5], and every σ > 0, Q_i ⪰ 0 with kernel exactly the span of those
κ_{y0} = O_{y0} + A_{0,y0} + A_{1,y0} + B_{10,y0} − B_{01,y0} that are nonzero (κ_{y0} is nonzero iff O_{y0} is
present). This is K_j ∩ V_i (i ≥ 2: the other known kernel vectors live at i ≤ 1).

*Proof (computer-assisted).* The certificate of Section 2 on each of the eleven families, with three variables (the
symbolic one of q0 and t, then z and σ), both parts eliminated. For each family: the vectors present, the split blocks
(x, y0), the pivots in the order taken (first the three-row part, then the gauge-fixed two-row part) with the terms of
numerator / denominator, the dimension of the kernel, and the time (`logs/run_edges_q.log`, `logs/run_edges_t.log`):

- [q0 = −1, t ≥ 4, i ≥ 6]: 2 vectors, no split block; pivots G_0 9/10, B_{000} 215/208; kernel dimension 0; 0.0 s
- [q0 = 0, t ≥ 4, i ≥ 6]: 8 vectors, split blocks [(0, 0)]; pivots S_{00} 33/19, G_0 9/10, G_1 3/3, B_{001} 285/285, B_{100} 1131/1129, W_{00} 1965/1756, U_{00} 4011/4379; kernel dimension 1; 0.5 s
- [q0 = 1, t ≥ 4, i ≥ 6]: 14 vectors, split blocks [(0, 0), (0, 1), (1, 0)]; pivots S_{01} 55/35, S_{10} 157/119, S_{00} 212/156, O_0 8/10, G_0 160/162, G_1 160/160, U_{01} 3195/2920, B_{101} 4867/5323, W_{01} 5418/6404, W_{00} 8355/7731, W_{10} 8968/8982, U_{00} 14099/14946; kernel dimension 2; 133.5 s
- [q0 = 2, t ≥ 4, i ≥ 6]: 16 vectors, split blocks [(0, 0), (0, 1), (1, 0), (1, 1)]; pivots S_{01} 118/83, S_{11} 359/348, S_{10} 443/354, S_{00} 449/355, O_0 8/10, O_1 5/13, G_0 115/119, G_1 114/115, W_{00} 4869/4442, W_{10} 7578/7596, W_{01} 10354/10354, W_{11} 11012/11674, U_{01} 19276/20180, U_{00} 20646/21682; kernel dimension 2; 424.9 s
- [q0 = 3, t ≥ 4, i ≥ 6]: 16 vectors, split blocks [(0, 0), (0, 1), (1, 0), (1, 1)]; pivots S_{10} 215/164, S_{01} 552/450, S_{00} 555/450, S_{11} 554/543, O_0 8/10, O_1 5/13, G_0 116/119, G_1 116/116, W_{00} 4883/4457, W_{10} 7592/7610, W_{01} 10367/10368, W_{11} 11026/11687, U_{01} 21830/22821, U_{00} 21982/23137; kernel dimension 2; 477.4 s
- [q0 ≥ 4, t = −2, i ≥ 5]: 1 vector, no split block; pivots B_{111} 20/20; kernel dimension 0; 0.0 s
- [q0 ≥ 4, t = −1, i ≥ 5]: 6 vectors, split blocks [(1, 1)]; pivots S_{11} 44/19, O_1 3/4, B_{110} 155/154, B_{011} 355/350, W_{11} 627/546; kernel dimension 1; 0.1 s
- [q0 ≥ 4, t = 0, i ≥ 5]: 13 vectors, split blocks [(0, 1), (1, 0), (1, 1)]; pivots S_{10} 79/55, S_{01} 213/163, S_{11} 247/163, G_1 1/1, O_0 22/25, O_1 15/22, W_{11} 1096/948, W_{01} 2189/1961, U_{01} 3200/3138, B_{010} 4446/4425, W_{10} 4404/4446; kernel dimension 2; 6.5 s
- [q0 ≥ 4, t = 1, i ≥ 5]: 16 vectors, split blocks [(0, 0), (0, 1), (1, 0), (1, 1)]; pivots S_{01} 83/56, S_{11} 348/220, S_{10} 348/284, S_{00} 347/273, G_1 2/4, G_0 1/4, O_0 82/84, O_1 75/82, W_{10} 3586/3248, W_{11} 5262/4886, W_{01} 6829/6339, W_{00} 7331/7334, U_{00} 12982/13656, U_{01} 14009/13931; kernel dimension 2; 148.1 s
- [q0 ≥ 4, t = 2, i ≥ 5]: 16 vectors, split blocks [(0, 0), (0, 1), (1, 0), (1, 1)]; pivots S_{10} 215/165, S_{00} 450/362, S_{01} 450/362, S_{11} 526/362, O_0 8/10, O_1 5/8, G_1 117/119, G_0 102/152, W_{10} 3498/3171, W_{11} 5654/5280, W_{01} 7960/7410, W_{00} 7843/7844, U_{01} 15637/15482, U_{00} 15896/16844; kernel dimension 2; 209.6 s
- [q0 ≥ 4, t = 3, i ≥ 5]: 16 vectors, split blocks [(0, 0), (0, 1), (1, 0), (1, 1)]; pivots S_{10} 215/165, S_{00} 555/453, S_{01} 555/453, S_{11} 644/453, O_0 8/10, O_1 5/8, G_1 117/119, G_0 108/152, W_{10} 3923/3569, W_{11} 6233/5831, W_{01} 9411/8791, W_{00} 9286/9287, U_{00} 19153/20068, U_{01} 19251/19165; kernel dimension 2; 355.5 s

Every pivot has numerator and denominator with positive coefficients only and a nonzero part free of the size
variables; every entry denominator is one-signed; the off-diagonal block vanishes identically; the known kernel vectors
are annihilated identically; the gauge determinant is 1. ∎

## 5. The families [q0 ≥ 4, t ≥ 4, i = 5] and [q0 ≥ 4, t ≥ 4, i = 6] (`prove_gen.py logs/mid_i5.pkl n,n,5`, `prove_gen.py logs/mid_i6.pkl n,n,6`)

**Theorem N5.** For every k, j with q0 = j − 5 ≥ 4 and t = k − 6 − j ≥ 4, and every σ > 0, Q_5 ⪰ 0 with kernel exactly
span(κ_0, κ_1) = K_j ∩ V_5.

**Theorem N6.** For every k, j with q0 = j − 6 ≥ 4 and t = k − 7 − j ≥ 4, and every σ > 0, Q_6 ⪰ 0 with kernel exactly
span(κ_0, κ_1) = K_j ∩ V_6.

*Proof (computer-assisted).* As for Theorems N1–N4 of `NOTES.md`: q0 = 4 + x, t = 4 + y, i = 5; all sixteen
vectors present, all four blocks split, parts 12 + 4, gauge U_{1,0}, U_{1,1}; entries of at most 1 174 terms in
(x, y, σ); the 9 sanity points (q0, t ∈ {4, 5, 7}) agree with `mlib.py`. The pivots, in the order taken, with the terms
of numerator / denominator (`logs/run_mid_i5.log`, 1 432 s): S_{10} 359/280, S_{00} 976/797, S_{01} 976/797,
S_{11} 976/886; G_1 9/12, G_0 5/15, O_0 182/189, O_1 147/215, W_{10} 8 293/7 965, W_{11} 12 673/13 277,
W_{01} 19 136/18 474, W_{00} 19 099/18 576, U_{00} 38 513/40 772, U_{01} 39 873/40 986 (the same term counts as for
i = 4 in `NOTES.md`). Every pivot has numerator and denominator with positive coefficients only and a nonzero part
free of x, y. Theorem N6 is the same with i = 6 (`logs/run_mid_i6.log`, in the extension `run_extra_mi.sh`, 1 078 s):
entries of at most 1 174 terms, and the same pivot order with exactly the same term counts, all with positive
coefficients only and a nonzero part free of x, y. ∎

## 6. The corollaries

Throughout, as in Corollary P1 of `NOTES.md`: the rescaling of the frame is a positive diagonal congruence
that commutes with S_u, so it keeps the isotypic parts and the inertia on each and maps kernels to kernels; on the
isotypic part of an irrep λ of index i the block is the tensor product of the form on the λ-part of V_i with the
identity of λ; and by Lemma S the λ-part of V_i is spanned by the S vectors when λ = (u − i − 1, i, 1) and by the other
vectors of the new basis when λ = (u − i, i). Since Q_i = Q₂ ⊕ Q₃ and K_j ∩ V_i lies in the two-row part, "Q_i ⪰ 0 with
kernel exactly K_j ∩ V_i" gives the statement on both parts, and for i ≥ 2 (where K_j ∩ V_i = span(κ_0, κ_1)) it says
that Q₃ is positive definite.

**Corollary T.** For every k ≥ 2, every j = −1, …, k and every σ > 0, on the isotypic parts of all S_{k−1}-irreps
(k − 2 − i, i, 1), the letter-1 block of (R∞)_j is ⪰ 0 with kernel exactly the part of K_j there.

*Proof.* Fix (k, j, i) with a three-row part. If k ≤ 40: Theorems MK, M′. If j ≤ 8 or j ≥ k − 8: Theorem MB. If
i ≤ 4: Theorem MS. Otherwise k ≥ 41, 9 ≤ j ≤ k − 9 and i ≥ 5 (and always q0 ≥ −1, t ≥ −2 when V_i ≠ 0). If q0 ≥ 4
and t ≥ 4: Theorem G3. If q0 = a ≤ 3 and t ≥ 4: then i = j − a ≥ 6, and Theorem E applies. If t = c ≤ 3 and q0 ≥ 4:
Theorem E. If both q0 ≤ 3 and t ≤ 3: Theorem F. `coverage_mi.py` checks this bookkeeping for every k ≤ 120. ∎

**Corollary P6.** For every k ≥ 2, every j and every σ > 0, on the isotypic part of the S_{k−1}-irrep (k − 1 − i, i),
the letter-1 block of (R∞)_j is ⪰ 0 with kernel exactly the part of K_j there whenever i ≤ 6, or q0 = j − i ≤ 3, or
t = k − 1 − j − i ≤ 3.

*Proof.* Outside the region k ≥ 41, 9 ≤ j ≤ k − 9, i ≥ 5 as for Corollary T (Theorems MK, M′, MB, MS). Inside it: if
q0 ≥ 4, t ≥ 4 and i ∈ {5, 6}, Theorems N5 and N6; if q0 = a ≤ 3 and t ≥ 4, then i ≥ 6 and Theorem E applies; if
t = c ≤ 3 and q0 ≥ 4, Theorem E; if both are ≤ 3, Theorem F. (Without Theorem N6, the same with i ≤ 5: Corollary P5.) ∎

**Theorem MF′.** For every k ≥ 2, every j and every σ > 0, the letter-1 block of (R∞)_j with squared weights
(σ, 1, …, 1) is ⪰ 0 with kernel exactly K_j on every isotypic part of S_{k−1} except possibly the two-row parts
(k − 1 − i, i) with i ≥ 7, j − i ≥ 4 and k − 1 − j − i ≥ 4. Such parts occur only for k ≥ 23 (k = q0 + t + 2i + 1), and
for k ≤ 55 they are covered by Theorems MK, M′ and M″; so for every k ≤ 55 and every j the whole block is ⪰ 0 with
kernel exactly K_j. The same holds for every weight vector with all squared weights but one equal (permutations and
homogeneity, Corollary M of `relative-anchors/code/marked_copy/NOTES.md`).

*Proof.* Corollaries T and P6 and Theorem M″. The bookkeeping is checked twice: `coverage_mi.py 120 32 6 51` in
`run_all_mi.sh` (without Theorem N6 and the certificates for k ≥ 51: every nonempty form with k ≤ 120 that is not
covered, by Theorem M for k ≤ 32, a certified member or a whole-form family, has k ≥ 51, 9 ≤ j ≤ k − 9, q0 ≥ 4, t ≥ 4
and i ≥ 6, and its three-row part lies in the family of Theorem G3), and `coverage_mi.py 120 32 7 56` in
`run_extra_mi.sh` (with them: the same with k ≥ 56 and i ≥ 7). ∎

## 7. Fixed k = 41, …, 55 (`cert_rest_mi.py 41 50`, `cert_rest_mi.py 51 55`)

For 41 ≤ k ≤ 55 every reduced form whose whole form lies in no proved family (the families of `NOTES.md` and of
Sections 4–5; Theorem G3 alone does not count) is certified by `cert_fixed.certify` (the form over Q(σ), the kernel
identities, a sign-definite gauge, every pivot of the gauge-fixed form positive on (0, ∞) by the sign test of
`signtest.py`). For 41 ≤ k ≤ 50 (`run_all_mi.sh`, without Theorem N6) these are exactly the forms with q0 ≥ 4, t ≥ 4
and i ≥ 6: 1 775 forms (121 at k = 41, rising to 240 at k = 50), all certified, in 430 s
(`logs/run_cert_rest_mi.log`). For 51 ≤ k ≤ 55 (`run_extra_mi.sh`, with Theorem N6) they are the forms with q0 ≥ 4,
t ≥ 4 and i ≥ 7: 1 282 forms (225 at k = 51, rising to 289 at k = 55), all certified, in 229 s
(`logs/run_cert_rest_mi2.log`).

**Theorem M″.** For 41 ≤ k ≤ 55, every j = −1, …, k and every σ > 0, the letter-1 block of (R∞)_j with squared weights
(σ, 1, …, 1) is ⪰ 0 with kernel exactly K_j; and, by Corollary M, the same for every weight vector on these k with all
squared weights but at most one equal.

*Proof.* Every nonempty (k, j, i) with 41 ≤ k ≤ 55 lies in a proved whole-form family or is certified one by one
(`coverage_mi.py`). ∎

## 8. Checks and controls

All exact. `run_all_mi.sh` reruns every step of Sections 3–8 except Theorem N6 and the certificates for 51 ≤ k ≤ 55
(`logs/run_all_mi.out`, one log per step; 64 minutes, ending with ALL PASS); `run_extra_mi.sh` reruns those (below).
- `famcheck_mi.py` (steps `famcheck_gen`, `famcheck_edgeq`): the forms of the generic family (whose three-row part is
  Theorem G3), of [q0 ≥ 4, t ≥ 4, i = 5] and of the eleven edge families, times N0, equal the numeric forms of
  `mlib.py` (exact; the Newton window asserted) entry by entry at every member whose symbolic sizes lie in
  q0 ∈ {6, 17, 29}, t ∈ {5, 19, 31}, i ∈ {7, 13, 19} and with k ≤ 84, at σ ∈ {3/11, 7/2, 40}: 393 (member, σ) points
  and 71 826 entries, with k up to 76. The certificates themselves compare only at the sizes min, min + 1, min + 3.
- `memcheck_mi.py 41 60 7`: at the 1 401 members with 41 ≤ k ≤ 60 and k + j + i ≡ 0 (mod 7) of the generic family
  [q0 ≥ 4, t ≥ 4, i ≥ 5] (828 members) and of the eleven edge families taken with i ≥ 5 (573 members), from the
  numeric orbit calculus: the known kernel vectors are annihilated, the off-diagonal block of the split vanishes, and
  the gauge-fixed form (and its three-row part alone) is positive definite by exact LDL at σ ∈ {10⁻⁶, 1/7, 1, 3, 10⁶}.
  So at these points Q_i ⪰ 0 with kernel exactly K_j ∩ V_i. For the generic family this includes the two-row part,
  which is not proved: there it is evidence only.
- `coverage_mi.py 120 32 6 51`: for every k ≤ 120, every j and every i with V_i ≠ 0: 3 909 forms with k ≤ 32
  (Theorem M), 3 753 certified one by one (1 978 with 33 ≤ k ≤ 40 from `NOTES.md`, 1 775 with 41 ≤ k ≤ 50), 68 120 in
  proved whole-form families; the remaining 84 525 all have k ≥ 51, 9 ≤ j ≤ k − 9, q0 ≥ 4, t ≥ 4 and i ≥ 6, and the
  three-row part of each lies in the family of Theorem G3. Always q0 ≥ −1 and t ≥ −2, and no family reported empty
  contains a nonempty form. (The lookup of the families is indexed; it is compared with the direct test on a sample.)
- Negative controls, B → (1 + 1/1000)B in the symbolic calculus (steps `neg_gen3`, `neg_edges_q`, `neg_edges_t`,
  `neg_i5`; a certificate fails when no remaining diagonal entry passes the sign test, and the smallest of them is
  named): the certificate of Theorem G3 fails in the three-row part (S_{10}); that of [q0 ≥ 4, t ≥ 4, i = 5] at S_{10};
  of the edge families, [q0 = −1, …] and [q0 = 0, …] fail in the two-row part (G_0), [q0 ≥ 4, t = −1, …] there too
  (O_1), and the other seven in the three-row part (S_{10}); the one
  that still passes is [q0 ≥ 4, t = −2, i ≥ 5], whose form is the single vector B_{111} (as for the three one-pivot
  families of Theorem F, the gap there exceeds the change of B). At the sampled members of `memcheck_mi.py`, with B
  scaled the gauge-fixed form is not positive definite at σ = 10⁻⁶ or 10⁶ for all 828 members of the generic family
  (its three-row part alone fails at all of them) and for 512 of the 573 edge members (all but the 61 members of
  [q0 ≥ 4, t = −2, i ≥ 5]); the true forms are positive definite there.

`run_extra_mi.sh` (`logs/run_extra_mi.out`; 22 minutes, ending with ALL PASS) adds: Theorem N6 (step `mid_i6`); the
certificates for 51 ≤ k ≤ 55 (`cert_rest_mi2`); `coverage_mi.py 120 32 7 56` (`coverage_mi2`: 3 909 forms with k ≤ 32,
5 035 certified one by one, 72 705 in proved whole-form families, and the remaining 78 658 all with k ≥ 56,
9 ≤ j ≤ k − 9, q0 ≥ 4, t ≥ 4, i ≥ 7, their three-row parts in the family of Theorem G3); the form of
[q0 ≥ 4, t ≥ 4, i = 6] against `mlib.py` at the far members (`famcheck_i6`: 27 points, 6 912 entries, k up to 73); and
its negative control (`neg_i6`: with B scaled the certificate fails in the three-row part, at S_{10}).

## 9. What is not covered, and where the generic family stopped

- **The two-row parts with i ≥ 7, q0 ≥ 4, t ≥ 4, k ≥ 56.** They need the generic family on the two-row part: 12
  vectors, kernel span(κ_0, κ_1), 10 pivots after the gauge, in the four variables (x, y, z, σ). One value of i at a
  time, the families [q0 ≥ 4, t ≥ 4, i = b] (15 to 25 minutes each for b = 2, …, 6) never finish the job.
- **The attempt** (`logs/dev_mi/probe_gen2.log`: `prove_gen.py ... n,n,n --imin=5 --parts=2`, stopped by a 330-second
  limit). The form took 3.5 s and the sanity check 3 s. Pivots in the order taken, with the terms of numerator /
  denominator and the time when chosen: G_1 15/18 (7.4 s), G_0 9/25 (7.6 s), O_0 494/504 (11.6 s), O_1 423/610
  (34.5 s). After O_1 the remaining entries had 1.2·10⁵ to 3.6·10⁵ terms and each row update took about 80 s. The fifth
  pivot, W_{10} with 73 712/69 508 terms, passed the sign test and was chosen at 256 s; the updates after it were
  stopped at 330 s. For comparison, with i fixed (b = 4 or 5) the same fifth pivot has 8 293/7 965 terms and is reached
  at 11 s, and the last five pivots then take about 24 minutes with up to 4.1·10⁴ terms. So at the fifth pivot the
  generic two-row part is about nine times larger and twenty times slower, and the remaining entries grow faster still;
  at this rate the elimination would take many hours at least, and it was not pursued.
- **The determinant route** (Proposition C of the README) would need det Q₂ (10 × 10 after the gauge) in four variables;
  it is the product of the ten pivots, so computing it is not cheaper than the elimination itself. The structure
  noted in Section 8 of `NOTES.md` (i enters polynomially through orbit sizes and rationally through
  u − c = t + i − δ) bounds the degree in i only over Q(q0, t, t + i, σ), again four variables.
- **Cost data.** The edge families with i symbolic are much cheaper than the generic family: three variables, entries of
  at most 798 terms (against 4 849), largest pivots 23 137 terms, 29 minutes for all eleven. The cost comes from having
  q0 and t symbolic together with i.
- Two or more marked copies, two weight values with both levels large, general weights: not attempted.
- **Stage 2** (the type I blocks with a letter 1 at a marked copy) was not started.

## Files

| File | Content |
|---|---|
| `prove_gen.py` | the family certificate with the split and chosen lower bounds of the symbolic sizes; `--parts=3` for the three-row part only (Section 2) |
| `coverage_mi.py` | coverage bookkeeping for every (k, j, i) with k ≤ KCHK, with whole-form and three-row families and the certified members |
| `cert_rest_mi.py` | fixed-k certificates (`cert_fixed.certify`) for the forms with KMIN ≤ k ≤ KMAX outside the proved whole-form families |
| `famcheck_mi.py` | the new family forms against `mlib.py` at far members (k up to 76) |
| `memcheck_mi.py` | exact eliminations at sampled members with 41 ≤ k ≤ 60, and negative controls there |
| `run_all_mi.sh` | runs the new steps, one Python process at a time under `guard.sh`: `bash run_all_mi.sh > logs/run_all_mi.out`, ending with ALL PASS |
| `run_extra_mi.sh` | the extension, run after `run_all_mi.sh`: Theorem N6, the certificates for 51 ≤ k ≤ 55, the coverage with them, a far-member check and the negative control: `bash run_extra_mi.sh > logs/run_extra_mi.out`, ending with ALL PASS |
| `logs/run_all_mi.out`, `logs/run_extra_mi.out`, `logs/run_*.log` of their steps, `logs/gen3.pkl`, `logs/edges_q.pkl`, `logs/edges_t.pkl`, `logs/mid_i5.pkl`, `logs/mid_i6.pkl`, `logs/members_41_50.pkl`, `logs/members_51_55.pkl` | the final runs and their results |
| `logs/dev_mi/` | development runs, not used by any proof: the first runs of the three-row parts and of the smaller edge families, a run of [q0 ≥ 4, t = 1, i ≥ 5] stopped by hand to free the machine for the final run, the attempt at the generic two-row part (`probe_gen2.log`), `table.awk`, which summarises a certificate log, and its outputs `etable*.txt` used for Section 4 |
| all other files | those of `NOTES.md`, unchanged (`logs/*.pkl` and `logs/frozen/` hold the earlier results read by the coverage) |
