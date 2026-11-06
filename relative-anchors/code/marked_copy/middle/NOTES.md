# The letter-1 block with one marked copy: the parts of low irrep index for every k

5 November 2026. Nothing here is refereed. In `relative-anchors/README.md` (Section 4b) Theorems N1 to N4 with Corollary P4 are summarized as Theorem MS, and Theorem M′ extends Theorem MK to k ≤ 40. Before running `run_all.sh` here, run the parent `run_all.sh` and then `prepare_frozen.sh`. The parts with irrep index i ≥ 5 are continued in `NOTES_MI.md` (run `run_all_mi.sh` and `run_extra_mi.sh` after the extensions here). This folder continues `marked_copy/` (the parent folder) (read its `NOTES.md` for the setting,
the reduction and Theorems M, B, F, I); `marked_copy/` (the parent folder) itself is not modified.

## Status

Squared weights s_0 = σ > 0 (copy 0 is marked), s_1 = … = s_{k−1} = 1; u = k − 1; the letter-1 block of (R∞)_j; K_j its
known kernel (Lemma N, Proposition K2, Theorem T of the README). For an irrep index i: q0 = j − i, t = k − 1 − j − i,
m = q0 + t; V_i the (H_i, χ_i)-twisted invariants and Q_i the reduced form on V_i (`marked_copy`, Section 2). The parts
with index i are the isotypic parts of the S_u-irreps (u − i, i) and (u − i − 1, i, 1).

| # | Statement | Status |
|---|---|---|
| 1 | **Lemma S (irrep split).** On V_i the basis O, G, U_{x,y0} = A_{x,y0} + B_{x,0,y0}, W_{x,y0} = B_{x,0,y0} + B_{x,1,y0}, S_{x,y0} = i·A_{x,y0} − (n + 1)·B_{x,0,y0} + (m − n)·B_{x,1,y0} (n = q0 − y0 − x; on the blocks with A_{x,y0} ≠ 0) separates the (u − i, i) part (at most 12 vectors) from the (u − i − 1, i, 1) part (the S, at most 4); the known kernel vectors lie in the first part | **proved** by hand (Section 2); in every family certificate the block-diagonality is also **verified identically** in the sizes and σ; **exact check** against the Lemma B block, k ≤ 7 (100 forms, 6 664 entries, σ ∈ {4, 1/9}, B scaled by 1 and 1001/1000), with the class sum of the transpositions acting by the two content sums |
| 2 | **Theorem N1.** For every member of the family [q0 ≥ 4, t ≥ 4, i = 1] and every σ > 0, Q_1 ⪰ 0 with kernel exactly span(κ_0, κ_1) = K_j ∩ V_1 | **proved**, computer-assisted, exact: 4 + 10 pivots, each a ratio of polynomials in (x, y, σ) with positive coefficients after q0 = 4 + x, t = 4 + y (largest about 3·10⁴ terms); about 9 minutes (Sections 3, 4) |
| 3 | **Corollary P1 (the parts i ≤ 1, every k).** For every k ≥ 2, every j = −1, …, k and every σ > 0, on the isotypic parts of the S_{k−1}-irreps (k − 1), (k − 2, 1), (k − 3, 1, 1) the letter-1 block of (R∞)_j is ⪰ 0, with kernel exactly the part of K_j there | **proved**, given Theorems M, F, I of `marked_copy` (Theorems F and I rerun here) and row 2 (Section 4); bookkeeping checked for k ≤ 80 |
| 4 | **Theorem M′ (fixed k = 33, …, 40).** For 33 ≤ k ≤ 40, every j and every σ > 0, the letter-1 block of (R∞)_j is ⪰ 0 with kernel exactly K_j (and, as in Corollary M, the same for every weight vector with all squared weights but one equal) | **proved**, computer-assisted, exact: the fixed-k certificates of `marked_copy` (`cert_fixed.certify`) for the 1 978 reduced forms with 33 ≤ k ≤ 40 not in a proved family, plus the families (Section 5) |
| 5 | **Theorems N2, N3, N4 and Corollary P4 (the parts i ≤ 4, every k).** For b = 2, 3, 4, every member of [q0 ≥ 4, t ≥ 4, i = b] and every σ > 0: Q_b ⪰ 0 with kernel exactly span(κ_0, κ_1). Hence for every k ≥ 2, every j = −1, …, k and every σ > 0, on the isotypic parts of all S_{k−1}-irreps (k − 1 − i, i) and (k − 2 − i, i, 1) with i ≤ 4, the letter-1 block of (R∞)_j is ⪰ 0 with kernel exactly the part of K_j there | **proved**, computer-assisted, exact (`run_extra.sh`, `run_extra3.sh`, `run_extra4.sh`, each ending ALL PASS): 4 + 10 pivots per family with positive coefficients after the shifts (largest about 4.1·10⁴ terms; 17, 15 and 20 minutes), with the band families of Theorem B and Theorem M for the edges (Section 6); bookkeeping checked for k ≤ 80 |
| 6 | The split certificate reproduces Theorem F (the 52 families with one symbolic size: 45 proved, 7 empty) | **exact**, a second elimination (other basis, other gauge, other pivot order) |
| 7 | Negative controls, B → (1 + 1/1000)B: the split certificate fails on 42 of the 45 nonempty one-size families (the same three one-pivot families still pass, as in `marked_copy`) and on [q0 ≥ 4, t ≥ 4, i = 1], [q0 ≥ 4, t ≥ 4, i = 2], [q0 ≥ 4, t ≥ 4, i = 3] and [q0 ≥ 4, t ≥ 4, i = 4] (at the first pivot of the (u − i − 1, i, 1) part); at sampled members with 33 ≤ k ≤ 40 of the new families the gauge-fixed split form with B scaled is not positive definite at σ = 10⁻⁶ or 10⁶ (exact LDL), while the true one is positive definite at σ ∈ {10⁻⁶, 1, 10⁶} | **exact** |
| 8 | k ≥ 41, 9 ≤ j ≤ k − 9 and i ≥ 5, except the members with both q0 ≤ 3 and t ≤ 3 (Theorem F): the generic family [q0 ≥ 4, t ≥ 4, i ≥ 2] (four variables) and the families with i symbolic and one of q0, t symbolic ([q0 = a, t ≥ 4, i ≥ 2], [q0 ≥ 4, t = c, i ≥ 2]) were not proved (Section 8); two or more marked copies; general weights; the type I blocks with a letter 1 at a marked copy (stage 2, not started) | **open** |

Nothing here says that Conjecture 12 of `single-excitation-range/` holds in general, or anything about the NPT (Werner
distillability) problem. Computer-assisted means exact rational arithmetic (python-flint `fmpq`, `fmpq_poly`,
`fmpq_mpoly`) and sign checks of coefficients. No floating point enters a proof.

## 1. What was open, and what is done here

`marked_copy` proves the statement (the letter-1 block of (R∞)_j with squared weights (σ, 1, …, 1) is ⪰ 0 with kernel
exactly K_j) for k ≤ 32 and every j (Theorem M), and for every k when j ≤ 8 or j ≥ k − 8 (Theorem B), and the part
i = 0 for every k and j (Theorem I). Open there: k ≥ 33, 9 ≤ j ≤ k − 9, parts i ≥ 1.

In that region the parts with q0 ≥ 4 and t ≥ 4 need families with two symbolic sizes ([q0 ≥ 4, t ≥ 4, i = b] for each
b ≥ 1, or the generic [q0 ≥ 4, t ≥ 4, i ≥ 2]); the parts with exactly one of q0 ≤ 3 (then i ≥ 6) and t ≤ 3 (then
i ≥ 5) need families with i symbolic and one more symbolic size ([q0 = a, t ≥ 4, i ≥ 2], [q0 ≥ 4, t = c, i ≥ 2]). The
band families of Theorem B do not reach them; the one-size families of Theorem F only cover the parts with both
q0 ≤ 3 and t ≤ 3 ([q0 = a, t = c, i ≥ 2]).

Done here: the irrep split (Section 2), which with a better gauge makes the two-size families with i fixed tractable
(Section 3); the family [q0 ≥ 4, t ≥ 4, i = 1] and hence every part i ≤ 1 for every k (Section 4); the families
[q0 ≥ 4, t ≥ 4, i = b] for b = 2, 3, 4, and hence every part i ≤ 4 for every k (Section 6); and, by fixed-k
certificates, every j and every i for 33 ≤ k ≤ 40 (Section 5). So stage 1 (every k and every j) is not finished: what
remains is listed in Section 8.

## 2. The irrep split (Lemma S)

Notation as in Section 2 of `marked_copy/NOTES.md`: q0 = j − i, t = k − 1 − j − i, m = q0 + t, u = k − 1, and the
sixteen twisted vectors O_{y0}, A_{x,y0}, B_{x,y,y0}, G_{z0} spanning V_i (the (H_i, χ_i)-twisted invariants).

**Lemma S.** Fix x, y0 ∈ {0, 1} with A_{x,y0} ≠ 0 (that is, 0 ≤ n ≤ m − 1 for n = q0 − y0 − x), and i ≥ 1. Put

    S_{x,y0} = i·A_{x,y0} − (n + 1)·B_{x,0,y0} + (m − n)·B_{x,1,y0},
    U_{x,y0} = A_{x,y0} + B_{x,0,y0},       W_{x,y0} = B_{x,0,y0} + B_{x,1,y0}.

Then S_{x,y0} lies in the isotypic part of (u − i − 1, i, 1), and U_{x,y0}, W_{x,y0} lie in that of (u − i, i).

*Proof.* The states with the row letter in U, letter x at the row site and y0 at the marked copy form the permutation
module of pairs (r, S), r ∈ U, S ⊂ U − r, |S| = n + i. The two equivariant maps (r, S) ↦ S and (r, S) ↦ S ∪ {r} go to
permutation modules of subsets of U, which contain only two-row irreps; so they vanish on the (u − i − 1, i, 1) part.
On the twisted vectors they act as follows (count the states): A_{x,y0} ↦ (m − n)·D_n and (n + 1)·D_{n+1},
B_{x,0,y0} ↦ 0 and i·D_{n+1}, B_{x,1,y0} ↦ −i·D_n and 0, where D_c is the product of the singlets with the Dicke vector
of weight c on the unpaired sites (D_n, D_{n+1} ≠ 0 since 0 ≤ n ≤ m − 1). The common kernel on span(A, B0, B1) is
spanned by S_{x,y0}; it contains the twisted vector of the (u − i − 1, i, 1) part, which is nonzero (A, B0, B1 are
linearly independent with disjoint supports and the (u − i, i) part has multiplicity at most 2 there, by the Kostka
number), so the two coincide. The (u − i, i) part is the orthogonal complement for the standard inner product, in which
|A|² : |B0|² : |B1|² = (m − n)/i : (m − n)/(n + 1) : 1 (from the norms of `mlib.norm_of`); a vector aA + bB0 + cB1 is
orthogonal to S_{x,y0} iff (m − n)(a − b + c) = 0, which U and W satisfy. ∎

The proofs below do not use Lemma S: on each family the certificate checks directly that the change of basis T
(O, G and the B vectors of blocks with A_{x,y0} = 0 unchanged; U, W, S on the other blocks, det = m + 1 + i > 0 per
block) gives Tᵀ Q T = Q₂ ⊕ Q₃ identically in the sizes and σ. Lemma S explains why, and `check_split.py` confirms it
against the Lemma B block and the class sum of the transpositions.

The known kernel vectors lie in the first part: κ_{y0} = O_{y0} + A_{0,y0} + A_{1,y0} + B_{1,0,y0} − B_{0,1,y0}
= O_{y0} + U_{0,y0} − W_{0,y0} + U_{1,y0} on split blocks.

## 3. The family certificate with the split (`prove_split.py`)

For a family (each of q0, t, i fixed or symbolic, symbolic meaning q0 = 4 + x, t = 4 + y, i = 2 + z with x, y, z ≥ 0;
σ > 0 always symbolic):
1. the reduced form Q over Q(sizes, σ) by the formal continuation of `famlib.py` (`marked_copy`, Section 5);
2. sanity: N0·Q equals the exact numeric form of `mlib.py` at 3^s sample points and two values of σ (mlib asserts its
   Newton window), and every entry denominator has one-signed coefficients after the shifts (so Q is continuous and
   equal to the true reduced form at every member and every σ > 0);
3. the change of basis T and the identity Tᵀ Q T = Q₂ ⊕ Q₃ (the block between the two parts vanishes identically);
4. the known kernel vectors, written in the new basis, lie in the first part and are annihilated identically;
5. gauge: delete one coordinate per kernel vector, with a sign-definite determinant. Here the deleted coordinate of
   κ_{y0} is U_{1,y0} (or A_{1,y0} when that block is not split), where κ_{y0} has coefficient 1 (gauge "R"; deleting
   O_{y0} as in `marked_copy` also works but gives larger pivots);
6. greedy symmetric elimination of Q₃ and of the gauge-fixed Q₂: every pivot is a ratio of two polynomials in
   (x, y, z, σ) whose coefficients all have one common sign and whose part free of the size variables is nonzero.

*Why this is a proof.* At every member (an integer point of the family) and every σ > 0, the entries of Q are the
entries of the true reduced form (step 2), Tᵀ Q T is a congruence by an invertible matrix (step 3), so Q ⪰ 0 with
kernel K iff Q₂ ⪰ 0 with kernel K and Q₃ ≻ 0 (step 4: K lies in the first part). The gauge coordinates give a
complement of K in the first part (step 5). A symmetric elimination with positive pivots in any order shows that the
corresponding symmetric permutation has all leading principal minors positive (each is the product of the pivots so
far, as an identity of rational functions that specialises at every point where the entry denominators and the
earlier pivots are positive), so it is positive definite (Sylvester).

## 4. The family [q0 ≥ 4, t ≥ 4, i = 1] and the parts i ≤ 1 (`prove_mid.py n,n,1`)

**Theorem N1.** For every k, j with q0 = j − 1 ≥ 4 and t = k − 2 − j ≥ 4, and every σ > 0, the reduced form Q_1 is
⪰ 0 with kernel exactly span(κ_0, κ_1) = K_j ∩ V_1.

*Proof (computer-assisted).* The certificate of Section 3 on the family with q0 = 4 + x, t = 4 + y symbolic and i = 1.
All sixteen twisted vectors are present and all four blocks split, so the parts have 12 and 4 vectors. Sanity at the 9
sample points (q0, t ∈ {4, 5, 7}) and σ ∈ {2/7, 5}; every entry denominator one-signed; the off-diagonal 12 × 4 block
vanishes identically; κ_0, κ_1 are annihilated identically; the gauge deletes U_{1,0}, U_{1,1} (determinant 1). The
pivots, in the order taken: part (u − 2, 1, 1): S_{10}, S_{00}, S_{01}, S_{11} (at most 975 terms); part (u − 1, 1):
G_1, G_0, O_0, O_1, W_{01}, W_{10}, W_{11}, W_{00}, U_{00}, U_{01} (the last two with about 2.8·10⁴ and 3.0·10⁴ terms in
numerator and denominator). Every pivot has numerator and denominator with positive coefficients only, and a nonzero
part free of x, y. ∎

**Corollary P1.** For every k ≥ 2, every j = −1, …, k and every σ > 0, on the isotypic parts of the S_{k−1}-irreps
(k − 1), (k − 2, 1) and (k − 3, 1, 1), the letter-1 block of (R∞)_j with squared weights (σ, 1, …, 1) is ⪰ 0, and its
kernel there is exactly the part of K_j there. The same holds for every weight vector with all squared weights but one
equal (by the argument of Corollary M of `marked_copy`; the isotypic parts are then those of the permutations of the
copies with equal weight).

*Proof.* The rescaling of `marked_copy` (Section 1 there) is a positive diagonal congruence that commutes with S_u, so
it keeps the isotypic parts, the inertia on each, and maps kernels to kernels. On the isotypic part of an irrep λ of
index i the block is the tensor product of the form on the λ-part of V_i with the identity of λ; the parts i ≤ 1 are
exactly the irreps listed ((u − 1, 0, 1) is not a partition). So it suffices that Q_0 and Q_1 are ⪰ 0 with kernels
exactly K_j ∩ V_0 and K_j ∩ V_1 (K_j is S_u-invariant). For i = 0 this is Theorem I. For i = 1: if q0 ≥ 4 and t ≥ 4,
Theorem N1; if exactly one of q0, t is ≥ 4, the member lies in a family of Theorem F ([q0 = a, t ≥ 4, i = 1] with
a ∈ {−1, …, 3}, or [q0 ≥ 4, t = c, i = 1] with c ∈ {−2, …, 3}); if both are ≤ 3 then k = q0 + t + 3 ≤ 9 and Theorem M
applies (q0 ≥ −1 and t ≥ −2 whenever V_1 ≠ 0). `coverage_mid.py` checks this bookkeeping (with all other proved
families) for every k ≤ 80. ∎

## 5. Fixed k = 33, …, 40 (`cert_rest.py`)

For 33 ≤ k ≤ 40 every reduced form not contained in a proved family (Theorems F, I, B, N1) is certified by
`cert_fixed.certify` of `marked_copy` (Section 3 there: the form over Q(σ), the kernel identities, a sign-definite
gauge, every pivot of the gauge-fixed form positive on (0, ∞) by the sign test of `signtest.py`): 1 978 forms (185 at
k = 33, rising to 313 at k = 40; about 6 minutes). With the families this covers every (k, j, i) with k ≤ 40
(`coverage_mid.py`: 1 978 forms certified one by one, the smallest k of an uncovered form is 41).

**Theorem M′.** For 33 ≤ k ≤ 40, every j = −1, …, k and every σ > 0, the letter-1 block of (R∞)_j with squared weights
(σ, 1, …, 1) is ⪰ 0 with kernel exactly K_j; and, as in Corollary M of `marked_copy`, the same for every weight vector
on these k with all squared weights but at most one equal.

## 6. The families [q0 ≥ 4, t ≥ 4, i = b] for b = 2, 3, 4 (`run_extra.sh`, `run_extra3.sh`, `run_extra4.sh`)

**Theorem N_b.** Let b ∈ {2, 3, 4}. For every k, j with q0 = j − b ≥ 4 and t = k − 1 − j − b ≥ 4, and every σ > 0, the
reduced form Q_b is ⪰ 0 with kernel exactly span(κ_0, κ_1) = K_j ∩ V_b.

*Proof (computer-assisted).* The certificate of Section 3 on the family with q0 = 4 + x, t = 4 + y symbolic and i = b,
exactly as for Theorem N1: all sixteen vectors present, all four blocks split, parts 12 + 4, gauge U_{1,0}, U_{1,1};
now all five Johnson nodes of the window enter (b + e ≥ 0 for e ≥ −2). The pivots, in the order taken, with the terms
of numerator / denominator:
- b = 2 (`logs/run_mid_i2.log`, 1 032 s): S_{10} 359/280, S_{00} 976/797, S_{01} 976/797, S_{11} 976/886; G_1 8/12,
  G_0 4/14, O_0 175/182, O_1 141/208, W_{11} 7 449/7 563, W_{10} 12 129/12 185, W_{01} 18 416/17 788,
  W_{00} 19 070/18 546, U_{00} 38 513/40 737, U_{01} 39 873/40 986.
- b = 3 (`logs/run_mid_i3.log`, 916 s): S_{10} 359/280, S_{00} 976/797, S_{01} 976/797, S_{11} 976/886; G_1 9/12,
  G_0 5/15, O_0 182/189, O_1 147/215, W_{11} 7 874/7 965, W_{10} 12 673/12 736, W_{01} 19 136/18 474,
  W_{00} 19 099/18 576, U_{00} 38 512/40 772, U_{01} 39 873/40 986.
- b = 4 (`logs/run_mid_i4.log`, 1 210 s): S_{10} 359/280, S_{00} 976/797, S_{01} 976/797, S_{11} 976/886; G_1 9/12,
  G_0 5/15, O_0 182/189, O_1 147/215, W_{10} 8 293/7 965, W_{11} 12 673/13 277, W_{01} 19 136/18 474,
  W_{00} 19 099/18 576, U_{00} 38 513/40 772, U_{01} 39 873/40 986.

Every pivot has numerator and denominator with positive coefficients only, and a nonzero part free of x, y. ∎

**Corollary P4.** For every k ≥ 2, every j = −1, …, k and every σ > 0, on the isotypic parts of all S_{k−1}-irreps
(k − 1 − i, i) and (k − 2 − i, i, 1) with i ≤ 4, the letter-1 block of (R∞)_j with squared weights (σ, 1, …, 1) is
⪰ 0, with kernel exactly the part of K_j there (and the same for every weight vector with all squared weights but one
equal).

*Proof.* As for Corollary P1, it suffices to treat Q_i for 2 ≤ i ≤ 4. If q0 ≥ 4 and t ≥ 4, Theorem N_i. If
q0 = a ≤ 3 and t ≥ 4, the member lies in the band family [q0 = a, t ≥ 4, i] of Theorem B (proved there for a + i ≤ 8,
and here a + i ≤ 7); if q0 ≥ 4 and t = c ≤ 3, in the band family [q0 ≥ 4, t = c, i] (proved for c + i ≤ 7, and
here c + i ≤ 7); if both are ≤ 3, k = q0 + t + 2i + 1 ≤ 15 and Theorem M applies.
`coverage_mid.py 80 32 band:5` checks the bookkeeping for k ≤ 80: every form not covered has 9 ≤ j ≤ k − 9,
i ≥ 5 and k ≥ 41. ∎

The band families stop covering the edges at i = 5 ([q0 ≥ 4, t = 3, i = 5] has c + i = 8), so beyond i = 4 the edge
families with i fixed would have to be added as well.


## 7. Checks and controls

All exact unless said otherwise; `run_all.sh`, `run_extra.sh`, `run_extra3.sh`, `run_extra4.sh` rerun them.
- `check_builder.py 6 6 16` (from `marked_copy`): the reduced forms of `mlib.py` against the Lemma B block of `j1lib.py`
  for k ≤ 6 (B scaled by 1 and 1001/1000), end-to-end kernels and positivity of the whole block for k ≤ 6, symbolic
  against numeric forms and window against all Johnson nodes for k ≤ 16.
- `check_split.py 7`: for k ≤ 7, every j and every i with a split block, σ ∈ {4, 1/9}, B scaled by 1 and by
  1001/1000: the vectors U, W, S built over the states from the explicit twisted vectors of `explicit.py` give
  v_aᵀỸv_b (Ỹ the rescaled Lemma B block) equal to the entries of Tᵀ Q T from `mlib.py` (6 664 entries, 100 forms), in
  particular the off-diagonal block vanishes; the class sum of the transpositions of S_u acts on every S as the content
  sum of (u − i − 1, i, 1) and on every vector of the first part as that of (u − i, i); the known kernel vectors lie in
  the first part and are annihilated.
- `prove_fam.py 1 1` and `prove_lowirr.py 0` (from `marked_copy`, rerun): Theorems F and I.
- `prove_mid.py` on the 52 one-size families with the split certificate: 45 proved, 7 empty, the same outcome as
  Theorem F with a different basis, gauge and pivot order.
- `coverage_mid.py 80 32 band:2` (in the extensions `band:3`, `band:4`, `band:5`, with the new families): for every
  k ≤ 80, every j and i with V_i ≠ 0: k ≤ 32, or a proved family contains (q0, t, i), or the form is certified by
  `cert_rest.py`; every uncovered form has 9 ≤ j ≤ k − 9, k ≥ 41 and i ≥ 2 (in the extensions i ≥ 3, 4, 5); always
  q0 ≥ −1, t ≥ −2; no family reported empty contains a nonempty form. After `run_extra4.sh`, 29 230 forms with k ≤ 80
  remain uncovered, all with 9 ≤ j ≤ k − 9, i ≥ 5, k ≥ 41.
- Negative controls, B → (1 + 1/1000)B: `prove_mid.py … --bscale=1001/1000` fails on 42 of the 45 nonempty one-size
  families (19 of them at the first pivot S_{10} of the (u − i − 1, i, 1) part, the other 23 at a pivot of the first
  part: G_0 in 14, O_0, O_1, A_{11}, B_{100}, B_{101}, B_{110} in the rest; the three families that still pass have a
  single pivot, as in `marked_copy`) and on [q0 ≥ 4, t ≥ 4, i = b] for b = 1, 2, 3, 4 (each at S_{10});
  `neg_mid.py 40` takes the members of [q0 ≥ 4, t ≥ 4, i = b] with 33 ≤ k ≤ 40 and k + j + i ≡ 0 (mod 7) (30, 28, 25
  and 23 members for b = 1, 2, 3, 4): the gauge-fixed split form (from `mlib.py`, exact) is positive definite at
  σ ∈ {10⁻⁶, 1, 10⁶} and, with B scaled, not positive definite at σ = 10⁻⁶ or 10⁶.

## 8. What is not covered

- **The middle band for k ≥ 41 and i ≥ 5.** The members there with q0 ≥ 4 and t ≥ 4 need the families
  [q0 ≥ 4, t ≥ 4, i = b] for every b ≥ 5 (or the generic [q0 ≥ 4, t ≥ 4, i ≥ 2]); those with exactly one of
  q0 ≤ 3 (then i ≥ 6) and t ≤ 3 (then i ≥ 5) need [q0 = a, t ≥ 4, i ≥ 2] and [q0 ≥ 4, t = c, i ≥ 2]. None of these is
  proved here. (The members with both q0 ≤ 3 and t ≤ 3 are covered by Theorem F.)
- **Cost data** (development runs, `logs/dev/`). [q0 ≥ 4, t ≥ 4, i = 1] with the split and the gauge of `marked_copy`
  (deleting O_{y0}) reached 13 of its 14 pivots in 12 minutes (largest 3.3·10⁴ terms) and was stopped by the limit;
  with the gauge deleting U_{1,y0} it finished in 9 minutes. Without the split, `marked_copy` reached 11 of 14 pivots in
  12 minutes. With the same gauge, the families with i = 2, 3, 4 took 17, 15 and 20 minutes (largest pivots about
  4·10⁴ terms). At i ≥ 2 all five Johnson nodes of the window enter (at i = 1 only four), which makes the entries
  larger. Each further b costs about a quarter of an hour, but beyond b = 4 the edge families with i = b fixed are
  needed as well. The generic family has one more variable than these; it was not attempted.
- **A determinant certificate** in the style of Proposition C would need a structure lemma bounding the degree in one
  size. In the calculus, i enters polynomially (degree ≤ 2) through the orbit sizes 2i, 2im, 2i(2i − 2), but also
  rationally through u − c = t + i − δ in the energies of Lemma E, so the argument of Proposition C (where τ entered only
  polynomially) does not apply as it stands.
- Two or more marked copies, two weight values with both levels large, general weights: not attempted.
- **Stage 2** (the type I blocks with a letter 1 at a marked copy) was not started, since stage 1 is not finished.

## Files

| File | Content |
|---|---|
| `splitlib.py` | Lemma S: the new basis U, W, S, the congruence Tᵀ Q T, the kernel vectors in the new basis |
| `prove_split.py` | the family certificate with the split (Section 3); `GAUGE=R` selects the gauge deleting U_{1,y0} |
| `prove_mid.py` | runs `prove_split` on a list of families, saves the outcomes (`logs/*.pkl`) |
| `cert_rest.py` | fixed-k certificates for the forms with KMIN ≤ k ≤ KMAX outside the proved families |
| `coverage_mid.py` | coverage bookkeeping for every (k, j, i) with k ≤ KCHK |
| `check_split.py` | the split against the Lemma B block, class sums, kernel vectors (k ≤ 7) |
| `neg_mid.py` | negative controls at explicit members of the new families |
| `mlib.py`, `famlib.py`, `prove_fam.py`, `prove_lowirr.py`, `prove_band.py`, `cert_fixed.py`, `explicit.py`, `check_builder.py`, `coverage_fam.py`, `coverage_irr.py`, `j1lib.py`, `flintfield.py`, `ufield.py`, `signtest.py`, `guard.sh` | copied unchanged from `marked_copy/` (the parent folder) |
| `run_all.sh` | reruns everything above, one Python process at a time under `guard.sh`: `bash run_all.sh > logs/run_all.out`, ending with ALL PASS |
| `run_extra.sh`, `run_extra3.sh`, `run_extra4.sh` | the extensions, run after `run_all.sh` and in this order: the families [q0 ≥ 4, t ≥ 4, i = b], b = 2, 3, 4, with their coverage (`band:3`, `band:4`, `band:5`) and negative controls: `bash run_extra.sh > logs/run_extra.out`, `bash run_extra3.sh > logs/run_extra3.out`, `bash run_extra4.sh > logs/run_extra4.out`, each ending with ALL PASS |
| `logs/run_all.out`, `logs/run_extra*.out`, `logs/run_*.log`, `logs/*.pkl` | the final runs (one log per step); `logs/frozen/` holds the band and Theorem F/I outcomes and the `run_all.out` of `marked_copy` (copied unchanged) |
| `logs/dev/` | development runs and probes, the templates of these notes and the helper that filled them in; not used by any proof |

