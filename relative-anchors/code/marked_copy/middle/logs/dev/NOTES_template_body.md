## 1. What was open, and what is done here

`frozen_e1s` proves the statement (the letter-1 block of (R∞)_j with squared weights (σ, 1, …, 1) is ⪰ 0 with kernel
exactly K_j) for k ≤ 32 and every j (Theorem M), and for every k when j ≤ 8 or j ≥ k − 8 (Theorem B), and the part
i = 0 for every k and j (Theorem I). Open there: k ≥ 33, 9 ≤ j ≤ k − 9, parts i ≥ 1.

In that region the parts with q0 ≥ 4 and t ≥ 4 need families with two symbolic sizes ([q0 ≥ 4, t ≥ 4, i = b] for each
b ≥ 1, or the generic [q0 ≥ 4, t ≥ 4, i ≥ 2]); the parts with exactly one of q0 ≤ 3 (then i ≥ 6) and t ≤ 3 (then
i ≥ 5) need families with i symbolic and one more symbolic size ([q0 = a, t ≥ 4, i ≥ 2], [q0 ≥ 4, t = c, i ≥ 2]). The
band families of Theorem B do not reach them; the one-size families of Theorem F only cover the parts with both
q0 ≤ 3 and t ≤ 3 ([q0 = a, t = c, i ≥ 2]).

Done here: the irrep split (Section 2), which with a better gauge makes the two-size families with i fixed tractable
(Section 3); the family [q0 ≥ 4, t ≥ 4, i = 1] and hence every part i ≤ 1 for every k (Section 4){{N2SHORT}}{{RESTSHORT}}.
So stage 1 (every k and every j) is not finished: what remains is listed in Section 8.

## 2. The irrep split (Lemma S)

Notation as in Section 2 of `frozen_e1s/NOTES.md`: q0 = j − i, t = k − 1 − j − i, m = q0 + t, u = k − 1, and the
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
1. the reduced form Q over Q(sizes, σ) by the formal continuation of `famlib.py` (`frozen_e1s`, Section 5);
2. sanity: N0·Q equals the exact numeric form of `mlib.py` at 3^s sample points and two values of σ (mlib asserts its
   Newton window), and every entry denominator has one-signed coefficients after the shifts (so Q is continuous and
   equal to the true reduced form at every member and every σ > 0);
3. the change of basis T and the identity Tᵀ Q T = Q₂ ⊕ Q₃ (the block between the two parts vanishes identically);
4. the known kernel vectors, written in the new basis, lie in the first part and are annihilated identically;
5. gauge: delete one coordinate per kernel vector, with a sign-definite determinant. Here the deleted coordinate of
   κ_{y0} is U_{1,y0} (or A_{1,y0} when that block is not split), where κ_{y0} has coefficient 1 (gauge "R"; deleting
   O_{y0} as in `frozen_e1s` also works but gives larger pivots);
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
equal (by the argument of Corollary M of `frozen_e1s`; the isotypic parts are then those of the permutations of the
copies with equal weight).

*Proof.* The rescaling of `frozen_e1s` (Section 1 there) is a positive diagonal congruence that commutes with S_u, so
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
`cert_fixed.certify` of `frozen_e1s` (Section 3 there: the form over Q(σ), the kernel identities, a sign-definite
gauge, every pivot of the gauge-fixed form positive on (0, ∞) by the sign test of `signtest.py`): 1 978 forms (185 at
k = 33, rising to 313 at k = 40; about 6 minutes). With the families this covers every (k, j, i) with k ≤ 40
(`coverage_mid.py`: 1 978 forms certified one by one, the smallest k of an uncovered form is 41).

**Theorem M′.** For 33 ≤ k ≤ 40, every j = −1, …, k and every σ > 0, the letter-1 block of (R∞)_j with squared weights
(σ, 1, …, 1) is ⪰ 0 with kernel exactly K_j; and, as in Corollary M of `frozen_e1s`, the same for every weight vector
on these k with all squared weights but at most one equal.

{{SEC6}}

## 7. Checks and controls

All exact unless said otherwise; `run_all.sh`, `run_extra.sh`{{B3SCRIPT}} rerun them.
- `check_builder.py 6 6 16` (from `frozen_e1s`): the reduced forms of `mlib.py` against the Lemma B block of `j1lib.py`
  for k ≤ 6 (B scaled by 1 and 1001/1000), end-to-end kernels and positivity of the whole block for k ≤ 6, symbolic
  against numeric forms and window against all Johnson nodes for k ≤ 16.
- `check_split.py 7`: for k ≤ 7, every j and every i with a split block, σ ∈ {4, 1/9}, B scaled by 1 and by
  1001/1000: the vectors U, W, S built over the states from the explicit twisted vectors of `explicit.py` give
  v_aᵀỸv_b (Ỹ the rescaled Lemma B block) equal to the entries of Tᵀ Q T from `mlib.py` (6 664 entries, 100 forms), in
  particular the off-diagonal block vanishes; the class sum of the transpositions of S_u acts on every S as the content
  sum of (u − i − 1, i, 1) and on every vector of the first part as that of (u − i, i); the known kernel vectors lie in
  the first part and are annihilated.
- `prove_fam.py 1 1` and `prove_lowirr.py 0` (from `frozen_e1s`, rerun): Theorems F and I.
- `prove_mid.py` on the 52 one-size families with the split certificate: 45 proved, 7 empty, the same outcome as
  Theorem F with a different basis, gauge and pivot order.
- `coverage_mid.py 80 32 band:2` (in the extensions `band:3`{{B3BAND}}, with the new families): for every k ≤ 80, every j
  and i with V_i ≠ 0: k ≤ 32, or a proved family contains
  (q0, t, i), or the form is certified by `cert_rest.py`; every uncovered form has 9 ≤ j ≤ k − 9, i ≥ 2 and k ≥ 41;
  always q0 ≥ −1, t ≥ −2; no family reported empty contains a nonempty form.
- Negative controls, B → (1 + 1/1000)B: `prove_mid.py … --bscale=1001/1000` fails on 42 of the 45 nonempty one-size
  families (19 of them at the first pivot S_{10} of the (u − i − 1, i, 1) part, the other 23 at a pivot of the first
  part: G_0 in 14, O_0, O_1, A_{11}, B_{100}, B_{101}, B_{110} in the rest; the three families that still pass have a
  single pivot, as in `frozen_e1s`) and on [q0 ≥ 4, t ≥ 4, i = 1] (at S_{10}){{NEG2B}}; `neg_mid.py 40` takes the 30
  members of [q0 ≥ 4, t ≥ 4, i = 1] with 33 ≤ k ≤ 40 and k + j + i ≡ 0 (mod 7){{NEG2C}}: the gauge-fixed split form
  (from `mlib.py`, exact) is positive definite at σ ∈ {10⁻⁶, 1, 10⁶} and, with B scaled, not positive definite at
  σ = 10⁻⁶ or 10⁶.

## 8. What is not covered

- **The middle band for k ≥ 41 and i ≥ {{OPENIMIN}}.** The members there with q0 ≥ 4 and t ≥ 4 need the families
  [q0 ≥ 4, t ≥ 4, i = b] for every b ≥ {{OPENIMIN}} (or the generic [q0 ≥ 4, t ≥ 4, i ≥ 2]); those with exactly one of
  q0 ≤ 3 (then i ≥ 6) and t ≤ 3 (then i ≥ 5) need [q0 = a, t ≥ 4, i ≥ 2] and [q0 ≥ 4, t = c, i ≥ 2]. None of these is
  proved here. (The members with both q0 ≤ 3 and t ≤ 3 are covered by Theorem F.)
- **Cost data** (development runs, `logs/dev/`). [q0 ≥ 4, t ≥ 4, i = 1] with the split and the gauge of `frozen_e1s`
  (deleting O_{y0}) reached 13 of its 14 pivots in 12 minutes (largest 3.3·10⁴ terms) and was stopped by the limit;
  with the gauge deleting U_{1,y0} it finished in 9 minutes. Without the split, `frozen_e1s` reached 11 of 14 pivots in
  12 minutes. {{COSTI2}} At i ≥ 2 all five Johnson nodes of the window enter (at i = 1 only four), which makes the
  entries larger. The generic family has one more variable than these; it was not attempted.
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
| `mlib.py`, `famlib.py`, `prove_fam.py`, `prove_lowirr.py`, `prove_band.py`, `cert_fixed.py`, `explicit.py`, `check_builder.py`, `coverage_fam.py`, `coverage_irr.py`, `j1lib.py`, `flintfield.py`, `ufield.py`, `signtest.py`, `guard.sh` | copied unchanged from `frozen_e1s/` |
| `run_all.sh` | reruns everything above, one Python process at a time under `guard.sh`: `bash run_all.sh > logs/run_all.out`, ending with ALL PASS |
| `run_extra.sh`{{B3SCRIPT}} | {{EXTRAFILE}} |
| `logs/run_all.out`, `logs/run_*.log`, `logs/*.pkl` | the final runs; `logs/frozen/` holds the band and Theorem F/I outcomes and the `run_all.out` of `frozen_e1s` (copied unchanged) |
| `logs/dev/`, `dev/` | development runs and probes, not used by any proof |
