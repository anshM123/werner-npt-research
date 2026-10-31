# Werner states at the NPT endpoint: multi-copy undistillability results

Ansh Mishra and Aryan Senthilkumar

This repository collects our results, as of 31 October 2026, on the question of whether NPT bound entanglement exists.
That is Problem 4 of Horodecki, Rudnicki and Życzkowski, PRX Quantum 3, 010101 (2022). The standard route to it is the
qutrit Werner state ρ ∝ I − F/2: is it undistillable from every number of copies? At the endpoint α = 1/2, n-copy
undistillability is equivalent (Costa Rico; Qi–Gupur–Wu–Guo; Song–Chen) to

    q_n(C) = Σ_J (−1/2)^{|J|} ‖Tr_J C‖_2² ≥ 0   for every operator C of rank at most two on (C^d)^{⊗n}.

**This repository does not solve the problem.** The all-copy statement is still open: we have neither a proof for
every n nor a counterexample. In our exploratory searches (local optimisation; that code is not part of this
repository) the minimum of q_n over rank-two C at α = 1/2 was always 0, and never negative, for n ≤ 11 copies at
d = 3, 7 at d = 4, 6 at d = 5 and 5 at d = 6. The scripted searches here cover the bosonic sector for n ≤ 7 and
local runs at three copies. What follows are partial results, each with its status.

## Contents

| Folder | Result | Status |
|---|---|---|
| [`first-order-bound/`](first-order-bound) | Σ_i ‖Tr_i Z‖² ≤ (n−1)‖Z‖² + ‖Z‖₁² for every operator Z on n factors (sharp). This settles Chattopadhyay–Sarkar's first-order assumption (2006). The paper (13 pp) also has sharp rank-two bounds, a graded form of the "no collective gain" conjecture, and the three-copy case of that conjecture in every dimension given Tomar's theorem (MacCormack and Bao proved it earlier for local dimension at most four). | Proved; checked independently; not refereed |
| [`two-copy-thickness/`](two-copy-thickness) | A margin linear in thickness for two copies: q_2(C) ≥ min_l τ_l ‖C‖²/(48 − 8√2), so m_2(V) ≥ min_l λ_3/37. This answers an open question of Tomar. For operators whose range and co-range lie in one plane, including all normal C, Tomar's inequality holds with the best possible constant 1 − 1/√2 > 1/4. (Tomar proved that constant on a different class of planes.) | Proved; checked independently; not refereed |
| [`copy-amplification/`](copy-amplification) | m-copy undistillability is a budget that can be spent across copies. Combined with Tomar's certificates, this gives the widest known windows valid for every n at fixed d (the only ones that stay above 1/d); at d = 3, 4 they are the best known for every n ≥ 6. | Proved, partly conditional on deposited (unrefereed) computer certificates; not refereed |
| [`bosonic-sector-notes/`](bosonic-sector-notes) | Operators whose range and co-range lie in the symmetric subspace, every n. The endpoint form has an exact SU(3) spectrum. It is nonnegative, indeed at least (1−2t)^n‖C‖², for every n and every rank, on operators whose range lies in one occupation sector z^l·Sym^{n−l}(z^⊥). Numerically, the minimum over rank-two operators is 0 (n ≤ 7), attained at qubit-supported operators. | Sector theorem proved; the full bosonic sector is open |
| [`three-copy-notes/`](three-copy-notes) | Three copies, rank-two projectors. Exact identities: a monogamy form, a triple-cross-product form and an entanglement form of the endpoint. The margin q_3(P_V) ≥ ¾ λ_3(V) is proved on diagonal planes in every dimension (proof corrected after an audit). The sharp three-copy margin is conjectured for all planes, with numerical evidence. A remark shows that crossed Cauchy–Schwarz is not tensor-stable. | Identities and diagonal theorem proved; general margin conjectured |
| [`contraction-notes/`](contraction-notes) | Passing the endpoint from n − 1 to n copies through contractions. An exact identity shows that an "aligned" class of rank-two operators (common backbone on the other copies, arbitrary tails) is nonnegative at n copies whenever n − 1 copies are. Nesting it gives an explicit class, thick at every copy, on which the endpoint inequality holds for every n with no unproven input. The zeros P ⊗ |a⟩⟨b| have nonnegative second variation for every n, given n − 1. A natural site-contraction inequality is false (an explicit three-copy operator). | Identity, nested class and second-variation statement proved; the counterexample is numerical |
| [`compression-notes/`](compression-notes) | A compression conjecture that would give every n: some copy and some direction w with q_n(C) ≥ q_n(C compressed to w^⊥). It comes with an exact formula and adversarial evidence up to five copies. Averages of contractions that do not depend on C can only certify t ≤ 1/d, the PPT range. The projector case is an exact monogamy relation for a purifying qubit, for every n. Translation-invariant MPS filters are described by a completely positive transfer channel; a search up to 33 copies finds no negative value. An explicit witness shows that two-pair Werner mixtures with a positive margin can be activated by one Werner copy. Operators whose range lies in a product of planes over two-copy blocks are nonnegative for every n and every rank. This gives another every-n class, thick at every copy. Given n − 1 copies, a filter on one copy can lower the form in at most one direction, and at a first zero below the endpoint the filter form has kernel exactly the identity. On the code space of a pure qutrit code with n = 2d − 1 the form is exactly A_n‖x‖² + B_n|Tr x|²; for the five-qutrit code this gives q_5 ≥ (71/288)‖C‖² for every operator whose range lies in the code space (co-range arbitrary). So products of planes over blocks of at most three copies and five-qutrit code spaces over five-copy blocks give an every-rank class with 3-dimensional factors. Finally, a reduction through cones of bond operators: if every local filtering of m + 1 Werner copies into dimension 3^m lies in the cone generated by filterings of m copies (KL_m), then the endpoint holds for every n on operators of bond dimension at most 3^m; for 2m copies, output dimension 2·3^(m−1) is enough. So a statement about 6 × 6 filterings of three copies would give the four-copy endpoint. Compressing all copies but one always gives a 6 × 6 filtering, so the endpoint for every n would follow from one such statement per number of copies (each implied by KL_m for m ≥ 2); the analogous one-copy statement at dimension 6 looks false, which matches the dimension bound. The hypotheses are numerically supported for m = 1, and under test for m = 2. Exact positive spectra for the antisymmetric subspace of two copies (smallest eigenvalue 1/4) and for the qutrit [[3,1,2]] code (1/2) give every-n, every-rank classes containing operators thick at every copy, with ranges of dimension up to 3^(n/2). The compressions to the antisymmetric subspace, the [[3,1,2]], [[4,2,2]] and [[5,3,2]] code spaces and the five-qutrit code space are even separable (explicit decompositions, the Gurvits–Barnum ball on syndrome blocks, or exact certificates by stabiliser product states). So a block of two, three, four or five copies where the range lies in such a subspace can be removed. Iterating, the endpoint holds for every n, with no unproven input, whenever the range lies in a product of such blocks and planes with at most three copies left free; this includes cases with up to eight copies. For any qutrit stabiliser code, the partial-transpose spectrum of the compression is an exact weight sum over cosets of the stabiliser group; the codes [[b, b−2, 2]] give PPT blocks exactly for b ≤ 6. A fixed-dimension sufficient condition: if every 3 × 3 local filtering of any number of Werner copies lies in the cone generated by filterings of one copy (Conjecture F3), the endpoint holds for every n; F3 passes adversarial tests at three and four copies. Two Werner copies have Schmidt number 2, and no one-step induction by local congruences of F_2 ⊗ ρ exists (a numerical dual certificate). | Propositions proved; the conjectures (QC, KL_m, F3, one-sided separability) and the large-n searches are numerical |

Each folder has its own README or notes, the LaTeX source and PDF where there is a paper, and verification scripts.
The scripts need Python 3 with NumPy and SciPy. copy-amplification also needs SymPy and Matplotlib, and two of its
scripts need PyTorch.

## What is open

- The endpoint for every n. `compression-notes/` states a compression conjecture that would give it, with evidence up to five copies.
- The endpoint for every n, even restricted to operators whose range and co-range lie in one plane. Section 4 of
  `three-copy-notes/NOTES.md` explains why the most natural inductive route fails.
- The endpoint for every n on the bosonic sector (range and co-range in Sym^n(C^3)). Section 4 of
  `bosonic-sector-notes/NOTES.md` isolates the remaining obstruction: interference between occupation blocks on one
  diagonal.
- Four copies at the exact endpoint. Tomar proves α ≤ 0.4979 for qutrits. `contraction-notes/` covers one class
  of operators that are thick at a copy. Section 9 of `compression-notes/NOTES.md` reduces it to a three-copy
  statement: every 6 × 6 local filtering of three Werner copies lies in the cone generated by filterings of two.
  Section 10 of the same notes proves the four-copy endpoint for every operator whose range at some pair of copies
  lies in a local-unitary image of the antisymmetric subspace, or at some triple in the [[3,1,2]] code space.
- Tomar's two-copy constant 1/4 for operators whose range and co-range differ. We prove 1/37 in general, and the
  sharp constant for the plane class.

## How to cite

Please cite the individual papers: each folder has a `CITATION.cff`. For the collection as a whole, use the
`CITATION.cff` at the top level.

## License

MIT, see `LICENSE`.
