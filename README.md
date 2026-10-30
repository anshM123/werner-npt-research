# Werner states at the NPT endpoint: multi-copy undistillability results

Ansh Mishra and Aryan Senthilkumar

This repository collects our results, as of 30 October 2026, on the question of whether NPT bound entanglement exists.
Concretely: is the qutrit Werner state ρ ∝ I − F/2 undistillable from every number of copies? (This is Problem 4 of
Horodecki, Rudnicki and Życzkowski, PRX Quantum 3, 010101 (2022).) At the endpoint α = 1/2, n-copy undistillability
is equivalent to

    q_n(C) = Σ_J (−1/2)^{|J|} ‖Tr_J C‖_2² ≥ 0   for every operator C of rank at most two on (C^d)^{⊗n}.

**This repository does not solve the problem.** The all-copy statement is still open: we have neither a proof for
every n nor a counterexample. Every numerical search we ran ends exactly at the product value, for n ≤ 11 copies at
d = 3, 7 at d = 4, 6 at d = 5 and 5 at d = 6. What follows are partial results, each with its status.

## Contents

| Folder | Result | Status |
|---|---|---|
| [`first-order-bound/`](first-order-bound) | Σ_i ‖Tr_i Z‖² ≤ (n−1)‖Z‖² + ‖Z‖₁² for every operator Z on n factors (sharp). This settles Chattopadhyay–Sarkar's first-order assumption (2006). The paper (13 pp) also has sharp rank-two bounds, a graded form of the "no collective gain" conjecture, and the three-copy case of that conjecture given Tomar's theorem. | Proved; checked independently; not refereed |
| [`two-copy-thickness/`](two-copy-thickness) | A margin linear in thickness for two copies: q_2(C) ≥ min_l τ_l ‖C‖²/(48 − 8√2), so m_2(V) ≥ min_l λ_3/37. This answers an open question of Tomar. For operators whose range and co-range lie in one plane, including all normal C, Tomar's inequality holds with the best possible constant 1 − 1/√2 > 1/4. (Tomar proved that constant on a different class of planes.) | Proved; checked independently; not refereed |
| [`copy-amplification/`](copy-amplification) | m-copy undistillability is a budget that can be spent across copies. Combined with Tomar's certificates, this gives the best known n-copy undistillability windows at fixed d. | Proved, partly conditional on published computer certificates; not refereed |
| [`three-copy-notes/`](three-copy-notes) | Three copies, rank-two projectors. Exact identities: a monogamy form, a triple-cross-product form and an entanglement form of the endpoint. The margin q_3(P_V) ≥ ¾ λ_3(V) is proved on diagonal planes in every dimension. The sharp three-copy margin is conjectured for all planes, with numerical evidence. A remark shows that crossed Cauchy–Schwarz is not tensor-stable. | Identities and diagonal theorem proved; general margin conjectured |

Each folder has its own README or notes, the LaTeX source and PDF where there is a paper, and verification scripts.
The scripts need Python 3 with NumPy; some also need SciPy.

## What is open

- The endpoint for every n, even restricted to operators whose range and co-range lie in one plane. Section 4 of
  `three-copy-notes/NOTES.md` explains why the most natural inductive route fails.
- Four copies at the exact endpoint. Tomar proves α ≤ 0.4979 for qutrits.
- Tomar's two-copy constant 1/4 for operators whose range and co-range differ. We prove 1/37 in general, and the
  sharp constant for the plane class.

## How to cite

Please cite the individual papers: each folder has a `CITATION.cff`. For the collection as a whole, use the
`CITATION.cff` at the top level.

## License

MIT, see `LICENSE`.
