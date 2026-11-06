# The letter-1 block with one marked copy: the parts of low irrep index for every k

5 November 2026. Nothing here is refereed. This folder continues `frozen_e1s/` (read its `NOTES.md` for the setting,
the reduction and Theorems M, B, F, I); `frozen_e1s/` itself is not modified.

## Status

Squared weights s_0 = σ > 0 (copy 0 is marked), s_1 = … = s_{k−1} = 1; u = k − 1; the letter-1 block of (R∞)_j; K_j its
known kernel (Lemma N, Proposition K2, Theorem T of the README). For an irrep index i: q0 = j − i, t = k − 1 − j − i,
m = q0 + t; V_i the (H_i, χ_i)-twisted invariants and Q_i the reduced form on V_i (`frozen_e1s`, Section 2). The parts
with index i are the isotypic parts of the S_u-irreps (u − i, i) and (u − i − 1, i, 1).

| # | Statement | Status |
|---|---|---|
| 1 | **Lemma S (irrep split).** On V_i the basis O, G, U_{x,y0} = A_{x,y0} + B_{x,0,y0}, W_{x,y0} = B_{x,0,y0} + B_{x,1,y0}, S_{x,y0} = i·A_{x,y0} − (n + 1)·B_{x,0,y0} + (m − n)·B_{x,1,y0} (n = q0 − y0 − x; on the blocks with A_{x,y0} ≠ 0) separates the (u − i, i) part (at most 12 vectors) from the (u − i − 1, i, 1) part (the S, at most 4); the known kernel vectors lie in the first part | **proved** by hand (Section 2); in every family certificate the block-diagonality is also **verified identically** in the sizes and σ; **exact check** against the Lemma B block, k ≤ 7 (100 forms, 6 664 entries, σ ∈ {4, 1/9}, B scaled by 1 and 1001/1000), with the class sum of the transpositions acting by the two content sums |
| 2 | **Theorem N1.** For every member of the family [q0 ≥ 4, t ≥ 4, i = 1] and every σ > 0, Q_1 ⪰ 0 with kernel exactly span(κ_0, κ_1) = K_j ∩ V_1 | **proved**, computer-assisted, exact: 4 + 10 pivots, each a ratio of polynomials in (x, y, σ) with positive coefficients after q0 = 4 + x, t = 4 + y (largest about 3·10⁴ terms); about 9 minutes (Sections 3, 4) |
| 3 | **Corollary P1 (the parts i ≤ 1, every k).** For every k ≥ 2, every j = −1, …, k and every σ > 0, on the isotypic parts of the S_{k−1}-irreps (k − 1), (k − 2, 1), (k − 3, 1, 1) the letter-1 block of (R∞)_j is ⪰ 0, with kernel exactly the part of K_j there | **proved**, given Theorems M, F, I of `frozen_e1s` (Theorems F and I rerun here) and row 2 (Section 4); bookkeeping checked for k ≤ 80 |
| 4 | **Theorem M′ (fixed k = 33, …, 40).** For 33 ≤ k ≤ 40, every j and every σ > 0, the letter-1 block of (R∞)_j is ⪰ 0 with kernel exactly K_j (and, as in Corollary M, the same for every weight vector with all squared weights but one equal) | **proved**, computer-assisted, exact: the fixed-k certificates of `frozen_e1s` (`cert_fixed.certify`) for the {{NREST}} reduced forms with 33 ≤ k ≤ 40 not in a proved family, plus the families (Section 5) |
| 5 | {{N2ROW}} | {{N2STATUS}} |
| 6 | The split certificate reproduces Theorem F (the 52 families with one symbolic size: 45 proved, 7 empty) | **exact**, a second elimination (other basis, other gauge, other pivot order) |
| 7 | Negative controls, B → (1 + 1/1000)B: the split certificate fails on 42 of the 45 nonempty one-size families (the same three one-pivot families still pass, as in `frozen_e1s`) and on [q0 ≥ 4, t ≥ 4, i = 1]{{NEG2}} (at the first pivot of the (u − i − 1, i, 1) part); at sampled members with 33 ≤ k ≤ 40 of the new families the gauge-fixed split form with B scaled is not positive definite at σ = 10⁻⁶ or 10⁶ (exact LDL), while the true one is positive definite at σ ∈ {10⁻⁶, 1, 10⁶} | **exact** |
| 8 | k ≥ 41, 9 ≤ j ≤ k − 9 and {{OPENI}}: the generic family [q0 ≥ 4, t ≥ 4, i ≥ 2] (four variables) and the families with i symbolic and one of q0, t symbolic ([q0 = a, t ≥ 4, i ≥ 2], [q0 ≥ 4, t = c, i ≥ 2]) were not proved (Section 8); two or more marked copies; general weights; the type I blocks with a letter 1 at a marked copy (stage 2, not started) | **open** |

Nothing here says that Conjecture 12 of `single-excitation-range/` holds in general, or anything about the NPT (Werner
distillability) problem. Computer-assisted means exact rational arithmetic (python-flint `fmpq`, `fmpq_poly`,
`fmpq_mpoly`) and sign checks of coefficients. No floating point enters a proof.

{{BODY}}
