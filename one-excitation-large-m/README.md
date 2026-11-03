# One free copy and one excitation: the endpoint for every n ≥ 17

Ansh Mishra and Aryan Senthilkumar, November 2026.

This folder proves the inequality (⋆) of `two-sector-class/README.md` (section "One free copy and one excitation")
for every m ≥ 16. That gives the endpoint for the whole one-free-copy, one-excitation class from n = m + 1 = 17 copies
on.

## Setting

Take m sector copies with planes K_s and e = e_1 ⊗ … ⊗ e_m with e_s ⊥ K_s. Put A = K_1 ⊗ … ⊗ K_m, and let
V = C^3 ⊗ (A ⊕ C e), with one free copy. For C of rank at most two with range and co-range in V, write
C = [[X, P], [Q, Y]] in blocks. The inequality is

    (⋆)   √q_1(Y) · (√Γ_X − √Σ_full) ≤ (2^m − 2) √(q̃(P) q̃(Q)),

with Γ_X, Σ_full, q_1 and q̃ as in `two-sector-class/README.md`. There, (⋆) is shown to imply q_n(C) ≥ 0 on the whole
class. By the compression of Corollary B the same holds for range in V with any co-range.

## Result

**Theorem G.**
- (⋆) holds for every m ≥ 16, for Y of either rank.
- (⋆) holds for every m ≥ 11 when rank Y ≤ 1.

Hence q_n(C) ≥ 0 for every rank-two C in the one-free-copy, one-excitation class when n ≥ 17, and when n ≥ 12 if the
excitation components give rank Y ≤ 1.

Together with Theorem E (range in C^3 ⊗ span{c, e}) and Theorem F (case 1a, m ≥ 4), this leaves the class open only
at 3 ≤ m ≤ 15. There we know numerically that (⋆) holds with room to spare. The largest ratios found are 4.06, 4.24,
4.39 and 4.49 at m = 3, 4, 5, 6, against the needed 6, 14, 30 and 62, and all of them lie in case 1a.

## Structure of the proof

1. **Recursion** (`proof/star19_energy_aware_recursion.md`). For arbitrary fields and either rank of Y, the best
   constant K_m in (⋆) satisfies K_m ≤ √2·K_{m−1} + B, provided a "merged-sector bridge" MB_B holds. Here K_1 = 0.
   The recursion uses:
   - a separable six-atom frame for 2I − F;
   - the exact splitting Γ̄ = Σ_j w_j Γ_{m−1}(X_j) and Σ_m = Σ_j w_j Σ_{m−1}(X_j);
   - atoms C_j of rank at most two with the same Y and Σ_j w_j p_j q_j ≤ 2pq;
   - a Minkowski triangle inequality.
2. **Bridge** (`proof/orb37_flat_margin_functional_bridge.md`). MB_B holds with B = 125. When rank Y ≤ 1 it holds with
   B = (163/10)√2 < 24.
   - **Rank Y ≤ 1.** A gauge (`proof/induction17_rank_one_Y_bridge.md`, `proof/orb28_combined_uniform_bridges.md`)
     and a functional bound with constant 163/10. That bound combines the sharp flat margins
     (`proof/induction16_qubit_merged_gap.md`, `proof/orb33_universal_merged_flat_margin.md`), a curvature bound
     (`proof/orb27_full_safe_functional.md`) and the pure weighted-kernel bounds
     (`proof/star20_weighted_qubit_kernel.md`).
   - **Rank Y = 2.** A Schur gauge and adjugate bound (`proof/star22_uniform_merger_bridge.md`). For small amplitudes
     the functional bound applies. For large amplitudes, a two-step energy path with the determinant bound of
     `proof/orb9_internal_concurrence.md` and the concurrence bound of
     `proof/orb26_arbitrary_gram_coherent_triangle.md`. The latter rests on the free rectangle
     (`proof/star24_free_rectangle_reduction.md`, `proof/star25_free_completion_rectangle.md`,
     `proof/induction24_free_completion_uniform_bridge.md`), on the coherent addition of separable Grams
     (`proof/orb24_coherent_separable_sum.md`), and on the component bounds of
     `proof/star21_concurrence_rectangle.md`.
3. **Arithmetic.** With K_1 = 0 and B = 125, K_m ≤ 2^m − 2 first holds at m = 16, since K_16 = 125(255 + 127√2) <
   65534 while K_15 > 32766. If K_m ≤ 2^m − 2, then √2·K_m + 125 ≤ 3·2^{m−1} + 122 ≤ 2^{m+1} − 2 for m ≥ 8. So the
   bound propagates to every m ≥ 16. With B = (163/10)√2 the first m is 11.

The one external input is the exactness of the homogeneous two-qubit concurrence dual (Eisert, Brandão and Audenaert,
New J. Phys. 9, 46 (2007), Eq. 66), used in the last step of `proof/orb24_coherent_separable_sum.md`.

The notes in `proof/` are our working notes, kept close to how the argument was found. Each one also records
statements outside this theorem. Only the steps listed above are part of the claim.

## Verification

An independent referee re-derived every step of the chain above from the definitions. It checked each step with its
own code (`code/`, outputs in `logs/`):
- **Conventions.** The definitions and Proposition M reproduce to 2·10⁻¹⁵. The exact m = 3 operator of
  `two-sector-class/code/check_star_constant.py` gives Γ_X·2^{40} = 1358550811983 and Σ_full·2^{40} = 1308487705385.
- **Adversarial searches.** Each lemma was searched with 12–40 restarts at internal dimensions 1–8. None was
  violated, and the margins are large. The bridge supremum found is about 4.0, against 125, or 23.05 when rank
  Y ≤ 1. The functional bound reaches 2.82 against 16.3, and the curvature bound 0.50 against 1. A summary is in
  `logs/results_final.txt`.
- **Exact arithmetic.** The rational certificates and the thresholds m = 16 and m = 11 were checked in exact
  arithmetic (`logs/stage6_arith.log`).
- **The theorem itself at small m.** Searches at m = 3–6 gave the ratios quoted under Result.

The referee found no error and no unproved step. It checked the concurrence dual numerically, with deviation below
7.5·10⁻⁸·Tr B on 52 states, and did not re-prove it. Two of the inputs are not load-bearing for B = 125: replacing
the double-residual bound of `orb26` by the elementary one of `star21` still gives 119.5 < 125. The coherent-addition
theorem of `orb24` is needed; without it the constant would be about 152. Not journal peer-reviewed.

The code uses Python 3 with NumPy, SymPy and PyTorch.
