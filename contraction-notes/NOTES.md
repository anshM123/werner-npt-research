# Contractions and the n-copy Werner endpoint

Ansh Mishra and Aryan Senthilkumar, October 2026.

These notes record small positive results and one negative one about passing the endpoint inequality from n − 1
copies to n copies. The main one is an explicit class, thick at every copy, on which the inequality holds for every n. They do not settle the all-copy problem.

## Setting

For an operator C on (C^d)^{⊗n} (d ≥ 3) and 0 ≤ t ≤ 1/2 put

    q_n(t, C) = Σ_J (−t)^{|J|} ‖Tr_J C‖_2² = ⟨C, B(C)⟩,   B = Π_s (1 − t E_s),   E_s X = I_s ⊗ Tr_s X,

the sum running over all subsets J of the copies. The Werner state ρ ∝ I − F/2 is n-copy undistillable if and only if
q_n(1/2, C) ≥ 0 for every C of rank at most two [CR23, Thm. 1; QGWG, Thm. V.1; SC, Thm. 1]. Two facts are used
throughout.

- **One-site recursion.** For any orthonormal basis |i⟩ of the first copy, with X_ij = ⟨i|_1 C |j⟩_1,

      q_n(t, C) = Σ_{i,j} q_{n−1}(t, X_ij) − t q_{n−1}(t, Σ_i X_ii).

  Every contraction ⟨w|_1 C|w'⟩_1 of a rank-two C has rank at most two.
- **Crossed pairing.** For rank-one operators the polarised form is

      q_{n−1}(t; |x⟩⟨y|, |u⟩⟨v|) = ⟨x ⊗ v| W |u ⊗ y⟩,   W = Π_s (I − t F_s) ≥ (1 − t)^{n−1},

  where F_s exchanges the s-th factors of two copies (the swap trick).

## 1. An aligned class

Fix orthonormal |0⟩, |1⟩, |2⟩ in the first copy. Take any vectors a, b, α_0, α_1, β_0, β_1 on the other n − 1 copies,
and put

    y_k = |k⟩ ⊗ a + |2⟩ ⊗ α_k,   z_k = |k⟩ ⊗ b + |2⟩ ⊗ β_k   (k = 0, 1),   C = Σ_k |y_k⟩⟨z_k|,   X = Σ_k |α_k⟩⟨β_k|.

The two range vectors share the "backbone" a on the other copies, and the two co-range vectors share b. The "tails"
α_k, β_k are arbitrary. With a = b = |c⟩^{⊗(n−1)} and α = β = 0 one gets P ⊗ |c⟩⟨c|^{⊗(n−1)}, a zero at t = 1/2.

**Theorem 1.** For every n ≥ 2 and 0 ≤ t ≤ 1/2,

    q_n(t, C) = (1 − t) q_{n−1}(t, X) + (2 − 4t) q_{n−1}(t, |a⟩⟨b|) + (1 − 4t²) Σ_k q_{n−1}(t, |a⟩⟨β_k|)
                + Σ_k ⟨v_k, W v_k⟩,   v_k = α_k ⊗ b − 2t a ⊗ β_k.

Hence q_n(t, C) ≥ 0 whenever q_{n−1}(t, ·) ≥ 0 on operators of rank at most two. The three rank-one terms are
nonnegative for every t ≤ 1/2, X has rank at most two, and W ≥ 0.

*Proof.* In the recursion above, X_kk' = δ_kk' |a⟩⟨b| for k, k' ∈ {0, 1}, X_2k = |α_k⟩⟨b|, X_k2 = |a⟩⟨β_k| and X_22 = X.
So

    q_n = 2q(A) + q(X) + Σ_k [q(|α_k⟩⟨b|) + q(|a⟩⟨β_k|)] − t q(2A + X),   A = |a⟩⟨b|,

with q = q_{n−1}(t, ·). Expand q(2A + X) = 4q(A) + 4 Re q(A, X) + q(X), and write Re q(A, X) = Σ_k Re ⟨a ⊗ β_k|W|α_k ⊗ b⟩
by the crossed pairing. Since q(|x⟩⟨y|) = ⟨x ⊗ y|W|x ⊗ y⟩, the terms q(|α_k⟩⟨b|) + 4t² q(|a⟩⟨β_k|) − 4t Re ⟨a ⊗ β_k|W|α_k ⊗ b⟩
combine into ⟨v_k, W v_k⟩, which leaves the stated remainder. ∎

**Corollary 2.** At t = 1/2, the endpoint inequality holds on the aligned class for four copies:
- in dimensions three and four, by the three-copy theorem [Zen; Liu; MB];
- in every dimension, given Tomar's computer-assisted three-copy theorem [T1].

On the class, n copies follow from the full (n − 1)-copy statement. The argument is one step only: five copies would
need the full four-copy statement.

**Corollary 3 (a nested class, every n).** Let K_1 be the operators of rank at most two on one copy. For n ≥ 2 let
K_n be the operators that, after a conjugation C ↦ UCU* by a product unitary U = U_1 ⊗ ⋯ ⊗ U_n (the same on range and
co-range) and a permutation of the copies, have the form above with tail operator X in K_{n−1}. Then q_n(t, C) ≥ 0 for every C in K_n, every n ≥ 1 and every 0 ≤ t ≤ 1/2.

*Proof.* Induction on n. For n = 1, q_1(t, X) ≥ (1 − 2t)‖X‖² because |Tr X|² ≤ 2‖X‖² for rank two. For n ≥ 2,
Theorem 1 gives q_n(t, C) ≥ (1 − t) q_{n−1}(t, X) ≥ 0. ∎

Corollary 3 uses no unproven input. Random members of K_n are thick at every copy: the reductions of the range and of
the co-range have rank three, and τ_s = Tr(P_s Q_s) = 3 > 2, at every copy we checked (`check_nested.py`, n ≤ 4). So
thin-copy removal reduces nothing here. K_n is a measure-zero family, but it contains P ⊗ |c⟩⟨c|^{⊗(n−1)}. By Theorem 1,
the only nonzero zeros of K_n at t = 1/2 are the operators P ⊗ |a⟩⟨b| (up to local unitaries and permutations), and
q_n > 0 on K_n ∖ {0} for t < 1/2, so the class has no new equality cases. K_n is new from n ≥ 4; for n ≤ 3 it follows
from the two- and three-copy theorems.

**Remarks.**
- **The class is not covered by earlier results.** A generic member has three-dimensional local supports on the
  first copy (so the copy is thick, and thin-copy removal [WZ, Thm. 7.2; T2, Lemma 6.1] does not apply). Its range
  contains no product vector unless c_0 α_0 + c_1 α_1 ∈ span(a) for some (c_0, c_1) ≠ 0, so the product-ray theorem
  [WZ, Thm. 7.3] does not apply either. Tomar's every-copy class [T2, Thms. 4.1-4.7] has GHZ backbones |0^k⟩, |1^k⟩
  with tails at single copies, which is a different family. Generic members of K_n are also outside the every-n classes
  of `two-sector-class/` and `compression-notes/`: for n ≥ 5 every two-copy reduction of their range and co-range has
  rank 9, while each of those classes forces a smaller rank on some pair of copies (at most 5, 6 or 3).
- **Aligned backbones are needed.** The same grouping does not work when the backbones differ (y_k = |k⟩ ⊗ a_k + |2⟩ ⊗ α_k
  with a_0 ≠ a_1). There is still an exact identity, but its part outside the squares reaches −‖C‖² in numerical
  searches, while q_n itself stays nonnegative.

## 2. Local stability of the zeros P ⊗ |a⟩⟨b|

At t = 1/2, q_1(P) = 0 for the rank-two projector P = |0⟩⟨0| + |1⟩⟨1|. So C_0 = P ⊗ |a⟩⟨b| is a zero of q_n for all
vectors a, b on the other copies, including entangled ones.

**Proposition 4.** Let n ≥ 2, and assume q_{n−1}(1/2, ·) ≥ 0 on operators of rank at most two. Then along every smooth
curve C(ε) = (Y_0 + εY_1 + …)(Z_0 + εZ_1 + …)^* of rank-two operators through C_0 = Y_0 Z_0^*, the second derivative of
q_n(1/2, C(ε)) at ε = 0 is nonnegative.

*Proof sketch.*
- The derivative of q_n at C_0 in direction H is 2Re⟨B(C_0), H⟩, with B(C_0) = −P^⊥ ⊗ B'(|a⟩⟨b|), where
  P^⊥ = I − P on the first copy and B' is the map B for the other n − 1 copies. B(C_0) is orthogonal to the tangent
  space, so the second variation (the ε² coefficient) is

      q_n(T) + 2 Re⟨B(C_0), Y_1 Z_1^*⟩,   T = Y_1 Z_0^* + Y_0 Z_1^*.

- Split the first-copy parts of the columns of Y_1 and Z_1 into the span of |0⟩, |1⟩ and its complement:
  y_{1,k} = Σ_i |i⟩ α_k^i and z_{1,k} = Σ_i |i⟩ β_k^i.
- The blocks |j⟩⟨k| with j ≥ 2 are traceless and mutually orthogonal, so

      q_n(T) = q_n(T_VV) + Σ_{k, j≥2} [q_{n−1}(|α_k^j⟩⟨b|) + q_{n−1}(|a⟩⟨β_k^j|)].

- The curvature term equals −2 Σ_{k, j≥2} Re⟨α_k^j ⊗ b|W|a ⊗ β_k^j⟩. Together with the last sum it gives
  Σ_{k, j≥2} ‖W^{1/2}(α_k^j ⊗ b − a ⊗ β_k^j)‖² ≥ 0.
- T_VV is supported on the qubit span{|0⟩, |1⟩} at the first copy. There B_1 acts as the projection onto traceless
  operators, so q_n(T_VV) = ½ Σ_σ q_{n−1}(Tr_1[σ T_VV]) over the three Pauli matrices σ.
- Each Tr_1[σ T_VV] = |α_σ⟩⟨b| + |a⟩⟨β_σ| has rank at most two, so these terms are ≥ 0 by assumption. ∎

The script `check_zero_family.py` evaluates the Hessian for random entangled a, b at n = 3, 4 and finds no negative
eigenvalue. An adversarial search over a, b (not part of the scripts here) found none either. The statement is only
second order.

Proposition 4 is the second-order part of a known result. Tomar [T2, Thm. 5.2] proves that, given EP_{n−1}(d), the
endpoint form is nonnegative on a whole neighbourhood of every product plane U ⊗ χ, with χ any vector on the other
copies, entangled or not, and arbitrary range vectors. The co-range of P ⊗ |a⟩⟨b| is U ⊗ b, so the proposition follows.
By [T1] this is unconditional at four copies. We include only a short direct computation of the second variation; its
split into a qubit part and a complement part parallels Tomar's Lemma B.9.

## 3. A contraction inequality that fails

A natural way to pass from n − 1 to n copies would be the site-contraction inequality

    (SCI)   R_n(C) ≥ (1 − t) min_{s, w, w'} R_{n−1}(⟨w|_s C|w'⟩_s),   R_k = q_k / ‖·‖².

Since contractions keep rank at most two, SCI for all n would give R_n ≥ (1 − 2t)(1 − t)^{n−1}, the conjectured sharp
curve, and in particular the endpoint for every n.
- SCI holds with equality at P ⊗ |c⟩⟨c|^{⊗(n−1)} for every t (a zero only at t = 1/2), and at the GHZ zeros for t = 1/2.
  The minimum runs over nonzero contractions.
- At t = 1/2, every zero we generated has a contraction that is a zero (60 random zeros at three copies, 12 at four;
  exploratory runs, not scripted here).
- **The rank-one case holds.** For C = |x⟩⟨y|, I − tF_s ≥ (1 − t)I gives
  q_n(C) ≥ (1 − t) Σ_{i,j} q_{n−1}(⟨i|_s C|j⟩_s).
- **Rank two fails at three copies.** The operator in `data/sci_counterexample_n3.npy` has rank two and singular
  values 0.98 ‖C‖ and 0.20 ‖C‖. At t = 1/2 it has R_3(C) = 0.101905, while the smallest value of R_2 over its contractions that
  we find is 0.207236, so SCI fails by 1.7 × 10^{−3}. Positivity itself is not in question: R_3(C) > 0.
  - The minimum over contractions is also verified by a decomposition: at each site, M_s = Q_s − 0.2072 H_s =
    P_s + Q_s^Γ with Q_s ⪰ 0 and P_s ≻ 0 (an SDP solution checked in double precision, with relative margins
    3.5·10⁻³, 1.0·10⁻³ and 1.2·10⁻⁴). So every contraction has R_2 ≥ 0.2072, and SCI fails by about 1.69·10⁻³.
  - The stored operator violates SCI for t ≥ 0.4; small perturbations of it violate SCI at t = 0.30 and 0.35 (checked
    the same way).

Two related observations from our searches:
- **A defect version survives, with a thin margin.** Subtract the conjectured sharp bound (1 − t)^{n−1}(‖C‖_2² − t‖C‖_1²)
  (see the remarks of `first-order-bound/`) from both levels:

      D_n(C)/‖C‖² ≥ (1 − t) min D_{n−1}(X)/‖X‖²,   D_k(X) = q_k(X) − (1 − t)^{k−1}(‖X‖² − t‖X‖_1²).

  - This is consistent with every search we ran. The worst ratio found at three copies is about 0.54, against the
    required 1/2. It is reached near the zeros of Section 2.
  - Appending a pure copy does not change the ratio, so its infimum can only decrease with n.
  - We have no proof.
- **The dephasing condition is false at two copies.** The condition "some basis of some copy has
  Σ_i q_{n−1}(X_ii) ≥ t q_{n−1}(Σ_i X_ii)" would close the induction through the recursion, together with EP_{n−1}.
  At t = 1/2 it fails for C = Y Z* with the Gaussian-integer factors in `code/check_dephasing_counterexample.py`:
  q_2(C) = 14193/4 > 0 (‖C‖² = 20969), but in every orthonormal basis the dephasing defect
  Σ_i q_1(X_ii) − ½ q_1(Tr_s C) is at most −65/12 at the first copy and at most −3/4 at the second.
  - *Certificate.* If H_1, H_2 are Hilbert–Schmidt orthonormal and span the traceless diagonal matrices of the basis,
    the defect equals Σ_j q_1(Tr_s[(H_j ⊗ I)C]) − q_1(Tr_s C)/6. Expand H_j = Σ_a r_{j,a} K_a in a basis K_a of
    traceless Hermitian matrices, with Gram matrix G. Then Σ_j r_jᵀ G r_j = 2, and since H_1² + H_2² = (2/3)I,
    Σ_j r_jᵀ D_c r_j = 0 for D_c[a, b] = Re Tr K_c(K_a K_b + K_b K_a)/2. So the defect is at most 2λ − q_1(Tr_s C)/6
    whenever λG − A + Σ_c y_c D_c ⪰ 0, where A[a, b] = Re q_1(Tr_s[(K_a ⊗ I)C], Tr_s[(K_b ⊗ I)C]). The script checks,
    in exact arithmetic, explicit integers λ and y_c for which twice this matrix is an integer matrix with nonnegative
    diagonal-dominance surpluses.
  - A direct numerical maximisation over bases gives about −619 and −565, far below zero.
  - This confirms the two-copy failure found in our searches of 30 October. The referee's min-over-C, max-over-basis
    searches had missed it; such searches drift to the zeros of q_2, where the defect is 0.

  So the dephasing condition is not a valid induction step. The endpoint itself is not in question: q_2(C) > 0.

## References

- [CR23] P. Costa Rico, New partial trace inequalities and distillability of Werner states, Lett. Math. Phys. 115, 47
  (2025); arXiv:2310.05726.
- [QGWG] S.-Y. Qi, G. Gupur, Y.-C. Wu, G.-P. Guo, Phys. Rev. A 110, 012406 (2024); arXiv:2402.18037.
- [SC] Z. Song, L. Chen, A partial-trace matrix inequality and Werner-state distillability, arXiv:2607.23416.
- [WZ] T. Wu, Q. Zou, Sharp Plücker geometry for three-copy Werner distillation, arXiv:2608.02647.
- [T1] T. S. Tomar, Undistillability of Werner states from three, four and five copies, Zenodo 23718035 (2026).
- [T2] T. S. Tomar, Equality cases and every-copy results for the Werner endpoint inequality, Zenodo 23718037 (2026).
- [MB] I. MacCormack, N. Bao, Three copies of the Werner state: the sharp curve in dimensions three and four and a
  bound in every dimension (2026), github.com/physical-superintelligence/oqp-batch-2026.
- [Zen] Anonymous, Three-copy undistillability of the boundary qutrit Werner state (2026), Zenodo 10.5281/zenodo.23134454.
- [Liu] A.-X. Liu, C.-F. Qiao, Three-copy undistillability of the boundary two-ququart Werner state (2026),
  github.com/GLiuax/d4-three-copy-werner.

## Code

Python 3 with NumPy. Outputs are in `logs/`.
- `code/werner.py`: shared helpers.
- `code/check_aligned.py`: Section 1. It checks the crossed pairing and the identity of Theorem 1 on random members,
  for n = 2, 3, 4 and t ∈ {0.1, 0.3, 0.45, 0.5}.
- `code/check_nested.py`: Corollary 3. It checks positivity and thickness of random members of K_n, for n = 2, 3, 4.
- `code/check_zero_family.py`: Section 2. It computes Hessians at P ⊗ |a⟩⟨b| for random entangled a, b, at n = 3, 4.
- `code/check_sci_counterexample.py`: Section 3. It re-checks the stored operator against SCI.
- `code/check_dephasing_counterexample.py`: Section 3. It recomputes q_2 and the dephasing certificate in exact
  arithmetic (SymPy) and maximises the defect over bases numerically.

## Referee

An independent referee re-derived Theorem 1, Corollaries 2 and 3, Proposition 4 and the SCI counterexample with its
own code. It checked:
- the identity of Theorem 1 at d = 3, 4 and n ≤ 5, to 2.4·10⁻¹⁵;
- W ≥ (1 − t)^{n−1}, exactly;
- K_n at n ≤ 6: positivity, thickness, and adversarial minima ≥ 0 against controls that go negative.

It found no mathematical error. Its corrections are applied above:
- Proposition 4 is a special case of [T2, Thm. 5.2];
- the dephasing remark was not reproduced (a later exact certificate shows that the dephasing condition does fail at
  two copies; see Section 3);
- the definition of K_n now fixes the unitary;
- the SCI minimum is now verified by a decomposition.

