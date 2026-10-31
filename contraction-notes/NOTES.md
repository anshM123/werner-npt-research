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

**Corollary.** At t = 1/2, the endpoint inequality holds on the aligned class for four copies:
- in dimensions three and four, by the three-copy theorem [Zen; Liu; MB];
- in every dimension, given Tomar's computer-assisted three-copy theorem [T1].

On the class, n copies follow from the full (n − 1)-copy statement. The argument is one step only: five copies would
need the full four-copy statement.

**Corollary 3 (a nested class, every n).** Let K_1 be the operators of rank at most two on one copy. For n ≥ 2 let
K_n be the operators that, after a unitary on each copy and a permutation of the copies, have the form above with
tail operator X in K_{n−1}. Then q_n(t, C) ≥ 0 for every C in K_n, every n ≥ 1 and every 0 ≤ t ≤ 1/2.

*Proof.* Induction on n. For n = 1, q_1(t, X) ≥ (1 − 2t)‖X‖² because |Tr X|² ≤ 2‖X‖² for rank two. For n ≥ 2,
Theorem 1 gives q_n(t, C) ≥ (1 − t) q_{n−1}(t, X) ≥ 0. ∎

Corollary 3 uses no unproven input. Random members of K_n are thick at every copy: the reductions of the range and of
the co-range have rank three, and τ_s = Tr(P_s Q_s) = 3 > 2, at every copy we checked (`check_nested.py`, n ≤ 4). So
thin-copy removal reduces nothing here. K_n is a measure-zero family, but it contains P ⊗ |c⟩⟨c|^{⊗(n−1)}.

**Remarks.**
- **The class is not covered by earlier results.** A generic member has three-dimensional local supports on the
  first copy (so the copy is thick, and thin-copy removal [WZ, Thm. 7.2; T2, Lemma 6.1] does not apply). Its range
  contains no product vector unless c_0 α_0 + c_1 α_1 ∈ span(a) for some (c_0, c_1) ≠ 0, so the product-ray theorem
  [WZ, Thm. 7.3] does not apply either. Tomar's every-copy class [T2, Thms. 4.1-4.7] has GHZ backbones |0^k⟩, |1^k⟩
  with tails at single copies, which is a different family.
- **Aligned backbones are needed.** The same grouping does not work when the backbones differ (y_k = |k⟩ ⊗ a_k + |2⟩ ⊗ α_k
  with a_0 ≠ a_1). There is still an exact identity, but its part outside the squares reaches −‖C‖² in numerical
  searches, while q_n itself stays nonnegative.

## 2. Local stability of the zeros P ⊗ |a⟩⟨b|

At t = 1/2, q_1(P) = 0 for the rank-two projector P = |0⟩⟨0| + |1⟩⟨1|. So C_0 = P ⊗ |a⟩⟨b| is a zero of q_n for all
vectors a, b on the other copies, including entangled ones.

**Proposition 2.** Let n ≥ 2, and assume q_{n−1}(1/2, ·) ≥ 0 on operators of rank at most two. Then along every smooth
curve C(ε) = (Y_0 + εY_1 + …)(Z_0 + εZ_1 + …)^* of rank-two operators through C_0 = Y_0 Z_0^*, the second derivative of
q_n(1/2, C(ε)) at ε = 0 is nonnegative.

*Proof sketch.*
- The gradient of q_n at C_0 is B(C_0) = −P^⊥ ⊗ B'(|a⟩⟨b|), where P^⊥ = I − P on the first copy and B' is the map B
  for the other n − 1 copies. It is orthogonal
  to the tangent space, so the second variation is

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

Tomar proves positivity near generic equality cases supported on two-dimensional subspaces of every copy [T2]. When a
and b are entangled, the zeros here have higher-rank reductions on the other copies.

## 3. A contraction inequality that fails

A natural way to pass from n − 1 to n copies would be the site-contraction inequality

    (SCI)   R_n(C) ≥ (1 − t) min_{s, w, w'} R_{n−1}(⟨w|_s C|w'⟩_s),   R_k = q_k / ‖·‖².

Since contractions keep rank at most two, SCI for all n would give R_n ≥ (1 − 2t)(1 − t)^{n−1}, the conjectured sharp
curve, and in particular the endpoint for every n.
- SCI is tight at the product zeros for every t, and at the GHZ zeros for t = 1/2.
- At t = 1/2, every zero we generated has a contraction that is a zero (60 random zeros at three copies, 12 at four;
  exploratory runs, not scripted here).
- **The rank-one case holds.** For C = |x⟩⟨y|, I − tF_s ≥ (1 − t)I gives
  q_n(C) ≥ (1 − t) Σ_{i,j} q_{n−1}(⟨i|_s C|j⟩_s).
- **Rank two fails at three copies.** The operator in `data/sci_counterexample_n3.npy` has rank two and singular
  values 0.98 ‖C‖ and 0.20 ‖C‖. At t = 1/2 it has R_3(C) = 0.101905, while the smallest value of R_2 over its contractions that
  we find is 0.207236, so SCI fails by 1.7 × 10^{−3}. Positivity itself is not in question: R_3(C) > 0.
  - This is a numerical finding, from many restarts of an exact alternating minimisation at every site
    (`check_sci_counterexample.py`), not a certificate.
  - Related operators violate SCI at every t we tried in [0.3, 1/2].

Two related observations from our searches:
- **A defect version survives, with a thin margin.** Subtract the conjectured sharp bound (1 − t)^{n−1}(‖C‖_2² − t‖C‖_1²)
  (see the remarks of `first-order-bound/`) from both levels:

      D_n(C)/‖C‖² ≥ (1 − t) min D_{n−1}(X)/‖X‖²,   D_k(X) = q_k(X) − (1 − t)^{k−1}(‖X‖² − t‖X‖_1²).

  - This is consistent with every search we ran. The worst ratio found at three copies is about 0.54, against the
    required 1/2. It is reached near the zeros of Section 2.
  - Appending a pure copy does not change the ratio, so its infimum can only decrease with n.
  - We have no proof.
- **The dephasing condition fails.** The condition "some basis of some copy has Σ_i q_{n−1}(X_ii) ≥ t q_{n−1}(Σ_i X_ii)"
  would close the induction through the recursion, but it fails already at two copies.

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
