# Three copies of the Werner endpoint: rank-two projectors

Ansh Mishra and Aryan Senthilkumar, 30 October 2026

These notes collect results on the three-copy endpoint form. They are the step after our two-copy thickness note,
in the direction of T. S. Tomar's Conjecture 7.1 (Zenodo 23718037). Results are labelled **proved** or
**conjectured**. Every identity is checked to machine precision by the scripts in `code/`.

## Setting

For an operator C on (C^d)^{⊗3} let

    q_3(C) = Σ_J (−1/2)^{|J|} ‖Tr_J C‖_2²,

with J running over subsets of the three sites. The qutrit Werner state at α = 1/2 is three-copy undistillable if and
only if q_3(C) ≥ 0 for every C of rank at most two. This is known, by Tomar's computer-assisted proof. Wu and Zou
(arXiv:2608.02647) proved it for positive semidefinite and normal C, and for every C whose range and co-range lie in
one plane.

Throughout, V is a plane (a two-dimensional subspace) with orthonormal basis u, v. P_V is its projector, and λ_3(V) is
the smallest, over the three sites, of the third largest eigenvalue of the one-site reduction of P_V. In this
normalisation, Tomar's Conjecture 7.1 at k = 3 asks for q_3(P_V) ≥ λ_3(V)/4.

## 1. Exact identities (proved)

Let α_l, β_l be the one-site marginals (site l) of uu* and vv*. Let τ_l = Tr_{other sites}(uv*), and ψ_l the one-site
marginal of (uu* + vv*)/2.

**(a) Monogamy form, any d.**

    q_3(P_V) = 3/2 − Σ_l C_l,   where   C_l = ¼‖α_l − β_l‖² + ‖τ_l‖².

Purify P_V/2 as ψ = (u|0⟩ + v|1⟩)/√2, with a maximally mixed qubit R. Then C_l = 2‖ψ_{lR} − ψ_l ⊗ I/2‖²_2 is the
Hilbert–Schmidt correlation between R and site l. So the projector case of the three-copy endpoint is a monogamy
inequality: R's correlations with the three sites add up to at most 3/2. The total correlation of R is 3/2. The
bound is tight both for product planes, where one site carries everything, and for span{|000⟩, |111⟩}, where each
site carries 1/2. For four or more sites the single-site sum can exceed 3/2: it reaches n/2.

**(b) Cross-product form, d = 3.** With the triple cross product

    (u ⊠ v)_{pqr} = Σ ε_{pii'} ε_{qjj'} ε_{rkk'} u_{ijk} v_{i'j'k'},

which depends only on V,

    q_3(P_V) = 5/2 − Σ_l Tr ψ_l² − ‖u ⊠ v‖².

**(c) Entanglement form, d = 3.** For x Haar-random in V,

    q_3(P_V) = ¾ ( E_{x∈V} Σ_l S_L(x_l) − ‖u ⊠ v‖² ),   with S_L(ρ) = 1 − Tr ρ².

So the projector case says that the cross-product norm of the plane is at most the average total single-site
entanglement of the states in it. For comparison, the two-copy analogue is q_2(P_V) = 3 E_{x∈V} S_L(x), which is
trivially nonnegative.

The identities follow from ‖Tr_l P_V‖² = ‖α_l‖² + ‖β_l‖² + 2‖τ_l‖², from the Lagrange identity
‖x ⊠_S y‖² = Σ_{R⊆S} (−1)^{|R|} ⟨x⊗y|F_R|x⊗y⟩ for cross products taken at a set S of sites, and from
E[xx* ⊗ xx*] = (P_V ⊗ P_V)(I + F)/6. Check: `code/check_identities.py`.

## 2. A margin on diagonal planes (proved, every d)

**Theorem.** Let V lie in span{|jjj⟩ : j = 1, …, d}, so that P_V = Σ_{jk} G_{jk} |jjj⟩⟨kkk| for a rank-two projector G
on C^d with diagonal g. Then

    q_3(P_V) = 3/2 − ¾ Σ_j g_j²   and   λ_3(V) = g_(3),

where g_(3) is the third largest entry of g. Moreover

    q_3(P_V) ≥ ¾ λ_3(V).

*Proof.* Every partial trace kills the off-diagonal terms of P_V. So all one- and two-site marginals are diagonal,
with squared norm Σ_j g_j², and the formula for q_3 follows. By Schur–Horn, g ranges over
{0 ≤ g_j ≤ 1, Σ g_j = 2}, and the claim becomes Σ g_j² + g_(3) ≤ 2. Order g_1 ≥ g_2 ≥ …, and put t = g_(3) and
r = Σ_{j≥3} g_j, so that Σ_{j≥3} g_j² ≤ t r.
- If g_1 + g_2 = 2 − r ≥ 1, then g_1² + g_2² ≤ 1 + (1 − r)².
  - For r ≤ 1/2 (where t ≤ r) the total is at most 2 − r + 2r² ≤ 2.
  - For r ≥ 1/2 (where t ≤ g_2 ≤ 1 − r) the total is at most 3 − 2r ≤ 2.
- If g_1 + g_2 = s < 1, the total is at most s²/2 + 3s/2 < 2. ∎

Equality holds, with λ_3 > 0, for g = (2/3, 2/3, 2/3, 0, …) (the plane span{GHZ, GHZ_ω}) and for
g = (1, 1/2, 1/2, 0, …) (the plane span{(|000⟩ − |111⟩)/√2, |222⟩}). In Tomar's units this is a ratio of 3 on this
family, against his conjectured 1. Check: `code/check_diagonal_family.py`.

## 3. Conjecture: the sharp three-copy margin for projectors

**Conjecture.** For every plane V ⊂ (C³)^{⊗3}, q_3(P_V) ≥ ¾ λ_3(V).

Evidence (numerical):
- **The two equality planes of Section 2 are local minima** of q_3/λ_3 = 3/4 over all planes. Local searches from
  perturbations at scales 10^{-3} to 0.15 never go lower (`code/check_local_minima.py`).
- **Symmetric planes** (V ⊂ Sym³(C³)): the minimum found is 0.75017.
- **Searches over all planes** stay at 0.758 or above.
- **Searches started next to thin planes** move away from them, so the bound appears linear down to λ_3 = 0.
- **Only the minimum over sites works.** The analogous bound with max_l λ_3^{(l)} fails: there are exact zeros of
  q_3(P_V) that are thick at two sites and thin at the third.

Equivalent forms, from Section 1:
- monogamy: Σ_l C_l ≤ (3/2)(1 − min_l λ_3(ψ_l));
- cross product: ‖u ⊠ v‖² + Σ_l Tr ψ_l² ≤ 5/2 − ¾ λ_3(V);
- entanglement: E_x Σ_l S_L(x_l) − ‖u ⊠ v‖² ≥ λ_3(V).

Remarks on a proof:
- **Wu and Zou's proof of q_3 ≥ 0 has three inequalities:** their Theorem 1.1, their Lemma 5.1, and an AM–GM step.
  The total slack of these three equals the left side of the conjecture, up to a factor.
- **No single step carries the margin.** At span{GHZ, GHZ_ω} the whole margin sits in the slack of Theorem 1.1. At
  span{(|000⟩ − |111⟩)/√2, |222⟩} it sits in the AM–GM step. A proof must bound the three slacks jointly.
- **The best simple bound is short of the target.** A correlation-complementarity argument (anticommuting
  σ_i^R ⊗ O operators, with traceless HS-unit O of operator norm ≤ √((d−1)/d)) gives Σ_l C_l ≤ 3(d − 1)/d. That is
  exactly 3/2 at d = 2, but only 2 at d = 3.

For Hermitian rank-two C with range in a plane (not only projectors), the numerical minimum of
2³ q_3(C)/(‖C‖² λ_3) is about 1.40, again above Tomar's conjectured 1.

## 4. Remark: crossed Cauchy–Schwarz is not tensor-stable

For every n, the endpoint for operators whose range and co-range lie in one plane (in particular every normal C) is
equivalent to the following. Tr[Z_n (X ⊗ X)] ≥ 0 for every Hermitian X of rank ≤ 2, where Z_n = (F − 1/2)^{⊗n}.
This is the same as a Cauchy–Schwarz inequality h(u,v)² ≤ h(u,u) h(v,v) for h(u,v) = ⟨u⊗v|Z_n|u⊗v⟩, and also the same
as asking that every two-dimensional local compression of Z_n be a positive combination of A ⊗ A.

If this property passed to tensor products, it would give the plane class for every n. **It does not.**
- **A qutrit counterexample.** The 9×9 operator Z1 in `data/` satisfies the property, with margin 4.9 × 10^{-3} in
  local searches, and F Z1 ≥ 0. But Z1 ⊗ (F − 1/2) violates it by about 1.6 × 10^{-2}. The violation is the same with
  the qubit Werner factor.
- **An example invariant under local unitaries.** Z = Z_W ⊗ Z_W + ¾ I ⊗ I on two qutrits satisfies the property,
  since Tr[Z(X⊗X)] = q_2(X) + ¾(Tr X)² ≥ 0, and F Z ≥ 0. But Z ⊗ Z_W violates it, with minimum about −1/4. This Z
  even satisfies full two-copy positivity: Tr[Z(C⊗C*)] ≥ 0 for every rank-two C.

So any all-copy argument for the plane class must use the specific Werner weights. Check:
`code/check_tensor_stability.py`. Its CCS margins come from local optimisation, so they are numerical, not certified.

## References

- T. S. Tomar, Equality cases and every-copy results for the Werner endpoint inequality, Zenodo 23718037 (2026),
  Conjecture 7.1 and Theorem 6.3.
- T. S. Tomar, Undistillability of Werner states from three, four and five copies, Zenodo 23718035 (2026).
- T. Wu, Q. Zou, Sharp Plücker geometry for three-copy Werner distillation, arXiv:2608.02647.
