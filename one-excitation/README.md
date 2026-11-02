# A thin operator plus one product excitation

Ansh Mishra and Aryan Senthilkumar, 2 November 2026

This note proves the endpoint inequality for every number of copies on a class of rank-two operators. The operators in
this class can have one-copy reductions of rank three at every copy. The class contains the lowest thick critical point
we found at two copies. It is a small class, and the argument does not extend to general operators.

## Setting

For an operator C on (C^3)^⊗n,

    q_n(C) = Σ_{J ⊆ {1..n}} (−1/2)^{|J|} ‖Tr_J C‖²,

and n-copy undistillability of the qutrit Werner state at the endpoint is q_n(C) ≥ 0 for every C of rank at most two.

Fix the following:
- a plane K_s ⊂ C^3 and a unit vector e_s ⊥ K_s at every copy;
- K = K_1 ⊗ … ⊗ K_n, a copy of (C^2)^⊗n;
- e = e_1 ⊗ … ⊗ e_n, so that e ⊥ K;
- V = K ⊕ C e.

**Proposition.** Let n ≥ 2. If C has rank at most two and both ran C and ran C* lie in V, then q_n(C) ≥ 0.

Without the e-component this is the thin case (range in a product of planes; Tomar, Lemma 2.6). With it, the
one-copy reduction of the range has rank three at every copy, unless the e-component is zero. A generic member has no
product vector in its range or co-range, so the product-ray theorem (Wu–Zou, Thm. 7.3, which holds for every n) does
not apply. Tomar's every-copy classes [T2, Thms. 4.1–4.7] are built on the two GHZ strings 0…0 and 1…1 with tails at
single copies. Here the qubit part is all of K, and only the single string e is added.

## Proof

Write C_K = P_K C P_K, κ = P_K C e, κ' = P_K C* e, γ = ⟨e|C|e⟩, so

    C = C_K + |κ⟩⟨e| + |e⟩⟨κ'| + γ |e⟩⟨e|.

**Step 1 (exact identity).**

    q_n(C) = q_n(C_K) + ‖κ‖² + ‖κ'‖² + 2^{−n}|γ|² + 2 (−1/2)^n Re(γ̄ Tr C_K).

A partial trace over a copy s pairs the copy-s factor of the range with that of the co-range. In |κ⟩⟨e| and |e⟩⟨κ'|
that pairing is between e_s and a vector of K_s, which gives zero. So these two pieces contribute only through J = ∅,
and they are orthogonal to the other pieces there. For |e⟩⟨e| every partial trace is |e_{J^c}⟩⟨e_{J^c}|, which gives
Σ_J (−1/2)^{|J|} = 2^{−n}. Against Tr_J C_K this is nonzero only for J = {1..n}, which gives the last term.

**Step 2 (the qubit part).** On M_2 the one-copy map X ↦ X − ½ Tr(X) I kills I and fixes the Pauli matrices. So for X
supported on K,

    q_n(X) = ‖Π X‖²,

where Π is the orthogonal projection onto the span of the 3^n Pauli strings with no identity factor. The rank-one bound
q_n(xy*) ≥ 2^{−n}‖x‖²‖y‖² holds because q_n(xy*) = ⟨x ⊗ y|(I − F/2)^⊗n|x ⊗ y⟩ and I − F/2 ≥ I/2.

**Step 3 (C_K of rank at most one).** Write C_K = xy*, possibly zero. Drop ‖κ‖² + ‖κ'‖² ≥ 0 and minimise over γ:

    q_n(C) ≥ q_n(xy*) − 2^{−n}|⟨y|x⟩|² ≥ 2^{−n}(‖x‖²‖y‖² − |⟨y|x⟩|²) ≥ 0.

**Step 4 (C_K of rank two).** Then rank C = rank C_K, which forces κ ∈ ran C_K, κ' ∈ ran C_K*, and
γ = κ'* C_K⁺ κ. Hence |γ| ≤ ‖κ‖‖κ'‖/s₂, with s₂ = s₂(C_K) the smaller nonzero singular value. Put g = |γ| and
t = Tr C_K. Using ‖κ‖² + ‖κ'‖² ≥ 2‖κ‖‖κ'‖ ≥ 2 s₂ g,

    q_n(C) ≥ q_n(C_K) + 2 s₂ g + 2^{−n} g² − 2^{1−n} g |t| ≥ q_n(C_K) − 2^{−n}(|t| − 2^n s₂)₊².

**Step 5 (qubit lemma).** For X of rank two on (C^2)^⊗n with n ≥ 2,

    2^n q_n(X) ≥ (|Tr X| − 2^n s₂(X))₊².

Take the singular value decomposition X = s₁ u₁v₁* + s₂ u₂v₂*, and let a_P = s₁⟨v₁|P|u₁⟩, b_P = s₂⟨v₂|P|u₂⟩ over
the full-weight Pauli strings P. Then 2^n q_n(X) = ‖a + b‖².
- By the rank-one bound, ‖a‖² = s₁² 2^n q_n(u₁v₁*) ≥ s₁².
- By completeness of all 4^n Pauli strings, ‖b‖² ≤ s₂² 2^n.
- So ‖a + b‖ ≥ s₁ − 2^{n/2} s₂.

Also |Tr X| ≤ s₁ + s₂. If |Tr X| > 2^n s₂, then 0 < |Tr X| − 2^n s₂ ≤ s₁ − (2^n − 1)s₂ ≤ s₁ − 2^{n/2}s₂, since
2^n − 1 ≥ 2^{n/2} for n ≥ 2. Squaring gives the lemma.

Steps 4 and 5 give q_n(C) ≥ 0 when C_K has rank two, and Step 3 covers the rest. ∎

## The two-copy critical point

The lowest thick critical point of R = q_2/‖·‖² on the rank-two variety that we found at two copies is in this class
(`orbit-notes` and `code/critical_n2.npy`). Its value is

    R = (9 − √41)/20 ≈ 0.12984,

and its singular-value ratio is (√41 − 4)/5. Three of four biased searches converged to it. In the eigenbasis of the
one-copy reduction of its range (one local unitary per copy, applied on both sides), it commutes with Z ⊗ Z, where
Z = diag(1, −1, 1). Its range and co-range lie in

    (K ⊗ K) ⊕ C|11⟩,    K = span{|0⟩, |2⟩},

and they differ by a rotation of about 30° inside K at each copy. It is a saddle on the rank-two variety (6 negative
Hessian directions, the smallest −0.58). The Proposition shows that no operator of this shape violates the endpoint
inequality, for any number of copies.

## Checks

`code/check_one_excitation.py` checks (`logs/check_one_excitation.log`, n = 2, 3, 4):
- the identity of Step 1 on 20 random members (error ≤ 2.2·10⁻¹⁶);
- the qubit lemma, by adversarial minimisation of the gap (minimum +2.7·10⁻¹⁸);
- the minimum of q_n/‖C‖² over the class (−3.5·10⁻¹⁶, attained on thin members);
- that the two-copy critical point lies in the class after the local unitary (weight outside 7·10⁻²⁸).

## References

- [T2] T. S. Tomar, *Equality cases and every-copy results for the Werner endpoint*, Zenodo (2026): Lemma 2.6 (range
  in a product of qubits), Theorems 4.1–4.7 (GHZ backbones with single-copy tails).
- [WZ] Wu and Zou, arXiv:2608.02647: Theorem 7.3 (product ray).
