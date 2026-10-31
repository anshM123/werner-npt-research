# Compressions, black-box averages and large-n checks for the Werner endpoint

Ansh Mishra and Aryan Senthilkumar, October 2026.

These notes collect several small results about the n-copy endpoint inequality, with code:
- a conjecture that would settle every n, with the evidence for it;
- a proposition explaining why arguments that do not adapt to the operator stop at the PPT boundary;
- an exact form of the projector case as a monogamy relation;
- a transfer-channel description of translation-invariant MPS filters, with a large-n search;
- an exact activation witness for two-pair Werner mixtures;
- an every-n class built from sums of two block-product operators.

**Nothing here settles the all-copy problem.**

## Setting

As in `contraction-notes/`, for an operator C on (C^d)^{⊗n},

    q_n(t, C) = Σ_J (−t)^{|J|} ‖Tr_J C‖_2² = ⟨C, B(C)⟩,   B = Π_s (1 − t E_s),   E_s X = I_s ⊗ Tr_s X,

and the qutrit Werner state ρ ∝ I − F/2 is n-copy undistillable iff q_n(1/2, C) ≥ 0 for every C of rank at most two
("EP_n"). We use the following tools.
- **One-site recursion.** With X_ij = ⟨i|_s C|j⟩_s in an orthonormal basis of copy s,

      q_n(C) = Σ_{ij} q_{n−1}(X_ij) − t q_{n−1}(Tr_s C).

- **Crossed pairing.** q_{n−1}(t; |x⟩⟨y|, |u⟩⟨v|) = ⟨x ⊗ v|W|u ⊗ y⟩ with W = Π_s (I − t F_s) ≥ 0.
- **Thin-copy removal** [WZ, Thm. 7.2; T2, Lemma 6.1]. If the rows and columns of C at some copy lie in a plane, then
  q_n(1/2, C) ≥ 0 follows from EP_{n−1}.

## 1. A compression conjecture

For a copy s and a unit vector w ∈ C³, let P = I − |w⟩⟨w| act on copy s, and put C_{w⊥} = P C P. Then C_{w⊥} has rank
at most two and is supported on the plane w^⊥ at copy s, so thin-copy removal gives q_n(1/2, C_{w⊥}) ≥ 0 from
EP_{n−1}. The definitions are:
- X_{uv} = ⟨u|_s C|v⟩_s and X_w = X_{ww};
- row_s(w) = Σ_l q_{n−1}(X_{wl}) and col_s(w) = Σ_l q_{n−1}(X_{lw}), summed over an orthonormal basis of copy s;
- q_{n−1}(·, ·) is the polarised form.

**Proposition 1.** For every operator C, copy s, unit w and t,

    Δ_s(w) := q_n(t, C) − q_n(t, C_{w⊥}) = row_s(w) + col_s(w) − (1 − t) q_{n−1}(X_w) − 2t Re q_{n−1}(Tr_s C, X_w).

The map w ↦ row_s(w) + col_s(w) − 2t Re q_{n−1}(Tr_s C, X_w) is a Hermitian form w^* M_s w, and Tr M_s = 2 q_n(t, C).

*Proof.* Take an orthonormal basis e_0, e_1, e_2 = w of copy s and write T = Tr_s C. The one-site recursion gives:
- q_n(C) = Σ_{ij} q(X_ij) − t q(T);
- q_n(C_{w⊥}) = Σ_{i,j ≤ 1} q(X_ij) − t q(T − X_w).

The blocks with i = 2 or j = 2 contribute row(w) + col(w) − q(X_w). Also q(T) − q(T − X_w) = 2 Re q(T, X_w) − q(X_w),
which gives the formula. Summing the quadratic part over an orthonormal basis gives 2 Σ_{ij} q(X_ij) − 2t q(T) =
2 q_n(C). ∎

**Conjecture QC.** At t = 1/2, for every n ≥ 2 and every C of rank at most two on (C³)^{⊗n}, there are a copy s and a
unit w with Δ_s(w) ≥ 0, that is q_n(C) ≥ q_n(C_{w⊥}).

QC at n copies and EP_{n−1} together give EP_n, so QC for all n gives the endpoint for all n. QC is tight on the known
zeros: if C is thin at copy s, then w ⊥ (its plane) gives C_{w⊥} = C. The trace identity shows QC is not a
consequence of an averaging: Σ_k Δ_s(e_k) = 2 q_n(C) − (1 − t) Σ_k q_{n−1}(X_{e_k}), and this can be negative while
q_n(C) > 0.

**Evidence** (`check_qc.py` and our search logs).
- Random rank-two operators at n = 2, 3, 4: max_{s,w} Δ_s(w) > 0 in every case.
- Derivative-free adversarial minimisation of max_{s,w} Δ_s(w)/‖C‖² on the exact formula: n = 2 (six runs, smallest
  value +6 × 10^{−3}) and n = 3 (three runs, smallest +0.19).
- Gradient (Danskin) adversarial minimisation: n = 3 (nine runs) and n = 4 (twenty-two runs). The value stays positive.
  It tends to 0 only as the adversary is driven to zeros of q_n, where Δ ≈ q_n → 0⁺.
- The ratio Δ/q_n at the optimum can be small (0.007 at n = 3, 0.03 at n = 4), so the margin is relative, not
  uniform. This happens only near thin zeros.
- Gradient adversaries penalised to keep every copy thick (smallest local eigenvalue at least about 0.02-0.1)
  find max Δ ≈ q_n/2 at n = 2 and n = 3. QC has a robust relative margin on thick operators.

**Canonical directions that fail.** A proof would be easier with a rule that picks w from C. Two natural rules fail
at two copies (`data/`):
- the eigenvector of Tr_{rest}(CC^* + C^*C) at copy s with the smallest eigenvalue: value −0.104 at an operator with
  q_2 = 0.289, where the best w gives +0.283;
- the top eigenvector of M_s: value −0.0077, where the best w gives +0.0062.

**Remark (first zeros).** Suppose EP_{n−1} holds and EP_n fails. Then [T2, Thm. 6.3] gives t* < 1/2 and C* of rank two
minimising q_n(t*, ·)/‖·‖² with value 0. The filtered operators (Q + aP)C*(Q + bP), with P = |w⟩⟨w| at copy s and
Q = I − P, have rank at most two. So the value at a = b = 1 is a double zero, and expanding gives
q_n(t*, C*_{w⊥}) = (1 − t*) q_{n−1}(t*, X_w), that is Δ_s(w) = −(1 − t*) q_{n−1}(t*, X_w) for every s and w. A version
of QC at t* therefore needs a direction with q_{n−1}(t*, ⟨w|_s C*|w⟩_s) = 0, which is a degenerate diagonal block. At
a first zero, every other compression strictly increases the form.

## 2. Averages that do not depend on C stop at t = 1/d

Thin-copy removal works by averaging q_{n−1} over rank-one contractions ⟨u|_s C|v⟩_s (which keep rank at most two)
with u, v inside the plane of C. The next proposition shows that such averages must depend on C to reach t = 1/2.

**Proposition 2.** Let μ be a probability measure on pairs (u, v) of unit vectors of C^d, invariant under
(u, v) ↦ (Uu, Uv) for every unitary U. For every operator C on n copies and every copy s,

    E_μ q_{n−1}(⟨u|_s C|v⟩_s) = a S + b q_{n−1}(Tr_s C),   S = Σ_{ij} q_{n−1}(X_ij),

with a > 0 and −1/d ≤ b/a ≤ 1.
- For u = v Haar distributed: (S + q(T))/(d(d+1)).
- For a Haar-random orthonormal pair: (dS − q(T))/((d+1)d(d−1)).

Since q_n(t, C) = S − t q_{n−1}(Tr_s C), a nonnegative combination of such averages equals q_n(t, ·) only when
t ≤ 1/d. For qutrits that is t ≤ 1/3, the PPT range.

*Proof.* Condition on c = |⟨u|v⟩|² and the phase φ of ⟨u|v⟩. Write v = √c e^{iφ} u + √(1−c) u′ with u′ a unit vector
orthogonal to u. Unitaries fixing u preserve μ and act transitively on the unit sphere of u^⊥, so given (u, c, φ) the
vector u′ is uniform there. In particular u′ and e^{iθ} u′ have the same law. Expanding q(⟨u|C|v⟩), the cross term
between ⟨u|C|u⟩ and ⟨u|C|u′⟩ is odd under u′ ↦ e^{iθ}u′ and averages to zero. So

    E_μ q = E[c] E_u q(⟨u|C|u⟩) + E[1 − c] E_{u ⊥ u′} q(⟨u|C|u′⟩),

where (u, u′) is a Haar-random orthonormal pair.
- The first average is (S + q(T))/(d(d+1)), by the fourth moments of a Haar vector.
- For the second, Σ_{i≠j} q(X_{e_i e_j}) = S − Σ_i q(X_{e_i}) for any orthonormal basis. Averaging over Haar bases
  gives (dS − q(T))/((d+1)d(d−1)).

Their ratios b/a are 1 and −1/d. ∎

At t = 1/2 the averaging in thin-copy removal runs over orthonormal pairs inside the plane of C, which is not
U(d)-invariant. Every argument that reaches t = 1/2 must use the operator in this way. The black-box impossibility
(activation of ρ by a PPT state) is a stronger form of the same fact.

## 3. The projector case as a monogamy relation, every n

Let V be a plane in (C^d)^{⊗n} with orthonormal basis v_0, v_1, and purify P_V/2 as ψ = (|0⟩_R v_0 + |1⟩_R v_1)/√2
with a qubit R. For a set J of copies put c_J = ‖ψ_{JR} − ψ_J ⊗ I/2‖_2². This is the correlation between R and the
copies in J, with c_∅ = 0 and c_{[n]} = 3/4.

**Proposition 3.**

    q_n(1/2, P_V) = (16/3) Σ_J [ (−1/2)^{|J|} + ½ (−1/2)^{n−|J|} ] c_J.

*Proof.* Put P(K) = Tr ψ_K² for sets K of copies.
- Since P_V = 2 Tr_R ψψ^*, ‖Tr_J P_V‖² = 4 P(J^c).
- Since ψ is pure, P(J^c) = Tr ψ_{JR}². Splitting ψ_{JR} into ψ_J ⊗ I/2 plus a part with zero partial trace over R
  (these are Hilbert-Schmidt orthogonal) gives P(J^c) = P(J)/2 + c_J.
- Let A = Σ_J (−1/2)^{|J|} P(J^c) and B = Σ_J (−1/2)^{|J|} P(J). The last identity gives A = B/2 + Σ_J (−1/2)^{|J|} c_J.
  Applied to complements, it gives B = A/2 + Σ_J (−1/2)^{n−|J|} c_J.
- Eliminating B gives (3/4)A = Σ_J [(−1/2)^{|J|} + ½(−1/2)^{n−|J|}] c_J, and q_n(P_V) = 4A. ∎

For small n:
- n = 2: q_2(P_V) = 4(3/4 − c_1 − c_2).
- n = 3: q_3(P_V) = 3/2 − 2(c_1 + c_2 + c_3). This is the monogamy form in `three-copy-notes/`.
- n = 4:

      q_4(P_V) = (16/3)[ 27/64 + (3/8) Σ_{|J|=2} c_J − (9/16) Σ_{|J|=1} c_J − (3/8) Σ_{|J|=3} c_J ].

At every n the endpoint for projectors says that a maximally mixed qubit cannot be correlated with the copies in the
alternating pattern above. `check_monogamy.py` checks the identity for random planes: n ≤ 4 at d = 3, and n ≤ 3 at
d = 4.

## 4. Translation-invariant MPS filters at large n

Take bond matrices A^x, B^x (x = 0, 1, 2) of size D, left vectors l_a, m_a (a = 0, 1) and right vectors r, r′:

    y_a = Σ_x (l_a^T A^{x_1} ⋯ A^{x_n} r) |x⟩,   z_a = Σ_x (m_a^T B^{x_1} ⋯ B^{x_n} r′) |x⟩,   C = Σ_a |y_a⟩⟨z_a|.

Define the one-copy transfer channel on operators on C^D ⊗ C^D,

    T(X) = Σ ⟨x x′|I − tF|u u′⟩ (A^u ⊗ B^{u′}) X (A^x ⊗ B^{x′})^*.

It is completely positive for t ≤ 1. In fact it is the local filtering of X ⊗ (I − tF) by v ⊗ |u⟩ ↦ A^u v on one
side and v ⊗ |u′⟩ ↦ B^{u′} v on the other.

**Proposition 4.** q_n(t, C) = ⟨ω| [T^n(|r r′⟩⟨r r′|)]^Γ |ω⟩, with ω = Σ_a conj(l_a) ⊗ m_a (Schmidt rank at most two)
and Γ the partial transpose on the second bond factor.

*Proof.* Expand ⟨y_a ⊗ z_b|W|y_b ⊗ z_a⟩ in the physical indices. Products of bond matrices collect into T^n, and the sum
over a, b with the crossed labels is ⟨ω|X^Γ|ω⟩. ∎

**Corollary.** Suppose T is primitive with Perron eigenvalue μ and right and left Perron eigen-operators X_*, Y_*
(both positive semidefinite). Then q_n/μ^n → ⟨Y_*, |rr′⟩⟨rr′|⟩ ⟨ω|X_*^Γ|ω⟩.
- If EP_n holds for every n, every such X_* is 2-undistillable as a state on C^D ⊗ C^D; for D = 2 it is PPT.
- Conversely, an NPT Perron operator at D = 2 would give counterexamples at every large n.

**Searches** (our logs; `check_mps.py` checks Proposition 4 and re-evaluates one search output).
- Translation-invariant MPS rank-two filters with D = 2, 3, 4. The right boundary vectors were also allowed to depend
  on a. The objective was a soft minimum of 2^n q_n/‖C‖² over n ∈ {3, …, 9, 12, 13, 16, 17, 24, 25, 32, 33}.
  - No negative value appeared.
  - Every candidate that looked negative in double precision is positive when re-evaluated with 60 digits
    (for example +2.4 × 10^{−11} at n = 33).
- At D = 2 we maximised the NPT-ness of the Perron operators X_*, Y_* over the bond matrices. The smallest
  partial-transpose eigenvalue was driven to about +10^{−12}, never below.

## 5. Activation of two-pair Werner mixtures

Let Π_S (S ∈ {s, a}²) be the products of the unnormalised symmetric and antisymmetric projectors on two qutrit pairs,
and σ_c = Σ_S c_S Π_S. Then ρ_W(1/2)^{⊗2} ∝ σ_c with c = (c_ss, c_sa, c_as, c_aa) = (1, 3, 3, 9). The first label
refers to pair 1.

**Proposition 5.** Let C_0 = |Φ⟩⟨Φ| ⊗ P. Here Φ is the maximally entangled state between pair 1 and a third, Werner
pair, and P is a rank-two projector on pair 2. Then the endpoint form of σ_c ⊗ ρ_W(1/2) satisfies

    Q(C_0)/‖C_0‖² = c_ss/2 − c_sa/6 + 3c_as/4 − c_aa/4.

*Proof.* The crossed-pairing form of a product operator is the product of the forms. The one-pair forms of P_sym and
P_anti are (‖X‖² ± |Tr X|²)/2, so:
- they take the values 3 and −1 on P;
- tensored with the Werner form, they take 1/3 and 1/2 on |Φ⟩⟨Φ|;
- ‖C_0‖² = 2. ∎

The expression vanishes at the Werner point, where C_0 is a zero of q_3 (an instance of P ⊗ |a⟩⟨b| with a = b = Φ
entangled across two copies). Wherever it is negative, a rank-two filter on three pairs distils σ_c ⊗ ρ_W.

On the slice c = (1, x, x, y) this means:
- **Activation** whenever y > 2 + 7x/3. The line passes through the Werner point (3, 9).
- **Positive margin is not enough.** At c = (1, 2, 2, 6.75), σ_c is one-copy 2-undistillable with margin 1/16 (local
  searches; `check_activation.py`), yet the witness value is −1/48. So two-pair Werner mixtures with a positive margin
  can be activated by one Werner copy. ρ_W^{⊗2} lies exactly on the boundary of this witness.
- **Outside the product cone, not activated** (numerical observations, not certificates). The one-copy
  2-undistillable part of the slice is x ≤ 3, y_lo(x) ≤ y ≤ 3 + 2x. Our searches found no activation by one or two
  further Werner copies at points outside the product cone K ⊗ K (y ≤ 3x, y ≥ 6x − 9) and below the witness line, for
  example (1, 1, 1, 4), (1, 2, 2, 0) and (1, 3, 3, 3). These are local searches, and they also missed activation at
  some points where the witness proves it, for example (1, 1, 1, 4.5). So this is weak evidence.
- **Conditional ratios are not enough.** "Every conditional anti/sym ratio is at most 3" does not imply one-copy
  2-undistillability: c = (1, 3, 3, 0) has two-pair minimum −0.42.

## 6. Sums of two block-product operators

**Proposition 6.** Partition the copies into blocks S_1, …, S_m and let

    C = Σ_{a=1,2} ⊗_j |α_a^{(j)}⟩⟨β_a^{(j)}|,

with arbitrary vectors α_a^{(j)}, β_a^{(j)} on block S_j. Let N^{(j)} be the 2 × 2 Gram matrix of q_{|S_j|}(1/2, ·) on
the two rank-one factors on block j. Then

    q_n(1/2, C) = Σ_{a,b} Π_j N^{(j)}_{ab}.

If EP holds for every block size that occurs, q_n(1/2, C) ≥ 0. In particular this holds for every n when all blocks
have at most two copies (by the two-copy theorem), and when all blocks have at most three copies (given the
three-copy theorem).

*Proof.* The form is a product over blocks, which gives the identity. Every element of span{|α_1⟩⟨β_1|,
|α_2⟩⟨β_2|} on a block has rank at most two, so N^{(j)} ⪰ 0 by EP for that block. The Hadamard product of positive
semidefinite matrices is positive semidefinite (Schur), so its entry sum is ≥ 0. ∎

When the vectors are entangled inside blocks of two copies, the reduced ranges are three-dimensional at every copy.
So this is another explicit class, thick at every copy, on which the endpoint inequality holds for every n with no
unproven input. Compare `contraction-notes/` (nested class) and:
- [T2], whose every-copy class has GHZ backbones |0⟩^{⊗k}, |1⟩^{⊗k} with tails at single copies and arbitrary
  co-range;
- [BGH, Thm. 8.3], which covers single tensor products, not sums.

The argument is tied to the two-term structure: a general rank-two operator expanded into many product terms gives
Gram matrices that are not positive semidefinite. `check_blocks.py` checks the identity and positivity for blocks
2 + 2 and 2 + 2 + 1.

## References

- [WZ] T. Wu, Q. Zou, Sharp Plücker geometry for three-copy Werner distillation, arXiv:2608.02647.
- [T2] T. S. Tomar, Equality cases and every-copy results for the Werner endpoint inequality, Zenodo 23718037 (2026).
- [BGH] K. Bharti, R. Gajjala, T. Haug, Two-copy nondistillability of Werner states: sharp partial-trace inequalities
  and finite-copy extensions, arXiv:2607.24479.
- [FGP] Fu, Gao, Park, A solution to 2-copy distillability of Werner states, arXiv:2607.21367.
- [SC] Z. Song, L. Chen, A partial-trace matrix inequality and Werner-state distillability, arXiv:2607.23416.

## Code

Python 3 with NumPy and SciPy; `check_mps.py` also uses mpmath. Outputs are in `logs/`.
- `code/werner.py`: shared helpers (as in `contraction-notes/`).
- `code/check_qc.py`: Section 1. It checks Proposition 1 and the trace identity, the two failing rules and the
  adversarial operators in `data/`.
- `code/check_averaging.py`: Section 2. It averages exactly over the four mutually unbiased bases of C³ (a 2-design).
- `code/check_monogamy.py`: Section 3.
- `code/check_mps.py`: Section 4. It checks Proposition 4 against dense q_n, shows a Perron operator, and re-evaluates a
  stored search output with 60 digits.
- `code/check_blocks.py`: Section 6.
- `code/check_activation.py`: Section 5. It checks the closed form and evaluates the two-pair minimum at points of the
  slice.
