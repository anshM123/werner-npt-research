# Compressions, black-box averages and large-n checks for the Werner endpoint

Ansh Mishra and Aryan Senthilkumar, October 2026.

These notes collect several small results about the n-copy endpoint inequality, with code:
- a conjecture that would settle every n, with the evidence for it;
- a proposition explaining why arguments that do not adapt to the operator stop at the PPT boundary;
- an exact form of the projector case as a monogamy relation;
- a transfer-channel description of translation-invariant MPS filters, with a large-n search;
- an exact activation witness for two-pair Werner mixtures;
- an every-n positivity statement for operators whose range lies in a product of planes over two-copy blocks;
- a lemma that a filter on one copy can lower the form in at most one direction, with a corollary about first zeros;
- an exact formula on the code spaces of pure qutrit codes, which shows the five-qutrit code space is far from the
  boundary;
- a reduction of the 2m-copy endpoint to a statement about m + 1 copies (a cone of bond operators that is invariant
  under the transfer channels), with numerical evidence for that statement at m = 1, 2. At m = 2 it would give the
  four-copy endpoint;
- exact positive spectra for the antisymmetric subspace of two copies and the qutrit [[3, 1, 2]] code. These give
  every-n, every-rank classes containing operators thick at every copy, with ranges of dimension up to 3^{n/2}.
- removal lemmas for blocks of two, three, four and five copies: the compressions to the antisymmetric subspace, the
  [[3, 1, 2]], [[4, 2, 2]] and [[5, 3, 2]] code spaces and the five-qutrit code space are separable. Iterating them gives
  the endpoint for every n whenever the range lies in a product of such blocks and planes, with at most three copies
  left free;
- an exact formula for the partial-transpose spectrum of the compression to any qutrit stabiliser code, as weight sums
  over cosets of the stabiliser group, and a separability test by stabiliser product states.

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
- Five copies (`check_qc_fast.py`). Fix a copy s and put Q[(a,b),(c,d)] = q_{n−1}(X_ab, X_cd) for the 3 × 3 blocks X_ab of
  C at s. Then q_n(C) and q_n(C_{w⊥}) are 9 × 9 expressions in Q and the projector (checked against dense q_n at n = 3,
  4), which makes the inner maximisation over w cheap. Four Danskin runs of 120 steps at n = 5 found no violation.
  The two random starts stay at +0.175 (q_5 ≈ 0.30). The two starts near a zero are driven to Δ ≈ 4 × 10^{−5} and
  2 × 10^{−5}, with q_5 ≈ 10^{−4} and 6 × 10^{−5}.
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

## 6. Ranges inside products of block planes

**Proposition 6.** Partition the copies into blocks S_1, …, S_m, each of at most two copies, and choose a plane V_j in
the space of block j (it may be entangled across the two copies). If the range of C lies in ⊗_j V_j, then
q_n(1/2, C) ≥ 0, whatever the rank and co-range of C. The same holds with co-range in place of range. With blocks of
at most three copies, it holds given the three-copy theorem.

*Proof.* Two facts combine:
- EP for k copies is equivalent to this: the one-sided compression (P_V ⊗ I) ρ^{⊗k} (P_V ⊗ I) to any plane V is PPT.
  Indeed, a vector has Schmidt rank at most two iff it lies in V ⊗ H for some plane V on the first side, and the
  partial transpose on the second side commutes with compressions on the first. So each block compression σ_j is
  PPT by the two-copy theorem.
- Tensor products of PPT operators are PPT, so ⊗_j σ_j, which is the compression of ρ^{⊗n} to ⊗_j V_j, is PPT.

Its partial transpose is the compression of B^{⊗n}, so ⟨ψ_C|B^{⊗n}|ψ_C⟩ ≥ 0 for every ψ_C in (⊗_j V_j) ⊗ H, that is
q_n(1/2, C) ≥ 0. ∎

**Special case (Schur form).** For C = Σ_{a=1,2} ⊗_j |α_a^{(j)}⟩⟨β_a^{(j)}| the range lies in
⊗_j span{α_1^{(j)}, α_2^{(j)}}, and directly q_n(1/2, C) = Σ_{a,b} Π_j N^{(j)}_{ab}. Here N^{(j)} is the positive
semidefinite 2 × 2 Gram matrix of q_{|S_j|} on the two rank-one factors, and the Schur product theorem applies.

When the planes are entangled inside two-copy blocks, the rank-two members of this class have three-dimensional
reduced ranges at every copy. So this is another explicit class, thick at every copy, on which the endpoint inequality
holds for every n with no unproven input. Compare `contraction-notes/` (nested class) and:
- [T2], whose every-copy class has GHZ backbones |0⟩^{⊗k}, |1⟩^{⊗k} with tails at single copies and arbitrary
  co-range;
- [BGH, Thm. 8.3], which covers single tensor products, not sums.

With single-copy blocks, Proposition 6 is the elementary fact that compressing a Werner copy to a plane on one side
gives a PPT state. The new content is that two-copy blocks with entangled planes can be used, thanks to the two-copy
theorem.

**Open.** If the block compression σ_j were separable (not just PPT), then a block where C is "thin" could be removed,
as single thin copies are. That would reduce n copies to n − 2 for operators thin at a two-copy block. A 2 × 9 PPT
state need not be separable, so this is open. Partial answers and evidence:
- **Planes thin at both copies.** If V lies in V_1 ⊗ V_2 for planes V_1, V_2 of the two copies, then σ_V is a local
  compression of σ_{V_1} ⊗ σ_{V_2}. Each factor is a PPT operator on C² ⊗ C³, hence separable, so σ_V is separable.
- **Planes in an antisymmetric subspace.** For planes inside a local-unitary image of Λ, see Section 10.
- **Numerical evidence (our logs).** For random planes, products alone fit σ_V to 10^{−8}, while a PPT entangled
  2 × 4 state (Horodecki's) embedded in 2 × 9 stays at 10^{−2}. A local search that minimises the smallest
  partial-transpose eigenvalue of σ_V drives it to 0. Every plane it found on that boundary was thin at both copies,
  so it is separable by the first point. The planes thick at a copy that we sampled had a positive margin.

`check_blocks.py` checks the Schur identity, and positivity for ranges in products of block planes with arbitrary
co-range (ranks 4 and 8).

## 7. At most one negative filter direction

For an operator C, a copy s and G ∈ M_3, write G_s C for G applied to copy s from the left.

**Proposition 7.** Let n ≥ 2 and 0 ≤ t ≤ 1/2, and assume EP_{n−1}(t). For every C of rank at most two and every copy
s, the Hermitian form Q_s(G) = q_n(t, G_s C) on M_3 has at most one negative eigenvalue. The same holds for right
filters, G ↦ q_n(t, C G_s).

*Proof.* If det G = 0, the range of G_s C at copy s lies in a plane, so q_n(t, G_s C) ≥ 0 by [T2, Lemma 6.1]. If Q_s
had two negative eigenvalues, it would be negative definite on a two-dimensional subspace L of M_3. On L, det is a
binary cubic form, so it has a nonzero zero over C. That zero is a singular G ≠ 0 with Q_s(G) < 0, a
contradiction. ∎

So if EP_{n−1} holds and q_n(C) < 0, then at every copy the filters G with q_n(t, G_s C) < 0 form, up to scalars, a
complex-hyperbolic ball in PGL(3). This ball contains I and avoids every singular matrix.

**Corollary (first zeros).** Suppose EP_{n−1} holds at t = 1/2 and EP_n fails there. Let t* < 1/2 and C* be the first
zero of [T2, Sec. 6]:
- q_n(t*, ·) ≥ 0 on operators of rank at most two;
- q_n(t*, C*) = 0;
- the local supports of the range and co-range of C* are three-dimensional at every copy.

Then for every copy s, Q_s is positive semidefinite with kernel exactly C·I, for left and for right filters.

*Proof.* G_s C* has rank at most two, so Q_s ≥ 0, and Q_s(I) = 0.
- Write the weight −t at each copy as (1 − 2t)·0 + 2t·(−1/2). This expresses q_{n−1}(t, X) as a positive combination
  of terms, one for each set S of copies. The term for S sums q_{|S|}(1/2, ·) over the contractions of X on the other
  copies, so it is nonnegative by EP_{|S|}(1/2), which follows from EP_{n−1}(1/2). The term for S = ∅ is
  (1 − 2t)^{n−1}‖X‖². So q_{n−1}(t*, X) ≥ (1 − 2t*)^{n−1}‖X‖² for X of rank at most two.
- If the kernel had dimension at least two, it would contain a nonzero singular G. Then G_s C* ≠ 0, because the range of C* at
  copy s is all of C³.
- [T2, Lemma 6.1] with this margin gives q_n(t*, G_s C*) ≥ (1 − 2t*)^n ‖G_s C*‖² > 0, a contradiction. ∎

**The bound is attained.** In the activation example of Section 5 at c = (1, 2, 2, 6.75), Q(C_0) < 0. There the
filter form at the Werner copy has exactly one negative eigenvalue (`check_signature.py`). The other factor σ_c is
one-copy 2-undistillable by the local searches of Section 5, which is what the proof needs.

The proposition does not rule out a counterexample. It says that a failure, if there is one, is one-dimensional in
the filter at each copy, and that a first zero is nondegenerate in the filter directions other than I.

## 8. Code spaces

Local searches are drawn towards the thin zeros, so we test a structured, highly entangled family directly:
operators supported on the code space of a qutrit error-correcting code. Let V : C^K → (C³)^{⊗n} be the encoding
isometry of a code that is pure up to distance d: for every code state, the reduced state on any set of fewer than d
copies is maximally mixed.

**Proposition 8.** For x ∈ M_K and every set L of copies:
- if |L^c| < d, then ‖Tr_L(V x V^*)‖² = 3^{−|L^c|} |Tr x|²;
- if |L| < d, then ‖Tr_L(V x V^*)‖² = 3^{−|L|} ‖x‖².

So when n = 2d − 1 (every set or its complement has fewer than d copies),

    q_n(1/2, V x V^*) = A_n ‖x‖² + B_n |Tr x|²,   A_n = Σ_{l<d} C(n,l) (−1/6)^l,   B_n = Σ_{l≥d} C(n,l) (−1/2)^l 3^{l−n}.

*Proof.* Write v_i = V e_i.
- If |L^c| < d, then Tr_L(V y y^* V^*) = ‖y‖² I/3^{|L^c|} for every y. By polarisation, Tr_L(V x V^*) = Tr(x) I/3^{|L^c|}.
- If |L| < d, the swap trick gives ‖Tr_L(V x V^*)‖² = Σ x̄_{ij} x_{kl} Tr[Tr_{L^c}(v_j v_l^*) Tr_{L^c}(v_k v_i^*)]. By the
  first case applied to L, Tr_{L^c}(v_j v_l^*) = δ_{jl} I/3^{|L|}. So the sum is 3^{−|L|} Σ |x_{ij}|². ∎

Consequences.
- **The five-qutrit code.** The stabilisers are the cyclic shifts of X Z Z^† X^† I. This code has dimension 3, and no
  nontrivial stabiliser element has weight below 3. Moreover Tr_{S^c}(v_i v_j^*) = δ_{ij} I/3^{|S|} for every set S
  of at most two copies and every pair of basis vectors (to 2 × 10^{−16}), so it is a pure [[5, 1, 3]]_3 code. Here
  A_5 = 4/9 and B_5 = −19/288.
  - Since |Tr x|² ≤ 3‖x‖², q_5(1/2, C) ≥ (A_5 + 3B_5)‖C‖² = (71/288)‖C‖² for every operator C on the code space, of
    any rank.
  - Equivalently, the two-sided compression of ρ^{⊗5} to the code space is PPT, with a margin.
  - **One-sided version (computed).** The partial transpose of the one-sided compression (V^* ⊗ I) ρ^{⊗5} (V ⊗ I) has
    spectrum in [71/288, 4/9] to 10 digits (`check_codes.py`, part (c)). So q_5(1/2, C) ≥ (71/288)‖C‖² for every C
    whose range lies in the code space, whatever its co-range and rank.
  - These operators are thick at every copy.
- **Three-dimensional blocks.** Section 6 needs only that each block compression be PPT. So Proposition 6 extends to
  partitions into blocks of at most three copies, with a plane in each, and blocks of five copies, with the code space
  (or an image of it under local unitaries) in each. If the range of C lies in the product of these block subspaces,
  then q_n(1/2, C) ≥ 0 for every rank. Unlike Section 6, the block factors here have dimension 3.
  - At one copy, no 3-dimensional subspace works.
  - At five copies, a random 3-dimensional subspace also gave a positive smallest eigenvalue (0.243 in one sample), so
    generic small subspaces behave like the code.
- **Larger codes.** For a pure [[2d − 1, 1, d]]_3 code, A_n + 3 min(B_n, 0) > 0 for every odd n from 3 to 21
  (`check_codes.py`), and the effective parameter −B_n/A_n tends to 0. We have only checked the five-qutrit code. Codes
  with these parameters may not exist for larger n, so the table describes what such a code would give, not a theorem
  for every n.

So operators that are highly entangled in the code sense sit far from the boundary. Every zero we know is thin.

## 9. Invariant cones: from 2m copies to m + 1

Proposition 4 does not use translation invariance. Write the two families of C = Σ_a |y_a⟩⟨z_a| as matrix product
states with site-dependent tensors, in a common copy order:

    y_a = Σ_x (l_a^T A_1^{x_1} ⋯ A_n^{x_n} r) |x⟩,   z_a = Σ_x (m_a^T B_1^{x_1} ⋯ B_n^{x_n} r′) |x⟩.

Here A_s^x is D_{s−1} × D_s and B_s^x is D′_{s−1} × D′_s. Successive singular value decompositions give such a form
with D_s equal to the rank of y as a matrix with rows indexed by (a, x_1, …, x_s), so D_s ≤ min(2·3^s, 3^{n−s}). Let
T_s be the transfer channel of Section 4 built from A_s and B_s, and put X_n = |r r′⟩⟨r r′| and X_{s−1} = T_s(X_s).
The proof of Proposition 4 gives

    q_n(1/2, C) = ⟨ω| X_0^Γ |ω⟩,   ω = Σ_a conj(l_a) ⊗ m_a.

Each T_s is a local filtering of X ⊗ ρ, where ρ = I − F/2 is one Werner copy. So the endpoint inequality follows
whenever the bond operators X_s stay in a cone of operators that are 2-undistillable. We use the cones generated by a
fixed number of copies.

**Definition.** For m ≥ 1, let K_m(D, D′) be the closed convex cone generated by the local filterings
(a ⊗ b) ρ^{⊗m} (a ⊗ b)^*, with a : (C³)^{⊗m} → C^D and b : (C³)^{⊗m} → C^{D′}. Rank-one filters give every product
operator, so K_m(D, D′) contains the separable cone. Contracting a copy with a product vector shows K_k ⊆ K_m for
k ≤ m. **L_m(D)** is the statement that every local filtering of ρ^{⊗(m+1)} with output dimensions at most D lies in
K_m, and **KL_m** = L_m(3^m). L_m(D) implies L_m(D′) for D′ ≤ D.

**Proposition 9.**
- (a) Assume EP_m and KL_m. Let C have rank at most two on n copies, for any n, with MPS forms as above whose bond
  dimensions are all at most 3^m. Then q_n(1/2, C) ≥ 0.
- (b) EP_m and L_m(2·3^{m−1}) imply EP_{2m}. Since EP_2 holds [BGH, FGP, SC], L_2(6) implies EP_4: every 6 × 6 local
  filtering of three Werner copies would have to lie in the cone generated by filterings of two.

*Proof.* We show X_s ∈ K_m(D_s, D′_s) by downward induction on s. For s = n this holds because X_n is a nonnegative
scalar. T_s is linear and continuous, so it is enough to check generators.
- A product operator |e⟩⟨e| ⊗ |f⟩⟨f| is sent to a local filtering of one copy of ρ. That lies in K_m, because ρ is a
  positive multiple of a filtering of ρ^{⊗m}: contract the other m − 1 copies with |0…0⟩ on both sides, which
  gives the factor (1/2)^{m−1}.
- A filtering (a ⊗ b) ρ^{⊗m} (a ⊗ b)^* is sent to a local filtering of ρ^{⊗(m+1)} with output dimensions D_{s−1},
  D′_{s−1} ≤ 3^m. That lies in K_m by KL_m.

Now let X = (a ⊗ b) ρ^{⊗m} (a ⊗ b)^* be a generator of K_m(D_0, D′_0). Then
⟨ω|X^Γ|ω⟩ = ⟨ω′|(ρ^{⊗m})^Γ|ω′⟩ with ω′ = (a^* ⊗ b^T) ω. This vector has Schmidt rank at most two, so the value is
nonnegative by EP_m. Separable generators also give nonnegative values, so ⟨ω|X_0^Γ|ω⟩ ≥ 0. This proves (a).

For (b), take n = 2m. Then D_s ≤ min(2·3^s, 3^{2m−s}).
- For s ≥ m, X_s is a filtering of 2m − s ≤ m copies, so it lies in K_m with no hypothesis.
- For s ≤ m − 1, the bond dimensions are at most 2·3^{m−1}. So the induction step above needs only L_m(2·3^{m−1}).
The last step is the same as in (a). ∎

Remarks.
- **Dual form.** W is in the dual cone of K_m(D, D′) iff W ⊗ ρ^{⊗m} is block-positive for the cut
  (C^D ⊗ (C³)^{⊗m}) | (C^{D′} ⊗ (C³)^{⊗m}). So KL_m says this: for operators W on C^{3^m} ⊗ C^{3^m},
  block-positivity of W ⊗ ρ^{⊗m} survives one more Werner copy.
- **Strength.** KL_m concerns m + 1 copies, and together with EP_m it implies EP_{m+1} (take output dimension two).
  The first case, KL_1, says that every 3 × 3 local filtering of ρ^{⊗2} lies in the cone generated by local
  filterings of ρ. Together with EP_1 (Cauchy–Schwarz), KL_1 would give q_n ≥ 0 for every n on all operators whose
  bond dimension is at most 3 at every cut, with arbitrary site-dependent tensors. That class contains operators that
  are thick at every copy.
- **The dimension bound is needed.** If ρ ⊗ ρ itself (output dimension 9) were in K_1, induction would put every
  ρ^{⊗n} in the cone generated by filterings of one copy, and EP_1 would give EP_n for all n. Numerically it is not:
  see below.
- **No induction in m.** A filter from 3^{m+1} to 3^m dimensions need not factor through 3^{m−1} dimensions at any
  cut, so KL_{m−1} does not give KL_m.

**Removal by compression.** The following identity is a common form of Proposition 9(b), of Propositions 11 and 12
in Section 10, and of thin-copy removal. For a set T of b copies, let S_T and S′_T be the range and co-range of C at
T: the spans of (I_T ⊗ ⟨r|) y_a and (I_T ⊗ ⟨r|) z_a, where C = Σ_a |y_a⟩⟨z_a| and r runs over a basis of the other
copies. Then

    q_n(1/2, C) = ⟨ψ_C| (τ_T ⊗ ρ^{⊗(n−b)})^Γ |ψ_C⟩,   τ_T = (P_{S_T} ⊗ P_{S′_T}) ρ^{⊗b} (P_{S_T} ⊗ P_{S′_T}).

Here ψ_C = Σ_a y_a ⊗ conj(z_a). Its second side lies in conj(S′_T), and moving that projector through the partial
transpose turns it into P_{S′_T}; compressing to conj(S′_T) instead gives a different, wrong value
(`check_removal.py`). As in Proposition 11, if τ_T lies in K_k (with K_0 the separable cone), then EP_{n−b+k}
implies q_n(1/2, C) ≥ 0. The dimensions of S_T and S′_T are at most 2·3^{n−b}. For n = 4 and b = 3, τ_T is a 6 × 6
local filtering of three copies, and k = 2 is L_2(6) for that one filtering. So the four-copy endpoint needs, for each
C, only one triple for which this holds.

Numerical picture (`check_removal.py` and our logs):
- For random C at n = 4, every triple gives a τ_T that is PPT with margin: the smallest normalised partial-transpose
  eigenvalue is about 0.5. Products alone fit τ_T to 10^{−8}, so k = 0 applies and EP_1 suffices for such C.
- At n = 5, random C give margins of about 0.8 for the best block of four copies.
- At n = 3 the analogous pair compressions are NPT, so this route does not reduce EP_3 to EP_1.
- Near the zeros of q_4 the triple compressions sit on the PPT boundary. In one near-minimiser, with
  q_4/‖C‖² = 3 × 10^{−6}, the triple that leaves out a nearly thin copy is separable (margin 0.2). In another, all four
  triples are slightly NPT (down to −0.025). There the fits by filterings of one or two copies are no worse than
  5 × 10^{−4}, the error floor for such targets (see the calibration below), so the test is inconclusive.
- For C = |y_1⟩⟨y_1| + |y_2⟩⟨y_2| with y_1 = Φ_12 ⊗ Φ_34, maximally entangled across copies on each side, all four
  triple compressions contain a full Werner copy and are NPT (−0.75). Products fit them only to 0.11, and filterings of
  one or two copies to about 10^{−3}. So k = 0 is not always available.

**Every n from the cone statements (proved reduction).** Take T to be all copies but one. Then S_T and S′_T have
dimension at most 6, so τ_T is a 6 × 6 local filtering of ρ^{⊗(n−1)}, and the identity with k = n − 2 gives:
- EP_{n−1} and L_{n−2}(6) imply EP_n, for every n ≥ 3.
- For m ≥ 2 we have 6 ≤ 3^m, and then L_m(6) follows from KL_m. Embed C^6 in C^{3^m}, apply KL_m to the embedded
  filtering, and compress each term of the decomposition back to C^6: filterings stay filterings and products stay
  products.
- EP_3 is known (Proposition 12(c) uses it). So the endpoint for every n would follow from KL_m for every m ≥ 2, one
  cone statement per number of copies. KL_1 is not needed for this chain.
- The chain cannot start at n = 3, where it would need L_1(6). The fits below suggest that L_1(6) is false. That is
  consistent with the dimension bound 3^m in KL_m, which is then sharp at m = 1, and it does not bear on L_2(6),
  since 6 ≤ 9.

**Evidence** (`check_cones.py`; numerical, not proofs). Every term of the fitting model lies in the cone, so a small
residual is an explicit approximate decomposition. Residuals are relative Frobenius norms.
- Part (a) checks the transfer identity against the dense q_n for n ≤ 5. The bond dimensions are at most 3, 6, 9 and
  18 for n = 2, 3, 4, 5.
- **KL_1** (parts (b), (e)).
  - Random 3 × 3 filterings of ρ^{⊗2} are PPT and fit by products alone (to 10^{−8}), so they are separable and test
    nothing. The informative targets are NPT. Part (b) drives λ_min(Y^Γ)/Tr Y down along a path of two-copy filterings
    and stops at −0.01 and −0.03. Neither products alone (0.02–0.08) nor a single filtering of ρ (0.10–0.14) fits
    these targets, so they are neither separable nor filterings of one copy. Six filterings of ρ plus 12 products fit
    them to 2–5 × 10^{−7} (1–2 × 10^{−8} with restarts, in our logs). The antisymmetric projector (1-distillable)
    stays at 0.43.
  - In our logs, an adversarial search that moves the filter to maximise the fit residual never got past
    1.1 × 10^{−4}, and its worst points refit to 4 × 10^{−8}.
  - A possible handle for a proof (our logs). NPT targets at −0.01 fit as a single filtering of ρ plus 12 products,
    to 1–5 × 10^{−8}. The filtering carries 44–57% of the trace. The remainder is full rank and PPT, with its smallest
    partial-transpose eigenvalue at 10^{−4}–10^{−3} of its trace, so it sits on the PPT boundary.
    - The filtering is not a contraction of the target's own filters: a a^* is 0.16–0.33 away from every A(I ⊗ |u⟩)
      and A(|u⟩ ⊗ I) (normalised Frobenius distance).
    - A proof of KL_1 along these lines would need a rule for (a, b), and a reason why the remainder is separable
      rather than only PPT.
  - Every element of K_1 has Schmidt number at most 2, so it satisfies Y ≤ 2 Tr_B(Y) ⊗ I. Over two-copy filterings,
    the smallest eigenvalue of 2 Tr_B(Y) ⊗ I − Y was driven to 0 (−3 × 10^{−16}), never below.
  - Filters of rank at most two need no test: the output lives on C² ⊗ C³, is PPT by EP_2, and so is separable. For
    full-rank filters the output is positive definite, so the range criterion says nothing. Positive definiteness and
    2-undistillability alone do not imply membership: the Tiles state plus a small multiple of the identity is
    outside K. So KL_1 has to use the two-copy structure.
- **KL_2** (parts (c), (f)). The KL_2 model is 96 filterings of ρ^{⊗2} plus 240 products. In our logs, a random
  element of the cone fits to about 10^{−6}.
  - Two states outside the cone stay far away. The first is the 2-distillable P_a ⊗ P (0.43). The second is the PPT
    entangled Tiles state tensored with |00⟩⟨00| (0.07). The Tiles state is provably outside: the generators of K_m
    form a compact set avoiding 0, so K_m consists of finite sums of generators, and each generator has range
    range(a) ⊗ range(b). So the range of every element of K_m is spanned by product vectors, which fails for Tiles.
    The fit therefore also separates the cone from PPT entangled states.
  - Random 9 × 9 filterings of ρ^{⊗3} are PPT, close to the boundary, and 700 products fit them to 10^{−7}–10^{−6}. So
    they are separable for practical purposes and test nothing.
  - NPT filterings, taken along a path that lowers λ_min(Y^Γ)/Tr Y and stopped at −0.003 and −0.016, are the real
    test. Products alone (0.02–0.10) and a single filtering of ρ^{⊗2} (0.35–0.43) do not fit them. The KL_2 model
    fits them only to 5–7 × 10^{−5}, well above in-cone controls. In our logs, column generation repeatedly adds the
    filterings of ρ^{⊗2} that are most negative against the current residual, then refits. On one target it brought
    the residual from 10^{−4} to 4.8 × 10^{−6}, still falling, and the oracle still found descent directions.
  - Calibration (part (f)): a single generator of K_2, which is in the cone exactly, fits a large model started at
    random only to 6 × 10^{−4}, while a one-term model fits it to 4 × 10^{−7}. So residuals of this size on near-extreme
    targets reflect the fitting method, not a distance to the cone.
  - With smaller term budgets all residuals were larger, controls included.
  - In our logs, two other tests were run:
    - Filterings that minimise λ_min(Y^Γ)/Tr Y converge to −1/15, the smallest value we found on generators of K_2.
      They turn out to have rank 9, with flat rank-3 reductions: they are filterings of a single copy, so they lie
      in K_2 trivially. So this test does not probe KL_2.
    - The 2-reduction test above, applied at 9 × 9, was not violated by two-copy or three-copy filterings (minimum 0
      to six digits).
  - We read this as no detected gap, but it is the weakest part of the evidence. L_2(6), which is what EP_4 needs, is
    implied by KL_2: a 6 × 6 filtering is a compression of a 9 × 9 one, and K_2 is closed under filtering.
  - A primal adversary for KL_2 (`check_lmD_adversary.py` with m = 2, D = 9) drives a model of 16 filterings and 80
    products to 2.8 × 10^{−2}. That is the size of the model, not a gap (`check_kl2_refit.py`):

    | Model size | Element of K_2(9) | Adversarial target |
    |---|---|---|
    | 16 filterings + 80 products | 1.1 × 10^{−3} | 1.6 × 10^{−2} |
    | 32 + 200 | 5.0 × 10^{−4} | 4.7 × 10^{−4} |
    | 48 + 400 | 3.2 × 10^{−4} | 1.3 × 10^{−4} |

    So for two copies the dimension 3^2 = 9 shows no gap.
- **Uncompressed** (part (d)). ρ ⊗ ρ as a 9 × 9 operator fits K_1(9, 9) only to 0.12, while a random element of
  K_1(9, 9) fits to 1–2 × 10^{−4} at the same budget (10^{−6} with more terms, in our logs). So the dimension bound in
  KL_m is needed.
- **L_1(6)** (`check_l1_six.py` and our logs). Output dimension 6 already seems to be too much for one copy.
  - With 16 filterings of ρ and 80 products, an element of K_1(6) built from 3 filterings and 10 products fits only to
    4 × 10^{−4}. This floor comes from the optimisation: it does not improve with more terms or with column generation.
  - A Schmidt-rank-4 state fits to 0.61.
  - Random 6 × 6 filterings of ρ^{⊗2} (partial-transpose minimum −0.001 to −0.002) fit to 1.2–1.7 × 10^{−3}.
    Filterings pushed to −0.01 fit to 4–7 × 10^{−3}, and those pushed to −0.03 to 0.9–1.6 × 10^{−2}.
  - Doubling the model leaves these 3–10 times above the floor. Six rounds of column generation move the −0.01 target
    only from 5.01 to 5.00 × 10^{−3}.
  - The pair compressions τ_T of random three-copy operators, which are what removal by compression would need at
    n = 3, fit to 0.7–1.2 × 10^{−2}.
  - The antisymmetric witness of [FGP, Thm. 2.4] does not detect Schmidt number 3 in these filterings (its maximum is
    0.3 < 1/2).

  - A primal adversary settles the comparison (`check_lm6_adversary.py`). It takes ascent steps on the filters, holding
    the fit fixed (envelope gradient), then refits. Within 60 steps it reaches residual 6 × 10^{−2}, at
    partial-transpose minimum −0.015. Its worst target refits from fresh starts to 6.1 × 10^{−2}, and to 5.9 × 10^{−2}
    with twice the terms: 150 times the floor.

  - The same adversary at other output dimensions (`check_lmD_adversary.py`) finds no gap for D = 4: the worst target
    refits to 2.5 × 10^{−4}, the floor for near-extreme targets. It finds a robust gap for D = 5: 1.3 × 10^{−2}, the
    same with twice the terms. So for one copy, L_1(D) holds for D ≤ 4 and fails for D ≥ 5.

  We read this as L_1(6) being false. The threshold for one copy lies between 4 and 5, slightly above 3^1. This is
  numerical, not a proof, and the floor depends on the target. A genuinely three-copy filtering below
  has partial-transpose minimum −0.022 and still fits K_1(6) to 6 × 10^{−5}. So only the strongly NPT filterings of two
  copies resist.
- **L_2(6)** (our logs). Random 6 × 6 filterings of ρ^{⊗3} are PPT (+0.005), and a gradient path started at random
  stalls at the PPT boundary. So we built three-copy targets from the two-copy filtering above at −0.030: tensor with a
  product vector on the third copy, then add δ times a random three-copy filter.
  - The two-copy cone K_2(6) fits its own random elements to 7 × 10^{−8}, and a single generator (δ = 0) to
    3.5 × 10^{−4}. That floor comes from random starts on an extreme ray.
  - δ = 0.1 (minimum −0.029) fits K_2(6) to 2.9 × 10^{−4}, and δ = 0.3 (−0.022) fits to 6 × 10^{−5}. For comparison,
    K_1(6) gives 9.3 × 10^{−3}, 6.7 × 10^{−3} and 6 × 10^{−5}.
  - The same adversary for L_2(6) finds nothing beyond the floor (`check_lm6_adversary.py`,
    `check_l2six_adversary_npt.py`).
    - From a random start it stays at 10^{−5}–10^{−4}, and its worst target refits to 4 × 10^{−8}.
    - Started from the δ = 0.3 target, it climbs to 1.4 × 10^{−3} while it tracks. But the worst target refits from
      fresh starts to 3.0 × 10^{−4} (2.9 × 10^{−4} with twice the terms), which is the floor for near-extreme NPT
      targets.
  - So nothing here argues against L_2(6): the two-copy cone absorbs what the one-copy cone cannot.

**Open.** KL_1, a two-copy statement about 3 × 3 operators, is the smallest open case. A proof would give an every-n
class of bond dimension 3. For four copies, L_2(6) is enough: a three-copy statement about 6 × 6 operators. For every
n, L_m(6) for all m ≥ 2 is enough (each follows from KL_m): one statement about 6 × 6 operators per number of copies.

## 10. Antisymmetric and code blocks

Section 6 needs only that the one-sided compression of each block be PPT, as Section 8 already used. Two small
blocks have this property with a margin, and both are thick at every copy.

**Proposition 10.**
- (a) Let Λ ⊂ C³ ⊗ C³ be the antisymmetric subspace of two copies on one side, spanned by |ij⟩ − |ji⟩ (dimension 3).
  The partial transpose of the one-sided compression (P_Λ ⊗ I) ρ^{⊗2} (P_Λ ⊗ I) has spectrum 1/4 (×8), 3/4 (×9) and
  1 (×10). So q_2(1/2, C) ≥ ¼‖C‖² for every operator C whose range lies in Λ, whatever its rank and co-range.
- (b) Let Q ⊂ (C³)^{⊗3} be the code space of the qutrit [[3, 1, 2]] code, spanned by Σ_k |k, k + j, k + 2j⟩ for
  j = 0, 1, 2 (indices mod 3). The partial transpose of the one-sided compression of ρ^{⊗3} to Q has spectrum
  1/2 (×56), 3/4 (×24) and 7/8 (×1). So q_3(1/2, C) ≥ ½‖C‖² when the range of C lies in Q.

*Proof.* Both spectra are exact rational computations (`check_blocks2.py`, part (a)). As in Section 6, the partial
transpose on the second side commutes with compressions on the first. So for a unit vector ψ in S ⊗ H, which is
ψ_C for an operator C with range in S, ⟨ψ|(ρ^{⊗b})^Γ|ψ⟩ is at least the smallest eigenvalue. ∎

The multiplicities in (a) are the dimensions 10, 8 + 1 and 8 of the pieces of C³ ⊗ C³ ⊗ C³ under U(3), which is
consistent with the symmetry: ρ^{⊗2} commutes with U^{⊗4}, and Λ is a U(3)-irreducible subspace.

**Consequence (every n, every rank).** Partition the copies into blocks. In each block choose one of the following:
- a plane (blocks of at most three copies);
- a local-unitary image of Λ (two copies);
- a local-unitary image of Q (three copies);
- a local-unitary image of the five-qutrit code space (five copies);
- a local-unitary image of the qutrit [[4, 2, 2]] code space (four copies, dimension 9; stabilisers X^{⊗4} and
  Z ⊗ Z² ⊗ Z ⊗ Z²), whose compression has partial-transpose spectrum in [1/3, 43/48] (`check_blocks2.py`, part (g))
  and is separable (Proposition 14).

If the range of C lies in the tensor product of the chosen subspaces, then q_n(1/2, C) ≥ 0, whatever the rank and the
co-range. The proof is that of Proposition 6.
- Local unitaries U_1 ⊗ ⋯ ⊗ U_b applied on both sides preserve ρ^{⊗b}, so each block compression is unitarily
  equivalent to that of the chosen subspace.
- Tensor products of PPT operators are PPT.

If there are no plane blocks, the margins multiply. For example, q_n(1/2, C) ≥ 4^{−n/2}‖C‖² when n is even and every
block is Λ. Generic operators in this class are thick at every copy. With Λ blocks the range has dimension 3^{n/2}.

Products of planes on single copies have larger ranges, of dimension 2^n, but they are thin at every copy and follow
from thin-copy removal. Among the classes in these notes that are thick at every copy, the Λ blocks give the largest
ranges.

**Larger blocks** (`check_blocks2.py`, part (b), and our logs; numerical).
- At two copies, Λ looks maximal. Adding |00⟩, Σ_i|ii⟩/√3 or (|01⟩ + |10⟩)/√2 to Λ makes the compression NPT. So do
  Sym²(C³) and random 4-dimensional subspaces; the best 4-dimensional subspace found had smallest eigenvalue −10^{−4}.
- At three copies, Q ⊕ C·ε is still PPT (smallest eigenvalue 0.041). Here ε is the totally antisymmetric vector.
- At three copies, 200 random 3-dimensional subspaces were all PPT (smallest 3.6 × 10^{−4}), and 95% of random
  4-dimensional ones. Sym³(C³) and the mixed-symmetry subspaces are NPT.
- Not every 3-dimensional subspace works: C³ ⊗ |00⟩ gives ρ ⊗ (positive), which is NPT.
- For the [[4, 2, 2]] code, the compression is block diagonal over its 9 syndrome sectors, but the 9 × 9 blocks have
  purity 0.0145, above the Gurvits–Barnum bound 1/80. Proposition 14 shows separability by another route.
- The codes [[b, b − 2, 2]] give PPT blocks only for b ≤ 6 (Proposition 14(d)).

**Removing an antisymmetric pair.** The compression to Λ is separable, not just PPT. This answers, for Λ, the
question left open at the end of Section 6, and it gives a removal lemma for a pair of copies. It works like
thin-copy removal [WZ, T2], but here C can be thick at every copy.

**Proposition 11.** Let C have rank at most two on n ≥ 3 copies. Suppose that for two copies s ≠ s′ and unitaries
U, U′, the range of C lies in (U ⊗ U′)Λ at copies s, s′, tensored with the full space of the other copies. If
EP_{n−2} holds, then q_n(1/2, C) ≥ 0. The same holds with the co-range in place of the range. Since EP_2 and EP_3
hold, the four- and five-copy endpoint inequalities hold for every such C.

*Proof.*
- **The compression.** Let a_z = Σ_{xy} ε_{zxy}|xy⟩/√2 (z = 0, 1, 2) be the basis of Λ. Expand
  ρ^{⊗2} = Σ_J (−1/2)^{|J|} F_J and use the determinant formula for ε_{zab} ε_{z′cd}. This gives
      σ_Λ = (P_Λ ⊗ I) ρ^{⊗2} (P_Λ ⊗ I) = ¾ I_Λ ⊗ P_anti + ½ (I_Λ ⊗ P_sym)(I + D)(I_Λ ⊗ P_sym).
  Here P_sym and P_anti are the symmetric and antisymmetric projectors of the second side's two copies. The operator
  D = |Φ′⟩⟨Φ′| ⊗ I has Φ′ = Σ_z |z⟩|z⟩ between Λ (in the basis a_z) and the second side's first copy.
- **Separability.** I + |Φ′⟩⟨Φ′| = Σ_{k=1}^{12} |u_k⟩⟨u_k| ⊗ |ū_k⟩⟨ū_k|, summed over the 12 vectors of four mutually
  unbiased bases of C³. This is the partial transpose of Σ_k |u_k⟩⟨u_k|^{⊗2} = I + F, which holds because these
  vectors form a 2-design. The term I_Λ ⊗ P_anti is a product, and the projector I_Λ ⊗ P_sym acts on the second side
  only. So σ_Λ is separable across Λ | (second side), with this explicit decomposition. Local unitaries preserve
  ρ^{⊗2}, so the same holds for (U ⊗ U′)Λ.
- **Removal.** Write the compression as Σ_k |e_k⟩⟨e_k| ⊗ |f_k⟩⟨f_k|. As in Section 6, the partial transpose commutes
  with the compression on the first side, and ψ_C lies in (U ⊗ U′)Λ ⊗ (rest). So
      q_n(C) = ⟨ψ_C|(σ ⊗ ρ^{⊗(n−2)})^Γ|ψ_C⟩ = Σ_k ⟨φ_k|(ρ^{⊗(n−2)})^Γ|φ_k⟩,   φ_k = (⟨e_k| ⊗ ⟨f̄_k|) ψ_C.
  Each φ_k has Schmidt rank at most two, so every term is nonnegative by EP_{n−2}. ∎

Remarks.
- A plane inside some (U ⊗ U′)Λ also has a separable compression, since compressions of separable operators are
  separable. So for these planes, the open question of Section 6 has a positive answer.
**Proposition 12 (removing a code triple, and iterating).**
- (a) The one-sided compression σ_Q of ρ^{⊗3} to the code space Q is separable. On the second side, the 27-dimensional
  space splits into nine syndrome sectors (X^a Z^b ⊗ I ⊗ I)Q, with a, b ∈ Z_3 and X, Z the qutrit shift and clock.
  σ_Q is block diagonal in these sectors. Take the basis |j⟩ (logical) ⊗ (X^a Z^b ⊗ I ⊗ I)|k_L⟩ (sector).
  - In the code sector the block is ½I + ⅛F, with F the swap. Since I + F = Σ_{k=1}^{12} |u_k u_k⟩⟨u_k u_k| over the
    12 vectors of four mutually unbiased bases (a 2-design), this equals ⅜I + ⅛ Σ_k |u_k u_k⟩⟨u_k u_k|.
  - In each of the other eight sectors the block is ½I + ¼Π_{ab}, where Π_{ab} is a rank-three projector. It is
    diagonal in a product of two of the mutually unbiased bases, so it is a sum of three product projectors.

  So σ_Q is separable, with 36 explicit product terms besides identities (`check_blocks2.py`, part (e)).
  A shorter proof: every block M has purity Tr M²/(Tr M)² ≤ 0.117 < 1/8, so it is separable by the theorem of
  Gurvits and Barnum [GB] (a state on C³ ⊗ C³ with Tr ρ² ≤ 1/8 is separable).
- (b) If the range of a rank-two C lies in (U ⊗ U′ ⊗ U″)Q at three copies, tensored with the rest, then EP_{n−3}
  implies q_n(1/2, C) ≥ 0. The proof is that of Proposition 11. Since EP_1, EP_2 and EP_3 hold, the four-, five- and
  six-copy endpoint inequalities hold for every such C.
- (c) Suppose the range of C lies in S_1 ⊗ ⋯ ⊗ S_m ⊗ (C³)^{⊗r}, where the S_i sit on disjoint blocks of copies. Each S_i
  is a plane at one copy, or a local-unitary image of Λ at two copies, of Q at three copies, of the [[4, 2, 2]] code
  space at four copies, of the [[5, 3, 2]] code space at five copies (Proposition 14) or of the five-qutrit code space
  at five copies (Proposition 13). Then EP_r
  implies
  q_n(1/2, C) ≥ 0. Some blocks may use the co-range instead of the range. Since EP_r holds for r ≤ 3, the endpoint
  inequality holds for every n, with no unproven input, whenever at most three copies are left outside the blocks.
  The vectors in the range may be entangled between the blocks and the remaining copies.

*Proof of (c).* Remove the blocks one at a time, as in Proposition 11.
- After contracting a block with e_k and f̄_k, each vector φ_k still has Schmidt rank at most two. Its range still lies
  in the product of the remaining blocks tensored with the rest.
- For a plane at one copy, the compression is a PPT operator on C² ⊗ C³, hence separable; this is thin-copy removal.
- For a block on the co-range side, use the compression on the second side. It is separable because ρ is invariant
  under exchanging the two sides.
- At the end r copies remain, and EP_r applies. ∎

**Proposition 13 (removing a five-qutrit code block).** The one-sided compression of ρ^{⊗5} to the code space of the
five-qutrit code of Section 8 is separable. So if the range of a rank-two C lies in a local-unitary image of that code
space at five copies, tensored with the rest, then EP_{n−5} implies q_n(1/2, C) ≥ 0. With at most three free copies
(Proposition 12(c)), this gives the endpoint at up to eight copies for such operators.

*Proof.* On the second side, the 243-dimensional space splits into the 81 syndrome sectors of the code (the joint
eigenspaces of the four stabiliser generators), each of dimension 3. The compression is block diagonal in them; the
blocks are 3 × 3 operators on (logical) ⊗ (sector). Every block has purity Tr M²/(Tr M)² ≤ 0.1135 < 1/8, so it is
separable by [GB], and so is their sum. The code-sector block is (4/9)I − (19/288)F, the coefficients of
Proposition 8 (`check_blocks2.py`, part (f)). Removal is as in Proposition 11. ∎

**Stabiliser codes in general.** Write qutrit Pauli operators, up to phase, as X^a Z^b with a, b ∈ F_3^b, and wt(P)
for the number of copies where P is not the identity. Let S be a stabiliser group on b copies with |S| = 3^{b−k}, N(S)
the Paulis that commute with S, and Q the code space (dimension 3^k). The Paulis E with a given syndrome form a coset
of N(S). They split into 9^k classes E·S.

**Proposition 14 (stabiliser code blocks).**
- (a) The one-sided compression σ_Q = (P_Q ⊗ I) ρ^{⊗b} (P_Q ⊗ I) is block diagonal over the syndrome sectors of the
  second side. In each sector, the partial transpose of its block is diagonal in a basis of maximally entangled
  vectors labelled by the classes E·S with that syndrome. The eigenvalue at E·S is

      p(E·S) = |S|^{−1} Σ_{s ∈ S} (−1/2)^{b − wt(Es)}.

- (b) Let L ⊇ S be a maximal set of commuting Paulis (|L| = 3^b). In a fixed sector, the uniform distribution on the
  3^k classes contained in a coset E_0 L is the spectrum of a separable block. So a sector block is separable whenever
  its weights p are a nonnegative combination of these uniform distributions.
- (c) For the [[4, 2, 2]] code and for the [[5, 3, 2]] code (stabilisers X^{⊗5} and Z ⊗ Z ⊗ Z ⊗ Z ⊗ Z², dimension 27),
  this holds in every sector, with exact rational coefficients. So both compressions are separable. If the range (or
  the co-range) of a rank-two C lies in a local-unitary image of one of these code spaces, tensored with the rest,
  then EP_{n−4} (respectively EP_{n−5}) implies q_n(1/2, C) ≥ 0. With at most three free copies, this gives the
  endpoint at up to seven (respectively eight) copies for such operators, and both blocks can be used in
  Proposition 12(c).
- (d) For the codes [[b, b − 2, 2]] with stabilisers X^{⊗b} and Z^{c_1} ⊗ ⋯ ⊗ Z^{c_b} (all c_i ≠ 0 and Σ c_i ≡ 0 mod 3),
  the compression is PPT if and only if b ≤ 6.

*Proof.*
- **(a)** One copy gives ρ = I − F/2 = (5/6) I ⊗ I − (1/6) Σ_{P ≠ I} P ⊗ P^†, summed over the eight non-identity Paulis,
  because F = (1/3) Σ_P P ⊗ P^†. So ρ^{⊗b} = Σ_P c(P) P ⊗ P^† with c(P) = (5/6)^b (−1/5)^{wt(P)}.
  - P_Q P P_Q is 0 unless P ∈ N(S), and then it equals P P_Q. So σ_Q = Σ_{P ∈ N(S)} c(P) P P_Q ⊗ P^†.
  - These P preserve the syndrome sectors on the second side, so σ_Q is block diagonal.
  - In a sector with isometry B, the partial transpose of the block is Σ_{P ∈ N(S)} c(P) G(P), where
    G(P) = (P|_Q) ⊗ conj(B^† P B). The unitaries G(P) commute, since the two commutation phases cancel, and
    G(s) is a scalar for s ∈ S.
  - Their joint eigenvectors are labelled by the characters χ_E(P) = ω^{⟨E, P⟩} (ω = e^{2πi/3}, ⟨·,·⟩ the symplectic
    form), for E with the given syndrome, one for each class E·S. Applying a logical Pauli on Q permutes them, and
    each is maximally entangled.
  - The eigenvalue at E·S is Σ_{P ∈ N(S)} c(P) χ_E(P). Poisson summation over N(S), whose symplectic complement is S,
    turns this into |S|^{−1} Σ_{s ∈ S} ĉ(Es). Here ĉ(E) = Σ_P c(P) χ_E(P) factorises over copies. On one copy,
    Σ_P (−1/5)^{wt(P)} χ_E(P) is −3/5 for E = I and 6/5 otherwise, so ĉ(E) = (−1/2)^{b − wt(E)}.
- **(b)** Let a be a common eigenvector of the operators P|_Q with P ∈ L (a stabiliser state of the code), and f a
  common eigenvector of the operators conj(B^† P B) with P ∈ L.
  - Averaging (a ⊗ f)(a ⊗ f)^* over conjugation by all G(P) is a mixture of local unitaries. So it gives a separable
    operator, diagonal in the basis of (a), with weights |⟨Φ_{E·S}|a ⊗ f⟩|².
  - Since a ⊗ f is an eigenvector of G(P) for every P ∈ L, these weights vanish outside the classes of one coset E_0 L.
  - Writing the basis vectors as logical Paulis applied to one of them shows the weights are equal to 3^{−k} on that
    coset.
- **(c)** For the [[4, 2, 2]] code there are 40 such sets L, so each sector has 360 cosets. A linear program finds
  a nonnegative combination of them equal to p in every sector. We then solve for the coefficients exactly on the
  support, and check positivity and the identity in rational arithmetic (`check_code_blocks.py`, part (b)).
  - For [[5, 3, 2]] (1120 sets L, 30240 cosets per sector) we use a margin instead of an exact LP. Put d = 3^k.
    Averaging the cosets through one class over all L gives (𝟙 + d e_i)/(d(d + 1)), because a nonzero logical vector
    lies in a fraction 1/(d + 1) of the sets L. So every r ≥ 0 with (d² + d) min r ≥ Σ r is a nonnegative combination
    of cosets.
  - The LP is solved for p − 10^{−3}·𝟙. Its coefficients are rounded to rationals, the remainder r is computed exactly,
    and r passes this test in every sector (`check_code_blocks_margin.py`; the same check also certifies [[4, 2, 2]]
    and the five-qutrit code).
  - Removal is as in Proposition 11.
- **(d)** S consists of the nine operators s_{x,y} = X^{x(1, …, 1)} Z^{y c}, and all eight non-identity ones have full
  weight.
  - On copy i, E_i s_{x,y} is the identity for exactly one pair (x, y), since c_i ≠ 0. So the classes E·S correspond
    to distributions of the b copies over the nine cells (x, y), and p = (1/9) Σ_{cells} (−1/2)^{N_cell}, where N_cell
    is the number of copies in that cell. Every distribution occurs.
  - Since (−1/2)^m ≥ 1 − (3/2)m for every integer m ≥ 0, with equality for m ∈ {0, 1}, we get p ≥ 1 − b/6, with
    equality when no cell holds two copies (possible for b ≤ 9).
  - For b ≥ 7, put one copy in each of seven cells and the rest in an eighth. Then 9p ≤ 7(−1/2) + 1 + 1 < 0. ∎

The formula in (a) reproduces the dense spectra of the [[3, 1, 2]] and [[4, 2, 2]] compressions to 5 × 10^{−15}
(`check_code_blocks.py`, part (a)). The weights for [[4, 2, 2]] are 1/3, 7/12, 17/24, 5/6 and 43/48, and the
compression is not inside the Gurvits–Barnum ball. The test in (b) also certifies the [[3, 1, 2]] code again. For
[[5, 3, 2]] the weights are 1/6, 5/12, 13/24, 2/3, 35/48, 19/24 and 85/96. We have not tested [[6, 4, 2]], whose
smallest weight is 0.

**NPT code blocks still give no counterexample** (numerical). For b ≥ 7 the compression is NPT, so operators whose
range lies in the code space form a slice of the b-copy problem where a negative value is not ruled out by (d).
- By (a), on this slice q_b(1/2, C) = Σ_E p(E) |⟨Φ_E|ψ_C⟩|². Applying logical Paulis to one Bell vector turns it
  into Σ_{u,v} p(E_0 · X̃^u Z̃^v) |Tr((X^u Z^v)^† M)|², where M is a rank-two operator on the 3^k-dimensional logical
  space and X̃_j, Z̃_j are logical Paulis.
- For [[7, 5, 2]] (c = (1, 1, 1, 1, 1, 2, 2); 243 logical dimensions), the weight −1/6 occurs on 4–5% of the
  classes of each sector. We minimised q/‖C‖² over rank-two M (Adam, eight restarts in each of the nine sectors).
  The smallest value found is +0.114, so this slice stays well away from zero (`search_npt_code_block.py`).

## References

- [WZ] T. Wu, Q. Zou, Sharp Plücker geometry for three-copy Werner distillation, arXiv:2608.02647.
- [T2] T. S. Tomar, Equality cases and every-copy results for the Werner endpoint inequality, Zenodo 23718037 (2026).
- [BGH] K. Bharti, R. Gajjala, T. Haug, Two-copy nondistillability of Werner states: sharp partial-trace inequalities
  and finite-copy extensions, arXiv:2607.24479.
- [FGP] Fu, Gao, Park, A solution to 2-copy distillability of Werner states, arXiv:2607.21367.
- [SC] Z. Song, L. Chen, A partial-trace matrix inequality and Werner-state distillability, arXiv:2607.23416.
- [GB] L. Gurvits, H. Barnum, Largest separable balls around the maximally mixed bipartite quantum state, Phys. Rev. A
  66, 062311 (2002).

## Code

Python 3 with NumPy and SciPy; `check_mps.py` also uses mpmath. Outputs are in `logs/`.
- `code/werner.py`: shared helpers (as in `contraction-notes/`).
- `code/check_qc.py`: Section 1. It checks Proposition 1 and the trace identity, the two failing rules and the
  adversarial operators in `data/`.
- `code/check_qc_fast.py`: the five-copy QC adversary of Section 1 (needs PyTorch; `python check_qc_fast.py 5 120`, about
  fifteen minutes).
- `code/check_averaging.py`: Section 2. It averages exactly over the four mutually unbiased bases of C³ (a 2-design).
- `code/check_monogamy.py`: Section 3.
- `code/check_mps.py`: Section 4. It checks Proposition 4 against dense q_n, shows a Perron operator, and re-evaluates a
  stored search output with 60 digits.
- `code/check_blocks.py`: Section 6.
- `code/check_activation.py`: Section 5. It checks the closed form and evaluates the two-pair minimum at points of the
  slice.
- `code/check_signature.py`: Section 7. It prints the eigenvalues of the filter form for random operators and for the
  activation example.
- `code/check_codes.py`: Section 8. It computes A_n and B_n exactly, builds the five-qutrit code, and checks the formula
  against the dense q_5.
- `code/check_cones.py`: Section 9. It checks the site-dependent transfer identity, then runs the decomposition fits
  for KL_1 and KL_2 with their controls, the uncompressed test and the Schmidt-number test. Parts (b)–(e) need PyTorch
  and take about half an hour.
- `code/check_l1_six.py`: Section 9, the L_1(6) fits (needs PyTorch, about half an hour).
- `code/check_lmD_adversary.py` (`python check_lmD_adversary.py m D 16 80 60 0`) and `code/check_kl2_refit.py`: the
  dimension threshold and the KL_2 adversary in Section 9 (needs PyTorch; `check_kl2_refit.py` reads the saved worst
  target `adv_lmD_m2_D9_s0.npy`).
- `code/check_lm6_adversary.py` (`python check_lm6_adversary.py m 16 80 60 0` for m = 1, 2) and
  `code/check_l2six_adversary_npt.py` (`python check_l2six_adversary_npt.py 16 80 60 0.3`): the primal adversaries for
  L_1(6) and L_2(6) in Section 9 (needs PyTorch).
- `code/check_removal.py`: Section 9, removal by compression. It checks the identity against the dense q_4, shows
  that the conjugated compression gives a different value, and prints the triple margins for random C.
- `code/check_blocks2.py`: Section 10. It computes the two spectra exactly (SymPy), tests larger candidate blocks and
  checks the product bound against the dense q_n. It also checks the separable decompositions of Propositions 11 and
  12 and the purity test of Proposition 13.
- `code/check_code_blocks.py`: Proposition 14. It checks the weight formula against dense spectra, builds the exact
  separability certificates for the [[3, 1, 2]] and [[4, 2, 2]] codes (needs SymPy), and computes the minimum weights
  of the [[b, b − 2, 2]] family.
- `code/check_code_blocks_margin.py` (with `code/stab_fast.py`): Proposition 14(c) by the margin method, for
  [[4, 2, 2]], [[5, 3, 2]] and the five-qutrit code. It takes about five minutes.
- `code/search_npt_code_block.py`: the search on the [[7, 5, 2]] slice after Proposition 14 (needs PyTorch;
  `python search_npt_code_block.py 7`, about ten minutes).
