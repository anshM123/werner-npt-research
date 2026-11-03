# Product-containing planes for every n, and a reduction to projectors

Ansh Mishra and Aryan Senthilkumar, November 2026.

These notes record four facts about the n-copy endpoint form, and some numerics:
- **Section 1.** The form is nonnegative on every operator whose range lies in a plane through a product vector. The
  co-range is arbitrary, and n is arbitrary. This is the product-ray theorem of Wu and Zou, which they state for three
  copies; their argument works for every n. We also show the form is positive definite there whenever the plane is
  thick at one copy.
- **Section 2.** The same holds with anchors built from maximally entangled qutrit pairs.
- **Section 3.** The endpoint for every number of copies is equivalent to its restriction to rank-two orthogonal
  projectors.
- **Section 4.** A negative witness forces the rank-two minimum to fall exponentially below the rank-one scale 2^{−N}.
- **Section 5 (numerical).** Which single vectors of a plane certify it in the Schur argument.

**None of this settles the all-copy problem.**

## Setting

For an operator X on (C^3)^{⊗n},

    q_n(X) = Σ_J (−½)^{|J|} ‖Tr_J X‖_2² = ⟨X, Φ_n(X)⟩,   Φ_n = (id − ½ Tr(·) I)^{⊗n},

and EP_n is the statement that q_n(C) ≥ 0 for every C of rank at most two. On one copy,
q_1(X) = ⟨vec X, M vec X⟩ with M = I − ½|Ω⟩⟨Ω| on C^3_row ⊗ C^3_col and Ω = Σ_a |aa⟩. Its blocks with respect to the row
factor are M_ab = ⟨a|M|b⟩ = δ_ab I − ½|a⟩⟨b|. On n copies the matrix is M^{⊗n}. Γ denotes the transpose on the row
factor and F the swap. Sections 1 and 2 work in any local dimensions ≥ 2.

## 1. Planes through a product vector

**Theorem 1.** Let p be a product vector and L a two-dimensional subspace containing p. Then q_n(X) ≥ 0 for every
operator X with range in L, for every n. The co-range is arbitrary, so X can be any element of L ⊗ (C^3)^{⊗n}*. The
same holds when the co-range lies in such a plane.

Wu and Zou prove this for three copies [WZ, Thm. 7.3], by a Schur complement and a completely copositive order
interval. Their argument does not use the number of copies, and we write it out for every n.

*Proof.* By local unitaries take p = |0⟩^{⊗n}, and write L = span{p, u} with u ⊥ p.

*One copy.* Put P = |0⟩⟨0|, Q = I − P, Z = P − Q and D = M_00 = I − ½P > 0. Define R_ab = M_a0 D^{−1} M_0b. A direct
computation gives

    R = ½ |vec Z⟩⟨vec Z| + P ⊗ Q,   A := M^Γ = I − ½F,   B := R^Γ = ½ (Z ⊗ Z) F + P ⊗ Q.

The three spaces |00⟩, span{|0i⟩, |i0⟩} (i ≠ 0) and Q ⊗ Q are invariant under A and B. On them:

| Invariant space | A − B | A + B |
|---|---|---|
| \|00⟩ | 0 | 1 |
| span{\|0i⟩, \|i0⟩} | [[0, 0], [0, 1]] | [[2, −1], [−1, 1]] |
| Q ⊗ Q | I − F | I |

So −A ≤ B ≤ A, with A > 0.

*n copies.* Compress M^{⊗n} to L on the row factor. In the basis p, u this is the block operator [[D^{⊗n}, K], [K*, H]]
on the co-range space, with K = ⟨p|M^{⊗n}|u⟩ and H = ⟨u|M^{⊗n}|u⟩. Since all blocks factorise over the copies, its
Schur complement is

    H − K* (D^{⊗n})^{−1} K = ⟨u| (M^{⊗n} − R^{⊗n}) |u⟩,

a partial expectation on the row factor. The matrix T = A^{−1/2} B A^{−1/2} is a Hermitian contraction, so

    (M^{⊗n} − R^{⊗n})^Γ = A^{⊗n} − B^{⊗n} = (A^{1/2})^{⊗n} (I − T^{⊗n}) (A^{1/2})^{⊗n} ≥ 0.

Whenever W^Γ ≥ 0, the operator ⟨u|W|u⟩ is positive semidefinite, because ⟨u ⊗ y|W|u ⊗ y⟩ = ⟨ū ⊗ y|W^Γ|ū ⊗ y⟩. So
the Schur complement is positive semidefinite, hence so is the compressed block operator, and that is the claim.
Finally q_n(X*) = q_n(X), which gives the co-range version. ∎

This is positivity of the whole difference A^{⊗n} − B^{⊗n}, not of the terms of an expansion.

**Strictness.** On the mixed blocks det(B − λA) = (λ − 1)(3λ + 1)/4, so T has eigenvalues 1 and −1/3 only. A tensor
product of these equals 1 only if every factor does. So the kernel of A^{⊗n} − B^{⊗n} is K^{⊗n}, with

    K = ker(A − B) = (|0⟩ ⊗ C^3) ⊕ Sym²(Q C^3).

Suppose a nonzero X = pv* + uw* has q_n(X) = 0. Then the Schur form vanishes in some direction z ≠ 0, so
ū ⊗ z ∈ K^{⊗n}. At each copy s, apply the vectors spanning K^⊥:
- **Vectors |i⟩|0⟩ (i ≠ 0).** Either u is locally |0⟩ at s, or z has no |0⟩ component there.
- **The vector |12⟩ − |21⟩.** In the second case, ⟨1|_s ū ⊗ ⟨2|_s z = ⟨2|_s ū ⊗ ⟨1|_s z. Then the local support of
  u at s lies in span{|0⟩, r_s} for a single vector r_s.

So the range of every zero has local dimension at most two at every copy. This gives:

**Corollary 2.** If L contains a product vector and has local dimension three at some copy, then q_n is positive
definite on {X : ran X ⊆ L}, with arbitrary co-range.

The zeros include P ⊗ |0…0⟩⟨b| for a qubit projector P, with b arbitrary and possibly entangled. A generic plane in
(C^3)^{⊗n}, n ≥ 2, contains no product vector, so Theorem 1 covers a thin set of planes. By Theorem 1, a negative rank-two
witness must have a completely entangled range and a completely entangled co-range.

## 2. Anchors with maximally entangled pairs

The proof of Theorem 1 runs block by block. Split the copies into blocks and let the anchor be p = ⊗_j p_j, with p_j a
vector on block j. For each block put

    V_j = p_j ⊗ I,   D_j = V_j* M_j V_j,   R_j = M_j V_j D_j^{−1} V_j* M_j,

where M_j is the endpoint matrix of the block. If D_j > 0 and −M_j^Γ ≤ R_j^Γ ≤ M_j^Γ for every block, then the Schur
complement factorises as before, and the tensor product of the contractions is again a contraction.

For a pair of qutrit copies with p_j = Ω_3 = Σ_a |aa⟩:
- D_j = 2I + ¼|Ω_3⟩⟨Ω_3| > 0.
- The 81 × 81 matrices 352 (M_j ± R_j)^Γ are integer matrices. Both are positive semidefinite: the signs of their
  characteristic polynomials alternate, checked exactly.
- (M_j − R_j)^Γ has a kernel of dimension 9.

**Theorem 3.** Let p be a tensor product of single-copy vectors and maximally entangled states on disjoint pairs of
copies. For every n, q_n(X) ≥ 0 for every X whose range lies in a plane containing p, with any co-range. The same holds
with range and co-range exchanged.

The zero analysis of Section 1 is not claimed for these anchors.

## 3. Reduction to rank-two projectors

**Proposition 4.** Suppose q_n(C) < 0 for some rank-two C on n copies. Then q_N(P) < 0 for some orthogonal rank-two
projector P on N > n copies. So EP_n for every n is equivalent to the same statement for rank-two projectors alone.

*Proof.* The proof has two steps.

*Flattening.* Write C = σ_1 E_1 + σ_2 E_2 with E_i = u_i v_i*, and put G_ij = q_n(E_i, E_j). Its diagonal is positive,
since q_n is positive on rank one. A combination of the E_i is negative, so det G < 0. The form tensorises, so

    D = E_1 ⊗ E_2 − E_2 ⊗ E_1   has   q_{2n}(D) = 2(G_11 G_22 − |G_12|²) < 0.

D = UV* for isometries U = [u_1 ⊗ u_2, u_2 ⊗ u_1] and V = [v_1 ⊗ v_2, −v_2 ⊗ v_1].

*Flags.* On k further copies let f_0 = |0^k⟩ and f_1 = |1^k⟩. Put Z = (U ⊗ f_0 + V ⊗ f_1)/√2 and P_k = ZZ*. Then
Z*Z = I_2, so P_k is a rank-two orthogonal projector. Its flag blocks are UU*/2, D/2, D*/2 and VV*/2. The off-diagonal
flag units have q_k(f_0 f_1*) = 1, and they are q-orthogonal to the rest. The diagonal ones give 2^{−k} and (−½)^k.
Hence

    q_{2n+k}(P_k) = ½ q_{2n}(D) + 2^{−k}/4 · [q_{2n}(UU*) + q_{2n}(VV*) + 2(−1)^k Re q_{2n}(UU*, VV*)].

Φ has operator norm one and ‖UU*‖² = ‖VV*‖² = 2, so the bracket is at most 8 in absolute value. Hence q(P_k) < 0 once
2^k > 8/|q_{2n}(D)|. ∎

The construction raises the number of copies, so it says nothing at a fixed n. Our four-copy theorem
(`four-replica-notes/`) already covers every projector at n = 4. Flag constructions of this kind are standard; we have
not seen the reduction stated for this problem.

For a projector P_V onto V = span{v_1, v_2}, the swap trick gives

    q_n(P_V) = Tr[(P_V ⊗ P_V) Π_s (F_s − ½)] = Tr[W Π_{Sym²V}] − ⟨ω|W|ω⟩,   W = Π_s (I − ½F_s) > 0,

where ω = (v_1 ⊗ v_2 − v_2 ⊗ v_1)/√2. So the all-copy problem asks whether, for every plane V, the W-weight of its
antisymmetric vector is at most the W-weight of its three-dimensional symmetric square.

## 4. Amplification: a negative witness grows at a fixed exponential rate

Iterating the flattening step of Proposition 4 gives a quantitative criterion. Let

    m_N = min { q_N(C) : rank C ≤ 2, ‖C‖_2 = 1 }.

**Proposition 5.** Suppose m_n < 0, and let G be the Gram matrix of Proposition 4 for a negative witness on n copies.
Put A = G_11 G_22 and B = |G_12|². Then B > A ≥ 4^{−n}, and along N = 2^k n (k ≥ 1)

    2^N m_N ≤ −[(4^n B)^{N/(2n)} − (4^n A)^{N/(2n)}].

So the right side tends to −∞ exponentially in N. Consequently the endpoint holds for every number of copies if and only
if 2^N max(0, −m_N) grows subexponentially in N, for instance if m_N ≥ −K·2^{−N} for one constant K.

*Proof.* Rank-one positivity gives q_n(xy*) ≥ 2^{−n}‖x‖²‖y‖², so A ≥ 4^{−n}, and det G < 0 gives B > A. The flattened
witness D_1 = E_1 ⊗ E_2 − E_2 ⊗ E_1 has the orthonormal rank-one pieces F_1 = E_1 ⊗ E_2 and F_2 = −E_2 ⊗ E_1. Since
q_{2n} is the tensor square of q_n, their Gram matrix is [[A, −B], [−B, A]]. Flattening (F_1, F_2) again squares every
entry: after k steps the Gram matrix is [[A^{2^{k−1}}, −B^{2^{k−1}}], [−B^{2^{k−1}}, A^{2^{k−1}}]], on N = 2^k n
copies. Its smallest eigenvalue A^{2^{k−1}} − B^{2^{k−1}} is the value of q_N at a unit rank-two operator, namely
(F_1 + F_2)/√2 at that level. Multiply by 2^N = (4^n)^{2^{k−1}}. ∎

**Corollary 6 (a constant suffices).** For unit vectors with u_1 ⊥ u_2 and v_1 ⊥ v_2 on N copies, put
W_N = Π_s (I − ½F_s) on two replicas and

    ρ_N = sup |⟨u_1 ⊗ v_2|W_N|u_2 ⊗ v_1⟩|² / (⟨u_1 ⊗ v_1|W_N|u_1 ⊗ v_1⟩ ⟨u_2 ⊗ v_2|W_N|u_2 ⊗ v_2⟩).

Then EP_N ⟺ ρ_N ≤ 1. Also ρ_{N+1} ≥ ρ_N and ρ_{2N} ≥ ρ_N². So the endpoint holds for every N if and only if
sup_N ρ_N < ∞: an orthogonal crossed Cauchy–Schwarz inequality with any constant independent of N would suffice.

*Proof.* The entries of G for E_i = u_i v_i* are the three W_N-expectations in ρ_N, by the crossed pairing
q(|x⟩⟨y|, |u⟩⟨v|) = ⟨x ⊗ v|W|u ⊗ y⟩. For C = σ_1 E_1 + σ_2 E_2, the phases of u_2 are free, so the worst case of
σ_1²G_11 + σ_2²G_22 + 2σ_1σ_2 Re G_12 is nonnegative for all σ exactly when |G_12|² ≤ G_11 G_22. Appending a common
product copy multiplies G_12 by ½ and G_11 G_22 by ¼. One flattening step squares the ratio, as in Proposition 5. ∎

The constant is attained: the GHZ plane span{|0^N⟩, |1^N⟩} (odd N) and the zeros of Section 1 give ratio 1. The
naive Cauchy–Schwarz bound loses an exponential factor (4^N on the GHZ plane), so a proof must keep the crossing.

The rank-one scale 2^{−N} is the natural one: the zeros of Section 1 are approached at that scale. This criterion is
in the spirit of the tube-radius equivalences that circulate for this problem. It does not prove positivity.

## 5. Universal anchors (numerical)

Section 1 used only one vector of the plane. Here is the general version. For a unit vector p let
R_p = M^{⊗n} V D^{−1} V* M^{⊗n}, with V = p ⊗ I and D = V* M^{⊗n} V, and put Y_p = (M^{⊗n} − R_p)^Γ. If Y_p ⪰ 0, the
Schur argument of Theorem 1 goes through for every plane containing p, so q_n ≥ 0 on all operators whose range plane
contains p, with any co-range. Call such p a universal anchor. Y_p always vanishes on p̄ ⊗ (C^3)^{⊗n}, since
(M^{⊗n} − R_p)(p ⊗ x) = 0.

**Numerics** (`code/check_anchors.py` and the n = 4 runs noted below). We find two kinds of behaviour.
- **Universal:** at n = 2, every vector we tried, including W = |01⟩ + |10⟩. At n = 3, Haar-random vectors are
  universal, with kernel exactly p̄ ⊗ (C^3)^{⊗n} and a gap of about 0.045. So are products and GHZ vectors.
- **Not universal:** at n = 3, the W state Σ_s |0…1_s…0⟩, and more generally every single excitation Σ_s b_1…x_s…b_n
  over a product, has smallest eigenvalue −0.0239.

At n = 4, random vectors stay universal (gap 0.045), W is not (−0.0149), and GHZ is universal with a small gap
(0.0014).

**Two exact facts.** Put A = (M^{⊗n})^Γ and T_p = A^{−1/2} R_p^Γ A^{−1/2}, so that p is universal exactly when
T_p ≤ 1.
- Anchors tensorise: T_{p⊗p′} = T_p ⊗ T_{p′}. Since T_{p′} always has the eigenvalue 1, a vector that is not universal
  stays so after tensoring with anything.
- Spectra in [−1, 1] are closed under tensor products; this is the block condition of Section 2. Spectra in [−3, 1]
  are preserved by appending product copies, because a product copy has spectrum {1, −1/3}.

**Anchors do not cover every plane.** Every vector of span{Σ_s |0…1_s…0⟩, Σ_s |0…2_s…0⟩} is locally unitarily
equivalent to W, so this plane contains no universal anchor at n = 3. Yet q_n ≥ 2^{1−n}‖X‖² on it, for every co-range,
in our computations up to n = 7.

## Code

Python 3 with NumPy and SymPy. `code/check_pair_anchor.py` uses python-flint for the exact characteristic polynomials;
without it, it falls back to floating point. Outputs are in `logs/`.
- `code/endpoint.py`: the form q_n.
- `code/check_product_anchor.py`: Theorem 1, Corollary 2 and Theorem 3 on examples. It prints the smallest eigenvalue
  of q_n on {X : ran X ⊆ L}, for random planes through |0…0⟩ (positive), for planes inside a product of local qubits
  (zero), and for planes through anchors with maximally entangled pairs at n = 3, 4 (positive).
- `code/check_pair_anchor.py`: the exact block condition of Section 2.
- `code/check_projector_reduction.py`: the flattening and flag identities of Proposition 4, on a random rank-two C, and
  the Gram identity behind Proposition 5.
- `code/check_anchors.py`: Section 5 at n = 2, 3.

## References

- [WZ] T. Wu, Q. Zou, Sharp Plücker geometry for three-copy Werner distillation, arXiv:2608.02647.
