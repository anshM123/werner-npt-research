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
- **Section 5.** Which single vectors of a plane certify it in the Schur argument ("universal anchors"). Numerics,
  and three exact facts for every n ≥ 3: the W-plane and every plane near it contain no anchor (Proposition 7); a
  single excitation is an anchor exactly when it lives on at most two copies (Proposition 8); and the endpoint still
  holds, with a margin, on all planes near the W-plane (Corollary 9). Also: when a tensor product is an anchor
  (Proposition 10), a second open family of anchor-free planes (Proposition 11), and strict anchors with explicit
  open balls for every n (Theorem 12). Finally, anchor-free planes all of whose vectors have maximally mixed one-copy
  marginals (Proposition 13), and an exact reduction of the endpoint to range planes of that balanced kind
  (Theorem 14): anchors alone cannot give every n.
- **Section 6.** Local structure at zero planes: spin-flip wedge zeros have nonnegative second variation for every n
  (Proposition 15); the GHZ plane and generic local-qubit planes have open neighbourhoods on which the endpoint holds,
  and m(L) ≤ 2σ(L) for every plane (Theorem 16); at projector-core zeros the second variation is equivalent to the
  endpoint on one copy fewer (Proposition 17); Kronecker forms of flattenings do not detect anchors (Proposition 18).

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

Two remarks sharpen this.
- **Orthogonality is not needed.** a·u_1 v_1* + b·u_2 v_2* has rank at most two for any vectors, so the same equivalence
  holds with the supremum over all unit vectors.
- **A subexponential bound suffices.** Tensoring two configurations multiplies their ratios, so ρ_{N+M} ≥ ρ_N ρ_M.
  By Fekete's lemma, lim (log ρ_N)/N = sup_N (log ρ_N)/N. Hence the endpoint holds for every N if and only if
  ρ_N = e^{o(N)}. For example, √ρ_{N+1} ≤ √ρ_N + c with a fixed c would suffice.

The constant is attained: the GHZ plane span{|0^N⟩, |1^N⟩} (odd N) and the zeros of Section 1 give ratio 1. The
naive Cauchy–Schwarz bound loses an exponential factor (4^N on the GHZ plane), so a proof must keep the crossing.

The rank-one scale 2^{−N} is the natural one: the zeros of Section 1 are approached at that scale. This criterion is
in the spirit of the tube-radius equivalences that circulate for this problem. It does not prove positivity.

## 5. Universal anchors

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
(0.0014). At n = 5 (matrix-free Lanczos, tolerance 1e−5), two random real vectors (Gaussian entries, so uniform on
the real unit sphere; not complex Haar samples) are universal with gaps 0.0267 and 0.0268, a complex Haar-random
vector is universal with gap 0.0270, the GHZ vector Σ_a |aaaaa⟩ is universal with gap 0.0217, and the Dicke vector with two letters 1 is not
(−0.131).

**Two exact facts.** Put A = (M^{⊗n})^Γ and T_p = A^{−1/2} R_p^Γ A^{−1/2}, so that p is universal exactly when
T_p ≤ 1.
- Anchors tensorise: T_{p⊗p′} = T_p ⊗ T_{p′}. Since T_{p′} always has the eigenvalue 1, a vector that is not universal
  stays so after tensoring with anything.
- Spectra in [−1, 1] are closed under tensor products; this is the block condition of Section 2. Spectra in [−3, 1]
  are preserved by appending product copies, because a product copy has spectrum {1, −1/3}.

**Anchors do not cover every plane.** Every vector of span{Σ_s |0…1_s…0⟩, Σ_s |0…2_s…0⟩} is locally unitarily
equivalent to W, so this plane contains no universal anchor at n = 3. Yet q_n ≥ 2^{1−n}‖X‖² on it, for every co-range,
in our computations up to n = 7.

The rest of this section makes these observations exact for every n.

**Notation.** e_t is the string with the letter 1 at copy t and 0 elsewhere, and r_S the string with the letter 1 on a
pair S of copies. Vectors of row ⊗ column are written |row⟩|column⟩. Put c_n = 2^{1−n} and W_n = Σ_t |e_t⟩ (not
normalised; Y_p does not depend on the normalisation of p). For a local unitary U^{⊗n}, Y_{U^{⊗n}p} is unitarily
equivalent to Y_p, because M^{⊗n} commutes with (U ⊗ Ū)^{⊗n}.

**Proposition 7 (no anchor near the W-plane).** Let n ≥ 3 and put

    a_n = Σ_{|S|=2} Σ_{t∈S} |r_S⟩|e_t⟩,   b_n = Σ_{|S|=2} Σ_{t∉S} |r_S⟩|e_t⟩.

Then the Gram matrix of Y_{W_n} on (a_n, b_n) is

    c_n (n − 1)(n − 2) · [[0, −1], [−1, 2(n + 1)]].

It has negative determinant. With z = 2(n + 1)a_n + b_n this gives λ_min(Y_{W_n}) ≤ −β_n, where
β_n = 2^{3−n}(n + 1)(n − 2)/(n[8(n + 1)² + n − 2]). At n = 3 the matrix is [[0, −½], [−½, 4]] and β_3 = 4/387.

*Proof.* Every entry of M^{⊗n} is a product over the copies of δ_{ab}δ_{yx} − ½δ_{ay}δ_{bx}.
- **The rows of MV_p.** Let p = Σ_t w_t|e_t⟩. Going through the copies one by one, the row (r_S, e_t) of MV_p is
  −c_n w_{S∖t}[t ∈ S] times the vacuum column ⟨0…0|, and has no other entry. D_p preserves the number of excitations
  (Lemma 1 of `relative-anchors/`) and equals c_n Σ_t w_t² on the vacuum. Hence, on these states,

      R_p^Γ[(S, t), (S′, t′)] = (c_n/Σ_u w_u²) w_{S′∖t} w_{S∖t′} [t ∈ S′][t′ ∈ S].

- **A = (M^{⊗n})^Γ on these states.** It is c_n I on the states (S, t) with t ∈ S. On the states with t ∉ S it is
  block diagonal over the three-sets S ∪ {t}, with blocks c_n(3I + J) (J the all-ones matrix). It does not couple the
  two families.
- **Counting** (w_t = 1). ‖a_n‖² = n(n − 1) and ‖b_n‖² = n(n − 1)(n − 2)/2. Fix a state of a_n. Then n states of a_n
  and n − 2 states of b_n have a nonzero R^Γ entry with it. A state of b_n has 2(n − 2) such partners in b_n. So
  G_A = c_n diag(n(n − 1), 3n(n − 1)(n − 2)) and G_{R^Γ} = c_n [[n(n − 1), (n − 1)(n − 2)], [(n − 1)(n − 2), (n − 1)(n − 2)²]].
  Their difference is the claim. ∎

**Consequences.**
- Every unit vector of the W-plane L_n = span{Σ_s |0…1_s…0⟩, Σ_s |0…2_s…0⟩} is U^{⊗n}W_n/√n for a one-copy unitary U
  fixing |0⟩. So λ_min(Y_p) ≤ −β_n on all of L_n: for every n ≥ 3 the W-plane contains no universal anchor.
- p ↦ Y_p is continuous, so the anchors form a closed set and the planes containing no anchor form an open set. By
  the above this open set is nonempty for every n ≥ 3. So vectors that are not anchors, and planes without an anchor,
  have positive measure. Such planes cannot be confined to finitely many proper algebraic families.
- Quantitatively, ‖M^{⊗n}‖ = 1 and D_p ⪰ 2^{−n} give ‖Y_p − Y_q‖ ≤ K_n‖p − q‖ for unit p, q, with
  K_n = 2·3^n(2^n + 4^n). (Use ‖D_p − D_q‖ ≤ 2‖p − q‖, ‖D_p^{−1} − D_q^{−1}‖ ≤ 2·4^n‖p − q‖ and ‖X^Γ‖ ≤ 3^n‖X‖.) So
  no plane L′ with η(L_n, L′) < β_n/(2K_n) contains an anchor. Here η(L, L′) is the smallest ‖U − U′‖ over
  isometries U, U′ : C² → (C³)^{⊗n} onto L and L′.

**Proposition 8 (anchors among single excitations).** Let p ≠ 0 be a single excitation over a product vector, that is
p = Σ_t x_t with x_t a single excitation at copy t. Then p is a universal anchor if and only if x_t ≠ 0 for at most two
copies. Quantitatively, local unitaries fixing the product vector bring p to Σ_t w_t|e_t⟩ with w_t ≥ 0 and
Σ_t w_t² = 1, and for any three distinct copies i, j, k,

    λ_min(Y_p) ≤ −4c_n π/(4(n − 1) + π),   π = w_i² w_j² w_k².

*Proof.* If p lives on at most two copies, it is a two-copy vector tensored with a product vector. It is then an
anchor by Proposition 10 of `single-excitation-range/` (its exact two-copy certificate, with product copies appended as
in the second exact fact above). Conversely, put a(w) = Σ_{|S|=2} Σ_{t∈S} w_{S∖t}|r_S⟩|e_t⟩ and b = |r_{ij}⟩|e_k⟩. The
formulas in the proof of Proposition 7 give ⟨a, Y_p a⟩ = 0, ‖a‖² = n − 1, ⟨a, Y_p b⟩ = −2c_n w_iw_jw_k,
⟨b, Y_p b⟩ = 4c_n and b ⊥ a. If Y_p ⪰ 0, then ⟨a, Y_p a⟩ = 0 forces Y_p a = 0, which fails when w_iw_jw_k ≠ 0. The vector
a + (w_iw_jw_k/2)b gives the bound. ∎

**Planes inside the single excitations.** Let L be a plane of single excitations over a product vector, with
orthonormal frame U and one-copy blocks U_s : C² → C², so that Σ_s A_s = I with A_s = U_s^*U_s. A vector Ug has
excitation U_s g at copy s. By Proposition 8, L contains an anchor iff some nonzero vector of L lives on two copies,
that is, iff max_{i<j} λ_max(A_i + A_j) = 1. These are exactly the planes of Proposition 10 of
`single-excitation-range/`. Equivalently: with r the number of rank-two blocks, m the number of rank-one blocks and M
the largest number of rank-one blocks sharing a kernel line, the fewest copies a nonzero vector of L can occupy is
r + m − M, and L contains an anchor iff r + m − M ≤ 2. If κ = 1 − max_{i<j} λ_max(A_i + A_j) > 0, every unit vector of L
has three weights with x_1x_2x_3 ≥ ρ = κ²/(n(n − 2)²) (x_s = ‖U_s g‖²), so λ_min(Y_p) ≤ −4c_nρ/(4(n − 1) + ρ) on all of L.
For the remaining single-excitation planes the relative anchors of `relative-anchors/` take over.

**Corollary 9 (the endpoint near the W-plane).** For every n ≥ 3 and every operator C whose range lies in a plane L′
with η(L_n, L′) = η, with any co-range,

    q_n(C) ≥ (μ_n − 2η)‖C‖_F²,   μ_n = min{2^{1−n}, 1/(2n)}.

So the endpoint holds on an explicit open set of planes around the W-plane, and for η < β_n/(2K_n) these planes contain
no anchor.

*Proof.* For the W-plane itself, the proof of Theorem 2 of `single-excitation-range/` gives explicit sector bounds:
with h = n − c holes, the hole form is at least (4h − 2)/n per unit norm, so q_n ≥ (2h − 1)/(n2^h)‖X_c‖² for
1 ≤ h ≤ n − 1; the sector c = n gives (n − 1)/(2n); and the vacuum sector is exactly 2^{1−n}. The minimum over these is
μ_n (for n ≥ 3, (2h − 1)/2^h is smallest at h = 1 or h = n − 1, and (2n − 3)/(n2^{n−1}) ≥ 2^{1−n}). For n ≥ 4, μ_n = 2^{1−n}
is attained by the vacuum sector. For a plane L′ with frame U′, the compressed form (U′^* ⊗ I)M^{⊗n}(U′ ⊗ I) moves by at
most 2‖U − U′‖, since ‖M^{⊗n}‖ = 1. ∎

At n = 3 the exact minimum on the W-plane is 1/4 (`single-excitation-range/code/check_wplanes.py`), so there the margin
is 1/4 − 2η. The radii are small (β_3/(2K_3) = 1/752328); the point is that both sets are open.

**Proposition 10 (anchors of tensor products).** Let p = p_1 ⊗ ⋯ ⊗ p_m be a product over disjoint blocks of copies,
put a_j = max{0, −λ_min(T_{p_j})} and order them a_(1) ≥ a_(2) ≥ ⋯ (a missing second value counts as 0). Then p is a
universal anchor if and only if every p_j is one and a_(1)a_(2) ≤ 1.

*Proof.* T_p = ⊗_j T_{p_j}, and every T_{p_j} has the eigenvalue 1. Necessity: take the eigenvalue 1 in all factors
but one, or the two most negative eigenvalues in two factors and 1 elsewhere. Sufficiency: a positive eigenvalue of
T_p uses an even number of negative factors; pair their magnitudes in decreasing order. Each pair is at most
a_(1)a_(2) ≤ 1 (note a_(2)² ≤ a_(1)a_(2)), and the positive factors are at most 1. ∎

For b = (|10⟩ + |01⟩)/√2, which is an anchor with λ_min(T_b) = −11/9 (`single-excitation-range/`, Proposition 10),
b ⊗ b is not an anchor, while b ⊗ |0…0⟩ is (the vacuum has a = 1/3).

**Proposition 11 (a second family of anchor-free planes).** Let n ≥ 4, let b be as above, let E be any two-dimensional
subspace of the antisymmetric space Λ²C³ ⊂ C³ ⊗ C³, and put L_n = b ⊗ E ⊗ |0⟩^{⊗(n−4)} (copies 1–2, 3–4, then the
rest). Then every unit p ∈ L_n has λ_min(Y_p) ≤ −(5/98)·2^{4−n}, and every C with range in L_n, with any co-range,
has q_n(C) ≥ 2^{1−n}‖C‖_F², with equality attained.

*Proof.*
- **Witness.** With A_2 = (M^{⊗2})^Γ and B_b = R_b^Γ, an exact computation shows that
  v = 2|00⟩|11⟩ + 3|01⟩|01⟩ + |01⟩|10⟩ + |10⟩|01⟩ + 3|10⟩|10⟩ + 2|11⟩|00⟩ (row ⊗ column) satisfies B_b v = −(11/9)A_2 v,
  ‖v‖² = 28 and ⟨v, A_2 v⟩ = 9. Every unit a ∈ Λ²C³ has two equal Schmidt coefficients, so it is locally unitarily
  equivalent to b and has a transformed witness v_a. Since D and R factor over products of row vectors,
  Y_{b⊗a} = A_2 ⊗ A_2 − B_b ⊗ B_a, and ⟨v ⊗ v_a, Y_{b⊗a}(v ⊗ v_a)⟩ = 81 − 121 = −40 with ‖v ⊗ v_a‖² = 784. Each further
  vacuum copy, with the witness |0⟩|0⟩, halves the form and keeps the norm.
- **Margin.** The compression of M^{⊗n} to L_n ⊗ (columns) is the tensor product of D_b, the compression of M^{⊗2} to
  E ⊗ C⁹, and D_{|0⟩} = ½ on each further copy. For normalised b, D_b has spectrum {1/2 (three times), 3/4 (five times),
  1}. All two-dimensional E are equivalent under U ⊗ U (Λ²C³ is the conjugate defining representation of U(3) times the
  determinant), and for E_0 = span{|01⟩ − |10⟩, |02⟩ − |20⟩} the compression has spectrum {1/4 (three times), 1/2 (twice),
  3/4 (nine times), 1 (four times)}. The minimum is ½ · ¼ · 2^{−(n−4)} = 2^{1−n}, attained by products of minimisers. ∎

At four copies every vector of L_4 factors across the cut {1,2}|{3,4}, while every vector of the W-plane has Schmidt
rank two across every two-two cut, so the two families are not related by invertible maps on single copies (even with
permutations). Their neighbourhoods (as in Corollary 9) are further open sets of anchor-free planes on which the
endpoint holds.

**Theorem 12 (strict anchors for every n).** Let g_2 = (|00⟩ + |11⟩ + |22⟩)/√3, g_3 = (|000⟩ + |111⟩ + |222⟩)/√3, and let
p_n = g_2^{⊗n/2} for even n and p_n = g_3 ⊗ g_2^{⊗(n−3)/2} for odd n ≥ 3. Then the kernel of Y_{p_n} is exactly
p̄_n ⊗ (C³)^{⊗n}, and Y_{p_n} ⪰ δ_n on its orthogonal complement, with δ_n = 2^{−n}/216.

*Proof.* For m = 2, 3 put A_m = (M^{⊗m})^Γ, B_m = R_{g_m}^Γ and P_m = (I − |g_m⟩⟨g_m|) ⊗ I. Exact elimination
(`code/check_anchor_families.py`) gives A_m − B_m − P_m/64 ⪰ 0 with nullity exactly 3^m, and B_m + (215/216)A_m ≻ 0.
The nonzero eigenvalues of I − T_{g_m} = A_m^{−1/2}(A_m − B_m)A_m^{−1/2} are those of (A_m − B_m)^{1/2}A_m^{−1}(A_m − B_m)^{1/2},
hence at least (1/64)/λ_max(A_m) = (1/64)(2/3)^m ≥ 1/216. So T_{g_m} has the eigenvalue 1 with multiplicity exactly
3^m, and all its other eigenvalues lie in [−215/216, 215/216]. By tensorisation the same holds for T_{p_n}, with
eigenvalue 1 of multiplicity 3^n (a product of eigenvalues equals 1 only if every factor is 1). Congruence by A^{1/2},
with λ_min(A) = 2^{−n}, gives the gap. ∎

Consequences, with K_n as in the consequences of Proposition 7.
- **Open sets of anchors.** The lifted matrix Ŷ_p = Y_p + |p̄⟩⟨p̄| ⊗ I is ⪰ 0 exactly at anchors, ≻ 0 exactly at strict
  anchors, and (K_n + 2)-Lipschitz in unit p. So every unit q with min_θ ‖q − e^{iθ}p_n‖ < 2^{−n}/(432(K_n + 2)) is a strict
  anchor with gap at least δ_n/2.
- **Margins through strict anchors.** If p is a unit vector with Y_p ⪰ g > 0 off its forced kernel (g ≤ 1), then
  q_n(C) ≥ g·2^{−n−1}‖C‖_F² for every C whose range lies in a plane containing p. (In an orthonormal row frame (p, u) the
  compressed form is [[D_p, X], [X^*, D_u]] with 2^{−n} ⪯ D_p, D_u ⪯ I and D_u − X^*D_p^{−1}X ⪰ g; after normalising the
  diagonal blocks the off-diagonal block has norm at most √(1 − g), and 1 − √(1 − g) ≥ g/2.) Every plane through p_n
  has margin at least 4^{−n}/432.
- **Measure.** For every n ≥ 3, both the anchors and the non-anchors contain nonempty open sets, and so do the planes
  with and without an anchor. The boundary of the set of anchors lies in the real algebraic set det Ŷ_p = 0, which is
  not everything because of p_n; so almost every anchor is strict. Whether anchors become typical as n grows is open.

**Remark (certifying a given plane).** f(p) = λ_min(Ŷ_p) is ≥ 0 exactly at anchors, and H(L) = max_{p∈L, |p|=1} f(p)
satisfies |H(L) − H(L′)| ≤ (K_n + 2)η(L, L′). A mesh of the projective line of L with rigorous eigenvalue bounds therefore
certifies either H(L) < 0 (no anchor near L) or H(L) > 0 (a strict anchor, and q_n ≥ H(L)·2^{−n−2}‖C‖² on all nearby
planes) whenever H(L) ≠ 0. H(L) < 0 says nothing against the endpoint: the W-plane and the planes of Proposition 11
have H < 0 and a positive margin.

**Proposition 13 (anchor-free planes with maximally mixed marginals).** Let X|j⟩ = |j + 1⟩ and Z|j⟩ = ω^j|j⟩ on C³
(ω = e^{2πi/3}, labels in F_3), let c(a, b) = (a, b, a + b, a + 2b) ∈ F_3⁴, and let C(a, b) be the 48-letter word that
repeats each coordinate of c(a, b) twelve times. For t ≥ 0 define the isometry from three copies into n = 51 + t copies

    ℰ_t v = (1/3) Σ_{a,b∈F_3} (X^aZ^b)^{⊗3} v ⊗ |C(a, b)⟩ ⊗ |a⟩^{⊗t},

and let L_n = ℰ_t(span{w_1, w_2}) with w_ℓ = (|ℓ00⟩ + |0ℓ0⟩ + |00ℓ⟩)/√3 (the three-copy W-plane). Then every unit vector
of L_n has all one-copy marginals equal to I/3, and λ_min(Y_p) ≤ −1087/198144 for every unit p ∈ L_n.

*Proof.*
- **Marginals.** Distinct codewords differ in at least three groups, so the nine branches stay orthogonal after any
  one copy is traced out, and the marginal is the average of the branch marginals. On a core copy this is the twirl
  (1/9) Σ_{a,b} X^aZ^b ρ (X^aZ^b)^* = (Tr ρ) I/3. On a flag copy each symbol occurs in three of the nine branches.
- **A compression.** Let y be the column flag word with six letters 1 and then six letters 2 in each group of twelve,
  and the letter 1 on the appended flags, and fix the row flags to 0…0. Copy by copy, with the entries
  δ_{ab}δ_{yx} − ½δ_{ay}δ_{bx}:
  - the column flag string y spans an invariant block of D_{ℰv}: two different codewords differ on a whole group,
    where y takes both values 1 and 2, so the cross terms vanish. On the block, D_{ℰv} = (1/9)(D_v + N_v) with
    N_v = Σ_{(a,b)≠0} ε_{ab} D_{U_{ab}v} and 0 ≤ ε_{ab} ≤ 2^{−18}, since a nonzero codeword matches y at 18 of the 48
    flags and the zero codeword at none (U_{ab} = (X^aZ^b)^{⊗3});
  - the rows of M V_{ℰv} with row flags 0…0 and column flags y receive only the branch (0, 0), with coefficient 1/3;
  - on these states A contributes the factor 1 at every flag.

  So, with J the embedding of the three-copy row-column space at these flags,

      J^* Y_{ℰv} J = A_3 − [M_3 V_v (D_v + N_v)^{−1} V_v^* M_3]^Γ = Y_v + (ΔR)^Γ,
      ΔR = M_3 V_v [D_v^{−1} − (D_v + N_v)^{−1}] V_v^* M_3,   0 ⪯ ΔR,   ‖ΔR‖ ≤ 8 · 8ε · 8 = 512ε,

  using D_v ⪰ I/8 for unit v, ‖N_v‖ ≤ 8ε (ε = 2^{−18}) and ‖M_3‖ = 1.
- **The W witness survives.** For unit v in the W-plane, Proposition 7 gives (after the local unitary) a vector z_v
  with ‖z_v‖² = 387 and ⟨z_v, Y_v z_v⟩ = −4. Its coefficient matrix is 8J − 7I, with singular values 17, 7, 7, so
  ‖(|z_v⟩⟨z_v|)^Γ‖_1 = 31² = 961. Trace duality bounds the change of the form by 512ε · 961, so
  λ_min(Y_{ℰv}) ≤ (−4 + 961/512)/387 = −1087/198144. ∎

So no property of the one-copy marginals of a vector (full local rank, or even maximally mixed marginals) makes it an
anchor.

**Theorem 14 (a balanced reduction of the endpoint).** For k ≥ 1 and even r ≥ 2 put N = k + 4r and

    ℰ_{k,r} v = (1/3) Σ_{a,b∈F_3} (X^aZ^b)^{⊗k} v ⊗ |c(a, b)^{[r]}⟩,

where c^{[r]} repeats each coordinate of c r times. This is an isometry from k copies into N copies, and every unit
vector in its image has all one-copy marginals I/3 (as in Proposition 13). Then EP_k holds if and only if, for every
even r ≥ 2, q_N(C) ≥ 0 for every C of rank at most two whose range lies in the image of ℰ_{k,r}, with any co-range.

*Proof.*
- **(⇒)** On one group of r flag copies, compress M^{⊗r} to the rows span{|b^r⟩} (a qutrit). Nonconstant column strings
  give a diagonal block with entries 2^{−(number of matches)}. On constant columns the block is
  σ_r = 2^{−r}|Ω⟩⟨Ω| + Σ_{a≠b} |ab⟩⟨ab| (here r even is used: the cross entries are (−½)^r). It is separable, because
  E_ξ[|ξ⟩⟨ξ| ⊗ |ξ̄⟩⟨ξ̄|] = |Ω⟩⟨Ω| + Σ_{a≠b} |ab⟩⟨ab| for ξ with independent unit phases (the 27 choices of cube roots
  suffice). Tensoring the four groups and compressing the rows to the nine codewords keeps separability between the
  codeword label and the flag columns: the compression is Σ_i |α_i⟩⟨α_i| ⊗ β_i with β_i ⪰ 0. Writing C = ℰ_{k,r}C̃, this
  gives q_N(C) = (1/9) Σ_i Σ_f q_k(C_{i,f}), where C_{i,f} is the slice at flag column f of
  (Σ_λ ᾱ_{iλ} U_λ) C̃ (I ⊗ β_i^{1/2}). Every slice has rank at most two, so EP_k gives q_N(C) ≥ 0.
- **(⇐)** Let y be the column flag word with r/2 letters 1 and then r/2 letters 2 in each group, and for a k-copy
  operator C put C′ = ℰ_{k,r} C J_y^* with J_y z = z ⊗ |y⟩. Then C′ has the same rank and norm, its range lies in the
  image of ℰ_{k,r}, and the same flag computation gives the exact identity

      q_N(C′) = (1/9)[ q_k(C) + 2^{−3r/2} Σ_{(a,b)≠(0,0)} q_k((X^aZ^b)^{⊗k} C) ].

  Since |q_k| ≤ ‖·‖_F², a C with q_k(C) = −d‖C‖_F² < 0 gives q_N(C′) ≤ (1/9)(−d + 8·2^{−3r/2})‖C‖_F² < 0 once
  2^{3r/2} > 8/d. ∎

**Consequence.** The endpoint for every number of copies is equivalent to the endpoint on rank-two operators whose range
lies in one of these balanced encoded spaces, all of whose vectors have maximally mixed one-copy marginals. By
Proposition 13 such ranges need not contain an anchor, although every vector in them has maximally mixed one-copy
marginals. So a proof of the endpoint for every n cannot come from anchors alone, and in particular not from showing
that vectors of full local rank are anchors.

## 6. Local structure at zero planes

q_n vanishes on many rank-two operators, so a proof of the endpoint has to understand q_n near its zeros. This section
records what can be said locally. Throughout, Q = span{|0⟩, |1⟩} at each copy, J|0⟩ = |1⟩, J|1⟩ = −|0⟩, Θv = J^{⊗n}v̄ on
Q^{⊗n}, and for a plane L with orthogonal projector Π_L we write ρ_ν(L) = Tr_{other copies} Π_L (a 3 × 3 matrix of trace
2), τ(L) = max_ν λ_min(ρ_ν(L)), σ(L) = Σ_ν λ_min(ρ_ν(L)), and m(L) = min q_n(C)/‖C‖_F² over C ≠ 0 with range in L (any
co-range).

**Proposition 15 (spin-flip wedges).** For every n and independent a, b ∈ Q^{⊗n}, the rank-two operator
C_0 = a(Θb)^* − b(Θa)^* has q_n(C_0) = 0, is a critical point of q_n on the manifold of rank-two operators, and q_n has
nonnegative second variation there.

*Proof sketch.* After the column change C ↦ CJ^{⊗n}, C_0 becomes the skew matrix ab^T − ba^T, and on qubit rows and
columns M^{⊗n} becomes the tensor product of the local symmetric projections, which kills antisymmetric tensors; so
q_n(C_0) = 0. For a set S of copies at which row and column carry the letter 2, the residual M^{⊗n}C_0 has the block
(−½)^{|S|} 𝒫(Tr_S C_0). It vanishes for |S| even, and otherwise has row and column outside Q^{⊗n}, so the residual is
orthogonal to the tangent space and C_0 is critical. In the second variation the masks decouple. Even masks contribute
squared norms. For an odd mask of size s, a symplectic basis of Θ_S (2^{s−1} pairs, each a one-qubit spin flip) writes
the contribution as an average of models with one effective qutrit. Each such model is nonnegative on rank-two
operators: the qubit projection I − ½|Ω_2⟩⟨Ω_2| is separable (one half the sum of the projectors onto |01⟩, |10⟩, |+−⟩,
|−+⟩, |+i +i⟩, |−i −i⟩), and the one-qutrit endpoint form is nonnegative on rank-two operators. ∎

**Theorem 16 (positive neighbourhoods of local-qubit planes).**
- (a) At the GHZ zero C_* (a = |0…0⟩, b = |1…1⟩), the kernel of the second variation is exactly the tangent space of the
  smooth family of zeros U K (J^{⊗n})^* U^* (U a product of local frames of two-dimensional subspaces, K skew of rank
  two), of real dimension 4·2^n − 6 + 4n, and the transverse second variation is at least 2^{1−n}. The range compression
  of L_* = span{|0…0⟩, |1…1⟩} has the single zero C_*, and its other eigenvalues are at least 2^{1−n}. Hence (a
  Morse–Bott estimate near the zero, compactness away from it) every plane in some neighbourhood of L_* satisfies the
  endpoint, with any co-range.
- (b) The same holds for every plane inside a fixed product of local two-dimensional subspaces, outside a proper real
  algebraic subset: two polynomial rank conditions, both satisfied at L_*.
- (c) For every plane L, m(L) ≤ 2σ(L). Near the planes in (a) and (b), also m(L) ≥ c τ(L) for some c > 0, so there the
  margin vanishes exactly on planes inside a product of local two-dimensional subspaces.

The radii in (a) and (b) and the constant c are not explicit.

*Proof of the upper bound in (c).* If σ(L) ≥ ½ the bound is trivial (m ≤ ‖M^{⊗n}‖ = 1). Otherwise let Q_ν be the top two
eigenvectors of ρ_ν(L), P = ⊗_ν Q_ν, U a frame of L and A = U^*PU. The union bound gives ‖I − A‖ ≤ σ(L), so
U′ = PUA^{−1/2} is a frame of a local-qubit plane L′ with ‖U − U′‖² ≤ 2(1 − √(1 − σ)) ≤ 2σ. Let C′ = U′X be a unit wedge
zero of L′ (Proposition 15) and C = UX. Then C − C′ has co-range inside that of C′, so it is tangent at the critical zero
C′, and q_n(C) = q_n(C − C′) ≤ ‖U − U′‖² ≤ 2σ(L). ∎

**Proposition 17 (projector-core zeros).** For vectors a, b on k copies and P = |0⟩⟨0| + |1⟩⟨1| on one more copy,
Z = P ⊗ ab^* is a critical zero of q_{k+1} on rank-two operators, with residual −|2⟩⟨2| ⊗ M^{⊗k}(ab^*). Along the rank-two
curve C(t) = Z + t|0⟩⟨1| ⊗ (xb^* + ay^*),

    q_{k+1}(C(t)) = t² q_k(xb^* + ay^*).

Every rank-two operator on k copies has the form xb^* + ay^*, so nonnegative second variation at all these zeros on
k + 1 copies is equivalent to EP_k: a local statement at all zeros already contains the whole problem. In the other
direction, suppose a has a range margin h > 0, that is q_k(X) ≥ h‖X‖² for every X of rank at most two whose range lies in
a plane containing a. (For the strict anchors of Theorem 12 one can take h = 4^{−k}/432.) Then the second variation at Z
has kernel exactly the tangent space of the family P_{Q′} ⊗ (rank one) and transverse gap min{h, 2^{−k}}. So every plane
near Q ⊗ span{a} satisfies the endpoint, with margin comparable to 2 − λ_max(Tr_{first copy} Π_L). ∎

**Proposition 18 (Kronecker data do not see anchors).** At n = 3 let g = Σ_a |aaa⟩, h = |012⟩ + |120⟩ + |201⟩, w_ℓ the
unnormalised vectors of Proposition 13, A = span{g, h}, and B_ε = span{w_1 + εg, w_2 + εh} with ε = 10^{−8}.
- A contains the strict anchor g (Theorem 12), and every vector of A has one-copy marginals I/3.
- B_ε lies within ε of the W-plane, so it contains no anchor and q_3 ≥ (1/4 − 2ε)‖C‖_F² on it (Proposition 7 and its
  consequences).
- Yet at every copy the stacked 6 × 9 coefficient matrix of both planes has rank 6 (for B_ε a 6 × 6 minor equals ε³), so
  every flattened pencil of either plane has the Kronecker form [I_3, zI_3, 0].

A genericity argument gives anchored and anchor-free planes with identical flattened Kronecker forms at every
bipartition for every odd n ≥ 3. So the Kronecker forms of the flattenings cannot decide whether a plane contains an
anchor. ∎

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
- `code/check_anchors_matrixfree.py`: the smallest eigenvalue of Y_p off its kernel, matrix-free (Lanczos), used for the
  n = 5 values (`python code/check_anchors_matrixfree.py 5 rand 1`, etc.; log in `logs/check_anchors_matrixfree_n5.log`).
- `code/check_anchor_families.py`: Theorem 12 and Proposition 11 in exact arithmetic (python-flint): the two
  certificates for g_2 and g_3 (with their ranks and component sizes), the spectrum of D_b, the witness v, and the
  spectrum of the compression to E_0.
- `code/check_local_structure.py`: Section 6 numerically: the projector-core identities (k = 1, 2), the wedge zeros
  (value, criticality, nonnegative second variation and its kernel dimension at n = 1, 2) and m(L) ≤ 2σ(L) on 300 planes.
- `code/check_kronecker_example.py`: Proposition 18 (exact ranks and the minor ε³; anchors numerically) and the range
  compressions of the GHZ plane and of a generic local-qubit plane (n = 2, 3).
- `code/check_balanced_encoding.py`: Proposition 13 and Theorem 14: the flag factors of the proofs (exactly, for all
  flag strings when r = 2, 4 and on samples when r = 12), the reduction identity by direct partial traces at
  (k, r) = (1, 2) and (2, 2), the marginals of encoded vectors, and the numbers 387, 961 and −1087/198144.
- `code/check_anchor_obstruction.py`: Propositions 7 and 8 in exact rational arithmetic for n = 3, …, 7: the rows of
  MV_p live on the vacuum column, the Gram matrix on (a_n, b_n), and the weighted identities at rational weights.

## References

- [WZ] T. Wu, Q. Zou, Sharp Plücker geometry for three-copy Werner distillation, arXiv:2608.02647.
