# Product-containing planes for every n, and a reduction to projectors

Ansh Mishra and Aryan Senthilkumar, November 2026.

These notes record four facts about the n-copy endpoint form, and some numerics:
- **Section 1.** The form is nonnegative on every operator whose range lies in a plane through a product vector. The
  co-range is arbitrary, and n is arbitrary. This is the product-ray theorem of Wu and Zou, which they state for three
  copies; their argument works for every n. We also show the form is positive definite there whenever the plane is
  thick at one copy.
- **Section 2.** The same holds with anchors built from maximally entangled qutrit pairs.
- **Section 3.** The endpoint for every number of copies is equivalent to its restriction to rank-two orthogonal
  projectors. Proposition 4′ does this with k + r copies instead of 2k + r, through an exact identity, and gives a
  rate: positivity on projectors at N copies implies q_k ≥ −8·2^{−(N−k)}‖·‖_F² on rank two at k copies.
- **Section 4.** A negative witness forces the rank-two minimum to fall exponentially below the rank-one scale 2^{−N}.
- **Section 5.** Which single vectors of a plane certify it in the Schur argument ("universal anchors"). Numerics,
  and three exact facts for every n ≥ 3: the W-plane and every plane near it contain no anchor (Proposition 7); a
  single excitation is an anchor exactly when it lives on at most two copies (Proposition 8); and the endpoint still
  holds, with a margin, on all planes near the W-plane (Corollary 9). Also: when a tensor product is an anchor
  (Proposition 10), a second open family of anchor-free planes (Proposition 11), and strict anchors with explicit
  open balls for every n (Theorem 12). Finally, anchor-free planes all of whose vectors have maximally mixed one-copy
  marginals (Proposition 13), and an exact reduction of the endpoint to range planes of that balanced kind
  (Theorem 14): anchors alone cannot give every n. Theorem 14′ balances the range and the co-range at once, the
  balancing keeps every anchor obstruction (Remark), and Theorem 14″ reduces the endpoint for every n to rank-two
  orthogonal projectors onto planes all of whose vectors have maximally mixed one-copy marginals.
- **Section 6.** Local structure at zero planes: spin-flip wedge zeros have nonnegative second variation for every n
  (Proposition 15); the GHZ plane and generic local-qubit planes have open neighbourhoods on which the endpoint holds,
  and m(L) ≤ 2σ(L) for every plane (Theorem 16); at projector-core zeros the second variation is equivalent to the
  endpoint on one copy fewer (Proposition 17); Kronecker forms of flattenings do not detect anchors (Proposition 18).
- **Section 7.** For balanced rank-two projectors, the standard weight-enumerator relations of quantum codes alone give
  q_n > 0 for n ≤ 8 (Proposition 19, exact certificates). They cannot give it at n = 9, and two-copy balance does not
  help (an exact negative enumerator at n = 10).
- **Section 8.** Zeros decide the sign. The endpoint at n holds if and only if every rank-two zero of q_n is a critical
  point, and if and only if the zeros form a set of codimension two (Proposition 20). So it follows from a
  classification of zeros: if every zero has the qubit form of Section 6 (Conjecture ZC_n, true at two copies), the
  endpoint holds at n (Corollary 21). Off the endpoint, EP_n follows if the critical-point equations of
  q_n(t, ·) have no rank-two solution for 0 < t < ½ (Proposition 22). A random-start census finds only zeros of the
  qubit form for n ≤ 6. It is enough that every zero be thin at some copy (Proposition 23), with no assumption on
  fewer copies.
- **Section 9.** Wu and Zou's sharp rank-two bound for positive operators appears to extend to every rank-two
  operator and the whole Werner line: q_n(t, C) ≥ (1 − t)^{n−1}(‖C‖_F² − t‖C‖_1²), at the endpoint
  2^n q_n(C) ≥ (σ_1 − σ_2)² (Conjecture TN_n). It is attained for every ratio of singular values, implies the endpoint
  and the sharp curve, and forces every zero to have equal singular values (all known zeros do); random-start
  minimisation finds no violation for n ≤ 5 or with other local dimensions, and it is proved for every n when a
  product vector is a common left and right singular vector (Proposition 26). At two copies it holds for every
  rank-two operator in all local dimensions, by the refined partial-trace inequality of `first-order-bound/`
  (Proposition 28). Two natural inductions on it fail
  (slice by slice, and through the cross term of two orthogonal pure states). Its qualitative part already suffices:
  if every rank-two zero has equal singular values, the endpoint holds (Corollary 25); all known zeros do. In a
  purified form with a reference qubit, its one-copy case is Uhlmann's theorem. A sharp merging budget (Conjecture FB:
  merging two sites loses at most the factor (1 − 2t)/2 of the gap) would give it for every n, and at any t in
  (1/3, ½) NPT bound entanglement; FB is proved with two sites (Proposition 28) and for every number of sites on the
  class of Proposition 26 (Proposition 29), and is open from three sites on.

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

**Proposition 4′ (a dilation without flattening).** Let C = l_1r_1^* + l_2r_2^* on k copies. On r further copies let
f_0 = |0^r⟩ and f_1 = |1^r⟩, and put w_i = l_i ⊗ f_0 + r_i ⊗ f_1, P = w_1w_1^* + w_2w_2^*, D_0 = Σ_i l_il_i^* and
D_1 = Σ_i r_ir_i^*. Then

    q_{k+r}(P) = 2 q_k(C) + 2^{−r} q_k(D_0 + (−1)^r D_1).

Let C have singular values σ_1 ≥ σ_2 > 0 and singular vectors u_i, v_i. Put λ = 2σ_1, let α_i ≥ β_i > 0 be the roots of
t² − λt + σ_i², and take l_i = √α_i u_i, r_i = √β_i v_i. Then w_1 ⊥ w_2 and |w_i|² = λ, so P = λΠ for a rank-two
orthogonal projector Π, and

    q_{k+r}(Π) ≤ 2q_k(C)/λ² + 4·2^{−r}.

Consequently, if q_N(Π) ≥ 0 for every rank-two orthogonal projector Π on N copies, then q_k(C) ≥ −8·2^{−(N−k)}‖C‖_F²
for every k < N and every C of rank at most two.

*Proof.* The flag blocks of P are D_0, C, C^* and D_1. On the flag copies, the entry of M^{⊗r} between |f_a⟩|f_b⟩ and
|f_c⟩|f_d⟩ (row, column) is 1 if (a, b) = (c, d) with a ≠ b, 2^{−r} if (a, b) = (c, d) with a = b, (−½)^r between
(0, 0) and (1, 1), and 0 otherwise. So q_{k+r}(P) = q_k(C) + q_k(C^*) + 2^{−r}[q_k(D_0) + q_k(D_1) + 2(−1)^r q_k(D_0, D_1)],
and q_k(C^*) = q_k(C) because M is real and invariant under the swap. For the bound, |q_k(X)| ≤ ‖X‖_F²,
‖D_0 ± D_1‖_F ≤ Tr D_0 + Tr D_1 = 2λ, and ‖C‖_F² ≥ σ_1² = λ²/4. The last statement follows with r = N − k, since
λ² ≤ 4‖C‖_F². ∎

Compared with Proposition 4, this uses k + r copies instead of 2k + r and gives a rate. Since P ⪰ 0, it also shows that
the endpoint for every n is equivalent to q_n(P) ≥ 0 for every positive semidefinite P of rank at most two. Writing
P = aa^* + bb^* and rescaling a and b separately, this is the scalar inequality

    q_n(aa^*, bb^*) ≥ −√(q_n(aa^*) q_n(bb^*))   for all vectors a, b,

where q_n(X, Y) is the polarised form. In terms of the marginals ρ^a_S = Tr_{S^c} aa^*, it reads
q_n(aa^*, bb^*) = Σ_S (−½)^{n−|S|} Tr[ρ^a_S ρ^b_S].

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

**Remark (added 9 November 2026: a negative margin never decays).** Padding with GHZ states is cheaper than flattening.
For g_r = (|0^r⟩ + |1^r⟩ + |2^r⟩)/√3 every partial trace of g_rg_r^* other than the empty and the full one has squared
norm 1/3, so q_r(g_rg_r^*) = (2/3)(1 + (−½)^r) + 2^{−r}/3 ≥ 0.6 for r ≥ 2. Since C ⊗ g_rg_r^* has the rank and the norm
of C and q is multiplicative, a witness with q_n(C) = −δ‖C‖_F² gives m_N ≤ −0.6δ for every N ≥ n + 2. Hence the
endpoint holds for every number of copies if and only if m_N ≥ −ε_N along some unbounded set of N with ε_N → 0; the
bound m_N ≥ −K·2^{−N} above is far more than needed. The same holds on the Werner line (q_r(t, g_rg_r^*) ≥ 2/3 for even
r) and inside the class of real balanced projectors, since g_r is balanced (`real-balanced-reduction/`, Corollary
"a vanishing error is enough"). This reformulates the problem; it proves nothing.

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

**Theorem 14′ (range and co-range balanced at once).** Index the eight nonzero linear forms ℓ on F_3². For λ ∈ F_3² let
c_λ ∈ F_3^{8r} list the values ℓ·λ, each repeated r times, and let d_λ = c_λ + (1, …, 1). With U_λ = (X^aZ^b)^{⊗k} for
λ = (a, b), put

    E_r v = (1/3) Σ_λ U_λ v ⊗ |c_λ⟩,   F_r v = (1/3) Σ_λ U_λ v ⊗ |d_λ⟩.

These are isometries from k copies into N = k + 8r copies, every vector in either image has all one-copy marginals
I/3, and for every operator C on k copies

    q_N(E_r C F_r^*) = (1/9) [ q_k(C) + 2^{−3r} Σ_{ν≠0} q_k(C U_ν) ].

So EP_k holds if and only if, for every r ≥ 1, q_N ≥ 0 on every operator of rank at most two whose range lies in the
image of E_r and whose co-range lies in the image of F_r.

*Proof.* The marginals are as in Proposition 13: each nonzero ℓ takes each value at three labels λ, and two different
codewords differ at 6r places. The 81 flag pair states |c_λ⟩ ⊗ |d_μ⟩ (row, column) are mutually q-orthogonal. A
diagonal entry is 2^{−a}, where a is the number of places at which the row and column flags agree: a = 0 when λ = μ,
and a = 3r when λ ≠ μ, because ℓ·(λ − μ) = 1 for exactly three of the eight ℓ. An entry between two different pair
states can be nonzero only if the row and column flags of the first state agree at every place where the two row
codewords, or the two column codewords, differ. That is 6r places, while a pair state has at most 3r agreements.
Finally q_k(U_λ C U_μ^*) = q_k(C U_μ^* U_λ) = q_k(C U_{λ−μ}), since q_k is invariant under X ↦ UXU^* for local unitaries
U and phases do not matter. An operator with range in the image of E_r and co-range in the image of F_r is E_r C F_r^*
with C = E_r^* X F_r of the same rank. For the converse, q_k(C) = −d‖C‖_F² < 0 gives
q_N(E_rCF_r^*) ≤ (−d + 8·2^{−3r})‖C‖_F²/9 < 0 once 2^{3r} > 8/d. ∎

**Remark (balancing keeps the anchor obstructions).** Let p be a unit vector on k copies and P = E_r p. The column
flag word y = 1…1 spans an invariant block of D_P, on which D_P = (D_p + N)/9 with 0 ⪯ N ⪯ 8·2^{−3r} I, and the rows
of MV_P with row flags 0…0 and column flags y come from the branch λ = 0 alone (the flag computation is the one in the
proof above). Hence the compression of Y_P to these flags is Y_p + (ΔR)^Γ with 0 ⪯ ΔR and
‖ΔR‖ ≤ ζ = 8·4^k·2^{−3r}/(1 + 8·2^k·2^{−3r}), using D_p ⪰ 2^{−k}I. A unit witness w with ⟨w, Y_p w⟩ = −β < 0 therefore
still gives a negative value once t(w)ζ < β, where t(w) ≤ 3^k is the squared trace norm of the coefficient matrix of
w (trace duality). So every non-anchor stays a non-anchor after encoding with r large: the balanced test spaces of
Theorems 14 and 14′ keep all the anchor obstructions.

**Theorem 14″ (balanced projectors suffice).** Keep the notation of Theorem 14′, let r ≥ 2 be even, let E_r^{(s)} be
E_r with every codeword shifted by s(1, …, 1) (s ∈ F_3), and put

    H_r v = (1/√3) Σ_{s∈F_3} E_r^{(s)} v ⊗ |s, s⟩,

an isometry from k copies into N = k + 8r + 2 copies. Every vector in its image has all one-copy marginals I/3. For
every operator Π on k copies, with ε = 2^{−3r}, η = 2^{−2r−2} and ζ = 2^{−8r−2},

    q_N(H_r Π H_r^*) = (2/27) q_k(Π) + ((2ε + η)/27) Σ_{ν≠0} q_k(Π U_ν) + (ζ/81) q_k(Σ_λ U_λ Π U_λ^*).

If Π is a rank-two orthogonal projector, so is H_rΠH_r^*, and the last two terms together are at most 2^{−2r−1} in
absolute value. With Proposition 4′: the endpoint holds for every number of copies if and only if q_N(Π) ≥ 0 for every
N and every rank-two orthogonal projector Π onto a plane all of whose vectors have maximally mixed one-copy marginals.

*Proof.* The 27 branch words α = (λ, s) have the codewords (c_λ + s(1, …, 1), s, s) ∈ F_3^{8r+2}, and
H_rΠH_r^* = (1/27) Σ_{α,β} U_λΠU_μ^* ⊗ |code_α⟩⟨code_β|. Two different codewords differ at 8r + 2, 5r + 2 or 6r places.
Among the 729 flag pair states, those with α ≠ β are isolated, by the argument of Theorem 14′ (they have at most 3r or
2r + 2 agreements), and their diagonal entries are 1 (s ≠ t, λ = μ: no agreement), ε (s ≠ t, λ ≠ μ: 3r agreements)
and η (s = t, λ ≠ μ: 2r + 2 agreements). The pair states with α = β have row and column flags equal; they do not
couple to the others, and among themselves the entries are ζ(−1)^{dist(α, α′)} = ζ, because all distances are even
when r is even. Counting: the 54 pairs of the first kind give (2/27)q_k(Π); the 432 of the second kind and the 216 of
the third give ((2ε + η)/27) Σ_ν q_k(ΠU_ν), each ν occurring equally often; and the α = β block gives
(ζ/729)q_k(3Σ_λ U_λΠU_λ^*). The marginals are as before: tracing out the last two flags removes the cross terms in s.
For the bound, |q_k(ΠU_ν)| ≤ 2, ‖Σ_λ U_λΠU_λ^*‖_F² ≤ 162, ε ≤ η and ζ ≤ 2^{−6r}η. ∎

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

## 7. Weight enumerators of balanced projectors

For a rank-two orthogonal projector Π on n copies and the 9^n Weyl operators E = ⊗_ν X^{a_ν}Z^{b_ν} (the weight of E is
the number of ν with (a_ν, b_ν) ≠ (0, 0)) put

    A_j = ¼ Σ_{wt E = j} |Tr ΠE|²,   B_j = ½ Σ_{wt E = j} Tr(ΠEΠE^*).

Φ multiplies a Weyl operator of weight j by (−½)^{n−j}, so

    q_n(Π) = (4/3^n) Σ_j (−½)^{n−j} A_j.

The usual relations of quantum codes hold: A_0 = B_0 = 1, Σ_j A_j = 3^n/2, B_j ≥ A_j ≥ 0, the MacWilliams identity
B(x, y) = (2/3^n) A(x + 8y, x − y) for A(x, y) = Σ_j A_j x^{n−j} y^j, and the coefficients of A(4x + 2y, x − y) are
nonnegative (this polynomial is (6^n/4) Tr[(Π ⊗ Π) ⊗_ν (xP_+ + yP_−)], with P_± the local symmetric and antisymmetric
projectors on two replicas). If every vector of ran Π has maximally mixed one-copy marginals, then ΠEΠ = 0 for every E
of weight one, so A_1 = B_1 = 0; the MacWilliams identity then gives Σ_j jA_j = 4n·3^{n−2}.

**Proposition 19.** Let Π be a rank-two orthogonal projector on n copies whose range consists of vectors with maximally
mixed one-copy marginals. Then q_3(Π) = 3/2 and

    q_4(Π) ≥ 11/12,  q_5(Π) ≥ 5/8,  q_6(Π) ≥ 9/16,  q_7(Π) ≥ 1239/5344,  q_8(Π) ≥ 6487/45408.

At n = 9 the relations above do not imply q_9(Π) ≥ 0: the linear program they define has optimum
−5515759/22425216, attained by an explicit rational enumerator.

*Proof.* At n = 3 the two moment equalities force A = (1, 0, 3/2, 11). For n = 4, …, 8 put F_n = A(−½, 1) and let S_j
be the coefficients of A(4x + 2y, x − y). On the affine space cut out by the equalities,

    F_4 = 297/16 + (9/4) A_2,
    F_5 = 1215/32 + (27/8)(B_2 − A_2),
    F_6 = 6561/64 + (81/16)(B_2 − A_2) + (1/8) S_6,
    F_7 = 2709693/21376 + (729/1336) A_5 + (67959/5344)(B_2 − A_2) + (891/2672)(B_4 − A_4) + (51/2672) S_6,
    F_8 = 14187069/60544 + (50301/30272) A_5 + (19683/60544) A_6 + (1400409/60544)(B_2 − A_2)
          + (62937/60544)(B_4 − A_4) + (6517/60544) S_8.

Every term after the constant is nonnegative, and q_n(Π) = (4/3^n) F_n. At n = 9 a dual certificate of the same kind
and a feasible enumerator with the same value give the stated optimum. ∎

**The bound at five copies is attained (added 9 November 2026).** Let S be the stabilizer group of the qutrit
five-qudit code [[5,1,3]]_3, generated by the cyclic shifts of X ⊗ Z ⊗ Z^{−1} ⊗ X^{−1} ⊗ I. Its code space is
three-dimensional, every vector in it has maximally mixed one- and two-copy marginals, and so every two-dimensional
subspace V has A_1 = B_1 = A_2 = B_2 = 0. Then F_5 = 1215/32, that is q_5(Π_V) = 5/8: the bound of Proposition 19 at
n = 5 is the minimum of q_5 over balanced rank-two projectors. This also follows from the compression of the
five-qutrit code computed in `compression-notes/` (A_5 = 4/9, B_5 = −19/288 there, so q_5 ≥ (A_5 + 2B_5)‖C‖_F² = (5/16)‖C‖_F²
for rank-two C on the code space, with equality on projectors); we had not connected the two. We checked q_5 = 5/8 by
direct partial traces on random subspaces (`code/global_reductions/perfect_code_check.py`). Gradient descent over balanced real planes from random
starts does not find these: it stops at q_5 ≈ 0.6336, on planes whose two-copy marginals are not maximally mixed. So the
hardest balanced planes at five copies are subspaces of a good quantum code.

More generally, for a two-dimensional subspace of a qutrit stabilizer code with one logical qutrit,
Π = P_code − |c⟩⟨c| with c an eigenvector of a logical Weyl class ℓ,

    q_N(t, Π) = (4/3^N) [ Σ_{E∈S} (1 − 3t)^{N−wt E} + ½ Σ_{E∈ℓ+S} (1 − 3t)^{N−wt E} ],

because |Tr ΠE|² is 4 on the stabilizer group, 1 on the cosets ℓ and −ℓ, and 0 elsewhere. This is the rank-two case of
the coset-weight formula of Proposition 14(a) in `compression-notes/` (checked against direct partial traces on random
five-qutrit codes, same script). This makes code subspaces cheap to scan at large N. A short simulated
annealing over codes (`code/global_reductions/stabilizer_anneal.py`, log `logs/stabilizer_anneal.log`) found no negative
value at the endpoint for N = 5, 7, 9; at five and seven copies it reaches the known zeros, and at nine copies it is too
weak to reach them, so it says little there. The most uniform codes are not witnesses either. For quantum MDS codes the weight enumerators are
fixed by the parameters; if the eight logical cosets have equal enumerators (as for the five-qutrit code), two-dimensional
subspaces of [[7,1,4]]_3 and [[9,1,5]]_3 codes, should such codes exist, have q_7/‖Π‖_F² = 21/64 and
q_9/‖Π‖_F² = 127/768 at the endpoint (`code/global_reductions/mds_code_values.py`, exact). So the negative enumerator of
Proposition 19 at n = 9 is not that of an MDS code. A code built to make the weight-three terms large does not give one either. The CSS
code with both classical codes equal to the ternary Reed–Muller code RM_3(2, 2) on the nine points of the affine plane
AG(2, 3) is a [[9, 3, 3]]_3 code whose normalizer has 96 elements of weight three (the lines of the plane) and whose
stabilizer sum Σ_{E∈S} (−½)^{9−wt E} = −28 is negative. Encoded rank-two operators C = Enc(C_L) have
q_9(C) = 3^{−9} Σ_{L} W_{L+S} |Tr(C_L L)|² over the 729 logical Weyl operators L, and minimising over rank-two C_L gives
q_9 ≥ 0.1497‖C‖_F² at the endpoint (`code/global_reductions/rm_code.py`, log `logs/rm_code.log`; floating point,
eight starts). So at nine copies, where the relaxation of Proposition 19 is negative, this balanced code stays clearly
positive.

**Stabilizer domination (added 9 November 2026).** Every stabilizer state is a graph state up to local Cliffords,
which do not change weights, and a Pauli acting on a graph state |g⟩ (adjacency Γ over F_3) is a Z-shift of it up to
phase. So the rank-two operators whose range and co-range are spanned by Pauli images of one stabilizer state are, up to
local unitaries, C = Σ_{i,j} c_{ij} Z^{x_i} G Z^{−y_j} with G = |g⟩⟨g|. In the Weyl expansion, Z^x G Z^{−y} is supported
on the coset Z^{y−x}T of the stabilizer group T, so q(C) splits over these cosets and is nonnegative for all such C
exactly when, for every shift d and every character v,

    |S(d, v)| ≤ S(d, 0),   S(d, v) = Σ_{a∈F_3^N} ω^{v·a} (1 − 3t)^{N − w_d(a)},   w_d(a) = |supp a ∪ supp(Γa + d)|.   (7.1)

Here S(d, 0) = 3^N q_N(t, Z^xGZ^{−y}) for any x, y with y − x = d, and S(d, v) is 3^N times the cross term between
Z^xGZ^{−y} and its conjugate by Z^v. In cut-rank form (r_J the F_3-rank of Γ restricted to rows J and columns J^c),
S(d, v) = 3^N Σ_J (−t)^{|J|} 3^{−r_J} [d_J ∈ col Γ_{J,J^c}] [v_{J^c} ∈ col Γ_{J^c,J}] ω^{θ_J}, a phase θ_J; for the plane
spanned by |g⟩ and Z^u|g⟩ this gives q/‖Π‖_F² = Σ_J (−t)^{|J|} 3^{−r_J} (1 + [u_{J^c} ∈ col Γ_{J^c,J}]).

**Conjecture SD_t.** For every N, every Γ, every d and every v ≠ 0, |S(d, v)| ≤ τ S(d, 0) with τ = t/(1 − t).

At t = ½ this is (7.1). Product states attain τ (a single-copy shift gives ratio exactly τ), and the GHZ and product
zeros give equality at the endpoint. SD_t implies the sharp bound of Section 9 on this class: for the plane,
q(Π) = 2q(G) + 2 Re S(0, u)/3^N ≥ 2(1 − τ) q(G) ≥ 2(1 − 2t)(1 − t)^{N−1} for Π the projector onto span{|g⟩, Z^u|g⟩}. Evidence (`code/global_reductions/
stabilizer_domination.py` and `graph_cut_rank.py`, log `logs/stabilizer_domination.log`): exhaustive over all Γ, d, v
at N = 2, 3, 4 for t = ½ (largest ratio exactly 1, thousands of equality cases), t = 0.45 and 0.4 (largest ratio exactly
τ = 9/11 and 2/3), and t = 0.35 at N = 3 (τ = 7/13); 300 random graphs with all d, v at N = 5 and t = ½ (largest ratio 1).
The class contains the minimisers of q_5 over balanced projectors (subspaces of the five-qutrit code) and every known
zero, so SD_t for every N would give the endpoint, and the sharp bound, on all of them. We have no proof beyond N = 4;
the cycle graphs are covered for every N by the transfer-matrix argument of `cyclic-code-planes/` (on their code spaces).
Two exact tools for an induction. Expanding each site factor in characters of (Γa + d)_i gives the dual form
S(d, v) = Σ_{s∈F_3^N} ω^{s·d} Π_i (3[s_i = 0 and (v + Γs)_i = 0] − t), in which d becomes a character and v a shift.
Deleting a vertex n (measure it in the Z basis: |g⟩ = 3^{−1/2} Σ_k |k⟩_n ⊗ Z^{kγ}|g′⟩ with γ = Γ_{n,·} and g′ the graph
state of Γ′ = Γ − n) and using the one-copy recursion for the slices gives, with d = (d′, d_n) and v = (v′, v_n),

    S_N(d, v) = Σ_{m∈F_3} [ ω^{m v_n} S′(d′ + mγ, v′) − t ω^{m d_n} S′(d′, v′ + mγ) ]                    (7.2)

(`code/global_reductions/graph_recursion.py`, log `logs/graph_recursion.log`: relative error below 10^{−13} on random
graphs). Bounding the two sums in (7.2) separately with SD_t on Γ′ does not close: the second sum lowers S_N(d, 0) by about
t S′(d′, 0), and an induction has to keep the joint values S′(d′ + mγ, v′ + m′γ) for all m, m′.

These are bounds for balanced projectors at a fixed n ≤ 8. Theorems 14′ and 14″ need N ≥ k + 18 copies, so
Proposition 19 gives the endpoint at no copy number. The negative optimum at n = 9 is the value of a relaxation; no
projector is constructed. Rank-two codes satisfy further constraints that these relations do not see, for instance
that Λ²(ran Π) is one-dimensional.

More purity does not rescue the relations (added 8 November 2026). If every vector of ran Π also has maximally mixed
two-copy marginals, then A_2 = B_2 = 0 as well, and still the relations do not imply q_10(Π) ≥ 0: the enumerator

    A = (1, 0, 0, 489/14, 315, 405, 201, 1467/14, 74601/14, 118019/7, 43920/7)

satisfies all of them and gives q_10 = −101/1792 (at n = 11 the value is −661705/3319296). Written for rank-three
projectors (A_j = 1/9 Σ|Tr ΠE|², B_j = 1/3 Σ Tr(ΠEΠE^*), B(x, y) = (3/3^n) A(x + 8y, x − y)), the same relations reach
−196723/472064 at n = 10, so per unit ‖Π‖_F² the rank-two relations are only about five times stronger. Weight
enumerators are two-replica data, while rank two is a three-replica condition. A floating-point linear program proposes
each vertex; the vertex is then recomputed from its active constraints and every relation is checked in exact arithmetic
(`code/global_reductions/check_enumerator_purity.py`).

## 8. Zeros decide the sign (added 8 November 2026)

Write X₂ for the set of operators of rank exactly two on (C^3)^{⊗n}. It is a connected real manifold of dimension
8·3^n − 8 (one orbit of GL × GL). Rank-one operators are strictly positive:

    q_n(ab^*) = ⟨ā ⊗ b̄| (I − ½F)^{⊗n} |ā ⊗ b̄⟩ ≥ 2^{−n} ‖a‖² ‖b‖²,

because M^Γ = I − ½F ≥ ½ I. So every nonzero zero of q_n of rank at most two lies in X₂. Let Z_n ⊂ X₂ be the set of
these zeros. The tangent space of X₂ at C is {XC + CY}, so C = AB^* is a critical point of q_n on X₂ exactly when
Φ_n(C)B = 0 and Φ_n(C)^*A = 0, that is, when CΦ_n(C)^* = 0 and Φ_n(C)^*C = 0.

**Proposition 20 (zeros decide the sign).** Fix n.
- (a) EP_n holds if and only if Z_n has codimension at least two in X₂.
- (b) EP_n holds if and only if every zero is critical: q_n(C) = 0 with C of rank two implies CΦ_n(C)^* = 0 and
  Φ_n(C)^*C = 0.

*Proof of (a).* Z_n is a closed semialgebraic subset of X₂, a finite union of smooth submanifolds of dimension at most
dim X₂ − 2. Two points of X₂ \ Z_n are joined by a path in X₂, and a generic perturbation of that path with the same
ends misses every one of these submanifolds, since 1 + (dim X₂ − 2) < dim X₂. So X₂ \ Z_n is connected. The form has no
zeros there, and it is positive at ab^* + εcd^* for small ε ≠ 0, so it is positive on all of X₂ \ Z_n. For the
converse (pointed out by an independent audit), let C be a zero and take a unit rank-one D = ub^* with b in the
co-range of C. The complex line C + wD stays in rank at most two, and q_n(C + wD) = 2 Re(w̄ β) + α|w|² with
α = q_n(D) ≥ 2^{−n}. If EP_n holds, C is a minimum, so β = 0 and the two real gradient components along this line have
independent differentials at C (their Hessian block is 2α I₂). Near C every zero solves these two equations, so Z_n
has codimension at least two there. ∎

*Proof of (b).* At a zero of a nonnegative function the derivative vanishes, which gives "only if". For "if" (this
short argument is due to an independent audit; an earlier version used real algebraic geometry), suppose
q_n(C) = −η < 0 for some C of rank two. With D = ub^* as above (b in the co-range of C, ‖u‖ = ‖b‖ = 1), every C + sD,
s real, has rank at most two, and q_n(C + sD) = αs² + 2βs − η with α = q_n(D) ≥ 2^{−n} > 0. This quadratic has two
distinct real roots. At a root, C + sD is not zero and not of rank one (rank-one operators have q_n > 0), so it lies in
X₂, and the derivative of q_n along the curve, 2(αs + β), is not zero. So C + sD is a zero that is not critical. ∎

So the endpoint is a statement about zeros only: there is no zero at which the form has a nonzero derivative. Every
zero in Section 6 is critical (Propositions 15 and 17), and by Proposition 17 the second variation at zeros already
contains the whole problem. Part (b) shows that the first variation does too. Neither part bounds how negative q_n
could be, and a proof that all zeros are critical cannot assume EP_n, since it must also exclude zeros at which q_n
changes sign.

**The known zeros.** For disjoint I, J with I ∪ J = [n] and J ≠ ∅, let 𝓑_J be the set of operators

    C = R_I ⊗ V K V^*,   K = (ab^T − ba^T)(J^{⊗|J|})^{−1},

with R_I of rank one on the copies in I (its vectors may be entangled), a, b ∈ (C^2)^{⊗|J|}, and V = ⊗_{s∈J} V_s a
product of isometries C^2 → C^3. They are zeros for every n: q_n factorises over tensor products, and q_{|J|}(VKV^*) = 0
by Proposition 15. An independent audit of this repository proves the converse for operators supported on a product of
local qubit planes, for every n. We checked the argument; in our words: on such operators Φ_n acts as the projection
Π^{⊗n} onto locally traceless operators (the parts of Φ_n(C) outside the qubit block are orthogonal to C), and right
multiplication by J^{⊗n} turns locally traceless into locally symmetric. So q_n(C) = 0 exactly when the polynomial
z^T CJ^{⊗n} z (z = z_1 ⊗ ⋯ ⊗ z_n, z_s ∈ C²) vanishes identically. With CJ^{⊗n} = f_1 g_1^T + f_2 g_2^T this says
f_1 g_1 + f_2 g_2 = 0 for multilinear polynomials, of degree one in each pair z_s. Put h = gcd(f_1, f_2), f_1 = ha,
f_2 = hb. Then a g_1 = −b g_2, so g_1 = bk and g_2 = −ak. Factors of multihomogeneous polynomials are multihomogeneous,
and every degree is 0 or 1, so the copies split into I, where h and k have degree one, and J, where a and b do; J is
nonempty because f_1 and f_2 are independent. Hence CJ^{⊗n} = (hk^T)_I ⊗ (ab^T − ba^T)_J, which is of the form above.

The same audit extends this to an arbitrary co-range. If only the range lies in a product of local qubit planes, then
on each copy the endpoint operator restricted to rows in the plane, I − ½|Ω_Q⟩⟨Ω_Q| on Q ⊗ C^3 with
Ω_Q = Σ_j q_j ⊗ q̄_j over an orthonormal basis of Q, is an orthogonal projection of rank five. So q_n(C) ≥ 0 for every
such C of any rank. For a fixed range plane span{h_I a_J, h_I b_J} with gcd(a_J, b_J) = 1, the zeros form a linear
space of dimension exactly 3^{|I|}, consisting of the operators
|h_I⟩⟨y_I| ⊗ K_J with y_I an arbitrary vector on the copies in I. We checked the kernel dimensions numerically (1 and 3
at two copies, 1 and 9 at three, for generic and factorised planes). So every zero whose range, or co-range, lies in a
product of local qubit planes is in some 𝓑_J.

**Conjecture ZC_n.** Every rank-two zero of q_n lies in some 𝓑_J.

At two copies this holds: Tomar's classification of equality cases places every two-copy zero in a product of local
planes [T2] (as quoted in `four-replica-notes/`; we have not re-checked [T2]), and the argument above then applies.
From three copies on it is open.

**Corollary 21.** ZC_n implies EP_n.

*Proof.* 𝓑_J has real dimension at most 4·3^i + 4·2^j + 4j − 10, with i = |I| and j = |J|: the rank-one factor gives
4·3^i − 2, the local planes 4j, the operators on a fixed product of planes whose product with J^{−1} is skew of rank
two 4·2^j − 6 (this condition depends only on the planes, since a change of basis in a plane multiplies J by a phase),
and the common scalar removes 2. The largest of these over j ≥ 1 is below 8·3^n − 10 = dim X₂ − 2 for every n ≥ 1
(for n ≥ 3 it is 4·3^{n−1} + 2 or 4·2^n + 4n − 6). Proposition 20(a) applies. ∎

**Proposition 22 (critical zeros along the Werner line).** For 0 ≤ t ≤ ½ put Φ_{n,t} = (id − t Tr(·) I)^{⊗n} and
q_n(t, C) = ⟨C, Φ_{n,t}(C)⟩, so that q_n = q_n(½, ·). Suppose that for no t with 0 < t < ½ is there a rank-two
C = AB^* with Φ_{n,t}(C)B = 0 and Φ_{n,t}(C)^*A = 0. Then EP_n holds. Conversely, if EP_k holds for every k ≤ n, there
is no such C.

*Proof.* Let m(t) be the least value of q_n(t, C) over operators of rank at most two with ‖C‖_F = 1; it is continuous
in t. For t < 1/3 the one-copy map id − t Tr(·) I is positive definite on Hilbert–Schmidt space (eigenvalues 1 and
1 − 3t), so m(t) > 0. Suppose m(½) < 0, and let t_* be the least t with m(t) ≤ 0; then 1/3 ≤ t_* < ½ and m(t_*) = 0. A
minimiser C_* at t_* has rank two, because q_n(t, ab^*) = ⟨ā ⊗ b̄|(I − tF)^{⊗n}|ā ⊗ b̄⟩ ≥ (1 − t)^n ‖a‖²‖b‖². As a
minimum with value zero it is a critical point of q_n(t_*, ·) on X₂, which is the excluded system. For the converse,
write id − t Tr(·) I = (1 − 2t) id + 2t (id − ½ Tr(·) I) at each copy and expand:

    q_n(t, C) = Σ_S (1 − 2t)^{n−|S|} (2t)^{|S|} Σ_{α,β} q_{|S|}(C^S_{αβ}),

where C^S_{αβ} are the blocks of C with respect to a basis of the copies outside S; each block has rank at most two. By
EP_{|S|} every term is nonnegative, and the term S = ∅ is (1 − 2t)^n ‖C‖² > 0. So q_n(t, ·) > 0 on nonzero operators of
rank at most two for t < ½, and a zero, critical or not, cannot exist. ∎

The system in Proposition 22 contains no inequality: it asks that the residual Φ_{n,t}(C) have range orthogonal to the
range of C and co-range orthogonal to the co-range of C, and that this be impossible for rank-two C strictly inside the
Werner line. At one copy it is immediate. P_A Φ_{1,t}(C) = 0 reads C = t Tr(C) P_A with P_A the range projector, so
Tr C = 2t Tr C, and for t ≠ ½ this forces C = 0. At t = ½ the solutions include every zero of Section 6 and 𝓑_J.

The same expansion localises what has to be excluded. Assume EP_{n−1}. Then every term with S ≠ [n] is nonnegative,
so q_n(t, C) ≥ (1 − 2t)^n ‖C‖² + (2t)^n q_n(C). Hence a zero of q_n(t, ·) with t < ½ has q_n(C) < 0, and it is thick at
every copy, because thin-copy removal gives q_n(C) ≥ 0 whenever the range or the co-range has local dimension at most
two at some copy. So, given EP_{n−1}, EP_n is equivalent to: for 0 < t < ½, no critical point of q_n(t, ·) on X₂ with
value zero is thick at every copy.

**Proposition 23 (thin zeros are enough).** Call C ∈ X₂ thick if at every copy both its range and its co-range have
local dimension three, and thin otherwise. Let n ≥ 2. If every zero of q_n in X₂ is thin (Conjecture NT_n), then
EP_n holds.

*Proof.* The operators whose range lies in Q ⊗ (C^3)^{⊗(n−1)} for some plane Q at a fixed copy form a set of real
dimension 4 + 4(2·3^{n−1} − 2) + 4(3^n − 2) + 8 = 20·3^{n−1} − 4 (the plane Q, a range plane inside Q ⊗ (C^3)^{⊗(n−1)}, a
co-range plane, and the 2 × 2 core), against dim X₂ = 24·3^{n−1} − 8; the same holds for the co-range. So the thin
operators form a closed semialgebraic set of codimension at least 4·3^{n−1} − 4 ≥ 8, and the thick ones a connected
open dense set (as in the proof of Proposition 20(a)). By NT_n, q_n has no zero on it, so it has one sign there. That
sign is positive: ab^* + εcd^* is thick for generic vectors and small ε, and q_n(ab^*) > 0. Thin operators are limits
of thick ones, so q_n ≥ 0 on them by continuity. ∎

Equivalently, under NT_n the zero set lies in the thin locus, of codimension at least 8, and Proposition 20(a)
applies. An earlier version of the proposition also assumed EP_{n−1}, for the last step of the proof; an independent
audit pointed out that continuity makes it unnecessary. So EP_n follows from NT_n alone, separately for each n ≥ 2
(EP_1 is the one-copy inequality). Conjecture ZC_n implies
NT_n, because an operator in 𝓑_J is thin at the copies in J. NT_n asks only that zeros be thin somewhere, not that
they have the qubit form. Nothing in the argument excludes a thick operator at which q_n is positive but small.

**No margin uniform in n.** Let θ(C) be the least eigenvalue of the one-copy marginals of the range and co-range
projectors of C, minimised over copies; it vanishes exactly on thin operators and is at most 2/3. A bound
q_n(C) ≥ c θ(C) ‖C‖_F² with c independent of n is false. Take the one-copy zero P = |0⟩⟨0| + |1⟩⟨1|. Theorem 14′ gives
q_N(E_r P F_r^*) = (1/9)·2^{−3r} Σ_{ν≠0} q_1(P U_ν) = (5/3)·2^{−3r} on N = 1 + 8r copies (the eight values q_1(P U_ν) are
3/2 twice and 2 six times), while ‖E_r P F_r^*‖_F² = 2 and θ = 2/3, because every vector in either image has maximally
mixed marginals. So any such c_N is at most (5/4)·2^{−3(N−1)/8}. In floating point (`code/thickness_margin.py`) the
smallest ratios q_n/(θ‖C‖_F²) found are 0.295, 0.234 and 0.190 at n = 2, 3, 4 (at two copies this is close to
1 − 1/√2 ≈ 0.293, the constant `two-copy-thickness/` finds for its measure). The minimisers are strongly thick
(θ between 0.35 and 0.65), unbalanced and far from Hermitian, and they are not near the zero families. So NT_n is a
statement about exact zeros. Thick operators come exponentially close to zero, and a proof cannot go through a margin
that is uniform in n.

**Below the endpoint (floating point).** On the Werner line q_n(t, ·) of Proposition 22, random-start minimisation
reproduces to six digits the sharp curve conjectured in `contraction-notes/`,

    min over rank-two C of q_n(t, C)/‖C‖_F² = (1 − 2t)(1 − t)^{n−1}, attained at P ⊗ |a⟩⟨a|^{⊗(n−1)},

at n = 2 and 3 for t = 0.34, 0.4, 0.45 and 0.49, and at n = 4 for t = 0.45 (`code/below_endpoint/`, logs in
`logs/below_endpoint/`). Below the endpoint the landscape is simple. Every local minimum found (52 runs: n = 3 at
t = 0.40, 0.45, 0.49 and n = 4 at t = 0.45) is thin: it is supported on local qubits at the copies of some J and is a
product of rank-one projectors elsewhere, so it is a deformation of a family 𝓑_J. The values are ordered by |J|; at
n = 3, t = 0.45 they are 1, 1.818 and 2.488 times the sharp value for |J| = 1, 2, 3. Thick critical points located
directly (by driving the invariant residual of the critical-point equations to zero; values 0.17 to 0.64 at n = 3,
t = 0.45) are all saddles. By Proposition 22 and the remark after it, EP_n follows from EP_{n−1} together with:

  (TLM_n) for 1/3 < t < ½, q_n(t, ·) has no local minimum on rank-two operators that is thick at every copy.

**An elementary bound at t = 1/3 for every n (from an independent audit).** Let Q_d(X) = X − Tr(X) I/d. For every
operator C of rank at most r < d on n local d-dimensional factors, ‖Q_d^{⊗n} C‖_F² ≥ (1 − r/d)^n ‖C‖_F². As
superoperators on Hilbert–Schmidt space,

    Q_d = (1/d) Σ_{i<j, θ} |A_{ijθ}⟩⟨A_{ijθ}| + (1 − 1/d) Σ_{i≠j} |E_{ij}⟩⟨E_{ij}|,

with A_{ijθ} = uv^*, u, v = (e_i ± e^{iθ} e_j)/√2, θ ∈ {0, π/2, π, 3π/2}, and E_{ij} = |i⟩⟨j|: every atom has rank one
and trace zero. Contracting C with a rank-one atom at one copy, (u^* ⊗ I) C (v ⊗ I), keeps the rank at most r, so by
induction ‖Q_d^{⊗n} C‖² ≥ (1 − r/d)^{n−1} ‖(Q_d ⊗ I) C‖² = (1 − r/d)^{n−1}(‖C‖² − ‖Tr_1 C‖²/d) ≥ (1 − r/d)^n ‖C‖², using
‖Tr_1 C‖_F ≤ ‖C‖_1 ≤ √r ‖C‖_F. For qutrits of rank two this gives q_n(1/3, C) ≥ 3^{−n} ‖C‖_F² (the sharp curve
predicts ½ (2/3)^n). Each map id − t Tr(·) I has norm at most one for 1/3 ≤ t ≤ ½, so telescoping gives
q_n(t, C)/‖C‖_F² ≥ 3^{−n} − 3n(t − 1/3): positivity on an interval that shrinks with n, which proves nothing about NPT
bound entanglement but confirms that a first zero on the Werner line cannot sit at t = 1/3. We checked the frame
identity for d = 2, 3, 4. The same frame proof applied to the rank-sensitive bound of Section 9 needs a one-step
inequality of the kind that fails there.

The descent at thick critical points is analysed in `orbit-notes/`: averaging over the squeezes of one copy cannot
find it, and at the low critical points most of the curvature comes from coupling between copies. At t = 1/3 the sharp
curve would say ‖P_⊥^{⊗n} C‖_F² ≥ ½ (2/3)^n ‖C‖_F² for rank-two C, with P_⊥ the projection onto traceless operators. That
is a quantitative form of the fact that no rank-two operator is a sum Σ_s I_s ⊗ X_s.

**A census of zeros (floating point).** `code/zero_census.py` minimises q_n(AB^*)/‖AB^*‖_F² from random complex
factors (L-BFGS) and classifies the end point by the rank of the one-copy marginal of range plus co-range at each copy,
and by the operator Schmidt coefficients across every cut.

| n | runs | largest end value | on local qubit planes | rank-one factor ⊗ qubit zero (factor thick at its copies) | other |
|---|---|---|---|---|---|
| 2 | 12 | 4·10^{−16} | 12 | 0 | 0 |
| 3 | 40 | 8·10^{−13} | 36 | 4 | 0 |
| 4 | 24 | 5·10^{−13} | 15 | 9 | 0 |
| 5 | 16 | 2·10^{−12} | 5 | 11 | 0 |
| 6 | 8 | 4·10^{−14} | 0 | 8 | 0 |

Every run ended at a zero of the conjectured form (to rounding). One five-copy run stopped at 2·10^{−12} with best cut
ratio 1.4·10^{−4}, just above the threshold the script uses; continuing it from the same start
(`code/zero_census_continue.py`) brings it to 2·10^{−13} with ratio 4·10^{−5}, a rank-one factor on four copies times
a qubit projector on the fifth, and it is counted that way above. None ended at another local minimum, and none went
negative. As n grows the runs mostly end on the largest family (|J| = 1). EP_n itself is proved only for n ≤ 3 (the
three-copy case is [T1]); at four and five copies [T1] reaches Werner parameters just short of the endpoint (0.4979
and 0.47 against ½), and the four-copy endpoint is proved here only for Hermitian operators. So at n = 4, 5, 6 the
census is evidence for EP_n as well as for ZC_n. It is evidence about where the zeros are, not a proof. (An earlier
version of this paragraph said EP_n was known for n ≤ 5; an independent audit that read [T1] corrected this.)

## 9. The sharp rank-two bound for every rank-two operator (added 8 November 2026)

For positive C of rank at most two, Wu and Zou prove 2^n q_n(C) ≥ 2 Tr C² − (Tr C)² at three copies [WZ, Thm. 1.3].
`OBSTRUCTIONS.md` ("The sharp rank-two bound") records what is known about this bound for more copies; for odd n it
would give the endpoint through the projector reductions (Proposition 4 and `real-balanced-reduction/`). Here we
observe that the same bound appears to hold for every rank-two operator, positive or not, once Tr C is replaced by
the trace norm, and along the whole Werner line. In this section C has rank at most two, with singular values
σ_1 ≥ σ_2 and trace norm ‖C‖_1 = σ_1 + σ_2, and q_n(t, ·) is the Werner-line form of Proposition 22.

**Conjecture TN_n.** For 0 ≤ t ≤ ½ and every C of rank at most two on (C^3)^{⊗n},

    q_n(t, C) ≥ (1 − t)^{n−1} (‖C‖_F² − t ‖C‖_1²).                                          (9.1)

At the endpoint this reads 2^n q_n(C) ≥ 2‖C‖_F² − ‖C‖_1² = (σ_1 − σ_2)². For positive C it is the bound of Wu and
Zou, so TN_3 holds at t = ½ for positive C.

The general form matters for an induction on n: the slices and partial traces met there are not positive, and (9.1)
is tight on a large set. What is known:

- At n = 1 it is |Tr C| ≤ ‖C‖_1.
- At n = 2 it holds for every rank-two operator, in all local dimensions (Proposition 28 below).
- TN_n implies EP_n and the sharp curve of Section 8: for rank two, ‖C‖_1² ≤ 2‖C‖_F², so the right side of (9.1) is at
  least (1 − 2t)(1 − t)^{n−1} ‖C‖_F². At t = 1/3 it reads (3/2)^n ‖P_⊥^{⊗n} C‖_F² ≥ ‖C‖_F² − σ_1σ_2.
- It is attained for every ratio σ_2/σ_1. If X ≥ 0 has rank two and a is a unit vector, then at
  C = X ⊗ |a⟩⟨a|^{⊗(n−1)} (with X at any copy) both sides equal (1 − t)^{n−1}(‖X‖_F² − t (Tr X)²). At the endpoint the
  equality set is larger: for any set J of copies and λ_1, λ_2 ≥ 0, the operator
  (λ_1|0^J⟩⟨0^J| − (−1)^{|J|} λ_2|1^J⟩⟨1^J|) ⊗ |a⟩⟨a|^{⊗(n−|J|)} has 2^n q_n = (λ_1 − λ_2)². For |J| = 1 this is the
  product family, and λ_1 = λ_2 gives the signed GHZ zeros of Section 6. Below the endpoint only |J| = 1 attains the
  bound: for |J| = m ≥ 2 the difference of the two sides is 2λ_1λ_2 t (1 − t)^{n−m}[(1 − t)^{m−1} − t^{m−1}] > 0.
- At the endpoint TN_n says that every rank-two zero has σ_1 = σ_2. All zeros in 𝓑_J do: a rank-two antisymmetric
  matrix ab^T − ba^T has two equal nonzero singular values, and the unitary (J^{⊗|J|})^{−1}, the isometries V and the
  rank-one factor R_I keep them equal. So do all 100 zeros of the census (s_2/s_1 = 1 to four digits, n = 2, …, 6).
  This qualitative consequence of TN_n already implies EP_n:

**Corollary 25.** If every rank-two zero of q_n has σ_1 = σ_2, then EP_n holds. Equivalently (given EP_n), q_n(C) > 0
for every rank-two C with σ_1(C) ≠ σ_2(C).

*Proof.* The rank-two operators with σ_1 = σ_2 are the positive multiples of rank-two partial isometries. They are
fixed by a range plane and a co-range plane (real dimension 4(D − 2) each, D = 3^n), a unitary between them (4) and the
scale (1): 8D − 11 in all, against dim X₂ = 8D − 8. So they form a closed semialgebraic set of codimension 3, and
Proposition 20(a) applies. ∎

Like Corollary 21 and Proposition 23, this asks only that the zeros lie in a small set; here the set is defined by the
singular values alone, and every zero known lies in it.

*A form without singular values.* For C = Σ_i a_i b_i^* one has the crossed form
q_n(t, C) = Σ_{i,j} ⟨a_i ⊗ b_j|(I − tF)^{⊗n}|a_j ⊗ b_i⟩, where F swaps the two factors at each copy. With τ = t/(1 − t)
and P_A the antisymmetric projector, I − tF = (1 − t)(I + 2τP_A); put G = τ^{−1}((I + 2τP_A)^{⊗n} − I), which is
positive semidefinite. Write C = σ_1 u_1 v_1^* + σ_2 u_2 v_2^* with orthonormal pairs (u_1, u_2) and (v_1, v_2). Since
⟨u_1 ⊗ v_2, u_2 ⊗ v_1⟩ = 0, the difference of the two sides of (9.1) is

    (1 − t)^n τ [σ_1² G_11 + σ_2² G_22 + 2σ_1σ_2 (Re G_× + 1)],

with G_ii = ⟨u_i ⊗ v_i|G|u_i ⊗ v_i⟩ and G_× = ⟨u_1 ⊗ v_2|G|u_2 ⊗ v_1⟩. The ratio σ_2/σ_1 and the phase of u_2 are free,
so TN_n is equivalent to

    |⟨u_1 ⊗ v_2|G|u_2 ⊗ v_1⟩| ≤ 1 + (⟨u_1 ⊗ v_1|G|u_1 ⊗ v_1⟩ ⟨u_2 ⊗ v_2|G|u_2 ⊗ v_2⟩)^{1/2}           (9.2)

for all orthonormal pairs. At the endpoint (τ = 1), EP_n is the same inequality with the right side replaced by the
larger ((1 + G_11)(1 + G_22))^{1/2}. At n = 1, G = 2P_A and the left side of (9.2) is |⟨u_1, v_1⟩⟨v_2, u_2⟩| ≤ 1.

In distillability language: let M = (U ⊗ V)^*(2I − F)^{⊗n}(U ⊗ V) be the compression of the n-copy Werner operator
to a 2 ⊗ 2 subspace (U, V isometries from C² into (C^3)^{⊗n}). EP_n says that M^Γ ≥ 0 for every such M. TN_n at the
endpoint says more: ⟨s|M^Γ|s⟩ ≥ |s|² − Conc(s) for every s ∈ C² ⊗ C², where Conc(s) = |⟨s̄|σ_y ⊗ σ_y|s⟩| = 2|det S| is
the concurrence of s = vec(S). The one-copy compression to a qubit plane, M = 2I − F, attains it. When M^Γ is real in
the magic basis this is the spectral condition λ_1 + λ_2 ≥ 2 on its two smallest eigenvalues (together with λ_1 ≥ 0);
in general it is not: the smallest value of λ_1 + λ_2 over two-copy compressions is 1.967 (floating point).

**Proposition 26 (a product singular vector).** Let a = a_1 ⊗ ⋯ ⊗ a_n be a unit product vector, and let
C = σ_1 e^{iφ} |a⟩⟨a| + σ_2 u v^* with σ_1, σ_2 ≥ 0 and unit vectors u, v orthogonal to a. Then (9.1) holds for C at
every t with 0 ≤ t ≤ ½. These are the rank-two operators having a product vector as a common left and right singular
vector; they include operators that are thick at every copy.

*Proof.* Use (9.2) with u_1 = e^{iφ}a, v_1 = a, u_2 = u, v_2 = v. Since a ⊗ a is symmetric at every copy and
G = Σ_{K≠∅} 2^{|K|} τ^{|K|−1} P_A^K (each P_A^K antisymmetric at the copies in K), G_11 = 0, so we must show
|G_×| ≤ 1. At each copy take an orthonormal basis whose first vector is a_s, and split u = Σ_{K≠∅} u^{(K)}, where
u^{(K)} is the part in which exactly the copies in K carry basis vectors orthogonal to a_s (K = ∅ is absent because
u ⊥ a); split v the same way. At one copy, for basis vectors x (of u) and y (of v),

    ⟨a_s ⊗ y|x ⊗ a_s⟩ = ⟨a_s, x⟩⟨y, a_s⟩,   ⟨a_s ⊗ y|P_A|x ⊗ a_s⟩ = ½(⟨a_s, x⟩⟨y, a_s⟩ − ⟨y, x⟩).

The first is 1 if x = y = a_s and 0 otherwise; the second is −½⟨y, x⟩ if x and y are both orthogonal to a_s, and 0
otherwise. Hence ⟨a ⊗ v|P_A^K|u ⊗ a⟩ = (−½)^{|K|}⟨v^{(K)}, u^{(K)}⟩, and

    G_× = e^{−iφ} Σ_{K≠∅} (−1)^{|K|} τ^{|K|−1} ⟨v^{(K)}, u^{(K)}⟩.

For t ≤ ½ we have τ ≤ 1, so |G_×| ≤ Σ_K ‖v^{(K)}‖ ‖u^{(K)}‖ ≤ 1 by Cauchy–Schwarz. ∎

Equality needs τ^{|K|−1} = 1 on the support: below the endpoint u and v live on single excitations, as in the product
family; at the endpoint any excitation sizes are allowed, as in the signed GHZ family. The class contains thick
operators (for instance σ_1|0^n⟩⟨0^n| + σ_2 uv^* with generic u, v ⊥ |0^n⟩), so this is a genuine all-n case of TN_n,
not a reduction to fewer copies. For positivity alone it is weaker than Wu and Zou's product-ray theorem, which needs a
product vector only in the range; the point is the sharp quantitative form.

*Purified form.* For C = Σ_{i=1,2} a_i b_i^* put α = Σ_i |i⟩ ⊗ a_i and β = Σ_i |i⟩ ⊗ b_i in C² ⊗ (C^3)^{⊗n}, with
a reference qubit R, and let ρ, σ be the projectors onto α and β. Then ‖Tr_K C‖_F² = Tr(ρ_{RK} σ_{RK}), where the
subscript means the marginal on R and the copies in K, so

    q_n(t, C) = Σ_K (−t)^{|K|} Tr(ρ_{RK} σ_{RK}) = ⟨α ⊗ β|F_R ⊗ ⊗_s (I − tF_s)|α ⊗ β⟩,

where F_R and F_s exchange the reference qubits and the copies s of the two factors.

Moreover ‖C‖_F² = Tr(ρ_R σ_R), Tr C = ⟨β, α⟩, and ‖C‖_1 = F(ρ_R, σ_R) = Tr|ρ_R^{1/2} σ_R^{1/2}|: the singular values of
C are the square roots of the eigenvalues of G_A G_B, with G_A and G_B the Gram matrices of the a_i and of the b_i, and
ρ_R, σ_R are their transposes. In this language TN_1 is Uhlmann's theorem |⟨β, α⟩| ≤ F(ρ_R, σ_R). For n copies the
right side of (9.1) depends only on the two reference marginals, while the left side compares the two purifications
on every set of copies. With balanced factors a_i = √σ_i u_i, b_i = √σ_i v_i both reference marginals equal diag(σ_1, σ_2).
Put W = ⊗_s (I − tF_s), A = ⟨α ⊗ β|I_R ⊗ W|α ⊗ β⟩ and B = ⟨α ⊗ β|P_A^R ⊗ W|α ⊗ β⟩, the part on the reference singlet,
whose vector has squared norm σ_1σ_2. Then q_n = A − 2B, rank-one positivity gives A ≥ (1 − t)^n (σ_1 + σ_2)², and
TN_n is exactly

    B ≤ (1 − t)^{n−1} σ_1σ_2 + ½ (A − (1 − t)^n (σ_1 + σ_2)²)

(this formulation is from an independent audit; we checked the algebra). The excess of the symmetric part has to be
kept: bounding B on its own is the cross-term route below, which fails.

*An operator form.* Let Φ′_t(Z) = Tr(Z) I − tZ on one copy (the map whose n-th tensor power has Choi matrix
((I − tF)^{⊗n})^Γ) and R(Z) = Tr(Z) I − Z the reduction map, so that Φ′_t = (1 − t) Tr(·) I + tR. In the purified form,
⟨α ⊗ β|F_R ⊗ ⊗_s(I − tF_s)|α ⊗ β⟩ = ⟨β|(id_R ⊗ Φ′_t^{⊗n})(αα^*)|β⟩, and a factorisation that is not balanced only
makes ‖α‖²‖β‖² larger than ‖C‖_1². So TN_n is equivalent to the operator inequality

    (id_2 ⊗ Φ′_t^{⊗n})(ψψ^*) ⪰ (1 − t)^{n−1} (ρ_R − t I_2) ⊗ I   for every unit ψ ∈ C² ⊗ (C³)^{⊗n},

with ρ_R the marginal of ψψ^* on the reference qubit. Expanding the tensor power in the two terms of Φ′_t and dividing by
t(1 − t)^{n−1} gives, with τ = t/(1 − t) and ρ_{RK} the marginal of ψψ^* on R and the copies in K,

    Σ_{K≠∅} τ^{|K|−1} (id ⊗ R^{⊗K})(ρ_{RK}) ⊗ I_{K^c} + (I_2 − ρ_R) ⊗ I ⪰ 0.                                     (9.8)

At n = 1 this is I − ψψ^* ⪰ 0. When ρ_R = I_2/2 and t = ½ the right side above is 0 and the inequality is the endpoint
itself in operator form; compare the operator hierarchy H_k(ψ) ⪰ 0 of Bharti, Gajjala and Haug for pure states with a
maximally mixed qubit marginal [BGH]. (9.8) is the quantitative version for an arbitrary reference marginal. Its terms
with K ≠ ∅ are not positive one by one (R is not 2-positive), and peeling one copy off (9.8) leaves exactly the slice
term of the failed slice induction below. We checked the identity behind (9.8) numerically for n ≤ 3.

*On local qubits.* For C supported on a product of local qubit planes the endpoint form is ‖Π^{⊗n}C‖_F², with Π the
projection onto traceless 2 × 2 matrices, so 2^n q_n(C) = Σ_P |Tr(PC)|² over the full-weight Pauli strings
P ∈ {X, Y, Z}^{⊗n}. There TN_n says that the full-weight Pauli correlations of a rank-two operator are at least
(σ_1 − σ_2)². For a state ρ = p_1 ψ_1ψ_1^* + p_2 ψ_2ψ_2^* it reads ‖p_1 T_1 + p_2 T_2‖² ≥ (p_1 − p_2)², with T_i ∈ R^{3^n}
the full-weight correlation vectors of the two orthogonal pure states: every positive combination of (1, T_1) and
(−1, T_2) is spacelike for the form −x_0² + ‖x⃗‖². For one pure state this is ‖T‖² ≥ 1, the rank-one case (ψ ⊗ ψ has
no weight on odd antisymmetry patterns). At odd n the GHZ pair (|0^n⟩ ± |1^n⟩)/√2 has ‖T_i‖² = 2^{n−1} and
⟨T_1, T_2⟩ = −2^{n−1}, so ⟨T_1, T_2⟩ + 1 = −((‖T_1‖² − 1)(‖T_2‖² − 1))^{1/2}: the pair is extremal for (9.2), and its
equal mixture is the GHZ zero. Random-start minimisation on three, four and five qubits finds no violation
(`logs/trace_norm/tn_dims_2x2x2*_t0.5.log`).

*Positive operators.* For Hermitian C, q_n(t, C) = Tr[(C ⊗ C) ⊗_s (F_s − t)], where F_s swaps copy s of the first
factor with copy s of the second, and 2‖C‖_F² − (Tr C)² = Tr[(C ⊗ C)(2F − I)] with F = ⊗_s F_s. The F_s commute; on
the range of Π_T = ⊗_{s∈T} P_A^{(s)} ⊗_{s∉T} P_S^{(s)} each 2F_s − I is 1 or −3. Hence

    2^n q_n(C) − (2‖C‖_F² − (Tr C)²) = Σ_T c_{|T|} Tr[(C ⊗ C) Π_T],   c_k = (−3)^k − 2(−1)^k + 1,          (9.3)

with c_0 = c_1 = 0, c_2 = 8, c_3 = −24, c_4 = 80; this is the two-replica identity used in `OBSTRUCTIONS.md`, with its
coefficients written out. Every trace on the right is nonnegative when C ≥ 0, and then ‖C‖_1 = Tr C.

**Proposition 24.** For every positive semidefinite C on C^{d_1} ⊗ C^{d_2}, of any rank, and every t ≥ 0,

    q_2(t, C) = (1 − t)(‖C‖_F² − t (Tr C)²) + 4t Tr[(C ⊗ C)(P_A ⊗ P_A)] ≥ (1 − t)(‖C‖_F² − t (Tr C)²).

In particular TN_2 holds for positive semidefinite C (at every t, and in any local dimensions). This two-copy case
is elementary; we record it because the identity is exact for every t and every rank.

*Proof.* Expand: q_2(t, C) − (1 − t)(‖C‖_F² − t (Tr C)²) = t[‖C‖_F² + (Tr C)² − ‖Tr_1 C‖_F² − ‖Tr_2 C‖_F²], and the
bracket equals 4 Tr[(C ⊗ C)(P_A ⊗ P_A)] by P_A = (I − F)/2 at each copy. ∎

From three copies on, the rank matters even for positive operators: the n = 3 case of (9.3) is negative for some
positive C of rank three (the smallest value found is −(8/9)(Tr C)²), while for rank two none was found.

The two-copy case does not need positivity. The inequality (9.4) below, a refinement of the partial-trace inequality
of Fraser, Huber, Pozsgay and Vona [FHPV], is already in `first-order-bound/` of this repository (with a
sum-of-squares identity at r = 2, and the endpoint consequence 4q_2 ≥ (σ_1 − σ_2)²). An independent audit re-derived it
by keeping a variance term that the balanced-decomposition argument of [FHPV] drops; we give that proof in our words,
and add the consequences along the Werner line and for merging.

**Proposition 28 (two copies, every operator).** Let C act on C^{d_1} ⊗ C^{d_2} and have rank at most r. Then

    ‖Tr_1 C‖_F² + ‖Tr_2 C‖_F² ≤ ‖C‖_F² + (1 − 1/r)‖C‖_1² + |Tr C|²/r.                                    (9.4)

Both coefficients are sharp: P ⊗ |a⟩⟨b| (P a projector of rank r, a ⊥ b) and X ⊗ |a⟩⟨a| (X ≥ 0 of rank at most r)
give equality. For r = 2 and 0 ≤ t ≤ ½,

    q_2(t, C) − (1 − t)(‖C‖_F² − t‖C‖_1²) ≥ t(½ − t)(‖C‖_1² − |Tr C|²).                                  (9.5)

So TN_2 holds for every rank-two operator in all local dimensions. At the endpoint 4q_2(C) ≥ (σ_1 − σ_2)², and every
rank-two zero of q_2 has σ_1 = σ_2.

*Proof.* From a singular value decomposition write C = Σ_{i=1}^r x_i y_i^* with x_i = √σ_i u_i and y_i = √σ_i v_i.
Replacing the columns (x_i) and (y_i) by the same unitary combination does not change C and keeps the two Gram
matrices equal; call the common one G = (g_ij), with g_ij = ⟨x_i, x_j⟩ = ⟨y_i, y_j⟩. Choose the combination so that
y_i^* x_i = Tr C/r for every i. This is possible because every square matrix is unitarily similar to one with constant
diagonal: the mean of the eigenvalues lies in the numerical range, so some unit vector has that expectation, and one
repeats on its orthogonal complement. Put δ_i = g_ii. Then Σ_i δ_i = Tr G = ‖C‖_1 and
‖C‖_F² = Tr G² = Σ_i δ_i² + 2Σ_{i<j}|g_ij|². With Φ(Z) = (Tr_1 Z, Tr_2 Z) we use two facts.
- ‖Φ(uv^*)‖² ≤ ‖u‖²‖v‖² + |v^*u|². In the crossed form the difference is ⟨u ⊗ v|(I − F_1)(I − F_2)|u ⊗ v⟩ ≥ 0.
- ⟨Φ(uv^*), Φ(st^*)⟩ = ⟨Φ(tv^*), Φ(su^*)⟩. For every set K of sites, ⟨Tr_K(uv^*), Tr_K(st^*)⟩ and
  ⟨Tr_{K^c}(tv^*), Tr_{K^c}(su^*)⟩ are the same contraction of the four vectors, and at two sites the complement
  exchanges the two single sites.

The diagonal terms give ‖Φ(x_iy_i^*)‖² ≤ δ_i² + |Tr C|²/r². For i ≠ j, crossing and then Cauchy–Schwarz give
|⟨Φ(x_iy_i^*), Φ(x_jy_j^*)⟩| = |⟨Φ(y_jy_i^*), Φ(x_jx_i^*)⟩| ≤ δ_iδ_j + |g_ij|², by the first fact applied to y_jy_i^*
and to x_jx_i^*. Adding up,

    ‖Φ(C)‖² ≤ ‖C‖_1² + |Tr C|²/r + 2Σ_{i<j}|g_ij|² = ‖C‖_F² + (1 − 1/r)‖C‖_1² + |Tr C|²/r − V,

with V = Σ_i δ_i² − ‖C‖_1²/r ≥ 0. For (9.5), write N = ‖C‖_F², S = ‖C‖_1, T = |Tr C|² and P = ‖Φ(C)‖². Then
q_2(t, C) − (1 − t)(N − tS²) = t[N − P + tT + (1 − t)S²], and (9.4) at r = 2 says P ≤ N + S²/2 + T/2. ∎

At r = 2 the proof is an exact sum-of-squares identity: the sum-of-squares lemma of `two-copy-thickness/` writes the
slack of (9.4) as ½(δ_1 − δ_2)² plus the four contraction defects 4‖P_A^{(1)}P_A^{(2)}(a ⊗ b)‖² (for (a, b) = (x_1, y_1),
(x_2, y_2), (x_1, x_2), (y_1, y_2)) plus ‖Φ(x_2x_1^* − y_2y_1^*)‖². When σ_1 = σ_2, G is a multiple of I_2 in every
gauge, the first term vanishes, and all of the slack sits in the defects (an observation of an independent audit).
This matters because the all-copy problem at a fixed t reduces to rank-two orthogonal projectors, where σ_1 = σ_2
(`real-balanced-reduction/`).

[FHPV] prove ‖Tr_1 C‖_F² + ‖Tr_2 C‖_F² ≤ r‖C‖_F² + |Tr C|²/r for rank at most r; at r = 2 this is the two-copy
endpoint q_2 ≥ 0 in every local dimension. (9.4) is sharper because ‖C‖_1² ≤ r‖C‖_F², and at the endpoint with r = 2
it is exactly TN_2. This section had stated TN_2 only for positive operators (Proposition 24) and for real skew ones
(the rectangular commutator inequality below); the endpoint case for every operator was already in
`first-order-bound/`. (9.5) is also the two-site case of the merging budget below. `code/trace_norm/refined_partial_trace.py` checks (9.4) on 4000 random operators of
rank up to six, its two equality families, and (9.5) on random rank-two operators (log
`logs/trace_norm/refined_partial_trace.log`).

*More sites.* The n-site version of (9.4) at rank two, stated as a conjecture in `first-order-bound/` and proved there
for normal C,

    Σ_s ‖Tr_s C‖_F² ≤ (n − 1)‖C‖_F² + ½‖C‖_1² + ½|Tr C|²,                                               (9.4′)

holds for rank one (also in `first-order-bound/`): there it follows from Σ_s ‖Tr_s(uv^*)‖² ≤ (n − 1)‖u‖²‖v‖² + |v^*u|², since (n − 1)I + F − Σ_s F_s
has the eigenvalue 2|T| − 1 + (−1)^{|T|} ≥ 0 on the sector antisymmetric exactly at the copies in T. For rank two we
find no violation: random-start minimisation (`code/trace_norm/multisite_refined.py`, log
`logs/trace_norm/multisite_refined.log`, 10 runs each) on (3,3,3), (3,3,3,3), (2,2,2), (2,2,2,2), (2,2,2,2,2), (4,4,4)
and (2,3,4) ends at equality to rounding in every case. Equality holds on the product family X ⊗ |a⟩⟨a|^{⊗(n−1)} with X ≥ 0.
TN_n to first order in t at t = 0 is the weaker Σ_s ‖Tr_s C‖² ≤ (n − 1)‖C‖² + ‖C‖_1², so (9.4′) is a separate
conjecture, not a consequence of TN_n. The proof above does not extend to it as it stands: crossing sends a single
site to the complement of that site.

*Merging two copies (a sharp budget below the endpoint).* Allow arbitrary local dimensions, write
Δ_d(t, C) = q_d(t, C) − (1 − t)^{m−1}(‖C‖_F² − t‖C‖_1²) for the gap in (9.1) on m sites of dimensions d = (d_1, …, d_m),
and let d′ = (d_1 d_2, d_3, …, d_m) be the same system with the first two sites merged into one.

**Conjecture FB(t).** For every rank-two C and all local dimensions, Δ_d(t, C) ≥ ((1 − 2t)/2) Δ_{d′}(t, C).

Merging the sites one at a time gives Δ_{(3,…,3)} ≥ ((1 − 2t)/2)^{n−1} Δ_{(3^n)}, and on a single site
Δ_{(3^n)}(t, C) = t(‖C‖_1² − |Tr C|²) ≥ 0. So FB(t) for the dimension patterns (3^k, 3, …, 3) implies TN_n(t) for every n,
and FB(t) at a single t in (1/3, ½) would therefore give n-copy undistillability of the Werner state at t for every n,
that is, NPT bound entanglement. With two sites and nothing else, FB(t) is exactly the two-copy endpoint inequality
2‖Tr_1 C‖² + 2‖Tr_2 C‖² ≤ 2‖C‖² + ‖C‖_1² + |Tr C|² in those dimensions, for every t, so it holds by Proposition 28:
Δ_{(d_1,d_2)} − ((1 − 2t)/2)Δ_{(d_1d_2)} = (t/2)(2‖C‖_F² + ‖C‖_1² + |Tr C|² − 2‖Tr_1 C‖_F² − 2‖Tr_2 C‖_F²) ≥ 0. The
constant is sharp:
C = P ⊗ |0⟩⟨1| gives Δ_{(3,3)} = t(2 − 4t) and Δ_{(9)} = 4t. Random-start minimisation of Δ_d − γΔ_{d′}
(`code/trace_norm/fusion_budget.py`, log `logs/trace_norm/fusion_budget.log`) finds no violation for γ = (1 − 2t)/2 at
(3,3,3) → (9,3) (t = 0.4 and 0.45), (9,3,3) → (27,3) and (3,3,3,3) → (9,3,3) (t = 0.4), with equality reached at
nontrivial points, while γ = 0.11 at t = 0.4 already fails (−4.8·10^{−3}). On matrix-product inputs with the first two of
n copies merged (`code/trace_norm/fusion_budget_mps.py`, logs `logs/trace_norm/fusion_budget_mps_*.log`) there is no
violation at n = 6, 8 (t = 0.4) and n = 10 (t = 0.45) with bond dimension two, or at n = 6 with bond dimension three;
every run ends on the equality line Δ_d = ((1 − 2t)/2) Δ_{d′}. At those equality points the operator is a
rank-two operator on qubit planes of the two merged copies, for example λ_1|00⟩⟨00| − λ_2|11⟩⟨11|, times a rank-one
product on the remaining copies: the budget is tight exactly where the two-copy endpoint inequality is tight on the merged
pair. From three sites on FB(t) is open: the other sites enter through Φ_t^{⊗(m−2)}, and the crossing step in the
proof of Proposition 28 does not survive them (see the crossed form below).

**Lemma 27 (merging identity).** With ρ = id − Tr(·) I on a site, for every C and every t,

    q_d(t, C) = (1 − t) q_{d′}(t, C) + t ⟨C, ρ_1 ⊗ ρ_2 ⊗ Φ_t^{⊗(m−2)} C⟩.

*Proof.* On the two copies of sites 1 and 2, (I − tF_1)(I − tF_2) = (1 − t)(I − tF_1F_2) + 4t P_A^{(1)} P_A^{(2)}: both
sides are diagonal in the symmetric/antisymmetric splitting at each site, with values (1 − t)², (1 − t)(1 + t) and
(1 + t)². The first term is the crossed form of the merged site, the second that of ρ ⊗ ρ (whose natural representation
has partial transpose (I − F_1)(I − F_2)). ∎

So FB(t) is equivalent to ½Δ_{d′} + t⟨C, ρ_1 ⊗ ρ_2 ⊗ Φ_t^{⊗(m−2)} C⟩ ≥ 0: merging costs exactly the two-site term, which
can be negative. Merging sites 1, …, k into one block B_k (k = 1, …, n − 1) and iterating gives an exact form of the whole
problem,

    q_n(t, C) = (1 − t)^{n−1}(‖C‖_F² − t|Tr C|²) + t Σ_{k=1}^{n−1} (1 − t)^{k−1} ⟨C, ρ_{B_k} ⊗ ρ_{k+1} ⊗ Φ_t^{⊗(n−k−1)} C⟩,

so TN_n(t) says that this weighted sum of two-site terms is at least −(1 − t)^{n−1}(‖C‖_1² − |Tr C|²), the Uhlmann gap of
the purified form at the scale of the product floor. We checked the lemma on random operators for dimension patterns
(3,3), (3,3,3), (2,3,4), (3,3,3,3). The individual two-site terms can be negative even when q_n is not small: for
C = a_1|0^N⟩⟨0^N| + a_2|1^N⟩⟨1^N| the weighted sum of the two-site terms is 2a_1a_2[(1 − t)^{N−1} − (−t)^{N−1}] ≥ 0
(zero at the endpoint for odd N), while for N = 3 and a_1 = a_2 = 1 the first term is −2t and the second 2 (independent
audit; `real-balanced-reduction/code/below_endpoint_check.py`). Further random-start tests of FB at t = 0.4 with spectator and merged dimensions
(3,3,9), (2,2,4), (3,9,3), (2,4,4), (4,2,3), and at t = 0.34, 0.37, 0.49 for (3,3,3), found no violation; every run ended
on the equality line.

For real skew operators there is a second proof of the two-site case. If J = xyᵀ − yxᵀ with x, y real on R^{d_1} ⊗ R^{d_2}, then Tr J = 0,
‖J‖_1² = 2‖J‖_F² = 4h with h = ‖x‖²‖y‖² − ⟨x, y⟩², and the two partial traces are, up to sign, XYᵀ − YXᵀ and XᵀY − YᵀX
for the matrices X, Y of x, y. An independent audit proves the rectangular commutator inequality

    ‖XYᵀ − YXᵀ‖_F² + ‖XᵀY − YᵀX‖_F² ≤ 4(‖X‖_F²‖Y‖_F² − ⟨X, Y⟩_F²)

(remove the component of Y along X, use a singular-value decomposition of X, and bound each off-diagonal pair by
2(s_i + s_j)²(Y_{ij}² + Y_{ji}²) ≤ 4‖X‖_F²(Y_{ij}² + Y_{ji}²)); we checked it on random rectangular pairs. This is exactly
TN_2 at the endpoint for real skew J, in all dimensions, with equality at x = (|00⟩ + |11⟩)/√2, y = (|01⟩ − |10⟩)/√2.
Real skew operators matter because of Proposition S of `whole-problem-routes/`: a negative real rank-two C becomes a
negative real skew rank-two operator after r added flag copies. Its error terms carry the factors (1 − t)^r and (−t)^r
on the Werner line (we redid the flag computation at general t), so that step works for every t ≤ ½; the dilation,
realification and balancing steps of `real-balanced-reduction/` also work at every fixed t ≤ ½ (its section "Below the
endpoint"). With three or more sites
the pairwise commutator bounds alone do not give FB or TN for real skew J.

*The merging budget in the crossed form.* The following formulation and Proposition 29 are from the same independent
audit as Proposition 28; we checked both. Put γ = (1 − 2t)/2 and

    W_b = (I − tF_1)(I − tF_2) − γ(I − tF_1F_2) = ½(I − tF_1F_2) + 4tP_A^{(1)}P_A^{(2)},   B = W_b ⊗ ⊗_{s≥3}(I − tF_s).

W_b has the eigenvalues (1 − t)/2, (1 + t)/2 and (1 + 7t)/2 on the sectors where both, one or neither of the first two
sites is symmetric, so B ≥ λI with λ = (1 − t)^{m−1}/2. For C = σ_1u_1v_1^* + σ_2u_2v_2^* with orthonormal pairs,

    Δ_d − γΔ_{d′} = σ_1²D_1 + σ_2²D_2 + 2σ_1σ_2(Re c + η),

with D_i = ⟨u_i ⊗ v_i|B|u_i ⊗ v_i⟩ − λ ≥ 0, c = ⟨u_1 ⊗ v_2|B|u_2 ⊗ v_1⟩ and η = t(1 − t)^{m−2}/2. Since the ratio
σ_2/σ_1 and the phase of u_2 are free, FB(t) is equivalent to

    |c| ≤ η + (D_1D_2)^{1/2}   for all orthonormal pairs (u_1, u_2), (v_1, v_2),                         (9.6)

which has the shape of (9.2). The crossing that proves Proposition 28 does not extend. In Σ_K (−t)^{|K|}⟨Tr_K(·), Tr_K(·)⟩
crossing sends K to its complement, so on the sites beyond the first two the crossed weights are (−t)^{|R|−|K|}. They
form the metric ⊗_s (D_s − t·id) with D_s(X) = Tr_s(X) I, which is negative (−t) already on a single off-diagonal matrix
unit, so Cauchy–Schwarz is not available. This shows that the direct lift fails, not that FB fails.

**Proposition 29 (the merging budget on the class of Proposition 26).** Let a be a unit product vector, u and v unit
vectors orthogonal to a, σ_1, σ_2 ≥ 0 and C = σ_1e^{iφ}|a⟩⟨a| + σ_2uv^*. Then FB(t) holds for C for every number of
sites, all local dimensions and 0 ≤ t ≤ ½.

*Proof.* Take a = |0…0⟩ and use (9.6) with u_1 = e^{iφ}a, v_1 = a, u_2 = u, v_2 = v. Since a ⊗ a is symmetric at every
site, D_1 = 0, and we must show |c| ≤ η. At one site, ⟨0 ⊗ y|(αI + βF)|x ⊗ 0⟩ = αδ_{x0}δ_{y0} + βδ_{xy} for basis vectors
x, y, so c = e^{−iφ} Σ_w v̄_w u_w b(w), summed over basis words w, where b(w) is a
product of local factors: on the first two sites (1 − t)/2 if both digits are zero, −t/2 if exactly one is, t/2 if
neither is; on every other site 1 − t or −t according as the digit is zero or not. The word 0…0 does not occur, because
u ⊥ a. If the first two digits are zero, another digit is not, and |b(w)| ≤ ((1 − t)/2) t (1 − t)^{m−3} = η; otherwise
|b(w)| ≤ (t/2)(1 − t)^{m−2} = η (we used t ≤ 1 − t). Hence |c| ≤ η Σ_w |u_w||v_w| ≤ η. ∎

The constant is attained for every m: (λ_1|00⟩⟨00| − λ_2|11⟩⟨11|) ⊗ |0⟩⟨0|^{⊗(m−2)} has
Δ_d = 2tλ_1λ_2(1 − 2t)(1 − t)^{m−2} and Δ_{d′} = 4tλ_1λ_2(1 − t)^{m−2}. Proposition 29 strengthens Proposition 26 from
TN to the merging budget on the same class; it does not reach entangled singular vectors, where D_1 and D_2 are both
positive and correlated with c. `code/trace_norm/refined_partial_trace.py` checks the eigenvalue floor and the identity
for Δ_d − γΔ_{d′} (to 5·10^{−16}), Proposition 29 on 300 random operators with two to four sites, and the sharp example.

*A weaker comparison would do.* FB(t) asks for more than is needed. At one t in (1/3, ½) it is enough that

    q_d(t, C) ≥ γ q_{d′}(t, C)                                                                          (9.7)

for every rank-two C and the patterns d = (3^k, 3, …, 3): merging repeatedly gives
q_n(t, C) ≥ γ^{n−1} q_{(3^n)}(t, C) ≥ (1 − 2t)γ^{n−1}‖C‖_F² > 0, with q_{(3^n)} the form on one site of dimension 3^n. At t = 2/5, for example, this would make
ρ = (5I − 2F)/39 undistillable from any number of copies, while its partial transpose has the eigenvalue −1/39. With two
sites (9.7) holds with room, q_2 − γq_1 ≥ γ‖C‖_F² + (t/2)(σ_1 − σ_2)². Even the approximate budget
Δ_d − γΔ_{d′} ≥ −γ(1 − t)^{m−2}‖C‖_F² would do: it is equivalent to q_d ≥ γq_{d′} + (t/2)(1 − t)^{m−2}(σ_1 − σ_2)², which
implies (9.7). By Lemma 27, (9.7) says ½q_{d′} + t⟨C, ρ_1 ⊗ ρ_2 ⊗ Φ_t^{⊗(m−2)}C⟩ ≥ 0. Real C suffice: by
Theorem R_t of `real-balanced-reduction/`, a complex negative witness at qutrit copies gives a real one on more copies,
and merging keeps C real (independent audit). A proof for fixed m, for product
singular vectors, or for finitely many merged dimensions does not complete this route; it needs an estimate uniform in
the number of sites.

*What does not work.*
- Slice by slice. On local qubit planes the form is exactly self-similar: with M = CJ^{⊗n}, q_n(C) = ‖P_S^{⊗n} M‖_F²,
  and peeling the last copy gives q_n(C) = 3 E_w q_{n−1}(M_w (J^{−1})^{⊗(n−1)}) with M_w = (I ⊗ w^T) M (I ⊗ w) of rank at most
  two and w uniform on the unit sphere of C². Applying TN_{n−1} to every M_w would need the one-step inequality
  ‖M‖_1² ≥ ‖M_{01} − M_{10}‖_F² + 12 E_w σ_1σ_2(M_w), with M_{ab} the blocks of the last copy. It fails: the minimum
  is −0.0499 ‖M‖_F² when the other factor is C² and −0.0570 ‖M‖_F² when it is C^4 (deterministic quadrature on the
  Bloch sphere; the minimisers give the same values at 64 × 128 nodes). The slices of an extremal operator are not
  extremal together.
- Through the cross term of two orthogonal pure states. For C = λ_1 ψ_1ψ_1^* + λ_2 ψ_2ψ_2^* with λ_i ≥ 0, (9.1) at the
  endpoint would follow from q_n(ψψ^*) ≥ 2^{−n} (rank-one positivity) if also ⟨ψ_1ψ_1^*, Φ_n(ψ_2ψ_2^*)⟩ ≥ −2^{−n}.
  That holds at n = 2 (it is ‖Tr_1 X‖_F² + ‖Tr_2 X‖_F² ≤ 1 for X = ψ_1ψ_2^*, which follows from
  (I − F_1)(I − F_2) ≥ 0), but not at n = 3, where the cross term reaches −4·2^{−3}, nor at n = 4 (−7·2^{−4} and
  +5·2^{−4}). So the bound uses the excess of q_n(ψ_iψ_i^*) over 2^{−n} for entangled ψ_i. This is the simplest
  case of the failure of uniform Cauchy–Schwarz certificates recorded in `OBSTRUCTIONS.md`.
- An operator-norm shortcut on local qubits. There the cross term of (9.2) is ⟨v_2|Y_1|u_2⟩ with Y_1 the full-weight
  part of v_1u_1^* (scaled by 2^n), so (9.2) would follow from ‖P_{v_1^⊥} Y_1 P_{u_1^⊥}‖ ≤ Σ_P |⟨v_1|P|u_1⟩|² for every
  rank-one operator. This holds at n = 1 and fails from n = 2 on (excess 0.25, 0.73, 1.70 at n = 2, 3, 4). If one of the
  two rank-one parts is a product projector, Y is a product of reflections and (9.2) does hold.
- Merging two copies into one site is not monotone either: q_2(t, P ⊗ |0⟩⟨1|) = 2 − 4t, less than (1 − t) times the
  one-site value 2 of the same operator on C^9. Yet (9.1) holds there, with room for 0 < t < ½, because ‖C‖_1 = 2
  while Tr C = 0. The trace norm, not the trace, is what survives merging (compare the fusion obstruction for the map
  route found by an independent audit, in `OBSTRUCTIONS.md`).

**Evidence (floating point).** `code/trace_norm/tn_bound.py` minimises the difference of the two sides of (9.1),
divided by ‖C‖_F², from random complex factors (L-BFGS):

| copies | t | runs | smallest difference found |
|---|---|---|---|
| 2 | 0.40, 0.45, 0.50 | 12 each | −2·10^{−16} |
| 3 | 0.42, 0.50 | 10 each | 8·10^{−17} |
| 4 | 0.45 / 0.50 | 6 / 8 | 2·10^{−15} / 1·10^{−15} |
| 5 | 0.50 | 4 | 3·10^{−16} |
| local dimensions (4,4), (2,4), (5,5), (4,4,2) | 0.50 | 8 each | −2·10^{−16} |
| qubits: (2,2,2), (2,2,2,2), (2,2,2,2,2) | 0.50 | 10 each | −2·10^{−16} |

Every run ends at equality to rounding, mostly on the product family above (below the endpoint) or on zeros with
σ_1 = σ_2 (at the endpoint). None goes below rounding. On matrix-product inputs (`code/trace_norm/tn_mps.py`, exact transfer-matrix
evaluation, objective divided by σ_1σ_2 so that the search stays away from rank one; logs `logs/trace_norm/mps_*.log`)
every run ends either on the product family with σ_1 = σ_2 (t = 0.45; n = 8, 12, 16 at bond dimension two and
n = 8, 12 at bond dimension three) or at a zero with σ_1 = σ_2 (t = ½; the same sizes, 40 runs in all); the differences are below 10^{−9}‖C‖_F² and none is
negative. This is evidence, not a proof. TN_2 holds for every rank-two operator (Proposition 28); for operators
that are not positive, TN_n is open for every n ≥ 3.

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
- `code/below_endpoint/` (`sharp_curve.py n runs seed t1,t2,...`, `landscape.py n runs seed t`): Section 8, the Werner line
  below the endpoint (logs in `logs/below_endpoint/`).
- `code/zero_census.py`: Section 8, the census of zeros (`python zero_census.py n runs seed`; logs
  `logs/zero_census_n2.log` to `_n6.log`, run with seeds 1, 7, 11, 13, 17); `code/zero_census_continue.py` continues one
  start (log `logs/zero_census_n5_run5_continued.log`).
- `code/global_reductions/` (`run_all.sh`): Proposition 4′ (the identity and the projector choice, k = 1, 2,
  r = 1, 2, 3); Theorem 14′ (the identity by direct partial traces on 9 and 10 copies, the 81-state flag compression
  exactly for r = 1, 2, 3, marginals); the Remark after it (the flag block structure at r = 1, over all 3^8 column
  flag words); Theorem 14″ (the 729-state flag compression exactly at r = 2, the identity contracted with core data at
  k = 1, 2, isometry and marginals); Section 7 (the enumerator identities on random projectors at n = 2, 3, the
  three-copy code, and the certificates for n = 4, …, 9 in exact rational arithmetic); `check_enumerator_purity.py` (the
  purity-two enumerators of Section 7, exact); `perfect_code_check.py` (Proposition 19 attained at n = 5 by the
  five-qutrit code, and the enumerator formula for code subspaces); `stabilizer_anneal.py N t steps restarts seed`;
  `mds_code_values.py`; `rm_code.py t runs seed` (the Reed–Muller code of AG(2, 3)); `stabilizer_domination.py N t all|sample count seed`
  and `graph_cut_rank.py check|search N ...` (Conjecture SD_t); `graph_recursion.py` (the recursion (7.2)); `two_replica_balanced.py t d N ...` (the copy-symmetric
  two-replica relaxation on real balanced projectors at any t; it reproduces Proposition 19 at t = ½).
- `code/trace_norm/`: Section 9. `tn_bound.py t runs seed d_1 ... d_n` (the random-start test of (9.1); logs
  `logs/trace_norm/tn_*.log`), `check_identities.py` (crossed form, the Hermitian formula, (9.3), the bookkeeping behind
  (9.2), the odd-n equality family), `two_copy_psd.py` (Proposition 24 on random positive operators of every rank),
  `psd_rank.py n r runs` (rank three breaks the positive case at n = 3), `check_purified.py` (the purified form), `cross_term.py n runs seed ±1` (the cross
  term of two orthogonal pure states), `one_step_qubit.py D runs seed nth nph` (the failing one-step inequality on
  local qubits), `tn_mps.py n chi t runs seed mode` (matrix-product inputs, any n), `product_singular_vector.py
  n t runs` (the identity in the proof of Proposition 26, and the bound on that class).
- `code/trace_norm/fusion_budget.py t gamma runs seed d_1 ... d_m`: Conjecture FB (merging the first two sites). `fusion_budget_mps.py n chi t gamma runs seed`: the
  same on matrix-product inputs. `refined_partial_trace.py [seed]`: Propositions 28 and 29 and the crossed form of the
  merging budget. `multisite_refined.py runs seed d_1 ... d_n`: the n-site inequality (9.4′). `degree4_relaxation.py n t samples seed tn|q`: the
  two-replica relaxation of `OBSTRUCTIONS.md` (nonnegative at two copies, negative at three). `operator_form.py`: the identity behind
  the operator form (9.8).

## Local positivity at the signed GHZ zero, every n (added 8 November 2026)

Let C_n^0 = |0^n⟩⟨0^n| − (−1)^n |1^n⟩⟨1^n|, a zero of q_n for every n, and write rank-two operators as C = AB^* with A, B
of size 3^n × 2. An independent audit of this repository computes the real Hessian of F_n(A, B) = 2^n q_n(AB^*) at
(A, B) = ([|0^n⟩, |1^n⟩], [|0^n⟩, −(−1)^n |1^n⟩]) block by block for every n, and finds:
- it is positive semidefinite;
- its kernel has real dimension 4·2^n + 4n + 2;
- its smallest positive eigenvalue is 8 at n = 1 and 4 for n ≥ 2.

The kernel dimension equals the dimension of the "twisted-transpose" zero family K(u, v) = (uv^T − vu^T) J^{⊗n}
(J = [[0, −1], [1, 0]]) placed on arbitrary local qubit planes, which is 4·2^n − 6 + 4n, plus the 8 real gauge directions of
C = AB^*. So C_n^0 is a nondegenerate (Morse–Bott) minimum of the endpoint form on the rank-two variety, and q_n ≥ 0 in a
neighbourhood of it, for every n. The audit gives the explicit radius ‖C − C_n^0‖ ≤ 2^{−n}/200, inside which
q_n(C) ≥ 2^{−n−4} dist(C, zero family)²; we have not re-derived that estimate. We checked the three spectral statements by
computing the full Hessian numerically for n = 1, …, 4 (`code/ghz_zero_hessian.py`, `logs/ghz_zero_hessian.log`). This is a
local statement at one zero family and says nothing about thick operators.

## References

- [BCR] J. Bochnak, M. Coste, M.-F. Roy, Real Algebraic Geometry, Springer (1998).
- [BGH] K. Bharti, R. Gajjala, T. Haug, Two-copy nondistillability of Werner states: sharp partial-trace inequalities and
  finite-copy extensions, arXiv:2607.24479.
- [FHPV] T. C. Fraser, F. Huber, B. Pozsgay, I. Vona, On the two-copy distillability of Werner states and a new partial
  trace inequality, arXiv:2607.24309.
- [T1] T. S. Tomar, Undistillability of Werner states from three, four and five copies, Zenodo 23718035 (2026).
- [T2] T. S. Tomar, Equality cases and every-copy results for the Werner endpoint inequality, Zenodo 23718037 (2026).
- [WZ] T. Wu, Q. Zou, Sharp Plücker geometry for three-copy Werner distillation, arXiv:2608.02647.
