# Two orthogonal product sectors

Ansh Mishra and Aryan Senthilkumar, 2 November 2026

This note proves the endpoint inequality for every number of copies on a class of rank-two operators. Their range lies
in the sum of two product subspaces that are orthogonal at every copy; the co-range can be arbitrary (or the other way
round). Operators in this class can have one-copy reductions of rank three at every copy. The class contains the lowest
thick critical point we found at two copies. It is a small class, and the argument does not extend to general
operators.

Later sections free one or two copies entirely. They prove the inequality when the two sectors carry rank-one blocks
(Theorems A and A₂), and for every rank-two operator when the planes are spread out over the sector copies (Theorem C,
an exact reduction to at most three copies of larger dimension). These cover a four-copy class that is thick at every
copy, some open cases of Tomar's §8, and (with three sectors) the rank-two operators on C^3 ⊗ span{1…1, 2…2, 3…3} for
every number of copies. Collapsing blocks of copies one at a time (Theorem D) gives, for every n, classes that are thick
at every copy and need only the two-copy theorem. It is enough that one side (range or co-range) lies in such a product.

## Setting

For an operator C on (C^3)^⊗n,

    q_n(C) = Σ_{J ⊆ {1..n}} (−1/2)^{|J|} ‖Tr_J C‖²,

and n-copy undistillability of the qutrit Werner state at the endpoint is q_n(C) ≥ 0 for every C of rank at most two.

At every copy s, let A_s and B_s be nonzero subspaces of C^3 with A_s ⊥ B_s. Put A = A_1 ⊗ … ⊗ A_n and
B = B_1 ⊗ … ⊗ B_n.

**Theorem.** For every n, if C has rank at most two and both ran C and ran C* lie in A ⊕ B, then q_n(C) ≥ 0.

**Corollary (one side is enough).** For every n, if C has rank at most two and ran C lies in V = A ⊕ B, with the
co-range arbitrary, then

    q_n(C) ≥ q_n(P_V C P_V) ≥ 0.

The same holds with ran C* in V and the range arbitrary, since q_n(C*) = q_n(C). The first inequality holds for every
rank.

The main example is **one excitation**: A_s = K_s is a plane, B_s = C e_s with e_s ⊥ K_s, so A ⊕ B = K ⊕ C e. Without
the e-component this is the thin case (range in a product of planes; Tomar, Lemma 2.6). With it, the one-copy
reduction of the range generically has rank three at every copy.

How the class relates to known ones (write a and b for the numbers of copies at which A_s, resp. B_s, is a plane):
- **Product ray.** When a + b ≥ 2, a generic member has no product vector in its range or co-range. Then the
  product-ray theorem (Wu–Zou, Thm. 7.3, stated for three copies; its proof extends to every n) does not apply. When
  a + b ≤ 1, every two-dimensional range meets A or B, whose vectors are products.
- **Tomar's every-copy classes [T2, Thms. 4.1–4.7].** These are built on the two GHZ strings 0…0 and 1…1, with tails at
  single copies. Up to local unitaries, our two-sided class with a ≤ 1 and b ≤ 1 lies inside [T2, Thm. 4.1].
- **Few copies.** For n ≤ 3 the Theorem also follows from the known two- and three-copy results.
- **What is new.** The two-sided class with a ≥ 2 or b ≥ 2 (for example one excitation) at n ≥ 4, the one-sided
  Corollary, and the free-copy results below (Theorems A, A₂ and C), which are new at n ≥ 4.
- **An open case of [T2, §8].** Let T_0 and T_1 be disjoint sets of copies, and take
  A = ⊗_{s∈T_0} span{0,2} ⊗ ⊗_{s∉T_0} span{0} and B = ⊗_{s∈T_1} span{1,2} ⊗ ⊗_{s∉T_1} span{1}. These are orthogonal
  at every copy. So the Corollary covers every rank-two operator whose range lies in
  span{0…0, 2_i 0…0 (i ∈ T_0), 1…1, 2_j 1…1 (j ∈ T_1)} and whose co-range is arbitrary: the two GHZ backbones of
  [T2] with their tails on disjoint sets of copies. [T2, §8] lists range in V_0 + V_1 with an arbitrary co-range as
  open, except for operators whose range is spanned by two vectors of V_0 ∪ V_1 or whose co-range is orthogonal to
  all partner strings. What stays open in that case is the ranges in which a tail of the 0…0 backbone and a tail of
  the 1…1 backbone sit at the same copy, outside Tomar's conditions (partly covered by Theorem C below).
- **Credits.** The co-range splitting in the proof of the Corollary is the argument of [T2, Lemma 4.8]. The rank-one
  bound below is [Song–Chen, Prop. 8]. Wu–Zou [WZ, App. C.2] treat locally orthogonal product codes at three copies,
  which contains the case where every copy has (dim A_s, dim B_s) = (1,1).

## Proof

**Notation.**
- Every copy has (dim A_s, dim B_s) equal to (1,1), (2,1) or (1,2).
- Let a be the number of copies with dim A_s = 2, and b the number with dim B_s = 2. Then a + b ≤ n.
- The block P_A C P_A is a fixed rank-one projector at each copy where A_s is a line. So it is an operator X on
  (C^2)^⊗a.
  Likewise P_B C P_B is an operator Y on (C^2)^⊗b.
- Put P = P_A C P_B and Q = P_B C P_A.
- On qubits put Q_k(Z) = 2^{−k} Σ_σ |Tr σZ|², where σ runs over the 3^k Pauli strings with no identity factor, and
  Q_0(z) = |z|².

**Step 1 (exact identity).**

    2^n q_n(C) = 2^a Q_a(X) + 2^b Q_b(Y) + 2^n (‖P‖² + ‖Q‖²) + 2 (−1)^n Re(Tr X · conj Tr Y).

A partial trace over copy s pairs the copy-s factor of the range with that of the co-range.
- For P and Q these factors lie in A_s and B_s, which are orthogonal. So these blocks contribute only through J = ∅,
  where they are orthogonal to everything else.
- The partial traces of the two diagonal blocks live in orthogonal subspaces unless J = {1..n}. That gives the last
  term.
- For a diagonal block, a copy where its subspace is a line contributes the factor 1 − 1/2 = 1/2.
- On a qubit the one-copy map X ↦ X − ½ Tr(X) I kills I and fixes the Pauli matrices, so the copies of dimension two
  give exactly Q_a (resp. Q_b).

**Step 2 (two facts on qubits).**
- *Rank-one bound* [Song–Chen, Prop. 8]. 2^k Q_k(xy*) ≥ ‖x‖²‖y‖² ≥ |Tr xy*|². Indeed q_k(xy*) = ⟨x ⊗ y|(I − F/2)^⊗k|x ⊗ y⟩, and
  I − F/2 ≥ I/2.
- *Qubit lemma.* For k ≥ 1 and Z of rank two on (C^2)^⊗k,

      2^k Q_k(Z) ≥ (|Tr Z| − 2^k s₂(Z))₊²,

  where s₁ ≥ s₂ are the singular values.

  For k = 1: 2Q_1(Z) − (|Tr Z| − 2s₂)² = 2[s₁² − (|Tr Z| − s₂)²] whenever |Tr Z| > 2s₂, and |Tr Z| ≤ s₁ + s₂.

  For k ≥ 2: write Z = s₁u₁v₁* + s₂u₂v₂*. Then 2^k Q_k(Z) = ‖α + β‖² with α_σ = s₁⟨v₁|σ|u₁⟩ and
  β_σ = s₂⟨v₂|σ|u₂⟩.
  - The rank-one bound gives ‖α‖ ≥ s₁.
  - Completeness of all 4^k Pauli strings gives ‖β‖ ≤ 2^{k/2} s₂.
  - So ‖α + β‖ ≥ s₁ − 2^{k/2}s₂ ≥ |Tr Z| − 2^k s₂, using |Tr Z| ≤ s₁ + s₂ and 2^k − 1 ≥ 2^{k/2}.

**Step 3 (cases).** Put x = |Tr X| and y = |Tr Y|.

- *rank X ≤ 1 and rank Y ≤ 1.* By the rank-one bound, 2^n q_n ≥ x² + y² − 2xy ≥ 0.
- *rank X = 2.* (The case rank Y = 2 is symmetric.)
  - Then rank C = rank X, so Y = Q X⁺ P (Schur complement), and y = |Tr(X⁺PQ)| ≤ ‖P‖‖Q‖/s₂(X).
  - Hence 2^n(‖P‖² + ‖Q‖²) ≥ 2^{n+1}‖P‖‖Q‖ ≥ 2^{n+1−a} a′ y ≥ 2a′y, where a′ = 2^a s₂(X).
  - If rank Y ≤ 1, the qubit lemma for X and the rank-one bound for Y give

        2^n q_n ≥ (x − a′)₊² + y² + 2a′y − 2xy ≥ 0.

  - If rank Y = 2, then also X = P Y⁺ Q and x ≤ ‖P‖‖Q‖/s₂(Y). Averaging the two bounds, with b′ = 2^b s₂(Y) and
    n − a ≥ b ≥ 1, n − b ≥ a ≥ 1, the off-diagonal part is at least 2a′y + 2b′x. So

        2^n q_n ≥ (x − a′)₊² + (y − b′)₊² + 2a′y + 2b′x − 2xy ≥ 0.

    (If x ≥ a′ and y ≥ b′, use (x − a′)² + (y − b′)² ≥ 2(x − a′)(y − b′). Otherwise one of 2a′y, 2b′x is at least
    2xy.) ∎

**Proof of the Corollary.** Let ran C ⊆ V and put X_A = P_A C and X_B = P_B C.

At every copy, split the co-range factor of X_A into its parts in A_s and in A_s^⊥. This writes X_A = Σ_T X_A^(T),
where T is the set of copies at which the co-range factor lies in A_s^⊥. Three facts about these pieces:
- A partial trace over a copy in T vanishes on X_A^(T).
- Two pieces with T ≠ T′ have orthogonal partial traces, since they differ at a copy that is not traced. So
  q_n(X_A) = Σ_T q_n(X_A^(T)).
- Every piece has q_n ≥ 0, for any rank. At the copies outside T, its range and co-range lie in A_s, of dimension at
  most two, where the one-copy form is positive semidefinite. At the copies in T no trace survives.

Two pieces are special:
- The piece T = ∅ is P_A C P_A.
- The piece T = {1..n} has no surviving trace, so its q_n is its squared norm. That is at least ‖P_A C P_B‖², because
  B ⊆ ⊗_s A_s^⊥.

The same holds for X_B. The partial traces of X_A and X_B live in orthogonal ranges, so their only cross term is the
full trace. There only P_A C P_A and P_B C P_B contribute. Comparing with Step 1 for P_V C P_V gives

    q_n(C) ≥ q_n(P_V C P_V).

Since P_V C P_V has rank at most two, the Theorem gives q_n(P_V C P_V) ≥ 0. ∎

**Other local dimensions.** The argument uses d = 3 only through dim A_s + dim B_s ≤ 3. The form q_n is the same for
qudit Werner states at the endpoint, and an operator supported on a plane gives the same one-copy form in any dimension.
So the Theorem and the Corollary hold for every local dimension d ≥ 3 when every copy has
(dim A_s, dim B_s) ∈ {(1,1), (2,1), (1,2)}.

For d ≥ 4 a copy can also carry two orthogonal planes. Then the last case of Step 3 can have a = n, and the case split
above does not cover it. The numerics stay nonnegative there (d = 4, 5, n ≤ 3, adversarial minimisation), but we have
no proof.

## One free copy

Keep the sectors at copies 1, …, m (m = n − 1) and leave copy 0 unconstrained. So A = A_1 ⊗ … ⊗ A_m and
B = B_1 ⊗ … ⊗ B_m as before, V = C^3 ⊗ (A ⊕ B), and C = [[X, P], [Q, Y]] in blocks with respect to
C^3 ⊗ A ⊕ C^3 ⊗ B. Let a and b count the sector copies at which A_s, resp. B_s, is a plane. For a Pauli string μ on the
planar factors of A (identity on the one-dimensional ones) put X̂_μ = Tr_A[(1 ⊗ σ_μ) X], a 3 × 3 matrix, and
X̂_0 = Tr_A X; define Ŷ_ν and Ŷ_0 in the same way. Write

    q_1(Z) = ‖Z‖² − ½|Tr Z|²,   B_1(Z, W) = Tr(Z^* W) − ½ conj(Tr Z) Tr W,   q̃(P) = ‖P‖² − ½‖Tr_0 P‖²,

where Tr_0 is the trace over copy 0. The argument of Step 1 gives, with copy 0 carried along,

    2^m q_n(C) = Σ_{μ ∈ {1,2,3}^a} q_1(X̂_μ) + Σ_{ν ∈ {1,2,3}^b} q_1(Ŷ_ν) + 2^m (q̃(P) + q̃(Q)) + 2(−1)^m Re B_1(X̂_0, Ŷ_0).

The first sum is 2^m q_n(X). It is nonnegative when X has rank at most two. At every sector copy the range and
co-range of X lie in the plane or line A_s. On a plane K the separable decomposition of 1 − |ω_K⟩⟨ω_K| removes the copy,
where ω_K is the normalised maximally entangled vector of K ⊗ K̄ (the thin-copy reduction). On a line the one-copy form
is ½. Doing this at each sector copy leaves the one-copy inequality at copy 0. The same holds for Y. Also q̃(P) ≥ 0
when P has rank at most two, since q̃(P) is a sum of q_1 over 3 × 3 blocks of rank at most two. So only the last term
can be negative.

**Lemma L.** For square matrices Z and W of the same (arbitrary) size, |B_1(Z, W)|² ≤ ψ(Z) ψ(W), where
ψ(Z) = ‖Z‖₁² − ½|Tr Z|². The bound is attained at W = Z of rank one.

*Proof.* First let z and w have rank one. For every λ > 0 and every phase θ, the matrix λz − e^{iθ}w/λ has rank at most
two, so q_1 of it is nonnegative. This gives |B_1(z, w)|² ≤ q_1(z) q_1(w), and q_1 = ψ in rank one. In general, take
singular value decompositions Z = Σ_k s_k u_k v_k^* and W = Σ_l t_l u′_l v′_l^*. Expanding B_1 bilinearly and using the
rank-one case term by term gives |B_1(Z, W)| ≤ S_Z S_W with S_Z = Σ_k s_k (1 − ½|c_k|²)^{1/2}, c_k = ⟨v_k, u_k⟩. Put
p_k = s_k/‖Z‖₁. Concavity of the square root and Σ_k p_k|c_k|² ≥ |Σ_k p_k c_k|² give
S_Z² ≤ ‖Z‖₁²(1 − ½|Σ_k p_k c_k|²) = ψ(Z), since Σ_k s_k c_k = Tr Z. ∎

**Lemma R.** If X = xy^* has rank one, then Σ_μ q_1(X̂_μ) ≥ ψ(X̂_0).

*Proof.* Write F_0 and F_s for the swaps of the two replicas at copy 0 and at copy s, and F_A for the product of the
F_s over the planar copies of A. By the replica form ‖Tr_J xy^*‖² = ⟨x ⊗ y|F_J|x ⊗ y⟩,
Σ_μ q_1(X̂_μ) = ⟨x ⊗ y| W_0 ⊗ ⊗_s (2 − F_s) |x ⊗ y⟩, with the product over the planar copies of A and
W_0 = 1 − ½F_0 ≥ 0; one-dimensional copies contribute a factor 1. For f_s = ±1 we have Π_s(2 − f_s) ≥ 2 − Π_s f_s, so
⊗_s(2 − F_s) ≥ 2 − F_A. Moreover W_0 ⊗ (2 − F_A) = 1 − ½F_0 ⊗ F_A + 4P_0^− ⊗ P_A^− ≥ 1 − ½F_0 ⊗ F_A, where the P^− are
antisymmetric projectors. The expectation of the right side is ‖x‖²‖y‖² − ½|⟨y, x⟩|². Finally ‖X̂_0‖₁ ≤ ‖x‖‖y‖ and
Tr X̂_0 = ⟨y, x⟩. ∎

**Theorem A.** For every n ≥ 2, if ran C and ran C^* lie in V and rank X ≤ 1, rank Y ≤ 1, rank P ≤ 2 and rank Q ≤ 2,
then q_n(C) ≥ 0. Here C itself may have rank up to six.

*Proof.* By the identity, Lemma R (the X-sum is 0 if X = 0, and then so is X̂_0), q̃ ≥ 0 and Lemma L,

    2^m q_n(C) ≥ ψ(X̂_0) + ψ(Ŷ_0) − 2 (ψ(X̂_0) ψ(Ŷ_0))^{1/2} ≥ 0.  ∎

**Corollary B.** For every n ≥ 2, if C has rank at most two and its range is spanned by one vector of C^3 ⊗ A and one
vector of C^3 ⊗ B, then q_n(C) ≥ 0, whatever the co-range. The same holds with range and co-range exchanged, since
q_n(C^*) = q_n(C).

*Proof.* Follow the proof of the Corollary above, splitting the co-range only at the sector copies. Each piece is thin
at the sector copies outside T, so its form reduces to q_1 at copy 0 on pieces of rank at most two, which is
nonnegative. Two steps change. The piece T = {1..m} keeps the trace over copy 0, so its q_n is q̃ of it. Splitting its
co-range ⊗_s A_s^⊥ = B ⊕ R gives parts with orthogonal copy-0 traces, so its q_n is q̃(P_A C P_B) plus q̃ of a part of
rank at most two, hence at least q̃(P). The cross terms are J = {1..m} and J = {0..m}, and only P_A C P_A and P_B C P_B
contribute to them. This gives q_n(C) ≥ q_n(P_V C P_V). Unlike in the Corollary, this step needs rank C ≤ 2, since copy
0 carries the indefinite one-copy form. The blocks of P_V C P_V satisfy the hypotheses of Theorem A. ∎

The same argument covers the block-triangular members of rank at most two (P = 0 or Q = 0). There
rank X + rank Y ≤ 2, and if one of X and Y vanishes, so does the cross term.

**Scope.** At n = 2 this class is all of C^3 ⊗ C^3, since with one sector copy of dimensions (2, 1) we have
A_1 ⊕ B_1 = C^3. The general rank-two member therefore contains the full two-copy theorem, and we prove only the cases
above. For n ≤ 3 the rank-two cases (Corollary B and the block-triangular members) also follow from the known two- and
three-copy results, so for rank two the new content is at n ≥ 4. Theorem A at ranks three to six is not covered by
them. In the
setting of [T2, §8] with one shared tail copy, ranges spanned by one vector of each backbone space are already covered
there. The general case, with a diagonal block of rank two, is covered by Theorem C below when a, b ≤ 2 (using [T1]), and
is open for max(a, b) ≥ 3. For one excitation, Theorem E below covers ranges in C^3 ⊗ span{c, e}, and the rest reduces
to one explicit inequality (⋆). There the q̃ terms can vanish on projector-like off-diagonal blocks, so the
Schur-complement step of Step 3 has nothing to pay with.

## One free copy and one excitation

Here the sector copies carry one excitation: A = K_1 ⊗ … ⊗ K_m with planes K_s, and B = C e with e = e_1 ⊗ … ⊗ e_m,
e_s ⊥ K_s. This is the case a = m, b = 0 of "One free copy", and it is open in general for m ≥ 3. Treat K = A ⊕ C e
(dimension 2^m + 1) as one copy, and write q_{0,K} for the two-copy form on C^3 ⊗ K. Put
Σ_full = Σ_{μ∈{1,2,3}^m} q_1(X̂_μ), which is 2^m q_n(X), and Γ_X = 2 q_{0,K}(X) = 2 Σ_{α,β} q_1(X_{αβ}) − q_1(X̂_0), where the
X_{αβ} are the 3 × 3 blocks of X in a basis of A. Since b = 0, Y is a 3 × 3 matrix.

**Proposition M.** For ρ, τ > 0 let D_1 = 1_A ⊕ (−1)^{m+1}√(ρ/τ) 1_B, D_2 = 1_A ⊕ √(ρτ) 1_B and Z = (1 ⊗ D_1) C (1 ⊗ D_2).
Then for every operator C on C^3 ⊗ (A ⊕ C e),

    2^m q_n(C) − (2/ρ) q_{0,K}(Z) = Σ_full − Γ_X/ρ + (1 − ρ) q_1(Y) + (2^m − 2τ) q̃(P) + (2^m − 2/τ) q̃(Q).

*Proof.* The blocks of Z are X, √(ρτ) P, (−1)^{m+1}√(ρ/τ) Q and (−1)^{m+1} ρ Y. In q_{0,K}(Z) the trace over K gives
X̂_0 + (−1)^{m+1} ρ Y, and its q_1 contains the cross term 2(−1)^{m+1} ρ Re B_1(X̂_0, Y). Subtracting from the identity of
"One free copy" (where, with b = 0, the Y-sum is q_1(Y)) gives the formula. At ρ = τ = 1, Z = UC with U = P_A + (−1)^{m−1} P_B, and the formula is the one-group collapse plus the
term Σ_full − Γ_X, which appears because all m planar copies share the group (see "Two planes on one side cannot share a
group"); it vanishes for m = 1. ∎

Since Z has rank at most that of C, the two-copy theorem on C^3 ⊗ K gives, for C of rank at most two and all ρ, τ > 0,
that 2^m q_n(C) is at least the right-hand side.

**Lemma R′.** If X has rank at most two and ran X ⊆ C^3 ⊗ c for one vector c ∈ A (entangled or not), then
Σ_full ≥ Γ_X. The same holds if ran X* ⊆ C^3 ⊗ c.

*Proof.* Write X = (1 ⊗ |c⟩) Z with Z : C^3 ⊗ A → C^3 of rank at most two. Then X̂_μ = L(σ_μ c), where L(w) = Z(1 ⊗ |w⟩)
is a 3 × 3 matrix of rank at most two. So w ↦ q_1(L(w)) is a Hermitian form ⟨w|𝔐|w⟩ with 𝔐 ≥ 0. With weights
ω_μ = [μ has no identity factor] − 2^{1−m} + [μ = 0], we have Σ_full − Γ_X = Σ_μ ω_μ q_1(X̂_μ) = Tr(𝔐 Δ_c), where
Δ_c = Σ_μ ω_μ σ_μ |c⟩⟨c| σ_μ. Since Σ_{μ full} |⟨v|σ_μ|c⟩|² = ⟨c ⊗ v|⊗_s(2 − F_s)|c ⊗ v⟩ and Σ_μ |⟨v|σ_μ|c⟩|² = 2^m‖c‖²‖v‖²,
we get ⟨v|Δ_c|v⟩ = ⟨c ⊗ v|⊗_s(2 − F_s) − 2 + F_A|c ⊗ v⟩ ≥ 0 by the operator inequality in the proof of Lemma R. So
Δ_c ≥ 0 and Tr(𝔐 Δ_c) ≥ 0. The co-range version follows by taking adjoints. ∎

Without the hypothesis on the range the inequality can fail at rank two. The gap Σ_full − Γ_X reaches −2‖X‖², for
example at X = |0⟩⟨1| ⊗ 2^{−1/2} 1_{K_1} ⊗ |k⟩⟨k′| ⊗ |k⟩⟨k| ⊗ … with k ⊥ k′. It fails only on a small open set; random
rank-two X satisfy it.

**Theorem E.** For every m, every C of rank at most two whose range lies in C^3 ⊗ span{c, e} for some c ∈ A has
q_n(C) ≥ 0, whatever its co-range.

*Proof.* First let the co-range lie in V = C^3 ⊗ (A ⊕ C e). Then ran X ⊆ C^3 ⊗ c, so Lemma R′ gives Σ_full ≥ Γ_X.
Proposition M with ρ = τ = 1 then gives 2^m q_n(C) ≥ (Σ_full − Γ_X) + (2^m − 2)(q̃(P) + q̃(Q)) ≥ 0, because the
3 × 3 blocks of P and Q have rank at most two. For an arbitrary co-range, the proof of Corollary B gives
q_n(C) ≥ q_n(P_V C P_V), and P_V C P_V satisfies the hypotheses with co-range in V. ∎

For m ≤ 2 this is also in Theorem C (with [T1] when m = 2; Theorem E needs only the two-copy theorem). If c is a product
vector, span{c, e} is a two-sector space of lines and the one-sided Theorem D covers it. More generally, if c has
two-dimensional local support at two or fewer sector copies, the one-sided Theorem C covers the class (with [T1] when
there are two). So the new cases are those where c has two-dimensional local support at three or more sector copies.
Generic members are thick at every copy.

**The rest of the class reduces to one inequality.** In Proposition M, the bound is separable in ρ and τ. Since
Γ_X ≥ 0 (the two-copy theorem for X) and q_1(Y), q̃(P), q̃(Q) ≥ 0, taking the supremum over ρ, τ > 0 gives

    2^m q_n(C) ≥ Σ_full + q_1(Y) + 2^m (q̃(P) + q̃(Q)) − 2√(Γ_X q_1(Y)) − 4√(q̃(P) q̃(Q)).

Write s = √Σ_full, g = √Γ_X and y = √q_1(Y). If g ≤ s the right side is at least (s − y)² ≥ 0. If g > s it is at least
(y − s)² − 2y(g − s) + (2^{m+1} − 4)√(q̃(P) q̃(Q)). So the endpoint holds on the whole class (one free copy, one
excitation, every rank-two C with range and co-range in V) as soon as

    (⋆)   √q_1(Y) · (√Γ_X − √Σ_full) ≤ (2^m − 2) √(q̃(P) q̃(Q))   for every C of rank at most two in the class.

We have not proved (⋆). It holds when q_1(Y) = 0, and also when q̃(P) q̃(Q) = 0. For example, if q̃(Q) = 0 and Q ≠ 0,
every block of Q is a multiple of one rank-two projector. That forces the co-range of C into C^3 ⊗ span{w, e}, so Lemma
R′ applies to C*. If Q = 0, then rank X ≤ 1 whenever q_1(Y) > 0. So (⋆) is needed only where q_1(Y), q̃(P) and q̃(Q) are
all positive.

Numerically the largest ratio of the two sides is 4, attained, at m = 2, 3, 4, where the needed constants are 2, 6 and
14. An explicit operator attains 4 at every m (strings in orthonormal bases of the planes): C* = |R_1⟩⟨L_1| + |R_2⟩⟨L_2|
with R_1 = |0⟩|0…0⟩ + |1⟩|e⟩, L_1 = |1⟩|0…0⟩ + |0⟩|e⟩, R_2 = −|0⟩|110…0⟩ + |1⟩|e⟩ and L_2 = |1⟩|110…0⟩ + |0⟩|e⟩. It has
Σ_full = 0, Γ_X = 4, q_1(Y) = 4 and q̃(P) = q̃(Q) = 1. Adding a sector copy with a product factor leaves all five
quantities unchanged, so the supremum does not decrease with m. The maximisers found share a shape: Y has rank one and
trace zero, X has rank two with equal singular values, and Σ_full = 0, X̂_0 = 0, Γ_X = 2‖X‖². (⋆) with constant 4 would
settle every m ≥ 3. By the compression of Corollary B it would also cover range in V with any co-range. At m = 2 the needed constant 2 is too small, and that case stays with Theorem C. Fixed choices
ρ = τ = 1 do not suffice (the bound reaches −2‖C‖²), and a version of (⋆) with squares instead of square roots is false.

**Partial results on (⋆).** For every operator Z on C^3 ⊗ A, Σ_full(Z) ≤ 2^m‖Z‖², since q_1 ≤ ‖·‖² and the Pauli strings
on A satisfy Σ_μ ‖Ẑ_μ‖² = 2^m‖Z‖². If Z has rank at most two, also Γ(Z) ≤ 2‖Z‖². Indeed
2‖Z‖² − Γ(Z) = ‖Tr_0 Z‖² + ‖Tr_A Z‖² − ½|Tr Z|². Here Tr_A Z is 3 × 3, and Tr_0 Z has rank at most 6, since each rank-one
term contributes rank at most 3. So |Tr Z|² ≤ 3‖Tr_A Z‖² and |Tr Z|² ≤ 6‖Tr_0 Z‖², and 1/3 + 1/6 = 1/2. Both bounds are
attained at rank one, and the second fails at rank three for m ≥ 3.

*Lemma T.* Let X = X_0 + E, where every operator in span{X_0, E} has rank at most two and Γ(X_0) ≤ Σ_full(X_0). Then

    √Γ_X − √Σ_full(X) ≤ √Γ(E) + √Σ_full(E) ≤ (√2 + 2^{m/2}) ‖E‖.

*Proof.* On rank-two operators, Σ_full ≥ 0 (thin copies) and Γ ≥ 0 (two-copy theorem), so both forms are positive
semidefinite on the span and their square roots are seminorms there. So
√Γ(X) ≤ √Γ(X_0) + √Γ(E) ≤ √Σ_full(X_0) + √Γ(E) ≤ √Σ_full(X) + √Σ_full(E) + √Γ(E). ∎

*Theorem F (case 1a).* For m ≥ 4, q_n(C) ≥ 0 for every C of rank at most two with range and co-range in V such that the range and
the co-range each contain a nonzero vector of C^3 ⊗ A.

*Proof.* Rank one is covered by the rank-one bound, so let C = r_1 l_1* + r_2 l_2* have rank two. The condition says
that the e-components of r_1 and r_2 are linearly dependent (they may vanish), and likewise for l_1 and l_2. If either
pair vanishes, C is block-triangular with Y = 0, a case covered in "One free copy". Otherwise the e-components are
λ_k u and κ_k v, and Y = (Σ_k λ_k κ̄_k) uv*. If Y = 0, the left side of (⋆) vanishes. If Y ≠ 0, a change of basis
r ↦ rG, l ↦ lG^{−*} moves the e-components into r_1 and l_1 alone. This gives Y = uv*, P = x_1 v*, Q = u y_1* and
X = x_1 y_1* + x_2 y_2*, where x_k and y_k are the A-parts of r_k and l_k. Then q_1(Y) ≤ ‖u‖²‖v‖², q̃(P) ≥ ½‖x_1‖²‖v‖² and
q̃(Q) ≥ ½‖u‖²‖y_1‖². Lemma T with X_0 = x_2 y_2* (rank one, so Γ ≤ Σ_full by the operator inequality in Lemma R) and
E = x_1 y_1* gives (⋆) with constant 2√2 + 2^{1+m/2}. This is at most 2^m − 2 exactly when m ≥ 4. ∎

Theorem F uses only the two-copy theorem. It is not contained in Theorems A, C, D, E or Corollary B. Generic members
have rank X = 2. Neither their range nor their co-range contains a vector of C^3 ⊗ e or lies in C^3 ⊗ span{c, e}. They
are thick at every copy, and their range has two-copy reductions of rank 5 on every pair of sector copies. The new cases
are at n ≥ 5. A proof of √Γ_X − √Σ_full(X) ≤ 2σ_2(X) for rank-two X would give (⋆) with constant 4 in case 1a, since
‖x_1‖‖y_1‖ ≥ σ_2(X), and so would extend Theorem F to m = 3 (n = 4). This bound holds numerically with sharp constant 2 at
m = 2, 3, 4, attained at the maximiser shape.

For m ≥ 4, the remaining cases of (⋆) are those where Y has rank two (the e-components are independent on both sides),
and those where the e-components are parallel on one side only. At m = 3 (n = 4), case 1a also remains, since
2√2 + 2^{5/2} ≈ 8.49 exceeds 6. In the remaining cases P or Q can have 3 × 3 blocks of rank two, so q̃ can vanish on it
while it is nonzero, at the form c ⊗ Π. If the e-components are parallel on the range side, Q still has rank one and
q̃(Q) ≥ ½‖Q‖², so only P can degenerate; if Y has rank two, both can. Lemma T bounds the left side of (⋆) linearly in
‖E‖, which does not control these degenerations. `code/check_case_1a.py` (`logs/check_case_1a.log`) checks:
- the two bounds and Lemma T, on random two-piece X;
- in case 1a, the largest ratio in (⋆): 4.0 at m = 3 and 4, against the constants 8.49 and 10.83 that the proof gives and
  the 6 and 14 that are needed;
- adversarial minima of q_n over case 1a, which are 0 up to rounding. They are attained at the block-triangular corner
  Y, P, Q → 0. Away from it, with σ_2(X) ≥ σ_1(X)/2 and ‖Y‖, ‖P‖, ‖Q‖ ≥ 0.2‖C‖, the minimum of q_n/‖C‖² is about 0.18 at
  m = 4 (referee's check).

An independent referee re-derived these results with its own code at m = 1 to 4. It checked:
- the norm bounds (with rank-three controls);
- Lemma T adversarially (with rank-four controls);
- the gauge of Theorem F, including the case Y = 0;
- the constant;
- the novelty of the class.

It found no mathematical error, and its fixes are applied above.

Checks: `code/check_free_copy_excitation.py` (`logs/check_free_copy_excitation.log`, m = 3, 4) covers:
- Proposition M, to 4·10⁻¹⁵;
- Lemma R′, adversarial minimum ≥ −2·10⁻¹⁶ against a generic control at −2;
- Theorem E with arbitrary co-range, minimum ≥ 3·10⁻¹⁷ against a rank-three control at −0.5;
- the largest ratio in (⋆).

**Two free copies.** Everything above holds verbatim with two free copies F (H_F = C^9), with q_1 replaced by the
two-copy form q_F and K = A ⊕ C e as a third copy. Proposition M becomes an identity with q_F and q_{F∪{K}}. Lemma R′ holds
because w ↦ q_F(L(w)) is positive semidefinite by the two-copy theorem; the rest of its proof is unchanged. This gives:

**Theorem E₂.** For every m, every C of rank at most two whose range lies in C^9 ⊗ span{c, e} for some c ∈ A has
q_n(C) ≥ 0, whatever its co-range.

The main term uses the three-copy theorem in arbitrary dimensions [T1] on C^3 ⊗ C^3 ⊗ K, the q̃ terms use the two-copy
theorem, and the co-range compression is the one of Corollary B₂. With two free copies, Theorem C covers c with
two-dimensional local support at one sector copy at most. So the new cases are those where c has it at two or more.
Checked by `code/check_free_copies_excitation.py` (`logs/check_free_copies_excitation.log`, m = 2, 3, so n = 4, 5):
- Proposition M with q_F to 2·10⁻¹⁵;
- Lemma R′ with q_F, adversarial minimum ≥ 0;
- Theorem E₂ with arbitrary co-range, minimum ≥ 4·10⁻¹⁷ against a rank-three control at −0.5.

This paragraph has not had a separate referee.

**General two sectors.** The argument works for any two-sector structure, with planes on both sides. Let Σ_X and Σ_Y
be the sums of q_1(X̂_μ) and q_1(Ŷ_ν) over Pauli strings with no identity factor on the planar copies of A, respectively
of B. Put Γ_X = 2q_{0,A}(X) and Γ_Y = 2q_{0,B}(Y), with A and B each read as one copy, and K = A ⊕ B. The proof of
Proposition M gives

    2^m q_n(C) − (2/ρ) q_{0,K}(Z) = (Σ_X − Γ_X/ρ) + (Σ_Y − ρΓ_Y) + (2^m − 2τ) q̃(P) + (2^m − 2/τ) q̃(Q),

with D_1 = 1_A ⊕ (−1)^{m+1}√(ρ/τ) 1_B and D_2 = 1_A ⊕ √(ρτ) 1_B as before. (When B is a line this is Proposition M,
since then Σ_Y = Γ_Y = q_1(Y).) Lemma R′ holds on each side.

**Theorem E′.** For every two-sector structure, every C of rank at most two whose range lies in C^3 ⊗ span{c, c′} for some
c ∈ A and c′ ∈ B has q_n(C) ≥ 0, whatever its co-range.

*Proof.* As for Theorem E, with ρ = τ = 1: Lemma R′ gives Σ_X ≥ Γ_X and Σ_Y ≥ Γ_Y, the q̃ terms are nonnegative, and the
compression of Corollary B handles an arbitrary co-range. ∎

Theorem E′ complements Corollary B. There the range is spanned by one vector of C^3 ⊗ A and one of C^3 ⊗ B, so the
diagonal blocks have rank one. Here X and Y may have rank two, but their ranges are products between copy 0 and the
sector. With one free copy, Theorem C covers max(a, b) ≤ 2, so Theorem E′ is new when c or c′ has two-dimensional local
support at three or more sector copies.

Checked by `code/check_free_copy_general_sectors.py` (`logs/check_free_copy_general_sectors.log`) on four structures,
(a, b) = (2, 1), (2, 2), (1, 3) and (3, 0) at n = 4, 5: the identity holds to 6·10⁻¹⁵, and the adversarial minima are
≥ 2·10⁻¹⁶ against rank-three controls at −0.5. This paragraph has not had a separate referee.

An independent referee re-derived this section and checked it with its own code at m = 1 to 4. It checked Proposition M
against the one-free-copy identity and the sign of D_1, every step of Lemma R′ (with rank-three controls), the
compression for an arbitrary co-range piece by piece, and the reduction to (⋆), and it found the explicit maximiser
above. It also found no gap at the edge cases. It found no mathematical error, and its fixes are applied above.

## Two free copies

The rank-one argument does not need Lemma L. It needs only the endpoint inequality for one more copy, where that copy
may have any dimension. This contains the f-copy inequality used below, since q_{F ∪ {K}}(Z ⊗ |e⟩⟨e|) = ½ q_F(Z). This
lets us free a second copy.

Let F be a set of f free copies and keep the sector copies 1, …, m as above, so n = f + m. Write H_F = (C^3)^{⊗f},
q_F for the f-copy form on operators on H_F, B_F(Z, W) = Σ_{J ⊆ F} (−½)^{|J|} Tr[(Tr_J Z)^* Tr_J W] for its sesquilinear
form, and q̃_F(P) = Σ_{J ⊆ F} (−½)^{|J|} ‖Tr_J P‖². The identity of the previous section holds verbatim with q_1, B_1, q̃
replaced by q_F, B_F, q̃_F; now the X̂_μ are operators on H_F. The proof of Lemma R also goes through as far as the
operator inequality, since W_F = ⊗_{c∈F}(1 − ½F_c) ≥ 0. Write x = Σ_α x_α ⊗ α and y = Σ_β y_β ⊗ β over an orthonormal
basis {α} of A (likewise x′, y′ over an orthonormal basis of B). Then

    Σ_μ q_F(X̂_μ) ≥ ⟨x ⊗ y| W_F ⊗ (2 − F_A) |x ⊗ y⟩ = 2 Σ_{α,β} q_F(x_α y_β^*) − q_F(X̂_0).

**Theorem A₂.** For every n ≥ 3, let F be any two copies, and let the other m = n − 2 copies be sector copies as above.
Suppose ran C and ran C^* lie in H_F ⊗ (A ⊕ B). Take blocks with respect to H_F ⊗ A ⊕ H_F ⊗ B. If rank X ≤ 1,
rank Y ≤ 1, rank P ≤ 2 and rank Q ≤ 2, then q_n(C) ≥ 0. Here C may have rank up to six.

*Proof.* Write X = xy^* and Y = x′y′^* (the case X = 0 or Y = 0 is covered by putting x = 0 or x′ = 0). Bound the
X-sum and the Y-sum by the inequality above, and put e^{iφ} = (−1)^{m+1}:

    2^m q_n(C) ≥ 2 q_{F ∪ {K}}(𝒞) + 2^m (q̃_F(P) + q̃_F(Q)),   𝒞 = xy^* ⊕ e^{iφ} x′y′^*.

Here 𝒞 is read as an operator on H_F ⊗ K, where K = A ⊕ B counts as one more copy with the endpoint form there. Indeed

    q_{F ∪ {K}}(𝒞) = Σ_{α,β} q_F(x_α y_β^*) + Σ_{α′,β′} q_F(x′_{α′} y′_{β′}^*) − ½ q_F(X̂_0 + e^{iφ} Ŷ_0).

With this φ, q_F(X̂_0 + e^{iφ} Ŷ_0) = q_F(X̂_0) + q_F(Ŷ_0) − 2(−1)^m Re B_F(X̂_0, Ŷ_0). So 2 q_{F ∪ {K}}(𝒞) is exactly the
bound for the X-sum, the Y-sum and the cross term together.

The operator 𝒞 has rank at most two. So q_{F ∪ {K}}(𝒞) ≥ 0 by the three-copy endpoint theorem with factors C^3, C^3
and K, which holds for arbitrary factor dimensions [T1, Thm. 1.2]. Each 9 × 9 block of P has rank at most two, so
q̃_F(P) ≥ 0 by the two-copy theorem; likewise for Q. ∎

With f = 1 the same proof uses the two-copy theorem for C^3 ⊗ K, which holds in arbitrary dimensions [FHPV, BGH].
This gives a second proof of Theorem A. In general, the analogue of Theorem A with f free copies follows from the
endpoint inequality for f + 1 copies, one of them of arbitrary dimension. If a cited theorem is stated for one common
local dimension, unequal factors follow by embedding each factor isometrically into the largest one; q_n and the rank
are unchanged.

**Corollary B₂.** For every n ≥ 3 with two free copies: if C has rank at most two and its range is spanned by one vector
of H_F ⊗ A and one vector of H_F ⊗ B, then q_n(C) ≥ 0, whatever the co-range.

*Proof.* As for Corollary B, with the co-range split only at the sector copies. Here P_A C = ug^* and P_B C = vh^* have
rank at most one, so every piece reduces at the sector copies to q_F of rank-one operators, which is nonnegative. The
piece T = {1..m} has q_n equal to q̃_F of it; splitting its co-range ⊗_s A_s^⊥ = B ⊕ R gives q̃_F(P_A C P_B) plus a
nonnegative term. The cross terms are J_F ∪ {1..m}, J_F ⊆ F, and only P_A C P_A and P_B C P_B contribute. So
q_n(C) ≥ q_n(P_V C P_V); this step holds for every rank-two C with ran C ⊆ V, by the two-copy theorem. The blocks of
P_V C P_V have rank at most one, so Theorem A₂ applies. ∎

**Status.** Theorem A₂ and Corollary B₂ rest on Tomar's computer-assisted three-copy theorem [T1]. For rank two they
are new at n ≥ 4. At four copies, [T1] reaches only α ≤ 0.4979. The Hermitian members, and those with range and co-range
in one plane, are also covered by the certificates in `four-replica-notes/`. Theorem A₂ at ranks three to six is not a
special case of the known rank-two results, even at n = 3, where it is derived from [T1].

## Grouping the sector copies

Without a rank condition on the blocks there is still an exact reduction, as long as the planes are spread out. Split
the sector copies into g groups G_1, …, G_g so that each group contains at most one copy with a planar A_s and at most
one copy with a planar B_s. Then A_i = ⊗_{s∈G_i} A_s and B_i = ⊗_{s∈G_i} B_s have dimension at most two. Treat
K_i = A_i ⊕ B_i (of dimension at most four) as one copy, with the endpoint form there. Here K_i ⊆ ⊗_{s∈G_i} C^3 and
A ⊕ B ⊆ K_1 ⊗ … ⊗ K_g. An operator on H_F ⊗ (A ⊕ B) is read as an operator on H_F ⊗ K_1 ⊗ … ⊗ K_g, and Tr_{K_i} is
the trace over the copies in G_i.

Write X = P_A C P_A and Y = P_B C P_B, and let X̂_0 and Ŷ_0 be their traces over all sector copies (operators on H_F).
In a group the one-copy forms multiply out to 2^{1−|G_i|} times the endpoint form on A_i. A line contributes ½, and
A_i is isometric to its planar factor if it has one, and is a line (form ½) otherwise. The same holds for B_i. The
off-diagonal blocks have no surviving trace over a sector copy or over a K_i. The two diagonal blocks meet only in the
terms that trace out every sector copy (every K_i), with any J ⊆ F. As in Step 1,

    q_n(C) = q_n(X) + q_n(Y) + q̃_F(P) + q̃_F(Q) + 2(−½)^m Re B_F(X̂_0, Ŷ_0),
    q_{F ∪ {K_1, …, K_g}}(UC) = 2^{m−g} (q_n(X) + q_n(Y)) + q̃_F(P) + q̃_F(Q) + 2(−½)^g (−1)^{m−g} Re B_F(X̂_0, Ŷ_0),

with U = P_A + (−1)^{m−g} P_B. Since 2^{g−m}(−½)^g(−1)^{m−g} = (−½)^m, this gives, for every operator C on
H_F ⊗ (A ⊕ B),

    q_n(C) = 2^{g−m} q_{F ∪ {K_1, …, K_g}}(UC) + (1 − 2^{g−m}) (q̃_F(P) + q̃_F(Q)).

**Theorem C.** Let a and b count the sector copies with planar A_s, resp. planar B_s. If f + max(a, b, 1) ≤ 3, then
q_n(C) ≥ 0 for every n and every C of rank at most two with ran C and ran C^* in H_F ⊗ (A ⊕ B). (For
f + max(a, b, 1) = 3 this uses [T1].)

*Proof.* Pair each planar-A copy with a different planar-B copy. Together with the unpaired planar copies, this gives
g = max(a, b, 1) groups of the required kind (if a = b = 0, one group of lines); the lines can go anywhere. UC has rank
at most two. So the first term is nonnegative by the endpoint theorem for f + g ≤ 3 copies, g of them of dimension at
most four, which holds in arbitrary dimensions [T1, Thm. 1.2]; for f + g ≤ 2 the two-copy theorem suffices [FHPV, BGH].
Next, q̃_F(P) = Σ_{α,β} q_F(⟨α|P|β⟩) over orthonormal bases of A and B, and each ⟨α|P|β⟩ has rank at most two. So
q̃_F(P) ≥ 0 (trivially for f ≤ 1, by the two-copy theorem for f = 2), and likewise q̃_F(Q) ≥ 0. Finally
1 − 2^{g−m} ≥ 0, since g ≤ m. ∎

No rank condition on X, Y, P or Q is needed.

**One-sided version.** If f + max(a, b, 1) ≤ 3, C has rank at most two and ran C ⊆ V = H_F ⊗ (A ⊕ B), then
q_n(C) ≥ q_n(P_V C P_V) ≥ 0 for any co-range. The first inequality is the proof of Corollary B; its pieces are
nonnegative because q_F ≥ 0 on rank-two operators on H_F. The second is Theorem C.

The cases covered are:
- f = 2 with at most one plane on each side;
- f = 1 with at most two planes on each side;
- f = 0 with up to three, which the Theorem above already covers for every a and b.

**Example (four copies).** Take F = {1, 2}. At copy 3 put A_3 = L_3, a plane, and B_3 = C e_3 with e_3 ⊥ L_3. At copy 4
put A_4 = C e_4 and B_4 = L_4, a plane orthogonal to e_4. Then a = b = 1 and g = 1, and [T1] is used with factors C^3,
C^3 and K_1 ≅ C^4. Every rank-two C with range in C^9 ⊗ (L_3 ⊗ e_4 ⊕ e_3 ⊗ L_4), with any co-range, has
q_4(C) ≥ 0. For generic members the one-copy reductions of range and co-range have rank three at every copy.

**An open case of [T2, §8] with shared tails.** For the GHZ backbones of [T2, §8], with tails T_0 and T_1, take every
copy in T_0 ∩ T_1 as a free copy. At the copies in T_0 \ T_1 put A_s = span{0, 2} and B_s = span{1}. At the copies in
T_1 \ T_0 put A_s = span{0} and B_s = span{1, 2}. At all other copies put A_s = span{0} and B_s = span{1}, of type
(1, 1). Every backbone string then lies in H_F ⊗ A or in H_F ⊗ B. So the one-sided version covers every rank-two
operator whose range lies in span{0…0, 2_i 0…0 (i ∈ T_0), 1…1, 2_j 1…1 (j ∈ T_1)}, with arbitrary co-range, whenever

    |T_0 ∩ T_1| + max(|T_0 \ T_1|, |T_1 \ T_0|, 1) ≤ 3.

This includes tails of both backbones at the same copy, the case that the Corollary leaves open. Examples:
- one shared tail copy, with up to two further tails on each backbone;
- two shared tail copies, with up to one further tail on each.

**Two planes on one side cannot share a group.** For rank-one X the previous section gives
2^m q_n(X) ≥ 2 q_{F ∪ {A}}(X) for any number of planar factors, but for rank two this fails. With planar A_s and A_t in
one group, 2^m q_n(X) − 2 q_{F ∪ {A}}(X) is twice the form with coefficient −½ at the free copies and −1 at s and t.
For a ⊥ a′ in A_t, the operator X = |0…0⟩⟨1…1|_F ⊗ (1/√2) 1_{A_s} ⊗ |a⟩⟨a′| has rank two and gives −2‖X‖², the minimum
over all X. This limits the method, not the inequality: in these cases q_n itself stays nonnegative in our searches.

**Status.** Theorem C rests on the two-copy theorem and on Tomar's computer-assisted three-copy theorem [T1]. For
n ≤ 3 it is contained in [T1] and the two-copy theorem; it is new at n ≥ 4. Tomar merges copies in other ways (a
linear program over grouped partial-trace norms in [T1], and a merging step inside the proof of [T2, Thm. 4.6]); we
did not find this identity there.

Checks:
- `code/check_one_plane_each.py` checks one group, `code/check_grouped_sectors.py` two groups, and
  `code/check_tomar8_shared_tails.py` the [T2, §8] shapes; the logs are `logs/check_one_plane_each.log` (which also has
  the §8 runs) and `logs/check_grouped_sectors.log`.
- The identity holds to 2.5·10⁻¹⁵ for random C of any rank.
- Adversarial minima of q_n/‖C‖² over rank-two members are ≥ −2.6·10⁻¹⁶ at n = 3, 4, 5 (2 to 12 L-BFGS starts per
  configuration; the counts are in the logs).
- For the [T2, §8] shapes with shared tails at n = 4 and arbitrary co-range, including {0,1,2}/{0,3} (two groups) and
  {0,1,2}/{0,1,3}, the minima are ≥ −1.6·10⁻¹⁶.
- `code/check_group_gap.py` (`logs/check_group_gap.log`, f = 0 and 1) checks the closed form of the gap, its smallest
  eigenvalue −2 and the rank-two example above; over rank one the gap is ≥ −7·10⁻¹⁶.

**More than two sectors.** The same bookkeeping works for r sectors S_i = ⊗_s V_{i,s} (i = 1, …, r), with
V_{1,s}, …, V_{r,s} mutually orthogonal at every sector copy and grouped so that each sector has at most one planar
factor in each group, as long as there are phases ε_i with ε_i conj(ε_j) = (−1)^{m−g} for all i ≠ j. Then the identity
above holds with U = Σ_i ε_i P_{S_i}, and with Σ_{i≠j} q̃_F(X_ij) in place of q̃_F(P) + q̃_F(Q), where X_ij is the block
of C from H_F ⊗ S_j to H_F ⊗ S_i. For r = 2 such phases always exist (this is U above). For r = 3 every V_{i,s} is a
line. The conditions on the three pairs force ε_1 = ε_2 = ε_3 and (−1)^{m−g} = 1, so the phases exist exactly when
m − g is even, and then U = 1. So if m − g is even, then for every operator C on H_F ⊗ (S_1 ⊕ S_2 ⊕ S_3),

    q_n(C) = 2^{g−m} q_{F ∪ {K_1, …, K_g}}(C) + (1 − 2^{g−m}) Σ_{i ≠ j} q_F(X_ij),

where X_ij is read as an operator on H_F (S_i and S_j are lines) and each K_k ≅ C^3 is spanned by the three sector
strings on G_k. Take g = 1 for odd m and g = 2 for even m. So for every m and every C of rank at most two with range
and co-range in C^3 ⊗ span{1…1, 2…2, 3…3}, q_{m+1}(C) ≥ 0. Here i…i is the product of the i-th vectors of an
orthonormal basis chosen at each sector copy, the same basis for range and co-range. The merged operator is C itself,
so its form is nonnegative by the two-copy theorem when m is odd and by [T1] at d = 3 when m is even. Each X_ij has rank
at most two, so q_1(X_ij) ≥ 0, and 1 − 2^{g−m} ≥ 0. For odd m the same argument with two free copies gives
q_{m+2}(C) ≥ 0 for range and co-range in C^9 ⊗ span{1…1, 2…2, 3…3}. It uses [T1] for the merged operator and the
two-copy theorem for the blocks X_ij.

Generic members are thick at every copy, and Theorem C covers them only when f + m ≤ 3. For n ≤ 3 these statements
follow from the two- and three-copy theorems, so they are new from four copies on. They are also outside
[T2, Thms. 4.1–4.7], which need the range (or co-range) in V_0 + V_1 up to local unitaries. For m ≥ 3, contracting the
free copy of a generic range would put three product vectors that are orthogonal at every copy into V_0 + V_1 on the m
sector copies, and there are none. With no free copy and three copies, the class is the same-code product code of
Wu–Zou [WZ, Remark C.3], and the identity reproduces their closed form.

Checks: `code/check_three_sectors.py` (`logs/check_three_sectors.log`) gives the identity to 8·10⁻¹⁶ and adversarial
minima of q_n/‖C‖² ≥ −1.3·10⁻¹⁶ at n = 4, 5, 6; controls with m − g odd fail at 3·10⁻².

**Products of blocks.** Both identities hold whatever the free copies are, since their derivation never uses the
dimensions of the free copies. So they can be applied to one block of copies at a time. Split the copies into blocks, and let the range and co-range of C lie in S_1 ⊗ … ⊗ S_k, where each
S_i lives on its own block and is one of the following:
- C^3 on a single copy (a free copy), with weight g_i = 1;
- a two-sector space A_i ⊕ B_i on the block, grouped as in Theorem C, with weight g_i = max(a_i, b_i, 1);
- a three-sector space span{1…1, 2…2, 3…3} on the block (in an orthonormal basis at each copy of the block), with
  weight g_i = 1 if the block has odd size and g_i = 2 if it has even size.

**Theorem D.** If g_1 + … + g_k ≤ 3, then q_n(C) ≥ 0 for every C of rank at most two with range and co-range in
S_1 ⊗ … ⊗ S_k. If the weights add up to at most 2, the proof uses only the two-copy theorem.

*Proof.* Induction on the number of blocks that are not single copies; single copies may have any dimension. Collapse
one such block, treating every other copy as free. The identity writes q_n(C) as 2^{g_i−m_i} times the form of UC
(U = 1 for a three-sector block), with that block replaced by g_i copies of dimension at most four, plus
(1 − 2^{g_i−m_i}) ≥ 0 times a sum of forms of operators of rank at most two on the other copies. These operators are the
blocks ⟨α|P|β⟩ and ⟨β|Q|α⟩ of the q̃ terms, or the X_ij. Their ranges and co-ranges lie in the product of the other
S_j. So do those of UC, with the new copies counted as single copies. Both are problems of the same kind, with one block
fewer to collapse and weight sum at most 3. When only single copies are left, there are at most three of them, each of
dimension at most four. The form is then nonnegative: by [T1, Thm. 1.2] for three copies, by the two-copy theorem
[FHPV, BGH] for two, and trivially for one or none. If the weights add up to at most 2, at most two copies are ever
left. ∎

Examples that use only the two-copy theorem and are thick at every copy:
- for even n, C^3 ⊗ span{1^{n−1}, 2^{n−1}, 3^{n−1}} (one free copy and an odd three-sector block);
- for odd n ≥ 3, (L_1 ⊗ e_2 ⊕ e_1 ⊗ L_2) ⊗ span{1^{n−2}, 2^{n−2}, 3^{n−2}}, with L_s a plane and e_s ⊥ L_s;
- at four copies, (L_1 ⊗ e_2 ⊕ e_1 ⊗ L_2) ⊗ (L_3 ⊗ e_4 ⊕ e_3 ⊗ L_4);
- at n = m + m′ copies with m, m′ odd, span{1^m, 2^m, 3^m} ⊗ span{1^{m′}, 2^{m′}, 3^{m′}}.

The first example is the last one with m′ = 1, and was already in "More than two sectors". The second and third were
already covered with [T1]: the second by the two-free-copy case of "More than two sectors", the third by Theorem C with
f = 2, as in the four-copy example there. Here they need only the two-copy theorem. The last with m, m′ ≥ 3 is new.
With [T1] the budget is 3: for example three odd three-sector blocks, or a free copy and two odd ones.

**One side is enough.** The co-range can be left arbitrary. Take a three-sector block of m copies and let R be the
other copies. For C with range in H_R ⊗ span{1…1, 2…2, 3…3}, expand the co-range factor on the block in the local
bases. A partial trace over copies J of the block keeps |j…j⟩⟨w| only if the string w equals j on J. Two such pieces
have orthogonal partial traces unless every copy of the block is traced, and then only the pieces with w = j…j survive,
and these belong to CQ.
This gives the exact identity

    q_n(C) = q_n(CQ) + Σ_j Σ_{w not constant} 2^{−|A(j,w)|} q_R(⟨j…j|C|w⟩),

where Q projects the co-range onto H_R ⊗ span{1…1, 2…2, 3…3}, A(j, w) is the set of copies where w_s = j, and
⟨j…j|C|w⟩ is the operator on H_R obtained by contracting the block. For a two-sector block, the proof of Corollary B
writes q_n(C) − q_n(CQ) as a sum of nonnegative multiples of terms q_R(Y), using the separable form of the one-copy form
at the planar copies. In both cases each Y has rank at most two, and its range lies in the product of the other blocks'
spaces. So q_n(C) ≥ q_n(CQ) once these terms are known to be nonnegative. The compression is the co-range splitting of
[T2, Lemma 4.8], with the local bases of the block in place of span{0, 1, 2}.

**Theorem D (one-sided).** If the weights add up to at most 3, then q_n(C) ≥ 0 for every C of rank at most two with
range in S_1 ⊗ … ⊗ S_k, whatever the co-range. The same holds with range and co-range exchanged.

*Proof.* Induction on the number of blocks that are not free copies. Compress the co-range one block at a time as
above. Each extra term is a one-sided problem on the other blocks, with weight sum at most 3 − g_i, so it is nonnegative
by induction. When only free copies are left, there are at most three of them, and the endpoint theorem covers them.
After the last compression the operator is CP_V = P_V C P_V, with V = S_1 ⊗ … ⊗ S_k. It has rank at most two, and its
range and co-range lie in V, so Theorem D applies. If the weights add up to at most 2, only the two-copy theorem is used.
The statement with range and co-range exchanged follows from q_n(C*) = q_n(C). ∎

For example, for every even n, every C of rank at most two whose range lies in C^3 ⊗ span{1^{n−1}, 2^{n−1}, 3^{n−1}} has
q_n(C) ≥ 0, whatever its co-range in (C^3)^{⊗n}. This uses only the two-copy theorem.

For C with range and co-range in one plane inside S_1 ⊗ … ⊗ S_k (for instance a Hermitian C of rank two), and blocks
that are free copies or three-sector blocks, the budget is 4. The collapses keep C itself (U = 1), and every copy left at
the end is a qutrit. So the last step is the four-copy theorem of `four-replica-notes/` for such C on (C^3)^{⊗4}, which
is computer-assisted and not refereed. Each remainder has weight at most 3 and is covered by two-sided Theorem D.
Two-sector blocks are excluded, because their collapse can turn C into UC with U ≠ 1 and creates copies of dimension
four. Blocks of size two do not shrink, so this adds new cases only from five copies on.

Checks:
- `code/check_one_sided_three_sectors.py` (`logs/check_one_sided_three_sectors.log`): the identity to 6·10⁻¹⁶ and
  adversarial minima ≥ −3.7·10⁻¹⁷ at n = 3, 4, 5.
- `code/check_product_blocks_one_sided.py` (`logs/check_product_blocks_one_sided.log`): mixed products up to n = 6 with
  arbitrary co-range. The minima are ≥ 0 up to rounding; rank-three members reach −0.125 to −0.5.

**Where the content lies.** The copies of type (1, 1) in Theorem C are thin: their local support A_s ⊕ B_s is a plane.
The thin-copy reduction removes them and stays inside the class, since a contraction at such a copy maps A ⊕ B into
the same kind of space on the other copies. After they are removed, f + a + b ≤ f + 2g ≤ 6 − f copies are left. For f = 0, Theorem C is contained in the original
Theorem. Its new cases have f ≥ 1, so they live on at most five thick copies, and on four or five once the known results
for up to three copies are set aside. Three-sector blocks are different: every copy of such a block is
thick. So Theorem D gives classes that are thick at every copy for every n. Both theorems are exact reductions to at most
three copies. What differs is thickness, which by [T2, §6] is where a first failure would have to sit.

Checks:
- `code/check_ghz_blocks.py` (`logs/check_ghz_blocks.log`): the two-stage identity for two three-sector blocks, to
  10⁻¹⁵ (n ≤ 6), and adversarial minima.
- `code/check_product_blocks.py` (`logs/check_product_blocks.log`): mixed products of free copies, two-sector groups and
  three-sector blocks up to n = 6. The adversarial minima of q_n/‖C‖² over rank two are ≥ 0 up to rounding, while
  rank-three members of the same supports reach −0.125 to −0.5, so the search finds negative values when they exist.

**An LOCC reading.** For odd m, project m copies of the endpoint Werner state locally, on both sides, onto the
repetition code span{1…1, 2…2, 3…3}. By the identity of "More than two sectors" with no free copy, the result is (up to
normalisation) the two-qutrit state ρ_m whose partial transpose is

    2^{1−m} (1 − ½|Ω⟩⟨Ω|) + (1 − 2^{1−m}) Σ_{i≠j} |ij⟩⟨ij|,   Ω = Σ_i |ii⟩.

It has full rank and is NPT for every m: the smallest eigenvalue of its partial transpose is −2^{−m}. Its Werner twirl
has parameter a = −(1 + 2^{1−m})/(3 + 2^{1−m}), which tends to the PPT boundary −1/3 (a = −5/13 for m = 3). The endpoint
inequality on products of k blocks of size m is exactly k-copy undistillability of ρ_m, so Theorem D says that ρ_m is
undistillable for up to three copies. This does not make the general problem easier. ρ_m^Γ is a mixture of the one-copy
endpoint operator and a separable operator, so k-copy undistillability of ρ_m follows from the endpoint inequality for
up to k copies, and it in turn implies undistillability of the Werner twirl. Still, undistillability of ρ_m for every k,
for a single m, would already give NPT bound entanglement. Checked by `code/check_repetition_state.py`
(`logs/check_repetition_state.log`, m = 1, 3, 5): the compression identity holds exactly. This remark has not had a
separate referee; it is a direct reading of the refereed identity.

## The two-copy critical point

The lowest thick critical point of R = q_2/‖·‖² on the rank-two variety that we found at two copies is in this class
(`code/critical_n2.npy`). Its value is

    R = (9 − √41)/20 ≈ 0.12984,

and its singular-value ratio is (√41 − 4)/5. Three of four biased searches converged to it.

In the eigenbasis of the one-copy reduction of its range (one local unitary per copy, applied on both sides), it
commutes with Z ⊗ Z, where Z = diag(1, −1, 1). Its range and co-range lie in

    (K ⊗ K) ⊕ C|11⟩,    K = span{|0⟩, |2⟩},

which is the one-excitation case. They differ by a rotation of about 30° inside K at each copy.

It is a saddle on the rank-two variety: 6 negative Hessian directions, the smallest −0.58. By the Theorem, no operator
of this shape violates the endpoint inequality, for any number of copies.

## Checks

- `code/check_one_excitation.py` (`logs/check_one_excitation.log`, n = 2, 3, 4), one-excitation case:
  - the identity on 20 random members (error ≤ 2.2·10⁻¹⁶);
  - the qubit lemma, by adversarial minimisation of the gap (minimum +2.7·10⁻¹⁸);
  - the minimum of q_n/‖C‖² over the class (≥ −3.5·10⁻¹⁶);
  - that the two-copy critical point lies in the class (weight outside 7·10⁻²⁸).
- `code/check_two_sectors.py` (`logs/check_two_sectors.log`), general case: random orthogonal subspaces of random
  dimensions (1,1), (2,1), (1,2) in random bases, at n = 2, 3, 4.
  - The identity of Step 1 holds to 4.4·10⁻¹⁵.
  - The minimum of q_n/‖C‖² over the class is ≥ −3.3·10⁻¹⁶.
- `code/check_one_sided.py` (`logs/check_one_sided.log`), the Corollary, with ran C in V and the co-range arbitrary,
  at n = 2, 3, 4.
  - The difference q_n(C) − q_n(P_V C P_V) is at least 0.17 ‖C‖² on random members.
  - By adversarial minimisation it is at least −8.9·10⁻¹⁶ ‖C‖²; zero is attained when the co-range lies in V.
- `code/check_other_dimensions.py` (`logs/check_other_dimensions.log`): d = 4, 5, including copies with two orthogonal
  planes; the minimum of q_n/‖C‖² over the one-sided class is ≥ −5·10⁻¹⁷.
- `code/check_free_copy.py` (`logs/check_free_copy_n3.log`, `logs/check_free_copy_n4.log`), one free copy, random
  sector dimensions and bases:
  - the identity, to 4.8·10⁻¹⁵;
  - Lemma L: the ratio |B_1|²/(ψψ) never exceeds 1 (random and adversarial); it equals 1 at W = Z of rank one;
  - Theorem A, adversarial minimum of q_n/‖C‖² ≥ −2.6·10⁻¹⁶;
  - Corollary B, adversarial minimum ≥ −2.6·10⁻¹⁶, and q_n(C) − q_n(P_V C P_V) ≥ −2.8·10⁻¹⁷.

- `code/check_free_copies.py` (`logs/check_free_copies_f2_m2.log`, `logs/check_free_copies_f2_m3.log`,
  `logs/check_free_copies_f1_m3.log`), f free copies and m sector copies (n = 4 and 5, including planar sectors on both
  sides). It checks:
  - the identity with q_F, to 3.9·10⁻¹⁵;
  - the inequality Σ_μ q_F(X̂_μ) ≥ 2Σ q_F(x_α y_β^*) − q_F(X̂_0) (gap ≥ 0, equal to 0 when A has at most one planar
    factor);
  - the reduction 2^m q_n(C) ≥ 2 q_{F∪{K}}(𝒞) + 2^m(q̃_F(P) + q̃_F(Q)), with slack ≥ −1.5·10⁻¹⁵;
  - Theorem A₂ and Corollary B₂ (with f = 1 also Theorem A and Corollary B), adversarial minima ≥ −2.1·10⁻¹⁶.

- `code/check_one_plane_each.py`, `code/check_grouped_sectors.py`, `code/check_tomar8_shared_tails.py`
  (`logs/check_one_plane_each.log`, `logs/check_grouped_sectors.log`): the grouping identity for random C of any rank,
  to 2.5·10⁻¹⁵, and adversarial minima for Theorem C and the [T2, §8] shapes; see "Grouping the sector copies".
- `code/check_group_gap.py` (`logs/check_group_gap.log`): the gap when two planes share a group.
- `code/check_three_sectors.py` (`logs/check_three_sectors.log`): three sectors, the identity and adversarial minima.
- `code/check_ghz_blocks.py`, `code/check_product_blocks.py` (`logs/check_ghz_blocks.log`,
  `logs/check_product_blocks.log`): products of blocks (Theorem D).
- `code/check_one_sided_three_sectors.py`, `code/check_product_blocks_one_sided.py`
  (`logs/check_one_sided_three_sectors.log`, `logs/check_product_blocks_one_sided.log`): the one-sided version.

An independent referee re-derived every step of the sections before "One free copy" and re-ran the checks with
separately written code, up to n = 5. It found no mathematical error. Its wording and attribution fixes are applied
above. A second independent referee did the same for "One free copy". It re-derived every step and wrote its own checks
at n = 2 to 5, including Theorem A at ranks three to six and the pieces of the proof of Corollary B. It found no
mathematical error, and its fixes (the proof of Corollary B, the scope statement, hypotheses and wording) are applied
above. It also confirmed that the rank conditions matter: with rank X = 2, rank Y = 1 and P = Q = 0 the form goes
negative. A third independent referee did the same for "Two free copies" at n = 4 and 5. It checked the identity,
the replica form, the reduction (which holds only with the stated phase), Theorem A₂, Corollary B₂ and control runs
showing that the rank hypotheses are needed. It found no mathematical error, and its fixes are applied above. A
fourth independent referee did the same for "Grouping the sector copies" at n = 1 to 5 (652 identity tests, every
valid g, ranks one to full, controls with the wrong sign of U and with two planes in one group). It found no
mathematical error; its fixes (the derivation of the identity, the embedding, the one-sided statement, the example's
notation, the §8 subspaces and the gap remark) are applied above. A fifth independent referee re-derived the identity of "More than two sectors" and repeated the checks with its own
code at n = 1 to 6: 414 identity tests for three sectors and 200 for two sectors and d = 4, adversarial minima at
n = 4, 5, 6, and controls showing that the optimiser finds negative values at rank three and at t = 0.55. It found no
mathematical error, and its fixes (U for general r, the parity condition, the bases, the theorems used, novelty) are
applied above. A sixth independent referee checked "Products of blocks" and "Where the content lies" with its own code at n = 2
to 6. It ran the induction to its leaves for 17 configurations of blocks: the identity is exact, the coefficients are
nonnegative, and the leaves have rank at most two on at most three copies of dimension at most four. It also ran
adversarial minima for the examples and for weight-3 products, with rank-three and t = 0.55 controls, and checked the
thin-copy reduction. It found no mathematical error, and its fixes are applied above. A seventh independent referee checked "One side is enough" with its
own code at n = 3 to 6. It found the one-sided identities exact, ran the one-sided induction end to end, and ran
adversarial minima for one-sided and Hermitian weight-4 products, with rank-three and t = 0.55 controls. It found no
mathematical error, and its fixes are applied above.

## References

- [T1] T. S. Tomar, *Undistillability of Werner states from three, four and five copies*, Zenodo 23718035 (2026):
  Theorem 1.2 (three copies at the endpoint, arbitrary factor dimensions; computer-assisted).
- [FHPV] Fraser, Huber, Pozsgay and Vona, arXiv:2607.24309; [BGH] Bharti, Gajjala and Haug, arXiv:2607.24479: the
  two-copy endpoint in arbitrary factor dimensions.
- [T2] T. S. Tomar, *Equality cases and every-copy results for the Werner endpoint*, Zenodo (2026):
  - Lemma 2.6 (range in a product of qubits);
  - Theorems 4.1–4.7 (GHZ backbones with single-copy tails);
  - Lemma 4.8 (co-range compression for range in V);
  - Section 8 (open cases).
- [WZ] Wu and Zou, arXiv:2608.02647:
  - Theorem 7.3 (product ray, three copies);
  - Appendix C.2 (locally orthogonal product codes, three copies).
- [SC] Song and Chen, arXiv:2607.23416: Proposition 8 (rank-one bound).
