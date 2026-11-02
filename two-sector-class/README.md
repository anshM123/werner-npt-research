# Two orthogonal product sectors

Ansh Mishra and Aryan Senthilkumar, 2 November 2026

This note proves the endpoint inequality for every number of copies on a class of rank-two operators. Their range lies
in the sum of two product subspaces that are orthogonal at every copy; the co-range can be arbitrary (or the other way
round). Operators in this class can have one-copy reductions of rank three at every copy. The class contains the lowest
thick critical point we found at two copies. It is a small class, and the argument does not extend to general
operators.

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
- **What is new.** The two-sided class with a ≥ 2 or b ≥ 2 (for example one excitation) at n ≥ 4, and the one-sided
  Corollary.
- **An open case of [T2, §8].** Let T_0 and T_1 be disjoint sets of copies, and take
  A = ⊗_{s∈T_0} span{0,2} ⊗ ⊗_{s∉T_0} span{0} and B = ⊗_{s∈T_1} span{1,2} ⊗ ⊗_{s∉T_1} span{1}. These are orthogonal
  at every copy. So the Corollary covers every rank-two operator whose range lies in
  span{0…0, 2_i 0…0 (i ∈ T_0), 1…1, 2_j 1…1 (j ∈ T_1)} and whose co-range is arbitrary: the two GHZ backbones of
  [T2] with their tails on disjoint sets of copies. [T2, §8] lists range in V_0 + V_1 with an arbitrary co-range as
  open, except for operators whose range is spanned by two vectors of V_0 ∪ V_1 or whose co-range is orthogonal to
  all partner strings. What stays open in that case is the ranges in which a tail of the 0…0 backbone and a tail of
  the 1…1 backbone sit at the same copy, outside Tomar's conditions.
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

The first sum is 2^m q_n(X). It is nonnegative when X has rank at most two. At every sector copy the range of X lies in
the plane or line A_s, and there the separable decomposition of 1 − |ω_K⟩⟨ω_K| (the thin-copy reduction) removes the
copy. Doing this at each sector copy leaves the one-copy inequality at copy 0. The same holds for Y. Also q̃(P) ≥ 0
when P has rank at most two, since q̃(P) is a sum of q_1 over 3 × 3 blocks of rank at most two. So only the last term
can be negative.

**Lemma L.** For square matrices Z and W of any size, |B_1(Z, W)|² ≤ ψ(Z) ψ(W), where ψ(Z) = ‖Z‖₁² − ½|Tr Z|².

*Proof.* First let z and w have rank one. For every λ > 0 and every phase θ, the matrix λz − e^{iθ}w/λ has rank at most
two, so q_1 of it is nonnegative. This gives |B_1(z, w)|² ≤ q_1(z) q_1(w), and q_1 = ψ in rank one. In general, take
singular value decompositions Z = Σ_k s_k u_k v_k^* and W = Σ_l t_l u′_l v′_l^*. Expanding B_1 bilinearly and using the
rank-one case term by term gives |B_1(Z, W)| ≤ S_Z S_W with S_Z = Σ_k s_k (1 − ½|c_k|²)^{1/2}, c_k = ⟨v_k, u_k⟩. Put p_k = s_k/‖Z‖₁. Concavity of the
square root and Σ_k p_k|c_k|² ≥ |Σ_k p_k c_k|² give S_Z² ≤ ‖Z‖₁²(1 − ½|Σ_k p_k c_k|²) = ψ(Z), since
Σ_k s_k c_k = Tr Z. ∎

**Lemma R.** If X = xy^* has rank one, then Σ_μ q_1(X̂_μ) ≥ ψ(X̂_0).

*Proof.* Write F for the swap of the two replicas. In replica form,
Σ_μ q_1(X̂_μ) = ⟨x ⊗ y| W_0 ⊗ ⊗_s (2 − F_s) |x ⊗ y⟩, with the product over the planar copies of A and
W_0 = 1 − ½F_0 ≥ 0; one-dimensional copies contribute a factor 1. For f_s = ±1 we have Π_s(2 − f_s) ≥ 2 − Π_s f_s, so
⊗_s(2 − F_s) ≥ 2 − F_A. Moreover W_0 ⊗ (2 − F_A) = 1 − ½F_0 ⊗ F_A + 4P_0^− ⊗ P_A^− ≥ 1 − ½F_0 ⊗ F_A, where the P^− are
antisymmetric projectors. The expectation of the right side is ‖x‖²‖y‖² − ½|⟨y, x⟩|². Finally ‖X̂_0‖₁ ≤ ‖x‖‖y‖ and
Tr X̂_0 = ⟨y, x⟩. ∎

**Theorem A.** For every n, if rank X ≤ 1, rank Y ≤ 1, rank P ≤ 2 and rank Q ≤ 2, then q_n(C) ≥ 0. Here C itself may
have rank up to six.

*Proof.* By the identity, Lemma R (the X-sum is 0 if X = 0, and then so is X̂_0), q̃ ≥ 0 and Lemma L,

    2^m q_n(C) ≥ ψ(X̂_0) + ψ(Ŷ_0) − 2 (ψ(X̂_0) ψ(Ŷ_0))^{1/2} ≥ 0.  ∎

**Corollary B.** For every n, if C has rank at most two and its range is spanned by one vector of C^3 ⊗ A and one
vector of C^3 ⊗ B, then q_n(C) ≥ 0, whatever the co-range.

*Proof.* The proof of the Corollary above goes through with the co-range split only at the sector copies. Each piece is
thin at the sector copies outside T, so its form reduces to q_1 at copy 0 on pieces of rank at most two, which is
nonnegative. This gives q_n(C) ≥ q_n(P_V C P_V). The blocks of P_V C P_V satisfy the hypotheses of Theorem A. ∎

The same argument covers the block-triangular members of rank at most two (P = 0 or Q = 0). There
rank X + rank Y ≤ 2, and if one of X and Y vanishes, so does the cross term.

**Scope.** At n = 2 this class is all of C^3 ⊗ C^3, since with one sector copy of dimensions (2, 1) we have
A_1 ⊕ B_1 = C^3. The general rank-two member therefore contains the full two-copy theorem, and we prove only the cases
above. For n ≤ 3 they also follow from the known two- and three-copy results, so the new content is at n ≥ 4. In the
setting of [T2, §8] with one shared tail copy, ranges spanned by one vector of each backbone space are already covered
there. The general case, with a diagonal block of rank two, is open. There the q̃ terms can vanish on projector-like
off-diagonal blocks, so the Schur-complement step of Step 3 has nothing to pay with.

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
  - Lemma L, with largest ratio |B_1|²/(ψψ) = 0.9957 (random and adversarial);
  - Theorem A, adversarial minimum of q_n/‖C‖² ≥ −2.6·10⁻¹⁶;
  - Corollary B, adversarial minimum ≥ −2.6·10⁻¹⁶, and q_n(C) − q_n(C P_V) ≥ −2.8·10⁻¹⁷.

  The free-copy section was re-derived step by step and checked with separately written code. It has not had an
  outside referee.

An independent referee re-derived every step and re-ran the checks with separately written code, up to n = 5. It found
no mathematical error. Its wording and attribution fixes are applied above.

## References

- [T2] T. S. Tomar, *Equality cases and every-copy results for the Werner endpoint*, Zenodo (2026):
  - Lemma 2.6 (range in a product of qubits);
  - Theorems 4.1–4.7 (GHZ backbones with single-copy tails);
  - Lemma 4.8 (co-range compression for range in V);
  - Section 8 (open cases).
- [WZ] Wu and Zou, arXiv:2608.02647:
  - Theorem 7.3 (product ray, three copies);
  - Appendix C.2 (locally orthogonal product codes, three copies).
- [SC] Song and Chen, arXiv:2607.23416: Proposition 8 (rank-one bound).
