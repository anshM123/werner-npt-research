# Ranges in the one-excitation space

Ansh Mishra and Aryan Senthilkumar, 3 November 2026

We prove the endpoint inequality for every number of copies on one class of rank-two operators. The range lies among
the single excitations over a product vector. The co-range has at most two excitations. The class is thick at every
copy. It consists of exactly the planes where the anchor method of `support-and-projector-notes/` (Section 5) breaks
down, because single excitations over a product are not anchors.

## Setting

For an operator X on (C^3)^{⊗n} let

    q_n(X) = Σ_J (−1/2)^{|J|} ‖Tr_J X‖_F²,

with J running over subsets of the n copies. The Werner endpoint asks whether q_n(X) ≥ 0 for every X of rank at most
two. That question is still open, and this note does not settle it.

Fix a product vector b = b_1 ⊗ ⋯ ⊗ b_n. By a local unitary we may take b_s = |0⟩ at every copy. Write H_c for the span
of the basis strings with exactly c nonzero entries. These are the vectors with c "excitations" over b, and H_0 is
spanned by b itself.

## Result

**Theorem.** Let X have rank at most two, with range(X) ⊆ H_1 and range(X^*) ⊆ H_0 ⊕ H_1 ⊕ H_2. Then q_n(X) ≥ 0 for
every n. Since q_n(X^*) = q_n(X), the same holds with range and co-range exchanged.

Remarks.
- **Thick at every copy.** Take the range spanned by two generic single-excitation vectors w_k = Σ_s |0…x_{k,s}…0⟩ with
  x_{k,s} ⊥ |0⟩. Then the range marginal at copy s is span{|0⟩, x_{1,s}, x_{2,s}}, of dimension three.
- **No product vector in the range plane.** A product vector in H_1 is excited at only one copy. So for n ≥ 2 the
  range plane of a generic X contains no product vector, and the product-ray theorem does not apply.
- **The top sector is free.** The sector c = n (every copy excited on the co-range side) is nonnegative for every rank,
  by the same formula. The open cases are the sectors 3 ≤ c ≤ n − 1, so the first open case is n = 4, c = 3.

## Proof

**1. Charge splitting.** Let U_φ = ⊗_s (|0⟩⟨0| + e^{iφ}(|1⟩⟨1| + |2⟩⟨2|)). The form q_n is invariant under X ↦ UXU^* for
local unitaries U, with the same U on both sides. If the range lies in H_1, the piece X P_c (P_c the projector onto
H_c) picks up the phase e^{i(1−c)φ}. Averaging over φ removes the cross terms:

    q_n(X) = Σ_c q_n(X P_c).

Each X P_c has rank at most two. The piece c = 0 has rank at most one, and rank-one operators have q_n ≥ 0. So it is
enough to treat a single sector c ≥ 1.

**2. The sector formula.** Write the entries of X : H_c → H_1 as X[(s, z); σ]. The row (s, z) is the string with value
z ∈ {1, 2} at copy s and 0 elsewhere, and the column σ is a string of weight c.

The partial trace over J keeps only entries whose row and column agree on J, and this fails as soon as J meets the
support of σ outside s. That leaves two kinds of term:
- entries with s ∉ J and J ∩ supp σ = ∅;
- entries with s ∈ J, σ_s = z and J ∖ {s} disjoint from supp σ. These land on the vacuum row.

Summing the signs over J gives

    2^{n−c+1} q_n(X) = 4·Off + 2·Diag + Σ_ρ ( |T_ρ|² − 2 Σ_{s free in ρ} |τ_{ρ,s}|² ),

where:
- Off is the sum of |X[(s,z);σ]|² over entries with s ∉ supp σ, and Diag is the same sum over entries with s ∈ supp σ;
- ρ runs over strings of weight c − 1;
- τ_{ρ,s} = Σ_{z=1,2} X[(s,z); ρ + z@s], where ρ + z@s puts the value z at the free copy s;
- T_ρ = Σ_s τ_{ρ,s}, summed over the free copies s of ρ.

Grouping the terms by ρ gives 2^{n−c+1} q_n(X) = Σ_ρ F_c(X^ρ). Here X^ρ is the submatrix with rows (s, z), s free in ρ,
and columns ρ + z′@t, t free in ρ. It is a block matrix with 2 × 2 blocks Y_st. Each off-diagonal entry of X lies in
exactly c of these submatrices, and each diagonal entry in exactly one. With the off-diagonal weight split evenly,

    F_c(Y) = (4/c) Σ_{s≠t} ‖Y_st‖² + 2 Σ_s ‖Y_ss‖² + |Σ_s Tr Y_ss|² − 2 Σ_s |Tr Y_ss|².

For c = 2 this is F_2(Y) = 2‖Y‖² + |Tr Y|² − 2 Σ_s |Tr Y_ss|². For c = 1 there is one string ρ (the empty one), and
F_1 = F_2 + 2 Σ_{s≠t} ‖Y_st‖².

**3. Lemma (rank two, diagonal traces).** Let Y have rank at most two and be made of square blocks Y_st, of any size
and any number m. Then

    Σ_s |Tr Y_ss|² ≤ ‖Y‖² + ½ |Tr Y|².

*Proof.*
- **Gauge.** Write Y = U V^* with U, V having two columns, and use the GL(2) freedom to make U^*U = I_2. Then
  ‖Y‖² = ‖V‖², and Tr Y_ss = ⟨V_s, U_s⟩ is the Hilbert–Schmidt pairing of the s-th blocks.
- **Reduction to a matrix.** For fixed U, the slack ‖V‖² + ½|Σ_s φ_s(V)|² − Σ_s |φ_s(V)|², with φ_s(V) = ⟨U_s, V_s⟩,
  is a quadratic form in V.
  - The functionals φ_s act on different blocks, so their Riesz vectors are orthogonal. Their norms u_s = ‖U_s‖_F
    satisfy Σ_s u_s² = Tr U^*U = 2.
  - On the orthogonal complement of these vectors the form is ‖V‖² ≥ 0.
  - On their span it is the matrix M = I − diag(u_s²) + ½ u uᵀ.
- **Positivity of M.** If every u_s² ≤ 1, then M ⪰ I − diag(u_s²) ⪰ 0. Otherwise, since the u_s² add up to 2, exactly
  one of them is larger than 1, say u_0² = 1 + δ with 0 < δ ≤ 1.
  - If δ = 1 (all other u_s vanish, including the case m = 1), M is diagonal apart from u_0, and M_00 = −1 + 1 = 0.
  - If δ < 1, D = I − diag(u_s²) has exactly one negative eigenvalue and the others are positive. Adding the positive
    rank-one term ½uuᵀ cannot lower any eigenvalue, so by interlacing λ_2(M) ≥ λ_2(D) > 0. Hence M ⪰ 0 exactly when
    det M = det D · (1 + ½ uᵀD^{−1}u) ≥ 0, that is (as det D < 0) when 1 + ½ uᵀD^{−1}u ≤ 0. That condition reads

        Σ_{s≠0} g(u_s²) ≤ g(1 − δ),   g(x) = x/(1 − x).

    It holds because g is convex with g(0) = 0, hence superadditive, and Σ_{s≠0} u_s² = 1 − δ. ∎

The rank hypothesis cannot be dropped: the rank-three matrix Y = diag(I_2, −|e_1⟩⟨e_1|), with two blocks, has
Σ_s |Tr Y_ss|² − ‖Y‖² − ½|Tr Y|² = 5 − 3 − ½ = 1.5 > 0.

Rank two enters only through Σ_s u_s² = 2, which allows at most one "bad" block. Equality holds, for example, for
Y = |1⟩⟨1| ⊗ I_2 and for Y = [[I, I], [−I, −I]].

**4. Conclusion.** Every X^ρ is a submatrix of X, so it has rank at most two. For c = 2 the Lemma gives
F_2(X^ρ) ≥ 0 for every ρ. For c = 1, F_1 ≥ F_2 ≥ 0. Together with step 1 this proves the theorem. ∎

## Why the sectors c ≥ 3 are different

The even split gives each string the weight 4/c < 2 on its off-diagonal blocks, and the Lemma needs 2. This loss is
real. The rank-two matrix Y = [[I, I], [−I, −I]] has F_c(Y)/‖Y‖² = 2/c − 1, which is negative for every c ≥ 3. So a
proof for c ≥ 3 has to share the off-diagonal weight unevenly, or use the rank-two coupling between different strings.

Numerically the full sector form stays nonnegative. Adversarial minimisation at n = 4, 5 with co-range in H_3 gives
minima of order 10^{−6} to 10^{−8}, approached at thin operators whose range sits on one or two copies. The referee's
independent searches found minimum 0 (to machine precision) for co-range in H_3 at n = 3, …, 6, in H_0 + ⋯ + H_3 at
n = 4, 5, and in H_4 at n = 4, 5.

Two partial facts hold for every sector.
- If the range of X is spread over at most three copies, every sector reduces to the cases already proved. The copies
  outside the range support carry fixed labels or holes and decouple, which leaves an effective problem on at most
  three copies, where c ≤ 2 or c equals the number of copies. So the first open case needs four copies, with the range
  supported on all four. (This sub-case is not new as a class: such a range is thin at every other copy, so Tomar's
  thin-copy removal and the three-copy theorem also cover it. The sector argument gives a direct proof.)
- For c = n − 1 the sector inequality takes a compact form. Use the gauge Σ_s U_s^*U_s = I and write F_u for the
  column vectors whose hole is at copy u, as a vector in (C^2)^{⊗([n]∖u)} ⊗ C^2. The inequality is then

      Σ_{s<t} ‖K_s F_t − K_t F_s‖² ≤ Σ_u ⟨F_u| (2I + 2U_u^*U_u) |F_u⟩,

  where K_s contracts the label at copy s and the C^2 index with the matrix U_s. In numerical tests the best
  constant is exactly 1 for n = 2, …, 5, so the inequality is tight.

## Verification

The scripts need Python 3 with NumPy and SciPy. Each one ends with an assertion and prints PASS.
- `code/check_sector_formula.py` compares the sector formula with direct partial traces for n = 2, …, 5 and every c,
  on random rank-two and full-rank X. It also checks the charge splitting q(X) = Σ_c q(X P_c). The largest relative
  error is 1.5 × 10^{−15}.
- `code/check_lemma.py`:
  - adversarial minimisation of the Lemma's slack for block sizes 1, 2, 3 and m ≤ 8 (minimum +9 × 10^{−14});
  - the equality cases;
  - the matrix step for random and adversarial u;
  - the explicit failure of the even split for c ≥ 3.
- `code/check_class.py` brute-forces q_n/‖X‖² on the class of the theorem for n = 2, …, 5, with random and adversarial
  starts (minimum +1.6 × 10^{−9}). It also runs the open sector c = 3 at n = 4, 5, as evidence only.

The logs are in `logs/`.
