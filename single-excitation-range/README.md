# Ranges in the one-excitation space

Ansh Mishra and Aryan Senthilkumar, 3 November 2026

We prove the endpoint inequality on rank-two operators whose range lies among the single excitations over a product
vector. These are exactly the planes where the anchor method of `support-and-projector-notes/` (Section 5) breaks down,
because single excitations over a product are not anchors.
- In Theorem 1 the range plane is arbitrary in that space and the co-range has at most two excitations.
- In Theorem 2 the range is a weighted W-plane (this includes the W-plane itself) and the co-range is arbitrary.
- Theorem 3 (added 4 November 2026) settles the co-range sector with n − 1 excitations, for every range plane in that
  space and every n. Together with Theorem 1 this gives the endpoint at four copies for every such range plane and
  every co-range (Corollary 4). The sectors 3 ≤ c ≤ n − 2 remain open; the first open case is n = 5, c = 3.

Both classes are thick at every copy.

## Setting

For an operator X on (C^3)^{⊗n} let

    q_n(X) = Σ_J (−1/2)^{|J|} ‖Tr_J X‖_F²,

with J running over subsets of the n copies. The Werner endpoint asks whether q_n(X) ≥ 0 for every X of rank at most
two. That question is still open, and this note does not settle it.

Fix a product vector b = b_1 ⊗ ⋯ ⊗ b_n. By a local unitary we may take b_s = |0⟩ at every copy. Write H_c for the span
of the basis strings with exactly c nonzero entries. These are the vectors with c "excitations" over b, and H_0 is
spanned by b itself.

## Result

**Theorem 1.** Let X have rank at most two, with range(X) ⊆ H_1 and range(X^*) ⊆ H_0 ⊕ H_1 ⊕ H_2. Then q_n(X) ≥ 0 for
every n. Since q_n(X^*) = q_n(X), the same holds with range and co-range exchanged.

Remarks.
- **Thick at every copy.** Take the range spanned by two generic single-excitation vectors w_k = Σ_s |0…x_{k,s}…0⟩ with
  x_{k,s} ⊥ |0⟩. Then the range marginal at copy s is span{|0⟩, x_{1,s}, x_{2,s}}, of dimension three.
- **No product vector in the range plane.** A product vector in H_1 is excited at only one copy. So for n ≥ 2 the
  range plane of a generic X contains no product vector, and the product-ray theorem does not apply.
- **The top sector is free.** The sector c = n (every copy excited on the co-range side) is nonnegative for every rank,
  by the same formula. For a general range plane in H_1 the open cases are the sectors 3 ≤ c ≤ n − 1, so the first
  open case is n = 4, c = 3. Theorem 2 closes all of them for weighted W-planes.

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

## Theorem 2: W-planes, with any co-range

For one family of range planes inside H_1, every co-range sector can be handled. Among them is the W-plane span{W_1, W_2},
W_z = Σ_s |0…0 z_s 0…0⟩, the standard example of a plane on which single-vector anchors fail
(`support-and-projector-notes/`, Section 5).

**Theorem 2.** Let λ_1, …, λ_n > 0 with max_s λ_s² ≤ 2 min_s λ_s², and let W_s be a unitary on span{|1⟩, |2⟩} at copy
s. Put w_k = Σ_s λ_s |0…(W_s e_k)_s…0⟩ for k = 1, 2, where e_1 = |1⟩ and e_2 = |2⟩. Then q_n(C) ≥ 0 for every C of rank
at most two whose range lies in span{w_1, w_2}, with arbitrary co-range, for every n. In particular this holds for the
W-plane (all λ_s equal, W_s = I).

These planes are thick at every copy and contain no product vector. For the W-plane the exact minimum of q_n/‖C‖² over
all co-ranges is 2^{1−n} for n = 3, 4, 5 (`code/check_wplanes.py`), in line with the numerics of Section 5 there; the
referee found 0, 1/4, 1/4, 1/8, 1/16, 1/32 for n = 1, …, 6, attained in the sector c = 0. The rotations W_s do not
change the minimum. For n ≤ 3 the statement also follows from the two- and three-copy theorems; it is new from n = 4.

*Proof.*
- **Reduction.** A product unitary that acts as W_s^{−1} on span{|1⟩, |2⟩} at each copy leaves q_n unchanged and
  turns W_s into I. Rescaling all λ_s together does not change the plane, so assume Σ_s λ_s² = 1. In the notation of
  the first proof, the range factor is U with U_s = λ_s I. So U^*U = I, and A_s := U_s^*U_s = λ_s² I.
- **Sectors.** By step 1 of the first proof it is enough to treat X = C P_c for each c. The case c = 0 has rank one.
- **Hole form.** Fix c ≥ 1 and put h = n − c.
  - For a column σ of weight c let H(σ) be its set of holes (copies where σ is 0). Write X[(s,z); σ] = (U_s f_σ)_z with
    f_σ ∈ C². Then Off_σ = f_σ^* A_{H(σ)} f_σ and Diag_σ = f_σ^* (I − A_{H(σ)}) f_σ, where A_H := Σ_{u∈H} A_u = λ_H² I.
  - The sector formula therefore reads

        2^{n−c+1} q_n(X) = Σ_σ f_σ^* (2I + 2A_{H(σ)}) f_σ − Σ_ρ ( 2 Σ_{s free in ρ} |τ_{ρ,s}|² − |T_ρ|² ).

  - Collect the columns with a given hole set H into one vector F_H ∈ (C²)^{⊗([n]∖H)} ⊗ C²_V: the labels at the
    excited copies index the first factor, and the value f_σ is the last one.
  - Then τ_{ρ,s} = λ_s (⟨Φ|_{s,V} F_H)(labels of ρ), with H = H(ρ) ∖ {s} and Φ = Σ_{z=1,2} |z⟩_s |z⟩_V. Each pair
    (ρ, s) with s free in ρ corresponds to exactly one triple (H, s ∉ H, labels).
- **h ≥ 1.** In the formula, drop the terms +|T_ρ|² ≥ 0; this can only lower the right-hand side. What is subtracted
  is then 2 Σ_ρ Σ_s |τ_{ρ,s}|², and the regrouping above turns it into an equality:

      2 Σ_ρ Σ_s |τ_{ρ,s}|² = 2 Σ_H Σ_{s∉H} λ_s² ‖⟨Φ|_{s,V} F_H‖² ≤ 2 Σ_H (1 − λ_H² + max_{s∉H} λ_s²) ‖F_H‖².

  The inequality is the lemma below, applied with |Φ⟩⟨Φ| = 2P (P a maximally entangled projector). The positive sum
  is Σ_H (2 + 2λ_H²) ‖F_H‖². The difference is at least Σ_H (4λ_H² − 2 max_{s∉H} λ_s²) ‖F_H‖². This is ≥ 0, because
  λ_H² ≥ min λ² ≥ ½ max λ² when H ≠ ∅.
- **h = 0** (c = n). Every ρ has one free copy, so its bracket is |τ_{ρ,s}|² − 2|τ_{ρ,s}|² = −|τ_{ρ,s}|². The net
  subtracted amount is Σ_s λ_s² ‖⟨Φ|_{s,V} F_∅‖² = 2⟨F_∅| Σ_s λ_s² P_s |F_∅⟩ ≤ 2‖F_∅‖², since Σ_s λ_s² P_s ⪯ I. That is
  the positive sum. ∎

**Lemma (weighted central spin).** Let P_s be the projector onto a maximally entangled state of (qubit s) ⊗ V, and let
w_s ≥ 0. Then λ_max(Σ_s w_s P_s) ≤ (Σ_s w_s + max_s w_s)/2.

*Proof.* Every maximally entangled state of qubit s and V has the form (I ⊗ Y)|ψ⁻⟩ with Y unitary. Since
(Y ⊗ Y)|ψ⁻⟩ = det(Y)|ψ⁻⟩, this equals det(Y)·(Y^{−1} ⊗ I)|ψ⁻⟩, a unitary on qubit s alone applied to the singlet.
Conjugating by the product of these satellite unitaries reduces to singlets. For a set S of qubits, Σ_{s∈S} P_s = |S|/4 − **S**_S·**S**_V, with
**S**_S the total spin of S. Its largest eigenvalue is attained at total spin |S|/2 coupled down to |S|/2 − ½, and it
equals (|S| + 1)/2. For general weights, write w = Σ_k (θ_k − θ_{k+1}) 1_{S_k} over the level sets
S_k = {s : w_s ≥ θ_k}, with θ_1 > θ_2 > ⋯ the distinct values and the last θ equal to 0. The largest eigenvalue is
subadditive, so it is at most Σ_k (θ_k − θ_{k+1})(|S_k| + 1)/2 = (Σ_s w_s + θ_1)/2. ∎

Equal weights attain the bound. For the W-plane the proof gives (n − h + 1)/(n + h) < 1 for the ratio of the two sums
in every sector with h ≥ 1. The exact values are smaller, for example (n + 2)/(2n + 2) at h = 1, computed for n ≤ 8.

The weight condition comes from dropping −|T_ρ|². Numerically it is not needed. The exact sector ratios stay at most
0.84 for strongly unequal weights (one heavy copy, one light copy, geometric weights, two heavy copies; n = 4, 5, 6).
Brute force over all co-ranges outside the condition also gives positive minima. The argument can keep part of
−|T_ρ|². Per string, one coefficient may drop below 2 if the others rise: Σ_s c_s|τ_s|² ≥ 2Σ_s|τ_s|² − |Σ_s τ_s|²
exactly when diag(c − 2) + J ⪰ 0. Combined with the lemma, this certifies every weight vector we tested, two heavy
copies included, provided the coefficients are chosen by a small semidefinite program. The simple explicit rules we
tried fail. So the general weighted case is left open until an explicit rule is found.

## Theorem 3: the sector c = n − 1, for every n

Added 4 November 2026. In the gauge Σ_s U_s^*U_s = I, let F_u be the part of the co-range whose hole is at copy u, as
a vector in C^2 ⊗ (C^2)^{⊗([n]∖u)} (the reference index, then the letters of the other copies). K_s contracts the letter
at copy s and the reference index with U_s. With A_u = U_u^*U_u, the sector c = n − 1 reads

    Σ_{s<t} ‖K_s F_t − K_t F_s‖² ≤ Σ_u ⟨F_u| (2I + 2Ā_u) |F_u⟩.

This is the hole form with h = 1 (checked against direct partial traces by `code/check_wplanes.py`, part (c)).

**Theorem 3.** For every n ≥ 2, all U_1, …, U_n with Σ_s U_s^*U_s = I and all F,

    Σ_{s<t} ‖K_s F_t − K_t F_s‖² ≤ 2 Σ_u ‖F_u‖².

So the co-range sector c = n − 1 of q_n is nonnegative for every range plane in H_1, and the constant is sharp.

*Proof.*
1. Split the reference index: F_u = Σ_r e_r ⊗ F_u[r]. Put w_{s,r} = U_s e_r and let κ_{s,r} contract the letter at copy
   s against w_{s,r}. Then K_s F_t = Σ_r κ_{s,r} F_t[r]. So the left side is ‖𝒟F‖², where 𝒟F = Σ_r 𝒟_r F[r] and
   (𝒟_r G)_{st} = κ_{s,r} G_t − κ_{t,r} G_s.
2. Fix r and drop it from the notation. Put Γ = Σ_t w_t ⊗_t G_t ∈ (C^2)^{⊗n}, where w_t fills the hole of G_t. For
   s ≠ t, ⟨w_s ⊗_s G_s, w_t ⊗_t G_t⟩ = ⟨κ_t G_s, κ_s G_t⟩, so the cross terms of ‖𝒟_r G‖² and ‖Γ‖² cancel:

       ‖𝒟_r G‖² + ‖Γ‖² = Σ_t [ Σ_{s≠t} ‖κ_s G_t‖² + ‖w_t‖² ‖G_t‖² ].

3. Since ‖κ_s G_t‖ ≤ ‖w_s‖ ‖G_t‖ and Σ_s ‖w_s‖² = ⟨e_r| Σ_s A_s |e_r⟩ = 1, this gives ‖𝒟_r G‖² ≤ Σ_t ‖G_t‖².
4. Hence ‖𝒟F‖ ≤ ‖F[1]‖ + ‖F[2]‖ ≤ √2 ‖F‖. ∎

The sector inequality follows, since Ā_u ⪰ 0. It is sharp. Take a real tight frame v_1, …, v_n of R²
(Σ_s v_s v_s^T = I) and unit vectors z_s. Put U_s = |z_s⟩⟨v_s| and F_t = (J v_t) ⊗ ⊗_{j≠t} z_j, with
J = [[0, 1], [−1, 0]]. Then ‖𝒟F‖² = 2‖F‖², and ⟨F_t|Ā_t|F_t⟩ = 0 because J v_t ⊥ v_t. So both sides of the sector
inequality are equal.

**Corollary 4.** If X has rank at most two and range(X) ⊆ H_1, then q_4(X) ≥ 0, whatever the co-range.

*Proof.* By the charge splitting (step 1 of the proof of Theorem 1), q_4(X) = Σ_c q_4(X P_c) over the co-range sectors
c = 0, …, 4. The sectors c ≤ 2 are covered by Theorem 1, c = 3 by Theorem 3, and c = 4 is the free top sector. ∎

Three copies were known (every sector there is c ≤ 2 or c = n). Four copies is the first case where an open sector had
to be closed. For n ≥ 5 the sectors 3 ≤ c ≤ n − 2 remain. The same argument applied to one reference index works in
every sector: with Γ inserting w_t at the holes, the form Σ_ρ [2‖τ_ρ‖² − |Σ_t τ_{ρ,t}|²] is at most ‖G‖². The two
reference indices cannot be recombined as in step 4, because for c ≤ n − 2 the per-string form 2I − J is indefinite.

## Why the sectors c ≥ 3 are different

The even split gives each string the weight 4/c < 2 on its off-diagonal blocks, and the Lemma needs 2. This loss is
real. The rank-two matrix Y = [[I, I], [−I, −I]] has F_c(Y)/‖Y‖² = 2/c − 1, which is negative for every c ≥ 3. So a
proof for c ≥ 3 has to share the off-diagonal weight unevenly, or use the rank-two coupling between different strings.
Theorem 2 does the latter for weighted W-planes. There all range blocks U_s are multiples of unitaries, and monogamy of
a central qubit controls all strings at once. For general range blocks this monogamy bound alone is too weak.

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
- For c = n − 1 the sector inequality takes a compact form, and Theorem 3 proves it.

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
- `code/check_class.py` brute-forces q_n/‖X‖² on the class of Theorem 1 for n = 2, …, 5, with random and adversarial
  starts (minimum +1.6 × 10^{−9}). It also runs the open sector c = 3 at n = 4, 5, as evidence only.
- `code/check_theorem3.py` checks Theorem 3:
  - the identity of step 2 (n = 2, …, 5, letters in C^2 and C^3; error 4 × 10^{−16});
  - λ_max(𝒟^*𝒟) ≤ 2 at random and adversarial U for n = 3, 4, 5 (largest 1.9999999997);
  - the sharp example (ratio exactly 2 for n = 3, 4, 5);
  - the one-reference-index bound in every sector (identity error 4 × 10^{−16}, ratio at most 1);
  - Corollary 4 end to end at n = 4 from the dense operator (minimum over random and minimised planes +3 × 10^{−10}).
- `code/check_wplanes.py` checks Theorem 2:
  - the central-spin value (m + 1)/2 for m ≤ 8;
  - the weighted bound (largest violation 2 × 10^{−15});
  - the hole form against direct partial traces (n = 3, 4, every h; error 10^{−15});
  - the exact minimum over all co-ranges for the W-plane and for weighted W-planes, n = 2, …, 5, with random rotations
    W_s (minimum 2^{1−n} for n ≥ 3).

The logs are in `logs/`.
