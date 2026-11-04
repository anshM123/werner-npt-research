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
- Proposition 5 (added 4 November 2026) settles every sector 1 ≤ c ≤ n − 1 for range planes in which no c copies are
  heavy: Σ_{t∈S} A_t ⪯ 2(n − c)/(3(n − c) − 1) · I for every set S of c copies. It is Theorem 3 applied to every set
  of c + 1 copies.
- Proposition 6 (added 4 November 2026) puts Theorems 1–3 and Proposition 5 into one family of certificates. Each
  certificate is a finite set of matrix inequalities and proves its sector. Every tested POVM for n ≤ 6 has one, but
  a counterexample at n = 7 shows that the family is not enough for all n.
- Proposition 7 (added 4 November 2026) adds a third level, the creation pairing on sets of c + 1 copies. It fixes the
  counterexample (n ≤ 10 tested), and in its balanced form it only asks to split a positive semidefinite part of each
  pattern among the pattern's children. Every POVM we tested has such a split. An explicit rule for the splits would
  close all sectors for all n; we do not have one.
- Proposition 8 (added 4 November 2026): in a second exact form of each sector, every one-copy correction has a
  single positive eigenvalue √det A_t. With the central-spin bound this closes every sector, for every n, when the
  rank-two parts are small compared with the vacuum of every hole set.
- Proposition 9 (added 4 November 2026): when only one copy has rank two, each sector is equivalent to one operator
  inequality on the other copies. It holds for every n when that copy carries a full direction (A_0 has eigenvalue 1),
  which includes the heavy-copy counterexample family, in a quantitative neighbourhood of that case, and whenever
  that copy is heavy in both directions (1/s_1² + 1/s_2² ≤ 2 + 8/h²).
- Proposition 10 (added 4 November 2026): every sector, for every n, when two copies together carry a full direction
  (A_0 + A_1 has eigenvalue 1). The proof is an anchor argument: such a plane contains (two-copy vector) ⊗ (product),
  and an exact computation shows that every such vector, of any Schmidt rank, is an anchor for every n.
- Lemma 11 and Conjecture 12 (added 4 November 2026): a relative version of the anchor argument reduces every sector
  of every single-excitation plane to one positive semidefinite statement about a single weighted W-type vector
  Σ_t w_t |1⟩_t. Copies of weight zero reduce to a finite test on the others (Lemma 13). Theorem 16 (added 4 November
  2026) proves it whenever at most three weights are nonzero, by an exact certificate over Q(a, b, c). So every
  sector is nonnegative for every n whenever A_0 + A_1 + A_2 has the eigenvalue 1, and the same holds for the
  corresponding tangent planes (Corollary 17). Larger supports hold in all numerical tests and at exact special
  points.

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

## Proposition 5: the open sectors when no c copies are heavy

Added 4 November 2026. Theorem 3 also controls the open sectors, by applying it to every set of c + 1 copies. We work
in the gauge Σ_t A_t = I and write the sector form as

    K̃_c = I + Vac − T + ½ α^*α,

where Vac = Σ_t Ā_t ⊗ P_{0,t} (vacuum weight), T = Σ_t β_t^* β_t, α = Σ_t β_t, and β_t contracts the reference index
and the letter at copy t with U_t (β_t F = K_t F in the notation of Theorem 3). So K̃_c ⪰ 0 is the sector inequality
L ≤ 2‖F‖² + 2⟨F|Vac|F⟩ with L = 2T − α^*α. Let Ex = Σ_t Ā_t ⊗ P_{1,t}; on co-range strings with excited set S it is
Ā_S = Σ_{t∈S} Ā_t, and Vac = I − Ex there. Put h = n − c.

**Proposition 5.** For every n, every 1 ≤ c ≤ n − 1 and every single-excitation range plane,

    K̃_c ⪰ (2h I − (3h − 1) Ex)/(h + 1) + (h − 1)/(2(h + 1)) · α^*α.

In particular, sector c is nonnegative whenever λ_max(Σ_{t∈S} A_t) ≤ 2h/(3h − 1) for every set S of c copies.

The threshold is 1 for h = 1 (that is Theorem 3), 4/5 for h = 2, 3/4 for h = 3, and it decreases to 2/3 as h grows.
Item 3 of Proposition 3 in `povm-minimum-principle/` gives the same conclusion under the threshold 2/3 for every h.
So the sufficient condition is strictly weaker here.

*Proof.*
1. For a string ρ with c − 1 letters, let τ_ρ = (τ_{ρ,t}) over its h + 1 free copies, τ_{ρ,t} = (β_t F)(ρ). Then
   L = Σ_ρ τ_ρ^* (2I − J) τ_ρ, with J the all-ones matrix.
2. On C^{h+1}, 2I − J = (2/(h + 1)) Σ_{s<t} (e_s − e_t)(e_s − e_t)^T − ((h − 1)/(h + 1)) J. The J term adds up to
   ‖αF‖², so
       L = (2/(h + 1)) Σ_ρ Σ_{s<t} |τ_{ρ,s} − τ_{ρ,t}|² − ((h − 1)/(h + 1)) ‖αF‖².
3. Group the pair terms by Q = supp(ρ) ∪ {s, t}, a set of c + 1 copies. For fixed Q they are exactly the terms of
   Theorem 3 on the copies of Q: the co-range strings with excited set Q − u play the role of F_u, and the copies
   outside Q are vacuum spectators. So the double sum is Σ_Q ‖𝒟_Q F‖².
4. On Q the range blocks satisfy Σ_{t∈Q} A_t = A_Q ⪯ I. If A_Q is invertible, replace U_t by U_t A_Q^{−1/2} and the
   reference index of F by (Ā_Q)^{1/2}. Theorem 3 then gives ‖𝒟_Q F‖² ≤ 2 Σ_{u∈Q} ⟨F_{Q−u}| Ā_Q |F_{Q−u}⟩, and by
   continuity this holds for every A_Q.
5. A string with excited set S lies in h of these sets, Q = S ∪ {u} with u a hole. Summing,
   Σ_Q ‖𝒟_Q F‖² ≤ 2 Σ_S ⟨F_S| hĀ_S + Ā_{H(S)} |F_S⟩ = 2 Σ_S ⟨F_S| I + (h − 1) Ā_S |F_S⟩, using Ā_H = I − Ā_S.
6. Hence L ≤ (4/(h + 1)) ⟨F| I + (h − 1) Ex |F⟩ − ((h − 1)/(h + 1)) ‖αF‖². Insert this into
   2K̃_c = 2I + 2 Vac − L = 4I − 2 Ex − L. ∎

The bound is sharp: at h = 1 it is Theorem 3, and the adversarial searches of `code/check_subsystem_bound.py` reach
zero in every tested sector. It does not cover all planes. If some c copies carry almost all of the weight in one
direction, the threshold fails, and so does the bound. An example is one heavy copy with A_0 = diag(1 − d, 1) and light
copies elsewhere. The sector stays nonnegative there numerically, but for a different reason (see the remarks below).

**Remarks (proved identities and numerical facts).**
- *Hard-core boson form.* Let a_r^† = Σ_t |U_t e_r⟩⟨0|_t create the range vector u_r and let G^r be the components of
  F. Then ⟨F|K̃_c|F⟩ = ‖F‖² + Σ_{r,r'} ⟨a_{r'}^† G^r, a_r^† G^{r'}⟩ − ½‖αF‖². Equivalently,
  K̃_c = I + Σ_{r,r'} |r⟩⟨r'| ⊗ [a_{r'}, a_r^†] + ½ α^*α, where the commutators are sums of one-copy terms.
- *Twisted identity.* With F̃ = (G^2, −G^1) and γ_t F = Σ_r |U_t e_r⟩⟨0|_t G^r (creation at copy t),
  L = (Σ_t Tr A_t)‖F‖² − s − Σ_t ‖γ_t F̃‖² − ‖Θ_+‖² + ½ ‖Σ_t γ_t F̃‖². Here s ≥ 0 is the one-copy Cauchy–Schwarz slack,
  ‖Θ_+‖² = ‖a_1^† G^1‖² + ‖a_2^† G^2‖² + ½‖a_1^† G^2 + a_2^† G^1‖², and the gauge is not used. For c = 1 every
  created pair receives exactly two terms, so ½‖Σ_t γ_t F̃‖² ≤ Σ_t ‖γ_t F̃‖² and L ≤ 2‖F‖². This is a short proof of the
  sector c = 1 (the content of the Lemma there).
- *The vacuum terms are needed.* The bound L ≤ 2‖F‖² holds for c = 1 and c = n − 1, but it fails in every sector
  2 ≤ c ≤ n − 2. The heavy-copy example above gives λ_min of I − T + ½α^*α from −7.8 × 10^{−3} to −1.6 × 10^{−2} at
  n = 5, 6, while K̃_c itself stays positive.
- *No low-degree certificate* (numerical). For a fixed co-range, the sector form is a Hermitian form in the range
  blocks. We asked for a sum-of-squares certificate of degree two in both, modulo the gauge. Such certificates exist for
  c = 1, c = n − 1 and c = n at n = 3, 4. None was found for n = 4, c = 2 (Theorem 1's sector) or n = 5, c = 3: an
  alternating-projection solver stalls at a relative residual of 3 × 10^{−2}. So the remaining sectors seem to need an
  argument of the Lemma's kind, not only identities and Cauchy–Schwarz.
- *Splitting the budget per string* (numerical). One could split each column's budget 2I + 2A_H among the strings ρ
  that contain it and apply the Lemma to each X^ρ. With splits diagonal in the letters, this provably fails at the
  uniform W-plane whenever c > 2n/3: by symmetry the even split is optimal there, and it gives Riesz weights
  2(n − h)/(n + h) > 1 on every block. It also fails numerically near POVMs with nearly unused copies (n = 6, c = 4). With
  metrics that mix the letters, the same per-string decomposition is feasible at the uniform W-plane for n = 5, c = 4
  (margin 0.1). The next section builds on this.

The sectors 3 ≤ c ≤ n − 2 therefore remain open for planes with a heavy set of c copies.

## Proposition 6: one certificate family for all three mechanisms

The proofs above use three different mechanisms:
- Theorem 1 uses the rank-two Lemma on each string ρ;
- Theorem 2 uses monogamy of the reference qubit;
- Theorem 3 uses a Lagrange identity on the whole sector.

All three fit into one family of certificates. Numerically this family suffices at every POVM we tested for
n ≤ 6, but not for all n: it fails at n = 7 (see the counterexample below).

Write the sector form with co-range vectors f_σ ∈ C², in the gauge Σ_t A_t = I:

    Q_c(f) = Σ_σ f_σ^* (2I + 2A_{H(σ)}) f_σ − Σ_ρ [2 Σ_t |τ_{ρ,t}|² − |Σ_t τ_{ρ,t}|²],

where τ_{ρ,t} = Σ_{z,r} U_t[z, r] f_{ρ+z@t}[r] and y_t ∈ C⁴ has entries conj(U_t[z, r]).

A **certificate** for sector c consists of two families.
- Numbers θ_Q ≥ 0, one for each set Q of c + 1 copies.
- Hermitian 4 × 4 matrices m_{S,t}, one for each set S of c copies and each t ∈ S. Each acts on the letter at copy t
  and the reference index; m^{(t)} denotes it extended by the identity on the other letters.

They must satisfy two conditions.
- **(L)** For every set R of c − 1 copies, with free copies F(R): ⊕_{t∈F(R)} m_{R∪t,t} ⪰ Y(N_R). Here N_R is the matrix on
  C^{F(R)} given by N_R = 2I − J − Σ_{s<t} θ_{R∪{s,t}} (e_s − e_t)(e_s − e_t)^T, and Y(N) has 4 × 4 blocks N_st |y_s⟩⟨y_t|.
- **(P)** For every set S of c copies:
  Σ_{t∈S} m_{S,t}^{(t)} + 2 Σ_{u∉S} θ_{S∪u} I ⊗ A_{S∪u} ⪯ I ⊗ (2I + 2A_{H(S)}).

**Proposition 6.** If sector c has a certificate, then Q_c ≥ 0 on that sector.

*Proof.* For a string ρ with support R, split 2I − J = 𝓛_R + N_R, where 𝓛_R is the weighted Laplacian with pair weights
θ_{R∪{s,t}}. As in steps 3–5 of Proposition 5, the 𝓛_R parts add up to Σ_Q θ_Q ‖𝒟_Q f‖². Theorem 3 on each Q bounds
this sum by 2 Σ_Q θ_Q Σ_{u∈Q} f_{Q−u}^* A_Q f_{Q−u}.

The N_R part of the string ρ is ⟨v_ρ| Y(N_R) |v_ρ⟩, where v_ρ collects the columns ρ + z@t. By (L) it is at most
Σ_t ⟨v_ρ| m_{R∪t,t} |v_ρ⟩. Summed over ρ and regrouped by columns, this is Σ_S Σ_{t∈S} f_S^* m_{S,t}^{(t)} f_S.

By (P), the two bounds together are at most Σ_σ f_σ^* (2I + 2A_H) f_σ. ∎

The family contains the earlier arguments.
- θ = 1 and m = 0 is Theorem 3 (c = n − 1).
- θ_Q = 2/(h + 1) and m = 0 is Proposition 5.
- θ = 0 and m_{S,t} = 2|y_t⟩⟨y_t| works whenever every A_t has rank at most one. Then (L) holds with slack |⊕y⟩⟨⊕y|,
  and (P) follows from 2A_S ⪯ 2I + 2A_H, because each y_t is a product of a letter and a reference vector.
- For weighted W-planes, the same choice reduces (P) to the central-spin bound. This is the proof of Theorem 2.

With θ = 0, (P) is a monogamy condition on the reference qubit and (L) is the Lemma's condition on each string.

**Conjecture 6 (withdrawn for n ≥ 7).** We first conjectured that every co-range sector of every single-excitation
range plane has a certificate. This is false for n ≥ 7 (counterexample below). The question remains open for n ≤ 6.
There, a certificate in the sector c = 3 at n = 5 for every POVM would give the endpoint on all single-excitation
range planes at n = 5, with any co-range.

*Evidence.* We solved the certificate as a semidefinite program, maximising a common slack, and re-checked every
inequality of the returned solution with independent eigenvalue computations
(`code/check_combined_certificate.py`). The tested cases were n = 5, 6, sectors 3 ≤ c ≤ n − 1, at W-planes, at one
or two heavy copies, at random POVMs, at near rank-one frames and with two nearly unused copies. The slack was positive
every time, ranging from 3 × 10^{−5} to 0.42.

A sweep over 60 structured POVMs at n = 5, c = 3 (random rank-one frames with noise, nearly unused copies,
nearly parallel directions, one or two heavy copies, near W-planes) found worst slack +5.5 × 10^{−5}
(`code/sweep_combined_certificate.py`). The smallest margins occur next to nearly unused copies.

Neither half of the certificate works alone. With θ = 0 the program is infeasible at c = n − 1, n = 6, for one heavy copy
A_0 = diag(1 − d, 1) and five light copies (d/5)|e_1⟩⟨e_1|. A rough count with symmetric, rank-one pieces explains why.
The light copies' pieces are products, so they add up with no monogamy gain, and the count closes only for at most four
light copies (slack ≈ 0 at n = 5, negative at n = 6). This explains the numbers but is not a proof of infeasibility.
With m = 0 the program is Proposition 5. With splits diagonal in the letters it fails as shown above.

**Counterexample (n = 7, c = 5).** Take one heavy copy A_0 = diag(1 − d, 1) and six light copies U_t = √(d/6)|1⟩⟨e_1|.
The program is convex and invariant under permuting the light copies, so a certificate exists if and only if a
symmetric one does, and the symmetric program has five unknown blocks. Its maximal slack is −6.69 × 10^{−5} at
d = 0.2. The full program gives the same value (−6.688 × 10^{−5}, re-checked independently). At n = 8, c = 6 the
slack is −4.5 × 10^{−4} already at d = 0.05.

The sector form itself is positive at the n = 7 point: λ_min = +0.140. So the counterexample is to the
certificate family, not to the endpoint. Along this family the slack falls steadily with n (d = 0.05: +4.1 × 10^{−3},
+1.4 × 10^{−3}, +1.8 × 10^{−4} at n = 5, 6, 7). The light copies are rank-one products along one direction, so their
pieces add up with no monogamy gain. Theorem 3 on sets of c + 1 copies also does not see the light copies
collectively. A proof of the open sectors for all n needs a further, collective ingredient for many light copies;
checks in `code/check_certificate_counterexample.py`.

A first-order computation locates that ingredient. Here it is derived by degenerate perturbation theory and checked
numerically, not written as a proof. Write N = n − 1 and treat the light copies as spins (letter 1 = up, hole = down).
- To first order in d, the minimum of the sector form along the family is d(1/2 + 2/N). Extrapolated numerics for
  n = 5, …, 9 give 1, 0.9, 0.833, 0.786, 0.75 for min/d, matching this to 10^{−5}.
- Two terms compete. One is the cost of the strings where the heavy copy is excited, K = 2N_up − S⁺S⁻. The other comes
  from eliminating the strings where the heavy copy is a hole, W = S⁻S⁺. Their sum is exactly K + W = N, because
  [S⁺, S⁻] = 2S_z.
- Bounded separately, they have maxima N and 2N − 2, reached on different spin multiplets. A certificate that bounds
  them separately therefore certifies at most d(6 − N)/(2N). That is positive for n ≤ 6, zero at n = 7 and negative from
  n = 8, which matches the slacks above.

The missing ingredient is thus a collective inequality for the light copies: the commutator identity
[S⁺, S⁻] = 2S_z, used across all strings at once.

So Proposition 6 remains a sufficient condition that we checked at many POVMs for n ≤ 6, and no more.

## Proposition 7: adding the creation level

Added 4 November 2026. We use the notation of the hard-core form (remarks after Proposition 5). The certificates of
Proposition 6 never use the creation pairing ⟨Θ, SΘ⟩ = Σ_{r,r'} ⟨a_{r'}^† G^r, a_r^† G^{r'}⟩: Theorem 3 throws it
away as slack. The counterexample needs exactly this pairing, since its term W = S⁻S⁺ is a creation term. So we add a
third level to the certificate.

Call a set S' of c + 1 copies a *child* of the sets S' − u, u ∈ S'. Let Π_{S'}(F) be the part of ⟨Θ, SΘ⟩ that lands
on strings with excited set S'. It involves only the co-range on the c + 1 sets S' − u.

**An exact family.** The identity T − α^*α = Vac − ⟨Θ, SΘ⟩ of the hard-core form gives, for every real μ,

    2K̃_c = Σ_S F_S^* (2I + (2 − 2μ) Ā_H) F_S − Σ_ρ τ_ρ^* [(2 − 2μ) I − (1 − 2μ) J] τ_ρ + 2μ Σ_{S'} Π_{S'}(F).

Here S runs over sets of c copies with holes H = H(S), ρ over strings with c − 1 letters, and τ_ρ is as in
Proposition 5. The value μ = 0 is the original form and μ = 1 the hard-core form.

**Three-level certificate.** Besides the data of Proposition 6 (θ_Q ≥ 0 and m_{S,t}), it has a real number μ and a
Hermitian matrix k_{S,u} for each set S of c copies and each hole u of S, acting on the strings with excited set S.
The conditions are:
- (L) as in Proposition 6, with 2I − J replaced by (2 − 2μ)I − (1 − 2μ)J;
- (C) for every child S': 2μ Π_{S'} + ⊕_{u∈S'} k_{S'−u,u} ⪰ 0;
- (P) as in Proposition 6, with Σ_{u∉S} k_{S,u} added on the left and 2A_H replaced by (2 − 2μ)A_H on the right.

**Proposition 7.** If sector c has a three-level certificate, then K̃_c ⪰ 0. If every condition holds with a common
slack λ, then λ_min(2K̃_c) ≥ (n + 1)λ.

*Proof.* Start from the exact family. The first two sums are handled as in the proof of Proposition 6, with the new
N_R. By (P), what is left over is at least Σ_S Σ_{u∉S} F_S^* k_{S,u} F_S. Regroup these terms by children (the pair
(S, u) belongs to the child S ∪ u) and add 2μ Π_{S'}. By (C) each child gives a nonnegative total. For the slack
count, each column lies in one set S, in c strings ρ and in h children, so it collects λ from n + 1 conditions. ∎

**The balanced case μ = 1/2.** Then the J terms cancel:

    2K̃_c = Σ_S F_S^* R_S F_S + Σ_{S'} Π_{S'}(F),   R_S = 2I + Ā_H − Σ_{t∈S} β_t^*β_t.

The pattern part R_S is positive semidefinite, because ‖Σ_{t∈S} β_t^*β_t‖ ≤ Σ_{t∈S} Tr A_t ≤ 2. Take θ = 0 and the
smallest metrics that (L) allows at μ = 1/2, m_{S,t} = |y_t⟩⟨y_t|. Then a three-level certificate is just a way to
split R_S among the h children of S:

    Σ_{u∉S} k_{S,u} ⪯ R_S   for every S,      Π_{S'} + ⊕_{u∈S'} k_{S'−u,u} ⪰ 0   for every child S'.

We call this a *balanced certificate*. For c = n − 1 there is a single child, and a balanced certificate is the
sector inequality itself (Theorem 3). The content is in h ≥ 2.

Each child form has a simple structure. For a string s' with excited set S' and letters z_v, put

    N(s') = Σ_{v∈S'} F(s' − v) y_{v,z_v}^T,

a 2 × 2 matrix. Here F(s' − v) ∈ C² is the co-range at the string s' with copy v emptied, and y_{v,z} ∈ C² is row z
of U_v. Then Π_{S'} = Σ_{s'} Tr(N N̄) = Σ_{s'} (‖N_sym‖² − ‖N_anti‖²). So on each child string the negative part has
rank one: it is ½|Σ_{v∈S'} det(F(s' − v), y_{v,z_v})|². Summed over all children, these parts give
½‖a_1^† G^2 − a_2^† G^1‖², the creation of the twisted co-range (G^2, −G^1) that also appears in the twisted
identity of Proposition 5.

*Evidence* (`code/threelevel_certificate.py`, `code/check_threelevel_certificate.py`). Every solution below is
re-checked with independent eigenvalue computations, together with the guard λ_min(2K̃_c) ≥ (n + 1)λ.
- *The heavy-copy counterexample.* Restrict to the strings where every excited light copy carries letter 1. This
  subspace is invariant and contains the minimiser (the restricted minima equal the full ones). There a balanced
  certificate exists for n = 5, …, 10, with slack +3 × 10^{−3} to +1.6 × 10^{−2}; μ = 1/2 and θ = 0 are optimal. The
  certificates of Proposition 6 fail there from n = 7 on.
- *General real POVMs at n = 5, c = 3* (random, near rank-one frame, two nearly unused copies, two heavy copies, one
  heavy copy, W-plane). A balanced certificate exists in every case, with slack 1.3 × 10^{−3} to 0.17. With μ free the
  slack is somewhat larger at random and near rank-one frames (optimal μ between 0.02 and 0.2).
- *The splits must be matrices.* At a near rank-one real frame (n = 5, c = 3), splits k_{S,u} = w_{S,u} R_S with the
  best scalar weights fail, while a matrix-valued balanced certificate exists.
- *Three numbers per split, at n ≤ 6 with c = 3.* Put E_S = 2Ā_S − Σ_{t∈S} β_t^*β_t, which is positive
  semidefinite (a sum of the one-copy terms 2Ā_t ⊗ I − β_t^*β_t). A balanced certificate of the form
  k_{S,u} = x_{S,u} Ā_u/Tr A_u + w_{S,u} E_S + y_{S,u} I exists at all six POVM types at n = 5, c = 3, at all six at
  n = 6, c = 3, and on the heavy-copy family for n = 5, …, 10. The child slack ranges from 7.6 × 10^{−4} to 0.5
  (`code/check_structured_split.py`; each solution is rebuilt from its numbers, re-checked, and repaired for solver
  tolerance by a shift of at most 6 × 10^{−7}). At n = 6, c = 4 the three-number family works at six of the seven
  test POVMs (random, two nearly unused copies, two heavy, one heavy, the heavy-copy family, the uniform W-plane;
  first-order solver SCS, slacks 1.1 × 10^{−3} to 0.47; the re-check shows solver errors up to 6.5 × 10^{−3}, below
  each slack). At the near rank-one frame it fails, with slack
  −2.6 × 10^{−2}, while that sector itself is strictly positive. A general matrix-valued balanced certificate at
  that frame was not solved: the SCS run was stopped after two hours. So at c = 4 the three numbers are not enough.
  The numbers vary from pair to pair. Closed-form choices we tried
  (even or trace-weighted splits of R_S, of E_S or of 2I + 2Ā_H − Σβ^*β, with or without moving vacuum between
  children, matrix-weighted splits) each fail at some tested POVM; for example the trace-weighted split of
  2I + 2Ā_H − Σβ^*β leaves a child eigenvalue −8.0 × 10^{−2} at n = 6, c = 4.

Does every sector have a balanced certificate? If so, the endpoint would hold on all single-excitation range planes,
for every n and every co-range. We know no counterexample. Conjecture 6 showed how evidence at small n can mislead,
so we state no conjecture. What is still missing is an explicit rule for the splits.

## Proposition 8: small rank-two parts

Added 4 November 2026. The sector form has a second exact expression. With Ex = Σ_t Ā_t ⊗ P_{1,t}, so that
I = Vac + Ex, the sector form 2K̃_c = 2I + 2 Vac − 2T + α^*α becomes

    2K̃_c = 4 Vac + α^*α − 2 Σ_t Δ_t^{exc},   Δ_t = β_t^*β_t − Ā_t ⊗ I.

Here Δ_t acts on the letter at copy t and the reference index, on the strings where t is excited. Write the singular
value decomposition U_t = Σ_k s_k |f_k⟩⟨g_k| with s_1 ≥ s_2. In the basis ḡ_k ⊗ f_l,

    Δ_t = s_1 s_2 (Φ_t^+ − Φ_t^−) − s_1² |ḡ_1 f_2⟩⟨ḡ_1 f_2| − s_2² |ḡ_2 f_1⟩⟨ḡ_2 f_1|,

with Φ_t^± the projectors onto (ḡ_1 ⊗ f_1 ± ḡ_2 ⊗ f_2)/√2. So Δ_t has exactly one positive eigenvalue, √det A_t, and
its eigenvector is maximally entangled between the letter qubit and the reference qubit. If A_t has rank one, then
Δ_t ⪯ 0. So the all-rank-one case is immediate: 2K̃_c ⪰ 4 Vac + α^*α.

**Proposition 8.** For every n and every 1 ≤ c ≤ n − 1, sector c is nonnegative if, for every set S of c copies,

    4 λ_min(Σ_{u∉S} A_u) ≥ Σ_{t∈S} √det A_t + max_{t∈S} √det A_t.

*Proof.* Drop α^*α ⪰ 0. The other two terms act on each pattern separately. On the strings with excited set S,
Vac = Ā_{H(S)} ⪰ λ_min(A_{H(S)}) I. The projectors Φ_t^+ are onto maximally entangled states of the letter qubit t
and the reference qubit, so the weighted central-spin Lemma of Theorem 2 gives

    2 Σ_{t∈S} Δ_t ⪯ 2 Σ_{t∈S} √det A_t Φ_t^+ ⪯ (Σ_{t∈S} √det A_t + max_{t∈S} √det A_t) I. ∎

The condition asks that the rank-two parts be small compared with the vacuum of every hole set. It is not implied
by Proposition 5. At n = 5, c = 3 there is a POVM with every A_t of rank two, a set of three copies with
λ_max(A_S) = 0.986 > 4/5, and margin 3.8 × 10^{−3} in the condition (`code/check_small_ranktwo.py`, part (e)). The
condition fails as soon as the holes of some S are parallel rank-one copies while S contains a rank-two copy, as in
the heavy-copy family. There the term α^*α, dropped here, is what keeps the sector positive.

## Proposition 9: one copy of rank two

Added 4 November 2026. Suppose copy 0 has a rank-two range block and the other N = n − 1 copies have rank one. After
rotations of the reference index and of the letters, U_0 = diag(s_1, s_2) and U_t = |1⟩⟨x̄_t| for t ≥ 1, with
Σ_t x_t x_t^* = I − Ā_0. A rank-one copy uses one letter. The strings where it carries the other letter form an
invariant block. On that block the form is 2⟨F, Ā_t F⟩ plus the same sector problem for the other copies, with one
excitation fewer. So in (a), (c) and (d) below we assume that every rank-one copy carries only its range letter.

Put y_t = Ā_0^{−1/2} x_t; every family of vectors y_t ∈ C² arises this way. Treat the rank-one copies as qubits. Let
H_k be the functions on k-element sets of them, σ_t^− the annihilation at copy t, and n_t the projector onto "copy t
is a hole". On V ⊗ H_k define

    X'_k = Σ_t y_t y_t^* ⊗ n_t,    b_k = Σ_t y_t^* ⊗ σ_t^− : V ⊗ H_k → H_{k−1},
    M = 4 X'_{c−1} + b_{c−1}^* b_{c−1},    K = ¼ b_c (I + X'_c)^{−1} b_c^*.

**Proposition 9.**
- (a) Sector c is nonnegative if and only if Tr_V (2 + M)^{−1} ⪯ ½ I + (2 + 4K)^{−1} on H_{c−1}.
- (b) If A_0 has an eigenvalue equal to 1, then every sector is nonnegative, for every n. In that case every other A_t
  is a multiple of one rank-one projector, so this needs no further hypothesis.
- (c) Sector c is nonnegative whenever ν(h + (c + 3)² μ/2) ≤ μ, where μ = s_min^{−2} − 1, ν = s_max^{−2} − 1 and
  h = n − c, whatever the directions x_t.
- (d) Sector c is nonnegative whenever 1/s_1² + 1/s_2² ≤ 2 + 8/h², whatever the directions x_t. At n = 5, c = 3
  it is enough that s_1², s_2² ≥ ½: copy 0 is heavy in both directions.

Part (b) contains the heavy-copy family A_0 = diag(1 − d, 1), with any weights and phases on the light copies. That
is the family where the certificates of Proposition 6 fail from n = 7 on.

*Proof.*
- (a) Start from the identity of Proposition 8. There Δ_t = 0 for the light copies, so
  2K̃_c = 4 Vac + α^*α − 2Δ_0^{exc}. Conjugate the reference index by Ā_0^{−1/2}. Write the strings where copy 0 is
  excited as F_h = (I ⊗ Ā_0^{−1/2})G, with G = (G_1, G_2) indexed by the letter of copy 0, and the strings where copy 0
  is a hole as F_l = Ā_0^{−1/2}H. Then

      ⟨F, 2K̃_c F⟩ = 2‖G‖² − 2‖g‖² + Σ_z ⟨G_z, M G_z⟩ + 4⟨H, (I + X'_c)H⟩ + ‖g + b_c H‖²,   g = G_1[1] + G_2[2].

  The minimum over H is ⟨g, (I + K)^{−1} g⟩. For fixed g, the minimum of Σ_z ⟨G_z, (2 + M)G_z⟩ is ⟨g, R^{−1} g⟩
  with R = Tr_V (2 + M)^{−1}. So the sector is nonnegative iff R^{−1} − 2 + (I + K)^{−1} ⪰ 0, that is, iff
  R ⪯ (I + K)(I + 2K)^{−1} = ½ I + (2 + 4K)^{−1}.
- (b) If s_max = 1, then Λ = Σ_t y_t y_t^* = Ā_0^{−1} − I has rank one, so y_t = η_t e for one unit vector e. Then
  M = ee^* ⊗ m with m = 4D + B^*B, where D = Σ_t |η_t|² n_t and B = Σ_t η̄_t σ_t^−, and K = ¼ B (1 + D_c)^{−1} B^*. In
  the basis (e, e^⊥), Tr_V (2 + M)^{−1} = ½ + (2 + m)^{−1}, so (a) asks for m ⪰ 4K. On H_{c−1} the hard-core
  commutator is
      B B^* − B^* B = Σ_t |η_t|² (n_t − (1 − n_t)),
  because terms at different copies commute. Hence 4K ⪯ BB^* = B^*B + D − U with U = Σ_t |η_t|² (1 − n_t), and
  m − 4K ⪰ 3D + U ⪰ (Σ_t |η_t|²) I. For general data, induct on n through the invariant blocks above. Removing a
  light copy and renormalising keeps the eigenvalue 1 of A_0, because the removed A_t is orthogonal to its
  eigenvector.
- (c) Work in the eigenbasis (e_1, e_2) of Λ, with Λ_11 = μ ≥ Λ_22 = ν, and write M in 2 × 2 blocks.
  - Tr_V (2 + M)^{−1} = (2 + S_1)^{−1} + (2 + S_2)^{−1} with Schur complements S_1, S_2 ⪰ 0. Since
    (2 + S_2)^{−1} ⪯ ½, it suffices that S_1 ⪰ 4K.
  - Next, S_1 ⪰ M_11 − ½ M_12 M_21, and 4K ⪯ b_1 b_1^* + b_2 b_2^*, where b_r is the e_r part of b_c.
  - The commutator in direction e_1 gives M_11 − 4K ⪰ Λ_11 − b_2 b_2^*.
  - Finally b_2 b_2^* ⪯ hΛ_22 and ‖M_12‖ ≤ (c + 3)(Λ_11 Λ_22)^{1/2}, both by Cauchy–Schwarz over the created or
    annihilated copy.
- (d) Put g(x) = x/(2 + x), which is operator monotone. Since Tr_V (2 + Y)^{−1} = I − ½ Tr_V g(Y) and
  ½ + (2 + Z)^{−1} = I − ½ g(Z), the inequality in (a) reads Tr_V g(M) ⪰ g(4K). It follows from
  Tr_V g(N) ⪰ g(W), where N = 2X'_{c−1} + b_{c−1}^*b_{c−1} ⪯ M and W = b_c b_c^* ⪰ 4K.
  - The hard-core commutator in matrix form gives N = Λ ⊗ I + C^Γ, where C = 𝐛𝐛^* ⪰ 0 with 𝐛 = Σ_t ȳ_t ⊗ σ_t^−,
    Γ is the partial transpose on V, and Tr_V C = W.
  - With κ = ‖N‖, the chord bound g(x) ≥ x/(2 + κ) on [0, κ] gives Tr_V g(N) ⪰ (W + τ)/(2 + κ), where
    τ = Tr Λ = s_1^{−2} + s_2^{−2} − 2.
  - On a qubit, C^Γ = I ⊗ W − (J^{−1} ⊗ I) C (J ⊗ I) ⪯ I ⊗ W, so κ ≤ λ_max(Λ) + ‖W‖.
  - For each eigenvalue w of W, (w + τ)/(2 + κ) ≥ w/(2 + w) is the condition w² − (κ − τ)w + 2τ ≥ 0. It holds when
    (‖W‖ − λ_min(Λ))² ≤ 8τ, and ‖W‖ ≤ hτ by Cauchy–Schwarz. So τ ≤ 8/h² is enough. ∎

Part (d) also covers configurations missed by Propositions 5 and 8 and by (b) and (c), for example s_1² = 0.9,
s_2² = 0.5 at n = 5, c = 3. The same identities reduce the sector with h holes to a statement about one positive
operator C on C² ⊗ H and a 2 × 2 matrix Λ ⪰ 0: if Λ ⊗ I + C^Γ ⪰ 0 and ‖C‖ ≤ h Tr Λ, then
Tr_V g(Λ ⊗ I + C^Γ) ⪰ g(Tr_V C). This abstract statement is false once ‖C‖ may reach 2.82 Tr Λ: there are
instances with dim H = 3, certified in 60-digit arithmetic, so it cannot give h ≥ 3. (An earlier version of this
page said that no violation was found for h ≤ 8; those searches had only drifted to the trivial point C, Λ → 0.)
For h = 2, about 3000 local searches with dim H ≤ 12 found no violation. The physical data have more structure
than these hypotheses, so the counterexamples say nothing about the sectors themselves.

*The elimination works for any other copies.* The first step of the proof of (a) never uses that the other copies
have rank one. Let 2K̃'_k be the sector forms of the other n − 1 copies after the same change of reference index (blocks
U_t A_0^{−1/2}; these are not normalised, since Σ_t A'_t = A_0^{−1} − I, so 2K̃'_k = 2Ā'_tot + 2 Vac' − 2T' + α'^*α'), and
let α' be their annihilation from c to c − 1 excitations. The same congruence gives

    ⟨F, 2K̃_c F⟩ = 2‖G‖² − 2‖g‖² + Σ_z ⟨G_z, 2K̃'_{c−1} G_z⟩ + ⟨H, (4 + 2K̃'_c − α'^*α')H⟩ + ‖g + α'H‖².

Eliminating H and then G as in (a):

- (a′) If some copy has rank two and the other copies satisfy 2K̃'_{c−1} ⪰ 0 and 2K̃'_c ⪰ 0, then sector c is
  nonnegative if and only if

      Tr_V (2 + 2K̃'_{c−1})^{−1} + α' (4 + 2K̃'_c + α'^*α')^{−1} α'^* ⪯ I   on H'_{c−1}.

The hypotheses are the same statement for n − 1 copies, so (a′) turns the sector inequality into an induction
over copies: one rank-two copy is removed at a time, at the price of a quantitative condition on two adjacent sectors
of the rest. For rank-one other copies, (a′) is (a). We have no proof of (a′) in general (`code/check_one_ranktwo.py`,
part (e), checks the congruence at random POVMs with every copy of rank two).

Part (c) covers the opposite end from Proposition 5: copy 0 nearly saturates one direction. For example, at n = 5,
c = 3 with s_min² = ½ it needs s_max² ≥ 20/21, while Proposition 5 needs λ_max(A_S) ≤ 4/5.

For general light directions the case is exactly the single inequality (a). Random and adversarial searches up to
N = 7 found no violation, but we have no proof. The inequality uses dim V = 2: its analogue with a d-dimensional
reference space fails for d ≥ 3 as y → 0, where the left side tends to d/2 and the right side to 1. A separable split
M ⪰ ρ ⊗ 4K (ρ a density matrix on V) and operator subadditivity of (2 + ·)^{−1} both fail on explicit instances
(c = 3, N = 3 and 4), so a proof has to use the hard-core structure.

## Proposition 10: two copies carrying a full direction

Added 4 November 2026. Proposition 9 (b) needs one copy with A_0 having the eigenvalue 1. Here two copies share that
role. The proof does not use the sector forms. It goes back to the anchor method of `support-and-projector-notes/`
(the Schur complement argument of Wu and Zou, Sections 1 and 5 there).

Notation. On one copy let M = I − ½|Ω⟩⟨Ω| on C³_row ⊗ C³_col, with Ω = Σ_a |aa⟩, so that
q_n(X) = ⟨vec X, M^{⊗n} vec X⟩. Let Γ be the transpose on the row factor. Then A = (M^{⊗n})^Γ = (I − ½F)^{⊗n}, where F
swaps row and column at each copy, so A ≻ 0. For a row vector p put
- V_p = p ⊗ I;
- D_p = V_p^* M^{⊗n} V_p;
- R_p = M^{⊗n} V_p D_p^{−1} V_p^* M^{⊗n};
- T_p = A^{−1/2} R_p^Γ A^{−1/2}.

D_p is the form of the rank-one operators |p⟩⟨v|, so D_p ⪰ 2^{−n}‖p‖² I.

**Anchor lemma.** If T_p ⪯ I, then q_n(X) ≥ 0 for every X whose range lies in a plane containing p, whatever the
co-range.

*Proof.* Write X = |p⟩⟨v| + |u⟩⟨w|. Minimising q_n over v leaves (u ⊗ w̄)^*(M^{⊗n} − R_p)(u ⊗ w̄). For product vectors
(x ⊗ y)^* W (x ⊗ y) = (x̄ ⊗ y)^* W^Γ (x̄ ⊗ y), so this equals (ū ⊗ w̄)^*(A − R_p^Γ)(ū ⊗ w̄). It is ≥ 0 when T_p ⪯ I. ∎

All four objects factor over the copies, so T_{p′⊗p″} = T_{p′} ⊗ T_{p″}. For the vacuum on one copy,
spec T_{|0⟩} = {1, −1/3}.

**Proposition 10.**
- (a) Let φ be any vector on two copies and b a product vector on the other n − 2 copies. Then q_n(X) ≥ 0 for
  every X whose range lies in a plane containing φ ⊗ b, with any co-range, for every n.
- (b) Take a single-excitation range plane in the gauge Σ_t A_t = I. If A_0 + A_1 has the eigenvalue 1, then
  q_n ≥ 0 for every co-range, so every sector is nonnegative, for every n.

*Proof.*
- (a) Local unitaries turn φ ⊗ b into a multiple of φ_{β,γ} ⊗ |0⟩^{⊗(n−2)}, where φ_{β,γ} = |00⟩ + β|11⟩ + γ|22⟩
  is the Schmidt form of φ (β, γ ≥ 0). Local unitaries leave q_n unchanged, and the co-range is arbitrary anyway.
  - For all β, γ ≥ 0 we have −3 ⪯ T_{φ_{β,γ}} ⪯ I (and D ≻ 0, by the rank-one bound). This is an exact
    computation, described below; the boundary β = 0 or γ = 0 follows by continuity.
  - Since every one-copy vector is locally unitarily equivalent to |0⟩, T_p is unitarily equivalent to
    T_{φ_{β,γ}} ⊗ T_{|0⟩}^{⊗(n−2)}. Its eigenvalues are λ(−1/3)^j with λ ∈ [−3, 1] and j ≥ 0. All of them are at
    most 1.
  - The anchor lemma applies.
- (b) Let g be a unit eigenvector for the eigenvalue 1.
  - Then Σ_{t≥2} g^*A_t g = 0, so U_t g = 0 for t ≥ 2.
  - The plane therefore contains u_g = (|U_0 g⟩_0 |0⟩_1 + |0⟩_0 |U_1 g⟩_1) ⊗ |0…0⟩. Its two-copy factor has Schmidt
    rank at most two, because its copy-0 marginal lives on span{|0⟩, U_0 g}.
  - Apply (a). ∎

*The exact computation* (`code/check_pair_anchor_rank3_exact.py`, a few seconds). The matrices A − R^Γ and
R^Γ + 3A on the 81-dimensional space have entries in Q(β, γ). The script splits them into the connected components
of their nonzero pattern (19 components, the largest of size 15). In each component it runs symmetric Gaussian
elimination over Q(β, γ). Every pivot it uses is a ratio of two polynomials whose coefficients all have the same
sign, so it is positive for all β, γ > 0. The ranks are 72 and 81, so the kernel of A − R^Γ is exactly
φ̄ ⊗ (C³)^{⊗2}. The script rejects the perturbed condition with R^Γ replaced by 1.001 R^Γ.

This contains Schmidt rank one (the product-ray theorem), Schmidt rank two, and the maximally entangled pair
(`support-and-projector-notes/`). An earlier version of this section treated Schmidt rank at most two with the
one-parameter family φ_u = |10⟩ + u|01⟩. That computation, by exact real-root isolation of the elementary
symmetric functions of 27 torus blocks over Q(u), is kept as `code/check_pair_anchor_exact.py`.

Remarks.
- (b) contains Proposition 9 (b). There U_1 g = 0, so the anchor is a product vector. The bath of the other copies
  never enters: its weights, letters and phases are arbitrary.
- The bound T ⪯ I is tight: λ_max(T_φ) = 1 for every φ. For φ_u = |10⟩ + u|01⟩ the smallest eigenvalue is −11/9,
  at u = 1; for random two-copy vectors it lies between about −1.05 and −0.8. So this argument says nothing when the
  top eigenvalue of A_0 + A_1 is 1 − ε.
- Three copies carrying a full direction are different. Every single excitation on three copies with three nonzero
  weights has λ_max(T) > 1, so the anchor lemma itself does not apply there.

## Relative anchors and a one-vector conjecture

For a single-excitation plane both the second range vector and the co-range are restricted (to H_1 and to H_c), so
the anchor argument needs much less than T_p ⪯ I. For p ∈ H_1 and 1 ≤ c ≤ n let
- D_c be the compression of D_p to the column space H_c;
- R_c = M^{⊗n} V_p D_c^{−1} V_p^* M^{⊗n};
- Y_c(p) be the compression of (M^{⊗n} − R_c)^Γ to H_1 ⊗ H_c (row index in H_1, column index in H_c).

**Lemma 11.** If Y_c(p) ⪰ 0, then sector c is nonnegative for every single-excitation range plane that contains p.

*Proof.* Run the proof of the anchor lemma with v, w ∈ H_c and u ∈ H_1. Since H_1 and H_c have real bases, ū ∈ H_1 and
w̄ ∈ H_c, so only the compression Y_c(p) is used. ∎

Every vector of a single-excitation plane lies in H_1. Local unitaries that fix |0⟩ bring every p ∈ H_1 to
Σ_t w_t |1⟩_t with w_t ≥ 0. Lemma 11 therefore reduces the whole class to a statement about one vector.

**Conjecture 12.** Y_c(Σ_t w_t |1⟩_t) ⪰ 0 for every n, every 1 ≤ c ≤ n and all weights w_t ≥ 0.

Conjecture 12 implies that every sector of every single-excitation range plane is nonnegative for every n. It has
n − 1 real parameters, where the original problem has a plane and a co-range.
- **It holds when at most two weights are nonzero.** By charge conservation, R_c and R_p agree on H_1 ⊗ H_c: M^{⊗n}
  conserves q(row) − q(column) at every copy, so D_p is block diagonal in the column charge, and only the block H_c
  reaches these entries. Hence Y_c(p) is a compression of A − R_p^Γ, which is ⪰ 0 when T_p ⪯ I. For two nonzero
  weights, p = (w_0|10⟩ + w_1|01⟩) ⊗ |0…0⟩ up to the order of the copies, and Proposition 10 (a) gives T_p ⪯ I.
- Its kernel contains p̄ ⊗ H_c. For full support it contains n − 1 further vectors when c = n − 1.
- Numerics, three independent implementations: no violation found.
  - Direct computation of Y_c covers n ≤ 5 with every c, n = 6 with c ≤ 4, n = 7 with every c and n = 8 with
    c = 3, 4, 5. The weights were random, near-degenerate, structured or chosen by Nelder–Mead minimisation of
    λ_min(Y_c).
  - By Lemma 13 (below), the group conditions cover every n at once. They hold at about 950 weight vectors on four
    copies (including a full grid of step 1/48), 1000 on five, 900 on six and 7 on seven. The smallest values were
    re-evaluated in 50-digit arithmetic: they are exact zeros, and the next eigenvalue is positive.
  - At equal weights the smallest nonzero eigenvalue of the group conditions is 2^{−k}(4k − 5)/(5k − 4) for
    k = 3, …, 7. It tends to zero only at the boundary of the orthant, like (w_min/w_max)².

**Zero weights reduce to a finite test.** Let p = p_G ⊗ |0⟩^{⊗m}, where p_G is a single excitation on a group G of k
copies and the other m copies carry weight zero. On the k copies of G put A = (M^{⊗k})^Γ, B = R_{p_G}^Γ and Y = A − B.
Write G1_j for the states with row in H_1(G) and column in H_j(G), and G0_j for those with row = vacuum and column in
H_j(G).

**Lemma 13 (bath reduction).** Y_c(p) ⪰ 0 for every m ≥ 0 and every c if and only if
- (R∞) Y compressed to G1_j ⊕ G0_{j+1} is ⪰ 0 for −1 ≤ j ≤ k (G1_{−1} = G0_{k+1} = ∅, so j = −1 is the single state
  row = column = vacuum), and
- (R2) 3A + B compressed to G0_j is ⪰ 0 for 0 ≤ j ≤ k.

*Proof.*
1. **Factorisation.** By the charge argument above, Y_c(p) is the compression of A_G ⊗ A_0^{⊗m} − B_G ⊗ B_0^{⊗m} to
   H_1 ⊗ H_c. Here A_0 = I − ½F and B_0 = R_{|0⟩}^Γ on one bath copy. Both conserve the charge q(row) + q(column) of
   each copy and its letters. On the states that occur they are:
   - row and column vacuum: A_0 = B_0 = ½;
   - "c" = (row 0, column i) and "r" = (row i, column 0), in the basis (c, r): A_0 = [[1, −½], [−½, 1]] and
     B_0 = [[1, −½], [−½, 0]];
   - row and column both excited: A_0 = I − ½F and B_0 = ½F.
2. **Blocks.** A state of H_1 ⊗ H_c has one row excitation. So Y_c(p) splits into blocks labelled by the bath copies
   of each type. Hole copies contribute a common factor 2^{−N_0}.
3. **One bath copy excited in row and column.** The row excitation is there, and the other charged bath copies are
   "c". The block is A^{00} ⊗ (I − ½F) − B^{00} ⊗ ½F on G0_j. It is ½Y^{00} on symmetric letter pairs and
   ½(3A^{00} + B^{00}) on the antisymmetric one.
4. **No such copy.** Let N copies carry one column excitation each. The row excitation sits in G (state G1_j ⊗ c…c)
   or at one of these copies s (state G0_{j+1} ⊗ r_s). The tensor products give the block
   F_N = [[Y^{11}, −½Y^{10} ⊗ 1^T], [−½Y^{01} ⊗ 1, A^{00} ⊗ I_N + ¼Y^{00} ⊗ (J_N − I_N)]].
   - On 1^⊥ the lower block is ¼(3A^{00} + B^{00}).
   - On the symmetric vector, conjugation by diag(I, 2/√N) and then diag(I, −I) gives
     Y|_{G1_j⊕G0_{j+1}} + (1/N)(0 ⊕ (3A^{00} + B^{00})).
5. **Conclusion.** (R∞) and (R2) make every block ⪰ 0. Conversely, letting N → ∞ forces (R∞) for j ≥ 0, the blocks
   of step 3 force it for j = −1, and N = 2 forces (R2). ∎

So copies of weight zero never hurt once the group passes a finite test. That test does not depend on n.

**Proposition 14 (balanced triple).** Suppose A_0 + A_1 + A_2 has the eigenvalue 1, with a unit eigenvector g such
that ‖U_0 g‖ = ‖U_1 g‖ = ‖U_2 g‖. Then every sector is nonnegative, for every n.

*Proof.*
- As in Proposition 10 (b), U_t g = 0 for t ≥ 3. So the plane contains a vector that local unitaries fixing |0⟩ turn
  into a multiple of p_G ⊗ |0⟩^{⊗(n−3)}, with p_G = |100⟩ + |010⟩ + |001⟩.
- An exact computation over Q checks (R∞) and (R2) for this p_G, using symmetric Gaussian elimination
  (`code/check_balanced_triple_exact.py`). The state sets have sizes 1, 12, 48, 80, 48 and 1, 6, 12, 8.
- Lemma 13 and Lemma 11 give the result. ∎

For general three-copy weights, (R∞) and (R2) hold at every weight vector tested. That is a grid of the two-parameter
family plus random points, about 285 vectors. They are not proved there.

**Tangent planes.** Every plane L ⊆ span{|0…0⟩} ⊕ H_1 contains a nonzero p ∈ H_1, since the vacuum amplitude is a
linear functional on L.

*Lemma 11′.* Suppose that for every k the compression of (M^{⊗n} − R_p)^Γ to (H_0 ⊗ H_k) ⊕ (H_1 ⊗ H_{k−1}) is ⪰ 0.
Then q_n(X) ≥ 0 for every X whose range lies in such a plane containing p, with any co-range.

*Proof.* Use the proof of Lemma 11 with arbitrary v, w. Write u = u_0 + u_1 with u_0 ∈ H_0 and u_1 ∈ H_1. Since
(M^{⊗n} − R_p)^Γ conserves q(row) + q(column), the form (ū ⊗ w̄)^*(M^{⊗n} − R_p)^Γ(ū ⊗ w̄) splits into the pieces
ū_0 ⊗ w̄_k + ū_1 ⊗ w̄_{k−1}. ∎

For a vector with no zero weight, these compressions are exactly the conditions (R∞) of Lemma 13 with the whole
system as the group. Zero weights reduce in the same way.

*Lemma 13′.* For p = p_G ⊗ |0⟩^{⊗m}, the compressions of Lemma 11′ are ⪰ 0 for every m ≥ 0 and every k if and only
if (R∞) and (R2) hold for p_G.

*Proof.* Follow the proof of Lemma 13. Now a state may also carry no row excitation at all, G0_{j+1} ⊗ c…c. With N
charged bath copies, the block on G1_j ⊕ (G0_{j+1} ⊗ C^N) ⊕ G0_{j+1} is

    [[Y^{11}, −½Y^{10} ⊗ 1^T, Y^{10}], [−½Y^{01} ⊗ 1, A^{00} ⊗ I_N + ¼Y^{00} ⊗ (J_N − I_N), −½Y^{00} ⊗ 1],
     [Y^{01}, −½Y^{00} ⊗ 1^T, Y^{00}]].

- On 1^⊥ it is ¼(3A^{00} + B^{00}).
- On the symmetric vector, write the three components as z_1, z_2, z_3. The substitution z_3 ↦ z_3 − (√N/2) z_2 turns
  the block into Y|_{G1_j ⊕ G0_{j+1}} ⊕ ¼(3A^{00} + B^{00}).
- For N = 0 the block is Y|_{G1_j ⊕ G0_{j+1}}.
- The blocks with a bath copy excited in row and column are those of Lemma 13. ∎

**Corollary 15.** If Conjecture 12 holds for every n, then q_n ≥ 0 on every tangent plane (range in
span{|0…0⟩} ⊕ H_1), with any co-range, for every n.

*Proof.* By Lemma 13, Conjecture 12 for every n is equivalent to (R∞) and (R2) for every group vector with positive
weights. Up to local unitaries fixing |0⟩ and the order of the copies, the vector p ∈ L ∩ H_1 of a tangent plane L is
some p_G ⊗ |0⟩^{⊗m}. Apply Lemmas 13′ and 11′. ∎

Two cases hold without the conjecture: tangent planes containing a single excitation on at most two copies
(Proposition 10), and those containing a balanced three-copy single excitation (Proposition 14 and Lemma 13′).
Theorem 16 below extends this to every single excitation on three copies.
Numerically the compressions of Lemma 11′ are ⪰ 0 for n = 3 and 4, at every k, for random, structured and minimised
weights. For tangent planes this would replace the chain analysis of `tangent-class/`.

**Remark ((R2) always holds).** On G0_j both row strings are the vacuum. So B is there a principal submatrix of
R_{p_G} = M V D^{−1} V^* M ⪰ 0, and 3A + B ⪰ 3A ≻ 0. Lemmas 13 and 13′ therefore need only (R∞).

**Theorem 16 (three copies).** Conjecture 12 holds for every p ∈ H_1 supported on at most three copies: Y_c(p) ⪰ 0
for every n and every c.

*Proof.*
- **At most two copies.** This is Proposition 10, since Y_c(p) is a compression of A − R_p^Γ.
- **Three copies.** By Lemma 13 and the remark, it suffices to prove (R∞) for p_G = a|100⟩ + b|010⟩ + c|001⟩ with
  a, b, c > 0. The boundary of the orthant follows by continuity, because D_{p_G} stays invertible there.
  `code/check_group_exact.py "a,b,c"` proves (R∞) over the field Q(a, b, c):
  - It checks exactly that each (R∞) set splits into torus blocks.
  - In each block it runs symmetric Gaussian elimination with diagonal pivoting. Each pivot it uses is a ratio of two
    polynomials whose coefficients all have the same sign, so it is positive on the open orthant. When no nonzero
    diagonal entry is left, the rest of the block is identically zero.
  - The only denominators are determinants of blocks of D_{p_G}, which are positive on the orthant, and products of
    pivots. So the elimination commutes with evaluation at every point. At every point each block is congruent to a
    diagonal matrix with nonnegative entries. ∎

The run takes a few seconds with python-flint. An independent implementation in sympy
(`code/check_group_exact_sympy.py`) proves the same for the weights (1, v, u) in about 45 s, with the same ranks.
Since the matrices are homogeneous of degree 0 in the weights, that covers the orthant as well. The flint version
rejects the perturbed condition with B replaced by 1.001 B.

**Corollary 17.**
- (a) If A_0 + A_1 + A_2 has the eigenvalue 1, then every sector is nonnegative, for every n.
- (b) Every tangent plane whose single excitation is supported on at most three copies satisfies the endpoint, with
  any co-range, for every n.

*Proof.*
- (a) As in Proposition 10 (b), the plane contains a single excitation supported on copies 0, 1 and 2. Apply
  Theorem 16 and Lemma 11.
- (b) Apply Lemmas 13′ and 11′ with the certified (R∞). ∎

Corollary 17 (a) contains Propositions 9 (b), 10 (b) and 14.

**Larger groups.** The same certificate proves (R∞) in several further cases:
- at the equal-weight points for k = 4, 5 and 6 copies;
- on the lines (1, 1, 1, u) and (1, 1, 1, 1, u) for all u > 0;
- at the rational point (1, 3/7, 11/5, 1/1000, 5/4) for k = 5.

So the conclusions of Corollary 17 hold for single excitations of these shapes. The six-copy point takes about
eight minutes. For the whole four-copy family the certificate is not complete. In a run with squared weights, eleven
of the largest torus blocks (sizes 12 to 28, three parameters) did not finish.

## What remains for single-excitation ranges

The open part is the sectors 3 ≤ c ≤ n − 2 outside the cases above (W-planes, Propositions 5, 8, 9 and 10).
- **The shortest route is Conjecture 12.** It concerns one vector Σ_t w_t |1⟩_t and would settle every sector of every
  single-excitation plane, and with Corollary 15 every tangent plane, for every n. It is proved for at most three
  nonzero weights (Theorem 16). The first open case is four nonzero weights. By Lemma 13 that is a three-parameter
  family of finite matrix inequalities on four copies, independent of n. There the certificate of Theorem 16 has
  finished on all but eleven torus blocks. Results for every number of copies (closed forms, the condition
  (R∞)_{k−1} on most blocks, and equal weights up to 14 copies) are in `relative-anchors/`.

Three further precise statements would each make progress, and the first two would together settle the class for
every n:
- **One rank-two copy.** The inequality of Proposition 9 (a) for light directions that are not parallel, outside the
  regions of (b)–(d). The first case is n = 5, c = 3 (four light copies), for example s_max² = 0.9, s_min² = 0.2. Our random and adversarial searches up to N = 6 light copies found no
  violation.
- **Several rank-two copies.** An inductive argument for Proposition 9 (a′). Removing one rank-two copy leaves a
  quantitative condition on two adjacent sectors of the other copies. An all-n proof would need a statement about
  those sectors that survives adding a copy.
- **A rule for the splits.** An explicit choice of the matrices k_{S,u} in the balanced certificate of
  Proposition 7. Numerically, three numbers per pair suffice at n ≤ 6 with c = 3.

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
- `code/check_subsystem_bound.py` checks Proposition 5:
  - the pair decomposition of step 2 (error 3 × 10^{−16});
  - Theorem 3 on subsystems, step 4 (largest eigenvalue of 𝒟_Q^*𝒟_Q − 2Ā_Q P_Q at most 2 × 10^{−15});
  - the operator bound at random, heavy-copy, W-plane, near-extreme and adversarial POVMs, for n = 3, …, 6 and
    1 ≤ c ≤ n − 1 (smallest eigenvalue of the difference −2 × 10^{−16});
  - the corollary at random POVMs that satisfy the weight condition (n = 5, 6).
- `code/check_hardcore_form.py` checks the remarks after Proposition 5:
  - the hard-core form and the twisted identity (n = 3, 4, every c; error 4 × 10^{−16});
  - the c = 1 bound;
  - the failure of the vacuum-free bound in the heavy-copy example (n = 5, 6, every 2 ≤ c ≤ n − 2).
- `code/combined_certificate.py` sets up the certificate of Conjecture 6 as a semidefinite program (needs cvxpy with
  the CLARABEL solver). `code/check_combined_certificate.py` solves it at structured POVMs for n = 5, 6 and
  3 ≤ c ≤ n − 1, re-checks every inequality of each solution with independent eigenvalue computations, and shows
  that θ = 0 alone is infeasible at n = 6, c = 5 with one heavy copy. It also checks the quantitative form of
  Proposition 6: a certificate with common slack λ forces λ_min(Q_c) ≥ (c + 1)λ, since each column lies in c strings
  and one pattern. `code/sweep_combined_certificate.py` runs the
  60-POVM sweep at n = 5, c = 3. This is evidence, not proof. `code/check_certificate_counterexample.py` checks the
  counterexample at n = 7, 8: a symmetric reduction that matches the full program at n = 5, 6, the negative
  slacks, and the positive sector minimum.
- `code/threelevel_certificate.py` sets up the three-level certificate of Proposition 7 (cvxpy with CLARABEL or SCS),
  for the restricted heavy-copy model and for general range blocks. `code/check_threelevel_certificate.py` checks
  the identity T − α^*α = Vac − ⟨Θ, SΘ⟩ and its split over children (error 2 × 10^{−16}), the certificate on the
  heavy-copy family for n = 5, …, 10, general real POVMs at n = 5, c = 3 (free μ and balanced), the need for
  matrix-valued splits, and the child-form formula Π_{S'} = Σ_{s'} Tr(N N̄) with its rank-one antisymmetric part
  (complex POVMs, error 7 × 10^{−15}). Every solution is re-checked, with the guard λ_min(2K̃_c) ≥ (n + 1)λ. This is evidence,
  not proof.
- `code/check_structured_split.py` solves the three-number balanced certificates of Proposition 7 (n = 5, 6 and the
  heavy-copy family up to n = 10), rebuilds and re-checks every solution, and shows a closed-form split failing.
- `code/check_small_ranktwo.py` checks Proposition 8: the identity 2K̃_c = 4 Vac + α^*α − 2 Σ_t Δ_t^{exc} (complex
  POVMs, n = 3, 4, 5, every c; error 2 × 10^{−15}), the spectrum of Δ_t, the central-spin bound on every pattern, the
  corollary at 60 random POVMs, and an example outside the threshold of Proposition 5.
- `code/check_one_ranktwo.py` checks Proposition 9: the reduction (the sign of λ_min(2K̃_c), computed directly,
  matches the inequality in (a) in all 80 random cases), the commutator identity and m − 4K ⪰ (Σ|η_t|²) I of (b),
  λ_min(2K̃_c) ≥ 0 directly at 50 POVMs with an eigenvalue-1 copy (N ≤ 7, every c, rotated frames), 120 random
  configurations under the condition of (c), the failure of the separable split and of subadditivity, and the
  general congruence of (a′) at random POVMs with every copy of rank two (relative error 10^{−15}). Part (f) checks the chain
  of inequalities in the proof of (d) and λ_min(2K̃_c) ≥ 0 at 150 non-parallel configurations satisfying its condition.
- `code/check_pair_anchor_rank3_exact.py` is the exact proof behind Proposition 10 (a) for every two-copy vector
  (python-flint, a few seconds; the option `sanity` must fail). `code/check_pair_anchor_exact.py` is the earlier
  exact proof for Schmidt rank at most two (sympy, about 2 s). It checks D_{φ_u} ≻ 0,
  A − R^Γ ⪰ 0 and R^Γ + 3A ⪰ 0 for all u ≥ 0 on the 27 torus blocks, and the block structure itself.
  `code/check_pair_anchor.py` checks the same facts numerically with an independent implementation. It covers spec T_{|0⟩} = {1, −1/3}, the
  values u ∈ [0, 1] on a grid of 21 points and u = 1.5, 3, 10 (largest eigenvalue 1, smallest −11/9), random complex two-copy vectors of Schmidt rank two, and
  λ_min(2K̃_c) ≥ 0 directly at random configurations with A_0 + A_1 having eigenvalue 1 (n = 5, 6, every c).
- `code/relative_anchor.py` implements Y_c(p). It checks the Schur step of Lemma 11 as an identity (error 3 × 10^{−16}),
  that the full R_p agrees with R_c on H_1 ⊗ H_c, and Y_c ⪰ 0 with kernel containing p̄ ⊗ H_c at random weights
  (n ≤ 5, every c). `code/adversarial_relative.py` minimises λ_min(Y_c) over the weights. These two are evidence for
  Conjecture 12, not proof.
- `code/check_bath_reduction.py` checks Lemma 13. For random group vectors (k = 2, 3) and m = 1, 2, 3 bath copies, the
  eigenvalues of all blocks of the lemma reproduce the spectrum of Y_c(p) computed directly (every c; deviation
  ≤ 10^{−14}). It also checks Y_c ⪰ 0 directly for the balanced triple at n = 3, …, 6.
  `code/check_balanced_triple_exact.py` is the exact proof of (R∞) and (R2) for Proposition 14. It works over Q and
  rejects the false variant A − (11/10)B. `code/tangent_relative.py` gives the evidence for the tangent compressions
  (n = 3, 4). `code/check_tangent_bath.py` checks Lemma 13′: the eigenvalues of its blocks reproduce the spectrum of
  (M^{⊗n} − R_p)^Γ on (H_0 ⊕ H_1) ⊗ (all columns), computed directly, for k = 2, 3 and m = 1, 2, 3 (deviation
  ≤ 2 × 10^{−14}).
- `code/check_group_exact.py` (needs python-flint; `code/flintfield.py` is a small exact rational-function field on
  top of it) proves (R∞) and (R2) exactly for a group vector with rational or symbolic weights. `"a,b,c"` is the
  proof of Theorem 16 (a few seconds). The option `bscale=1001/1000` is the sanity test that must fail.
  `code/check_group_exact_sympy.py` is an independent implementation in sympy. Logs of the runs at
  (a, b, c), (1, 1, 1, 1), (1, 1, 1, u), (1, 1, 1, 1, 1), (1, 1, 1, 1, u), (1, 3/7, 11/5, 1/1000, 5/4) and
  (1, 1, 1, 1, 1, 1) are in `logs/`, with the sympy run at (1, v, u) and the failing sanity run.
- `code/check_wplanes.py` checks Theorem 2:
  - the central-spin value (m + 1)/2 for m ≤ 8;
  - the weighted bound (largest violation 2 × 10^{−15});
  - the hole form against direct partial traces (n = 3, 4, every h; error 10^{−15});
  - the exact minimum over all co-ranges for the W-plane and for weighted W-planes, n = 2, …, 5, with random rotations
    W_s (minimum 2^{1−n} for n ≥ 3).

The logs are in `logs/`.
