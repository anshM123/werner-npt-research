# (R∞)₁, the letter-1 block at two weight values, for every number of copies k

Ansh Mishra and Aryan Senthilkumar, 5 November 2026. Nothing here is refereed.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Theorem 2V.** For every k ≥ 4 and every weight vector whose squared weights take at most two values (u₁ copies of q₁ and u₂ copies of q₂, u₁ + u₂ = k, any u₁, u₂ ≥ 0, any q₁, q₂ > 0), the letter-1 block of (R∞)₁ (states P_{r,c} = (1@r, 1@c) and Q_{cc′} = (vac, 1@{c, c′}), size (3k² − k)/2) is ⪰ 0 with kernel exactly span{p ⊗ \|1@c⟩ : c ∈ [k]} | **proved**, computer-assisted, exact, symbolic in u₁, u₂ and in q₂/q₁ (Sections 2–4); 22 certificates, each a few seconds |
| 2 | The case u₁ = 1 (one copy with its own weight, all other weights equal), every k ≥ 4 | **proved** (contained in row 1; the regimes with one level of size 1) |
| 3 | One level (all weights equal), the letter-1 block of (R∞)₁ for every k ≥ 4: a second proof of this case of Theorem E1 of the README, by a different reduction | **proved**, computer-assisted (the regimes with u₂ = 0) |
| 4 | **Corollary.** For every k ≥ 4 and every weight vector with at most two distinct weights, all of (R∞)₁ is ⪰ 0 with kernel exactly p ⊗ H₁ | **proved**, given Theorems EQ, 2L, D₁ and Proposition T2′ of the README (themselves computer-assisted and not refereed) together with row 1 (Section 5) |
| 5 | Lemma 1: closed forms of all entries of the block in the rescaled frame, for all positive weights | **proved** (by hand, from Lemmas B and G and the Woodbury identity); **exact check** entry by entry against the Lemma B builder `j1lib.py` at 6 two-level weight vectors, with the B-part scaled by 1 and by 1 + 1/1000 (`check_entries.py`) |
| 6 | Lemmas 2 and 3: the isotypic decomposition under S_{u₁} × S_{u₂}, the test vectors, the pattern sums, the vacuum elimination | **proved**; **exact checks**: reduced forms against explicit projections of the `j1lib.py` block at 14 size pairs (k ≤ 7, both scalings; `check_reduce.py`) and at (5, 4), (6, 3) (k = 9; `check_reduce_large.py`); against `ylib.py` (built from the definitions, not from Lemma B) identically in q₂/q₁ at 7 size pairs (`check_ylib_sym.py`); the symbolic forms of the certificates against `j1lib.py` at 108 points, with the irrep on either level (`check_sym_points.py`); dimension counts and coverage of all (u₁, u₂) with 4 ≤ k ≤ 60 (`coverage.py`) |
| 7 | The conclusion on the whole block at 54 two-level rational weight vectors, k = 4, …, 7 | **exact check** (`check_block.py`; not a proof) |
| 8 | Negative controls: with B → (1 + 1/1000)B the certificates of the two parts that carry the kernel (inv, std) fail in every two-level regime, and so does the one of the (u − 2, 2) part at u₂ = 1, and at 7 explicit two-level rational points the scaled gauge-fixed block is not positive definite while the true one is | **exact** (`cert.py … 1001/1000`, `neg.py`) |
| 9 | Observation: the certificates of the remaining parts ((u − 2, 2) with the other level of size ≥ 2, (u − 2, 1, 1), and standard ⊗ standard) still pass with B → (1 + 1/1000)B | **exact** (certificate runs, `run_all.sh` steps `obs_scaled_*`); not used in any proof |
| 10 | Three or more weight values; general weights for k ≥ 5; Conjecture L1 (the Laplacian structure of the column-group elimination) | **open**; Theorem 2V is the block statement, not Conjecture L1 |

Nothing here says anything about Conjecture 12 of `single-excitation-range/` in general, or about the NPT (Werner
distillability) problem. What is new: the letter-1 block of (R∞)₁ for every k at two weight values (row 1), and with
the README results the whole condition (R∞)₁ for every k at two weight values (row 4). For k = 4 the block is already
covered for all weights by Theorem A4; at equal weights by Theorem E1.

Notation as in `relative-anchors/README.md` (Sections 1, 2) and `middle_j1/NOTES.md` (Section 4). Computer-assisted
means: exact rational-function arithmetic in python-flint and sign checks of coefficients; no floating point enters
any proof or check.

## 1. Setting and statement

The group has k copies. Put s_c = w_c² and L = |s| = Σ_c s_c. The letter-1 block of (R∞)₁ has the states
P_{r,c} = (1@r, 1@c) (r, c ∈ [k]) and Q_{ab} = (vac, 1@{a, b}) (a ≠ b), (3k² − k)/2 in all. We work in units 2^{k−1}
and in the rescaled frame of ylib: Ỹ = D Y D with D = diag(μ), μ(P_{r,c}) = w_r w_c, μ(Q_{ab}) = w_a w_b. Since D is
positive diagonal, Ỹ and Y have the same inertia and ker Ỹ = D⁻¹ ker Y. In this frame the kernel vector p ⊗ |1@c⟩ =
Σ_r w_r P_{r,c} becomes a multiple of the indicator 1_c of the column group {P_{r,c} : r ∈ [k]}; that these are in the
kernel is the identity Σ_r w_r 𝒯_{(r,1)} = E of Lemma N.

Two weight values: the copies form two levels U₁, U₂ with |U_i| = u_i and s_c = q_i on U_i, so L = u₁q₁ + u₂q₂. In the
scripts the levels are numbered 0 and 1, with sizes u0 = u₁, u1 = u₂ and squared weights Q0 = 1, Q1 = q = q₂/q₁. The
rescaled block is homogeneous in (q₁, q₂), so q₁ = 1 is no restriction.

**Theorem 2V.** Let k ≥ 4 and let the squared weights take at most two values (u₁ copies of q₁ > 0 and u₂ copies of
q₂ > 0, u₁ + u₂ = k, u₁, u₂ ≥ 0). Then the letter-1 block of (R∞)₁ is ⪰ 0, and its kernel is exactly
span{p ⊗ |1@c⟩ : c ∈ [k]}.

(For k = 3 the block has the extra kernel vector of Proposition K2; there the triples are all of [k], C₃ = 0 on them,
and the Woodbury form of Lemma G used below does not exist.)

## 2. Closed forms of the entries

Notation: d_c = 2(L − s_c), γ₀ = Σ_c s_c/d_c, σ = 2 + γ₀;
N₁[c, c] = Σ_{u≠c} s_u/(2(L − s_c − s_u));
Â = the k × k matrix with Â[c, c] = (1 + ½N₁[c, c])/s_c and Â[c, c′] = 1/(4(L − s_c − s_{c′})) for c ≠ c′, and
Ǎ = Â⁻¹;
for pairs P = {a, b}: N₂[P, P] = Σ_{u∉P} s_u/(2(L − s_a − s_b − s_u)), and Ŵ the matrix on pairs with
Ŵ[P, P] = (1 + ½N₂[P, P])/(s_a s_b), Ŵ[{c, a}, {c, b}] = 1/(4 s_c (L − s_a − s_b − s_c)) for a ≠ b, and 0 between
disjoint pairs.

**Lemma 1 (closed forms; k ≥ 4, all positive weights).**
1. *P–P.* Ỹ[(r₁, c₁), (r₂, c₂)] = (½ − [r₁ = r₂])[c₂ = r₁][c₁ = r₂] s_{r₁}s_{r₂} + [r₁ = r₂][c₁ = c₂](2 − [r₁ = c₁]) s_{r₁}s_{c₁}
   − β((r₂, c₁), (r₁, c₂)), where
   - β((r, c), (r′, c′)) = s_r s_{r′} γ₀/(2σ) if c = r and c′ = r′;
   - β = 2 s_r s_{r′} s_{c′}/(σ d_{c′}) if c = r and c′ ≠ r′ (and symmetrically);
   - β = 4 s_r s_{r′} ([c = c′] s_c/d_c − s_c s_{c′}/(d_c d_{c′} σ)) if c ≠ r and c′ ≠ r′.
2. *P–Q.* Ỹ[(r, c), {a, b}] = −s_r Ǎ[c, c*] if r ∈ {a, b}, with c* the other element of the pair, and
   Ỹ[(r, c), {a, b}] = s_r (s_a Ǎ[c, b] + s_b Ǎ[c, a])/(L − s_a − s_b) if r ∉ {a, b}.
3. *Q–Q.* Ỹ_QQ = 2Ŵ⁻¹, and Ŵ ≻ 0.
4. *Scaled B-part.* If every term of Lemma B that contains an E⁻¹ is multiplied by b (this is `bscale` of `j1lib.py`),
   the β-term in item 1 is multiplied by b; in item 2 the entry becomes A + b(entry − A) with
   A = −[r ∈ {a, b}][c = c*] s_a s_b; and Ỹ_QQ = 2D₂ + b(2Ŵ⁻¹ − 2D₂) with D₂ = diag(s_a s_b).

*Proof.* All three formulas of Lemma B are used with j = 1, in units 2^{k−1}.
- *E₁ and 𝒯.* On the letter-1 strings |1@c⟩ of H₁, C|1@c⟩ = 2(L − s_c)|1@c⟩ and J⁺J⁻|1@c⟩ = w_c w with
  w = Σ_t w_t|1@t⟩. So E₁ = Δ + ½wwᵀ with Δ = diag(d_c), and by Sherman–Morrison E₁⁻¹ = Δ⁻¹ − Δ⁻¹wwᵀΔ⁻¹/σ. In
  particular ⟨w, E₁⁻¹w⟩ = 2γ₀/σ and E₁⁻¹w = (2/σ)Δ⁻¹w. Next, 𝒯_{(1,r)}|1@c⟩ = ½[c = r] w + 2[c ≠ r] w_r |1@c⟩: the
  term ½J⁺σ⁻ needs the letter 1 at r, and the term 2w_r|0⟩⟨0|_r acts on a zero at r.
- *Item 1.* For the states (1@r₁, 1@c₁) and (1@r₂, 1@c₂), the swap term of Lemma B is
  ⟨1@c₁|σ⁺_{(1,r₂)}σ⁻_{(1,r₁)}|1@c₂⟩ = [c₂ = r₁][c₁ = r₂], the diagonal term is 2 − n_{r₁} = 2 − [r₁ = c₁] on |1@c₁⟩,
  and the B-term is ⟨𝒯_{(1,r₂)}|1@c₁⟩, E₁⁻¹ 𝒯_{(1,r₁)}|1@c₂⟩⟩. With the two facts above this inner product is, in the
  three cases, γ₀/(2σ), w_{r′} w_{c′} · 2/(σ d_{c′}) and 4 w_r w_{r′} (δ_{cc′}/d_c − w_c w_{c′}/(d_c d_{c′} σ)) (with
  (r, c) = (r₂, c₁), (r′, c′) = (r₁, c₂)). Multiplying by w_{r₁}w_{c₁}w_{r₂}w_{c₂} gives β.
- *Item 2.* Lemma B gives −⟨1@c|σ⁻_{(1,r)} − J⁻E₂⁻¹𝒯_{(1,r)}|P⟩. On the pair strings E₂ = C₂ + ½J⁺J⁻ with
  C₂|{a, b}⟩ = 2(L − s_a − s_b)|{a, b}⟩, and (I + ½J⁻C₂⁻¹J⁺)J⁻ = J⁻C₂⁻¹E₂ gives J⁻E₂⁻¹ = A⁻¹J⁻C₂⁻¹ with
  A = I + ½N₁, N₁ = J⁻C₂⁻¹J⁺ (this is Lemma L1.1 of `middle_j1`). If r ∈ P = {r, c*}, then σ⁻|P⟩ = |1@c*⟩ and
  𝒯_{(1,r)}|P⟩ = ½J⁺|1@c*⟩, so J⁻C₂⁻¹𝒯|P⟩ = ½N₁e_{c*} and A⁻¹·½N₁ = I − A⁻¹; the entry is
  −δ_{c c*} + δ_{c c*} − A⁻¹[c, c*]. If r ∉ P, then σ⁻|P⟩ = 0, 𝒯_{(1,r)}|P⟩ = 2w_r|P⟩ and
  J⁻C₂⁻¹|P⟩ = (w_a e_b + w_b e_a)/(2(L − s_a − s_b)). Multiply by w_r w_c w_a w_b. With D₁ = diag(w),
  w_c w_{c′} A⁻¹[c, c′] = (D₁A⁻¹D₁)[c, c′] = (D₁⁻¹AD₁⁻¹)⁻¹[c, c′], and D₁⁻¹AD₁⁻¹ = Â because
  N₁[c, c′] = w_c w_{c′}/(2(L − s_c − s_{c′})) for c ≠ c′.
- *Item 3.* By Lemma G the vacuum block is 2(I + ½N₂)⁻¹ with N₂ = J⁻C₃⁻¹J⁺ and C₃ = 2(L − s_a − s_b − s_c) > 0 on the
  triples (k ≥ 4). Its diagonal is N₂[P, P] above, the entry between {c, a} and {c, b} is w_a w_b/(2(L − s_a − s_b − s_c)),
  and disjoint pairs do not couple. Rescaling by D_Q = diag(w_a w_b) gives 2(D_Q⁻¹(I + ½N₂)D_Q⁻¹)⁻¹ = 2Ŵ⁻¹. Since N₂ ⪰ 0,
  Ŵ ≻ 0.
- *Item 4.* The terms without E⁻¹ (the A-part) are the swap and diagonal terms in item 1, −⟨1@c|σ⁻|P⟩ in item 2, and
  2I in item 3. ∎

## 3. The symmetry reduction

Let G = S_{u₁} × S_{u₂} act by permuting the copies inside each level. It preserves the weights, so it permutes the
states and commutes with D and Ỹ. Let V_i be the permutation module of U_i: V_i = triv ⊕ std_i if u_i ≥ 2 and V_i = triv
if u_i = 1. The P-states carry V ⊗ V with V = V₁ ⊕ V₂ (row index ⊗ column index); the Q-states carry the pairs inside
U₁, the pairs inside U₂ and V₁ ⊗ V₂. With std ⊗ std = triv ⊕ std ⊕ (u − 2, 2) ⊕ (u − 2, 1, 1) (for u ≥ 4; at u = 3 the
part (u − 2, 2) is absent and at u = 2 only triv is left) and pairs(U_i) = triv ⊕ std_i ⊕ (u_i − 2, 2), the irreps of G
that occur, and their multiplicities (P-states + Q-states), are:

| part | irrep of G | multiplicity (u_i ≥ 3, other level ≥ 2) | small sizes |
|---|---|---|---|
| inv | (triv, triv) | 6 + 3 | at u_i = 1 one P- and one Q-vector of that level drop out |
| std | (std₁, triv) | 5 + 2 | at u₁ = 2: 4 + 1; at u₂ = 1: no change; at u₂ = 0: 3 + 1 |
| two | ((u₁ − 2, 2), triv) | 1 + 1 | needs u₁ ≥ 4 |
| anti | ((u₁ − 2, 1, 1), triv) | 1 + 0 | needs u₁ ≥ 3 |
| mix | (std₁, std₂) | 2 + 1 | needs u₁, u₂ ≥ 2 |

and the same with the levels exchanged. The kernel span{1_c} carries V = V₁ ⊕ V₂: one trivial vector per level and one
copy of each std_i. So the expected kernel is 2-dimensional in the part inv (1-dimensional at one level), 1-dimensional
in std, and 0 in two, anti and mix. Lemma 2 checks this decomposition.

**Test vectors.** For each part we fix a seed vector of the irrep, with a few named copies, and list its images under
explicit equivariant maps. Rescaled coordinates X[r, c] on the P-states and φ{a, b} on the Q-states. With
f = e_a − e_b for named copies a, b of U₁:
- *inv* (no named copies): T1 = 1[r, c ∈ U₁, r ≠ c], T2 = 1[r = c ∈ U₁], T3 = 1[r ∈ U₁, c ∈ U₂], T4 = 1[r ∈ U₂, c ∈ U₁],
  T5 = 1[r, c ∈ U₂, r ≠ c], T6 = 1[r = c ∈ U₂]; T7 = 1[a, b ∈ U₁], T8 = 1[a, b in different levels],
  T9 = 1[a, b ∈ U₂]. Kernel vectors T1 + T2 + T4 and T3 + T5 + T6.
- *std*: S1 = f_r [r, c ∈ U₁, r ≠ c], S2 = f_c [r, c ∈ U₁, r ≠ c], S3 = f_r [r = c], S4 = f_r [c ∈ U₂],
  S5 = f_c [r ∈ U₂]; S6 = (f_a + f_b)[a, b ∈ U₁], S7 = f_a [b ∈ U₂] + f_b [a ∈ U₂]. Kernel vector S2 + S3 + S5. At u₁ = 2,
  S2 = −S1 and S6 = 0, and the kernel vector is −S1 + S3 + S5.
- *two*: h = e_{ac} + e_{bd} − e_{ad} − e_{bc} on the pairs of four named copies of U₁; H1 = h_{rc} (r ≠ c), H2 = h.
- *anti*: A = f₁ ∧ f₂ with f₁ = e_a − e_b, f₂ = e_b − e_c (three named copies of U₁); A1 = A[r, c].
- *mix*: f = e_a − e_b on U₁ and g = e_c − e_d on U₂; B1 = f_r g_c, B2 = g_r f_c, B3 = f_a g_b + f_b g_a.

**Lemma 2 (the reduction).** For each part, write v₁, …, v_m for the P-test vectors and b₁, …, b_n for the Q-test
vectors that exist at the given sizes, and R = (⟨x, Ỹ y⟩) over them.
1. The test vectors are images φ(seed) of the seed under linearly independent equivariant maps φ, and together they
   exhaust the state space: Σ_parts dim(irrep) · (number of test vectors) = k² for the P-states and k(k − 1)/2 for the
   Q-states. Hence the maps of each part form a basis of Hom_G(irrep, states).
2. Consequently Ỹ ⪰ 0 if and only if R ⪰ 0 for every part, and dim ker Ỹ = Σ_parts dim(irrep) · dim ker R.
3. Ỹ_QP maps the span of the v_i into the span of the b_l, so the Schur complement of the vacuum states, restricted to
   the part, is S = R_PP − R_PQ R_QQ⁻¹ R_QP; and R_QQ = 2 G_Q R_W⁻¹ G_Q with G_Q the Gram matrix of the b_l and R_W the
   reduced form of Ŵ, so S = R_PP − ½ R_PQ G_Q⁻¹ R_W G_Q⁻¹ R_QP. R_QQ ≻ 0, so R ⪰ 0 iff S ⪰ 0, with kernels in
   bijection (the Q-coordinates of a kernel vector are determined by its P-coordinates).

*Proof.*
1. Linear independence of the images is equivalent to that of the maps (an irrep is generated by any nonzero vector),
   and it holds since the Gram determinants of the test vectors are positive on every regime of Section 4
   (`cert.py` checks det G_P > 0 and det G_Q > 0 symbolically, by the same coefficient-sign test as the minors). The
   dimension identity: for u₁, u₂ ≥ 4, with dim std = u − 1, dim (u − 2, 2) = u(u − 3)/2, dim (u − 2, 1, 1) =
   (u − 1)(u − 2)/2, the P-count is 6 + 5(u₁ − 1) + 5(u₂ − 1) + (u₁² − 3u₁ + 1) + (u₂² − 3u₂ + 1) + 2(u₁ − 1)(u₂ − 1) =
   (u₁ + u₂)², and the Q-count is 3 + 2(u₁ − 1) + 2(u₂ − 1) + u₁(u₁ − 3)/2 + u₂(u₂ − 3)/2 + (u₁ − 1)(u₂ − 1) =
   k(k − 1)/2. `coverage.py` checks the identity, with the basis actually used at each size (small sizes included),
   for every (u₁, u₂) with 4 ≤ k ≤ 60, and `check_reduce.py` checks independence at its 14 size pairs. Independent
   equivariant images that fill the space must be bases of the multiplicity spaces.
2. The isotypic component of an irrep λ is λ ⊗ Hom_G(λ, states), and for a basis φ₁, …, φ_m of Hom_G(λ, states) the
   form (x, y) ↦ ⟨φ_i(x), Ỹ φ_j(y)⟩ is G-invariant on λ, hence equal to ⟨x, y⟩ Q_λ[i, j] (the irreps of symmetric
   groups and their products are absolutely irreducible). So R = ‖seed‖² Q_λ, and Ỹ on the component is congruent to
   ⟨·,·⟩ ⊗ Q_λ.
3. Ỹ_QP ∘ φ_i is an equivariant map from λ to the Q-states, hence a combination of the maps of the b_l (item 1). In the
   basis b, Ỹ_QQ restricted to their span has the matrix G_Q⁻¹R_QQ; this gives the formula for S. R_QQ is the reduced
   form of the positive definite Ỹ_QQ on independent vectors. The formula for R_QQ is Lemma 1, item 3, reduced in the
   same way (the reduced form of an invertible invariant operator O⁻¹ is G (reduced form of O)⁻¹ G). ∎

**Lemma 3 (computing R by pattern sums).** Every entry of R_PP, R_PQ, R_W and of the Gram matrices is a sum over pairs
of states of (coefficient of the first test vector) · (entry) · (coefficient of the second). Write each index of the two
states as a named copy of the seed or as an anonymous copy of level 1 or 2 (distinct from the named ones), and group the
terms by joint pattern: which indices coincide, and the level and name of each. The number of realisations of a
pattern is Π_i (u_i − n_i)(u_i − n_i − 1)⋯(u_i − n_i − m_i + 1), with n_i the number of named copies and m_i the number
of distinct anonymous copies of level i in the pattern. By Lemma 1 the entry depends only on the levels of the indices
and on which of them coincide; the sums over the other copies in γ₀, N₁ and N₂ are explicit in u₁ and u₂, and Ǎ is
obtained from the invariant 2 × 2 block of Â on (1_{U₁}, 1_{U₂}) and its two eigenvalues Â[c, c] − Â[c, c′] on the
standard parts. Hence every reduced form is a rational function of (q₁, q₂, u₁, u₂), equal to the true reduced form at
every admissible integer size. (`l1sym.py`; the Q-states are enumerated as ordered pairs with weight ½.)

All divisions in this computation are by quantities that are positive at admissible sizes (L − s_c − s_{c′} − s_u and
the like for k ≥ 4, the eigenvalues of the positive definite Â, the reduced forms of the positive definite Ŵ, orbit
sizes), so the rational functions specialise correctly at every admissible integer size.

## 4. Proof of Theorem 2V

1. By Lemmas 2 and 3 it suffices to show, for every part and every pair of sizes with k ≥ 4, that the reduced Schur
   complement S is ⪰ 0 with kernel exactly the span of the listed kernel vectors of the part. Then the kernel of the
   block has dimension 2 · 1 + (u₁ − 1) + (u₂ − 1) = k (with the obvious changes at u_i ≤ 1) and contains
   span{p ⊗ |1@c⟩}, which has dimension k; so it is that span.
2. *Regimes.* Each level size is either fixed (u = 0, 1 or 2) or symbolic, u = o + x with an offset o ∈ {3, 4} and
   x ≥ 0 real. Offsets of at least 3 are needed for the coefficient-sign test: a factor such as L − 3q₁ =
   (u₁ − 3)q₁ + u₂q₂ from the triples inside U₁ is one-signed only if u₁ − 3 ≥ 0. The certified regimes, with the
   irrep on the first level (the regimes of the parts on the second level follow by exchanging the levels, since every
   certificate holds for all q = q₂/q₁ > 0):

   | part | (u₁, u₂) regimes |
   |---|---|
   | inv | (≥ 3, 1), (2, 2), (≥ 3, 2), (≥ 3, ≥ 3), (≥ 4, 0) |
   | std | (≥ 3, 1), (2, 2), (2, ≥ 3), (≥ 3, 2), (≥ 3, ≥ 3), (≥ 4, 0) |
   | two | (≥ 4, 1), (≥ 4, 2), (≥ 4, ≥ 3), (≥ 4, 0) |
   | anti | (≥ 3, 1), (≥ 3, 2), (≥ 3, ≥ 3), (≥ 4, 0) |
   | mix | (2, 2), (2, ≥ 3), (≥ 3, ≥ 3) |

   `coverage.py` checks that every part at every (u₁, u₂) with u₁, u₂ ≥ 0, 4 ≤ k ≤ 60, lies in one of them; the
   regimes are uniform in the sizes, so the check is only bookkeeping.
3. *Certificate* (`cert.py PART U1 U2`, over Q(q, x₁, x₂)). The script builds S, checks that it is symmetric and that
   S κ = 0 identically for the kernel vectors κ of the part, and deletes the gauge coordinates T2, T6 (inv) and S3 (std),
   the test vectors supported on the diagonal states P_{c,c}. Each kernel vector has coefficient 1 on its own gauge
   coordinate and 0 on the other one, so S is ⪰ 0 with kernel exactly span(κ) if and only if the gauge-fixed S is
   positive definite. Every leading principal minor of the gauge-fixed S is a ratio of two polynomials in (q, x₁, x₂)
   with only positive coefficients, each with a monomial free of the x_i. So every minor is positive for q > 0 and
   x_i ≥ 0, and the gauge-fixed S is positive definite on the whole regime (Sylvester). No factorisation was needed:
   all numerators and denominators are one-signed as they stand. The largest minors:

   | regime | gauge-fixed size | last minor: numerator degree / terms | denominator degree / terms |
   |---|---|---|---|
   | inv (≥ 3, ≥ 3) | 4 | 65 / 7 100 | 56 / 4 493 |
   | std (≥ 3, ≥ 3) | 4 | 56 / 4 781 | 52 / 3 653 |
   | mix (≥ 3, ≥ 3) | 2 | 30 / 680 | 28 / 680 |
   | two (≥ 4, ≥ 3) | 1 | 18 / 220 | 18 / 220 |
   | anti (≥ 3, ≥ 3) | 1 | 6 / 20 | 6 / 20 |
   | inv (≥ 3, 1) | 3 | 20 / 155 | 15 / 135 |
   | std (≥ 3, 1) | 4 | 25 / 293 | 22 / 275 |
   | two (≥ 4, 1) | 1 | 9 / 54 | 9 / 54 |
   | anti (≥ 3, 1) | 1 | 3 / 10 | 3 / 10 |

   The other regimes are smaller (at fixed sizes the minors are polynomials in q alone); all logs are
   `logs/run_cert_*.log`. Each certificate takes at most about 5 seconds. ∎

So the proof uses Lemma 1 (from Lemmas B and G of the README), elementary representation theory of the symmetric
group, and exact arithmetic with rational functions. The case u₁ = 1 (one copy with its own weight, all other
weights equal) is the set of regimes in which one level has size 1: inv (≥ 3, 1), std (≥ 3, 1), two (≥ 4, 1),
anti (≥ 3, 1).

The one-level regimes (≥ 4, 0) give a second proof of the letter-1 block of (R∞)₁ at equal weights (Theorem E1 of the
README for j = 1), by a different reduction (the S_k-parts here are the two-level parts with an empty second level).

## 5. Consequence for (R∞)₁

(R∞)₁ is the compression of Y to G1₁ ⊕ G0₂. Its torus blocks are: the letter-1 block; the type I blocks with |T| = 1
and |T| = 2 (a state of G1₁ ⊕ G0₂ carries at most two letters in row and column together, so |T| ≤ 2); and the type II
blocks. For a weight vector with at most two distinct weights:
- the letter-1 block is ⪰ 0 with kernel exactly p ⊗ (letter-1 part of H₁) (Theorem 2V, k ≥ 4);
- each |T| = 1 block is ⪰ 0 with kernel exactly p ⊗ |2@t⟩ (Theorems EQ and 2L of the README: the weights on
  U = [k] ∖ {t} take at most two values, and the weight of t is arbitrary);
- each |T| = 2 block is positive definite (Theorem D₁, all weights);
- each type II block is positive definite (Proposition T2′, all weights).

**Corollary.** For every k ≥ 4 and every weight vector with at most two distinct weights, (R∞)₁ ⪰ 0 with kernel
exactly p ⊗ H₁. (k ≤ 3 is proved for all weights in `single-excitation-range/`.)

This rests on the README theorems as stated there; they are computer-assisted and not refereed, like Theorem 2V.

## 6. Validation, controls and observations

All checks are exact (rational arithmetic or identities of rational functions).
- *Builders.* `check_ylib.py` (copied unchanged from `middle_j1`): `j1lib.py` equals `ylib.Builder` on the letter-1
  and type I blocks of (R∞)₁ at general rational weights, k = 4, 5, 6.
- *Lemma 1, entry by entry* (`check_entries.py`): every P–P and P–Q entry and the whole vacuum block of the `j1lib.py`
  block equal the closed forms, and Ŵ is symmetric and positive definite, at 6 two-level weight vectors
  ((u₁, u₂) = (2, 2), (1, 3), (3, 1), (3, 2), (2, 4), (4, 2)), with the B-part scaled by 1 and by 1001/1000. Since the
  copies are concrete there, every equality pattern of up to four indices on both levels occurs.
- *Reduced forms against explicit projections* (`check_reduce.py`): at 14 size pairs ((2, 2), (3, 1), (1, 3), (3, 2),
  (2, 3), (4, 1), (1, 4), (4, 2), (3, 3), (5, 1), (4, 0), (5, 0), (4, 3), (5, 2)), both scalings, every part: R_PP, R_PQ,
  R_QQ and S from `l1sym.py` (exact rational arithmetic) equal the projections of the `j1lib.py` block on the explicit
  test vectors (S from the Schur complement of the whole vacuum block); the Gram matrices agree and are nonsingular;
  S κ = 0; the dimension counts hold. `check_reduce_large.py` repeats the comparison of R_PP, R_PQ, R_QQ, S and the
  Gram matrices at (5, 4) and (6, 3) (k = 9, block size 117), with the irrep on either level. (Every joint pattern of
  the pattern sums with nonzero count occurs at one of these sizes: the inv vectors have no named copies and at most
  four anonymous ones, which all occur at (4, ·) and (·, 4); the std vectors have at most one anonymous index per
  state; two, anti and mix have none.)
- *Against ylib over Q(q)* (`check_ylib_sym.py`): at (2, 2), (3, 1), (1, 3), (3, 2), (2, 3), (4, 1), (1, 4) the block
  built by `ylib.Builder` from the definitions (M, V_p, D, R and the partial transpose; not Lemma B) over Q(q), its
  kernel vectors 1_c, and the projections R_PP, R_PQ, R_QQ, S agree with `l1sym.py` identically in q.
- *The certified objects* (`check_sym_points.py`): for every regime, the symbolic S of `cert.py` evaluated at rational
  points (q = 1/9, 4, 25/49; sizes with x ∈ {0, 1}, k ≤ 7) equals the projection of the `j1lib.py` block, both with
  level 1 placed first and with level 2 placed first among the concrete copies; the kernel vectors and the leading
  minors agree and the minors are positive (108 comparisons). At larger sizes (x up to 9, k up to 21) the symbolic S
  equals the same reduction computed in rational arithmetic (consistency of the symbolic path).
- *Coverage* (`coverage.py`): 7 010 (size pair, part) cases with 4 ≤ k ≤ 60, all covered, all dimension counts right.
- *The conclusion directly* (`check_block.py`, exact checks at points, not a proof): at 54 two-level weight vectors
  (k = 4, 5, 6 with every split and weight ratios 1/40, 1/3, 7/2, 30; k = 7 with three splits), the indicators 1_c are
  in the kernel of the `j1lib.py` block, the block with the diagonal states P_{c,c} deleted is positive definite (exact
  LDL), and all reduced parts (irreps on either level) are positive definite, in agreement.
- *Negative controls* (B → (1 + 1/1000)B):
  - the certificate fails, with a minor having a factor of mixed signs, for inv in the regimes (≥ 3, 1), (2, 2),
    (≥ 3, 2), (≥ 3, ≥ 3), for std in (≥ 3, 1), (2, 2), (2, ≥ 3), (≥ 3, 2), (≥ 3, ≥ 3), and for two in (≥ 4, 1)
    (`run_all.sh`, steps `neg_cert_*`);
  - at 7 explicit points the gauge-fixed scaled block from `j1lib.py` is not positive definite while the true one is,
    and the reduced forms give the same verdicts and name the failing parts (`neg.py`): (u₁, u₂) = (2, 2), (3, 1),
    (1, 3), (1, 4), (4, 1), (3, 2), (2, 3) with weight ratios 1/10, 10 or 1/30; four of these points have a level of
    size 1, so the control also applies to the case u₁ = 1;
  - the failure depends on the gauge: with the states of row 0 deleted (the gauge of `middle_j1/check_l1.py`) the
    scaled block is still positive definite at w = (1, 1, 1/10, 1/10) and fails at (1, 1, 1/30, 1/30), the point
    reported there; with the diagonal states deleted (the gauge here) it fails at both (`neg.py`).
- *Observation* (not used): the certificates of two in (≥ 4, 2), (≥ 4, ≥ 3), of anti in all three two-level regimes
  and of mix in all three regimes still pass with the scaled B-part (`run_all.sh`, steps `obs_scaled_*`), so these
  parts stay positive definite when B is enlarged by 0.1%; the tight parts are the two that carry the kernel, and the
  part two when the other level is a single copy.

## 7. What is not covered

- Three or more weight values. The same reduction applies (more parts and larger forms: for three levels the
  gauge-fixed invariant part has size 9 in five variables); we have not tried it. For (R∞)₁ it would also need the
  |T| = 1 blocks at three U-weight values, which are open (`middle_j1_t1`, row 9).
- General weights for k ≥ 5 (no symmetry), and the letter-1 blocks of (R∞)_j for 2 ≤ j ≤ k − 2.
- Conjecture L1 (that each stage of the column-group elimination is a Laplacian) is a stronger structural statement;
  Theorem 2V does not prove it, even at two weight values.

## Files

| File | Content |
|---|---|
| `l1sym.py` | The closed forms (Lemma 1, class `Model`), the test vectors of every part, the pattern sums (Lemma 3), the reduced forms and the vacuum elimination (`reduce_part`), the bases at small sizes (`keep_for`) and the kernel vectors; generic over the number type (Fractions or flint rational functions) |
| `cert.py` | The certificates of Section 4 (`python cert.py PART U1 U2 [bscale]`, with U = an integer or `o+`) |
| `check_entries.py` | Lemma 1 entry by entry against `j1lib.py` |
| `check_reduce.py` | Reduced forms against explicit projections of the `j1lib.py` block; Gram matrices, kernels, dimension counts |
| `check_reduce_large.py` | The same comparison at k = 9 |
| `check_ylib_sym.py` | Reduced forms against `ylib.py` over Q(q), identically in q |
| `check_sym_points.py` | The symbolic forms used by the certificates against `j1lib.py` at rational points (both orders of the levels) and against the rational-arithmetic reduction at larger sizes |
| `coverage.py` | Coverage of all size pairs by the regimes; dimension counts for k ≤ 60 |
| `check_block.py` | The conclusion on the whole block at 54 weight vectors; agreement with the reduced forms |
| `neg.py` | Negative controls at explicit points; the gauge dependence of the control |
| `j1lib.py`, `ylib.py`, `flintfield.py`, `check_ylib.py`, `guard.sh` | Copied unchanged from `relative-anchors/code/middle_j1/` (`flintfield.py` also equals the copy in `middle_j1_t1/`) and `relative-anchors/code/no_letter1/guard.sh` (2.5 GB limit) |
| `run_all.sh` | Reruns everything, one Python process at a time under `guard.sh`; PASS/FAIL per step and ALL PASS at the end; about three minutes |
| `logs/run_all.out`, `logs/run_*.log` | Output of a complete rerun in a fresh copy of this folder (ALL PASS) |
| `logs/dev/` | Logs of the development runs (superseded by the final run) |
