# The letter-1 block of (R∞)₂ at two weight values: exact certificates for every k from 5 to 10

5 November 2026. Nothing here is refereed.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Theorem L.** For every k with 5 ≤ k ≤ 10 and every weight vector on k copies whose weights take at most two distinct values, the letter-1 block of (R∞)₂ (states P_{r,X} = (1@r, 1@X), X a pair, and Q_T = (vac, 1@T), T a triple; k·k(k − 1)/2 + k(k − 1)(k − 2)/6 states) is ⪰ 0 with kernel exactly span{p ⊗ \|1@X⟩ : X a pair} | **proved**, computer-assisted, exact: one certificate per size pair (u₀, u₁), u₀ ≥ u₁, u₀ + u₁ = k, each valid for every weight ratio q > 0 (Sections 2–4) |
| 2 | **Corollary.** For 5 ≤ k ≤ 10 and every weight vector with at most two distinct weights, (R∞)₂ ⪰ 0 with kernel exactly p ⊗ H₂ | **proved**, given row 1 and Theorems F1b, F2b (`code/middle_j2/`, valid for 5 ≤ k ≤ 26), DJ and Proposition T2′ of the README (all computer-assisted, not refereed) (Section 5) |
| 3 | Lemma A: closed forms of all entries of the block in the rescaled frame, all positive weights, k ≥ 5 | **proved** (by hand from Lemmas B and G of the README, Section 2); **exact check** entry by entry against the Lemma B builder `j1lib.py` (8 weight vectors, k = 5, 6, 7, general and two-value, B-part scaled by 1 and by 1001/1000) and against `ylib.py`, the builder from the definitions (3 weight vectors, k = 5, 6) |
| 4 | Lemma V: the vacuum block is 2Ŵ₃⁻¹ with an explicit Ŵ₃ (k ≥ 5); so it is positive definite and its Schur complement needs no inverse | **proved** (Lemma G with Woodbury); identity **checked exactly** at the 8 weight vectors of row 3 |
| 5 | Lemma Y: the reduction by Young symmetrizers at fixed level sizes, and the gauge | **proved** (standard representation theory of S_{u₀} × S_{u₁}, Section 3); **exact checks** in every certificate run: dimension counts Σ dim(λ)·mult = number of states, for the P-states, the kernel and the Q-states; **numerical check** (floating point): the spectrum of the Schur complement equals the union of the spectra of the parts with multiplicities, at the 8 two-level size pairs with k = 5, 6, 7 |
| 6 | The certified objects (the gauge-fixed forms of every part, as saved by the certificate runs) equal the projections of the Schur complement of the block built by `j1lib.py`, exactly, at rational points: every size pair with 5 ≤ k ≤ 10 (three weight ratios for k ≤ 6, one for k ≥ 7); and of the block built by `ylib.py` at 9 points (8 size pairs, k = 5, 6, 7) | **exact check** at points (`check_reduced.py`, `check_ylib3.py`) |
| 7 | End to end at points: the indicators of the column groups are in the kernel, the nullity is C(k, 2), and the block with one state per column group deleted is positive definite | **exact check** at 8 weight vectors, k = 5, 6, 7 (`check_builder.py`); not a proof |
| 8 | Negative controls: with B → (1 + 1/1000)B every two-level certificate with k = 5, 6, 7 fails (an explicit rational q where a leading minor of a reduced form is ≤ 0), and at an explicit two-level rational point for each of these size pairs the scaled gauge-fixed block built by `j1lib.py` is not positive definite while the true one is | **exact** (`run_all.sh` steps `neg_cert_*`, `neg_points.py`) |
| 9 | Observation: in every certificate of row 1, every leading principal minor of every gauge-fixed part is a ratio of two polynomials in q with only nonnegative coefficients (the Sturm fallback of the sign test was never used) | **exact** observation (logs); an analogue, for these reduced forms, of the positivity pattern of Section 3g of the README (which concerns principal minors in the rescaled state basis); not used beyond the sign test |
| 10 | The same for every k (level sizes symbolic, as Theorem 2V for (R∞)₁); k > 10; three or more weight values; general weights | **open**; not attempted symbolically. The fixed-size method extends to larger k at a cost that grows by a factor of about 3 to 4 per step in k (Section 7) |
| 11 | A second check written separately (Section 8), using only `j1lib.py`, `ylib.py` and `r2lib.py`: at a fresh two-value point for each of the 27 size pairs (levels at random positions, weight ratios not used above, some with an overall scale) the block built by `j1lib.py` has the C(k, 2) vectors p ⊗ \|1@X⟩ in its kernel and is positive definite without one state per column group (exact certificate), so it is ⪰ 0 with exactly that kernel there; the same for the block built by `ylib.py` at 5 points with non-square ratios q (k = 5, 6); `j1lib.py` equals `ylib.py` entry by entry at 2 points; the saved certified forms are positive definite at 31 values of q from 10⁻⁸ to 10⁸ | **exact check** at points (`crosscheck.py`, `run_crosscheck.sh`); not a proof |

Nothing here says anything about Conjecture 12 of `single-excitation-range/` in general, or about the NPT (Werner
distillability) problem. Computer-assisted means: exact rational-function arithmetic in one variable (python-flint
fmpq_poly) and sign checks of coefficients; floating point enters only the spectral check of row 5, which is not used in
any proof.

## 1. Setting and statement

Notation as in `relative-anchors/README.md` (Sections 1, 2) of the repository `werner-npt-research`. The group has k
copies, s_c = w_c², L = |s|. For a set Z of copies write π_Z = Π_{z∈Z} s_z, σ_Z = Σ_{z∈Z} s_z, w_Z = Π_{z∈Z} w_z, and
|Z⟩ for the column string with the letter 1 exactly on Z. (R∞)₂ is Y compressed to G1₂ ⊕ G0₃. Its letter-1 block has the
states P_{r,X} = (1@r, 1@X) (r ∈ [k], X a pair) and Q_T = (vac, 1@T) (T a triple). We use the entries of Lemma B in units
2^{k−2} and the rescaled frame Ỹ = DYD of `ylib.py`: D = diag(μ), μ(P_{r,X}) = w_r w_X, μ(Q_T) = w_T. Positive diagonal
rescaling keeps inertia, and the kernel vector p ⊗ |1@X⟩ = Σ_r w_r P_{r,X} becomes (1/w_X) times the indicator 1_X of the
column group {P_{r,X} : r ∈ [k]}.

Two weight values: u₀ copies with squared weight 1 and u₁ copies with squared weight q > 0 (in the scripts, level 0 is
copies 0, …, u₀ − 1).

**Theorem L.** Let 5 ≤ k ≤ 10 and let the weights take at most two distinct values. Then the letter-1 block of (R∞)₂ is
⪰ 0 and its kernel is exactly span{p ⊗ |1@X⟩ : X a pair of copies}.

The cases reduce to the certified ones as follows (`coverage.py`). A permutation of the copies permutes the states. Ỹ
is homogeneous of degree 3 in the squared weights (Lemma A; checked exactly in `check_builder.py`), so dividing all
squared weights by a constant multiplies Ỹ by a positive number. Hence a weight vector with m copies of squared weight α
and k − m of β gives, up to these operations, the size pair (u₀, u₁) = (max(m, k − m), min(m, k − m)) at
q = (squared weight of the smaller level)/(squared weight of the larger level) > 0, and equal weights are the pair
(k, 0).

## 2. Closed forms (Lemma A) and the vacuum block (Lemma V)

Let D_w = diag(w_Z) on the letter-1 strings, j̄⁺|Z⟩ = Σ_{t∉Z} |Z + t⟩ and j̄⁻ = (j̄⁺)ᵀ (unweighted), and for a copy r let
t_r be the operator t_r|X⟩ = ½ Σ_{t ∉ X − r} |X − r + t⟩ if r ∈ X and t_r|X⟩ = 2|X⟩ if r ∉ X. On the strings with m
letters 1 put

    Ê_m = D_w⁻¹ E_m D_w⁻¹:  Ê_m[Z, Z] = 2(L − σ_Z)/π_Z + ½ Σ_{t∈Z} 1/π_{Z−t},   Ê_m[Z, Z′] = 1/(2π_{Z∩Z′}) if |Z ∩ Z′| = m − 1,

and 0 otherwise.

**Lemma A (k ≥ 5, all positive weights).** In units 2^{k−2} and in the rescaled frame:
1. Ỹ[P_{r,X′}, P_{r′,X}] = (½ − [r = r′])[r ∈ X][X′ = X − r + r′] s_{r′}π_X + [r = r′][X = X′](2 − [r ∈ X]) s_r π_X
   − s_r s_{r′} ⟨t_{r′}|X′⟩, Ê₂⁻¹ t_r|X⟩⟩ (the condition X′ = X − r + r′ includes r′ ∉ X − r);
2. Ỹ[P_{r,X′}, Q_T] = −[r ∈ T][X′ = T − r] π_T + s_r ⟨X′| j̄⁻ Ê₃⁻¹ t_r |T⟩;
3. Ỹ[Q_{T′}, Q_T] = 2[T = T′] π_T − ⟨T′| j̄⁻ Ê₄⁻¹ j̄⁺ |T⟩.
With every term of Lemma B that contains an E⁻¹ multiplied by b (`bscale`), the three E⁻¹-terms are multiplied by b.

*Proof.* On the letter-1 strings, J⁻D_w⁻¹ = D_w⁻¹j̄⁻ and D_w⁻¹J⁺ = j̄⁺D_w⁻¹ (⟨Y|J⁻|Z⟩ = w_t when Z = Y + t, and
w_t/w_Z = 1/w_Y), and 𝒯_{(1,r)} = ½J⁺σ⁻_{(1,r)} + 2w_r|0⟩⟨0|_r = w_r D_w t_r D_w⁻¹ (the first term moves the letter at r
to any empty copy t, with weight w_t, and w_X w_t/w_{X−r+t} = w_r). Hence E_m = D_w Ê_m D_w with Ê_m = D_w⁻¹CD_w⁻¹ +
½ j̄⁺D_w⁻²j̄⁻, which gives the entries above, and:
- 𝒯_{(1,r′)}ᵀ E₂⁻¹ 𝒯_{(1,r)} = w_r w_{r′} D_w⁻¹ t_{r′}ᵀ Ê₂⁻¹ t_r D_w⁻¹; multiplied by μ(P_{r,X′})μ(P_{r′,X}) =
  w_r w_{r′} w_{X′} w_X this is item 1's B-term (Lemma B, first formula, with the row letter a = (1, r) acting on the
  column string X and the column letter b = (1, r′) on the row string X′: the partial transpose). The swap term
  σ_b⁺σ_a⁻|X⟩ = |X − r + r′⟩ needs r ∈ X and r′ ∉ X − r, and w_r w_{X′} w_{r′} w_X = s_{r′}π_X for X′ = X − r + r′; the
  diagonal term is 2 − n_r = 2 − [r ∈ X].
- J⁻E₃⁻¹𝒯_{(1,r)} = w_r D_w⁻¹ j̄⁻ Ê₃⁻¹ t_r D_w⁻¹, and −⟨X′|σ⁻_{(1,r)}|T⟩ = −[r ∈ T][X′ = T − r]; multiplied by
  w_r w_{X′} w_T this is item 2 (Lemma B, second formula).
- J⁻E₄⁻¹J⁺ = D_w⁻¹ j̄⁻ Ê₄⁻¹ j̄⁺ D_w⁻¹; multiplied by w_{T′}w_T this is item 3 (Lemma B, third formula). ∎

All entries are rational in the squared weights and homogeneous of degree 3 (Ê_m has degree 1 − m, and t_r, j̄^± are
constant).

**Lemma V (k ≥ 5).** Ỹ_QQ = 2Ŵ₃⁻¹ with

    Ŵ₃[T, T] = (1 + ½ Σ_{t∉T} s_t/(2(L − σ_T − s_t)))/π_T,   Ŵ₃[T, T′] = 1/(4π_{T∩T′}(L − σ_{T∪T′})) if |T ∩ T′| = 2,

and 0 otherwise. So Ỹ_QQ ≻ 0, and the Schur complement onto the P-states is S = Ỹ_PP − ½ Ỹ_PQ Ŵ₃ Ỹ_QP.

*Proof.* Lemma G of the README: for k ≥ 5 the energies C₄ = 2(L − σ_U) of the strings with four letters 1 are positive,
and I − ½J⁻E₄⁻¹J⁺ = (I + ½N₃)⁻¹ with N₃ = J⁻C₄⁻¹J⁺ (Woodbury). N₃[T, T] = Σ_{t∉T} s_t/(2(L − σ_T − s_t)) and, for
|T ∩ T′| = 2, N₃[T, T′] = w_t w_{t′}/(2(L − σ_{T∪T′})) with {t} = T′ ∖ T, {t′} = T ∖ T′. Rescaling,
Ỹ_QQ = 2D_T(I + ½N₃)⁻¹D_T = 2(D_T⁻¹(I + ½N₃)D_T⁻¹)⁻¹ with D_T = diag(w_T), and D_T⁻¹(I + ½N₃)D_T⁻¹ = Ŵ₃. N₃ ⪰ 0, so
Ŵ₃ ≻ 0. ∎

Since the vectors 1_X (on the P-states, zero on the Q-states) lie in the kernel of Ỹ (Lemma N of the README; checked as an
identity in q in every certificate), Ỹ_QP 1_X = 0 and S 1_X = 0, and Ỹ ⪰ 0 with kernel exactly span{1_X} if and only if
S ⪰ 0 with kernel exactly span{1_X} (a kernel vector of Ỹ is (x, −Ỹ_QQ⁻¹Ỹ_QP x) with x ∈ ker S).

`l2block.py` implements Lemma A for any number field (rationals, or rational functions of q). `check_builder.py`
compares it entry by entry with the rescaled `j1lib.py` block at 8 weight vectors (k = 5, 6, 7; general weights, and
two values with the levels interleaved), with b = 1 and b = 1001/1000, and checks Lemma V and the homogeneity there.
`check_ylib3.py` compares it with `ylib.py` (built from M, V_p, D_p, R_p and the partial transpose, not from Lemma B) at
3 weight vectors.

## 3. The reduction (Lemma Y) and the gauge

G = S_{u₀} × S_{u₁} permutes the copies inside each level. It permutes the states and commutes with Ỹ and with S. The
P-states span a permutation module M_P, the Q-states M_Q, and K = span{1_X} ⊂ M_P is G-stable.

For an irrep λ = (λ₀, λ₁) of G fill the Young diagrams of λ₀ and λ₁ row by row with the copies of level 0 and level 1,
and let e_λ = (Σ_{ρ∈R} ρ)(Σ_{κ∈C} sgn(κ)κ) be the Young symmetrizer (R, C the row and column groups). Since
e_λ x e_λ ∈ C e_λ for every x in the group algebra (von Neumann's lemma), e_λ acts as zero on every irrep other than λ
and with rank one on λ. So on any permutation module M, e_λ M has dimension mult(λ, M), and its vectors are s₀ ⊗ c for a
fixed vector s₀ of λ. If Ỹ (or S) is G-invariant and symmetric, then on the λ-isotypic component it is ⟨·,·⟩_λ ⊗ Q_λ
(irreps of symmetric groups are absolutely irreducible and self-dual), and for any basis B of e_λ M, BᵀSB =
|s₀|² CᵀQ_λC with C invertible.

**Lemma Y.** Let W ⊂ M_P be spanned by the states P_{r,X} with r ∉ X and by P_{a,X} − P_{b,X} for X = {a, b}. Then W is a
G-stable complement of K in M_P, and e_λM_P = e_λK ⊕ e_λW. Suppose S1_X = 0 for every pair X. Then S ⪰ 0 with kernel
exactly K if and only if, for every λ that occurs, the form S_λ = BᵀSB on a basis B of e_λW is positive definite.

*Proof.* K ∩ W = 0: a vector Σ_X c_X 1_X of K has the coordinates (c_X, c_X) on the two states P_{a,X}, P_{b,X}
(X = {a, b}), while every vector of W has opposite coordinates there; so all c_X = 0. dim W = (k − 2)C(k, 2) + C(k, 2)
= k·C(k, 2) − C(k, 2), so M_P = K ⊕ W. Since e_λ is linear, e_λM_P = e_λK + e_λW, and
the dimensions add (mult(λ, M_P) = mult(λ, K) + mult(λ, W)). If S1_X = 0, then in a basis of e_λM_P consisting of a basis
of e_λK followed by B the reduced form is diag(0, S_λ), so Q_λ ⪰ 0 with kernel of dimension mult(λ, K) if and only if
S_λ ≻ 0. Summing over λ, ker S has dimension Σ dim(λ) mult(λ, K) = dim K and contains K. ∎

*Vacuum elimination in the reduced form.* Ỹ_QP maps e_λM_P into e_λM_Q. With B_Q a basis of e_λM_Q, G_Q = B_QᵀB_Q and
R_W = B_QᵀŴ₃B_Q, the reduced form of S on B is BᵀỸ_PPB − ½R_QWᵀG_Q⁻¹R_WG_Q⁻¹R_QW with R_QW = B_QᵀỸ_QPB (as in
Lemma 2 of `middle_j1_l1/NOTES.md`). All matrices here are small; no inverse of a block of states is needed.

`young.py` builds e_λ (row symmetrization by orbit sums of the row group, column antisymmetrization by enumeration) and
reduced row echelon bases of e_λK, e_λW and e_λM_Q in exact rational arithmetic. Only irreps with at most three boxes
outside the first rows can occur (states have at most three indices); every certificate checks
Σ_λ dim(λ)·mult(λ, M_P) = k·C(k, 2), Σ_λ dim(λ)·mult(λ, K) = C(k, 2) and Σ_λ dim(λ)·mult(λ, M_Q) = C(k, 3), so no
irrep is missed.

## 4. Certificates and the proof of Theorem L

`cert_fixed.py U0 U1` works over Q(q) with exact univariate rational functions (`ufield.py`, python-flint fmpq_poly):
1. builds Ỹ_PP and Ỹ_PQ (Lemma A; Ê₂⁻¹ and Ê₃⁻¹ by exact Gauss–Jordan over Q(q)) and Ŵ₃ (Lemma V);
2. checks Ỹ_PP1_X = 0 and Ỹ_QP1_X = 0 for every pair X, identically in q;
3. for every irrep λ: the bases of Section 3, the dimension counts, the reduced kernel identities, and S_λ;
4. computes every leading principal minor of every S_λ (by symmetric elimination) and checks that it is > 0 for every
   q > 0: its numerator and denominator are polynomials in q with all coefficients ≥ 0 (and not zero). A Sturm count of
   the roots in (0, ∞) is the fallback if a polynomial has coefficients of both signs; it was never needed.

The entries of S_λ, as rational functions of q in lowest terms, are defined at every q > 0 and equal the true reduced
form there: Ê₂ and Ê₃ are positive definite for every q > 0 (they are congruent to E₂, E₃ ≻ 0), so their inverses over
Q(q) specialise correctly, and the other denominators (L − σ_Z, the entries of Ŵ₃, the constant Gram matrix G_Q) are
positive or constant. A leading minor is a polynomial in the entries, so its value at q equals the rational function
computed over Q(q), evaluated at q.

By Sylvester's criterion each S_λ is positive definite for every q > 0, so by Lemmas V and Y the block is ⪰ 0 with
kernel exactly span{p ⊗ |1@X⟩} for every q > 0 and the given sizes. `run_all.sh` runs the certificate for every pair
(u₀, u₁) with u₀ ≥ u₁, u₀ + u₁ = k, 5 ≤ k ≤ 10, and `coverage.py` checks that every two-value weight vector with
5 ≤ k ≤ 10 reduces to one of them (Section 1). ∎

Sizes, from the final run (`logs/run_cert_k*_*.log`): the number of parts with a nonempty gauge-fixed form, the largest
gauge-fixed size, the number of leading minors checked, the largest degrees in q of a numerator and of a denominator,
and the time. In all 27 certificates all 808 minors have numerator and denominator with only nonnegative coefficients
(at one level, u₁ = 0, they are positive rational numbers).

| k | (u₀, u₁) | parts | largest part | minors | num. degree | den. degree | time |
|---|---|---|---|---|---|---|---|
| 5 | (5, 0) | 5 | 3 | 9 | 0 | 0 | 0.1 s |
| 5 | (4, 1) | 4 | 7 | 17 | 34 | 32 | 0.1 s |
| 5 | (3, 2) | 6 | 8 | 27 | 43 | 37 | 0.2 s |
| 6 | (6, 0) | 6 | 3 | 10 | 0 | 0 | 0.2 s |
| 6 | (5, 1) | 5 | 7 | 19 | 35 | 33 | 0.4 s |
| 6 | (4, 2) | 8 | 9 | 34 | 51 | 45 | 0.6 s |
| 6 | (3, 3) | 8 | 8 | 36 | 54 | 40 | 0.6 s |
| 7 | (7, 0) | 6 | 3 | 10 | 0 | 0 | 0.6 s |
| 7 | (6, 1) | 6 | 7 | 20 | 38 | 36 | 1.4 s |
| 7 | (5, 2) | 9 | 9 | 36 | 60 | 54 | 2.0 s |
| 7 | (4, 3) | 10 | 9 | 43 | 65 | 58 | 2.4 s |
| 8 | (8, 0) | 6 | 3 | 10 | 0 | 0 | 1.5 s |
| 8 | (7, 1) | 6 | 7 | 20 | 34 | 32 | 5.8 s |
| 8 | (6, 2) | 10 | 9 | 37 | 54 | 48 | 7.9 s |
| 8 | (5, 3) | 11 | 9 | 45 | 60 | 53 | 8.7 s |
| 8 | (4, 4) | 12 | 9 | 50 | 69 | 51 | 8.5 s |
| 9 | (9, 0) | 6 | 3 | 10 | 0 | 0 | 3.8 s |
| 9 | (8, 1) | 6 | 7 | 20 | 40 | 38 | 22 s |
| 9 | (7, 2) | 10 | 9 | 37 | 62 | 56 | 31 s |
| 9 | (6, 3) | 12 | 9 | 46 | 70 | 63 | 34 s |
| 9 | (5, 4) | 13 | 9 | 52 | 79 | 64 | 34 s |
| 10 | (10, 0) | 6 | 3 | 10 | 0 | 0 | 9 s |
| 10 | (9, 1) | 6 | 7 | 20 | 38 | 36 | 77 s |
| 10 | (8, 2) | 10 | 9 | 37 | 59 | 53 | 99 s |
| 10 | (7, 3) | 12 | 9 | 46 | 67 | 60 | 113 s |
| 10 | (6, 4) | 14 | 9 | 53 | 76 | 61 | 116 s |
| 10 | (5, 5) | 14 | 9 | 54 | 78 | 60 | 122 s |

Most of the time goes into the exact inverse of Ê₃ over Q(q) (size C(k, 3)); the reductions take a few seconds.

## 5. Consequence for (R∞)₂

For a weight vector with at most two distinct weights and 5 ≤ k ≤ 10, every torus block of (R∞)₂ is now
covered: the letter-1 block by Theorem L; the type I blocks with |T| = 2 and |T| = 1 by Theorems F1b and F2b of
`code/middle_j2/NOTES.md` (they hold at two weight values for every k ≤ 26); the type I blocks with |T| = 3 by Theorem DJ;
the type II blocks by Proposition T2′. The kernels are the parts of p ⊗ H₂ that each block contains (p ⊗ |1@X⟩,
p ⊗ |2@t 1@c⟩, p ⊗ |2@t 2@τ⟩; the other blocks are positive definite).

**Corollary.** For 5 ≤ k ≤ 10 and every weight vector whose weights take at most two distinct values, (R∞)₂ ⪰ 0 with
kernel exactly p ⊗ H₂. At k = 4 all of (R∞)₂ is Theorem B4 (every weight vector), whose kernel is
p ⊗ H₂ ⊕ span(vac ⊗ J⁻|1111⟩).

This rests on the theorems of the README and of `code/middle_j2/` as stated there; they are computer-assisted and not
refereed, like Theorem L.

## 6. Checks and controls

All exact unless said otherwise.
- *Builders* (`check_builder.py`): Lemma A against `j1lib.py` entry by entry at 8 weight vectors (k = 5: three, k = 6:
  three, k = 7: two; general weights and two values with the levels interleaved), both scalings; Lemma V; homogeneity
  of degree 3; and, end to end, the indicators 1_X in the kernel, nullity C(k, 2), and the block with the states
  P_{min X, X} deleted positive definite (exact LDL). `check_ylib3.py` (a): 2^{k−2} times the `ylib.py` block equals
  Lemma A at 3 weight vectors (k = 5, 6).
- *The certified objects* (`check_reduced.py`): every certificate saves the forms S_λ it certifies (coefficient lists,
  `logs/cert_k*_*.pkl`). For every size pair with 5 ≤ k ≤ 10 they are evaluated at rational q (q = 1/9, 4, 49/25 for
  k ≤ 6; q = 1/9 for k ≥ 7; q = 1 at one level) and compared with the projections BᵀS_refB of the Schur complement
  S_ref = Ỹ_PP − Ỹ_PQỸ_QQ⁻¹Ỹ_QP of the rescaled `j1lib.py` block (exact inverse of its vacuum block, no Woodbury, flint
  fmpq_mat): equal for every part; the kernel vectors e_λK are annihilated by S_ref; all pivots positive.
  `check_ylib3.py` (b) does the same with the `ylib.py` block at 9 points: every size pair with k = 5, 6 (two weight
  ratios at (3, 2)) and (4, 3).
- *Completeness of the parts* (floating point, not used in a proof): at q = 1/4 and every two-level pair with k = 5, 6, 7
  the spectrum of S_ref (k·C(k, 2) eigenvalues, C(k, 2) of them zero) equals the union over λ of the eigenvalues of
  G⁻¹R on the basis e_λK + e_λW, each with multiplicity dim λ, to 10⁻⁹ relative.
- *Negative controls* (B → (1 + 1/1000)B, i.e. `bscale` 1001/1000 in `j1lib.py` and in Lemma A):
  - `cert_fixed.py U0 U1 1001/1000` fails for every two-level pair with k = 5, 6, 7, and names a rational q at which a
    leading minor of a reduced form is ≤ 0 (`run_all.sh` steps `neg_cert_*`, expected FAIL);
  - `neg_points.py`: for each of these pairs, at an explicit two-level rational point, the scaled block built by
    `j1lib.py` with the states P_{min X, X} deleted is not positive definite while the true one is (exact LDL), and the
    reduced forms at that point name the failing parts.

## 7. What is not covered

- *Every k.* The certificates are at fixed sizes. A proof for every k needs, as in Theorem 2V, explicit test vectors
  for each part with the level sizes symbolic, pattern sums over joint equality patterns of up to six indices, and
  closed forms of Ê₂⁻¹ and Ê₃⁻¹ (or of their Woodbury forms on singles and pairs) as functions of the patterns and of
  the level sizes. The parts that occur are those of G with at most three boxes outside the first rows (sixteen
  families with two levels; with one level (u), (u − 1, 1), (u − 2, 2), (u − 2, 1, 1), (u − 3, 3), (u − 3, 2, 1), of
  multiplicities 2, 4, 3, 2, 1, 1 in M_P). This was not attempted here. For orientation, the parts at (u₀, u₁) = (5, 4)
  (`logs/run_cert_k9_5_4.log`), with their multiplicity in M_P, in K, the gauge-fixed size and the number of Q-vectors:

  | λ | dim | M_P | K | gauge-fixed | Q |
  |---|---|---|---|---|---|
  | ((5), (4)) | 1 | 10 | 3 | 7 | 4 |
  | ((5), (3, 1)) | 3 | 11 | 2 | 9 | 3 |
  | ((5), (2, 2)) | 2 | 4 | 1 | 3 | 1 |
  | ((5), (2, 1, 1)) | 3 | 3 | 0 | 3 | 0 |
  | ((4, 1), (4)) | 4 | 11 | 2 | 9 | 3 |
  | ((4, 1), (3, 1)) | 12 | 8 | 1 | 7 | 2 |
  | ((4, 1), (2, 2)) | 8 | 2 | 0 | 2 | 1 |
  | ((4, 1), (2, 1, 1)) | 12 | 1 | 0 | 1 | 0 |
  | ((3, 2), (4)) | 5 | 5 | 1 | 4 | 2 |
  | ((3, 2), (3, 1)) | 15 | 2 | 0 | 2 | 1 |
  | ((3, 1, 1), (4)) | 6 | 3 | 0 | 3 | 0 |
  | ((3, 1, 1), (3, 1)) | 18 | 1 | 0 | 1 | 0 |
  | ((2, 2, 1), (4)) | 5 | 1 | 0 | 1 | 0 |

  With both level sizes symbolic, the invariant part (gauge-fixed size 7, four Q-vectors) and the two standard parts
  (size 9) would be certificates in three variables (q, x₀, x₁) of the size of the regimes of `code/middle_j2/` that did
  not finish (rows 5b and 10b there). A partial result with one level of fixed size m and the other symbolic (two
  variables) looks more within reach, but needs the same pattern-sum machinery as Theorem 2V, here with three indices
  per state and the inverses of Ê₂ and Ê₃ (or their Woodbury forms on singles and pairs) as functions of the patterns.
- *Larger k.* The cost of a certificate grows by a factor of about 3 to 4 per step in k, mostly the exact inverse of
  Ê₃ over Q(q) (a Woodbury form on pairs, as in Lemma V, would make it smaller); the comparison with `j1lib.py` (about
  100 s per size pair at k = 10) grows as fast. k = 11 is within reach of the same scripts; it was not run.
- Three or more weight values, and general weights.

## Files

| File | Content |
|---|---|
| `l2block.py` | Lemma A and Lemma V for any number field: the blocks Ỹ_PP, Ỹ_PQ, Ỹ_QQ, Ŵ₃ |
| `ufield.py` | exact univariate rational functions over Q (python-flint), Sturm counts, the sign test |
| `young.py` | Young symmetrizers of S_{u₀} × S_{u₁}, the bases of e_λK, e_λW, e_λM_Q, hook dimensions |
| `fixedred.py` | the reduced gauge-fixed forms S_λ (Lemma Y with the vacuum elimination of Lemma V), symmetric elimination |
| `cert_fixed.py` | the certificate for one size pair (Section 4); with bscale ≠ 1 the negative control |
| `coverage.py` | coverage of all two-value weight vectors with 5 ≤ k ≤ 10 by the certified pairs |
| `saved.py` | loads the forms saved by the certificates |
| `check_builder.py` | Lemma A against `j1lib.py`, Lemma V, homogeneity, end-to-end checks at points |
| `check_reduced.py` | the certified forms against projections of the `j1lib.py` block; spectral completeness |
| `check_ylib3.py` | Lemma A and the certified forms against the `ylib.py` block |
| `neg_points.py` | negative controls at explicit points |
| `j1lib.py`, `ylib.py`, `flintfield.py`, `r2lib.py`, `guard.sh` | copied unchanged from `relative-anchors/code/middle_j2/` (`guard.sh` from `no_letter1/`, 2.5 GB limit) |
| `run_all.sh` | reruns everything, one Python process at a time under `guard.sh`; PASS/FAIL per step, ALL PASS at the end (`bash run_all.sh > logs/run_all.out`; KMAX = 10 by default; about 37 minutes) |
| `crosscheck.py`, `run_crosscheck.sh` | the second check of Section 8 (`bash run_crosscheck.sh > logs/run_crosscheck.out`, after `run_all.sh`, whose saved forms it reads) |
| `logs/run_all.out`, `logs/run_*.log`, `logs/cert_k*.pkl` | the final run (the `.pkl` files are not in the repository; `run_all.sh` writes them) |
| `logs/run_crosscheck.out`, `logs/run_crosscheck_*.log` | the second check |
| `logs/dev/` | development runs (superseded by the final run) |

## 8. A second check

`crosscheck.py` was written separately and uses none of `l2block.py`, `young.py`, `fixedred.py`, `ufield.py`,
`cert_fixed.py` or `saved.py`: only the builders `j1lib.py` (Lemma B) and `ylib.py` (from the definitions), `r2lib.py`
(states and kernel vectors), python-flint and numpy.

- *An exact test of positive definiteness.* For a rational symmetric matrix A let A_int be A times the lcm of its
  denominators. A floating-point Cholesky factor L of A gives Z, a rounding of L⁻¹ to integers row by row (one power of
  two per row); Z is lower triangular with nonzero diagonal, so B = Z A_int Zᵀ is congruent to A_int. If B, rescaled by
  these powers of two, is strictly diagonally dominant with positive diagonal (checked in exact integer arithmetic),
  then B, and hence A, is positive definite. Floating point only chooses Z. As negative controls the test refuses the
  full block at k = 5, 6 (it is singular), the gauge-fixed block minus 1.01 λ_min vvᵀ/|v|² for a rationalized
  eigenvector v of its smallest eigenvalue (then vᵀAv < 0 exactly), and the gauge-fixed block with B → (1 + 1/1000)B
  at a two-value point where that block is not positive definite.
- *End to end from Lemma B.* For each of the 27 size pairs with 5 ≤ k ≤ 10, one two-value weight vector with the levels
  at random positions, a weight ratio not used in Sections 2–6 (q between 1/144 and 64) and, for some pairs, an overall
  scale ≠ 1. The block built by `j1lib.py` is symmetric, the C(k, 2) vectors p ⊗ |1@X⟩ (entries w_r on (1@r, X)) are
  in its kernel (exact), and the block without the states (1@0, X) is positive definite (the test above). Every
  vector x splits uniquely as x = v + y with v in the span of the kernel vectors and y zero on the deleted states, and
  xᵀYx = yᵀYy; so the block is ⪰ 0 with kernel exactly span{p ⊗ |1@X⟩} at these 27 points.
- *From the definitions.* The same at 5 points with non-square ratios q (k = 5, 6) for the block built by `ylib.py`
  (in the rescaled frame the kernel vectors are the indicators of the column groups); and `j1lib.py`, rescaled by μ,
  equals 2^{k−2} times `ylib.py` entry by entry at one point each for k = 5, 6.
- *The saved forms.* For all 27 pairs the forms saved by the certificates are symmetric, their irrep dimensions are the
  hook-length dimensions, Σ_λ dim(λ)·(size of the part) = (k − 1)C(k, 2) = dim W, and every part is positive definite
  (exact elimination) at 31 values of q from 10⁻⁸ to 10⁸.

`run_crosscheck.sh` runs all of this, one Python process at a time; the output is `logs/run_crosscheck.out`.
