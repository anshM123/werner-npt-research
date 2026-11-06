# The letter-1 block of (R∞)₂ at two weight values, with one level size symbolic

5 November 2026. Nothing here is refereed.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Theorem L₁₃.** Let k ≥ 11 and let the squared weights take two values, on two levels of sizes a ≥ b with 2 ≤ b ≤ 13 (a = k − b arbitrary, any weight ratio q > 0). Then the letter-1 block of (R∞)₂ (states P_{r,X} = (1@r, 1@X), X a pair; Q_T = (vac, 1@T), T a triple) is ⪰ 0 with kernel exactly span{p ⊗ \|1@X⟩ : X a pair} | **proved**, computer-assisted, exact, symbolic in the size of the other level and in q (Sections 2–5); 12 size regimes × 16 families of parts, every pivot a ratio of polynomials with coefficients of one sign |
| 2 | **Corollary.** The letter-1 block of (R∞)₂ is ⪰ 0 with kernel exactly span{p ⊗ \|1@X⟩} for every k ≥ 5 and every weight vector whose weights take at most two values and whose smaller level has at most 13 copies; in particular for every two-value weight vector with 5 ≤ k ≤ 27 | **proved**, given row 1 and, for the other cases, Theorem L (5 ≤ k ≤ 10, `code/middle_j2_l1/`), Theorem E1 (one level) and Theorem MB (a level of size 1, `code/marked_copy/`) of the README (all computer-assisted, not refereed) (Section 5) |
| 2b | **Theorem L₄₀.** The statement of row 1 for 2 ≤ b ≤ 40. Hence the letter-1 block of (R∞)₂ holds at two weight values for every k ≥ 5 whose smaller level has at most 40 copies, in particular for 5 ≤ k ≤ 81 | **proved**, computer-assisted, exact: the regimes (b+, b) for 14 ≤ b ≤ 40 (`run_ext30.sh`, `run_ext40.sh`; 27 × 16 certificates), forms, coverage for k ≤ 81 and negative controls (`run_ext_checks.sh`) (Section 6b) |
| 3 | **Corollary.** (R∞)₂ ⪰ 0 with kernel exactly p ⊗ H₂ (a) for every two-value weight vector with 5 ≤ k ≤ 26, and (b) for every k ≥ 5 and every two-value weight vector whose smaller level has at most 12 copies | **proved**, given row 2 and Theorems F1b, F2b (`code/middle_j2/`), D₂/DJ and Proposition T2′ of the README (Section 5) |
| 4 | Lemma 1 (the P–P entries with Ê₂⁻¹ in Woodbury form through the singles), Lemma 2 (Ỹ_PQ = (I ⊗ Â₂⁻¹)K with an inverse-free K), Lemma 3 (the congruence S′ = (I ⊗ Â₂)S(I ⊗ Â₂) has the same kernel) | **proved** by hand (Section 2); **exact check** entry by entry against the Lemma B builder (`fastj1.py`, the fast copy of `j1lib.py`) at 17 two-level points (15 size pairs) with k = 5, 6, 7, 8, 11, 12 (`check_points.py`) |
| 5 | Lemma 4 (the reduction: tableau test vectors, bases of the parts, products through the reduced bases, pattern sums with falling-factorial counts) | **proved** (Section 3); **exact check**: at the same 17 points every reduced form (Gram matrices and the whole reduced S′ on the K- and W-vectors, for each of the 16 families) equals the projection of the explicit S′ built from the Lemma B block; the dimension counts add up (`check_points.py`; at k = 11, 12 with the bases of the certified regimes) |
| 6 | The certified symbolic forms (whole W-blocks) equal the same reduction in exact rational arithmetic in 603 comparisons (q ∈ {1/9, 4, 49/25}, three values of the symbolic size, every family) in the regimes (9+, 2), (8+, 3), (7+, 4), (6+, 5), (13+, 13) | **exact check** (`check_forms.py`) |
| 7 | Coverage: every two-value configuration with 11 ≤ k ≤ 80 whose smaller level has at most 13 copies lies in a certified regime (or has a level of size 0 or 1: Theorems E1, MB), with all dimension counts right; the configurations with both levels ≥ 14 (smallest k = 28) are listed as outside the statement | **exact** bookkeeping (`coverage.py`) |
| 8 | Negative controls: with B → (1 + 1/1000)B the reduced forms of the parts inv, mix, std₀, std₁ (and two₁ at some sizes) are not positive definite at q = 1/10⁴ at (u₀, u₁) = (9, 2), (8, 3), (7, 4), (6, 5), (12, 9), (21, 20), while all parts are positive definite for the true B; at (3, 2) and (4, 2) these scaled reduced forms equal the projections of the scaled explicit block | **exact** (`neg.py`) |
| 9 | Both level sizes symbolic (both levels ≥ 14; after row 2b needed only when both levels are ≥ 41, k ≥ 82) | **open**. Attempted in three variables (q, x₀, x₁) at (6+, 6+): the part anti₀ passes; the part mix passed 4 of its 7 pivots (pivot 4: numerator degree 80 with 12 299 terms) and was stopped (`logs/dev/dev_cert2_*.log`). The parts inv, std₀, std₁ (7 and 9 pivots) were not attempted in three variables. A determinant-and-inertia certificate (Proposition C of the README) or a larger M are the next steps |
| 10 | Three or more weight values; general weights | **open**; not attempted |

Nothing here says anything about Conjecture 12 of `single-excitation-range/` in general, or about the NPT (Werner
distillability) problem. Computer-assisted means: exact arithmetic with rational functions in python-flint and sign
checks of coefficients; no floating point is used anywhere.

## 1. Setting and statement

Notation as in `relative-anchors/README.md`, Sections 1 and 2, and `code/middle_j2_l1/NOTES.md`, Sections 1 and 2
(repository `werner-npt-research`). Level 0 consists of u₀ copies with squared weight 1, level 1 of u₁ copies with
squared weight q > 0; L = u₀ + qu₁. For a set Z of copies, π_Z, σ_Z are the product and sum of its squared weights. We
use Lemma A of `middle_j2_l1` (entries of Lemma B in units 2^{k−2}, rescaled frame Ỹ = DYD) and Lemma V there
(Ỹ_QQ = 2Ŵ₃⁻¹ for k ≥ 5). The kernel vectors p ⊗ |1@X⟩ become multiples of the indicators 1_X of the column groups; let
K = span{1_X} ⊂ M_P (M_P, M_Q: the spaces of P- and Q-states).

**Theorem L₁₃.** Let k ≥ 11 and let the squared weights take two values, on two levels of sizes a ≥ b with
2 ≤ b ≤ 13. Then for every weight ratio the letter-1 block of (R∞)₂ is ⪰ 0 with kernel exactly
span{p ⊗ |1@X⟩ : X a pair}.

By homogeneity (Ỹ is homogeneous of degree 3 in the squared weights) and a permutation of the copies, the level of size
b can be taken as level 1 (squared weight q, any q > 0) and the other as level 0 (squared weight 1). Since a ≥ b and
a + b ≥ 11, u₀ = a ≥ o(b) := max(b, 11 − b), so the configuration lies in the regime (u₀ ≥ o(b), u₁ = b).

## 2. Closed forms

On the letter-1 strings with c letters, Ê_c = Δ_c + ½ j̄⁺Λ_{c−1}j̄⁻ with Δ_c = diag(2(L − σ_Z)/π_Z) and
Λ_{c−1} = diag(1/π_Y) on the (c − 1)-sets (this is the formula for Ê_c in Lemma A of `middle_j2_l1`).

**Lemma 1 (P–P).** Let N = j̄⁻Δ₂⁻¹j̄⁺ and Ψ = diag(s_c) + ½N on the singles (N[c, c′] = 1/Δ₂({c, c′}) for c ≠ c′,
N[c, c] = Σ_{t≠c} 1/Δ₂({c, t})). Then Ê₂⁻¹ = Δ₂⁻¹ − ½Δ₂⁻¹j̄⁺Ψ⁻¹j̄⁻Δ₂⁻¹, and the B-term of Lemma A item 1,
s_r s_{r′}𝔅 with 𝔅 = ⟨t_{r′}X′, Ê₂⁻¹t_rX⟩ (row state (r, X′), column state (r′, X)), is
- r ∈ X = {r, c}, r′ ∈ X′ = {r′, c′}: 𝔅 = ¼N[c′, c] − ⅛(NΨ⁻¹N)[c′, c];
- r ∈ X = {r, c}, r′ ∉ X′: 𝔅 = ([c ∈ X′] − ½Σ_{y∈X′}(Ψ⁻¹N)[y, c])/Δ₂(X′), and symmetrically;
- r ∉ X, r′ ∉ X′: 𝔅 = 4[X = X′]/Δ₂(X) − 2Σ_{x∈X, y∈X′}Ψ⁻¹[y, x]/(Δ₂(X)Δ₂(X′)).

*Proof.* Woodbury with Ê₂ = Δ₂ + ½j̄⁺Λ₁j̄⁻ and Λ₁⁻¹ = diag(s). Write t_r = ½j̄⁺σ_r + 2(I − n_r) (σ_r removes r):
t_r|X⟩ = ½j̄⁺e_c if X = {r, c} and 2|X⟩ if r ∉ X, and j̄⁻Δ₂⁻¹t_r|X⟩ = ½Ne_c, respectively 2Σ_{x∈X}e_x/Δ₂(X). Insert
into ⟨h′, Δ₂⁻¹h⟩ − ½(j̄⁻Δ₂⁻¹h′)ᵀΨ⁻¹(j̄⁻Δ₂⁻¹h). ∎

At two weight values N, Ψ and their products and inverse are *two-level operators* M[c, c′] = [c = c′]a_{l(c)} +
b_{l(c)l(c′)}; they form an algebra with (a, b)(a′, b′) = (aa′, diag(a)b′ + b diag(a′) + b diag(u)b′) and
(a, b)⁻¹ = (1/a, −(diag(a) + b diag(u))⁻¹ b diag(1/a)). Here a_l > 0 when u_l ≥ 2 (a_l is the eigenvalue on the standard
part of level l), and diag(a) + b diag(u) is the invariant 2 × 2 part of Ψ ≻ 0. So every P–P entry is an explicit
rational function of q, u₀, u₁ and of the equality pattern and levels of its six indices.

**Lemma 2 (P–Q without inverses on triples).** Let Â₂ = D⁻¹(I + ½J⁻C₃⁻¹J⁺)D⁻¹ on the pairs:
Â₂[Z, Z] = (1 + ½Σ_{t∉Z} s_t/(2(L − σ_Z − s_t)))/π_Z, Â₂[Z, Z′] = 1/(4π_{Z∩Z′}(L − σ_{Z∪Z′})) if |Z ∩ Z′| = 1, else 0
(the analogue of Ŵ₃ of Lemma V one level down; Â₂ ≻ 0). Let K[(r, Z), T] = −s_r[r ∈ T][Z = T − r] +
s_r[r ∉ T][Z ⊂ T] s_{T∖Z}/(L − σ_T). Then Ỹ_PQ = (I ⊗ Â₂⁻¹)K.

*Proof.* (I + ½j̄⁻Δ₃⁻¹j̄⁺Λ₂)j̄⁻ = j̄⁻Δ₃⁻¹Ê₃ gives j̄⁻Ê₃⁻¹ = Â₂⁻¹Λ₂j̄⁻Δ₃⁻¹, since Λ₂⁻¹ + ½j̄⁻Δ₃⁻¹j̄⁺ = Λ₂⁻¹Â₂Λ₂⁻¹. For
r ∉ T, t_r|T⟩ = 2|T⟩ and Lemma A item 2 gives s_r Σ_{t∈T} s_t Â₂⁻¹[X′, T − t]/(L − σ_T). For r ∈ T,
t_r|T⟩ = ½j̄⁺|T − r⟩ and ½Λ₂j̄⁻Δ₃⁻¹j̄⁺ = Â₂Λ₂⁻¹ − I, so the E⁻¹-term is s_rπ_{T−r}[X′ = T − r] − s_rÂ₂⁻¹[X′, T − r];
the first part cancels the A-part −[r ∈ T][X′ = T − r]π_T. ∎

**Lemma 3 (congruence).** With S = Ỹ_PP − ½Ỹ_PQŴ₃Ỹ_QP (the Schur complement of the vacuum states, Lemma V) and
M = I ⊗ Â₂ (positive definite),

    S′ := MSM = (I ⊗ Â₂)Ỹ_PP(I ⊗ Â₂) − ½ K Ŵ₃ Kᵀ.

S′ ⪰ 0 iff S ⪰ 0, and ker S′ = M⁻¹ker S. Since M⁻¹(1 ⊗ v) = 1 ⊗ Â₂⁻¹v and K = 1 ⊗ M_pairs, M⁻¹ maps K onto K; so
ker S′ = K iff ker S = K. Hence the block is ⪰ 0 with kernel exactly K iff S′ is (Lemma V and `middle_j2_l1`, Section 2).
Every operator in S′ is explicit; the only inverse left is the two-level inverse Ψ⁻¹ on the singles. ∎

Checks (`check_points.py`, exact, at every run point): every P–P entry of the Lemma B block equals the formula of
Lemma 1, (I ⊗ Â₂)Ỹ_PQ = K and Ỹ_QQŴ₃ = 2I, with Â₂, K, Ŵ₃ built from their formulas on the concrete copies.

## 3. The reduction

G = S_{u₀} × S_{u₁} acts on the copies of each level and commutes with every operator above. Only irreps
λ = (λ₀, λ₁) with at most three boxes outside the first rows occur in M_P and M_Q, sixteen families: inv;
std_l = (u_l − 1, 1) on level l; two_l = (u_l − 2, 2); anti_l = (u_l − 2, 1, 1); three_l = (u_l − 3, 3);
tto_l = (u_l − 3, 2, 1); mix = (std₀, std₁); twostd_l = (two_l, std_{1−l}); antistd_l = (anti_l, std_{1−l}) (l = 0, 1;
(u − 3, 1, 1, 1) does not occur, since the P-states of one level are the Young modules of shapes (u − 2, 1, 1) and
(u − 3, 2, 1)).

*Tableaux and test vectors.* For each family fix a tableau T on each level whose first row contains all copies but those
of the other rows; the copies in columns of length ≥ 2 are *named* (at most six), the others anonymous. R_T is the row
group, C_T the column group (signed sum), e′_λ = C_T R_T. As for the Young symmetrizer e_λ = R_TC_T (von Neumann's lemma,
`middle_j2_l1` Section 3), e′_λ x e′_λ ∈ C e′_λ, so on any G-module V, e′_λV has dimension mult(λ, V) and consists of
vectors s₀ ⊗ c for one fixed s₀ ∈ λ. Since R_TV is spanned by the indicators 1_O of the R_T-orbits O of basis
vectors, e′_λV is spanned by the vectors C_T1_O. W := span{P_{r,X} : r ∉ X} + span{P_{a,X} − P_{b,X} : X = {a, b}} is a
G-stable complement of K (Lemma Y of `middle_j2_l1`), and
- e′_λK is spanned by C_Tι(1_{O′}), O′ an R_T-orbit of pairs, ι(e_X) = 1_X;
- e′_λW is spanned by C_T1_O (O an orbit of states with r ∉ X) and C_T(1_{O₁} − 1_{O₂}) (O₁, O₂ the two orbits of
  states (x₁, {x₁, x₂}) and (x₂, {x₁, x₂}) when x₁, x₂ lie in different R_T-orbits of copies);
- e′_λM_Q is spanned by C_T1_O, O an R_T-orbit of triples.

`r2v.choose` takes a maximal independent subset of each list (exact Gram ranks at a sample size of the regime); B_P is
the K-part followed by the W-part, B_Q the Q-part.

**Lemma 4 (reduction).** Let R_X = B_PᵀXB_P etc. denote reduced forms, G_P = B_PᵀB_P, G_Q = B_QᵀB_Q. If B_P and B_Q are
bases of e′_λM_P and e′_λM_Q, then

    B_PᵀS′B_P = R_{I⊗Â₂} G_P⁻¹ R_{Ỹ_PP} G_P⁻¹ R_{I⊗Â₂} − ½ R_K G_Q⁻¹ R_{Ŵ₃} G_Q⁻¹ R_Kᵀ,

its rows on the K-part vanish, and the block is ⪰ 0 with kernel exactly K iff, for every λ, the W-block S′_λ of this
matrix is positive definite.

*Proof.* For equivariant X, Y and z ∈ e′_λM_P = span(B_P): z = B_PG_P⁻¹B_Pᵀz, and YB_P ⊂ e′_λM_P (Y commutes with
e′_λ); so B_PᵀXYB_P = R_XG_P⁻¹R_Y, and the same through M_Q. S′ vanishes on K, so the K-rows and K-columns vanish, and
the reduced matrix is diag(0, S′_λ). On the λ-isotypic component S′ is ⟨·,·⟩_λ ⊗ Q_λ and the reduced matrix is
|s₀|²CᵀQ_λC with C invertible; so Q_λ ⪰ 0 with kernel of dimension mult(λ, K) iff S′_λ ≻ 0. Summing over λ, ker S′ has
dimension dim K and contains K (Lemma Y of `middle_j2_l1`, with e′_λ in place of e_λ). ∎

*Pattern sums.* Each reduced matrix is a sum over pairs of states of (coefficient) × (entry) × (coefficient). The
coefficients of C_T1_O at a state depend only on the state's indices relative to the named copies (which named copy, or
an anonymous copy of level 0 or 1, and which indices coincide); the entries of Ỹ_PP (Lemma 1), I ⊗ Â₂, K and Ŵ₃ depend
only on the equality pattern and levels of the six indices of the two states. So each reduced entry is
Σ over joint patterns of count × coefficients × entry, where a joint pattern fixes for every index a named copy or an
anonymous one (anonymous indices of the two states may coincide), and count = Π_l (u_l − n_l)(u_l − n_l − 1)⋯
(u_l − n_l − m_l + 1) with n_l named and m_l distinct anonymous copies of level l. States are enumerated as ordered
tuples with the factors ½ (P) and 1/6 (Q). The sums inside the entries over the other copies (in N, Ψ, Â₂, Ŵ₃) are
taken by levels, (u_l − n_l(Z)) times the value for a copy of level l. (`r2v.joint_sum`, `r2v.assemble`.)

*Specialisation.* The certificates compute over Q(q, x) with u₀ = o + x (o ≥ 6) and u₁ = m. A joint pattern that needs
more than m copies of level 1 has count exactly 0 and is skipped. Every other denominator is positive at every real
x ≥ 0 and q > 0: L − σ_Z ≥ u₀ − 4 > 0 for the sets Z (at most four copies) of a pattern, s_l, Δ₂, the two-level
quantities a_l (u_l ≥ 2) and det(diag(a) + b diag(u)), det G_P and det G_Q (checked positive on the regime by the
certificate), and the pivots of the elimination (checked positive). So the symbolic forms specialise to the true
reduced forms at every integer size of the regime.

## 4. Certificates

`cert.py FAMILY U0 U1` (and `batch.py`, all sixteen families for one regime), over Q(q, x) in python-flint:
1. the bases at the sample size (o + 3, m);
2. det G_P and det G_Q pass the sign test (every irreducible factor of numerator and denominator has coefficients of one
   sign and a monomial free of x), so the vectors are bases of the multiplicity spaces on the whole regime;
3. the reduced S′ (K- and W-part) is symmetric and its K-rows vanish identically;
4. symmetric elimination of the W-block S′_λ: every pivot (the ratio of consecutive leading principal minors) passes
   the sign test with sign +1, so S′_λ ≻ 0 for every q > 0 and every real x ≥ 0 (Sylvester).

Regimes: (o(m)+, m) for 2 ≤ m ≤ 13, o(m) = max(m, 11 − m): (9+, 2), (8+, 3), (7+, 4), (6+, 5), (6+, 6), (7+, 7), …,
(13+, 13). In all 12 × 16 runs every pivot passes the sign test. Sizes (from `logs/run_cert_*.log`; parts with an
empty W-part at small m print sizes 0 there): the gauge-fixed parts have size at most 9 (std₀, std₁; inv 7 or 6, mix 7
or 5, two 4 or 3, anti 3 or 2, twostd 2, three, tto, antistd 1), with at most 4 Q-vectors; the largest pivot numerator has 2 143 terms (degree 64 in (q, x)). Each regime
(sixteen certificates) takes 10 to 70 seconds.

**Proof of Theorem L₁₃.** By Section 1 every configuration is in a regime (o(m)+, m) with m = b. In it, Lemmas 3 and 4
reduce the statement to the positive definiteness of the sixteen W-blocks, which the certificates prove for every q > 0
and every size u₀ ≥ o(m); the bases are complete at every size by the dimension count of `coverage.py` (Σ_λ dim λ ·
(number of vectors) equals k·C(k, 2), C(k, 2) and C(k, 3) for the P-states, K and the Q-states: an identity between
polynomials of degree at most 3 in u₀ (all shapes are partitions for u₀ ≥ 6), checked at every u₀ with k ≤ 80, at
least 50 values per regime). ∎

## 5. Coverage and consequences

`coverage.py 13 none 60` (and with 80): every two-value configuration with 11 ≤ k ≤ 60 (80) whose smaller level has at
most 13 copies is in a certified regime, or has a level of size 0 (Theorem E1) or 1 (Theorem MB, j = 2 ≤ 8), and all
dimension counts are right; the configurations with both levels ≥ 14 (289 with k ≤ 60, the smallest at k = 28) are
listed as outside the statement.

**Corollary (row 2).** For k ≥ 11 the cases are: smaller level 0 (equal weights): Theorem E1; 1: Theorem MB; 2 to 13:
Theorem L₁₃. For 5 ≤ k ≤ 10: Theorem L. If k ≤ 27, the smaller level has at most 13 copies.

**Corollary (row 3).** At two weight values the torus blocks of (R∞)₂ are the letter-1 block (row 2), the type I blocks
with |T| = 2 and |T| = 1 (Theorems F1b and F2b of `code/middle_j2/`: every two-value configuration with k ≤ 26, and for
larger k whenever one U-level has at most 12 copies, which holds if the smaller level has at most 12 copies), those with
|T| = 3 (Theorem DJ) and the type II blocks (Proposition T2′). Each is ⪰ 0 with kernel exactly the part of p ⊗ H₂ it
contains. This rests on those theorems as stated there; they are computer-assisted and not refereed, like Theorem L₁₃.

## 6. Checks and controls

All exact.
- `check_points.py U0 U1 Q… [sample=…]`: the explicit block from `fastj1.py` (the Lemma B builder, same output as
  `j1lib.py`), rescaled; the identities of Lemmas 1 and 2 and Lemma V entry by entry; S_ref = Ỹ_PP − Ỹ_PQỸ_QQ⁻¹Ỹ_QP
  with an exact inverse, S′_ref = (I ⊗ Â₂)S_ref(I ⊗ Â₂); for each family the test vectors realised on concrete copies,
  the Gram matrices and the whole reduced S′ (K- and W-part) equal `reduce_family` at the same numbers, the K-rows
  vanish, and the dimension counts hold. Run at (3, 2), (2, 3), (4, 2), (3, 3), (2, 4), (5, 2), (4, 3), (3, 4), (6, 2),
  (5, 3), (4, 4) (k = 5 to 8, 13 points), and inside the certified regimes, with the bases of the certificates, at
  (9, 2), (7, 4), (6, 5) (k = 11) and (8, 4) (k = 12).
- `check_forms.py`: the forms saved by the certificates (regimes (9+, 2), (8+, 3), (7+, 4), (6+, 5), (13+, 13)),
  evaluated exactly at q ∈ {1/9, 4, 49/25} and x ∈ {0, 1, 5}, equal the reduction in rational arithmetic with the
  same bases (603 comparisons of whole W-blocks). With the previous item this ties the certified objects to the Lemma B
  block at k = 11, 12. (In the first complete run this step failed because the saved coefficient strings were
  evaluated with float division; after the fix the evaluation is exact; `logs/dev/first_run/`.)
- `neg.py`: B → (1 + 1/1000)B (every E⁻¹-term of Lemma B times 1001/1000). The scaled Schur complement is reduced
  without the congruence (`reduce_scaled`: R_{Ỹ_PQ} = G_PR_{I⊗Â₂}⁻¹R_K, reduced inverses as GR⁻¹G); at (3, 2) and
  (4, 2) it equals the projection of the scaled explicit block, for b = 1 and 1001/1000. Inside the regimes, at
  (9, 2), (8, 3), (7, 4), (6, 5), (12, 9), (21, 20) and q = 1/10⁴, all parts are positive definite for the true B and
  the parts inv, mix, std₀, std₁ (and two₁ from (7, 4) on) are not for the scaled one.
- `coverage.py`: Section 5.

## 6b. The regimes (b+, b) for 14 ≤ b ≤ 40 (`run_ext30.sh`, `run_ext40.sh`, `run_ext_checks.sh`)

**Theorem L₄₀.** The statement of Theorem L₁₃ holds for 2 ≤ b ≤ 40.

*Proof.* For 14 ≤ b ≤ 40 the offset is o(b) = b, and `batch.py "b+,b"` runs the sixteen certificates of Section 4
in the regime (b + x, b), x ≥ 0 real. All 27 × 16 pass (`logs/ext_cert_*.log`; each regime about a minute). The
bases are complete at every size by the dimension counts of `coverage.py 40 none 81`, which also checks that every
two-value configuration with 11 ≤ k ≤ 81 lies in a certified regime or has a level of size 0 or 1. ∎

Checks (`run_ext_checks.sh`, ALL PASS): the saved forms of the regimes (14+, 14), (20+, 20), (27+, 27), (33+, 33),
(40+, 40) equal the rational-arithmetic reduction at the points of `check_forms.py`; `neg_ext.py`: at (u₀, u₁) =
(30, 14), (25, 20), (41, 27), (50, 33), (41, 40) every W-block is positive definite for the true B, and with
B → (1 + 1/1000)B the parts inv, mix, std₀, std₁ and two₁ are not, at q = 1/10⁴.

## 7. What is not covered

- Both levels larger than 40 (k ≥ 82; before Section 6b: larger than 13, k ≥ 28). The reduction of Section 3 is symbolic in both sizes (the code takes u₁ = o + x₁
  as well), but the certificates in three variables (q, x₀, x₁) grow fast: at (6+, 6+) the part mix passed 4 of its 7
  pivots (the fourth with 12 299 terms) before it was stopped; the parts inv, std₀, std₁ have 7, 9, 9 pivots. A
  determinant-and-inertia certificate (Proposition C of the README) or a larger fixed range M (each regime takes under a
  minute, so M can be raised at a linear cost, as Section 6b does) are the obvious next steps. Theorem L₄₀ is not
  claimed beyond b ≤ 40.
- The exclusions of Theorems F1b and F2b (both U-levels larger than 40 after `middle_j2/run_fext.sh`), so (R∞)₂ at
  two weight values for k ≥ 82 with both levels larger than 40.
- Three or more weight values, and general weights.

## Files

| File | Content |
|---|---|
| `r2v.py` | the closed forms (two-level operators, Lemma 1, Â₂, K, Ŵ₃), the sixteen families, tableau test vectors, bases, pattern sums, the reduced S′ (`reduce_family`), the scaled reduction for the controls (`reduce_scaled`) |
| `cert.py` | certificate for one family and one regime |
| `batch.py` | the sixteen certificates of one or more regimes in one process; saves the forms (`logs/forms_*.pkl`) |
| `run_ext30.sh`, `run_ext40.sh`, `run_ext_checks.sh`, `neg_ext.py` | Section 6b: the regimes (b+, b) for 14 ≤ b ≤ 30 and 31 ≤ b ≤ 40 (`logs/ext_cert_*.log`, `logs/run_ext30.out`, `logs/run_ext40.out`), then the forms, the coverage for k ≤ 81 and the negative controls (`logs/run_ext_checks.out`) |
| `check_points.py` | the reduction against the explicit Lemma B block (Sections 2, 3) |
| `check_forms.py` | the saved symbolic forms against the rational-arithmetic reduction |
| `neg.py` | negative controls |
| `coverage.py` | coverage of all configurations, dimension counts |
| `fastj1.py`, `j1lib.py` | the Lemma B builders (copied unchanged from the scratch folder `fastj1/` and from `code/middle_j2_l1/`; they are identical copies of `j1lib.py`) |
| `flintfield.py`, `guard.sh`, `l2block.py`, `r2lib.py` | copied unchanged from `code/middle_j2_l1/` (`l2block.py` is not used by the scripts here) |
| `run_all.sh` | reruns everything, one Python process at a time under `guard.sh` (2.5 GB); PASS/FAIL per step, ALL PASS at the end; output `logs/run_all.out`, one log `logs/run_<step>.log` per step |
| `logs/dev/` | development runs, including the stopped three-variable attempts of row 9; `logs/dev/first_run/`: the first complete run (all certificates PASS; `forms` FAIL from float evaluation of the saved strings, `coverage` FAIL because it counted the configurations outside the statement; both scripts fixed, then the full rerun in `logs/run_all.out`) |
