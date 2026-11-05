# The blocks without a letter 1 of every (R∞)_j, for every number of copies

Ansh Mishra and Aryan Senthilkumar, 5 November 2026. Nothing here is refereed.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Theorem DJ.** For every j ≥ 1, every k ≥ j + 1, every positive weight vector and every set T of j + 1 copies, the type I block of (R∞)_j in which the copies of T carry the letter 2 (size j + 2, no letter 1) is positive definite (Section 7) | **proved**, computer-assisted: closed form, vacuum elimination and the kernel representations by hand (Sections 2–4); the last steps are two exact coefficient-sign certificates, `cert_pair.py` (7 092 terms, one AM–GM repair) and `cert_row.py` (357 055 terms, all coefficients ≥ 0), each a few seconds |
| 2 | Lemma 1 (closed form of the block, every j and k) | **proved** by hand from Lemma B; **exact check**: equal to 2^{k−j} times the ylib block, entry by entry, for 51 blocks (all (k, j) with 1 ≤ j ≤ 6, j + 1 ≤ k ≤ 7, and (k, j) = (8, 1), (8, 3), (8, 5), (8, 7)); the A-part and the scaled block A − (1001/1000)B also agree (`check_closed.py`) |
| 3 | Lemma 2 (vacuum elimination; the formulas for F_ii, F_il and H_il) | **proved** (algebra); **exact check** at 12 rational configurations, every ordered pair (`check_reduction.py`) |
| 4 | Lemma 3 (kernel representations of D_il and N_il) | **proved** (multilinear expansion); **exact check** with the very objects certified in rows 5 and 6, at 9 rational configurations with 1 ≤ \|U\| ≤ 4 (`check_reduction.py`) |
| 5 | Proposition A (inequality (i)): F_il < 0 for all i ≠ l | **proved**: Lemma 3 and `cert_pair.py` |
| 6 | Proposition B (row inequality): a_l F_ii + (x_T − a_i) F_il > 0 for all i ≠ l | **proved**: Lemma 3 and `cert_row.py` |
| 7 | The conclusions at rational weight vectors: block positive definite, F_il < 0, strict diagonal dominance, Proposition B; \|T\| = 2, …, 8, \|U\| ≤ 12, random, corner, tiny, huge and spread weights (161 vectors); blocks taken directly from ylib positive definite for k ≤ 6 (24 vectors) | **exact check** (`check_pd.py`) |
| 8 | Negative control: with B replaced by (1 + 1/1000)B the block is indefinite at 12 explicit rational points (j = 1, 2, 3, 5; k ≤ 9), next to the corner where one copy of T carries almost all the weight, and Propositions A and B fail there for the scaled block; for k ≤ 5 the scaled block from ylib equals the scaled closed form | **exact check** (`negative_control.py`) |
| 9 | Margins: F_il/(a_i a_l) ≤ −0.16, (a_l F_ii + (x_T − a_i)F_il)/(a_l(c − a_i)) ≥ 0.80 | **numerical** (`numerics.py`, 6 000 random configurations, \|T\| ≤ 14, \|U\| ≤ 40) |
| 10 | Corollary: with Proposition T2′ of `relative-anchors/README.md` (type II blocks), every block of every (R∞)_j, j ≥ 1, that contains no letter 1 is positive definite, for every k ≥ max(j, 2) and every positive weight vector | **proved**, given Proposition T2′ (not re-proved here) |
| 11 | Separate implementations, written from the formulas of Sections 3 to 6 and not from `certlib.py` or `djlib.py`: both kernel certificates (the same polynomials: 7 092 terms with the single negative coefficient at p q e¹⁰ s₁² s₂², and 357 055 terms, all nonnegative), the two U = ∅ numerators (6 and 20 terms, all positive), and the conclusions (block positive definite, F_il < 0, strict diagonal dominance) at 60 further random rational configurations with \|T\| ≤ 12, \|U\| ≤ 15 | **exact check** (`independent_cert.py`, `independent_e2e.py`) |

Theorem DJ contains Proposition U0 (k = j + 1), Theorem D₁ (j = 1) and Theorem D₂ (j = 2) of the README, with a
different proof. Nothing here says anything about Conjecture 12 of `single-excitation-range/` in general, or about the
NPT (Werner distillability) problem. The letter-1 blocks and the type I blocks that contain a letter 1 of the middle
conditions (R∞)_j, 1 ≤ j ≤ k − 2, are untouched and remain open for k ≥ 5 (Conjectures L and L1 for j = 1).

## 1. Setting and notation

Notation as in `relative-anchors/README.md`, Sections 1 and 2 (Lemmas 1, B, N, G). There are k copies with weights
w_t > 0, s_t = w_t², and we normalise |s| = 1, so x_t = s_t. Fix j ≥ 1 and a set T = {t_0, …, t_j} of j + 1 ≤ k
copies; put U = [k] ∖ T.

A type I state of (R∞)_j without a letter 1 has j + 1 letters 2 on distinct copies (the row letter and the j column
letters, or the j + 1 column letters of a vacuum state). So the type I blocks of (R∞)_j without a letter 1 are exactly
the blocks indexed by a set T with |T| = j + 1, and such a block has the j + 2 states

- X_i = (2@t_i, 2 on T ∖ t_i), i = 0, …, j (row: the letter 2 at t_i; column: the letter 2 on T ∖ t_i, zeros elsewhere),
- X_vac = (vac, 2 on T).

Put a_i = x_{t_i} > 0, x_T = Σ_i a_i, e = x_U = 1 − x_T ≥ 0 (e = 0 exactly when k = j + 1), c = 2 − x_T, and for
u ∈ U

    t_u = c − 2x_u = 2 − 2x_u − x_T,    Φ(d) = Σ_{u∈U} x_u/(t_u + d)   (Φ = 0 if U = ∅).

Since 0 < x_u ≤ e and c − 2e = x_T, every t_u lies in [x_T, c). In homogeneous form (no normalisation) c = x_T + 2e,
and every quantity below is homogeneous of degree 0.

## 2. The closed form

**Lemma 1.** In units 2^{k−j} (that is, for M = 2^{k−j} Y restricted to the block), in the basis X_0, …, X_j, X_vac:

- M_ii = 2;
- for i ≠ l, M_il = (1 + G_i)(1 + G_l)/(2 + G_i + G_l + Φ(a_i + a_l)), with G_i = a_i/(c − a_i + a_l) and
  G_l = a_l/(c − a_l + a_i) (they depend on the pair; c − a_i + a_l = 2 − 2a_i − x_S with x_S = x_T − a_i − a_l);
- M_{i,vac} = −A_i, A_i = 2(1 + g_i)/(2 + g_i + Φ(a_i)), g_i = a_i/(c − a_i);
- M_{vac,vac} = 2/Γ, Γ = 1 + Φ(0)/2.

*Proof.* Write Ψ = Σ_i y_i e_{a_i} ⊗ ψ_i + q vac ⊗ φ with a_i = (2 at t_i), ψ_i = |2 on T ∖ t_i⟩ and φ = |2 on T⟩, and
evaluate the quadratic form 2^Z⟨Ψ, YΨ⟩ of Lemma B (Z = k − j) term by term.
- *Swap term.* For i ≠ l, σ_{a_l}⁻ψ_i = σ_{a_i}⁻ψ_l = |2 on T ∖ {t_i, t_l}⟩, so each ordered pair gives ½ y_i y_l; for
  i = l, σ_{a_i}⁻ψ_i = 0. So the A-part is ½ off the diagonal.
- *Diagonal term.* ψ_i has no letter at t_i, so it gives 2y_i².
- *B-term.* 𝒯_{a_i}ψ_i = 0, since both parts of 𝒯_{(t,2)} = ½J⁺σ⁻_{(t,2)} − w_t|1⟩⟨2|_t need a letter 2 at t_i. For
  l ≠ i put S = T ∖ {t_i, t_l}. Then 𝒯_{a_i}ψ_l = Σ_{u∉S} (½w_u − w_{t_i}[u = t_i]) |1@u, 2 on S⟩ =: v_i, a vector on
  the strings |1@u, 2 on S⟩ (u ∈ U ∪ {t_i, t_l}), which lie in H_j. On these strings
  E_j = C + ½J⁺J⁻ = diag(d_u) + ½w′w′ᵀ with d_u = 2(1 − x_u − x_S) + x_S = 2 − 2x_u − x_S, and w′ = (w_u)_{u∉S}.
  The term is −Σ_{i≠l} y_i y_l ⟨v_l, E_j⁻¹ v_i⟩ (here 𝒯_{a_l}ψ_i = v_l). Put Δ = diag(d_u) and
  γ = Σ_{u∉S} x_u/d_u = ⟨w′, Δ⁻¹w′⟩. Since x_{t_i}/d_{t_i} = G_i, x_{t_l}/d_{t_l} = G_l and d_u = t_u + a_i + a_l for
  u ∈ U, γ = G_i + G_l + Φ(a_i + a_l). Moreover ⟨w′, Δ⁻¹v_i⟩ = ½γ − G_i and ⟨v_l, Δ⁻¹v_i⟩ = ¼γ − ½(G_i + G_l), so
  Sherman–Morrison gives
  ⟨v_l, E_j⁻¹ v_i⟩ = ¼γ − ½(G_i + G_l) − ½(½γ − G_i)(½γ − G_l)/(1 + ½γ), and a short computation gives
  ½ − ⟨v_l, E_j⁻¹v_i⟩ = ½(1 + G_i)(1 + G_l)/(1 + ½γ). This is M_il.
- *Vacuum coupling.* The term is −2Σ_i y_i q (⟨ψ_i, σ_{a_i}⁻φ⟩ − ⟨ψ_i, J⁻E_{j+1}⁻¹𝒯_{a_i}φ⟩). Here σ_{a_i}⁻φ = ψ_i,
  𝒯_{a_i}φ = Σ_{u∈U∪{t_i}} (½w_u − w_{t_i}[u = t_i]) |1@u, 2 on T ∖ t_i⟩, E_{j+1} = diag(2 − 2x_u − (x_T − a_i)) +
  ½w″w″ᵀ on these strings, and J⁻ maps |1@u, 2 on T ∖ t_i⟩ to w_u ψ_i. With γ″ = g_i + Φ(a_i) (the u = t_i term is
  a_i/(2 − a_i − x_T) = g_i, and 2 − 2x_u − x_T + a_i = t_u + a_i), Sherman–Morrison gives
  ⟨w″, E⁻¹𝒯_{a_i}φ⟩ = (½γ″ − g_i)/(1 + ½γ″), so the bracket is 1 − (½γ″ − g_i)/(1 + ½γ″) = (1 + g_i)/(1 + ½γ″), and
  M_{i,vac} = −(1 + g_i)/(1 + ½(g_i + Φ(a_i))) = −A_i.
- *Vacuum diagonal.* J⁺φ = Σ_{u∈U} w_u |1@u, 2 on T⟩ and E_{j+2} = diag(t_u) + ½w_Uw_Uᵀ on these strings, so
  2(1 − ½⟨w_U, E⁻¹w_U⟩) = 2(1 − ½Φ(0)/(1 + ½Φ(0))) = 2/Γ. (If U = ∅ there is no J-term and the entry is 2.) ∎

`check_closed.py` compares 2^{k−j} times the block computed by ylib (the exact builder, all column strings) with this
closed form, entry by entry, for 51 blocks: all (k, j) with 1 ≤ j ≤ 6 and j + 1 ≤ k ≤ 7 (up to three sets T each),
and (k, j) = (8, 1), (8, 3), (8, 5), (8, 7). It also checks that the A-part of the block (2 on the diagonal, ½ between
X_i and X_l, −1 between X_i and X_vac, 2 at (vac, vac)) is ylib's A, so that the scaled block of the negative controls
is A − (1001/1000)(A − M).

## 3. The vacuum elimination

**Lemma 2.** Write M = [[C, −A], [−Aᵀ, 2/Γ]] and let F = C − (Γ/2)AAᵀ be the Schur complement of the vacuum entry.
Then M ≻ 0 if and only if F ≻ 0. For a pair i ≠ l put p = a_i, q = a_l, r = x_T − p − q ≥ 0 and

    X = 2 + Φ(0),  Y = 2 + G_i + G_l + Φ(p + q),  Z = 2 + g_i + Φ(p),  W = 2 + g_l + Φ(q),
    δ = (1 + g_i)(1 + g_l),  ε = (1 + G_i)(1 + G_l),  D = δXY − εZW,
    N = 2qZ²WY − q(1 + g_i)²XWY − (q + r)δXYZ + (q + r)εZ²W.

Then F_ii = 2 − X(1 + g_i)²/Z², F_il = −D/(ZWY), and

    H_il := a_l F_ii + (x_T − a_i) F_il = N/(Z²WY).

*Proof.* 2/Γ > 0, so M ≻ 0 iff its Schur complement F is ≻ 0. Since Γ = X/2, A_i = 2(1 + g_i)/Z, A_l = 2(1 + g_l)/W
and C_il = ε/Y, we get F_ii = 2 − (X/4)A_i² and F_il = ε/Y − δX/(ZW) = −D/(ZWY). Finally
q F_ii − (q + r) D/(ZWY) = N/(Z²WY), with a_l = q and x_T − a_i = q + r. ∎

So F_il < 0 iff D > 0, and H_il > 0 iff N > 0 (X, Y, Z, W > 0).

## 4. The measure and the two kernels

Let μ = Σ_{u∈U} x_u δ_{t_u}: a positive measure of mass e with support in [x_T, c), and Φ(d) = ∫ dμ(t)/(t + d). The
other copies enter the block only through μ. Suppose U ≠ ∅, so e > 0. Since ∫ dμ = e, every constant can be written as
an integral against μ, and

    (X, Y, Z, W) = ∫ ψ(t) dμ(t),   ψ(t) = (2/e + 1/t, (2 + G_i + G_l)/e + 1/(t + p + q), (2 + g_i)/e + 1/(t + p),
                                            (2 + g_l)/e + 1/(t + q)).

**Lemma 3.** (a) D = ½ Σ_{u,v∈U} x_u x_v Dk(t_u, t_v), with the symmetric kernel

    Dk(t, t′) = δ(X(t)Y(t′) + X(t′)Y(t)) − ε(Z(t)W(t′) + Z(t′)W(t)),

where X(t), …, W(t) are the components of ψ(t).

(b) Put α₁ = c − p, α₂ = c − q, β₁ = c − p + q, β₂ = c − q + p (so g_i = p/α₁, g_l = q/α₂, G_i = p/β₁, G_l = q/β₂) and
Λ = α₁²α₂β₁β₂. With den(t) = t(t + p)(t + q)(t + p + q) and the slot polynomials

    X̂(t) = e den(t)(2/e + 1/t),  Ẑ(t) = e den(t) α₁((2 + g_i)/e + 1/(t + p)),  Ŵ(t) = e den(t) α₂((2 + g_l)/e + 1/(t + q)),
    Ŷ(t) = e den(t) β₁β₂((2 + G_i + G_l)/e + 1/(t + p + q)),  M̂(t) = den(t),

let P(t₁, t₂, t₃, t₄) be the sum over the 24 bijections σ of {1, 2, 3, 4} onto the slots of

    2q Ẑ(t_σ1)Ẑ(t_σ2)Ŵ(t_σ3)Ŷ(t_σ4) − q c² X̂(t_σ1)Ŵ(t_σ2)Ŷ(t_σ3)M̂(t_σ4)
      − (q + r) c² X̂(t_σ1)Ŷ(t_σ2)Ẑ(t_σ3)M̂(t_σ4) + (q + r)(β₁ + p)(β₂ + q) Ẑ(t_σ1)Ẑ(t_σ2)Ŵ(t_σ3)M̂(t_σ4).

Then

    Λ N = Σ_{u1,u2,u3,u4 ∈ U} x_u1 x_u2 x_u3 x_u4 P(t_u1, t_u2, t_u3, t_u4) / (24 e⁴ den(t_u1) den(t_u2) den(t_u3) den(t_u4)).

*Proof.* (a) D is a bilinear expression in the four integrals; expand the products and symmetrise. (b) The
denominators of N are α₁ (in Z), α₂ (in W), β₁β₂ (in Y), α₁α₂ (in δ), β₁β₂ (in ε) and α₁² (in (1 + g_i)²); multiplying by
Λ and putting Z̃ = α₁Z, W̃ = α₂W, Ỹ = β₁β₂Y gives
Λ N = 2qZ̃²W̃Ỹ − qc²XW̃Ỹ − (q + r)c²XỸZ̃ + (q + r)(β₁ + p)(β₂ + q)Z̃²W̃ (using 1 + g_i = c/α₁, 1 + g_l = c/α₂,
1 + G_i = (β₁ + p)/β₁, 1 + G_l = (β₂ + q)/β₂). Each of X, Z̃, W̃, Ỹ is ∫ of a per-atom function; the three terms of
degree 3 are multiplied by 1 = ∫ dμ/e. Expanding the products of four integrals into a fourfold sum over atoms,
symmetrising over the four slots (the factor 1/24) and clearing the denominators e den(t) of the slot functions gives
the formula. ∎

Neither (a) nor (b) uses the special form of μ (mass x_u at the position c − 2x_u); they hold for every finite
measure of mass e.

## 5. Inequality (i)

**Proposition A.** For every j ≥ 1, k ≥ j + 1, all positive weights and all i ≠ l, D > 0, hence F_il < 0.

*Proof.* (A) U = ∅. Then Φ = 0, c = x_T, and D = 2δ(2 + G_i + G_l) − ε(2 + g_i)(2 + g_l) is a rational function of
(p, q, r). In exact arithmetic (`cert_pair.py`) its numerator has 6 terms, all with positive coefficients, one of them
a monomial in (p, q) alone, and its denominator is a product of four linear forms with positive coefficients. So D > 0
for p, q > 0, r ≥ 0.

(B) U ≠ ∅. By Lemma 3(a) it suffices that Dk(t, t′) > 0 for t, t′ ∈ [x_T, c). Substitute
t = x_T + 2e s₁/(1 + s₁), t′ = x_T + 2e s₂/(1 + s₂); as s₁, s₂ range over [0, ∞) this covers [x_T, c)². In exact
arithmetic Dk = N_D/Den, where Den is a product of polynomials with only positive coefficients (each containing a
monomial in (p, q, e) alone, so Den > 0), and N_D, a polynomial in (p, q, r, e, s₁, s₂) with 7 092 terms, has exactly
one negative coefficient: −8 p q e¹⁰ s₁² s₂². The coefficients of p q e¹⁰ s₁³ s₂ and p q e¹⁰ s₁ s₂³ are 8, so
N_D = N′ + 4 p q e¹⁰ s₁ s₂ (s₁ − s₂)², where N′ has only nonnegative coefficients (7 091 terms) and a positive monomial
in (p, q, e) alone. So N_D > 0 for p, q, e > 0 and r, s₁, s₂ ≥ 0, and Dk > 0 on [x_T, c)². Then
D = ½ Σ_{u,v} x_u x_v Dk(t_u, t_v) > 0. ∎

## 6. The row inequality

**Proposition B.** For every j ≥ 1, k ≥ j + 1, all positive weights and all i ≠ l,
H_il = a_l F_ii + (x_T − a_i) F_il > 0.

*Proof.* By Lemma 2 it suffices that N > 0.

(A) U = ∅. H_il is a rational function of (p, q, r) (X = 2, c = x_T). In exact arithmetic (`cert_row.py`) its numerator
has 20 terms, all with positive coefficients, including a monomial in (p, q) alone, and its denominator is a product of
positive-coefficient factors. So H_il > 0.

(B) U ≠ ∅. Substitute t_m = x_T + s_m (s_m ≥ 0, m = 1, …, 4) in P. In exact arithmetic (`cert_row.py`, about 3 s)
P is a polynomial in (p, q, r, e, s₁, s₂, s₃, s₄) with 357 055 terms, total degree 22 (degrees 20, 22, 21, 9, 4, 4, 4,
4 in the separate variables), **all coefficients nonnegative**, and a positive monomial in (p, q, e) alone. So P > 0
whenever p, q, e > 0 and r, s_m ≥ 0, that is, at every 4-tuple of positions in [x_T, ∞). The factors e den(t) and Λ
are products of positive-coefficient polynomials. By Lemma 3(b), Λ N > 0. ∎

## 7. The theorem

**Theorem DJ.** For every j ≥ 1, every k ≥ j + 1, every positive weight vector and every set T of j + 1 copies, the
type I block of (R∞)_j in which the copies of T carry the letter 2 is positive definite.

*Proof.* By Lemma 1 the block is 2^{−(k−j)} M, and by Lemma 2 it is positive definite iff F is. Fix i. By
Proposition A, |F_il| = −F_il for l ≠ i. By Proposition B, a_l F_ii > (x_T − a_i)|F_il| for every l ≠ i. Summing over
l ≠ i and using Σ_{l≠i} a_l = x_T − a_i > 0 (j ≥ 1),

    (x_T − a_i) F_ii > (x_T − a_i) Σ_{l≠i} |F_il|,   so   F_ii > Σ_{l≠i} |F_il|.

So F is a symmetric matrix with positive diagonal that is strictly diagonally dominant; by Gershgorin's theorem F ≻ 0.
Hence the block is positive definite. ∎

(F is a Z-matrix with positive row sums, so it is a nonsingular M-matrix; this is the route suggested by the earlier
numerical observations, with the row condition replaced by the pairwise inequality of Proposition B.)

**Corollary.** For every j ≥ 1, every k ≥ max(j, 2) and every positive weight vector, every block of (R∞)_j that
contains no letter 1 is positive definite: the type I ones by Theorem DJ, the type II ones by Proposition T2′ of the README.

## 8. Remarks

- *What the proof uses.* The other copies enter only through the measure μ, and only through the representation
  (X, Y, Z, W) = ∫ψ dμ: Proposition A holds for every positive measure of mass e on [x_T, c), and Proposition B for every
  positive measure of mass e on [x_T, ∞). No moment bounds, no vertices and no use of the link between masses and
  positions are needed. The joint structure of the four values Φ(0), Φ(a_i), Φ(a_l), Φ(a_i + a_l) (which the independent
  interval bounds lose) is kept exactly by expanding over pairs and 4-tuples of atoms.
- *Sharpness.* The block degenerates only where one copy of T carries all the weight: there F_ii → 0, and so does
  H_il. Numerically (`numerics.py`) F_il/(a_i a_l) ≤ −0.16 and H_il/(a_l(c − a_i)) ≥ 0.80 over 6 000 configurations
  (≥ 1.5 when U = ∅); the smallest eigenvalue of M in these samples is of order 10⁻⁴, near that corner. The negative controls sit there:
  with B → (1 + 1/1000)B the block is indefinite at s = (10⁴ or 10⁵, 1, 1.1, …, ½, 0.7, …) and Propositions A and B fail.
- *The pairwise form of the row condition.* Proposition B is the termwise version of the row-sum condition F1 > 0:
  it bounds each |F_il| by a_l F_ii/(x_T − a_i), and these bounds sum to F_ii. It holds with the margin above, so the
  sup over the free weight (which made the earlier row-sum reduction delicate) is not needed.
- *Explorations not used in the proof* (`explore/`, floating point and symbolic, no claim rests on them): the minimum of
  D over U-configurations is in all samples at the limit of many tiny U-weights, and the minimum of H_il at one U-copy
  or at that limit (`ex1.py`, `ex3.py`); moving every atom of μ to c does not increase D, by a derivative condition
  K = K₀ + eK₁ ≥ 0 whose first part has only positive coefficients (`sym1.py`) and whose second part is positive in
  floating point (`ex2.py`, `ex4.py`, `ex5.py`) but has negative coefficients in the substitutions tried (`sym1.py`,
  `sym2.py`); a single atom at any
  position gives an all-positive numerator for H_il (`sym3.py`), and so does a two-atom measure with arbitrary masses
  and positions (`sym4.py`, 40 648 terms); the 4-kernel was first checked in floating point (`ex6.py`) and then built
  exactly (`k4test.py`, identical to `cert_row.py`); the 2-kernel of D with the substitution t = x_T + s has 4 negative
  coefficients, the compact substitution one (`sym5.py`).

## 9. Files and checks

- `djlib.py`: the closed form of Lemma 1 (exact or floating point; `bscale` multiplies the B-part, everything that
  contains an E⁻¹, for negative controls), the vacuum elimination, and the ylib states of the block.
- `ylib.py`, `flintfield.py`: copies of `relative-anchors/code/ylib.py` and `flintfield.py`.
- `certlib.py`: helpers of the certificates (coefficient signs, positive monomials, positive factorisations of
  denominators, the AM–GM repair).
- `check_closed.py`: row 2 of the table.
- `cert_pair.py`: Proposition A (both cases). `cert_row.py`: Proposition B (both cases).
- `check_reduction.py`: rows 3 and 4: the formulas of Lemma 2 against the eliminated closed-form block at 12 rational
  configurations (|T| = 2, …, 5, |U| = 0, …, 4, every ordered pair), and the identities of Lemma 3 evaluated with the
  very rational function Dk of `cert_pair.py` and the very polynomial P of `cert_row.py` (built by the same functions),
  at 9 configurations with |U| = 1, …, 4 (two ordered pairs each); the symmetry of P at points.
- `check_pd.py`: row 7. `negative_control.py`: row 8. `numerics.py`: row 9 (floating point).
- `guard.sh`: a copy of `relative-anchors/code/middle_j1/guard.sh` with the memory limit set to 2 500 000 KB (kills its
  Python child above that or after a time limit).
- `independent_cert.py`, `independent_e2e.py`: row 11 (separate implementations of the certificates and of the
  closed form).
- `run_all.sh`: reruns every step, one Python process at a time under `guard.sh`, prints PASS or FAIL per step and ALL
  PASS at the end; `logs/run_all.out` and `logs/*.log` are from a complete rerun in a fresh copy of this folder.
- `explore/`: the explorations of Section 8 and their logs.

## 10. What is left

Theorem DJ settles every block without a letter 1 of the middle conditions. In the middle conditions (R∞)_j,
1 ≤ j ≤ k − 2, the letter-1 blocks and the type I blocks that contain a letter 1 remain open for k ≥ 5 (for j = 1 they
are reduced to Conjectures L and L1 in `relative-anchors/code/middle_j1/`); so do the type I blocks of (R∞)_{k−1} with
at least three copies free of the letter 2. Their entries contain pair and triple functions of the weights, so the
other copies do not enter through a single measure, and the kernel argument used here does not apply as it stands.
