# (R∞)₁, the type I blocks with |T| = 1 for general weights: two sub-blocks for every k, and a two-copy decomposition

Ansh Mishra and Aryan Senthilkumar, 5 November 2026. Nothing here is refereed.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Proposition Q (P- and Q-states, every k).** For every k ≥ 2, every copy t and every positive weight vector, the principal sub-block of the type I block T = {t} of (R∞)₁ on the P- and Q-states is ⪰ 0 with kernel exactly p ⊗ \|2@t⟩. More precisely X₂₂ = Π − Y_PQ Y_QQ⁻¹ Y_PQᵀ ⪰ c Π with c = (2 + 2γ′ + γ′²)/(2 + γ′)² ≥ ½ (γ′ as in Section 1). The proof is an exact identity of Theorem B0 type in which the copy t enters as an extra "phantom" copy (Section 2) | **proved** by hand; **exact checks**: the reduction against `j1lib.py` at 18 rational weight vectors (k = 3, …, 8), the identity as a quadratic form at 32 rational points (\|U\| = 1, …, 8) and symbolically for \|U\| = 2, 3, the conclusion at the 18 points, the sub-block from `j1lib.py` at 21 further points (k = 2, …, 8); negative controls at 4 points (`check_q.py`) |
| 2 | **Proposition RP (R- and P-states, every k).** For every k ≥ 2, t and positive weight vector, the principal sub-block on the R- and P-states is ⪰ 0 with kernel exactly p ⊗ \|2@t⟩ (Section 3) | **proved** by hand, with one exact one-variable polynomial certificate (Bernstein coefficients 2, 4, 2, 1); **exact checks** of the closed form, the bounds and the conclusion at 21 rational weight vectors (k = 2, …, 8); negative controls at 3 points (`check_rp.py`) |
| 3 | **Lemma K (the P-states see the R- and Q-states through k combinations).** Y_PQ = Y_RPᵀ K with an explicit K (diagonal on U plus the row t), so the coupling of the P-states is Y_RPᵀ [I K] and Φ = C − [I K]ᵀ (¼ Λ P₀ N⁻¹ P₀ Λ) [I K] (Section 3) | **proved** by hand; **exact check** against `j1lib.py` at 18 points (`check_q.py` (e)) |
| 4 | **Conjecture TC (two-copy decomposition).** For every k ≥ 3, t and positive weight vector, the block with the state P_t deleted is a sum of positive definite matrices each supported on the states of at most two copies (the copy t contributes R_t, a copy c ≠ t contributes R_c, P_c, Q_c). It implies the block statement (⪰ 0 with kernel exactly p ⊗ \|2@t⟩) (Section 4) | **conjecture**. **Exact certificates** (rational, every term checked positive definite by exact LDL, exact sum) at 25 rational weight vectors, k = 3, …, 10, including both degenerate limits, margins down to 6·10⁻⁶ (`cert_pairs.py`); at 3 points with B → (1 + 1/1000)B the scaled block is not positive definite and the construction fails. **Numerical**: in the cone with margin ratio ≥ 0.72 at 48 weight vectors, k = 3, …, 10; random positive definite matrices of the same shapes are outside the cone (`num_pairs.py`) |
| 5 | The block itself for k ≥ 5 and general weights; the principal sub-block on the R- and Q-states (C ≻ 0, part of Conjecture L) for every k; a proof of Conjecture TC | **open**. In the congruent form C ≻ 0 ⟺ A₁NA₁ ≻ ½HA₂Hᵀ the ratio stays near ¼ as σ → 0 but tends to 1 as σ → 1 (50-digit numerics, Section 5), so C needs an exact argument near σ → 1 |
| 6 | Simpler certificate shapes fail: Φ, Φ′ and the inverse-free form X are not H-matrices in general; the block norm comparison of X fails; the block is not a sum of positive semidefinite terms supported on two of the three state groups R, P, Q; so the sub-block results 1 and 2 cannot be glued by a fixed split of the shared diagonal blocks | **numerical** (`explore_structure.py`, Section 5); no claim rests on it |

What is new and proved: Propositions Q and RP and Lemma K, for every k and all positive weights. They are statements about
two of the three principal sub-blocks of the |T| = 1 block (necessary conditions for it), and a factorisation. They do
not prove the block itself for general weights and k ≥ 5, which stays open, as does Conjecture L of `middle_j1`. Nothing
here says anything about Conjecture 12 of `single-excitation-range/` in general, or about the NPT (Werner
distillability) problem. Conjecture TC is a conjecture; its exact checks are checks at points.

Notation as in `relative-anchors/README.md` (Sections 1, 2), `relative-anchors/code/middle_j1/NOTES.md` (Section 3)
and `relative-anchors/code/middle_j1_t1/NOTES.md` (Section 1). "Exact" means rational arithmetic (Python Fractions or
sympy); floating point enters only the searches and the surveys labelled numerical.

## 1. Setting and the closed forms used

Normalise |s| = 1, so x_c = s_c = w_c², and let σ = x_t, U = [k] ∖ {t}, N_t = 1 and N_c = 2 for c ∈ U. The block has the
states R_c = (2@t, 1@c) (c ∈ [k]), P_r = (1@r, 2@t) (r ∈ [k]) and Q_c = (vac, 2@t 1@c) (c ∈ U). In units 2^{k−1}
(Lemma T1.3 of `middle_j1/NOTES.md`, implemented in `tlib.py`):
- D_R = diag(N);
- Y_RP = ½ Λ P₀, with P₀ = I − wwᵀ, Λ = diag(τ_c/(1 − x_c)), τ_c = τ for c ∈ U and τ_t = τ − 2,
  τ = (1 + g_t)/(1 + γ/2), g_c = x_c/(2(1 − x_c)), γ = Σ_c g_c;
- Π = diag(N) − (N∘w)(N∘w)ᵀ/(2 − σ);
- Y_PQ[r, c] = τ′ (N_r w_r w_c − (2 − σ) δ_rc)/(2 − σ − 2x_c) for c ∈ U, with τ′ = 1/(1 + γ′/2),
  γ′ = Σ_{u∈U} x_u/(2 − σ − 2x_u);
- Y_QQ = 2(I + ½N′)⁻¹ with N′[c, c] = Σ_{u∈U∖c} x_u/(2 − σ − 2x_c − 2x_u), N′[c, c′] = w_c w_c′/(2 − σ − 2x_c − 2x_c′);
- Y_RQ = −(I + ½N₁)⁻¹ H (N₁ the pair matrix of Theorem B0, H as in Lemma T1.3).

C denotes the block on the R- and Q-states and Φ = C − Bᵀ diag(1/N) B the P-first reduction of Lemma P
(`middle_j1_t1/NOTES.md`), where B is the coupling of the P-states to the R- and Q-states. The kernel vector
p ⊗ |2@t⟩ is w on the P-states.

**Phantom weights.** Put y_c = 2x_c/(2 − σ) for c ∈ U and y₀ = σ/(2 − σ) for an extra index 0 (the "phantom copy",
which stands for the copy t). Then y₀ + Σ_U y_c = 1, and for c, u ∈ U

    2 − σ − 2x_c = (2 − σ)(1 − y_c),      2 − σ − 2x_c − 2x_u = (2 − σ) δ̃_cu,   δ̃_cu = 1 − y_c − y_u ≥ y₀ > 0.

Write λ̃_c = 1/(1 − y_c), g̃_c = y_c/(2(1 − y_c)) = (λ̃_c − 1)/2; then γ′ = Σ_{c∈U} g̃_c, τ′ = 2/(2 + γ′),
N′[c, c] = Σ_{u∈U∖c} y_u/(2δ̃_cu) and N′[c, c′] = w̃_c w̃_c′/(2δ̃_cc′) with w̃_c = (2/(2 − σ))^{1/2} w_c (so w̃_c² = y_c).
So the Q-states carry the column block of Theorem B0 for the copies U with the weights y, except that the total weight
of U is 1 − y₀ < 1: the phantom copy 0 takes the rest.

## 2. Proposition Q: the P- and Q-states

**Proposition Q.** For every k ≥ 2, every t and every positive weight vector,

    X₂₂ := Π − Y_PQ Y_QQ⁻¹ Y_PQᵀ ⪰ c Π,     c = (2 + 2γ′ + γ′²)/(2 + γ′)² ≥ ½.

Hence the principal sub-block of the |T| = 1 block on the P- and Q-states is ⪰ 0 with kernel exactly p ⊗ |2@t⟩.

*Step 1: the reduction.* For π ∈ R^k put π̃_r = π_r/w_r, m = Σ_r y_r π̃_r (with y_t := y₀) and ζ_r = π̃_r − m, so that
Σ_{a∈{0}∪U} y_a ζ_a = 0 (ζ₀ := ζ_t). Then

    πᵀ Π π = (2 − σ) S(ζ),      πᵀ Y_PQ Y_QQ⁻¹ Y_PQᵀ π = (2 − σ) L(ζ)/(2 + γ′)²,

where S(ζ) = Σ_{a∈{0}∪U} y_a ζ_a² and

    L(ζ) = Σ_{c∈U} y_c λ̃_c² ζ_c² + ¼ Σ_{c<u∈U} (y_c y_u/δ̃_cu) (λ̃_c ζ_c + λ̃_u ζ_u)².

Proof. Since N_r x_r/(2 − σ) = y_r, πᵀΠπ = (2 − σ)(Σ y_r π̃_r² − (Σ y_r π̃_r)²) = (2 − σ) S(ζ). Next,
(Y_PQᵀ π)_c = τ′(w_c Σ_r N_r w_r π_r − (2 − σ) π_c)/(2 − σ − 2x_c) = −τ′ λ̃_c w_c ζ_c, because Σ_r N_r w_r π_r = (2 − σ) m.
With Y_QQ⁻¹ = ½(I + ½N′) and x_c = (2 − σ) y_c/2 the Gram form is
½τ′²(Σ_U x_c λ̃_c² ζ_c² + ½ Σ_{c<u} x_c x_u (λ̃_c ζ_c + λ̃_u ζ_u)²/((2 − σ) δ̃_cu)) = (2 − σ)(τ′²/4) L(ζ). ∎

So X₂₂ ⪰ cΠ is the inequality L ≤ (1 − c)(2 + γ′)² S on the constraint subspace Σ y_a ζ_a = 0. This is inequality (C0)
of Theorem B0 for the copies {0} ∪ U, except that the phantom copy 0 appears only in the constraint and in S, not in L
and not in γ′.

*Step 2: the identity.* On the subspace Σ_{a∈{0}∪U} y_a ζ_a = 0, with d_ab = ζ_a − ζ_b,

    (2 + 2γ′) S − L = R₁ + R₂ + R₃ + R₄,
    R₁ = Σ_{r<s∈U} 2 y_r y_s (γ′ − g̃_r − g̃_s) d_rs² + Σ_{r∈U} y_r y₀ (1 + 2γ′ − 2g̃_r) d_r0²,
    R₂ = ½ Σ_{l∈U} y_l λ̃_l² Σ_{a∈U∖l} y_a y₀ d_a0²,
    R₃ = ¼ Σ_{r<s∈U} Σ_{u∉{r,s}} y_r y_s y_u (λ̃_s d_ru − λ̃_r d_su)²,
    R₄ = ¼ Σ_{r<s∈U} (y_r y_s/δ̃_rs) Σ_{u<v; u,v∉{r,s}} y_u y_v (λ̃_r + λ̃_s)² d_uv²,

where u and v range over {0} ∪ U (the phantom copy included). Every coefficient is ≥ 0: γ′ − g̃_r − g̃_s is a sum of
the other g̃_u, and 1 + 2γ′ − 2g̃_r > 0.

Proof (the three steps of Theorem B0, with the phantom copy kept in the constraint).
- *Variance identity.* For every r, the constraint and Σ y = 1 give ζ_r = Σ_{a≠r} y_a d_ra, so λ̃_r ζ_r is the mean of
  the d_ra with weights y_a/(1 − y_r), and
  λ̃_r Σ_{a≠r} y_a d_ra² − λ̃_r² ζ_r² = λ̃_r² Σ_{a<b; a,b≠r} y_a y_b d_ab².
  Multiply by y_r and sum over r ∈ U (not over the phantom copy): E − D = G, with D the first sum of L,
  E = Σ_{r∈U} y_r λ̃_r Σ_{a≠r} y_a d_ra² and G = Σ_{r∈U} y_r λ̃_r² Σ_{a<b; a,b≠r} y_a y_b d_ab².
- *Cauchy–Schwarz.* For r ≠ s in U, using λ̃ y = λ̃ − 1, λ̃_r ζ_r + λ̃_s ζ_s = Σ_{u∉{r,s}} y_u (λ̃_s d_ru + λ̃_r d_su), a
  sum with total weight δ̃_rs = Σ_{u∉{r,s}} y_u. Lagrange's identity and 2a² + 2b² − (a + b)² = (a − b)² give
  ¼(y_r y_s/δ̃_rs)(λ̃_r ζ_r + λ̃_s ζ_s)² = ½ y_r y_s Σ_u y_u (λ̃_s² d_ru² + λ̃_r² d_su²) − (the R₃ and R₄ terms of the pair
  r, s), since the differences of the summands are (λ̃_r + λ̃_s) d_vu. Summed over r < s, the first part counts every
  term y_l λ̃_l² y_a y_b d_ab² of G once when a, b ∈ U and one half when one of a, b is the phantom copy; so the second
  sum T of L equals G − R₂ − R₃ − R₄.
- *The coefficient.* E = Σ_{r<s∈U} y_r y_s (λ̃_r + λ̃_s) d_rs² + Σ_{r∈U} y_r y₀ λ̃_r d_r0² and, by the constraint,
  S = Σ_{a<b} y_a y_b d_ab² over all pairs of {0} ∪ U. With λ̃ = 1 + 2g̃, (2 + 2γ′)S − E = R₁.

So (2 + 2γ′)S − L = (2 + 2γ′)S − D − T = R₁ + (E − D − G) + R₂ + R₃ + R₄ = R₁ + R₂ + R₃ + R₄. ∎

*Step 3: conclusion.* By Steps 1 and 2, πᵀX₂₂π = (2 − σ)(S − L/(2 + γ′)²) ≥ (2 − σ) S (1 − (2 + 2γ′)/(2 + γ′)²) =
c πᵀΠπ, and c ≥ ½ because (2 + γ′)² − 2(2 + 2γ′) = γ′² ≥ 0. If πᵀX₂₂π = 0 then S(ζ) = 0, so ζ = 0, that is
π ∈ span(w). Since Y_QQ ≻ 0
and Y_PQᵀ w = 0, the sub-block on the P- and Q-states is ⪰ 0 and its kernel is {(π, −Y_QQ⁻¹Y_PQᵀπ) : π ∈ ker X₂₂} =
span((w, 0)). ∎

Remarks.
- For k = 2 (U = {c}) the identity reads (2 + 2γ′)S − L = y_c y₀ d_c0².
- Numerically the best constant is about 5/9 (smallest ratio 0.557 in 1 500 samples, `explore_structure.py`); the
  proved constant c tends to ½ only as γ′ → 0, that is when the copy t carries almost all the weight.
- Negative controls (`check_q.py` (f)): with the B-part scaled by 1001/1000 (`j1lib.py`, bscale) the sub-block is
  negative on (w, 0) at every point (−4·10⁻³ to −1.25·10⁻² at the four points), and with P_t deleted it is not positive
  definite at four explicit rational points with a small weight on one copy, while the true sub-block is.

## 3. Proposition RP and Lemma K

**Proposition RP.** For every k ≥ 2, every t and every positive weight vector, the principal sub-block of the
|T| = 1 block on the R- and P-states is ⪰ 0 with kernel exactly p ⊗ |2@t⟩.

*Proof.* D_R ≻ 0 and Y_RP w = 0, so it suffices that Π_R = Π − Y_RPᵀ D_R⁻¹ Y_RP ⪰ 0 with kernel span(w).
1. *Closed form.* For π ∈ R^k put e_r = π_r/w_r − Σ_s w_s π_s; then Σ_r x_r e_r = 0, that is
   Σ_{c∈U} x_c e_c = −σ e_t, and (P₀π)_r = w_r e_r. Since Πw = 0 and Σ_r N_r x_r e_r = −σ e_t,
   πᵀΠπ = 2Σ_U x_c e_c² + 2σ(1 − σ) e_t²/(2 − σ), while πᵀY_RPᵀD_R⁻¹Y_RPπ = ¼Σ_r Λ_r² x_r e_r²/N_r. So

       πᵀ Π_R π = Σ_{c∈U} a_c x_c e_c² + b_t e_t²,   a_c = 2 − τ² λ_c²/8,
       b_t = 2σ(1 − σ)/(2 − σ) − (2 − τ)² σ/(4(1 − σ)²),   λ_c = 1/(1 − x_c).

2. *Two bounds.* Since γ ≥ g_t + g_c, τλ_c ≤ φ(g_c) := (1 + g_t)(1 + 2g_c)/(1 + (g_t + g_c)/2), which increases in
   g_c; and g_c ≤ x_c/(2σ) ≤ (1 − σ)/(2σ) because 1 − x_c ≥ σ. With g_t = σ/(2(1 − σ)) this gives
   τλ_c ≤ 2(2 − σ)/(1 + 2σ − 2σ²), hence a_c ≥ A(σ) := 2 − (2 − σ)²/(2(1 + 2σ − 2σ²)²) > 0 for 0 < σ < 1. Next
   2 − τ = 2(1 + G)/(2 + g_t + G) with G = Σ_U g_c ≤ (1 − σ)/(2σ), increasing in G, so
   2 − τ ≤ 2(1 − σ²)/(1 + 2σ − 2σ²) and b_t ≥ −σ br(σ), br(σ) = (1 + σ)²/(1 + 2σ − 2σ²)² − 2(1 − σ)/(2 − σ).
3. *Cauchy–Schwarz.* σ² e_t² = (Σ_U x_c e_c)² ≤ (Σ_U x_c/a_c)(Σ_U a_c x_c e_c²) ≤ ((1 − σ)/A) Σ_U a_c x_c e_c².
4. *One variable.* F(σ) = σA(σ) − (1 − σ) br(σ) = 3σ p(σ)/(2(2 − σ)(2σ² − 2σ − 1)²) with
   p(σ) = 5σ³ − 12σ² + 6σ + 2, whose Bernstein coefficients on [0, 1] are 2, 4, 2, 1; so F > 0 on (0, 1).
5. *Conclusion.* If br ≤ 0 then b_t ≥ 0 and πᵀΠ_Rπ ≥ Σ_U a_c x_c e_c². If br > 0, choose θ < 1 with
   θσA/(1 − σ) ≥ br (possible by step 4); steps 2 and 3 give
   πᵀΠ_Rπ ≥ (1 − θ)Σ_U a_c x_c e_c² + σ e_t² (θσA/(1 − σ) − br) ≥ (1 − θ)Σ_U a_c x_c e_c². In both cases πᵀΠ_Rπ = 0
   forces e_c = 0 on U, then e_t = 0 by the constraint, so π ∈ span(w). ∎

`check_rp.py` checks the closed form of step 1 (as quadratic forms, by polarisation), the bounds of step 2, the final
inequality and the positivity of the sub-block with P_t deleted exactly at 21 rational weight vectors (k = 2, …, 8),
and step 4 symbolically. Negative control: with the B-part scaled by 1001/1000 the sub-block with P_t deleted is not
positive definite at three rational points next to the corner σ → 0 with one dominant copy of U (where a_c → 0).

**Lemma K.** For c ∈ U put

    K[c, c] = −(2τ′/τ)(2 − σ)(1 − x_c)/(2 − σ − 2x_c),   K[t, c] = 2 w_t w_c τ′ (1 − σ)/((2 − τ)(2 − σ − 2x_c)),

and K[r, c] = 0 for the other r. Then Y_PQ = Y_RPᵀ K. Hence the coupling of the P-states to the R- and Q-states is
B = Y_RPᵀ [I K]: the P-states see (ρ, q) only through the k combinations ρ̂ = ρ + Kq, and

    Φ = C − [I K]ᵀ Ω [I K],   Ω = Y_RP N⁻¹ Y_RPᵀ = ¼ Λ P₀ N⁻¹ P₀ Λ   (a diagonal matrix plus rank two, no pair sums).

*Proof.* (Y_RPᵀK)[r, c] = Y_RP[c, r] K[c, c] + Y_RP[t, r] K[t, c]. Here ½Λ_c K[c, c] = −τ′(2 − σ)/(2 − σ − 2x_c) and
½Λ_t K[t, c] = −w_t w_c τ′/(2 − σ − 2x_c), so the sum is
−τ′((2 − σ)(δ_rc − w_r w_c) + w_t w_c(δ_rt − w_r w_t))/(2 − σ − 2x_c) = τ′(N_r w_r w_c − (2 − σ)δ_rc)/(2 − σ − 2x_c),
using σ = w_t² and N_r = 2 − δ_rt. ∎

So, in the coordinates (π, ρ̂, q), the P-states couple only to ρ̂, through ½P₀Λ, exactly as the rows of Theorem B0
couple to its columns; the block is ⪰ 0 with kernel exactly p ⊗ |2@t⟩ if and only if the q-block Ĉ_QQ of C in the
coordinates (ρ̂, q) is positive definite and S ≻ Ω, where S is the Schur complement of C onto ρ̂ (q eliminated at
fixed ρ̂). Unlike Theorem B0, S⁻¹ = [I K] C⁻¹ [I K]ᵀ is not explicit (numerically it is not 2(I + ½N′) plus a matrix
of small rank), so this does not reduce the block to an inequality of type (C0). Lemma K is checked exactly against
`j1lib.py` at 18 rational weight vectors (`check_q.py` (e)); in Proposition Q it appears as (Y_PQᵀπ)_c = −τ′λ̃_c w_c ζ_c.

## 4. Conjecture TC: a decomposition over pairs of copies

Let G be the |T| = 1 block with the state P_t deleted (the kernel vector has P_t-coordinate w_t ≠ 0, so the block is
⪰ 0 with kernel exactly p ⊗ |2@t⟩ if and only if G ≻ 0). Group the states of G by copies: G_t = {R_t} and
G_c = {R_c, P_c, Q_c} for c ∈ U.

**Conjecture TC.** For every k ≥ 3, every t and every positive weight vector,

    G = Σ_{a<b} Z_ab + Σ_a S_a,

with every Z_ab positive definite and supported on the states of G_a ∪ G_b, and every S_a positive definite and
supported on G_a.

Conjecture TC implies G ≻ 0, hence the block statement for that k. The property does not depend on the basis chosen
inside each group (it is invariant under congruence by block-diagonal matrices with one block per copy, in particular
under the rescaling of `ylib`), but it does depend on the grouping by copies. It is a strong restriction: a positive
definite matrix in general position is not of this form.

Evidence.
- **Exact certificates** (`cert_pairs.py`): a floating-point semidefinite program finds the terms (maximising the
  margin eps with G − eps diag(G) in the cone); the off-diagonal blocks of every Z_ab are then set to the exact blocks
  of G, the diagonal blocks are rounded to rationals after a shift by a quarter of the margin, the S_a are computed
  exactly as the remainders, and every Z_ab and S_a is checked positive definite by exact LDL with the exact sum
  equal to G. This succeeds at 25 rational weight vectors, k = 3, …, 10 (random, spread over powers of 3, one small
  weight, one dominant weight, the shape of the near-failure below, and both degenerate limits of
  `middle_j1_t1/NOTES.md`, Section 5: small weight on t with two dominant copies of U, and t carrying almost all the
  weight). The smallest SDP margin is 6.1·10⁻⁶ (k = 4, w = (19/10, 1/100, 6/5, 8/5), t = 1).
- **Negative controls**: with the B-part scaled by 1001/1000, at three rational points the scaled G is not positive
  definite (exact LDL), so it is outside the cone; the SDP margin is negative (about −2·10⁻³) and the construction fails.
- **Numerical** (`num_pairs.py`, double precision): 48 weight vectors, k = 3, …, 10, including both limits: the margin is
  positive everywhere, and its ratio to the smallest eigenvalue of the normalised G is 1.000 (to solver accuracy) except
  when t carries most of the weight, where it is between 0.72 and 0.99. Random positive definite matrices of the same
  sizes (Wishart plus 0.02 I, normalised) are outside the cone (9 of 9, margins between −0.4 and −1.4). The solver
  (SCS) reports reduced accuracy at a few of the degenerate points; the survey is indicative only. In a first run
  without the floor 10⁻⁶ on the weights the margin was −5·10⁻⁸ at a point with x_t ≈ 10⁻⁷, where the smallest
  eigenvalue is 6·10⁻⁸, below the solver accuracy; at the same shape with x_t = 10⁻⁴, …, 10⁻⁶ the ratio is 1.000
  (`logs/dev_num_pairs.log`), and a rational point of that shape with x_t ≈ 10⁻⁴ has an exact certificate.
- **Numerical, other forms** (`explore_structure.py`, item 6): the undeleted block, with G_t = {R_t, P_t}, lies in the
  closed two-copy cone to solver accuracy (it is singular, so the margin is 0). The inverse-free congruent form
  diag(A₁, I, A₂)ᵀ M diag(A₁, I, A₂) (A₁ = I + ½N₁, A₂ = I + ½N′, entries polynomial in the pair kernels) does not: its
  margin is clearly negative at 5 of 12 points (down to −4.4·10⁻²; the undeleted block has |margin| ≤ 2.3·10⁻⁹ at
  the same points). The congruence by A₁ and A₂ mixes the copies, and the property is lost.

What a proof would need: an explicit choice of the terms. For a pair (c, u) of copies of U the off-diagonal block of
Z_cu is fixed (it is the block of G between the two copies); only the split of the diagonal blocks is free. If the split
can be chosen as an explicit function of (x_t, x_c, x_u) and finitely many symmetric sums of the weights, positivity of
each term is a finite-dimensional statement, and the measure-kernel method of Theorem DJ could apply to it. We did not
find such a split. The entries of G contain (I + ½N₁)⁻¹ and (I + ½N′)⁻¹, and the inverse-free form, whose entries are
explicit, is not in the cone (previous item), so a split has to be found for G itself (or for another form obtained
from G by a congruence that does not mix the copies).

## 5. Numerical remarks (no claim rests on them)

From `explore_structure.py` (double precision, 1 500 random and degenerate weight vectors, k = 3, …, 8):
- the smallest ratio X₂₂ : Π on w^⊥ is 0.557 (Proposition Q proves ≥ ½);
- the comparison matrix (H-matrix test) is not positive definite for Φ at 245, for Φ′ at 1 488 and for X with P_t
  deleted (rescaled) at 402 of the 1 500 samples; the block norm comparison of X with the 2 × 2 blocks (η_c, π_c) fails
  at 653; so certificates of diagonal-dominance type, in these bases, cannot work for general weights;
- the gauge-fixed block is not a sum of positive semidefinite terms supported on two of the three state groups R, P, Q
  (SDP margin negative at 7 of 13 points), so Propositions Q and RP cannot be combined by a fixed split of the shared
  diagonal blocks; the grouping by copies (Conjecture TC) is the one that works;
- the block C on the R- and Q-states (the remaining two-group sub-block; Conjecture L of `middle_j1` asks C ≻ 0) has
  normalised smallest eigenvalue down to 6·10⁻⁵ at σ ≈ 0.9998 in these samples. A better measure is the congruent
  form: C ≻ 0 if and only if A₁NA₁ ≻ ½HA₂Hᵀ (congruence by diag(A₁, A₂), then the Q-states eliminated). In 50-digit
  arithmetic (`explore_rq.py`) the largest generalised eigenvalue r of ½HA₂Hᵀ against A₁NA₁ is 0.2500003 at k = 3,
  σ = 10⁻⁶, so the soft direction of C as σ → 0 is harmless there, but r reaches 0.99996 at σ = 0.99995 (over 300
  rational weight vectors, k = 3, …, 7). So C needs an exact argument only near σ → 1. (In double precision r came
  out above 1 at σ = 10⁻⁶, k = 3: a cancellation artefact.)

A structural remark, not used: in the limit w_t → 0 the coupling of the R_U-states to the P-states is −½ times that
of the Q-states, and the substitution q ↦ q − ½ρ_U decouples the block into the letter-1 block of Theorem B0 on the
copies U, the block 2I − ½(I + ½N′)⁻¹ on the R_U-states and a 2 × 2 block on (R_t, P_t). This is the precise sense in
which "the spectator disappears" as σ → 0. The degenerate limit σ → 0 with two dominant copies of U is the limit in
which Theorem B0 itself degenerates (δ → 0 for those two copies).

## 6. Files

- `j1lib.py`, `ylib.py`, `flintfield.py`, `tlib.py`, `check_j1_ylib.py`: copies from `relative-anchors/code/middle_j1_t1/`
  (unchanged): the Lemma B builder, the exact builder of `relative-anchors/code`, the rational-function field, the
  closed forms of Lemma T1.3 and Φ, and the check that `j1lib.py` equals `ylib.Builder` on six |T| = 1 blocks.
- `nlib.py`: the same closed forms in floating point (numpy), for the searches and surveys.
- `qlib.py`: the phantom weights and the forms S, L, R₁, …, R₄ of Proposition Q (exact).
- `check_q.py`: Proposition Q (reduction, identity numerically exact and symbolic, conclusion, sub-block from
  `j1lib.py`, negative controls) and Lemma K.
- `check_rp.py`: Proposition RP.
- `cert_pairs.py`: the exact certificates and negative controls of Conjecture TC.
- `num_pairs.py`, `explore_structure.py`, `explore_rq.py`: floating-point surveys (Sections 4 and 5; `explore_rq.py` in
  50-digit arithmetic).
- `guard.sh` (copy of `relative-anchors/code/no_letter1/guard.sh`): runs one Python process, killed above 2 500 000 KB
  or after a time limit. `run_all.sh`: reruns every step above, one Python process at a time, PASS or FAIL per step and
  ALL PASS at the end; its output is `logs/run_all.out`, the step logs are `logs/run_*.log` (both from a complete
  rerun in a fresh copy of this folder). `logs/dev_*.log` are logs of
  development runs of the same scripts quoted in the text.
