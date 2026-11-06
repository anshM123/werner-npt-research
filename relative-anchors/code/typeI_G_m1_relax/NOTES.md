# Conjecture M (one copy carrying the letter 2): a pair-function-free lower bound for G, the reduced finite system, Conjecture R at equal U-weights, and its local condition for all weights

Ansh Mishra and Aryan Senthilkumar, 6 November 2026. Nothing here is refereed.

Notation as in `relative-anchors/README.md` (Lemmas B, N; Theorems T, C, U3₁, EU, TL),
`relative-anchors/code/typeI_u2/NOTES.md` (Part B: the form F, the splitting X = X′ + a wᵀ + w cᵀ, formulas (F1)–(F3),
Lemma H) and `relative-anchors/code/typeI_G_m1/NOTES.md`. Throughout m = |T| = 1, T = {t}, σ = s_t, w_t = √σ,
U the other u copies, w = w_U, S = Σ_{h∈U} s_h, K = σ + S/2, κ = (2σ + S)/(σ + S). Forms are in units of 2Y. G is the
Schur complement on the letter-2 states (ρ, ζ) of the type I block in the gauge x_t = 0 (the P_{t,h} deleted);
G_z = Dg G Dg with ζ_h = w_t w_h z_h.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Lemma R1.** For every u ≥ 2 and all positive weights: after the exact elimination of q, the hard part X′ of Theorem C's splitting can be relaxed entrywise by Lemma H, and Φ_q(X′ + a wᵀ + w cᵀ, ρ, ζ) ≥ Ψ(a, c, ρ, ζ) for all hard X′, a ⊥ w, c, ρ, ζ, with Ψ the explicit form of Section 2. Hence G(ρ, ζ) ≥ inf_{a ⊥ w, c} Ψ(a, c, ρ, ζ); if Ψ ≻ 0 then G ⪰ G_rel ≻ 0 | **proved** (by hand, from the closed form of the block, (F1)–(F3) and Lemma H of Theorem C); every identity used is rechecked in exact arithmetic (row 3) |
| 2 | **Corollary R1.** At every weight vector where Ψ ≻ 0, G ≻ 0, so the type I block of (R∞)_{k−1} with T = {t} is ⪰ 0 with kernel exactly span{p ⊗ f_h : h ∈ U} | **proved** (Lemma R1 with Theorem C) |
| 3 | The identities behind Lemma R1, against the closed form `cf.TypeI` of the whole block: (I1) the q-elimination formula with the explicit coupling ℋ; (I2), (I3) = (F2), (F3); (I4) Lemma H on the hard space; (I5) the inequality Φ_q ≥ Ψ at random points; (I6) G − G_rel ⪰ 0 | **exact check**, 15 rational configurations, u = 2…6 (`check_identity.py`, PASS) |
| 4 | Ψ ≻ 0 (all leading principal minors positive, exact integer arithmetic) at 192 rational weight vectors with u = 3…10 (24 per u, eight families), and at 8 (u = 12), 8 (u = 16) and 4 (u = 20) further ones. By Corollary R1, G ≻ 0 at each of these 212 weight vectors. The runs for u = 24 and u = 32 were stopped by hand for time and are not claimed | **exact check** (pointwise; not a proof for any u; `certify_points.py`, Section 4) |
| 5 | The relaxation is almost lossless: in the frame (ρ√(1 + σ/S), merged z√(1 + S/σ), z ⊥ w), the smallest value of λ_min(Ĝ_rel) over the samples is within 3·10⁻⁴ of that of λ_min(Ĝ) (u ≤ 8, 24 to 40 weight vectors per u, σ and the U-weights over e^{±8}, some U-weights over e^{±10}) | **numerical** (`relax1.py`, `relax2.py`) |
| 6 | In that frame λ_min(Ĝ) ≥ 0.38 for every sample (u ≤ 12), the soft direction being ρ as S/σ → 0 (eigenvalue ≈ S/σ before the ρ-rescaling) and the merged direction as S/σ → ∞ | **numerical** (`ex_frame.py`, `ex_soft_dirs.py`) |
| 7 | **Conjecture R.** For every u ≥ 3 and all positive weights, Ψ ≻ 0 on {a ⊥ w} × R^U × R × R^U. It implies G ≻ 0 for m = 1 and every u, i.e. every type I block of (R∞)_{k−1} with exactly one copy carrying the letter 2, for every k | **open**; numerically true with the margins of row 5 |
| 8 | **The reduced finite system** (Section 7): Ψ = Σ_h v_hᵀ D(s_h) v_h + yᵀ C y with v_h = (ah_h, ch_h, zh_h), y = (ρ, F_1, …, F_30), F_j = Σ_h φ_j(s_h) x_h, the φ_j taken from ten fixed functions of s per variable class; D(s) is 3 × 3 (73 terms), C is 31 × 31 (4952 terms), all entries polynomials in IG = 1/(2 + g), IGT = 1/(2 + g̃), u and ten moment sums M_k = Σ_h f_k(s_h), with coefficients rational in (s, σ, S) | **derived** by a symbolic program (`reduce_cas.py`); **exact check**: T ᵀ(D + U C Uᵀ)T = Ψ of `psi_exact.py` at 24 random rational weight vectors, u = 3…8 (`check_reduced.py`, two seeds, PASS) |
| 8a | **Lemma F.** If D(s_h) ≻ 0 for every h and B_F = Σ_h U_hᵀ D(s_h)⁻¹ U_h ≻ 0, then Ψ ≻ 0 on {a ⊥ w} iff N = C + (0 ⊕ B_F⁻¹), with the row and column of the constraint functional Σ_h s_h ah_h deleted, is ≻ 0 (a fixed 30 × 30 condition). For u ≤ 9 (rank B_F ≤ 3u < 30) and for repeated U-weights B_F is singular and the pseudo-inverse version applies | **proved** (linear algebra, Section 7) |
| 8b | Conditions (i) D_h ≻ 0 and (ii) N ≻ 0 at rational weight vectors with u = 10, 11, 12: (i) at 16/16, B_F ≻ 0 at 15/16, (ii) at all 15 where it applies; the verdict agrees with the direct certificate of Ψ at all 15 | **exact check** (`certify_reduced.py ... both`) |
| 8c | **Theorem R_EU (Conjecture R at equal U-weights).** If all U-weights are equal, then Ψ ≻ 0 on {a ⊥ w} for every real u ≥ 3 and all σ, s > 0. Hence G ≻ 0 there (a second proof of the conclusion G ≻ 0 of Theorem EU, through Lemma R1) | **proved**, computer-assisted, symbolic in u (Section 8): the six leading minors of D(s) and of the invariant 3 × 3 block, with u = 3 + x, are ratios of polynomials in (σ, s, x) with only positive coefficients (`equal_R.py 3`, PASS) |
| 8d | Negative controls: with the couplings between (a, c) and (ρ, ζ) multiplied by 1.001, Ψ is indefinite at explicit points with S ≫ σ; with the coupling m to Q multiplied by 1.001, at explicit points with S ≪ σ; without the projection Π, at points with S ≳ σ; condition (ii) with B_F replaced by 1.001 B_F fails at explicit points (Section 9) | **exact** (`controls_R.py`, `certify_reduced.py ... controlB`) |
| 9 | Route 1 on the relaxed system: in the scaled coordinates (a_h/w_h, c_h/w_h, ρ, z_h) many entries of Ψ change sign with the weights, so no sign change independent of the weights makes it a Z-matrix; without the constraint a ⊥ w it is indefinite | **numerical** (`signs1.py`, u = 4, 200 weight vectors) |
| 10 | Choice of the representative of the hard-part couplings: with both couplings projected by Π = I − wwᵀ/S (the form Ψ) the bound is within 10⁻⁴; with Theorem C's unprojected 𝒢 the bound loses about 30 % at u = 8; with the unprojected full coupling it is indefinite | **numerical** (`relax2.py`) |
| 8e | The local block D(s_h) depends on the other copies only through the four bounded sums g = Σ s_k/d_k, N = Σ s_k/d_k², g0 = Σ s_k/d0_k, N0 = Σ s_k/d0_k² (u cancels identically) | **derived** (`D_natural.py`, symbolic); **exact check** against the reduced system, 132 blocks (`check_D.py`) |
| 8f | **Theorem D (condition (i) for all weights).** For every u ≥ 1, all positive weights and every h ∈ U, D(s_h) ≻ 0 | **proved**, computer-assisted, exact (Section 11): monotonicity in N, N0 by one-signed coefficients; the three leading minors of the extreme block on the region (t, ξ, θ, θ0) ∈ (0,1) × (0,1] × [0,1]² certified factor by factor (positive coefficients or tensor Bernstein coefficients, one corner rule) (`prove_D.py`, PASS); negative control: a 1/10 diagonal shift breaks it at an explicit point |
| 8g | Measurements for condition (ii): it has no slack (Section 9), so per-copy Cauchy–Schwarz alone cannot prove it; all degeneracy sits in the global 3 × 3 block on (γ, ξ, ρ); on the perpendicular subspace the global part takes at most 31 % of the local part (u = 3; 15 % at u = 4, less for larger u); the plan "Theorem D + per-copy bound on the perpendicular part + exact 3 × 3 condition" holds at all 32 weight vectors tested; C has rank 16 with 7 negative eigenvalues and a moment-dependent kernel | **numerical** (Section 12; `ex_globperp.py`, `ex_perpeps.py`, `ex_plan.py`, `ex_inertiaC.py`, `ex_kerC.py`) |
| 8h | **Lemmas P1, P2.** (A_κ) Ψ_pp ⪰ κD_pp on p and (B_κ) Ψ_gg − Ψ_gp(κD_pp)⁻¹Ψ_gpᵀ ≻ 0 imply Ψ ≻ 0; the coupling term equals (1/κ)(B0 − B1ᵀB2⁻¹B1) with three sums over copies of explicit 3 × 3 matrices | **proved** (linear algebra, Section 13); identity of P2 **exact check** at 16 weight vectors (`check_Blambda.py`) |
| 8i | (A_1/2) and (B_1/2) at 56 rational weight vectors, u = 3, 4, 5, 6, 8, 10, 12 (so Conjecture R there, by P1); control: (A_99/100) fails at 14 of them | **exact check** (`check_plan.py`, PASS) |
| 11 | (Z) and (D) of Conjecture M themselves, and Conjecture R, for general weights and every u ≥ 4; in particular condition (ii) of Lemma F for all weights, (A_κ) and (B_κ) of Section 13 for all weights, and Conjecture R for U-weights with two or three distinct values (symbolic in the weights and level sizes) | **open** (Sections 10, 12 and 13 explain what blocked them) |

| 12 | The exact checks rerun separately with other seeds: `check_identity.py` (seed 777; I1–I6 at every trial), `certify_points.py` (u = 3…8, 12 weight vectors each, seed 5: Ψ ≻ 0 at all 72), `check_reduced.py` (seeds 99 and 11), `equal_R.py 3`, `controls_R.py`, `certify_reduced.py` in mode `both` (u = 10, 11, 12, seed 3); Theorem D from scratch: `reduce_cas.py` and `D_natural.py` rerun (the same `D_natural.pkl`), `prove_D.py`, `check_D.py` (seed 23), `check_D_homog.py` (seed 4); Lemmas P1 and P2: `check_plan.py` (u = 3…6, seed 29) and `check_Blambda.py` (seed 43) | **exact**, all as stated (`logs_rerun_*.log`, `logs_rerun2_*.log`, `logs_rerun3_*.log`) |
Nothing here says anything about Conjecture 12 or the NPT problem in general. Conjecture M and Conjecture R are open.

## 1. What was measured first

- `ex_conjM.py` (floating point) reproduces Conjecture M at random weights for u = 3, 4, 6, 8: all off-diagonal entries
  of G_z are ≤ 0 and the smallest normalised (D) margin is 4.6·10⁻³, 1.5·10⁻³, 9.8·10⁻⁴, 7.9·10⁻⁴.
- `ex_soft_dirs.py`: in the frame (ρ, merged z·√(1 + S/σ), z ⊥ w) of `typeI_G_u3/`, the smallest eigenvalue of Ĝ is
  ≈ S/σ as S/σ → 0, with eigenvector e_ρ, and stays near 0.49 as S/σ → ∞. So G has one soft direction in each
  degenerate regime: ρ when the U-copies are light, the merged direction ζ ∝ w when they are heavy.
- `ex_frame.py`: rescaling ρ by √(1 + σ/S) as well, λ_min(Ĝ) ≥ 0.382 over 60 weight vectors for each u = 3, 4, 6, 8 and
  20 for u = 12 (families: spread, wide spread, one U-weight varied over e^{±10}, geometric; σ over e^{±8}). The
  perpendicular block is near 0.72–1, as in `typeI_G_u3/NOTES.md` Section 4.

So in this frame the problem is well conditioned; what is missing is a way to evaluate G without the inverse of the
(C_u) block.

## 2. Lemma R1: relaxing the hard part after eliminating q

**The form.** In the gauge x_t = 0 the block is (`typeI_G_m1/equal_u.py`, docstring; `cf.py`)

    Φ = ½‖X‖² + Σ_h X_hh² + ½(tr X + ρ)² + (3/2)ρ² + |ζ|² + 2q[−tr X + wᵀXw/(2K) − (2Kρ + w_t w·ζ)/(σ + S)] + 2q²
        − 𝔅_U(X) − 2 Σ_{r∈U} ⟨𝒯_{(t,2)}ψ_r, W 𝒯_{(r,1)} y⟩,

with y = ρ e_t + ζ on U ∪ {t}, W = E_∅⁻¹ (E_∅ = diag(d⁰) + ½w̃w̃ᵀ on U ∪ {t}, d⁰_x = s_x + (σ + S)/2, w̃ = (w_t, w)),
𝒯_{(t,2)}ψ = ½(w·ψ)e_t − ½w_tψ and 𝒯_{(r,1)}y = c_r e_r + ½w_r y, c_r = ½Y0 + w_rζ_r, Y0 = w_tρ + w·ζ.

**Step 1 (q).** Φ_qq = 2, so with ℓ(X) = tr X − wᵀXw/(2K) and m(ρ, ζ) = κρ + w_t(w·ζ)/(σ + S),

    Φ_q(X, ρ, ζ) := min_q Φ = F(X) + ⟨ℋ(ρ, ζ), X⟩ + 𝒜(ρ, ζ),   𝒜 = 2ρ² + |ζ|² − ½m²,

where F is the (C_u) form of Theorem C (q eliminated) and ⟨ℋ, X⟩ = ρ tr X − ℓ(X)m − 2Σ_r⟨𝒯_{(t,2)}ψ_r, W𝒯_{(r,1)}y⟩.
Reading off the coefficient of X_rh (with Σ_r w_r𝒯_{(t,2)}ψ_r = ½(wᵀXw)e_t − ½w_t Xᵀw):

    ℋ_rh = (ρ − m)δ_rh + m w_rw_h/(2K) − c_r W_tr w_h + w_t c_r W_hr − ½(Wy)_t w_rw_h + ½w_t w_r (Wy)_h.

**Step 2 (splitting).** For X = X′ + Ω, Ω = a wᵀ + w cᵀ, a ⊥ w, X′w = 0 = wᵀX′, Theorem C gives
F(X) = F(X′) + σ f(a, c) + σ⟨𝒢(a, c), X′⟩ ((F1)–(F3) of `typeI_u2/NOTES.md`, with
𝒢 = 2E⁻¹diag(a∘w) + 2diag(c∘w)E⁻¹ + p(E⁻¹ − I/K), p = w·c, E = diag(s_h + S/2 + σ) + ½wwᵀ).

**Step 3 (the representative).** X′ = ΠX′Π with Π = I − wwᵀ/S, so ⟨M, X′⟩ = ⟨ΠMΠ, X′⟩ for every M, and every term of
the form x wᵀ or w yᵀ is annihilated. Hence ΠℋΠ = Πℋ̂Π with

    ℋ̂ = w_t [ −(Y0/(σ + S)) diag(ρ⁰) + diag(w∘ζ/d⁰_U) − (1/(2 + g̃)) ((Y0/2) η̃ + w∘ζ∘η̃) η̃ᵀ ],

ρ⁰_h = s_h/d⁰_h, η̃ = w/d⁰_U, g̃ = σ/d⁰_t + Σ_h ρ⁰_h. (The diagonal uses (ρ − m) + w_tY0/(2d⁰_r) = −w_tY0ρ⁰_r/(σ + S): the
order-one couplings ρ tr X′ and −ℓ(X′)m of the trace direction cancel against the cross term up to ρ⁰_r. Before q
is eliminated, q itself couples to tr X′ at order one, and the same relaxation then gives an indefinite form,
`relax3.py`, Section 5.) In z-coordinates Y0 = w_t(ρ + Σ_h s_h z_h) and w∘ζ = w_t(s∘z), so ℋ̂ carries an overall
factor σ, like σ𝒢: the hard part couples to everything else at order σ.

**Step 4 (relaxation).** Lemma H gives F(X′) ≥ Σ_{r,h} λ_rh X′_rh² with λ_rh = 5/12 (r ≠ h) and λ_hh = Λ_h = 3/2 − ρ_h,
ρ_h = s_h/(s_h + S/2 + σ). Completing the square entry by entry (minimising over all real X′_rh, a relaxation of the
constraints):

    Φ_q(X′ + Ω, ρ, ζ) ≥ Ψ(a, c, ρ, ζ) := σ f(a, c) + ⟨ℋ(ρ, ζ), Ω⟩ + 𝒜(ρ, ζ) − ¼ Σ_{r,h} R_rh²/λ_rh,
    R = Π(σ𝒢(a, c) + ℋ(ρ, ζ))Π.

**Lemma R1.** For every u ≥ 2, all positive weights, all hard X′, a ⊥ w, c ∈ R^U, ρ ∈ R and ζ ∈ R^U, the inequality above
holds. Consequently G(ρ, ζ) ≥ inf_{a⊥w, c} Ψ(a, c, ρ, ζ) for all (ρ, ζ), and if Ψ ≻ 0 on {a ⊥ w} × R^U × R × R^U, then
G ⪰ G_rel ≻ 0, where G_rel is the Schur complement of Ψ onto (ρ, ζ).

*Proof.* Steps 1–4: Step 1 is the closed form with q eliminated and ℋ read off from it; Step 2 is (F1)–(F3) of Theorem C;
Step 3 uses only X′ = ΠX′Π; Step 4 is Lemma H and λ_rh X² + R X ≥ −R²/(4λ_rh). Every X ∈ R^{U×U} is X′ + Ω for a
unique such triple (X′, a, c), and G(ρ, ζ) = min_X Φ_q(X, ρ, ζ) because the (C_u) block is positive definite
(Theorem C). ∎

Remark. The (a, c) block of Ψ is not Theorem C's Φ1 (the representative of the coupling to X′ is the projected one,
Π𝒢Π, not 𝒢), so its positivity is not covered by Theorem C; it is part of the statement Ψ ≻ 0 checked below.

**Why the off-diagonal constant matters.** Lemma H's constant 5/12 replaces ½ − ρ_rρ_h/(2 + g); with the latter, 1/λ_rh
would be a pair function again. With a constant off-diagonal weight, Σ R²/λ = (12/5)‖R‖_F² + Σ_h (1/Λ_h − 12/5)R_hh², and
every matrix that enters (E⁻¹ = D⁻¹ − ηηᵀ/(2 + g), W, Π, 𝒢, ℋ̂) is diagonal plus a correction of rank at most two, with
entries rational in s_h, σ and global sums over U. So Ψ is a sum over h ∈ U of local forms in (a_h, c_h, ζ_h) plus a
form in a bounded number of global linear functionals, and (Woodbury) G_rel is diagonal plus bounded rank in the
z-coordinates. This removes the pair-function obstruction of `typeI_G_m1/NOTES.md` Section 4, at the price of an
inequality. We have not written the reduced finite system out.

**Exact checks** (`check_identity.py`, `logs_check_identity.log`, PASS). For u = 2, …, 6 and three random rational
weight vectors each, against `cf.TypeI` over Q: (I1) Φ_q(X, ρ, ζ) = Φ_q(X, 0, 0) + ⟨ℋ, X⟩ + 𝒜 for a random full X;
(I2) Φ_q(Ω, 0, 0) = σf(a, c); (I3) the cross term = σ⟨𝒢, X′⟩; (I4) the LDLᵀ pivots of F − Σλ_rhX′² on the hard
space are positive; (I5) Ψ ≤ Φ_q(X′ + Ω, ρ, ζ) at random rational points; (I6) G − G_rel ⪰ 0 (positive pivots for
u ≥ 3; for u = 2 the hard space has dimension one, G − G_rel has rank one, and all principal minors are ≥ 0). In
floating point (`check_forms_float.py`) the same identities hold to 10⁻¹⁴ at 20 random configurations, and the
exact Ψ agrees with an independent floating-point construction by polarisation (`check_psi_float.py`).

## 3. How tight is the relaxation (floating point)

`relax1.py` builds the whole block from `cf.py`, eliminates q exactly, writes X in the basis (hard part, a wᵀ, w cᵀ),
and replaces the hard block by the diagonal weights, with the coupling projected onto the hard space. Smallest
eigenvalue in the frame of Section 1, exact G against relaxed G (40 weight vectors per u):

| u | 3 | 4 | 6 | 8 |
|---|---|---|---|---|
| exact | 0.38208 | 0.38228 | 0.38209 | 0.38221 |
| relaxed, Lemma H weights | 0.38208 | 0.38227 | 0.38208 | 0.38221 |

(with the tight weights ½ − ρ_rρ_h/(2 + g) the numbers agree to the same accuracy). `relax2.py` builds Ψ from the
formulas of Section 2 and compares representatives (24 weight vectors per u; exact / both couplings projected (Ψ) /
𝒢 unprojected and ℋ projected / both unprojected):

| u | exact | Ψ | 𝒢 unprojected | unprojected |
|---|---|---|---|---|
| 3 | 0.38224 | 0.38223 | 0.38216 | −0.86 |
| 5 | 0.38801 | 0.38779 | 0.36136 | −2.04 |
| 8 | 0.46935 | 0.46935 | 0.33238 | −2.24 |

So the projection of both couplings is what makes the bound tight; the u-dependence of the second column is the reason
not to reuse Theorem C's representative.

## 4. Exact pointwise certificates for larger u

`certify_points.py u_list n_list seed` builds Ψ exactly over Q (`psi_exact.py`; rational w_t, w_h, so all squared
weights are rational) and computes all its leading principal minors by fraction-free elimination over Z. Families,
cycled: log-uniform; one U-weight tiny (10⁻⁴ to 10⁻²); all U-weights small (S ≪ σ); all large (S ≫ σ); one U-weight
dominant; geometric spread; two tiny U-weights; σ tiny. Weights are rounded to three significant digits.

Final log lines (`logs_certify_small.log`, `logs_certify_large.log`; run by `run_certify.sh`, one process at a time
under `guard.sh`):

    u = 9: 24/24 weight vectors with all pivots of Psi positive; log2 of the smallest letter-2 pivot ratio -22.2; 6.5 s
    u = 10: 24/24 weight vectors with all pivots of Psi positive; log2 of the smallest letter-2 pivot ratio -21.6; 12.1 s
    PASS
    u = 12: 8/8 weight vectors with all pivots of Psi positive; log2 of the smallest letter-2 pivot ratio -27.0; 11.3 s
    u = 16: 8/8 weight vectors with all pivots of Psi positive; log2 of the smallest letter-2 pivot ratio -25.8; 59.8 s
    u = 20: 4/4 weight vectors with all pivots of Psi positive; log2 of the smallest letter-2 pivot ratio -24.8; 100.9 s
    stopped by hand during u = 24 (time budget); u = 24 and u = 32 not completed

(All eight lines for u = 3, …, 10 read 24/24.) Each line proves G ≻ 0 at those weight vectors only. The previous exact
checks (`typeI_G_m1/`, 280 points with u ≤ 10) needed the whole block of size u² + u + 2; Ψ has size 3u, which is
what makes u = 12 to 20 reachable here. The fraction-free elimination grows quickly with u (u = 24 did not finish four
points in about 25 minutes). The smallest letter-2 pivot of Ψ, relative to its diagonal entry, is about 2⁻²² to 2⁻²⁸ in
these families (degenerate weights); exact arithmetic resolves it.

## 5. Routes of the task, as measured

1. **Supersolutions (a Z-matrix basis).** Not found. For C itself the earlier measurement stands (frustrated cycles
   in the natural and in the symmetric/antisymmetric bases). On the relaxed system Ψ, which no longer contains the
   hard part, the scaled coordinates (a_h/w_h, c_h/w_h, ρ, z_h) give mixed signs over 200 weight vectors (u = 4) in
   most (a, c)–(a, c) and (c, z) entries; only the letter-2 block has a fixed Z-pattern (`signs1.py`). Without the
   constraint a ⊥ w, Ψ is indefinite (normalised smallest eigenvalue −0.11), so that constraint, which links all
   coordinates, cannot be dropped as Theorem C does for Φ1. We did not pursue this route further.
2. **Elimination in the other order.** The decisive point turned out to be the order of q and the hard part, not of
   the letter-2 states. Once q is eliminated exactly, the couplings of the hard part X′ to the letter-2 states are of
   order σ (Step 3), and relaxing X′ entrywise loses almost nothing (Section 3). This is Lemma R1. It reduces the
   problem to Conjecture R, but does not prove it. Control (`relax3.py`, 24 weight vectors per u): with q kept as a
   variable, the same entrywise relaxation (valid there too, since Φ(X′, 0, 0, 0) = F(X′) + ½(tr X′)²) gives an
   indefinite form, with normalised smallest eigenvalue −0.37, −1.78, −3.07, −4.44 for u = 2, 3, 4, 6, growing with
   u; this is consistent with the failure reported in Section B5 of `typeI_u2/NOTES.md`.
3. **Reduction of the pair functions.** Lemma R1 removes the pair functions (row 8). A reduction to Theorems EU and
   TL does not follow: those theorems are about G, not G_rel, and an argument by extreme points of the moment set
   would need weight distributions with more atoms than two (one per global sum that enters Ψ). Not pursued.
4. **Other.** None completed. A small adversarial search (`adv.py`: Nelder–Mead on the U log-weights in [−9, 9],
   σ = 1, six starts for u = 3, 4, 5) found no value of λ_min(G_rel) in the frame below 0.382; the searches ran to the
   boundary (all U-weights at e^{±9}), so the infimum seems to be approached as S/σ → 0 or ∞. The value near 0.382
   seen in every sweep is close to (3 − √5)/2 ≈ 0.38197; we have not derived it.

## 6. What a proof of Conjecture R needs, and what we expect

Ψ is well conditioned in the frame of Section 1 (smallest eigenvalue of G_rel ≈ 0.38), its two degenerate regimes are
one-dimensional (ρ as S/σ → 0; the merged pair (ζ ∝ w, c ∝ w) as S/σ → ∞), and the perpendicular part is stiff
(0.72–1) and weakly coupled (≤ 2·10⁻² in floating point, `typeI_G_u3/NOTES.md`). A natural plan is to keep the
global variables (ρ, the merged components of c and ζ) exact, bound the local variables (a, c ⊥ w, ζ ⊥ w) per
coordinate as in Section B4 of `typeI_u2/NOTES.md`, and close with a 3 × 3 matrix of global sums. The soft-direction
cancellations in both regimes are exact to first order in Ψ (they come from the exactly kept soft part), so per-
coordinate bounds must not touch the merged pair. We did not carry this out.

## 7. The reduced finite system (D + U C Uᵀ) and Lemma F

**Coordinates.** a_h = w_h ah_h, c_h = w_h ch_h, ζ_h = w_t w_h zh_h (h ∈ U), and ρ. In these coordinates every entry of Ψ
is rational in the squared weights. Put d = s + S/2 + σ, d0 = s + (σ + S)/2, dL = s + 3S/2 + 3σ (so Λ = dL/(2d)), and
d0_t = (3σ + S)/2.

**The decomposition** (`reduce_cas.py`). With v_h = (ah_h, ch_h, zh_h) and y = (ρ, F_1, …, F_30),

    Ψ = Σ_h v_hᵀ D(s_h) v_h + yᵀ C y,      F_j = Σ_h φ_j(s_h) x_h,

where (x, φ_j) runs over
- x = ah and x = ch: φ ∈ {1, s, s², 1/d, 1/d², 1/d³, 1/d⁴, 1/d0, 1/d0², 1/dL};
- x = zh: φ ∈ {1, s, s², 1/d, 1/d², 1/d0, 1/d0², 1/d0³, 1/d0⁴, 1/dL}.

So the columns of U are e_ρ and the vectors (φ_j(s_h))_h placed in the block of x. The constraint a ⊥ w is
F_{j0} = Σ_h s_h ah_h = 0. D(s) is a symmetric 3 × 3 matrix whose entries are polynomials in IG = 1/(2 + g),
IGT = 1/(2 + g̃), u, M_0 = Σ 1/d0_h, M_1 = Σ 1/d_h, M_2 = Σ 1/d0_h², M_3 = Σ 1/d_h², with coefficients rational in
(s, σ, S) (73 terms). C is a symmetric 31 × 31 matrix whose entries are polynomials in IG, IGT, u and the ten moment
sums M_k = Σ_h f_k(s_h), f_k ∈ {1/d0, 1/d, 1/d0², 1/d², 1/d³, 1/d0³, s², 1/dL, 1/d0⁴, 1/d⁴}, with coefficients rational
in (σ, S) (4952 terms). Since s/d = 1 − (S/2 + σ)/d, g = u − (S/2 + σ)M_1 and g̃ = σ/d0_t + u − ((σ + S)/2)M_0.

**Derivation.** A small symbolic program works with a generic index: every vector is an expression in s and the local
variables at that index, every matrix is a diagonal plus a sum of outer products, and a sum over U of a product of
entries is split into a local quadratic part (into D), a part linear in one local variable (a functional F_j, after
partial fractions in s over Q(σ, S)) and a constant part (a moment). It evaluates σf, 𝒜, ⟨ℋ, Ω⟩ (with the scaled matrix
of ℋ_full) and the relaxation term with R = Π(σ𝒢 + ℋ̂)Π written as a double centring of the scaled matrix. It does not
use `psi_exact.py`; the two constructions share only the formulas of Section 2, and the second uses ℋ̂ (after the
projection) where the first uses ℋ_full.

**Exact check** (`check_reduced.py`, PASS for two seeds): at 24 random rational weight vectors with u = 3, …, 8 the
matrix Tᵀ(blockdiag D + U C Uᵀ)T equals the matrix of `psi_exact.py` exactly, where T is the change of coordinates
(α, c, ρ, ζ) → (ah, ch, zh, ρ) (ah = Bm α / w, ch = c / w, zh = ζ/(w_t w)). This also checks ΠℋfullΠ = Πℋ̂Π.

**Lemma F.** Suppose D(s_h) ≻ 0 for every h and B_F := Σ_h U_hᵀ D(s_h)⁻¹ U_h ≻ 0, where U_h ∈ R^{3×30} holds the values
φ_j(s_h) in the row of the class of F_j. Let N = C + (0 ⊕ B_F⁻¹), a 31 × 31 matrix on y = (ρ, F), and let N_0 be N with
the row and the column of F_{j0} deleted. Then Ψ ≻ 0 on {a ⊥ w} if and only if N_0 ≻ 0.

*Proof.* For fixed y with y_{j0} = 0, the minimum of Σ_h v_hᵀ D_h v_h over the local vectors with F(v) = y_F is
y_Fᵀ B_F⁻¹ y_F (least norm in the D-metric; B_F ≻ 0 means that the 30 functionals are independent on the local space).
So the minimum of Ψ over this fibre is yᵀ N y. If y = 0 the fibre is ker F, on which Ψ = Σ v_hᵀ D_h v_h > 0 unless
v = 0. ∎

So condition (i) is "D(s_h) ≻ 0 for every h" and condition (ii) is "N_0 ≻ 0", a fixed-size condition whose entries are
rational in the 13 global quantities (σ, S, u, M_0, …, M_9) and in B_F = Σ_h U_hᵀ D_h⁻¹ U_h, itself a sum over h of
explicit rational functions of s_h with the global quantities as parameters. B_F has rank at most 3u, so it is singular
for u ≤ 9, and also whenever the U-weights take fewer than ten distinct values; then the same statement holds with the
pseudo-inverse of B_F and y_F restricted to its range.

**Exact checks of (i), (ii)** (`certify_reduced.py 10,11,12 8,4,4 31 both`): at 16 rational weight vectors (i) holds at
all 16, B_F ≻ 0 at 15 (one u = 10 vector, family 0, has B_F singular), (ii) holds at all 15, and the verdict agrees with
the direct certificate of Ψ at all 15. Final log lines:

    u = 10: (i) D_h > 0 for all h: 8/8; B_F > 0: 7/8; (ii) N > 0: 7/8, agreement with the direct certificate 7/8; 186.3 s
    u = 11: (i) D_h > 0 for all h: 4/4; B_F > 0: 4/4; (ii) N > 0: 4/4, agreement with the direct certificate 4/4; 150.0 s
    u = 12: (i) D_h > 0 for all h: 4/4; B_F > 0: 4/4; (ii) N > 0: 4/4, agreement with the direct certificate 4/4; 202.6 s

(The last word of that log is FAIL only because that version of the script counted the singular-B_F vector as a
failure; the script was then changed to report such vectors as "route not applicable".)

**Conditioning.** The partial-fraction basis is badly conditioned: the ten functions per class are nearly dependent
for spread weights, and C has large entries that cancel. In floating point the assembled matrix and N are meaningless
(a first attempt gave "eigenvalues" of order −10⁹ at points where Ψ is certified positive; that script was discarded).
In exact arithmetic the reduced route was slower than the direct certificate at u = 20 (stopped after about seven
minutes on the first vector, `logs_certify_reduced_time.log`), because B_F⁻¹ has very large entries. A better
conditioned basis (for example the quantities ρ_h, ρ⁰_h, η_h² of Theorem C instead of the powers 1/d^k) would be needed
before (ii) can be bounded by hand.

## 8. Theorem R_EU: Conjecture R at equal U-weights, every u

**Theorem R_EU.** Let all U-weights be equal, s_h = s, and let u ≥ 3 be real. Then Ψ ≻ 0 on {a ⊥ w} for all σ, s > 0.
Consequently (Lemma R1, Corollary R1) G ≻ 0, and the type I block of (R∞)_{k−1} with |T| = 1 is ⪰ 0 with kernel
exactly span{p ⊗ f_h}. (The conclusion G ≻ 0 is already contained in Theorem EU. Theorem R_EU is a second proof of it
through the relaxation, and shows that Lemma R1 loses no positivity at equal weights.)

*Proof (computer-assisted, symbolic in u; `equal_R.py 3`).* At equal weights Ψ commutes with the permutations of U.
- Standard part: vectors with Σ_h v_h = 0 have F(v) = 0 (each φ_j(s_h) is the same number) and satisfy the constraint;
  on them Ψ = Σ_h v_hᵀ D(s) v_h, and they are Ψ-orthogonal to the invariant vectors. So this part is ≻ 0 iff D(s) ≻ 0.
- Invariant part: v_h = (ah, ch, zh) for all h; the constraint u s ah = 0 forces ah = 0, and Ψ is a 3 × 3 form Ψ_inv
  on (ch, zh, ρ), with y_F = u φ_j(s)·ch or u φ_j(s)·zh.
- All entries are rational in (σ, s, u), with S = us, M_k = u f_k(s), IG = 1/(2 + us/d), IGT = 1/(2 + σ/d0_t + us/d0).
  With u = 3 + x, the three leading principal minors of D(s) and of Ψ_inv are ratios whose numerators have only
  positive coefficients as polynomials in (σ, s, x) (15, 30, 99 and 39, 130, 126 terms; degrees 10, 17, 31 and 16, 33,
  33) and whose denominators are products of polynomials with only positive coefficients and a positive constant. By
  Sylvester's criterion both parts are positive definite for σ, s > 0 and x ≥ 0. ∎

With u = 1 + x the minors of Ψ_inv are still subtraction-free, but the second and third minors of D(s) have coefficients
of both signs, so the certificate does not cover 1 < u < 3 (`logs_equal_R_u1.log`; not needed here).

## 9. Negative controls

`controls_R.py` (exact; eight rational weight vectors, u = 4 and 6; `logs_controls_R.log`):
- (C0) Ψ ≻ 0 at all eight.
- (C1) The couplings between (a, c) and (ρ, ζ) multiplied by 1.001: Ψ is indefinite at the two vectors with S ≫ σ
  (w_t = 1/100, U-weights near 1), and stays positive at the other six.
- (C4) The coupling m of the letter-2 states to Q multiplied by 1.001: Ψ is indefinite at the four vectors with S ≪ σ
  (U-weights of order 10⁻² and 10⁻³, w_t = 1), where 𝒜_ρρ = 2 − κ²/2 → 0. So each degenerate regime has a constant with
  no slack.
- (C2) Without the projection Π: Ψ is indefinite at the four vectors with S ≳ σ.
- (C3) Lemma H's diagonal weight Λ_h replaced by 5/12 (weaker but still valid): Ψ stays positive at all eight, so the
  relaxation has room on the diagonal.

`certify_reduced.py 10 8 31 controlB` (`logs_certify_reduced_controlB.log`): with B_F replaced by 1.001 B_F in
condition (ii), N_0 fails to be positive definite at three of the seven vectors where the route applies (all U-weights
large; one U-weight dominant; σ tiny), and it holds unperturbed at all seven (row 8b). (The log line for family 0 is the
vector with B_F singular, where the route does not apply.)

## 10. What blocked a proof for all weights, and the classes with two or three values

- **Condition (i).** (This was open when this section was written; it is now proved for all weights, Theorem D in
  Section 11.)
- **Condition (ii).** It is a 30 × 30 condition. Bounding its entries one by one is hopeless in the present basis,
  because of the cancellations described in Section 7.
- **Two or three U-weight values, symbolically.** At L distinct values, Ψ reduces exactly to the standard parts
  D(s_i) ≻ 0 and an invariant block of size 3L (minus the constraint) whose entries are rational in
  (σ, s_1, …, s_L, u_1, …, u_L). At equal weights (L = 1, three variables) the minors already reach degree 33. For L = 2
  (a 6 × 6 block in five variables) and L = 3 (9 × 9 in seven) the determinants would have degree well above 100; we
  did not attempt them in the time. So Conjecture R for three distinct values, which would be new beyond Theorem TL, is
  still open. Update: by Theorem D (Section 11) the standard parts D(s_i) ≻ 0 are now proved for every L, so only the
  invariant block remains. For one free U-weight and the others equal (L = 2 with u_1 = 1; three variables after
  normalisation, plus the level size) the invariant block is 6 × 6; the entries of C, after the moments are written for
  two levels, have common denominators of degree around 50 to 80, so the determinant would have degree of several
  hundred. We did not attempt it, and that class is covered for G by Theorem TL anyway. For three values, uniformly in
  the level sizes at fixed weight values, a rough count gives degree several hundred in three variables; also not
  attempted.

## 11. Theorem D: condition (i) for all weights

**The local block in natural sums** (`D_natural.py`, exact). With c = S/2 + σ and c0 = (σ + S)/2, the identities
s/d = 1 − c/d and s/d0 = 1 − c0/d0 give u − c M_1 = g, M_1 − c M_3 = N, u − c0 M_0 = g0, M_0 − c0 M_2 = N0, where

    g = Σ_k s_k/d_k,   N = Σ_k s_k/d_k²,   g0 = Σ_k s_k/d0_k,   N0 = Σ_k s_k/d0_k²   (sums over all k ∈ U).

Substituting, u cancels identically (degree 0 in u, symbolic), and D(s) becomes a rational function of
(s, σ, S, g, N, g0, N0) only (with IG = 1/(2 + g), IGT = 1/(2 + σ/d0_t + g0)): the (ah, ch) block depends on (g, N),
the zh entry on (g0, N0), the cross entries on (g, g0). N enters only D_aa and D_cc, N0 only D_zz, each linearly.
Check (`check_D.py`, exact): at 24 rational weight vectors (u = 3, …, 8) the natural-sum form evaluated at the true
sums equals the local blocks of the reduced system, for all 132 blocks.

**Theorem D.** For every u ≥ 1, all positive weights and every h ∈ U, the local block D(s_h) of the reduced system
(Section 7) is positive definite.

*Proof (computer-assisted, exact; `prove_D.py`).*
- **Part A (monotonicity).** The coefficient of N in D_aa and in D_cc, and that of N0 in D_zz, is minus a polynomial with
  only positive coefficients (10 terms each) over a denominator with only positive coefficients, in the positive
  variables (s, σ, S, IG, IGT, g, g0). Since d_k ≥ c and d0_k ≥ c0, N ≤ g/c and N0 ≤ g0/c0. Hence
  D ⪰ Dt := D(N = g/c, N0 = g0/c0), the difference being diagonal and positive semidefinite.
- **Part B (the region).** Dt depends on (s, σ, S, g, g0). Normalise K = σ + S/2 = 1: every entry of D satisfies
  D(λs, λσ, λS, g, N/λ, g0, N0/λ) = λ² D(s, σ, S, g, N, g0, N0) (checked exactly, `check_D_homog.py`; g, g0, IG, IGT are
  scale invariant), so positivity is unchanged. Then t = σ ∈ (0, 1), S = 2(1 − t), b = c0 = 1 − t/2,
  s = S ξ with ξ ∈ (0, 1]. Because g′ = Σ_{k≠h} s_k/(1 + s_k) ∈ [0, S − s] and g0′ = Σ_{k≠h} s_k/(s_k + b) ∈
  [0, (S − s)/b], every configuration has g = s/(1 + s) + θ(S − s) and g0 = s/(s + b) + θ0(S − s)/b for some
  θ, θ0 ∈ [0, 1].
- **Part C (certificate).** The three leading principal minors of Dt are rational functions of (t, ξ, θ, θ0). Every
  irreducible factor of every numerator and denominator is shown to have a fixed sign on t ∈ (0, 1), ξ ∈ (0, 1],
  θ, θ0 ∈ [0, 1]: the monomial factors t, ξ trivially; the factor t − 1 by Bernstein coefficients with the corner rule at
  t = 1; every other factor (the largest has degrees (19, 20, 4, 2) in (t, ξ, θ, θ0) and 4641 terms) by its tensor
  Bernstein coefficients on [0, 1]⁴, all of one sign without subdivision. Multiplying the signs, the three minors are
  positive, so Dt ≻ 0 by Sylvester's criterion, and D ⪰ Dt ≻ 0. ∎

Last lines of `logs_prove_D.log`: `m3: certified sign 1 (3.5 s)`, `PASS`.

**Checks and controls** (`check_D.py`, PASS): at the 24 weight vectors above, every D(s_h) is positive definite in exact
arithmetic, the parameters (t, ξ, θ, θ0) computed from the weights lie in the certified domain, and D(s_h) − Dt(s_h) is
diagonal with nonnegative entries. Control: Dt − diag(Dt)/10 is not positive definite at t = 9999/10000, ξ = 1,
θ = θ0 = 0, so the certificate has a margin below 1/10 of the diagonal there (floating point over the region: smallest
Jacobi-normalised eigenvalue of Dt about 0.089, `ex_Dbox2.py`). A second control, D evaluated at the false bounds
N = 2g/c, N0 = 2g0/c0, stayed positive definite on the grid tried; so the bound N ≤ g/c is not the critical step.

**What Theorem D gives.** With Lemma F: for every weight vector with B_F ≻ 0, Conjecture R holds if and only if condition
(ii) holds. Theorem D does not by itself prove Conjecture R or G ≻ 0 anywhere new.

## 12. Condition (ii) and the direct route: what was measured

- **No slack.** Condition (ii) fails at explicit points when B_F is replaced by 1.001 B_F (Section 9), and Ψ fails when
  the soft couplings or the Q-coupling are scaled by 1.001. So a certificate must be exact in the soft directions; a
  per-copy Cauchy–Schwarz bound of the bounded-rank part (Ψ ⪰ Σ_h v_hᵀ(D_h − W_h)v_h) cannot work on its own, because it
  loses a constant factor in exactly those directions.
- **Where the soft directions are** (`ex_globperp.py`, floating point from exact matrices, u = 4, 6, 8, 16 weight vectors
  each). Split c = γ w/|w| + c⊥ and ζ = ξ w/|w| + ζ⊥ with c⊥, ζ⊥ ⊥ w (a ⊥ w already). After Jacobi scaling, the
  perpendicular block on (a, c⊥, ζ⊥) has smallest eigenvalue at least 0.49, the global 3 × 3 block on (γ, ξ, ρ) carries
  all the degeneracy (down to 4·10⁻¹⁰), and the coupling satisfies ‖C_gp‖² ≤ 0.67 λ_min(global) λ_min(perp). So a
  proof could keep (γ, ξ, ρ) exact (a 3 × 3 matrix of sums), bound the perpendicular block per copy (Theorem D is the
  local ingredient), and control the coupling by a Schur complement.
- **The rank of C** (`ex_inertiaC.py`, eigenvalues in 300-digit arithmetic from exact entries, u = 4 and 8, eight weight
  vectors): C has rank 16 (15 or 16 zero eigenvalues), with 7 negative and 8 or 9 positive eigenvalues, of size up to
  10²⁶ with cancellations. The kernel of C is not determined by (σ, S) alone: at two weight vectors with the same σ, S
  and u (and at an equal-weight vector with the same σ, S) the kernels differ (`ex_kerC.py`, exact). So no change of
  basis depending only on (σ, S) reduces the 31 functionals to 16; a well-conditioned rewrite of C has to use
  moment-dependent combinations, or (better) be restricted to the perpendicular subspace as below.
- **The perpendicular part is dominated by the local part** (`ex_perpeps.py`, exact data, 120-digit generalized
  eigenvalues, eight weight vectors per u). On the subspace p = {Σ s_h ah_h = Σ s_h ch_h = Σ s_h zh_h = 0, ρ = 0} let
  ε* = max (−yᵀCy)/(Σ_h v_hᵀD_hv_h). Then Ψ ⪰ (1 − ε*)Σ_h v_hᵀD_hv_h on p, and

  | u | 3 | 4 | 6 | 8 | 10 | 12 | 16 |
  |---|---|---|---|---|---|---|---|
  | max ε* | 0.31 | 0.15 | 0.11 | 0.066 | 0.016 | 0.12 | 0.008 |

  On p the soft-direction terms vanish (ρ = 0, p = w·c = 0, Y0 = 0, m = 0, so 𝒜 = |ζ|² is local), and the remaining
  global terms are products of a bounded number of natural sums such as Σ_h s_h²ah_h/d_h², Σ_h ρ_h ah_h. So a per-copy
  Cauchy–Schwarz bound of the global part on p has a factor of about 3 to spare.
- **The whole plan, tested** (`ex_plan.py`, same data). With g = (γ, ξ, ρ) (ch = γ·1, zh = ξ·1, ρ) and the subspace p:
  if Ψ_pp ⪰ κ D_pp, then Ψ ≻ 0 follows from the 3 × 3 condition Ψ_gg − Ψ_gp (κ D_pp)⁻¹ Ψ_gpᵀ ≻ 0. For κ = 0.65 and
  κ = 0.85 this 3 × 3 matrix is positive definite at all 32 weight vectors tried (u = 3, 4, 6, 8); its smallest
  normalised eigenvalue (down to 6·10⁻¹⁰) is the same for both κ, i.e. it is set by the soft directions and not by the
  coupling to p. So a proof of Conjecture R can be organised as: Theorem D (proved), a per-copy bound ε ≤ 0.35 on p (to
  do; measured ε* ≤ 0.31), and an exact 3 × 3 condition in natural sums that keeps the soft directions exact (to do).

## 13. The two-part plan for Conjecture R: lemmas, exact pointwise checks, and what did not close

Notation: V = {Σ_h s_h ah_h = 0} (the constraint a ⊥ w), g = (γ, ξ, ρ) with ch = γ·1, zh = ξ·1, and
p = {Σ s_h ah_h = Σ s_h ch_h = Σ s_h zh_h = 0, ρ = 0}; then V = g ⊕ p. D_pp is the local part Σ_h v_hᵀ D(s_h) v_h
restricted to p (positive definite by Theorem D).

**Lemma P1 (assembly).** If for some κ ∈ (0, 1]
- (A_κ) Ψ_pp − κ D_pp ⪰ 0 on p, and
- (B_κ) Ψ_gg − Ψ_gp (κ D_pp)⁻¹ Ψ_gpᵀ ≻ 0 (a 3 × 3 matrix),

then Ψ ≻ 0 on V, i.e. Conjecture R holds at that weight vector. *Proof.* (A_κ) gives Ψ_pp ⪰ κD_pp ≻ 0, hence
Ψ_pp⁻¹ ⪯ (κD_pp)⁻¹, and the Schur complement of Ψ onto g dominates the matrix of (B_κ). ∎

**Lemma P2 (per-copy form of the coupling term).** Let b_h be the 3 × 3 coupling of g to the local variables of copy h
(the full-space representative), and B0 = Σ_h b_hᵀD_h⁻¹b_h, B1 = Σ_h s_h D_h⁻¹b_h, B2 = Σ_h s_h² D_h⁻¹. Then
Ψ_gp D_pp⁻¹ Ψ_gpᵀ = B0 − B1ᵀ B2⁻¹ B1. *Proof.* Minimise Σ_h v_hᵀD_hv_h − 2Σ_h b_h(g)ᵀv_h under the three constraints
of p with a Lagrange multiplier λ ∈ R³ (one per class): the optimum is min_λ Σ_h (b_h − s_hλ)ᵀD_h⁻¹(b_h − s_hλ). ∎
So (B_κ) is a condition on Ψ_gg and on three sums over copies of explicit 3 × 3 matrices. The multiplier is essential:
with the full-space representative alone (λ = 0) the condition fails at all 24 weight vectors tried (u = 4, 6, 8;
`ex_Bcrude.py`), since that representative carries large components along the constraint directions.

**Exact checks.**
- `check_Blambda.py` (PASS): Lemma P2's identity holds exactly at 16 weight vectors (u = 3, 4, 6, 8).
- `check_plan.py` (PASS): (A_1/2) and (B_1/2) both hold in exact arithmetic at 56 rational weight vectors (8 each for
  u = 3, 4, 5, 6, 8, 10, 12; all eight families, including the degenerate ones), so Conjecture R holds there by Lemma
  P1 (also known from the direct certificates). Control: (A_99/100) fails at 2, 4, 3, 2, 1, 1, 1 of the 8 vectors for
  u = 3, 4, 5, 6, 8, 10, 12, so the global part on p is genuinely negative and (A) is not the statement D ≻ 0 alone.
- `ctl_plan.py` (control for (B)): with the coupling term inflated, (B_κ) fails at 3 of 8 vectors for κ = 1/100 and at
  5 (u = 4) and 4 (u = 6) of 8 for κ = 1/1000, while it holds at all 8 for κ = 1/20 and κ = 1/2. So the coupling to p is
  not negligible in (B); the margin in κ is between 1/100 and 1/20 at these vectors.

**What did not close (this session).**
- **(A_κ) for all weights.** Not proved.
  - A natural-functional version of the reduced system on p (`reduce_perp.py`: no partial fractions, the constraint
    functionals set to zero) has 95 functionals with complicated weights, too many for term-by-term Cauchy–Schwarz.
  - A direct per-copy lower bound was tried (`ex_Qlb.py`): T ≤ (12/5)‖R̂‖²_s with ‖A + B‖² ≤ (1 + δ)‖A‖² +
    (1 + 1/δ)‖B‖², the centring and the variances N − g²/S, N0 − g0²/S for the rank-one parts, and Cauchy–Schwarz for
    the coupling term (σ²IGT/2)Γ_cΓ_z. The resulting 3 × 3 blocks have negative diagonal entries for every tried
    (δ, θ). The loss is the dropped diagonal correction Σ_h (1/Λ_h − 12/5) s_h² R̂_hh² of the relaxation, which reduces
    the local part of T by a factor 2 to 3; a bound that keeps it needs the copy-by-copy pair kernel of T (products of
    natural sums such as Σ_h s_h² ah_h/d_h²), which we did not derive.
- **(B_κ) for all weights.** Not proved. Its soft directions require exact treatment, and by Lemma P2 it involves sums
  of the per-copy matrices D_h⁻¹, D_h⁻¹b_h, b_hᵀD_h⁻¹b_h, which are not functions of finitely many natural sums.
- So Conjecture R is still open for general weights; Lemma P1 together with Theorem D reduces it to (A_κ) and (B_κ) for
  one κ, both verified exactly at the 56 weight vectors above.

## 14. Files

| File | Content |
|---|---|
| `cf.py`, `flintfield.py`, `gzexact.py`, `cu_form.py` | copied from the repository (closed form of the block; rational functions; exact G_z; the (C_u) form), unchanged |
| `numlib.py` | floating-point block in the gauge x_t = 0 (scaled or not), G_z, the vector v |
| `forms.py` | floating-point formulas of Section 2: ℋ, 𝒜, 𝒢, f |
| `psi_exact.py` | exact Ψ over Q (python-flint), its letter-2 index, LDLᵀ pivots |
| `check_identity.py` | exact checks (I1)–(I6) of Lemma R1 (`logs_check_identity.log`, PASS) |
| `certify_points.py`, `run_certify.sh` | exact leading minors of Ψ at rational weight vectors (`logs_certify_small.log`, `logs_certify_large.log`) |
| `check_forms_float.py`, `check_psi_float.py` | floating-point checks of the formulas and of the exact construction |
| `relax1.py`, `relax2.py` | tightness of the relaxation, representatives (Section 3); `logs_relax1.log`, `logs_relax2.log` (the first eight lines of the latter repeat the former) |
| `relax3.py` | control: the relaxation before eliminating q is indefinite (`logs_relax3.log`) |
| `adv.py` | adversarial search for small λ_min(G_rel) (`logs_adv.log`) |
| `ex_conjM.py`, `ex_soft_dirs.py`, `ex_frame.py`, `signs1.py` | floating-point exploration (Sections 1, 5); no claim rests on them (logs `logs_ex_*.log`, `logs_signs1.log`) |
| `guard.sh` | runs one Python process with a memory cap of 1.0 GB and a time limit |
| `reduce_cas.py` | symbolic derivation of the reduced system (Section 7); writes `reduced.pkl` (SymPy) and `reduced_compact.pkl` (integers); `logs_reduce_cas.log` lists the 30 functionals and the 10 moment functions |
| `reduced_c.py` | exact evaluation of the reduced system from `reduced_compact.pkl` (D blocks, U, C, assembly) |
| `check_reduced.py` | exact check Tᵀ(D + U C Uᵀ)T = Ψ, u = 3…8 (`logs_check_reduced.log`, `logs_check_reduced_seed7.log`, PASS) |
| `certify_reduced.py` | conditions (i), (ii) of Lemma F at rational weight vectors; mode `both` compares with the direct certificate (`logs_certify_reduced_both.log`), mode `controlB` uses 1.001 B_F (`logs_certify_reduced_controlB.log`); `logs_certify_reduced_time.log` is the u = 20 run stopped by hand |
| `equal_R.py` | Theorem R_EU: symbolic minors at equal U-weights (`logs_equal_R_u3.log`, PASS; `logs_equal_R_u1.log`, see Section 8) |
| `controls_R.py` | negative controls C0–C4 (`logs_controls_R.log`) |
| `logs_rerun_*.log`, `logs_rerun2_*.log`, `logs_rerun3_*.log` | the exact checks rerun with other seeds (status row 12) |
| `*.pkl` | written by `reduce_cas.py` (about four minutes), `D_natural.py` and `equal_R.py`; not in the repository |
| `reduced_eval.py`, `show_D.py` | first (SymPy-based) loader of `reduced.pkl`, superseded by `reduced_c.py`; display of D(s) (`logs_show_D.log`) |
| `D_natural.py` → `D_natural.pkl` | D(s) in the natural sums g, N, g0, N0 (`logs_D_natural.log`: u cancels) |
| `prove_D.py`, `bern.py` | Theorem D (`logs_prove_D.log`, PASS); `bern.py` copied from `typeI_u2/` (tensor Bernstein certificates) |
| `check_D.py`, `check_D_homog.py` | exact checks and controls for Theorem D (`logs_check_D.log`, PASS); the homogeneity used in the normalisation (`logs_check_D_homog.log`, PASS) |
| `ex_Du.py`, `ex_Dbox.py`, `ex_Dbox2.py` | the u-independence at random points (exact) and the margin of D over the region (floating point) |
| `ex_globperp.py`, `ex_inertiaC.py`, `ex_kerC.py`, `ex_perpeps.py`, `ex_plan.py` | measurements of Section 12 (`logs_ex_globperp.log`, `logs_ex_inertiaC.log`, `logs_ex_kerC.log`, `logs_ex_perpeps.log`, `logs_ex_perpeps_large.log`, `logs_ex_plan.log`) |
| `check_plan.py`, `check_Blambda.py`, `ctl_plan.py` | exact checks and the control of Section 13 (`logs_check_plan.log`, `logs_check_plan_more.log`, `logs_check_Blambda.log`, PASS; `logs_ctl_plan.log`) |
| `reduce_perp.py` → `perp_compact.pkl` | the reduced system on p with natural functionals (95 of them; `logs_reduce_perp.log`); not used further |
| `ex_Qlb.py`, `ex_Bcrude.py` | the failed crude bounds of Section 13 (`logs_ex_Qlb.log`, `logs_ex_Bcrude.log`) |
