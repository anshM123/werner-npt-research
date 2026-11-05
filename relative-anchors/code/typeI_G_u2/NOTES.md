# G ≻ 0 for the type I blocks of (R∞)_{k−1} with |U| = 2

Ansh Mishra and Aryan Senthilkumar, 4 November 2026. Nothing here is refereed. Notation as in
`relative-anchors/README.md` (Lemmas B, N, Theorem C, Proposition U1) and `relative-anchors/code/typeI_u2/NOTES.md`
(Parts A and B). U = copies without a letter 2, u = |U|; T = copies carrying one letter 2, m = |T|; s_x = w_x²;
σ = s_T, S = s_U. All forms are in units of 2Y.

## Status

| Statement | Status |
|---|---|
| (S4′) structural form of G in the gauge w_Tᵀ X_T = 0: G_tt depends only on (s_t, σ, w_U), G_tt′ only on (s_t, s_t′, σ, w_U), for every m | **proved** (Section 1, read off the closed form `cf.py`; uses (S1), (S2) and Theorem C) |
| Closed forms of the 3 × 3 blocks G_tt, G_tt′ for u = 2, as rational functions of (s_t, s_t′, σ, s_3, s_4) | **exact** (symbolic elimination over Q(w_0, …, w_4), `closed_u2.py`); all entries even in every weight (checked); **exact** agreement with direct elimination of the whole block for m = 1, …, 5 (115 blocks, `check_closed.py`) |
| (B_2): λ_min(Ĝ_tt) > (2/5)(1 − s_t/σ) for every t, m, and all weights (Ĝ_tt ≻ 0 if m = 1) | **proved**, computer-assisted (exact Bernstein certificates, `prove_B.py`) |
| (A_2): ‖Ĝ_tt′‖_F < (2/5)√(s_t s_t′)/σ for t ≠ t′, every m ≥ 2, all weights | **proved**, computer-assisted (exact Bernstein certificate, `prove_A.py`) |
| **Theorem U2.** G ≻ 0 for u = 2, every m ≥ 1 and all positive weights. Hence every type I block of (R∞)_{k−1} with \|U\| = 2 is ⪰ 0 with kernel exactly span{p ⊗ f_h : h ∈ U} (every k ≥ 3) | **proved** (weighted block Gershgorin from (A_2), (B_2); then Theorem C and Section A3 of `typeI_u2/NOTES.md`) |
| Negative controls: (A) with H = 1/4, (B) with H = 4/5 are false (explicit rational witnesses) | **exact** (`controls.py`) |
| G ≻ 0 for u ≥ 3 | **open** (not attempted beyond the observations in Section 5) |

So the type I blocks of (R∞)_{k−1} are now settled for |U| ∈ {0, 1, 2} (Propositions U0, U1 and Theorem U2 with
Theorem C). The cases 3 ≤ |U| ≤ k − 1 with |T| ≥ 1 remain open. Nothing here says anything about Conjecture 12 or the
endpoint problem beyond these blocks.

## 1. The structure of G (every u, every m)

The kernel vectors κ_c = p ⊗ Σ_h c_h f_h (X_{r,h} = w_r c_h for all r, ρ = ζ = q = 0) satisfy Yκ_c = 0 (Lemma N). So the
form is unchanged when a κ_c is added, and the Schur complement onto y = (ρ_t, ζ_t)_t after eliminating the other
variables does not depend on which complement of span{κ_c} is used for them (provided the eliminated block is positive
definite there). Section A3 of `typeI_u2/NOTES.md` uses {P_{r0,h} = 0}; here we also use the gauge w_Tᵀ X_T = 0
(η = 0 in the notation of (S3)), which is a complement because κ_c has w_Tᵀ X_T = σ cᵀ.

In that gauge the eliminated variables are X_U (all u² entries), q, and X′_T with Σ_t w_t ψ_t = 0 (ψ_t = row t of X).
- By (S1), (S2) the X′_T block is the identity and X′_T does not couple to (X_U, q). The (X_U, q) block is exactly the
  (C_u) block of Theorem C, positive definite, depending only on (σ, w_U).
- Reading the closed form (`cf.py`, Section A1 of `typeI_u2/NOTES.md`):
  - y–y terms: ½(Σ_t ρ_t)² + (3/2)|ρ|² + Σ_t |ζ_t|² minus the pair terms ⟨𝒯_{(t′,2)}y_t, E_{T∖tt′}^{−1}𝒯_{(t,2)}y_{t′}⟩;
  - y–(X_U, q) terms: (Σ_t ρ_t) tr_U X, the q-row −2qΣ_t(2Kρ_t + w_t w_U·ζ_t)/(2K − s_t), and the terms of
    −2Σ_{r∈U,t}⟨𝒯_{(t,2)}ψ_r, E_{T∖t}^{−1}𝒯_{(r,1)}y_t⟩; each is a sum over t of a map of y_t that depends on
    (w_t, σ, w_U) only;
  - y–X′_T terms: for r ∈ T ∖ t, 𝒯_{(r,1)}y_t = w_r y_t, so Σ_{r∈T∖t} contributes ⟨𝒯_{(t,2)}(ξ − w_tψ_t), E_{T∖t}^{−1}y_t⟩
    with ξ = Σ_r w_rψ_r = 0 in this gauge. So ψ_t couples only to y_t, through ψ_t·ℓ_t(y_t) with
    ℓ_t = −ζ_t − 2𝒯_{(t,2)}ᵀE_{T∖t}^{−1}(𝒯_{(t,1)} − w_t)y_t, which depends on (w_t, σ, w_U) only.
- Minimising Σ_t(|ψ_t|² + ψ_t·ℓ_t) under Σ_t w_tψ_t = 0 gives −¼Σ_t|ℓ_t|² + |Σ_t w_tℓ_t|²/(4σ).

Hence, with C_t the coupling of y_t into (X_U, q) and C the (C_u) block,

    G_tt  = diag(2, I) − ¼ℓ_tᵀℓ_t + (s_t/(4σ)) ℓ_tᵀℓ_t − C_tᵀ C^{−1} C_t,
    G_tt′ = ½ e_ρe_ρᵀ − (pair term of (t, t′)) + (w_t w_t′/(4σ)) ℓ_tᵀℓ_t′ − C_tᵀ C^{−1} C_t′      (t ≠ t′),

where E_{T∖t} has d_x = s_x + σ + S/2 − s_t/2 (x ∈ U), K + s_t/2 (x = t), and E_{T∖tt′} has
d_x = s_x + σ + S/2 − (s_t + s_t′)/2. Every ingredient depends only on the stated weights, which proves (S4′) for
every m (this is the structural reason behind (S4) of Part A). In particular the closed forms computed from one
configuration with m = 3 (T = {t, t′, v}, s_v = σ − s_t − s_t′) are the blocks for every m; for m = 2 (s_v = 0) and for
G_tt with m = 1 they extend by continuity, and the exact checks of Section 2 confirm them directly for m = 1, …, 5.

## 2. Closed forms for u = 2 (`closed_u2.py`, `check_closed.py`)

The whole type I block with T = {0, 1, 2}, U = {3, 4} is built by `cf.TypeI` over the rational function field
Q(w_0, …, w_4) (python-flint), the rows P_{3,h} are deleted (gauge of Section A3), and all other P_{r,h} and Q are
eliminated by exact Gaussian elimination (15 s). In the coordinates ζ_{t,h} = w_t w_h z_{t,h} every entry has only even
exponents, so it is a rational function of s_0, …, s_4. The 3 × 3 blocks G_00 (coordinates ρ_0, z_{0,3}, z_{0,4}) and
G_01 have common denominators

    den_00 = σ (σ + 2S) U_6 (2K − s_0)² E_0²,       den_01 = σ (σ + 2S) U_6 (2K − s_0)(2K − s_1) E_0 E_1 E_01,

with K = σ + S/2, E_0, E_1 cubic and E_01 quartic irreducible factors (they come from the inverses of E_{T∖t} and
E_{T∖{0,1}}; we did not identify them further, their signs are certified), and U_6 a degree-6 factor that depends
only on (σ, s_3, s_4) (it does not involve θ after normalisation; it comes from the (C_u) block). Numerators have degree ≤ 18 (G_00) and ≤ 22 (G_01). G_00 depends on s_1, s_2 only
through s_1 + s_2 (checked symbolically), as (S4′) says.

`check_closed.py` compares, exactly, the blocks of G obtained by direct elimination of the whole block over Q (random
rational weights, so the z-scaling is rational) with the closed forms evaluated at (s_t, σ − s_t, 0, s_3, s_4) and
(s_t, s_t′, σ − s_t − s_t′, s_3, s_4): m = 1, …, 5, 115 blocks, all equal. It also checks in floating point that the
two norm formulas below equal the quantities of `sweep.py` of the repo (relative error 2·10⁻¹⁵).

## 3. The two inequalities

Per t use the frame of Part A: ρ_t, the w_U-component of ζ_t rescaled by c = √(1 + S/σ), and the component ⊥ w_U; call
the blocks in this frame Ĝ (a congruence of G by a block-diagonal invertible matrix, so G ≻ 0 ⟺ Ĝ ≻ 0). With
Z_t = diag(1, w_t w_3, w_t w_4) and W = diag(1, I + w_Uw_Uᵀ/σ):
- ‖Ĝ_tt′‖_F² = tr(Ω_t G_tt′ Ω_t′ G_tt′ᵀ), Ω_t = Z_t^{−1}WZ_t^{−1} = diag(1, (diag(1/s_U) + 11ᵀ/σ)/s_t) (z-coordinates);
- λ_min(Ĝ_tt) > h ⟺ N := G_tt − h·diag(1, s_t(diag(s_U) − s_Us_Uᵀ/(σ + S))) ≻ 0 (z-coordinates).

Everything is homogeneous of degree 0, so normalise σ + S = 1 and put s_t = θσ, s_t′ = φσ, s_3 = μλ, s_4 = μ(1 − λ),
with σ = 1 − μ. Actual configurations: θ ∈ (0, 1] (θ = 1 iff m = 1), θ + φ ≤ 1, μ ∈ (0, 1), λ ∈ (0, 1). H = 2/5.

**(B_2)** (`build_B.py`, `prove_B.py`). Ñ = den_00·N has polynomial entries; det Ñ = θ²μ²λ(1 − λ)(1 − μ)⁶ P_3 with P_3 of
degrees (27, 45, 32) in (θ, μ, λ).
1. den_00 > 0 on the domain: each irreducible factor is certified sign-definite on the closed cube by Bernstein
   coefficients (U_6 > 0, σ + 2S = 1 + μ, 2K − s_t > 0, and −E_0 > 0 in this normalisation), except σ = 1 − μ itself.
2. P_3 vanishes on the edge θ = 1, μ = 0 (m = 1 and S → 0, where Ĝ_tt degenerates, as in Proposition U1). Two blow-up
   charts cover the cube: (θ, μ) = (1 − e, ev) and (1 − μv, μ). In each, the substituted P_3 is divisible by e
   (resp. μ) exactly once, and the quotient (degrees (71, 45, 32), resp. (71, 27, 32)) has all tensor Bernstein
   coefficients positive on [0, 1]³ (root box, no subdivision). So det N > 0 on the whole domain.
3. N ≻ 0 at (θ, μ, λ) = (1/2, 1/2, 1/3) (exact LDLᵀ).
4. The domain (0, 1] × (0, 1) × (0, 1) is connected, N is continuous on it and det N ≠ 0, so the inertia of N is
   constant; by 3, N ≻ 0 everywhere. ∎

**(A_2)** (`build_A.py`, `prove_A.py`). With Õm_t = s_t s_3 s_4 σ Ω_t (polynomial) and G_tt′ = Nm/den_01,
(A_2) ⟺ FA := H² s_t² s_t′² (s_3s_4)² den_01² − tr(Õm_t Nm Õm_t′ Nmᵀ) ≥ 0. FA = θ²φ²μ⁴λ²(1 − λ)²(1 − μ)⁶ P_A with P_A of
degrees (16, 16, 38, 28). After the onto map θ = x, φ = (1 − x)y of the triangle, P_A has all tensor Bernstein
coefficients positive on [0, 1]⁴ (degrees (24, 16, 38, 28), root box). The factors of den_01 are certified
sign-definite on the closed domain as in (B_2) (σ = 1 − μ > 0 aside). Hence ‖Ĝ_tt′‖_F < (2/5)√(s_t s_t′)/σ. ∎

**Theorem U2.** For every k ≥ 3, every T with |U| = 2 and all positive weights, G ≻ 0.

*Proof.* For m = 1, G = G_tt ≻ 0 by (B_2) at θ = 1. For m ≥ 2 and x = (x_t) in the frame of Section 3,
xᵀĜx ≥ Σ_t λ_min(Ĝ_tt)|x_t|² − Σ_{t≠t′}‖Ĝ_tt′‖₂|x_t||x_t′|, and 2|x_t||x_t′| ≤ √(s_t′/s_t)|x_t|² + √(s_t/s_t′)|x_t′|² gives
xᵀĜx ≥ Σ_t |x_t|² (λ_min(Ĝ_tt) − Σ_{t′≠t}‖Ĝ_tt′‖₂√(s_t′/s_t)). By (A_2), Σ_{t′≠t}‖Ĝ_tt′‖₂√(s_t′/s_t) <
(2/5)Σ_{t′≠t} s_t′/σ = (2/5)(1 − s_t/σ), and by (B_2) this is < λ_min(Ĝ_tt). So Ĝ ≻ 0, hence G ≻ 0. With Theorem C and
the reduction of Section A3 of `typeI_u2/NOTES.md`, the type I block is ⪰ 0 with kernel exactly span{p ⊗ f_h}. ∎

The certificates need no subdivision, so the margins are comfortable: numerically (Part A of `typeI_u2/NOTES.md`,
`sweep.py`) sup ‖Ĝ_tt′‖₂ σ/√(s_t s_t′) ≈ 0.278 and inf λ_min(Ĝ_tt)/(1 − s_t/σ) = 3/4, the same extreme constants as
u = 1; for the Frobenius norm used in (A_2) the largest value found is ≈ 0.288, near θ + φ = 1, μ → 0, λ small
(`ex_findA.py`).
Negative controls (`controls.py`): with H = 1/4 a rational point with FA < 0 is exhibited, and with H = 4/5 a rational
point where N has a negative LDLᵀ pivot; so neither constant can be moved past the numerical extremes.

## 4. What is computer-assisted, what is exact, what is numerical

- Proved by hand: Section 1 (structure of G in the gauge w_Tᵀ X_T = 0), the Gershgorin step, the inertia argument.
- Exact symbolic computation: the closed forms (rational function field over Q, python-flint), the evenness, the
  dependence of G_00 on s_1 + s_2 only, the factorisations, the polynomial reductions and divisibility.
- Exact certificates: tensor Bernstein coefficients (integers) of P_3 in two charts, of P_A, and of the denominator
  factors; exact LDLᵀ at one point.
- Exact checks (not used in the proof, but confirming it): closed forms against direct elimination for m = 1, …, 5.
- Numerical only: the extreme constants 0.278 and 3/4 and the location of the witnesses used by the controls.

## 5. u ≥ 3 (open)

- For u = 2 the relaxation route of Section B5 of `typeI_u2/NOTES.md` is exact (the trace-free part of X′ is empty),
  so the proof above is the exact problem; for u ≥ 3 that relaxation would be needed, and it was not carried out.
- The closed-form route does not extend uniformly in u as it stands: numerically, G depends on w_U through more than
  one scalar functional (matching one block entry between a u = 2 and a u = 3 configuration leaves the other entries
  different, `ex_g2.py`), and the term C_tᵀC^{−1}C_t′ involves the inverse of the whole (C_u) block. For fixed u the same
  pipeline applies in principle (G_tt′ then depends on u + 2 normalised parameters).
- Observed structure that may help: the components of ζ_t orthogonal to w_U are stiff (diagonal ≈ 0.75) and couple to
  the rest weakly (≤ 0.002 in the samples of `ex_g1.py`; exactly 0 at equal U-weights by symmetry), while the
  (ρ_t, w_U-component) part is close to the u = 1 blocks with the U-copies merged into one copy of weight √S.

## Files

| File | Content |
|---|---|
| `run_all.sh` | reruns every check below, one Python process at a time, prints PASS/FAIL per script and ALL PASS |
| `cf.py`, `lemmab.py`, `fastdirect.py`, `check_cf.py` | copied from `typeI_u2/`: closed form of the type I block and its exact check against the Lemma B builder |
| `bern.py`, `test_bern.py` | copied from `typeI_u2/`: exact tensor Bernstein certificates and unit tests |
| `flintfield.py` | copied from `relative-anchors/code/`: rational function field on python-flint |
| `closed_u2.py` → `closed_u2.pkl` | exact symbolic Schur complement, u = 2, m = 3; blocks G_00, G_01 over (s_0, …, s_4) |
| `check_closed.py` | closed forms = direct exact elimination, m = 1, …, 5 (115 blocks); norm formulas in floating point |
| `build_B.py`, `prove_B.py` | (B_2): minors, denominator signs, blow-up charts, one-point check |
| `build_A.py`, `prove_A.py` | (A_2): the polynomial FA, denominator signs, triangle map, certificate |
| `controls.py` | negative controls H = 1/4 for (A), H = 4/5 for (B), with exact witnesses |
| `ex_g1.py`, `ex_g2.py`, `ex_findA.py`, `diag_*.py`, `inspect_u2.py`, `sym_try*.py` | exploration only; no claim rests on them |
