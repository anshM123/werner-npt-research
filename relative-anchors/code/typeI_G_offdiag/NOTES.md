# Toward the off-diagonal bound for G: an exact reorganisation and a two-condition reduction

Ansh Mishra and Aryan Senthilkumar, 6 November 2026. Nothing here is refereed.

Notation as in `relative-anchors/README.md` (Theorem C), `../typeI_G_em/NOTES.md` Section 1 (the formula for G) and
`../typeI_G_diag/NOTES.md` (Theorem DG, Lemma GG). T carries the letter 2 (m = |T| ≥ 2), U is free (u = |U|),
σ = s_T, S = s_U, K = σ + S/2, θ_t = s_t/σ, Ĝ = PᵀGP in the frame of `../typeI_G_u3/numG.py`.

By Theorem DG and Lemma GG, every type I block of (R∞)_{k−1} at every weight vector would follow from the off-diagonal
bound ‖Ĝ_tt′‖ ≤ (3/4)√(θ_tθ_t′). That bound is **not proved**. This note records what is proved toward it: an exact
identity that removes the largest cancellation, a reduction to two explicit conditions with room on both sides, the
closed form of one of them, and two lemmas for it.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | **Lemma Q** (common-mode shift, Section 1) | **proved**; exact check against direct elimination of the whole block |
| 2 | **Lemma GG′** (one constant, Section 2) | **proved** |
| 3 | **Reduction.** G ≻ 0 at every weight vector follows from (O) ‖R̂_raw,tt′‖ ≤ κ√(θ_tθ_t′) for t ≠ t′ and (S) Ŝ_free,tt ⪯ (3/4 − κ)θ_t I for every t, for any κ ∈ (0, 3/4) | **proved** (Lemma GG′ and Theorem DG) |
| 4 | Closed form of R_raw,tt′ through five statistics of the free weights (Section 3) | derived by hand; **exact check** against R_raw at rational points |
| 5 | **Lemma L** (−2/9 < L(s) < 1) and **Lemma B** (block reduction of (O)), Section 4 | **proved** |
| 6 | (O) with κ = 3/8 at rational weight vectors, physical and formal free measures; control κ = 1/4 fails | **exact check** only |
| 7 | (O) and (S) for all weights | **open**. Floating point: sup ‖R̂_raw,tt′‖/√(θθ′) = 0.2778 (approached as σ/S → ∞, the same for formal measures); sup λ_max(Ŝ_free,tt)/θ_t = 0.279 for physical weights, and 1/3 for formal measures as σ → 0. So κ = 3/8 leaves room on both sides |

## 1. Lemma Q

Work in the gauge Σ_t w_tψ_t = 0. Let C be the (C_u) block on z = (X_U, q), C_t the coupling of y_t = (ρ_t, ζ_t) to z,
e_q the q-coordinate, e_ρ the ρ-coordinate, and recall C_qq = 2.

**Lemma Q.** Put δ_t = C_t + ½Ce_qe_ρᵀ. Then δ_t = 0 when w_t = 0, and for t ≠ t′

    G_tt′ = −P̃_tt′ + (w_tw_t′/(4σ)) L_tᵀL_t′ − δ_tᵀC⁻¹δ_t′,

with an explicit −P̃_tt′ in which every term carries the factor w_tw_t′ (formula in `lemmaq_exact.py`). So G = 𝓡 − 𝓓ᵀC⁻¹𝓓
with 𝓡 free of C⁻¹, and the constant ½e_ρe_ρᵀ of G_tt′ cancels exactly.

*Proof.* Substitute z = z′ + ½(Σ_tρ_t)e_q. The quadratic term −¼C_qq(Σρ_t)² = −½(Σρ_t)² cancels the term ½(Σρ_t)² of
½(tr X + Σρ_t)², and minimising over z′ gives the formula. The vanishing at w_t = 0 is checked column by column: the
ρ-column of C_t is ½ tr X − 2qK/(2K − s_t) minus the B-term, which at w_t = 0 equals ¼w_UᵀX_Uw_U/K, and the
ζ-columns vanish because the q-coupling and the B-term carry the factor w_t. ∎

Without the shift, the base, gauge and free parts of Ĝ_tt′ are each up to about 500 times the total at small θ. This
is why term-by-term bounds fail.

## 2. Lemma GG′

**Lemma GG′.** Let G = R − 𝓓ᵀ𝓜⁻¹𝓓 with 𝓜 ≻ 0, S_tt = 𝓓_tᵀ𝓜⁻¹𝓓_t, and suppose, for some κ > 0,
(O) ‖R̂_tt′‖ ≤ κ√(θ_tθ_t′) for t ≠ t′ and (D) R̂_tt − Ŝ_tt/θ_t ≻ κ(1 − θ_t)I for every t. Then G ≻ 0.

*Proof.* Since Σθ_t = 1, Cauchy–Schwarz gives |𝓜^{−1/2}Σ_t𝓓_tx_t|² ≤ Σ_t|𝓜^{−1/2}𝓓_tx_t|²/θ_t. Hence
x*Ĝx ≥ Σ_t x_t*(R̂_tt − Ŝ_tt/θ_t)x_t + Σ_{t≠t′} x_t*R̂_tt′x_t′ ≥ κ[Σ_t|x_t|² − (Σ_t√θ_t|x_t|)²] ≥ 0, strictly for x ≠ 0. ∎

Take R = R_raw, the form left after eliminating ψ, q and the common mean p (the free direction X_U ∝ w_Uw_Uᵀ) only, and
S = S_free, the Gram of the rest of the free elimination. Since R̂_tt − Ŝ_tt/θ_t = Ĝ_tt − ((1 − θ_t)/θ_t)Ŝ_tt and
Ĝ_tt ≻ (3/4)(1 − θ_t)I (Theorem DG), condition (D) follows from (S) Ŝ_free,tt ⪯ (3/4 − κ)θ_tI. This proves row 3 of the
status table.

## 3. The closed form of R_raw,tt′

Normalise S = 1 and view the free squared weights as a probability measure on [0, 1]. Write Z_h = ζ_h/w_h for the
field of a letter-2 copy. Then R_raw,tt′ is an explicit bilinear form in (ρ, Z) and (ρ′, Z′). The free measure enters
only through five statistics: F = E[1/(K + s)], three further values of the same Stieltjes transform, and one divided
difference. The fields enter through five linear statistics each, and a local multiplier

    L(s) = (K + s)²/((k_a + s)(k_b + s)) − σ/(k₀ + s),   a = s_t, b = s_t′, k_c = K − c/2, k₀ = K − (a + b)/2.

The formula is in `rawcf.py`. `check_exact.py` compares it entry by entry with R_raw, computed in exact rational
arithmetic from the whole block (`exact_offdiag.py`, which builds G with `lemmaq_exact.py` and eliminates ψ, q and p
exactly).

## 4. Lemma L and Lemma B

**Lemma L.** For σ > 0, a, b ≥ 0 with a + b ≤ σ and every s ≥ 0: −2/9 < L(s) < 1.

*Proof.* Put d = K + s, A = a/2, B = b/2, S′ = A + B ≤ σ/2, P = AB. The identity
(1 − L)(d − A)(d − B)(d − S′) = d(σ − S′)(d − S′) + P(σ + d − S′) has a positive right side, so L < 1. For the lower
bound, (d − A)(d − B) ≤ (d − S′/2)² and σ ≤ d give, with q = S′/d < ½,
L ≥ 1/(1 − q/2)² − 1/(1 − q) = −q²/(4(1 − q)(1 − q/2)²). Its absolute value increases on (0, ½] and equals 2/9 at
q = ½. ∎

**Lemma B.** Split the letter-2 space of one copy into span(e_ρ) ⊕ span(1) ⊕ H₀ (H₀: mean-zero fields), and write R̂
in blocks A (2 × 2 on (ρ, mean)), B, C (couplings to H₀) and D = ¼E[L Z°Z′°] + D_fin (on H₀). Then
‖R̂‖ ≤ ‖[[‖A‖, ‖B‖], [‖C‖, ¼ + ‖D_fin‖]]‖. So (O) with κ = 3/8 holds as soon as ‖A‖ ≤ x, ‖B‖, ‖C‖ ≤ β and
‖D_fin‖ ≤ d with x ≤ 3/8, d ≤ 1/8 and (3/8 − x)(1/8 − d) ≥ β².

*Proof.* The norm of a block operator is at most the norm of the matrix of block norms. The local part of D is the
compression of multiplication by L/4, of norm < ¼ by Lemma L. A nonnegative 2 × 2 matrix has norm monotone in its
entries, and [[x, β], [β, ¼ + d]] has norm ≤ 3/8 exactly under the stated conditions. ∎

The couplings and D_fin are bounded through variances of explicit functions, hence through Popoviciu's inequality.
Floating-point maps give ‖A‖ ≤ 0.2778, couplings ≤ 0.08 and ‖D_fin‖ ≤ 0.024. A certificate over the remaining
parameters (σ/S, θ_t, θ_t′ and the five statistics, each in a proved interval) has not been completed.

## 5. Checks

| File | What it checks | Log |
|---|---|---|
| `lemmaq_exact.py 1 8` | Lemma Q and the closed form of −P̃ against direct elimination of the whole block in another gauge, at 8 rational points (85 blocks); δ_t = 0, P̃ = 0 and G_tt′ = 0 at 4 points with w_t = 0 | `logs/lemmaq_exact.log` |
| `lemmaL_check.py` | the identities of Lemma L, the closed form of its lower bound, the value 2/9, the partial fractions of L (exact, SymPy) | `logs/lemmaL_check.log` |
| `check_exact.py 1 8` | the closed form `rawcf.py` equals R_raw entry by entry at 8 rational points (62 blocks) | `logs/check_exact.log` |
| `exact_gg2.py 11 12 3/8` | (O) and (D) of Lemma GG′ with κ = 3/8, and G ≻ 0 itself, at 12 rational points (74 off-diagonal, 35 diagonal blocks): no failure | `logs/exact_gg2_k38.log` |
| `exact_gg2.py 11 12 1/4` | control: κ = 1/4 breaks (O) at 4 blocks | `logs/exact_gg2_k14.log` |
| `exact_formal.py 5 20 3/8` / `1/4` | (O) on formal free measures (masses unrelated to locations) and physical ones, exact: κ = 3/8 holds at 20 points; κ = 1/4 fails at 2 | `logs/exact_formal_k38.log`, `logs/exact_formal_k14.log` |

Run `bash run_all.sh` (Python with python-flint, NumPy and SymPy). The builders reuse `../typeI_G_u3/cf.py`.
