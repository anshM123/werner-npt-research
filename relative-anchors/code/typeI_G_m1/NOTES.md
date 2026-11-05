# Conjecture M (one copy carrying the letter 2): every u when the U-weights take at most two values, exact checks for u ≤ 10, and why no closed form in moments exists

Ansh Mishra and Aryan Senthilkumar, 5 November 2026. Nothing here is refereed.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | For m = 1, in the gauge x_t = 0 the eliminated block is exactly the (C_u) block C, and G = diag(2, I_u) − Bᵀ C⁻¹ B. So G_z[i,j] = −⟨b_i, C⁻¹ b_j⟩ for i ≠ j: (Z) says that the couplings b_i of the letter-2 coordinates are pairwise positively correlated in the C⁻¹ metric (Section 1) | **proved** (read off the closed form); gauge independence **exact check** |
| 2 | Conjecture M for u = 1, 2, 3 and all positive weights | **proved** (exact closed forms in the gauge x_t = 0, coefficient signs); a re-derivation, in a second gauge, of the u ≤ 3 cases of typeI_G_u3 (same term counts) |
| 3 | **Theorem EU.** Conjecture M for every u ≥ 1 and every σ > 0 when all U-weights are equal | **proved**, computer-assisted, symbolic in u (Section 2) |
| 4 | **Theorem TL.** Conjecture M for every u ≥ 1, every σ > 0 and every U-weight vector with at most two distinct values (u₁ copies of s₁ and u₂ copies of s₂; all u₁, u₂, all s₁, s₂ > 0) | **proved**, computer-assisted, symbolic in u₁, u₂ (Section 3). It contains one small (or one dominant) U-weight with all others equal, at every σ/S; the exact checks of row 7 find the smallest margins of (D) for σ ≪ S in the row of the smallest U-weight |
| 5 | **Corollary.** For every k, every type I block of (R∞)_{k−1} with \|T\| = 1 whose U-weights take at most two values is ⪰ 0 with kernel exactly span{p ⊗ f_h : h ∈ U} | **proved** (rows 3, 4 with Theorem C) |
| 6 | Equal U-weights, u → ∞ at fixed σ, S: the (ρ, z) block in the frame of typeI_G_u3 Section 4 converges to an explicit matrix that is positive definite for all σ, S > 0. The values 0.579, …, 0.448 listed there are reproduced from the exact formula; the limit at σ = 0.004, S = 1 is 0.4221 | **proved** (exact limit and sign certificate, Section 2b); answers the question of typeI_G_u3 Section 4, item 3, at equal U-weights |
| 7 | (Z), (D) and G ≻ 0 at 28 rational points for each u = 1, …, 10 (280 points; nine families, among them one or two tiny U-weights, S/σ from 10⁻⁸ to 10⁸, geometric spreads); gauge independence at 40 of them. The smallest margins (down to 3·10⁻¹²) occur for σ ≪ S, in the row of the smallest U-weight | **exact check** (Section 5; not a proof) |
| 8 | Three levels (two single U-copies of arbitrary weights and u₃ equal ones): at 12 rational weight triples (σ = 1), G_z[ρ,h] < 0, G_z[h,h′] < 0 across levels and (D) hold **for every u₃ ≥ 2** | **exact check** (pointwise in the weights, uniform in u₃; Section 5) |
| 9 | The symbolic reductions agree with direct exact elimination of the whole block (36 points for equal weights, u ≤ 12; 51 points for two levels; 36 for the three-level checks); two independent implementations of the invariant two-level reduction agree exactly | **exact check** |
| 10 | Negative controls: with β doubled or halved, G_z v has a negative entry at explicit rational points (u = 1, 3, 5); the symbolic (D) certificate at equal U-weights fails for these vectors | **exact** |
| 11 | G_z has no closed form in finitely many moment sums of the U-weights, as one might hope: the minimiser behind G_z v has a genuine pair structure for u ≥ 4 | **numerical** observation with a structural explanation (Section 4); not a theorem |
| 12 | For u = 4, 5 and σ/S ∈ {0.05, 1, 20}, the normalised ρ-row margin (G_z v)_ρ/(G_z[ρ,ρ] v₀) is smallest at equal U-weights (over the whole weight simplex); (G_z v)_ρ itself is not concave in the weights | **numerical** (Section 5); suggests that for the ρ-row Theorem EU is the worst case |
| 13 | Three levels symbolically (u₃ as a parameter and all weights free) | **not completed**: the exact elimination slowed down about six-fold per pivot and was stopped (Section 6) |
| 14 | Conjecture M for general U-weights and u ≥ 4, in particular (u, m) = (4, 1) at k = 5 | **open** |

Nothing here says anything about Conjecture 12 in general or about the NPT (Werner distillability) problem. For general
weights the type I blocks of (R∞)_{k−1} with |T| = 1 and |U| ≥ 4 remain open.

## Notation and the statement

Notation as in `relative-anchors/README.md` (Lemmas B, N; Theorem C; Theorem U3₁) and
`relative-anchors/code/typeI_G_u3/NOTES.md` (Sections 1, 4, 4b). One copy t carries the letter 2 (m = |T| = 1), U is
the set of the other u copies, s_x = w_x², σ = s_t, S = s_U. Forms are in units of 2Y. G_z is the Schur complement on
the letter-2 states in the coordinates (ρ, z_h), ζ_h = w_t w_h z_h. **Conjecture M** (typeI_G_u3, Section 4b):
**(Z)** every off-diagonal entry of G_z is negative; **(D)** G_z v > 0 entrywise for
v = (4σ(σ + S), 3(2σ + S), …, 3(2σ + S)). (Z) and (D) give G ≻ 0, hence (Theorem C) the type I block of (R∞)_{k−1}
is ⪰ 0 with kernel exactly span{p ⊗ f_h : h ∈ U}.

## 1. The formula for G when m = 1

Let T = {t}. The kernel vectors p ⊗ f_h have P_{t,h}-coordinate w_t c_h, so deleting the coordinates x_t = (P_{t,h})_h
is a gauge. In this gauge the eliminated variables are X_U ∈ R^{U×U} and q, and their block is exactly the (C_u) block
of Theorem C (the closed form of `cf.py` at X_T = 0, ρ = ζ = 0), which is positive definite. For m = 1 there is no
letter-2/letter-2 B-term, and at x_t = 0 the A-part on the letter-2 states is ½ρ² + (3/2)ρ² + |ζ|². Hence, as forms in
(ρ, ζ),

    G = diag(2, I_u) − Bᵀ C⁻¹ B,

where B is the coupling of (ρ, ζ) to (X_U, q), i.e. the bilinear part of

    ρ tr X − 2q(2Kρ + w_t w_U·ζ)/(2K − σ) − 2 Σ_{r∈U} ⟨𝒯_{(t,2)}ψ_r, E_∅⁻¹ 𝒯_{(r,1)} y⟩,   y = ρ e_t + ζ,  K = σ + S/2,

with E_∅ = diag(s_x + (σ + S)/2) + ½ w wᵀ on the one-hole strings of all k copies. This is the case m = 1 of Section 1
of `typeI_G_u2/NOTES.md` (for m = 1 the gauge w_Tᵀ X_T = 0 is x_t = 0, and the ℓ_t-terms there cancel). In
z-coordinates the first term stays diagonal, so for i ≠ j

    G_z[i, j] = −⟨b_i, C⁻¹ b_j⟩,   b_i = the coupling column of the coordinate i (ρ or z_h).

So (Z) says that the couplings are pairwise positively correlated in the C⁻¹ metric. If C were a nonsingular M-matrix
in a basis in which all b_i are nonnegative, (Z) would follow at once. In the coordinate basis this fails: at u = 3,
24 of the 100 entries of the scaled (C_u) block change sign over 200 random weight vectors (`ex_sep_run.py`).

`gzexact.py` computes G_z exactly over Q at rational weights in this gauge (u = 8 in about 0.1 s) and can use the gauge
P_{r0,h} instead; the two agree exactly wherever compared (`check_points.py`).

## 2. Theorem EU: equal U-weights, every u

**Theorem EU.** Let m = 1, u ≥ 1, σ > 0 and w_h = ω for all h ∈ U, s = ω². Then G_z is a symmetric Z-matrix and
G_z v > 0. Hence G ≻ 0.

*Proof (computer-assisted, symbolic in u; `equal_u.py`, `check_equal_u.py`).*
1. **Reduction.** The gauge x_t = 0 and the closed form are invariant under the permutations of U, and C is positive
   definite, so the minimiser for given letter-2 data is unique and depends equivariantly on the data.
   - For invariant data (ρ, z_h = z) it is invariant: X = s(a I + b(J − I)) and q (three unknowns).
   - For standard data ζ = w_t ω z f with Σ_h f_h = 0 (u ≥ 3) it is X = s(α diag f + β f 1ᵀ + γ 1 fᵀ), q = 0: the
     standard representation occurs exactly three times in R^{U×U} (once on the diagonal, twice off it) and not in q.
2. **The restricted forms, symbolically in u.** `equal_u.py` evaluates the closed form (formula in its docstring) on
   these families. Each sum over r, r′ ∈ U is split into r = r′ (u terms) and r ≠ r′ (u(u − 1) terms); vectors are
   combinations of e_t, 1_U, e_r, e_r′, f with ⟨1, 1⟩ = u, ⟨1, e_r⟩ = 1, ⟨1, f⟩ = 0, ⟨e_r, f⟩ = f_r, and
   Σ_r f_r = 0, Σ_r f_r² = |f|², Σ_{r≠r′} f_r f_r′ = −|f|². No odd power of w_t or ω survives (checked). Elimination
   gives the 2 × 2 matrix Q_inv(σ, s, u) of G_z on (ρ, z) and λ_std(σ, s, u), with Q_inv[0,0] = G_z[ρ,ρ],
   Q_inv[0,1] = Σ_h G_z[ρ,h] = u G_z[ρ,h], Q_inv[1,1] = Σ_{h,h′} G_z[h,h′], λ_std = G_z[h,h] − G_z[h,h′] (h ≠ h′).
3. **Validation.** At three rational weight pairs for each u = 1, …, 12 these quantities equal the values from exact
   elimination of the whole block (`gzexact.py`), exactly (λ_std for u ≥ 3).
4. **Certificates.** (G_z v)_ρ = Q_inv[0,0]v₀ + Q_inv[0,1]v₁, u(G_z v)_h = Q_inv[0,1]v₀ + Q_inv[1,1]v₁ and
   −u G_z[h,h′] = λ_std − Q_inv[1,1]/u. After u = 1 + x (u = 3 + x for the last one and for λ_std), every irreducible
   factor of the numerator and of the denominator of −Q_inv[0,1], (G_z v)_ρ, u(G_z v)_h, λ_std − Q_inv[1,1]/u and
   λ_std is a polynomial in (σ, s, x) with only positive coefficients (largest: degree 32, 153 terms). So these are
   positive for σ, s > 0 and every real u ≥ 1 (resp. u ≥ 3).
5. For u = 1 the reduction of step 1 degenerates (J − I = 0), and for u = 2 the z–z entry is not covered by step 4;
   both cases are row 2 of the status table. ∎

Denominators: σ + S, S + 2s + σ and S² + sS + 4σS + 4sσ + 2σ² = (2 + g₀)(3σ + S)(2s + σ + S)/4 come from E_∅
(g₀ = Σ_x s_x/d⁰_x, d⁰_x = s_x + (σ + S)/2); one further factor, of degree 5 in (σ, s), comes from the elimination (not
identified further).

### 2b. The limit u → ∞ at equal U-weights (`equal_u_limit.py`)

typeI_G_u3 (Section 4, item 3) observed in floating point that at equal U-weights the smallest eigenvalue of the
(ρ, z) block in the rescaled frame (ζ_t = c z_t w_U/|w_U| + ζ′_t, c = √(1 + S/σ)) decreases with u (σ = 0.004, S = 1:
0.579 at u = 1, …, 0.448 at u = 12) and asked whether it has a positive limit. That block is M(u) = D Q_inv D with
D = diag(1, c/(w_t√S)) and s = S/u. From the exact formula: the listed values are reproduced (0.5788, 0.5302, 0.5045,
0.4887, 0.4780, 0.4703, 0.4644, 0.4598 for u = 1, …, 8 and 0.4484 for u = 12), the limit M(∞) = lim_{u→∞} M(u)
exists, and

    M(∞)₀₀  = S(3S⁷ + 32S⁶σ + 144S⁵σ² + 352S⁴σ³ + 508S³σ⁴ + 424S²σ⁵ + 184Sσ⁶ + 32σ⁷)
              / (2(S + σ)²(S² + 2Sσ + 2σ²)(S² + 4Sσ + 2σ²)²),
    det M(∞) = S(S² + 6Sσ + 4σ²)(3S⁴ + 14S³σ + 24S²σ² + 16Sσ³ + 4σ⁴)
              / (2(S + σ)(S² + 2Sσ + 2σ²)(S² + 4Sσ + 2σ²)²)

(M(∞)₁₁ and M(∞)₀₁² in `logs/run_equal_u_limit.log`). All are ratios of polynomials with positive coefficients, so
M(∞) ≻ 0 for all σ, S > 0; at σ = 0.004, S = 1 its smallest eigenvalue is 0.4221. Since M(u) ≻ 0 for every u
(Theorem EU) and M(u) → M(∞) ≻ 0, the smallest eigenvalue is bounded below uniformly in u at fixed (σ, S). As S → 0
the determinant vanishes like S, consistent with the absence of a uniform gap there.

## 3. Theorem TL: U-weights with at most two values, every u

**Theorem TL.** Let m = 1, U = U₁ ⊔ U₂ with |U_i| = u_i, s_h = s_i on U_i (s₁, s₂ > 0 arbitrary), and σ > 0. Then
(Z) and (D) hold. Hence G ≻ 0.

*Proof (computer-assisted, symbolic in u₁, u₂; `two_level_fast.py`, `check_two_level.py`, `check_two_level_std.py`).*
If u₁ = 0 or u₂ = 0 this is Theorem EU, and u₁ = u₂ = 1 is u = 2; so let u₁, u₂ ≥ 1, not both 1.
1. **Invariant part.** The closed form is invariant under S_{u₁} × S_{u₂}. For invariant data (ρ, z_h = z_i on U_i)
   the minimiser has x_hh = p_i (h ∈ U_i), x_rh = o_ij (r ∈ U_i, h ∈ U_j, r ≠ h) and q: seven unknowns (six if u₁ = 1,
   when o₁₁ does not exist). The result is the 3 × 3 matrix Q3 with Q3[0,0] = G_z[ρ,ρ], Q3[0,i] = Σ_{h∈U_i} G_z[ρ,h] and
   Q3[i,j] = Σ_{h∈U_i, h′∈U_j} G_z[h,h′]. Case `one`: u₁ = 1, u₂ ≥ 2; case `gen`: u₁, u₂ ≥ 2; the case u₁ ≥ 2, u₂ = 1
   is `one` with the levels exchanged.
2. **Certificates for the invariant part** (`check_two_level.py`; u_i = 2 + x_i): the numerators of Q3[0,1], Q3[0,2],
   Q3[1,2] have only negative coefficients, those of the three rows of Q3(v₀, v₁, v₁)ᵀ only positive ones, and the
   denominators only positive ones (case `one`: about 2 000 terms each, in (σ, s₁, s₂, x₂); case `gen`: 31 000 to
   36 000 terms, in (σ, s₁, s₂, x₁, x₂)). By symmetry this gives G_z[ρ,h] < 0 for all h, G_z[h,h′] < 0 for h, h′ in
   different levels, and (D): (G_z v)_ρ is the first row and u_i (G_z v)_h (h ∈ U_i) is row i.
3. **Standard part** (pairs inside one level). For data ζ = w_t ω₁ z f with f supported on U₁, Σ f = 0, the minimiser
   is x = α diag f + β f 1ᵀ + γ 1 fᵀ on U₁ × U₁, x_rh = δ f_r on U₁ × U₂, x_hr = ε f_r on U₂ × U₁, zero on U₂ × U₂, and
   q = 0: the standard representation of S_{u₁} occurs three times in the U₁ × U₁ block if u₁ ≥ 3 (twice if u₁ = 2;
   then γ is dropped), once in each mixed block, and nowhere else. The Schur complement is λ₁ = G_z[h,h] − G_z[h,h′]
   (h ≠ h′ in U₁), and −u₁ G_z[h,h′] = λ₁ − Q3[1,1]/u₁. `check_two_level_std.py` certifies λ₁ − Q3[1,1]/u₁ > 0 by
   one-signed coefficients in the regimes (i) u₁ = 3 + x₁, u₂ = 2 + x₂ (about 247 000 terms); (ii) u₁ = 3 + x₁, u₂ = 1
   (Q3 from case `one` with the levels exchanged); (iii) u₁ = 2, u₂ = 2 + x₂. The case u₁ = 2, u₂ = 1 is u = 3. Pairs
   inside U₂ follow by exchanging the levels. ∎

**Validation.** Q3 equals exact elimination of the whole block at 12 (case `one`, u₂ = 2, 3, 5, 7) and 15 (case `gen`)
rational points, and λ₁ at 24 points (u₁ = 2, …, 5). The invariant reduction was implemented twice, with SymPy
expression trees (`two_level.py`) and with sparse polynomial rings (`two_level_fast.py`); for case `one` the two agree
exactly (`cmp_impl.py`).

## 4. Why there is no closed form in moment sums

A natural plan is to solve the (C_u) block against the source of v in closed form, with entries built from sums
Σ_h f(s_h), and then to reduce (D) to a few moment parameters as in Theorem C. This does not work as stated.

- **The crossed term.** The (C_u) form contains (1/(2 + g)) Σ_{r,r′} ρ_r ρ_{r′} X_{rr′} X_{r′r} (Theorem C, (F1)),
  ρ_h = s_h/d_h. It couples X_rh with X_hr through a coefficient depending on both s_r and s_h. Inverting these 2 × 2
  blocks produces pair functions such as 1/((2 + g)² − 4ρ_r²ρ_h²). The aggregates that couple the pairs (Xw, Xᵀw, Xη,
  Xᵀη with η = w/d) then satisfy linear equations whose coefficients are sums Σ_h s_h φ(s_h)/((2 + g)² − 4ρ_r²ρ_h²)
  with a kernel that does not factor in (r, h). Such sums are not functions of finitely many moments.
- **Numerically** (`ex_sep_run.py`, `ex_dep.py`): for the source of v the off-diagonal part of the scaled minimiser is
  not of the form x_rh = a_r + b_h (best additive fit leaves residuals 1.6·10⁻⁴ to 9·10⁻⁴ for u = 4, 5, 6, far above
  roundoff), and the minimiser is not symmetric (asymmetry up to 9·10⁻³). For the direction (1, 1/σ, …, 1/σ) the
  minimiser is symmetric, but again not additively separable for u ≥ 4. At fixed (σ, s_h, S), G_z[ρ,ρ] moves in the
  fourth digit when the other U-weights are redistributed.
- **What survives.** When the U-weights take finitely many values, every pair function takes finitely many values (one
  per pair of levels), and the reduction is a finite exact system with the multiplicities as parameters. This is how
  Sections 2 and 3 work. The number of unknowns grows like the square of the number of levels, and already three
  levels did not finish symbolically (Section 6).

So a proof for general weights has to treat the pair structure by inequalities (as Lemma H does for the hard part in
Theorem C) or by a structural argument. One possible route, which we have only begun to test (row 12): Theorem TL
covers all configurations with at most two distinct U-weights, among them the configurations with all but one weight
equal that are extremal in "equal variable" arguments for symmetric inequalities. A principle of that kind for the
entries of G_z v (at fixed σ, S and u, with the boundary handled by induction on u) would reduce Conjecture M to
Theorems EU and TL; we have no such principle.

## 5. Exact checks at rational points, controls, and floating-point observations

**Rational points** (`check_points.py 10 28`, `logs/run_check_points.log`). For u = 1, …, 10 and 28 rational weight
vectors per u (w rational, so s = w² is rational), cycling through nine families: log-uniform weights; one U-weight
tiny (s down to about 10⁻⁸ times the others); all U-weights shrunk or enlarged (S/σ from about 10⁻⁸ to 10⁸); one
U-weight dominant; equal U-weights; σ tiny; two U-weights tiny; geometric spreads. At every point (Z) and (D) hold and
all LDLᵀ pivots of G_z are positive, in exact arithmetic; at every seventh point the gauge P_{r0,h} gives the same G_z.
Smallest normalised margin min_i (G_z v)_i/(G_ii v_i) over the 28 points:

| u | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| min margin | 0.25 | 7.8·10⁻¹¹ | 1.1·10⁻¹¹ | 3.1·10⁻¹² | 1.6·10⁻⁹ | 2.2·10⁻¹⁰ | 3.2·10⁻¹¹ | 1.4·10⁻¹¹ | 3.9·10⁻¹⁰ | 2.3·10⁻¹⁰ |

`where_min.py` replays the same points for u ≤ 8: for every u ≥ 2 the smallest margin sits in the row z_h of the
smallest U-weight, and it occurs when σ is small compared with S (families "all large" and "σ tiny"), with that
U-weight between 10⁻⁵ and 10⁻¹ of S; it is not the family "one tiny" at moderate S/σ. Exact arithmetic resolves these
margins. The u = 10 points take about 3 minutes in total.

**Three levels, uniformly in u₃** (`three_level_points.py 12`). For 12 rational triples (w₁, w₂, w₃) with σ = 1
(including w₁ = 10⁻³ and w₃ down to 7·10⁻⁵), the invariant reduction of Section 3 with three levels (two single copies
and u₃ copies of weight s₃) is carried out exactly over Q(u₃). It agrees with direct elimination at u₃ = 2, 3, 5, and
with u₃ = 2 + x the numerators and denominators of G_z[ρ,h], of G_z[h,h′] for h, h′ in different levels and of the
rows of G_z v have one-signed coefficients. So at these weight triples (Z) (except inside the u₃ equal copies) and (D)
hold for every u₃ ≥ 2.

**Negative controls** (`controls.py`). With β doubled, v′ = (4σ(σ + S), 6(2σ + S), …), G_z v′ has a negative entry at
w_t = 1, w_U = (1/10, …) and (1/10, 1/20, 1/30, …) for u = 1, 3, 5; with β halved, at w_U = (10, …) and
(10, 13, 16, …). The original v passes at all these points. At equal U-weights the numerators of the (D) rows for
these perturbed vectors have mixed coefficients, so the certificate of Section 2 fails for them, as it must.

**Floating-point observations** (no claim rests on them). `ex_sep_run.py`: Section 4 and the sign count of Section 1.
`ex_ev.py`: for u = 4, 5, S = 1 and σ ∈ {0.05, 1, 20}, the normalised ρ-row margin (G_z v)_ρ/(G_z[ρ,ρ] v₀), minimised
over the weight simplex by Nelder–Mead from 12 starts, attains its minimum at equal U-weights (u = 4: 0.2974, 0.4491,
0.2811; u = 5: 0.2952, 0.4435, 0.2801), and no point of the two-level families does better. The unnormalised
(G_z v)_ρ is not concave along random segments of the simplex (midpoint defects down to −6·10⁻²), so concavity is not
the mechanism. The rows h are not covered by this test; their margins are the very small ones of the table above.

## 6. What was tried and did not finish

- **(4, 1) by a closed form.** In the gauge x_t = 0 the exact elimination for u = 3 (`closed_t.py 3`, ten pivots in
  four variables) takes about 75 s, with pivot times growing about 1.5- to 2-fold (the last one 33 s). For u = 4 there
  are 17 pivots in five variables; with the stopped runs of `typeI_G_u3/logs/abandoned/` this is far beyond the time
  here, and it was not started. Nothing is claimed for general weights at u = 4.
- **Three levels symbolically** (`three_level_fast.py`; two single U-copies of arbitrary weights and u₃ equal ones,
  invariant part, eleven pivots). A first run (σ symbolic, diagonal unknowns first) was stopped during its fourth
  pivot, a second run (σ = 1, q first) during its sixth; in the second run the pivot times were 0.3, 1.9, 16.7 and
  117 s, about six-fold growth per pivot. Logs in `logs/abandoned/`. The checker `check_three_level.py` is ready but
  was never run. The pointwise version (`three_level_points.py`, Section 5) replaces it.
- **The first two-level implementation** (SymPy expression trees) spent more than 12 minutes building the form of case
  `gen` and was stopped; its standard part (`two_level_std.py`) did not finish within 400 s (logs in
  `logs/abandoned/`). The polynomial-ring version builds the forms in under a second; the slowest elimination (case
  `gen`) takes about two minutes.
- **An M-matrix structure before elimination.** The (C_u) block is not a Z-matrix in the coordinate basis (Section 1),
  so the closure of M-matrices under Schur complements does not apply directly.

## 7. Files

| File | Content |
|---|---|
| `run_all.sh`, `guard.sh` | reruns every step, one Python process at a time (each killed above 2.8 GB or after its time limit); PASS/FAIL per step and ALL PASS; logs in `logs/run_*.log`, summary in `logs/run_all.out` (about 20 minutes) |
| `cf.py`, `flintfield.py`, `closed_gen.py` | copied from `typeI_G_u3/` (closed form of the type I block; rational functions on python-flint; the scaled symbolic build) |
| `gzexact.py` | exact G_z at rational weights, gauge x_t = 0 (or P_{r0,h}) |
| `closed_t.py u`, `show_cf.py u` | exact closed form of G_z in the gauge x_t = 0 (u ≤ 3); display in factored form |
| `check_m_small.py` | Conjecture M for u = 1, 2, 3 by coefficient signs; closed form = direct elimination |
| `equal_u.py`, `check_equal_u.py`, `equal_u_limit.py` | Theorem EU (reduction, validation, certificates) and the limit u → ∞ |
| `two_level_fast.py`, `check_two_level.py`, `check_two_level_std.py` | Theorem TL |
| `two_level.py`, `cmp_impl.py` | second implementation (expression trees) of case `one` and the exact comparison; `two_level_std.py`: expression-tree standard part (too slow, not used) |
| `three_level_points.py` | three levels at rational weights, uniformly in u₃ (exact checks) |
| `three_level_fast.py`, `check_three_level.py` | symbolic three-level reduction (did not finish; not run by `run_all.sh`) |
| `check_points.py`, `where_min.py` | exact rational-point checks of (Z), (D), pivots and gauge independence, u ≤ 10; where the smallest margins occur |
| `controls.py` | negative controls |
| `ex_sep_run.py`, `ex_ev.py`, `ex_sep.py`, `ex_csign.py`, `ex_struct.py`, `ex_dep.py`, `ex_v1.py`, `ex_v1b.py` | floating-point exploration (Sections 1, 4, 5); no claim rests on them |
| `*.pkl` | closed forms and reduced matrices written by the scripts (regenerated by `run_all.sh`) |
| `logs/`, `logs/abandoned/` | logs of a complete rerun from a fresh copy (`logs/run_all.out`: ALL PASS); the stopped runs |
