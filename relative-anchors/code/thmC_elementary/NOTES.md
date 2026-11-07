# An elementary proof of Theorem C, with an explicit margin (Theorem C′)

Ansh Mishra and Aryan Senthilkumar, 6 November 2026. Nothing here is refereed.

Notation as in `relative-anchors/README.md`, Theorem C: a type I block of (R∞)_{k−1}, T the copies carrying the letter
2, U the free copies, σ = Σ_{t∈T} s_t, S = Σ_{h∈U} s_h, and F(X) the form on X ∈ R^{U×U} obtained from the block on
(X_U, q) by eliminating q (step 1 of the proof of Theorem C).

**Theorem C′.** For every k, every type I block of (R∞)_{k−1} and every positive weight vector,

    F(X) ≥ σ/(12(σ + S)) ‖X‖_F²   for every X ∈ R^{U×U}.

Since the q-pivot is 2 > 0, the block on (X_U, q) is positive definite. This replaces steps 2–4 of the proof of
Theorem C (the coefficient certificates) by identities and elementary inequalities, and improves its margin
9·10⁻⁵ σ/(σ + S). Along the way: after the soft free variables are completed, the hard part satisfies H_c ⪰ I/2 on
symmetric hard matrices (the constant 1/2 is sharp) and H_c ≻ 4I/9 on antisymmetric ones.

The symmetric part is checked exactly against our closed form of the block, the antisymmetric part numerically,
and every scalar step exactly (Section 5).

## 1. Normalisation and splitting

The block is unchanged when every squared weight is multiplied by the same constant, so normalise S = 1 and write
x_i = s_i = w_i² (i ∈ U), Σ x_i = 1, E q = Σ_i x_i q_i, t = σ, K = t + 1/2, f_i = 1/(K + x_i), F = E f,
α = 1/(2 + F), r_i = x_i f_i, B = E r = 1 − KF, D = 2K + 1, τ = t/(t + 1). Then σ/(12(σ + S)) = τ/12.

F is invariant under X ↦ Xᵀ, so it splits into a symmetric and an antisymmetric branch, orthogonal both for F and for
the Frobenius norm. In each branch X is the sum of a hard part (Yw = 0) and a soft part:

    symmetric:      X = Y + w_iw_j(ℓ_i + ℓ_j + p),   E ℓ = 0,     ‖X‖² = ‖Y‖² + 2Eℓ² + p²;
    antisymmetric:  X = Z + w_iw_j(d_i − d_j),       E d = 0,     ‖X‖² = ‖Z‖² + 2Ed².

## 2. The symmetric branch: an exact sum of squares

With m_L = 2t(x + 1/2 + 2αx)f and c_p = t − 2t²αF − t²/(2K²), the symmetric form is

    F(X) = A₀(ℓ, p) + H₊(Y) + 2⟨Y, R(ℓ, p)⟩,
    A₀ = E[m_L ℓ²] + 2t²α(Eℓf)² + c_p p² − 8t²α p Eℓf,
    H₊(Y) = ‖Y‖²/2 + Σ_i K f_i Y_ii² + α Σ_{i,j} r_ir_j Y_ij²,
    ⟨Y, R⟩ = Σ_i A_iY_ii + α⟨Y, S⟩,   A_i = 2t r_iℓ_i + p t² r_i/K,   S_ij = a_ia_j(t(ℓ_i + ℓ_j)/K + p t²/K²),

with a_i = w_i r_i. Completing the diagonal square and the weighted matrix square gives the identity

    F(X) − ‖Y‖²/2 = Σ_i K f_i (Y_ii + A_i/(K f_i))² + α Σ_{i,j} r_ir_j (Y_ij + S_ij/(r_ir_j))²
                    + (t/K) Eℓ² + (tα/K²) E[r(ℓ + p t/K)²] + R* p²,
    R* = t(t² + t + 1/2)/(2K²(t + 1)) + t²B/(4DK³(D − B)) > 0.

*Derivation.* With l = ℓ + p t/(2K) the two completed squares cost (4t²/K)E(r l²) and
(2αt²/K²)(B E(r l²) + (E r l)²). Use E ℓ = 0, E(rℓ) = −K E(ℓf) (since r = 1 − Kf), the identity
2α − t/K − αtB/K² = α/(2K²), and, for the p coefficient, α = K/(D − B) and 4K⁴ − tD(tD + 1) = 1/4. ∎

Consequences:
- Minimising over (ℓ, p) gives H_c,+ ⪰ I/2 on the symmetric hard space. The constant is sharp: at four equal free
  weights the hard matrix with rows (0, 1, −1, 0), (1, 0, 0, −1), (−1, 0, 0, 1), (0, −1, 1, 0) has no soft coupling and
  quotient 1/2 + αr² → 1/2 as t → ∞.
- At Y = 0: A₀ ≥ (t/K)Eℓ² + R*p² ≥ (τ/2)(2Eℓ² + p²), since t/K ≥ τ and R* ≥ τ/2.

## 3. The antisymmetric branch

Here q, ℓ and p decouple. The soft form is A_D(d) = E[m_D d²] + 2t²α(E d f)² with m_D = 2t(x + 1/2 − 2αx)f ≥ tf, so
A_D ≥ tE(fd²) ≥ tEd²/(K + 1) ≥ (τ/3)·2Ed². On antisymmetric hard Z,

    H₋(Z) = ½‖Z‖² − α⟨Z, RZR⟩,   R = Π diag(r) Π on w^⊥.

1. **A covariance bound (Jung's theorem).** For every probability vector z, the k largest eigenvalues of
   diag(z) − zzᵀ sum to at most k/(k + 1): for a rank-k projection P the sum Tr P(diag(z) − zzᵀ) is the variance of
   points Pe_i (weights z_i) of diameter at most √2 in a k-dimensional space, and such a set lies in a ball of squared
   radius k/(k + 1). For the eigenvalues λ₁ ≥ λ₂ ≥ … of R this gives Σ_{i≤k} λ_i ≤ k/((k + 1)K + 1) (compression of the
   operator concave function s/(K + s), whose nonzero eigenvalues of Π diag(x) Π are those of diag(x) − xxᵀ) and
   Σ_{i≤k} λ_i ≤ kF/(k + 1) (comparison of the secular equations of the compressions of diag(r) orthogonally to √x and
   to √(r/F)).
2. **Hard bound.** The negative part of H₋ has largest eigenvalue αλ₁λ₂ ≤ s²/(4(2 + 3s/2)) with s = λ₁ + λ₂ ≤ 2/(3K + 1)
   and F ≥ 3s/2. The bound is increasing in s, so H₋ ⪰ h(t) I with h(t) = 1/2 − 1/((3K + 1)(6K + 5)).
3. **The d completion.** The half-source of Z to d is J_i = 2tαw_if_i(Zv_f)_i with v_f = (w_i(f_i − F))_i. Since
   A_D ≥ tE(fd²), eliminating d costs at most E(J²/(tf)) = 4tα² Σ_i x_i²f_i(Zv_f)_i² ≤ δ(t)‖Z‖² with
   δ(t) = t/(8K²(K + 1)³), using x²f ≤ 1/(K + 1) on [0, 1], α ≤ 1/2, ‖v_f‖² = Var(f) ≤ 1/(4K²(K + 1)²) and
   ‖Z‖_op² ≤ ‖Z‖²/2 for antisymmetric Z.
4. Hence H_c,− ⪰ (h − δ)I, and (h − δ − 4/9)·144K²(K + 1)³(3K + 1)(6K + 5) =
   8(t − 1)² + 11/2 + t/4 + 1171t² + 4099t³ + 5296t⁴ + 3388t⁵ + 1104t⁶ + 144t⁷ > 0. So H_c,− ≻ 4I/9. ∎

## 4. The margin

**Lemma (normalised shell).** If a form x*Ax + 2Re x*By + y*Hy has A, H ≻ 0 and H − B*A⁻¹B ⪰ cH, then it is at least
(1 − √(1 − c))(x*Ax + y*Hy). (‖A^{−1/2}BH^{−1/2}‖² ≤ 1 − c, then Cauchy–Schwarz.)

- Symmetric branch: H₊ ⪯ (3/2)I (diagonal coefficients 3/2 − r_i + αr_i² ≤ 3/2, off-diagonal ones 1/2 + αr_ir_j ≤ 3/2),
  so H_c,+ ⪰ I/2 ⪰ H₊/3 and the lemma with c = 1/3 gives F ≥ (1 − √(2/3))(A₀ + H₊) ≥ (1/6)((τ/2)‖X_soft‖² + ½‖Y‖²)
  ≥ (τ/12)‖X‖².
- Antisymmetric branch: H₋ ⪯ I/2 and H_c,− ⪰ 4I/9 ⪰ (8/9)H₋, so c = 8/9 and F ≥ (2/3)(A_D + H₋)
  ≥ (2/3)((τ/3)‖X_soft‖² + (4/9)‖Z‖²) ≥ (2τ/9)‖X‖².

Together F(X) ≥ (τ/12)‖X‖² = σ/(12(σ + S)) ‖X‖². ∎

## 5. Checks

| File | What it checks | Log |
|---|---|---|
| `symmetric_identity.py` | exact (python-flint): the formula of Section 2 equals our q-eliminated form F (built from `../typeI_G_u3/cf.py` by `freeblock.py`) on symmetric X, and the sum-of-squares identity holds, at 90 cases (ten rational unit weight vectors with u = 2, …, 5, three values of t, three random rational X each) | `logs/symmetric_identity.log` |
| `antisymmetric_bound.py` | floating point: the hard Schur operator of the antisymmetric branch at 200 configurations (u = 3, …, 7, t from e⁻⁶ to e⁶) has λ_min ≥ h − δ (closest gap 3·10⁻⁸, the sharp case) and λ_min − 4/9 ≥ 0.009 | `logs/antisymmetric_bound.log` |
| `margin.py` | floating point: λ_min(F) ≥ 6 · σ/(12(σ + S)) at 150 unnormalised weight vectors (the factor 6 is approached as σ/S → ∞, where the symmetric constant 1/2 is sharp), and F(X) = F(Xᵀ) | `logs/margin.log` |
| `scalars.py` | exact (sympy): every scalar identity and inequality of Sections 2–4 | `logs/scalars.log` |

Run `bash run_all.sh` (Python with NumPy, SymPy and python-flint).
