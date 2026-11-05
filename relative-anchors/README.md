# Relative anchors: the group conditions for every number of copies

Ansh Mishra and Aryan Senthilkumar, 4 November 2026

In `single-excitation-range/` the endpoint inequality on two classes of rank-two operators was reduced to finitely many
matrix inequalities. The classes are operators whose range lies in the single excitations over a product vector, and
operators whose range lies in the tangent space of a product vector; every co-range and every number of copies n is
allowed. The reduction is Lemmas 11, 11′, 13 and 13′ there. Pick a single excitation p_G = Σ_{t<k} w_t |1_t⟩ in the
range plane, with all w_t > 0 on a group of k copies. The inequalities are the group conditions (R∞)_j, j = −1, …, k, on
those k copies alone (Conjecture 12 there is equivalent to them for every k). They are proved there for k ≤ 3.

This folder collects what we can prove about the group conditions for every k.

**Results.**
- **Every k, every positive weight vector.**
  - Closed forms for all entries (Lemmas 1, B, N, G).
  - The two extreme conditions hold: (R∞)_{−1} (Proposition M) and (R∞)_k, the latter with kernel exactly p ⊗ H_k and
    an explicit gap (Theorem K).
  - (R∞)_{k−1} holds on its letter-1 block (Theorem T), on every block in which one copy carries the letter 2 in row
    and column (Proposition T2), on the block in which every copy carries one letter 2 (Proposition U0), and on every
    block in which all copies but one carry one letter 2 (Proposition U1). On every remaining type I block, the part
    on the copies without a letter 2 is positive definite (Theorem C), so each such block is equivalent to a single
    Schur-complement condition G ≻ 0 on its letter-2 states. That condition holds when exactly two copies are free
    of the letter 2 (Theorem U2), and for k = 4 also when three are (Theorem U3₁). So (R∞)_{k−1} holds completely
    for k ≤ 4, and for k ≥ 5 it is open only on the type I blocks with at least three such copies.
  - (R∞)_{k−2} has the explicit kernel vector vac ⊗ J⁻|1…1⟩ (Proposition K2). On its letter-1 block the vacuum row
    eliminates exactly and the soft directions are known exactly; the gap there is at most min_h s_h/(2(|s| − s_h))
    (Section 3e), so a proof has to be exact to first order in the smallest weight.
  - (R∞)_0 holds completely for k ≥ 3, with kernel exactly p ⊗ H_0 (Theorem B0, an algebraic proof).
  - In every (R∞)_j, every block in which one copy carries the letter 2 in row and column is positive definite, with
    an explicit lower bound (Proposition T2′).
- **Four copies, every positive weight vector** (Section 3f, exact computer-assisted certificates). (R∞)_1 holds
  (Theorem A4), and (R∞)_2 holds on every block except the four type I blocks in which exactly one copy carries the
  letter 2 (Theorem B4); this includes the letter-1 block of Section 3e. With the results above, all of (R∞) at four
  copies is reduced to those four blocks, which are images of each other under permutations of the copies.
- **Equal weights** (w_t = 1).
  - The letter-1 block of every (R∞)_j is ⪰ 0 for every k, with the kernel described exactly (Theorem E1). The proof is
    computer-assisted, symbolic in k.
  - The blocks with one copy carrying the letter 2 in row and column are positive definite for every k
    (Proposition E2).
  - The blocks in which some copies carry the letter 2 once are ⪰ 0 for every k (Theorem E3, whose last family is
    Proposition C, an exact determinant certificate). The proofs are computer-assisted, symbolic in k and |T|.
  - Hence all of (R∞) holds at equal weights, for every k.

**Corollary.** Conjecture 12 of `single-excitation-range/` holds for every single excitation Σ_{t∈G} x_t (each x_t a
single-copy excitation, all ‖x_t‖ equal), for every finite set G of copies. Hence, for every n:
- every sector of every single-excitation range plane that contains such a vector is nonnegative;
- q_n ≥ 0 on every tangent plane that contains such a vector, with any co-range.

(Local unitaries fixing |0⟩ bring such a vector to Σ_{t∈G} |1_t⟩.)

**Open.**
- For general weights: the type I blocks of (R∞)_{k−1} with at least three copies free of the letter 2 and k ≥ 5
  (|U| ≥ 3, |T| ≥ 1; the case k = 4 is Theorem U3₁), reduced to the condition G ≻ 0 of Theorem C; and the letter-1
  and type I blocks of (R∞)_j for 1 ≤ j ≤ k − 2 (for k ≥ 4). At k = 4 only the type I blocks of (R∞)_2 with one
  copy carrying the letter 2 remain (Section 3f).

Nothing here is refereed. The endpoint problem itself is open.

## 1. Setting

Notation is as in `single-excitation-range/README.md`. The group consists of copies 0, …, k − 1. Put s_t = w_t² and
|s| = Σ_t s_t. One copy carries M_1 = I − ½|Ω⟩⟨Ω| on C³_row ⊗ C³_col. Let A = (M_1^{⊗k})^Γ and B = R_{p_G}^Γ with
R_p = M V_p D_p^{−1} V_p^* M, and put Y = A − B.
- G1_j is the set of states with row in H_1 and column in H_j.
- G0_j is the set of states with row equal to the vacuum and column in H_j.
- (R∞)_j is Y compressed to G1_j ⊕ G0_{j+1}.

Y preserves the number of letters 1 in row and column together, and at each copy the number of letters 2. The torus
blocks are labelled by these numbers.
- The *letter-1 block* is the one with no letter 2.
- In a *type I* block the letter 2 occurs exactly once (in row or column) at each copy of a set T ≠ ∅.
- In a *type II* block some copy carries the letter 2 in both row and column.

**Operators on column strings.** Z = k − j is the number of zeros of a string in H_j.
- For a single excitation a = (letter ℓ at copy r): σ_a⁻ = |0⟩⟨ℓ|_r. Also n_r = 1 − |0⟩⟨0|_r.
- J⁺ = Σ_t w_t |1⟩⟨0|_t and J⁻ = (J⁺)^T.
- C = Σ_t s_t (2|0⟩⟨0|_t + |2⟩⟨2|_t) and E_c = C + ½J⁺J⁻ on H_c.
- 𝒯_a = ½J⁺σ_a⁻ + [ℓ = 1] w_r (2|0⟩⟨0|_r + |2⟩⟨2|_r) − [ℓ = 2] w_r |1⟩⟨2|_r.
- N_r = (3/2)|0⟩⟨0|_r + ½|1⟩⟨1|_r + |2⟩⟨2|_r and 𝒩 = Σ_t s_t N_t.

## 2. Closed forms

**Lemma 1.** The compression D_c of D_{p_G} to H_c satisfies 2^{k−c} D_c = E_c.

*Proof.* D_c[x, x′] = Σ_{t,t′} w_t w_{t′} M[(e_t, x), (e_{t′}, x′)], and M factors over the copies.
- **t = t′.** The entry is diagonal. Copy t contributes 1 − ½[x_t = 1]. Every other copy contributes 1, or ½ if
  x_s = 0.
- **t ≠ t′.** Copy t must have x_t = 1, x′_t = 0 and contributes −½. Copy t′ must have x_{t′} = 0, x′_{t′} = 1 and
  contributes −½. The other copies agree and contribute ½ per zero.

Multiplying by 2^{k−c}, the diagonal is Σ_t s_t (2[x_t = 0] + ½[x_t = 1] + [x_t = 2]), which is the diagonal of
C + ½J⁺J⁻. The off-diagonal entries are ½ w_t w_{t′} for moving a letter 1 from t′ to t, which are those of ½J⁺J⁻. ∎

So D_c ⪰ 2^{−(k−c)−1}|s| ≻ 0. This also follows from the rank-one bound.

**Lemma B (the entries of Y).** Let a = (ℓ at r) and b = (ℓ′ at r′) be single excitations. Then:
- 2^Z Y[(a, y′), (b, y)] = ⟨y′| (½ − [r = r′]) σ_b⁺σ_a⁻ + [a = b](2 − n_r) − 𝒯_b^T E_j^{−1} 𝒯_a |y⟩ for y, y′ ∈ H_j;
- 2^Z Y[(a, y′), (vac, y)] = −⟨y′| σ_a⁻ − J⁻ E_{j+1}^{−1} 𝒯_a |y⟩ for y′ ∈ H_j, y ∈ H_{j+1};
- 2^Z Y[(vac, y′), (vac, y)] = 2⟨y′| I − ½ J⁻ E_{j+2}^{−1} J⁺ |y⟩ for y, y′ ∈ H_{j+1} (no J-term if j + 2 > k).

*Proof.*
- **The A-part.** A[(a, y′), (b, y)] = M[(b, y′), (a, y)] is a product over the copies. Copy by copy it equals the stated
  boson expression: the swap term ½σ_b⁺σ_a⁻, minus the same at a common copy, plus the diagonal 2 − n_r.
- **The B-part.** B[(a, y′), (b, y)] = ⟨G_{(b,y′)}, D^{−1} G_{(a,y)}⟩ with G_{(a,y)}[z] = Σ_t w_t M[(a, y), (e_t, z)].
- **G on G1 states.** Computing G copy by copy gives three terms:
  - the column unchanged, with weight w_r 2^{−Z+[y_r = 0]} when ℓ = 1;
  - the letter at r replaced by 1, with weight −½ w_r;
  - the letter at r moved to an empty copy t as a 1, with weight +½ w_t.

  The last two combine to ½J⁺σ_a⁻ plus a diagonal remainder, so G_{(a,·)} = 2^{−Z} 𝒯_a on H_j.
- **G on G0 states.** G_{(vac,y)} = −2^{−Z(y)} J⁺|y⟩.
- With Lemma 1 the powers of two combine as stated. The exchange of a and b between the two G factors is the partial
  transpose. ∎

As a quadratic form, write Ψ = Σ_a e_a ⊗ ψ_a + vac ⊗ φ. Then

    2^Z ⟨Ψ, YΨ⟩ = Σ_{a,b} (½ − [r_a = r_b]) ⟨σ_b⁻ψ_a, σ_a⁻ψ_b⟩ + Σ_a ⟨ψ_a, (2 − n_{r_a})ψ_a⟩ − Σ_{a,b} ⟨𝒯_b ψ_a, E_j^{−1} 𝒯_a ψ_b⟩
                  − 2 Σ_a ⟨ψ_a, (σ_a⁻ − J⁻E_{j+1}^{−1}𝒯_a) φ⟩ + 2 ⟨φ, (I − ½J⁻E_{j+2}^{−1}J⁺) φ⟩.

The crossed pairs ⟨K_b ψ_a, L K_a ψ_b⟩ are the partial transpose. Such a sum equals ⟨V_sym, L V_sym⟩ − ⟨V_anti, L V_anti⟩
for V = Σ_{a,b} e_a ⊗ e_b ⊗ K_b ψ_a.

**Lemma N (normal order).** On every H_c,
- 𝒯_{(r,1)} = ½σ_{(r,1)}⁻J⁺ + w_r N_r;
- 𝒯_{(t,2)} = ½σ_{(t,2)}⁻J⁺ − ½w_t|1⟩⟨2|_t;
- E_c = 𝒩 + ½J⁻J⁺.

So on ker J⁺ every 𝒯_a acts locally, and Σ_r w_r 𝒯_{(r,1)} = E_c. The last identity is the kernel p ⊗ H_j of
(R∞)_j.

*Proof.* Use [J⁺, σ_{(r,1)}⁻] = w_r(|1⟩⟨1|_r − |0⟩⟨0|_r), [J⁺, σ_{(t,2)}⁻] = w_t|1⟩⟨2|_t and
J⁺J⁻ − J⁻J⁺ = Σ_t s_t(|1⟩⟨1|_t − |0⟩⟨0|_t). ∎

**Lemma G (the vacuum-row block).** On G0_{j+1}, Y = 2^{−(k−j−1)}(I − ½J⁻E_{j+2}^{−1}J⁺). If C is invertible on H_{j+2},
this equals 2^{−(k−j−1)}(I + ½J⁻C^{−1}J⁺)^{−1}, by the Woodbury identity with E = C + XX^*, X = J⁺/√2. ∎

## 3. The condition (R∞)_{k−1}, every k

**Theorem T.** For every k ≥ 2 and all w_t > 0, the letter-1 block of (R∞)_{k−1} is ⪰ 0. Its kernel is exactly
{p ⊗ f : f ∈ span(f_h)} ⊕ {u ⊗ J⁻|1…1⟩ : u ∈ span(e_r), u ⊥ p}, of dimension 2k − 1.

The block consists of P_{r,h} = (e_r, f_h), where f_h = 1…1 − e_h has its hole at h (r, h ∈ [k]), and of
Q = (vac, 1…1).

*Proof.*
1. **Closed form.** By Lemma B with Z = 1, unrescaled:
   - A[P_{r,h}, P_{r′,h′}] = ¼[(r,h) = (r′,h′), r ≠ h] + [r = h, r′ = h′](¾[r = r′] + ¼).
   - A[Q,Q] = 1 and A[Q, P_{r,h}] = −½[r = h].
   - B[P_{r,h}, P_{r′,h′}] = m_{r′,h}^T D^{−1} m_{r,h′}, B[Q, P_{r,h}] = −½w_r w_h/|s| and B[Q, Q] = 0.

   Here D = D_{k−1} = ½(diag(d) + ½ww^T) in the basis f_h, with d_h = s_h + ½|s|, and
   m_{r,h} = ¼(w_r f_h + w_h f_r) + ½[r = h] w_h f_h.
2. **Eliminate Q.** Write Ψ = Σ X[r,h] P_{r,h} + qQ. Since Y[Q,Q] = 1, the minimum over q is
   ⟨X, Y_PP X⟩ − c(X)², with c(X) = −½ tr X + ½ w^TXw/|s|.
   - The vectors X = wc^T (that is, p ⊗ f) lie in the kernel, so we may assume w^TX = 0. Then c(X) = −½ tr X.
   - The A-part minus c² is ¼‖X‖² + ½Σ_a X_aa².
3. **The B-term.** Expand the products of the m's and use w^TX = 0. With v = Xw and Z = X diag(w), the B-term is
   𝔅(X) = Σ_{a,b} D^{−1}_{ab}(v_a v_b/16 + Z_ab v_b/8 + v_a Z_ba/8 + Z_ab Z_ba/4).

   Here D^{−1} = 2 diag(1/d) − (2/(2 + G))(w/d)(w/d)^T with G = Σ_a s_a/d_a. So the reduced form is
   F(X) = ¼‖X‖² + ½ΣX_aa² − 𝔅(X).
4. **The second radical.** Split X = X′ + a w^T, with a = Xw/|s|. Then X′w = 0, w^TX′ = 0 and w^Ta = 0.
   - All cross terms vanish, using D^{−1}_{ab} d_a = 2δ_ab − (2/(2 + G)) w_a w_b/d_b.
   - F(a w^T) = ¼|s||a|² + ½Σ_a s_a a_a² − ½Σ_a d_a a_a² = 0.
   - So F(X) = F(X′).
5. **Positivity.** Now v = 0, and with ρ_a = s_a/d_a ≤ 2/3,
   F(X′) = ¼Σ_{a≠b} X′_ab² + Σ_a(¾ − ½ρ_a)X′_aa² + (1/(2(2 + G))) Σ_{a,b} ρ_a ρ_b X′_ab X′_ba.

   For a ≠ b put x = s_a/|s| and y = s_b/|s|. Then ρ_a = 2x/(1 + 2x), ρ_b = 2y/(1 + 2y), x + y ≤ 1, and
   (1 + 2x)(1 + 2y)(2 + ρ_a + ρ_b − 12ρ_aρ_b) = 2 + 6(x + y) − 32xy ≥ 2(1 − x − y)(1 + 4(x + y)) ≥ 0,
   using xy ≤ (x + y)²/4. Since G ≥ ρ_a + ρ_b, this gives ρ_aρ_b/(2(2 + G)) ≤ 1/24. Bound each crossed pair by
   2|X′_ab X′_ba| ≤ X′_ab² + X′_ba², and use ρ_a ≤ 2/3 on the diagonal:
   F(X′) ≥ (5/24)Σ_{a≠b} X′_ab² + (5/12)Σ_a X′_aa² ≥ (5/24)‖X′‖². ∎

   (The first version had 5/36 here. As a form on all X, F has smallest eigenvalue exactly 5/24 in the limit of two
   weights ½ and the rest 0, so the crossed-pair bound cannot give more.)

**Proposition T2.** Every type II block of (R∞)_{k−1} is positive definite.

*Proof.*
- **Closed form.** All states have row e_t^{(2)}. By Lemma B, Y = ¼I − ½𝒯^T E^{−1} 𝒯 with 𝒯 = 𝒯_{(t,2)}.
- **The action of 𝒯.** On the states y_h it gives 𝒯y_h = ½(w_h u_t − w_t u_h), with E = diag(d′) + ½w_{U′}w_{U′}^T on
  the one-hole strings, d′_x = s_x + s_{T′} + ½s_{U′}.
- **The estimate.** For ψ = Σ c_h y_h put g = (c·w_U)e_t − w_t c. Cauchy–Schwarz and E ⪰ diag(d′) give
  ⟨g, E^{−1}g⟩ ≤ 2|c|²(a/(a + 3b) + b/(a + b)), with a = s_U and b = s_t.
- **Conclusion.** The bracket equals 1 − 2ab/((a + 3b)(a + b)) < 1, so Y ≻ 0. ∎

**Proposition U0.** The type I block of (R∞)_{k−1} in which every copy carries one letter 2 is positive definite.

*Proof.*
- **Closed form** (Lemmas B and N; units of 2Y; S = |s|).
  - Diagonal entries are 2 on the states R_t and 2 at Q.
  - Between R_t and R_{t′} the entry is ½ + ¼κ(s_t, s_{t′}), with
    κ(a, b) = (S(a + b) + a² + b²)/(S² + S(a + b)/2 + ab − (a² + b²)/2).
  - Between R_t and Q the entry is −c_t, with c_t = 2S/(2S − s_t).
- **Eliminate Q.** This gives F with diagonal D_t = 2(S − s_t)(3S − s_t)/(2S − s_t)² and off-diagonal entries
  β(s_t, s_{t′}) < 0, where, with c = S − a − b,
  β(a, b) = −ab(4a² + 7ab + 7ac + 4b² + 7bc + 3c²)/[(a + 2b + 2c)(2a + b + 2c)(2a² + 8ab + 5ac + 2b² + 5bc + 2c²)].
- **Diagonal dominance.** The function 2b(3S − a)/(2S − a)² + β(a, b) is b times a ratio of polynomials with
  nonnegative coefficients and positive numerator. Summing over t′ gives Σ_{t′≠t}|F_tt′| < D_t, so F is strictly
  diagonally dominant with positive diagonal. ∎

**Remark (the other type I blocks: an exact reduction).** Let T ≠ ∅ be the copies carrying one letter 2, with
U = [k] ∖ T ≠ ∅, m = |T| and u = |U|. The block has u(k + m) + m + 1 states:
- P_{r,h}: row letter 1 at r ∈ [k]; column 2 on T and 1 on U, with a hole at h ∈ U;
- R_t: row letter 2 at t ∈ T; column with a hole at t, 2 on T ∖ t and 1 on U;
- Z_{t,h}: row letter 2 at t ∈ T; column 1 at t, 2 on T ∖ t and 1 on U, with a hole at h ∈ U;
- Q: vacuum row; column 2 on T and 1 on U.

On the one-hole strings of U put E_0 = diag_{h∈U}(s_h + ½s_U) + ½w_Uw_U^T and E = s_T I + E_0. The principal block on
the states P_{t,h} with t ∈ T is

    K = ½ I − ½ w_Tw_T^T ⊗ E^{−1},   K^{−1} = 2I + 2 w_Tw_T^T ⊗ E_0^{−1},   K ⪰ s_U/(2(2s_T + s_U)) · I.

Indeed K acts as ½E_0E^{−1} on the direction of w_T and as ½ on its complement, and E_0 ⪰ ½s_U. So these mu
coordinates can be eliminated exactly. The u kernel vectors p ⊗ f_h survive the elimination and are removed by
deleting one coordinate P_{r_0,h} (r_0 ∈ U fixed) for each h. The block is then ⪰ 0 iff an explicit matrix of size
u(k − 1) + m + 1 is ⪰ 0. We have not proved that it is in general; for u = 1 it is Proposition U1.

**Proposition U1.** For every k ≥ 2 and every positive weight vector, every type I block of (R∞)_{k−1} with |U| = 1
(all copies but one carry one letter 2) is ⪰ 0, with kernel exactly p ⊗ f_0.

*Proof (computer-assisted, `code/typeI_u1/`, exact arithmetic).* Let U = {0}, σ = s_0, T = {1, …, m} (k = m + 1) and
S_T = Σ_{t∈T} s_t.
1. **Reduction.** Put the forced kernel coordinate x_0 to zero and eliminate the P_{t,0} (t ∈ T) with the explicit K^{−1}
   of the Remark. What is left is an explicit form Φ(ρ, ζ, q) in the coordinates of R_t, Z_{t,0} and Q, and the block
   is ⪰ 0 with kernel exactly p ⊗ f_0 iff Φ ≻ 0. Its q-entry is Φ_qq = (3S_T² + 2S_Tσ + σ²)/(2(S_T + σ/2)²) > 0;
   eliminating q leaves a 2m × 2m matrix G with 2 × 2 blocks G_tt, G_tt′ in the coordinates (ρ_t, ζ_t).
2. **Closed forms.** Every block depends only on (s_t, s_t′, σ, S_T): the other copies enter Φ only through S_T. The
   blocks are explicit rational functions (homogeneous of degree 0), checked as identities of rational functions for
   m ≤ 4 with S_T free, and exactly at 25 weight vectors with m ≤ 7.
3. **Rescaling and a weighted Gershgorin bound.** Rescale ζ_t by c = √(1 + σ/S_T) and call the result Ĝ. If
   λ_min(Ĝ_tt) > Σ_{t′≠t} ‖Ĝ_tt′‖ √(s_t′/s_t) for every t, then Ĝ ≻ 0 (weighted block Gershgorin with d_t = √s_t). This
   follows from two inequalities with H = 2/5:
   - (A) ‖Ĝ_tt′‖_F ≤ H √(s_t s_t′)/S_T for t ≠ t′;
   - (B) λ_min(Ĝ_tt) > H (1 − s_t/S_T).
   Indeed, (A) gives Σ_{t′≠t} ‖Ĝ_tt′‖ √(s_t′/s_t) ≤ (H/S_T) Σ_{t′≠t} s_t′ = H(1 − s_t/S_T).
4. **Certificates.** After normalising S_T = 1, (A) and (B) are polynomial inequalities in (s_t, s_t′, σ) and (s_t, σ)
   (degrees up to 22). On compactified domains (the triangle s_t + s_t′ ≤ 1 mapped to the unit square; σ ∈ [0, 1] and
   σ = 1/x for σ ≥ 1), all tensor Bernstein coefficients are positive, except at one vertex (s_t = S_T, σ = 0) that lies
   outside the domain of (B). ∎

Numerically sup ‖Ĝ_tt′‖_F S_T/√(s_t s_t′) ≈ 0.283 and inf λ_min(Ĝ_tt)/(1 − s_t/S_T) = 3/4, so H = 2/5 has room on
both sides; the certificate fails, as it must, for H = 7/25 and for H = 4/5.

**Theorem C (the copies without a letter 2, every u).** Let k ≥ 2, let T ≠ ∅ be the copies carrying one letter 2 in a
type I block of (R∞)_{k−1}, U = [k] ∖ T, u = |U| ≥ 1, and let the weights be positive. Then the principal submatrix on
the states P_{r,h} (r, h ∈ U) and Q is positive definite. It depends on T only through σ = Σ_{t∈T} s_t. Consequently
the block is ⪰ 0 with kernel exactly span{p ⊗ f_h : h ∈ U} if and only if G ≻ 0, where G is the Schur complement on
the states R_t, Z_{t,h} (t ∈ T) after eliminating P_{r,h} (all r) and Q modulo the kernel. For u = 1, G ≻ 0 is
Proposition U1; for u ≥ 2 it is open.

*Proof (computer-assisted in the last step; `code/typeI_u2/`, `NOTES.md` there).*
1. **Structure.** In the closed form of the block (derived from Lemmas B and N, checked exactly against the builder),
   the states P_{t,h} (t ∈ T) couple to the rest only through ξ = Σ_t w_t X_{t·}, and the part with w_Tᵀ X_T = 0
   decouples with identity block. So the statement reduces to the block on (X_U, q), where X ∈ R^{U×U} holds the
   coordinates of P_{r,h}, r, h ∈ U. Eliminating q leaves a form F(X). It involves T only through σ, with
   E = diag(s_h + S/2 + σ) + ½w_U w_Uᵀ and S = Σ_{h∈U} s_h.
2. **Splitting.** Write X = X′ + a wᵀ + w cᵀ with X′w = 0 = wᵀX′ and a ⊥ w. Then:
   - F(X′) is the form of step 5 of Theorem T, with every s_h + S/2 replaced by s_h + S/2 + σ;
   - F(awᵀ + wcᵀ) = σ f(a, c), with f explicit;
   - the cross term is σ⟨𝒢, X′⟩, with 𝒢 explicit.

   Every coupling is a multiple of σ.
3. **Hard part.** The crossed-pair bound of Theorem T gives
   F(X′) ≥ (5/12) Σ_{r≠h} X′_rh² + Σ_h (3/2 − ρ_h) X′_hh², with ρ_h = s_h/(s_h + S/2 + σ).
   Minimising entry by entry gives F(X) ≥ σ Φ1(a, c) for an explicit quadratic form Φ1 on R^{2u}.
4. **Φ1 ≻ 0.** With t = σ/(σ + S/2):
   - The a-block is ⪰ (1 − t)I.
   - The off-diagonal entries of the c-block are ≤ 0, so a weighted Gershgorin bound with weights √s_h applies.
   - A Schur test with parameter 4/5 bounds the cross block.

   Exact identities such as Σ_{k≠h} s_kρ_k = S′ − Σ_{k≠h} ρ_k bound every sum over the other copies through three
   parameters on [0, 1]. This reduces everything to three rational inequalities, (I), (II), (III), on [0, 1]³. All
   their tensor Bernstein coefficients are positive; (III) needs the corner rule at one vertex. Negative controls
   with the constants moved (4/5 → 7/10 or 2, and (III) shifted by 1/100) fail as they should. ∎

The same certificates with margins give F(X) ≥ 9·10⁻⁵ (σ/|s|)‖X‖², uniformly in u; numerically the constant is about
0.42. The weighted block Gershgorin route of Proposition U1 for G ≻ 0 is numerically sound for u = 2…5 with the same
extreme constants (5/18 and 3/4) as for u = 1, but its blocks have no closed form yet.

**Theorem U2 (two copies free of the letter 2).** For every k ≥ 3, every type I block of (R∞)_{k−1} with |U| = 2 is
⪰ 0 with kernel exactly span{p ⊗ f_h : h ∈ U}, for all positive weights.

*Proof (computer-assisted; `code/typeI_G_u2/`, `NOTES.md` there).* By Theorem C it remains to show G ≻ 0.
1. **Structure, by hand, for every u and m.** Work in the gauge w_Tᵀ X_T = 0, which is a complement of the kernel.
   Eliminating X_T with that constraint, and then X_U and Q through the block of Theorem C, gives the following.
   G_tt depends only on (s_t, σ, w_U), and G_tt′ (t ≠ t′) only on (s_t, s_t′, σ, w_U). So closed forms computed for
   one configuration with three letter-2 copies give the blocks for every m. They were also compared exactly with
   direct elimination for m = 1, …, 5 (115 blocks).
2. **Closed forms for u = 2.** Exact symbolic elimination over Q(w_0, …, w_4) (python-flint). In the coordinates
   ζ_{t,h} = w_tw_h z_{t,h} every entry is a rational function of the squared weights.
3. **Weighted block Gershgorin, as in Proposition U1** (same frame and rescaling, H = 2/5):
   - (B_2) λ_min(Ĝ_tt) > (2/5)(1 − s_t/σ). Its determinant is positive by Bernstein certificates in two blow-up charts
     around the degenerate edge (m = 1, S → 0); the matrix is positive definite at one rational point; the inertia is
     constant on the connected domain.
   - (A_2) ‖Ĝ_tt′‖_F < (2/5)√(s_t s_t′)/σ. A Bernstein certificate on [0, 1]⁴, after mapping the triangle onto the
     square.

   Every denominator factor is certified sign-definite. With H = 1/4 for (A) or H = 4/5 for (B), exact rational
   witnesses show the statements false. ∎

For |U| ≥ 3 the same pipeline applies in principle for each fixed u, but G then depends on w_U through more than
one scalar functional, and no closed form uniform in u is known.

**Theorem U3₁ (k = 4).** For k = 4, every type I block of (R∞)₃ with |U| = 3 (so |T| = 1) is ⪰ 0 with kernel exactly
span{p ⊗ f_h : h ∈ U}, for all positive weights. Together with Theorem T, Propositions T2, U0, U1 and Theorems C and
U2, the whole condition (R∞)_{k−1} holds for k ≤ 4 and all positive weights.

*Proof (computer-assisted, exact; `code/typeI_G_u3/`, `NOTES.md` there).* By Theorem C it remains to show G ≻ 0, a
4 × 4 matrix on the letter-2 states (ρ, ζ_{t,h}), h ∈ U.
1. **Coordinates.** Put X_{r,h} = w_rw_h x_{r,h} and ζ_{t,h} = w_tw_h z_h. Modulo a common part, the states carry the
   parity vectors e_r + e_h (P_{r,h}), e_t + e_h (Z_{t,h}) and 0 (R_t, Q) in (Z/2)^k, so every entry of the rescaled
   block is a rational function of the squared weights s_i (checked: an odd exponent would raise an error). The Schur
   complement in these coordinates is G_z = DGD with D diagonal and invertible, so G ≻ 0 ⟺ G_z ≻ 0.
2. **Closed form.** Exact Gaussian elimination over Q(s_0, …, s_3) (python-flint) of the ten states other than the
   gauge rows and the letter-2 states. The eliminated block is positive definite on the open orthant (Theorem C), so
   the result is the Schur complement at every positive weight vector. It agrees with direct elimination over Q at
   12 rational weight vectors, in two different gauges.
3. **Certificate.** G_z = Ñ/den with den = 8 f₇ f₁₇ (s_0 + S)² e₄², where f₇, f₁₇, e₄ are irreducible of degrees 7,
   17, 4 and S = s_U. Every factor has only positive coefficients. Every leading principal minor of G_z is a product
   of powers of these factors times a core polynomial (degrees 34 to 56, up to 23 374 terms), and every core has only
   positive coefficients. Sylvester's criterion gives G_z ≻ 0 on the open orthant.
4. **A second certificate.** G_z is a symmetric Z-matrix: every off-diagonal numerator has only negative
   coefficients. With v = (4s_t(s_t + S), 3(2s_t + S), 3(2s_t + S), 3(2s_t + S)), every entry of den·G_zv has only
   positive coefficients. So diag(v)G_z diag(v) has nonpositive off-diagonal entries and positive row sums, and it
   is strictly diagonally dominant. ∎

Negative controls: G_z − W/10 (W the identity of the rescaled frame of Theorem U2) and G_z − e_ρe_ρᵀ/1000 are
indefinite at explicit rational points, so there is no room near S → 0.

**Remarks.**
- The same minor certificate proves G ≻ 0 for (u, m) = (1, 1), (2, 1), (1, 2), (2, 2) and (1, 3), independently of
  Proposition U1 and Theorem U2. We conjecture that it works for every (u, m): every leading principal minor of G_z,
  in the order (ρ_t, z_{t,·}) copy by copy, is a ratio of polynomials with only positive coefficients (Conjecture P in
  the notes). This would give all of (R∞)_{k−1}.
- For m = 1 the Z-matrix certificate with the same vector v holds numerically for u ≤ 8 (Conjecture M in the notes).
  A proof would settle every type I block of (R∞)_{k−1} with exactly one copy carrying the letter 2, for every k. For
  m ≥ 2 the Z-pattern fails.
- The comparison with the one-copy system obtained by merging the copies of U into one copy of weight √S, which
  Proposition U1 suggests, is not one-sided: the difference is indefinite, also at equal U-weights (notes,
  Section 4). So the constants of Proposition U1 do not transfer directly.
- For k = 5 the same pipeline ((u, m) = (4, 1) and (3, 2)) did not finish in our runs.

**Proposition K2.** The vector vac ⊗ J⁻|1…1⟩ lies in the kernel of (R∞)_{k−2}.

*Proof.*
- **Vacuum row.** E_k = ½|s| on 1…1 and J⁺J⁻|1…1⟩ = |s||1…1⟩, so (I − ½J⁻E_k^{−1}J⁺)J⁻|1…1⟩ = 0.
- **Rows (r, 1).** By Lemma N, 𝒯_{(r,1)}J⁻|1…1⟩ = E_{k−1}σ_r⁻|1…1⟩. So the coupling σ_r⁻ − J⁻E_{k−1}^{−1}𝒯_{(r,1)}
  annihilates J⁻|1…1⟩.
- **Rows (t, 2).** These give 0, since there is no letter 2. ∎

Theorem T and Proposition K2 account for every kernel vector beyond p ⊗ H_j that our numerics found at j = k − 1 and
j = k − 2.

## 3b. The condition (R∞)_0, every k

**Theorem B0.** For every k ≥ 3 and every positive weight vector, the condition (R∞)_0 holds, and its kernel is exactly
p ⊗ H_0. (For k = 2 the condition is part of Theorem 16 of `single-excitation-range/`; there the form (C0) below
degenerates, since δ_12 = 0.)

Normalise |s| = 1 and write x_t = s_t and λ_t = 1/(1 − x_t). Put γ = Σ_t x_t/(2(1 − x_t)) = ½Σ_t(λ_t − 1) ≥ ½. The
blocks of (R∞)_0 are:
- the letter-1 block, with states (1@r, vac) and (vac, 1@r);
- for each r, a 2 × 2 letter-2 block with states (2@r, vac) and (vac, 2@r).

**Closed form** (Lemmas B and G at j = 0, scaled by 2^{k−1}). Let |w⟩ = J⁺|vac⟩ = Σ_t w_t|1@t⟩ and u = E_1^{−1}w.
Sherman–Morrison gives u_r = w_r/(2(1 − x_r)(1 + γ/2)).
- The row block is I − ww^T on letter-1 rows and I on letter-2 rows.
- The column block is I − ½J⁻E_2^{−1}J⁺.
- The cross block is −½(I − K), with K[a, y] = ⟨u|𝒯_a|y⟩.

**Letter-2 blocks.** The block at r is [[1, −½(1 − κ_r)], [−½(1 − κ_r), η_r]] with
- η_r = 1/(1 + γ_r/2), where γ_r = Σ_{t≠r} x_t/(2 − x_r − 2x_t);
- 1 − κ_r = (1 + g_r)/(1 + γ/2), where g_r = x_r/(2(1 − x_r)).

It is ⪰ 0 iff (1 + γ/2)² ≥ ¼(1 + g_r)²(1 + γ_r/2).

*Proof.* Put a = 1 + g_r and G′ = Σ_{t≠r} g_t, so that 1 + γ/2 = ½(a + 1 + G′). With e = 1 − x_r and y = x_t ∈ [0, e],

    (1 + e)y/(2e(1 + e − 2y)) − y/(1 − y) − y/(2e) = −y(e − y)(1 + e − y)/(e(1 − y)(1 + e − 2y)) ≤ 0.

The first term is a·x_t/(2 − x_r − 2x_t) and the second is 2g_t. Summing over t ≠ r, where Σ_{t≠r} x_t = e, gives
aγ_r ≤ 2G′ + ½. Hence (a + 1 + G′)² − a²(1 + γ_r/2) ≥ (a + 1 + G′)² − a² − a(G′ + ¼) = a(7/4 + G′) + (1 + G′)² > 0, so
every letter-2 block is positive definite. ∎

**Letter-1 block.**
- **The cross block factors.** With τ = 1/(1 + γ/2), Λ = diag(λ_r) and P = I − ww^T, it equals −½τPΛ. This uses
  (I − K)^T w = 0, so (w, 0), which is p ⊗ H_0, is in the kernel.
- **Reduction to (C0).** Eliminate the rows (Lemma G gives the column block as (I + ½N)^{−1} with N = J⁻C^{−1}J⁺). The
  block is ⪰ 0 iff v^TΛ(I + ½N)Λv ≤ 4(1 + γ/2)²|v|² for all v ⊥ w.
- **The form of (C0) we prove.** With v_r = w_r ζ_r, δ_rt = 1 − x_r − x_t and W_rt = x_r x_t/δ_rt, it reads: for all
  ζ with Σ_r x_r ζ_r = 0 and S = Σ_r x_r ζ_r²,

      L := Σ_r x_r λ_r² ζ_r² + ¼ Σ_{r<t} W_rt (λ_r ζ_r + λ_t ζ_t)² ≤ (2 + γ)² S.               (C0)

*Proof of (C0).* We prove the stronger bound L ≤ (2 + 2γ)S; it also holds for complex ζ. Since γ ≥ ½ and
(2 + 2γ)/(2 + γ)² decreases for γ ≥ 0, it gives L ≤ (12/25)(2 + γ)²S. Write L = D + T, with D the first sum and T the
second, put d_rt = ζ_r − ζ_t and δ_rt = 1 − x_r − x_t, and let

    E = Σ_{r<t} x_r x_t (λ_r + λ_t)|d_rt|²,   G = Σ_{r<t} x_r x_t (Σ_{u∉{r,t}} x_u λ_u²) |d_rt|².

- **E − D = G.** The constraint gives, for each r, the weighted variance identity
  λ_r Σ_{t≠r} x_t|d_rt|² − λ_r²|ζ_r|² = λ_r² Σ_{t<u; t,u≠r} x_t x_u |d_tu|². Multiply by x_r and sum over r.
- **T ≤ G.** The constraint also gives λ_rζ_r + λ_tζ_t = Σ_{u∉{r,t}} x_u (λ_t d_ru + λ_r d_tu). Cauchy–Schwarz with total
  weight δ_rt, and then |a + b|² ≤ 2|a|² + 2|b|², give
  T ≤ ½ Σ_{r<t} Σ_{u∉{r,t}} x_r x_t x_u (λ_t²|d_ru|² + λ_r²|d_tu|²) = G, counting each triple with its three differences.
- **E ≤ (2 + 2γ)S.** λ_r + λ_t = 2 + x_rλ_r + x_tλ_t ≤ 2 + 2γ, and Σ_{r<t} x_r x_t |d_rt|² = S.

So L = D + T ≤ D + G = E ≤ (2 + 2γ)S. The three gaps add up to an exact identity on the constraint subspace:

    (2 + 2γ)S − L = Σ_{r<t} Σ_{u∉{r,t}} x_r x_t x_u λ_u |d_rt|² + ¼ Σ_{r<t} Σ_{u∉{r,t}} x_r x_t x_u |λ_t d_ru − λ_r d_tu|²
                    + ¼ Σ_{r<t} (x_r x_t/δ_rt) Σ_{u<v; u,v∉{r,t}} x_u x_v (λ_r + λ_t)² |d_uv|²,

with every coefficient nonnegative. ∎

**Kernel.** In units of 2^{k−1} the letter-1 block is [[P, X], [X^*, B]] with P = I − ww^T (|w| = 1), X = −½τPΛ and
B = (I + ½N)^{−1} ≻ 0. The bound gives XB^{−1}X^* ⪯ (12/25)P, so eliminating the columns leaves
P − XB^{−1}X^* ⪰ (13/25)P on the rows. Hence the kernel of the block is exactly (w, 0), that is p ⊗ H_0 (note X^*w = 0).
The letter-2 blocks are positive definite.

The worst case of (C0) itself is the ratio 4/9, approached as x → (½, ½, 0, …). (The first version of this proof used a
regime analysis and interval arithmetic, `code/c0_interval.py`; we keep it as an independent check.)

## 3c. The conditions (R∞)_{−1} and (R∞)_k, every k

**Proposition M.** For every k ≥ 2 and every positive weight vector, (R∞)_{−1} holds (for k = 1 the entry is 0).

*Proof.* It is the 1 × 1 block (vac, vac). By Lemma G it equals 2^{−k}(1 − ½w^TE_1^{−1}w), where on the letter-1 part
of H_1, E_1 = 2 diag(|s| − s_r) + ½ww^T. Sherman–Morrison gives w^TE_1^{−1}w = γ/(1 + γ/2) with
γ = Σ_r s_r/(2(|s| − s_r)), so the block is 2^{−k}/(1 + γ/2) > 0. ∎

**Theorem K.** For every k ≥ 2 and every positive weight vector, (R∞)_k ⪰ 0 and its kernel is exactly p ⊗ H_k. On the
orthogonal complement of p ⊗ H_k in G1_k,

    Y ⪰ ν(w) I,   ν(w) = min{11/48, ½(1 − max_t s_t/|s|)}.

At j = k there is no vacuum row and Z = 0, so Lemma B gives Y itself. Every column string lies in {1, 2}^k. For a set
T of copies let y_T be the string with the letter 2 on T and 1 elsewhere, and s_T = Σ_{t∈T} s_t.

*Closed form.* On H_k we have J⁺ = 0 and n_r = 1, so by Lemma N:
- E_k is diagonal, E_k|y_T⟩ = ½(|s| + s_T)|y_T⟩;
- 𝒯_{(r,1)} = w_r N_r is diagonal, where N_r = 1 if the column has the letter 2 at r and ½ if it has the letter 1;
- 𝒯_{(t,2)}|y⟩ = −½ w_t [y_t = 2] |y with t → 1⟩;
- σ_b⁺σ_a⁻ vanishes on H_k unless a and b sit at the same copy.

So the torus blocks are:
- for each T ⊆ [k] (∅ included) a type I block on P_r = (1@r, y_T), r ∈ [k], and Z_t = (2@t, y_{T∖t}), t ∈ T;
- for each B and t ∈ B a type II block with the single state (2@t, y_B).

Together they contain all 2k·2^k states. Put N_r = 1 for r ∈ T and ½ otherwise, v_r = w_rN_r, a = Σ_r s_rN_r =
½(|s| + s_T), a_t = a − s_t/2 and a_tu = a − (s_t + s_u)/2. Lemma B gives:
- on the P_r: P = diag(N) − vv^T/a;
- between P_r and Z_t: C_rt = (a/(2a_t))(w_tv_r/a − δ_rt);
- on the Z_t: B_tt = 1 and B_tu = −w_tw_u/(4a_tu) for t ≠ u;
- a type II block: (a_B − s_t)/(2a_B − s_t), with a_B = ½(|s| + s_B).

*Proof.*
1. **The P-block.** Since Σ_r N_rs_r = a, x^TPx = min_λ Σ_r N_r(x_r − λw_r)². So P ⪰ ½(I − ww^T/|s|), and the
   kernel of P is spanned by w.
2. **Completing the square.** Let L be the k × |T| matrix with columns L_t = −(a/(2a_t))e_t. Then C = PL, so the form
   on (x, z) is (x + Lz)^TP(x + Lz) + z^TFz with F = B − L^TPL. Its entries are
   F_tt = 1 − a(a − s_t)/(4a_t²) and F_tu = −w_tw_us_ts_u/(16a_ta_ua_tu).
3. **F ⪰ (11/16)I.** a(a − s_t) ≤ a_t², so F_tt ≥ ¾. Put x_t = s_t/a. Then Σ_{t∈T} x_t ≤ 1, and a_t, a_u, a_tu ≥ a/2,
   so |F_tu| ≤ ½(x_tx_u)^{3/2}. Hence Σ_{u≠t}|F_tu| ≤ ½x_t^{3/2}(Σ_{u≠t} x_u)^{3/2} ≤ ½[x_t(1 − x_t)]^{3/2} ≤ 1/16.
4. **Type I blocks.** By steps 1–3 the form is ⪰ 0, and it vanishes only at z = 0, x ∈ span(w), which is p ⊗ y_T.
   For the gap, let x_⊥ be the part of x orthogonal to w and y = x_⊥ + Π_{w⊥}Lz. Since s_t ≤ a we have ‖L‖ ≤ 1, so
   ‖x_⊥‖² + ‖z‖² ≤ 2‖y‖² + 3‖z‖². The form is at least ½‖y‖² + (11/16)‖z‖², hence at least (11/48)(‖x_⊥‖² + ‖z‖²).
5. **Type II blocks.** (a_B − s_t)/(2a_B − s_t) = ½(1 − s_t/(|s| + s_B − s_t)) ≥ ½(1 − s_t/|s|), since s_B ≥ s_t. ∎

In our samples the smallest type I eigenvalue off the kernel is about 0.44, so 11/48 is not sharp. This fits the
searches, where (R∞)_k never came close to binding.

## 3d. The type II blocks of every (R∞)_j

**Proposition T2′.** For every k ≥ 2, every j and every positive weight vector, every torus block of (R∞)_j in which a
copy t carries the letter 2 in row and column is positive definite. More precisely, let u = |2@t⟩ be the row of such a
block and x_t = s_t/|s|. The compression of Y to u ⊗ (all columns) satisfies

    Y|_{u ⊗ ·} ⪰ (4/9) 4^{−k} min{x_t, 1 − x_t} · I.

*Proof.*
1. **Positivity.** u is a product vector, so by Theorem 1 of `support-and-projector-notes/` the form q_k is ≥ 0 on every
   operator whose range lies in span{u, p}, whatever the co-range. So the compression of M^{⊗k} to span{u, p} ⊗ (all
   columns) is ⪰ 0. In the basis (u, p) write it as [[D_u, C], [C^*, D_p]]. Its Schur complement D_u − CD_p^{−1}C^* is the
   compression of M^{⊗k} − R_p to the row u. Since u is a real product vector, this equals the compression of
   Y = (M^{⊗k} − R_p)^Γ to the row u. A type II block has the single row u, so it is a principal submatrix.
2. **One-copy facts.** For a one-copy computational vector v, Y_v has kernel K_v = (v ⊗ C³) ⊕ Sym²(v^⊥) (for v = |0⟩,
   Y_v = |10⟩⟨10| + |20⟩⟨20| + (|12⟩ − |21⟩)(⟨12| − ⟨21|), row ⊗ column), and T_v has spectrum {1, −1/3}. By
   tensorisation every eigenvalue of T_u other than 1 is at most 1/9, so Y_u ⪰ (8/9)2^{−k} Π_{K^⊥} with K = ⊗_i K_{u_i}.
3. **The other Schur complement.** At the copy t (where u has the letter 2), K_{|2⟩}^⊥ is spanned by |02⟩, |12⟩ and
   (|01⟩ − |10⟩)/√2, and the row of p at copy t has density (|s| − s_t)|0⟩⟨0| + s_t|1⟩⟨1|. Hence, for every column
   vector z, ⟨p ⊗ z, Y_u(p ⊗ z)⟩ ≥ ε‖z‖² with ε = (4/9)2^{−k} min{s_t, |s| − s_t}. By the product identity this says
   D_p − C^*D_u^{−1}C ⪰ εI.
4. **Exchange.** Put K = D_u^{−1/2}CD_p^{−1/2}. Since D_p ⪯ |s| (‖M^{⊗k}‖ = 1), step 3 gives I − K^*K ⪰ (ε/|s|)I, hence
   I − KK^* ⪰ (ε/|s|)I, hence D_u − CD_p^{−1}C^* ⪰ (ε/|s|)D_u ⪰ (ε/|s|)2^{−k}I. ∎

This contains Propositions T2 (j = k − 1) and E2 (equal weights). Numerically the bound holds with a factor of at least
9 to spare (`code/check_typeII_allj.py`).

**What blocks the rest.** Two parts are still open:
- the type I blocks of (R∞)_{k−1} with at least three copies free of the letter 2: by Theorem C each is equivalent to
  G ≻ 0 on its letter-2 states (proved for two such copies, Theorem U2, and for k = 4, Theorem U3₁);
- the letter-1 and type I blocks of the conditions (R∞)_j with 1 ≤ j ≤ k − 2 (at k = 4 all but one orbit of
  blocks are certified, Section 3f).

In both, the method of Theorem T leaves soft directions: eigenvalues 0.01–0.1, and as small as 0.0015 for
(R∞)_{k−2} at spread weights. They lie in the J⁻ part of the G1 columns and tend to 0 at the boundary of the weight
simplex. So no uniform crossed-pair bound can work, and these directions have to be treated exactly. Section 3e
makes this precise for the letter-1 block of (R∞)_{k−2}.

## 3e. The letter-1 block of (R∞)_{k−2}: an exact reduction and a sharpness bound

The block has the states (1@r, holes P) with |P| = 2 and (vac, hole h). Write Ψ = Σ_r e_r ⊗ ψ_r + vac ⊗ φ and
X[r, P] = ψ_r(P). Proofs and checks are in `code/middle_k2/` (`NOTES.md` there); here are the statements.

- **Lemma R1 (the vacuum row).** 4⟨Ψ, YΨ⟩ = F(X) + 2|φ_⊥ − ½c(X)|², where F = Y11 − ½Y10Y01 is explicit, c(X) ⊥ w,
  and φ_⊥ is the part of φ orthogonal to w. So the block is ⪰ 0 with kernel exactly
  p ⊗ H_{k−2} ⊕ span(vac ⊗ J⁻|1…1⟩) iff F ≻ 0 on {Σ_r w_r ψ_r = 0}.
- **Lemma W (the only nonlocal object).** Let 𝒩|P⟩ = (½|s| + s_P)|P⟩ and M = I + ½J⁺𝒩^{−1}J⁻, a k × k matrix with
  M_hg = ½w_hw_g/(½|s| + s_h + s_g) for h ≠ g. Then E_{k−2}^{−1} = 𝒩^{−1} − ½𝒩^{−1}J⁻M^{−1}J⁺𝒩^{−1} and
  E_{k−2}^{−1}J⁻ = 𝒩^{−1}J⁻M^{−1}.
- **Proposition S (soft directions).** Fix a copy h and a ∈ R^k with a_h = 0 and a ⊥ w. Put
  Ψ(a, h) = Σ_r a_r e_r ⊗ J⁻|h⟩ and u = Σ_t a_t|{h, t}⟩. Then Ψ(a, h) ⊥ p ⊗ H_{k−2}, |Ψ(a, h)|² = (|s| − s_h)|a|², and

      F(Ψ(a, h)) = 2s_h|a|² − 4s_h²⟨u, E_{k−2}^{−1}u⟩ − ½|c|²,
      c = 4s_h [½⟨E_{k−1}^{−1}e_h, a⟩ w + (E_{k−1}^{−1}e_h) ∘ a ∘ w].

- **Corollary S (sharpness).** For every k ≥ 3 and all weights, the smallest eigenvalue of the block on the orthogonal
  complement of its kernel is at most min_h s_h/(2(|s| − s_h)).

*Proof of Corollary S.* Take v = (Ψ(a, h), ½c), which is orthogonal to the kernel because c ⊥ w. By Lemma R1 and
Proposition S, ⟨v, Yv⟩ = F(Ψ(a, h))/4 ≤ s_h|a|²/2, while |v|² ≥ (|s| − s_h)|a|². ∎

So the gap tends to 0 linearly when any weight does, and a proof has to be exact to first order in the smallest weight
along the directions Ψ(a, h). Numerically (k = 4, …, 7, 610 weight vectors) the block is ⪰ 0 with exactly the stated
kernel, its gap lies between 0.37 and 1 times the bound of Corollary S, and F ⪰ 0.45 on the complement of the soft
directions. As s_h → 0 the k − 2 lowest eigenvalues are (s_h/(2|s|))κ_i + o(s_h), with κ_i ∈ [0.90, 0.97] in our
examples. Positivity of the corresponding (k − 2) × (k − 2) first-order matrix, which depends only on the other k − 1
copies, is a necessary condition; we have not proved it.

## 3f. The middle conditions at four copies

At k = 4 the conditions (R∞)_{−1}, (R∞)_0, (R∞)_3 and (R∞)_4 hold for all weights (Proposition M, Theorem B0, Theorem
U3₁ with the results of Section 3, Theorem K), and the type II blocks of every (R∞)_j are positive definite
(Proposition T2′). What remains are the letter-1 and type I blocks of (R∞)_1 and (R∞)_2. Code, logs and notes are in
`code/middle_k4/`.

**Reductions.**
- *Permutations.* A permutation of the copies, applied to rows and columns at once and with the weights permuted
  along, maps every torus block of (R∞)_j to a torus block of (R∞)_j. So a certificate symbolic in all four weights,
  for one block of each orbit, covers the whole orbit. The representatives are T = ∅ (the letter-1 block), {0},
  {0, 1} and, for j = 2, {0, 1, 2}.
- *Coordinates.* In the rescaled basis of `code/ylib.py` (each state divided by ∏_t w_t^{N_t}, N_t the number of
  letters 1 at copy t in row and column) every entry is a rational function of the squared weights s_t. It is
  homogeneous of degree j + 1 − |T|, so s_0 = 1 loses nothing.
- *Gauge.* The known kernel vectors are p ⊗ y for every column y of the block (in these coordinates, the all-ones
  vector on the states (1@r, y)), and, for the letter-1 block of (R∞)_2, vac ⊗ J⁻|1111⟩ (Proposition K2). Delete the
  coordinates (1@1, y), and for that block also (vac, hole 0). Each known kernel vector is the only one that is
  nonzero on its deleted coordinate. So the block is ⪰ 0 with kernel exactly the span of the known vectors iff the
  remaining principal submatrix is positive definite.

**Certificate.** Let N = den·M, where M is the gauge-fixed block over Q(s_1, s_2, s_3) and den is the lcm of the
denominators of its entries.
- Fraction-free (Bareiss) elimination of N over Z[s_1, s_2, s_3] produces every leading principal minor exactly.
- The irreducible factors of den all have only positive coefficients. They are divided out exactly after every step.
- Then D_r(M) = core_r · ∏ f_i^{x_ir − rE_i}/cst^r.

In every certified block, every core has only positive coefficients. So every leading principal minor is positive on
the open orthant, and Sylvester's criterion applies. This is exact on the whole open orthant, as Corollary S
requires: the minors may tend to 0 at its boundary.

Checks:
- every D_r is compared with a direct determinant at a random integer point;
- the certified minors agree with minors of the block built directly over Q, at four rational points per block;
- the symbolic blocks agree exactly with the point builder, and with the dense construction from M, V_p and D_p in
  floating point;
- with B replaced by (1 + 1/1000)B, every block fails, with an explicit rational witness.

| j | T | blocks | size → gauge-fixed | kernel (exactly) | cores: largest degree, most terms | status |
|---|---|---|---|---|---|---|
| 1 | ∅ | 1 | 22 → 18 | p ⊗ \|1@c⟩ | 93, 45 335 | proved |
| 1 | {t} | 4 | 11 → 10 | p ⊗ \|2@t⟩ | 62, 27 536 | proved |
| 1 | {t, t′} | 6 | 3 | none | 22, 2 296 | proved |
| 2 | ∅ | 1 | 28 → 21 | p ⊗ (letter-1 part of H_2) ⊕ span(vac ⊗ J⁻\|1111⟩) | 121, 47 905 | proved |
| 2 | {t} | 4 | 21 → 18 | expected p ⊗ \|2@t 1@c⟩ | | open (positive definite after the gauge at 33 rational weight vectors, exact) |
| 2 | {t, t′} | 6 | 12 → 11 | p ⊗ \|2@t 2@t′⟩ | 79, 59 548 | proved |
| 2 | {t, t′, t″} | 4 | 4 | none | 31, 5 980 | proved |

**Theorem A4 ((R∞)_1 at four copies).** For k = 4 and every positive weight vector, (R∞)_1 holds, and its kernel is
exactly p ⊗ H_1.

**Theorem B4.** For k = 4 and every positive weight vector, every block of (R∞)_2 except the type I blocks with
|T| = 1 is ⪰ 0 with the kernel listed in the table. In particular, the letter-1 block of Section 3e holds at four
copies, with kernel exactly p ⊗ (letter-1 part of H_2) ⊕ span(vac ⊗ J⁻|1111⟩).

*Proofs (computer-assisted, exact).* Type II blocks by Proposition T2′; the others by the reductions and the
certificate above, one representative per orbit. ∎

If the four remaining blocks are ⪰ 0, then all of (R∞) holds at four copies. By Lemma 13 and the remark on (R2) in
`single-excitation-range/`, Conjecture 12 then holds for every single excitation supported on at most four copies.
We do not claim this yet.

## 4. Equal weights, every k

At w_t = 1 the conditions commute with S_k.

**Reduction.**
- **The representation.** On the letter-1 block of (R∞)_j it is C^k ⊗ M^{(k−j,j)} ⊕ M^{(k−j−1,j+1)}, where M is the
  permutation module of j-subsets.
- **Which irreps occur.** C^k ⊗ M^{(k−j,j)} is the permutation module of pairs (r, Y). Its two orbits are Young modules
  of types (k−j, j−1, 1) and (k−j−1, j, 1). So every irrep λ that occurs has at most three rows, with λ_3 ≤ 1.
- **The twisted subgroups.** For 0 ≤ i ≤ k/2 let H_i = (S_2 ≀ S_i) × S_{k−2i} (i disjoint pairs and the rest), with
  the character χ_i that is the sign on each pair swap.
- **Pieri.** Restricted to such λ, Ind_{H_i}^{S_k} χ_i contains exactly (k−i, i) and (k−i−1, i, 1), each once.
- **Consequence.** The block is ⪰ 0 iff, for every i, the form Q_i on the χ_i-twisted invariants is ⪰ 0. These are
  spanned by at most seven explicit vectors: products of singlets on the pairs with Dicke vectors on the rest, and
  the row excitation in the rest or in a pair.

**Computing Q_i without enumeration.** The entries come from Lemma B.
- Sums over (a, b) reduce to seven orbits of H_i with explicit sizes.
- J^± act on Dicke vectors by explicit factors.
- E_c^{−1} is applied by Newton interpolation on the Johnson-scheme eigenvalues
  λ_{c,i′} = (k − c − i′)(c + 1 − i′) of J⁻J⁺. Only the nodes i′ ∈ [i − 2, i + 2] can occur, and only genuine nodes
  are used.
- So each entry is a rational function of q0 = j − i, t = k − j − i and i.

**Theorem E1.** At equal weights, for every k ≥ 2 and every j = −1, …, k, the letter-1 block of (R∞)_j is ⪰ 0. Its
kernel is exactly:
- p ⊗ (letter-1 part of H_j);
- plus the vector of Proposition K2 when j = k − 2;
- plus the k − 1 vectors of Theorem T when j = k − 1.

*Proof (computer-assisted, `code/prove_families.py`, a few seconds).*
1. **Families.** The parameter space splits into 42 families, according to q0 ∈ {−1, 0, 1, ≥ 2},
   t ∈ {−1, …, 3, ≥ 4} and i ∈ {0, 1, ≥ 2}. Here ≥ means the parameter stays symbolic.
2. **The computation.** In each family Q_i is computed over Q(q0, t, i). The script checks Q_i κ = 0 for the kernel
   vector κ and runs symmetric Gaussian elimination.
3. **Positivity of the pivots.** Shift to nonnegative variables (q0 = 2 + x, t = 4 + y, i = 2 + z). Then every pivot is
   a ratio of polynomials whose coefficients all have one sign, with a nonzero constant term. So every pivot is
   positive on the whole family, and the elimination commutes with evaluation.
4. **Finite cases.** Families without a symbolic parameter have k ≤ 6. These, and every k ≤ 40, are also eliminated
   exactly (`code/num_all.py`).
5. **Kernel.** It follows from the nullities of the same eliminations. ∎

**Proposition E2.** At equal weights, every block in which one copy carries the letter 2 in row and column is positive
definite, for every k and j. Its reduced form is one-dimensional per i, with a closed Clebsch–Gordan formula; the
formula is checked against the explicit forms (`code/check_caseII.py`).

**Letter-2 blocks of case I.** For the type I blocks the symmetry is S_T × S_U, with τ = |T| and u = k − τ. The T-part
carries only the irreps (τ) and (τ − 1, 1): the S_T-invariant vectors ("sym", at most 10) and the sign twist on a pair
of T ("anti", at most 3). The U-part is the letter-1 module on u copies with the twists of the reduction above. The
entries come from the same orbit calculus, with τ as a fourth parameter; it enters only through 𝒩 (each letter-2 copy
adds 1) and through orbit sizes. The reduced forms agree with the independent builder for k = 3, 4 (598 entries), with
explicit vectors for k ≤ 10 (1003 forms), and are ⪰ 0 by exact elimination for every k ≤ 14 (`code/run2.py`) and
k ≤ 16 (`code/num_letter2.py`, 4998 forms).

**Theorem E3.** At equal weights, for every k, every j and every T ≠ ∅ with no copy carrying the letter 2 twice, every
reduced form of the case I block is ⪰ 0. Anti forms are positive definite; sym forms have kernel exactly the part of
p ⊗ H_j in the block. (Here q0 = j − τ − i and t = k − j − i.)

*Proof (computer-assisted).*
1. **Families** (`code/prove_letter2.py`). The parameters split into families q0 ∈ {−1, 0, 1, ≥ 2},
   t ∈ {−1, …, 3, ≥ 4}, i ∈ {0, 1, ≥ 2}, plus the line q0 = 2, t = −1; τ is always symbolic, and both modes are
   treated. In each family the script checks that the Newton window is exact, compares the symbolic form with the
   numeric calculus at sample points, checks the kernel vector, and eliminates with pivots that are ratios of
   one-signed polynomials after each parameter is shifted to its minimum. All anti families pass, and all sym
   families except the generic one [q0 ≥ 2, t ≥ 4, i ≥ 2].
2. **The generic sym family** (`code/prove_letter2_split2.py`). An explicit basis without poles, with one-signed
   determinant, splits the form by S_U-irreps. The (u − i − 1, i, 1) part has two positive pivots, for every k. The
   (u − i, i) part is a 7 × 7 form: its first four pivots are positive on the whole family. With one more parameter
   fixed the elimination closes: τ = 1 for all q0 ≥ 2, t ≥ 4, i ≥ 2 (`code/prove_tau_loop.py 1`), and t = 4, 5 for all
   q0 ≥ 2, i ≥ 2, τ ≥ 1 (`code/prove_fixed_loop.py t 4 5`).
3. **The (u − i, i) part of the generic family** is Proposition C below. (Before it was found, the corner
   {q0 ≥ 2, t ≥ 6, i ≥ 2, τ ≥ 2} was eliminated exactly for 14 ≤ k ≤ 38: 11 011 forms, each ⪰ 0 with nullity exactly 1,
   `code/num_generic.py 14 38 6 2`; this is now a consistency check.) ∎

**Proposition C.** On the whole generic family F = {q0 ≥ 2, t ≥ 4, i ≥ 2, τ ≥ 1}, the form Q8 on
span{X1, X2, X3, X4, G0, C1, D0, D1} (the (u − i, i) part) is ⪰ 0 with kernel exactly span(κ), κ = X1 + X2 + C1.

*Proof (computer-assisted, `code/corner/`, about 15 minutes, exact arithmetic throughout).*
1. **Reduction.** κ lies in the kernel identically and has coefficient 1 on C1, so in the basis
   X1, X2, X3, X4, G0, D0, D1, κ the form is Q7 ⊕ 0. It suffices to show Q7 ≻ 0 on F.
2. **Continuity.** The entries of Q7 are rational functions whose denominators are products of 17 explicit factors,
   each of which has one-signed coefficients after shifting the parameters to their minimum. So Q7 is continuous on
   the real box F_R = [2, ∞) × [4, ∞) × [2, ∞) × [1, ∞), which is connected, and equals the reduced form at every
   integer point of F (the Newton window is exact there).
3. **A structure lemma.** With ν = 2i + τ, the calculus depends on (i, τ) only through ν and through orbit sizes of
   degree at most 2. So every entry of Q7 is a polynomial of degree at most 2 in τ over Q(q0, t, ν), with a τ-free
   denominator (checked for all 49 entries), and A·det Q7 has degree at most 14 in τ for the common denominator A.
4. **The determinant.** Exact determinants over Q(q0, t, ν) at τ = 2, …, 16 and Lagrange interpolation give
   det Q7 = 8 τ³ i² (i + q0 + t + 1)² q0 t³ · P/A, with P an integer polynomial of 80 573 terms (total degree 73) and A a
   product of the 17 factors. This is cross-checked against direct determinants at 69 random points and against the
   independent numeric calculus at 6 more points (k up to 44).
5. **No zero on the box.** After q0 = 2 + x, t = 4 + y, i = 2 + z, τ = 1 + w, all coefficients of P are positive and
   its constant term is nonzero. So det Q7 > 0 on F_R.
6. **Conclusion.** Q7 ≻ 0 at (q0, t, i, τ) = (2, 4, 2, 1) by exact elimination. Since det Q7 never vanishes on the
   connected set F_R, the inertia of Q7 is constant there, so Q7 ≻ 0 on all of F_R. ∎

Together with Theorem E1, Proposition E2, Theorem E3, Proposition M and the remark that (R2) always holds, all group
conditions hold at equal weights, for every k. This proves the Corollary stated at the top.

## Code and checks

Python 3 with NumPy, SymPy and python-flint (`code/flintfield.py` is a small exact rational-function field on top of
it). Every script prints PASS; the logs are in `logs/`.

| Script | What it checks |
|---|---|
| `ylib.py`, `check_builder.py` | Exact builder of the (R∞)_j blocks, checked against a dense floating-point construction of M, V_p, D, R and Γ |
| `check_boson_Y.py` | Lemma B, exactly, k ≤ 4, every j, all letters |
| `check_normal_order.py` | Lemma N and Proposition K2, exactly, k ≤ 6 |
| `check_g0_inverse.py` | Lemma G, exactly, k ≤ 4 |
| `check_top_block.py` | Theorem T: closed form, both radicals, the reduction and the bound, exactly, k = 2, …, 6 |
| `check_typeII.py`, `check_typeI_u0.py` | Propositions T2 and U0 (closed forms, exact pivots, coefficient signs) |
| `eqlib.py`, `eqlib2.py`, `orbit.py`, `orbit3.py`, `symf.py` | Equal-weight reduced forms and the orbit calculus |
| `prove_families.py` | Theorem E1 for every k (42 symbolic families) |
| `num_all.py` | Exact elimination of every reduced letter-1 form for k ≤ 40 |
| `check_vs_builder.py`, `check_orbit_vs_eqlib.py` | Reduced letter-1 forms against the builder (k ≤ 5) and against explicit vectors (k ≤ 9) |
| `check_caseII.py` | Proposition E2 |
| `check2_vs_builder.py`, `run2.py` | Reduced letter-2 forms against the builder (k = 3, 4); exact elimination for every k ≤ 14 |
| `orbit2.py`, `check_orbit2.py` | Orbit calculus for the case I letter-2 forms (numeric and over Q(q0, t, i, τ)); agreement with explicit vectors for every case I form with k ≤ 10 |
| `prove_letter2.py`, `prove_letter2_split2.py`, `prove_letter2_split_tau.py`, `prove_tau_loop.py`, `prove_fixed_loop.py`, `sanity_generic.py` | Theorem E3: the symbolic families, the irrep split of the generic family, the fixed-parameter families τ = 1 and t = 4, 5, and the window and sample-point checks of the generic family |
| `num_generic.py`, `num_letter2.py` | Exact elimination of the corner forms for 14 ≤ k ≤ 38, and of every case I form for k ≤ 16 |
| `typeI_u1/` (`run_all.sh`) | Proposition U1: the closed forms against the reduced form (exactly at rational points, and as identities of rational functions for m ≤ 4), the polynomials of (A) and (B), their Bernstein certificates, end-to-end checks, negative controls |
| `typeI_u2/` (`run_all.sh`, notes in `typeI_u2/NOTES.md`) | Theorem C and the reduction before it: the closed form of every type I block (exactly against the builder, k = 3, 4, 5), the structure facts, the splitting formulas, the hard-part lemma, Φ1, every per-coordinate bound at 200 rational configurations, the chain at 32 points, exact pivots at 33 configurations, the Bernstein certificates of (I), (II), (III), negative controls |
| `typeI_G_u3/` (`run_all.sh`, notes in `typeI_G_u3/NOTES.md`) | Theorem U3₁: the exact Schur complement for u = 3, m = 1 over Q(s_0, …, s_3) and its agreement with direct elimination, the leading minors with factor stripping and expansion checks, the coefficient-sign certificate, the Z-matrix certificate, negative controls; the same minor certificate for (u, m) = (1, 1), (2, 1), (1, 2), (2, 2), (1, 3); numerical sweeps behind Conjectures P and M |
| `typeI_G_u2/` (`run_all.sh`, notes in `typeI_G_u2/NOTES.md`) | Theorem U2: the exact symbolic Schur complement for u = 2, its agreement with direct elimination (m = 1, …, 5), the polynomial reductions, the blow-up charts and Bernstein certificates of (A_2) and (B_2), negative controls |
| `middle_k4/` (`run_all.sh`, notes in `middle_k4/NOTES.md`) | Section 3f: the torus blocks of (R∞)_1 and (R∞)_2 at k = 4 over Q(s_1, s_2, s_3), the known kernels and the gauge, the Bareiss leading-minor certificates, minors against direct determinants at rational points, symbolic blocks against the point builder and the dense construction, negative controls |
| `middle_k2/` (`run_all.sh`, notes in `middle_k2/NOTES.md`) | Section 3e: the closed form of the letter-1 block of (R∞)_{k−2} against a direct construction (k = 4, 5), Lemmas R1, J, W, BJ and Proposition S exactly (k = 4, 5, 6), and a numerical sweep (k = 4, …, 7) |
| `corner/` (`run_all.sh`) | Proposition C: the node determinants and the interpolation (`compute_nodes.py`, `interpolate.py`), the coefficient signs (`shift_check.py`), the structure lemma, window exactness and agreement with the numeric calculus (`check_structure.py`), the determinant at random points and the base point (`check_det_points.py`), and the irrep split (`check_split.py`); logs in `corner/logs/` |
| `check_bottom_condition.py`, `c0_interval.py` | Theorem B0: closed form against a direct construction (k ≤ 5), the letter-2 blocks, the reduction to (C0), the bound chains at random points, and the two reduced inequalities in interval arithmetic |
| `check_bottom_algebraic.py` | Theorem B0, algebraic proof: the identity E − D = G, the bound T ≤ G and the sum-of-squares identity exactly at rational points; the letter-2 estimate in SymPy; the bound L ≤ (2 + 2γ)S at random complex points |
| `check_typeII_allj.py` | Proposition T2′: the one-copy facts exactly, and the lower bound for the row u = 2@t against a direct construction (k = 2, …, 5, random and extreme weights) |
| `check_top_condition.py` | Theorem K: closed form against a direct construction (k ≤ 5, every T, all 2k·2^k states), the identities C = PL and the entries of F exactly at rational weights, the dominance margin 11/16, and the gap at random and extreme weights |
| `check_typeI_reduction.py` | The Remark after Proposition U0 (the block K, its inverse and lower bound, block sizes; k = 3, 4, 5) and the constant 5/24 in Theorem T (the scalar lemma in SymPy, the smallest eigenvalue of the form) |
