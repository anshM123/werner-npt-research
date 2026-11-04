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
  - (R∞)_{k−1} holds on its letter-1 block (Theorem T), on every block in which one copy carries the letter 2 in row
    and column (Proposition T2), and on the block in which every copy carries one letter 2 (Proposition U0).
  - (R∞)_{k−2} has the explicit kernel vector vac ⊗ J⁻|1…1⟩ (Proposition K2).
- **Equal weights** (w_t = 1).
  - The letter-1 block of every (R∞)_j is ⪰ 0 for every k, with the kernel described exactly (Theorem E1). The proof is
    computer-assisted, symbolic in k.
  - The blocks with one copy carrying the letter 2 in row and column are positive definite for every k
    (Proposition E2).
  - All of (R∞) holds, by exact computation, for every k ≤ 14.

**Corollary.** Conjecture 12 of `single-excitation-range/` holds for every single excitation Σ_{t∈G} x_t (each x_t a
single-copy excitation, all ‖x_t‖ equal) with |G| ≤ 14. Hence, for every n:
- every sector of every single-excitation range plane that contains such a vector is nonnegative;
- q_n ≥ 0 on every tangent plane that contains such a vector, with any co-range.

(Local unitaries fixing |0⟩ bring such a vector to Σ_{t∈G} |1_t⟩.)

**Open.**
- For general weights: the type I blocks of (R∞)_{k−1} with 1 ≤ |T| ≤ k − 1 (for k ≥ 5), and (R∞)_j for 0 ≤ j ≤ k − 2
  (for k ≥ 4; k = 4 is being certified).
- At equal weights: the "case I" letter-2 blocks for k ≥ 15.

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

   Bound the crossed pairs by AM–GM, |ρ_aρ_b X′_ab X′_ba| ≤ (4/9)·½(X′_ab² + X′_ba²), and use 2 + G > 2. This gives
   F(X′) ≥ (5/36)‖X′‖². ∎

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

**Proposition K2.** The vector vac ⊗ J⁻|1…1⟩ lies in the kernel of (R∞)_{k−2}.

*Proof.*
- **Vacuum row.** E_k = ½|s| on 1…1 and J⁺J⁻|1…1⟩ = |s||1…1⟩, so (I − ½J⁻E_k^{−1}J⁺)J⁻|1…1⟩ = 0.
- **Rows (r, 1).** By Lemma N, 𝒯_{(r,1)}J⁻|1…1⟩ = E_{k−1}σ_r⁻|1…1⟩. So the coupling σ_r⁻ − J⁻E_{k−1}^{−1}𝒯_{(r,1)}
  annihilates J⁻|1…1⟩.
- **Rows (t, 2).** These give 0, since there is no letter 2. ∎

Theorem T and Proposition K2 account for every kernel vector beyond p ⊗ H_j that our numerics found at j = k − 1 and
j = k − 2.

**What blocks the rest.** Two parts are still open:
- the type I blocks with 1 ≤ |T| ≤ k − 1;
- the conditions (R∞)_j with j ≤ k − 2.

In both, the method of Theorem T leaves soft directions: eigenvalues 0.01–0.1, and as small as 0.0015 for
(R∞)_{k−2} at spread weights. They lie in the J⁻ part of the G1 columns and tend to 0 at the boundary of the weight
simplex. So no uniform crossed-pair bound can work, and these directions have to be treated exactly.

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

**Letter-2 blocks, k ≤ 14.** For the type I blocks the symmetry is S_T × S_U. They reduce to forms of size at most
10 × 10 (S_T-symmetric) and 3 × 3 (sign twist on a pair in T). These reduced forms are checked against the independent
builder for k = 3 and 4. They are ⪰ 0 by exact elimination for every k ≤ 14, every j and every pattern
(`code/run2.py`).

Together with Theorem E1, Proposition E2 and the remark that (R2) always holds, all group conditions hold at equal
weights for every k ≤ 14. This proves the Corollary stated at the top.

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
