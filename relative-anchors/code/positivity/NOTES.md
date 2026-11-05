# Why are the minors of the rescaled (R∞)_j blocks subtraction-free? Exact tests at k = 2, 3, 4, an obstruction, and the row-group structure

Ansh Mishra and Aryan Senthilkumar, 5 November 2026. Nothing here is refereed.

Notation as in `relative-anchors/README.md` (Sections 1–3e). Blocks are the torus blocks of (R∞)_j in the rescaled
basis of `relative-anchors/code/ylib.py` (x_X = Ψ_X/μ_X, μ_X = ∏_t w_t^{N_t(X)}, entries in Q(s), s_t = w_t²), with
the known kernel removed by deleting coordinates ("deletion gauge", as in `../middle_k4/k4lib.py`; Section 1). A rational
function is called *subtraction-free* (SF) if it is a ratio of two nonzero polynomials with positive coefficients. A
symmetric matrix is *SF-P* if every principal minor is SF. Nothing below says anything about Conjecture 12 in general
or about the NPT (Werner distillability) problem.

## Status

| # | Statement | Status |
|---|---|---|
| 1 | k = 2 and k = 3, every torus block of every (R∞)_j, j = −1, …, k (one per orbit of the copy permutations: 12 blocks at k = 2, 20 at k = 3, type II included): **every** principal minor of the gauge-fixed block is SF (46 and 916 minors) | **exact check** (`allminors.py`) |
| 2 | k = 4: every principal minor of the blocks j = 1, T = {0, 1} (7), j = 2, T = {0, 1, 2} (15) and j = 1, T = {0} (all 1023) is SF | **exact check** (`allminors.py`; 7 min for the 1023) |
| 3 | k = 4: random principal minors of all seven blocks, 12 random 7-subsets per block, each with its 7 nested minors (36 for the 3 × 3 block, 48 for the 4 × 4 block; 504 minors, small subsets can repeat): all SF. This includes the block j = 2, T = {0}, whose leading minors are certified in `../middle_k4/` | **exact check** of a random sample (`allminors.py … sub:12:7:2026`) |
| 4 | k = 4: every principal minor of every *row group* (states with one fixed row string) of every block: 447 minors, all SF | **exact check** (`rowgroup_minors.py`) |
| 5 | Lemma 1: SF-P is preserved by principal submatrices, Schur complements, inversion and diagonal congruence by SF functions; a positive definite matrix Σ_e c_e x_e x_eᵀ with fixed real vectors x_e and SF weights c_e is SF-P (Cauchy–Binet) | **proved** (elementary) |
| 6 | Corollary 1: in a fixed deletion gauge, SF-P of the gauge-fixed type I block of (R∞)_{k−1} implies Conjecture P (`typeI_G_u3/NOTES.md`) for the Schur complement taken in that gauge | **proved** (from Lemma 1) |
| 7 | **Proposition R** (every k, every j, every block): the row group of the vacuum row and of every letter-1 row is, up to a diagonal monomial congruence and a power of 2, the Schur complement of an explicit matrix with a fixed-vector representation (0/1 vectors, weights Laurent monomials or nonnegative linear forms in s). Hence every principal minor inside one such row group is SF | **proved** (uses Lemmas B and G of the README); the explicit lifts are also checked exactly at k = 2, 3, 4 (`rowgroup_lift.py`) |
| 8 | Letter-2 rows: the same Schur-complement identity holds with an explicit C_a (checked exactly at k = 3, 4), but in the coordinates forced by the lift C_a has the obstruction of Lemma 2, so the proof of Proposition R does not extend | identity **exact check**; structure **open** |
| 8b | Column groups (states with one fixed column string): the partial transpose does not act, and the group is the Schur complement of a Gram-type lift; the lift has a fixed-vector representation for 12 of 14 column groups at k = 3 and 13 of 15 at k = 4 (exceptions: columns with a letter 1 at the gauge copy) | **exact check** (`col_lift.py`); every k **open** |
| 8c | No obstructed pair (Lemma 2) shares its column: at k = 3 and 4 they lie inside one row group (26, 124) or differ in both row and column (27, 200) | **exact check** (`obstr_cols.py`) |
| 9 | Lemma 2 (ratio obstruction): if M = Σ_e c_e x_e x_eᵀ with finitely many fixed real vectors and c_e ≥ 0, then M_ab/M_aa is bounded on the orthant | **proved** |
| 10 | By Lemma 2, **no** representation with fixed vectors and nonnegative weights exists for: k = 2: the letter-1 block of (R∞)₁; k = 3: 6 of 20 blocks (the letter-1 blocks of (R∞)₀, (R∞)₁, (R∞)₂, the type I blocks (R∞)₁ with label (1,1,0,0) and (R∞)₂ with label (2,1,0,0), the type II block (R∞)₂ with label (1,2,0,0)); k = 4: 5 of 7 blocks, namely all except j = 1, T = {0, 1} and j = 2, T = {0, 1, 2} (so every large block at k = 4) | **proved** (exact certificate per block: a pair (a, b) and a monomial curve along which M_ab/M_aa → ∞; `newton_obstr.py`) |
| 11 | Lemma 3: the obstruction survives **every** rescaling of the coordinates by monomials (any real exponents) for all these blocks except the type II block at k = 3 | **proved** (exact width or Farkas certificates; `newton_obstr.py`, `rescale_farkas.py`) |
| 12 | A fixed-vector representation with weights s^m/den (every monomial coefficient matrix of den·M_g positive semidefinite) exists for: k = 2: 10 of 12 blocks; k = 3: 13 of 20; k = 4: the two blocks j = 1, T = {0, 1} and j = 2, T = {0, 1, 2}; and for the type II block (1,2,0,0) at k = 3 after the rescaling diag(s₁, s₂) | **exact check** (`mono_psd.py`, `mono_rescale.py`) |
| 13 | Pattern for j < k: at k = 3 and 4 a representation as in 12 exists exactly for the blocks without a letter 1 (n₁ = 0; at k = 2 also the letter-1 block of (R∞)₀). All blocks of (R∞)_k have one at k ≤ 3 | **exact check** at k ≤ 4; for every k a **conjecture** |
| 14 | Two blocks have no obstruction of Lemma 2 but no representation was found: k = 2, (R∞)₁, label (1,1,0); k = 3, (R∞)₂, label (1,1,1,0). Multipliers \|s\|^N (N ≤ 12 resp. 6), products of den factors and (s₀ + c s₁)^N (k = 2) fail; no first- or second-order jet obstruction exists | **open** (`polya_mult.py`, `linear_mult.py`, `jet_obstr.py`) |
| 15 | k = 2, letter-1 block of (R∞)₁: with D = diag(\|s\|/s₀, 1), D M D = s₀s₁[½(1,1)(1,1)ᵀ + diag((s₀² + s₁²)/(2Δ), ½)], Δ = s₀² + 4s₀s₁ + s₁² | **proved** (exact identity, `k2_rescaled.py`) |
| 16 | The natural whole-block lift (one copy of D′ per pair of rows, Woodbury for the antisymmetric part; Section 6) is exact and PSD, but at k = 3 it is obstructed against every monomial rescaling for the letter-1 blocks of (R∞)₀, (R∞)₁ and for (R∞)₁ label (1,1,0,0) | lift identity **exact check**; obstruction **proved** (Farkas certificates, `lift_test.py`) |
| 17 | A representation of the whole block (fixed vectors after a Schur complement, Conjecture S below) uniform in k; a proof of the group conditions for every k | **open**; not obtained |

**In short.** (1) The literal goal (a positive combination of squares of fixed integer linear forms
with SF weights) is **impossible** for every block that matters at k = 3 and 4, by an exact obstruction, and it stays
impossible after any monomial rescaling of the basis. (2) What does hold, for every k, is that each vacuum or letter-1
row group is a *Schur complement* of a fixed-vector matrix (Proposition R). Schur complements produce exactly this
kind of obstruction (unbounded entry ratios) and keep subtraction-free principal minors (Lemma 1). The obstructed
pairs inside one row group (about half of them) are explained this way. All the others differ in both row and column
(no obstructed pair shares its column), so they sit exactly on the entries created by the partial transpose. Column
groups, where the partial transpose does not act, have a Gram-type lift that is a fixed-vector matrix in 25 of 29
cases at k = 3, 4 (exact checks). So the positivity pattern is most likely of "Kron reduction" type, not of Gram type
(Conjecture S, Section 6). The natural whole-block lift (one copy of D′ per pair of rows) is exact but is itself
obstructed (row 16), so a whole-block explanation needs a different lift. (3) Every principal minor is SF on all
blocks at k = 2, 3, on the small blocks and the block j = 1, T = {0} at k = 4, and on a random sample of the large
k = 4 blocks.

## 1. Blocks and gauge

`build_blocks.py k` (k = 2, 3) builds one block per orbit of the copy permutations for every j = −1, …, k, exactly over
Q(s₀, …, s_{k−1}) with `ylib.Builder`. It checks symmetry and agreement with the point builder at 3 rational points.
The kernel is computed exactly at a rational point and the gauge is chosen greedily. Coordinates (1@1, y) go first,
then the other rows with a letter 1, then vacuum states. The kernel restricted to the deleted set must be invertible.
For the known kernels this is the gauge of `../middle_k4/k4lib.py`: (1@1, y) for every y, plus (vac, hole 0) for the
letter-1 block of (R∞)_{k−2}. For the letter-1 block of (R∞)_{k−1} (Theorem T) it also deletes (1@r, hole 0), r ≠ 1.
The gauge-fixed block is positive definite at 3 rational points for every block (exact LDL). At k = 4 the seven blocks
of `blocks_sym.pkl` from `../middle_k4/` (s₀ = 1) are used as they are.

The tests below are about the gauge-fixed matrix M_g. For blocks whose kernel vectors are constant in the rescaled
basis (all except the letter-1 blocks of (R∞)_{k−2} and (R∞)_{k−1}), a fixed-vector representation of the full block
restricts to M_g, and M_g is a principal submatrix of the full block. So every negative statement for M_g holds for
the full block too.

## 2. Lemma 1 (closure) and Corollary 1

**Lemma 1.** Let M be symmetric over Q(s) and positive definite on the open orthant.
- (a) If M is SF-P, so is every principal submatrix.
- (b) If M is SF-P, so is every Schur complement M/M[α]. This uses det (M/M[α])[β] = det M[α ∪ β] / det M[α].
- (c) If M is SF-P, so is M⁻¹, by Jacobi: det M⁻¹[β] = det M[β^c] / det M.
- (d) If M is SF-P, so is DMD for D diagonal with SF entries.
- (e) If M = Σ_e c_e x_e x_eᵀ with fixed real vectors x_e and SF weights c_e, then M is SF-P. By Cauchy–Binet,
  det M[β] = Σ_{|E| = |β|} ∏_{e∈E} c_e · det X[β, E]². This is a nonnegative constant combination of SF functions, and
  it is not identically zero because M is positive definite.

**Corollary 1.** Conjecture P of `typeI_G_u3/NOTES.md` concerns the leading minors of G_z. G_z is the Schur complement
of the gauge-fixed type I block of (R∞)_{k−1} onto its letter-2 states, after a diagonal monomial rescaling. Its leading
minors are therefore ratios of principal minors of the gauge-fixed block. So SF-P of that block (in the same deletion
gauge) implies Conjecture P, in every order. At k = 3 the type I blocks of (R∞)₂ are SF-P in our gauge (row 1 of the
status table). For the label (2,1,0,0) this is the gauge convention of Theorem C (deleted row in U). For (1,1,1,0) our
deleted row lies in T, so that case is a different gauge.

## 3. The obstruction to fixed-vector representations

**Lemma 2.** Let M(s) = Σ_e c_e(s) x_e x_eᵀ on the open orthant, with finitely many fixed real vectors x_e and
c_e(s) ≥ 0. Then |M_ab| ≤ K_ab M_aa with K_ab = max{|x_e(b)/x_e(a)| : x_e(a) ≠ 0}.

*Proof.* |M_ab| ≤ Σ_{e: x_e(a)x_e(b) ≠ 0} c_e |x_e(a)| |x_e(b)| ≤ K_ab Σ_e c_e x_e(a)² = K_ab M_aa. ∎

Write M_ab/M_aa = P/Q with P = n_ab d_aa and Q = d_ab n_aa. Here Q has positive coefficients: all denominators and all
diagonal numerators of every block have positive coefficients (checked). Then the ratio is bounded iff
conv(supp P) ⊆ conv(supp Q). `newton_obstr.py` certifies unboundedness exactly. It gives an integer direction l for
which the maximum of l·m over supp P is attained at a single monomial and exceeds the maximum over supp Q, so the ratio
grows like t^gap along s = t^l. Results (all certificates found among l ∈ {−3, …, 3}^r):

| k | blocks with an obstruction | examples |
|---|---|---|
| 2 | 1 of 12: letter-1 block of (R∞)₁ | M[(10,10),(00,11)]/M[(10,10),(10,10)] = Δ/(2s₀(s₀+s₁)) |
| 3 | 6 of 20: letter-1 blocks of (R∞)₀, (R∞)₁, (R∞)₂; (R∞)₁ label (1,1,0,0); (R∞)₂ labels (2,1,0,0) and (1,2,0,0) | M[(001,001),(001,010)]/M[(001,001),(001,001)] ~ s₁ as s₁ → ∞ (letter-1 block of (R∞)₁) |
| 4 | 5 of 7: j = 1, T = ∅ and {0}; j = 2, T = ∅, {0} and {0, 1} | M[(0001,2001),(0001,2010)]/M[(0001,2001),(0001,2001)] ~ s₁ (j = 2, T = {0}) |

Obstructed pairs occur both inside one row group and between row groups (`obstr_rows.py`, `obstr_cols.py`, ordered
pairs with an exact certificate):

| k | same row string | same column string | neither ("crossed" pairs) |
|---|---|---|---|
| 2 | 0 | 0 | 1 |
| 3 | 26 | 0 | 27 |
| 4 | 124 | 0 | 200 |

**No obstructed pair shares its column.** Every obstructed pair either lies in one row group or differs in both row and
column. On the latter pairs the B-entry ⟨MV′(b, y′), D′⁻¹MV′(a, y)⟩ pairs each state with the column of the other:
this is exactly where the partial transpose acts. Inside a row group, the typical pair consists of two states that
differ by moving one letter 1 between two copies. This is the hopping term ½J⁺J⁻ seen through E_j⁻¹.

**Lemma 3 (monomial rescalings).** Rescaling the coordinates by monomials, x_a → s^{δ_a} x_a with any real δ,
multiplies M_ab/M_aa by s^{δ_b − δ_a}. This translates conv(supp P). Hence:
- (width) if conv(supp P) is wider than conv(supp Q) in some direction, no rescaling makes this ratio bounded;
- (Farkas) the conditions (δ_b − δ_a)·l ≤ g_ab(l) = max_{supp Q} l·m − max_{supp P} l·m are necessary for every
  pair and every l. A nonnegative rational combination of them that sums to 0 ≤ (negative number) shows that no δ
  exists.

The inequality holds for every l, even when several monomials tie for the maximum. The tied part of P is a nonzero
face polynomial, which does not vanish at a generic base point.

`rescale_farkas.py` finds the support by linear programming. It then recomputes the multipliers exactly over Q and
checks every condition exactly. Exact certificates exist for 1/1 obstructed block at k = 2, 5/6 at k = 3 and 5/5 at
k = 4 (objective values −1, −1/2, …, −4/5). The exception is the k = 3 type II block (1,2,0,0). For it the LP finds the
rescaling diag(s₁, s₂), and after that every monomial coefficient matrix is PSD (`mono_rescale.py`), which is
row 12 of the table.

So for all blocks that matter, no basis obtained from the rescaled one by monomials carries a fixed-vector
representation. Non-monomial rescalings can help. At k = 2 the letter-1 block of (R∞)₁ becomes a fixed-vector matrix
after D = diag(|s|/s₀, 1) (row 15, `k2_rescaled.py`).

## 4. Where representations do exist

`mono_psd.py`: with den the lcm of the entry denominators, a representation M_g = Σ_e (p_e/den) x_e x_eᵀ with fixed
vectors and polynomials p_e with nonnegative coefficients exists iff every monomial coefficient matrix C_m of den·M_g
is PSD. For ⇐, an LDLᵀ of each C_m gives the vectors. For ⇒, collect the monomials. The test is exact: elimination with diagonal pivoting, and a
rational witness x with xᵀC_m x < 0 when the test fails. The multipliers |s| and |s|² change no verdict (a
multiplier with positive coefficients keeps PSD coefficient matrices PSD; on the failing blocks some C_m stays
non-PSD).

| k | blocks with every C_m PSD (q = 1) |
|---|---|
| 2 | 10 of 12: all except (R∞)₁ labels (1,1,0) and (2,0,0) |
| 3 | 13 of 20: (R∞)_{−1}; (R∞)₀ label (0,1,0,0); (R∞)₁ labels (0,2,0,0), (0,1,1,0); (R∞)₂ labels (0,2,1,0), (0,1,1,1); all 7 blocks of (R∞)₃ |
| 4 | the two small blocks j = 1, T = {0, 1} (455 monomials) and j = 2, T = {0, 1, 2} (969 monomials) |

So for j < k these are the blocks without a letter 1 (n₁ = 0), with one exception at k = 2 (the letter-1 block of
(R∞)₀; Theorem B0 needs k ≥ 3). These blocks are the easy ones (Propositions U0, T2-type). The blocks with letters 1
are obstructed at k = 3, 4, except (R∞)₂ label (1,1,1,0) at k = 3 (row 14 of the table).

For the two unobstructed blocks without a representation (row 14): no non-PSD C_m sits at a vertex of the Newton
polytope. A first-order (vertex) or second-order (C_{m₁} on ker C_v) violation would exclude every SF-weighted
representation (`jet_obstr.py`; the argument is in its docstring), but none was found. Pólya-type multipliers did
not help either. These two blocks are undecided.

## 5. Proposition R: every row group is a Schur complement of a fixed-vector matrix

Fix k, j, a torus block and a row string a. The *row group* is the set of states (a, y) of the block. By Lemma B of the
README, on a row group with a = (letter ℓ at r) and y, y′ ∈ H_j,

    2^Z Y[(a, y′), (a, y)] = ⟨y′| Δ_a − 𝒯_aᵀ E_j⁻¹ 𝒯_a |y⟩,   Δ_a = (2 − n_r) − ½ σ_a⁺σ_a⁻ (diagonal: 2, ½, 1 on y_r = 0, ℓ, other),

and by Lemma G the vacuum row gives 2^{k−j−1} Y = I − ½ J⁻E_{j+2}⁻¹J⁺ on H_{j+1}. So 2^Z Y_aa = 𝔹_a/E with
𝔹_a = [[Δ_a, 𝒯_aᵀ], [𝒯_a, E]]. The other Schur complement is 𝔹_a/Δ_a = E − 𝒯_a Δ_a⁻¹ 𝒯_aᵀ =: C_a. A short computation
(the cross terms vanish because σ_a⁻ kills |0⟩ and |2⟩ at r) gives:
- vacuum row (Δ = I, 𝒯 = J⁺/√2): C = Σ_t s_t(2|0⟩⟨0|_t + |2⟩⟨2|_t), diagonal;
- letter 1 at r: C_a = C^{(r)} + ½ J⁺_{(r)} (1 − P_{r,0}) J⁻_{(r)}. Here C^{(r)} is C without copy r, J_{(r)} is J without
  copy r, and P_{r,0} projects on strings with a 0 at r;
- letter 2 at r: C_a = C + ½ J⁺_{(r)} (1 − P_{r,0}) J⁻_{(r)} + Σ_{t≠r} w_r w_t (h_{t←r} + h_{r←t}), where h moves a letter 1.

`rowgroup_identity.py` checks all three, and Lemma B on the diagonal row blocks, exactly at two rational points for
k = 3 and 4 (with w rational, so s = w²): 2 104 entries of Lemma B and 50 + 50 operator identities. It also checks
that C_a ⪰ 0 for the letter-2 rows at these points.

**Proposition R.** For every k ≥ 2, every j, every torus block, and a = vacuum or a = letter 1 at r: in the coordinates
ψ′_y = s^{N(y)} x_{(a,y)} (a monomial rescaling of the ylib coordinates; N(y) = letters 1 in y) the row group is
2^{−Z} (resp. 2^{−(k−j−1)}) times the Schur complement onto ψ′ of

    letter 1 at r:  𝔹′_a = Σ_{y∈H_j} (Δ(y) s_r / s^{N(y)}) f_y f_yᵀ + Σ_{x∈H_j} (C^{(r)}_x / s^{N(x)}) e_x e_xᵀ
                           + ½ Σ_{u∈H_{j−1}, u_r ≠ 0} s^{−N(u)} g_u g_uᵀ,
                    f_y = e_{ψ′,y} + e_{z′,y} + [y_r = 1] Σ_{t: y_t = 0} e_{z′, y−1_r+1_t},   g_u = Σ_{t≠r, u_t=0} e_{z′, u+1_t};
    vacuum:         𝔹′ = Σ_{y∈H_{j+1}} s^{−N(y)} f_y f_yᵀ + Σ_{x∈H_{j+2}} (2C_x / s^{N(x)}) e_x e_xᵀ,   f_y = e_{φ′,y} + Σ_{t: y_t=0} e_{z′, y+1_t}.

Every vector has entries 0 or 1 and every weight is a Laurent monomial or a linear form with nonnegative coefficients.
Hence (Lemma 1 (b), (d), (e)) every principal minor of the row group, gauge-fixed or not, that is nonzero is SF.

*Proof.* Write 𝔹_a = L diag(Δ_a, C_a) Lᵀ with L = [[I, 0], [𝒯_aΔ_a⁻¹, I]]. The quadratic form is then
Σ_y Δ_a(y)(ψ + Δ_a⁻¹𝒯_aᵀz)_y² + ⟨z, C_a z⟩. Now ⟨z|𝒯_a|y⟩ is w_r c(y) for z = y, with c = 2, ½, 1 on y_r = 0, 1, 2, and
½w_t for z = y − 1_r + 1_t (y_r = 1, y_t = 0). Put ψ_y = (w_r/μ_y) ψ′_y and z_x = z′_x/μ_x. Since c(y)/Δ(y) = 1 and
½/Δ(y) = 1 when y_r = 1, the linear forms become (w_r/μ_y) f_yᵀ(ψ′, z′). The C_a part becomes the stated diagonal and
g-terms, because ⟨J⁺_{(r)}u, z⟩ = μ_u⁻¹ Σ_t z′_{u+1_t}. The vacuum row is the same computation with
ψ = φ, z = √2 z̃ and E_{j+2} = C + ½J⁺J⁻. Finally ψ = Ψ restricted to the row a, and x = Ψ/μ_X with μ_X = w_r μ_y,
so ψ′ = s^{N(y)} x. ∎

`rowgroup_lift.py` builds 𝔹′ from these formulas and compares its Schur complement with ylib's row blocks, exactly, for
every row, every j (j = 0, …, k for letter rows, j = −1, …, k − 1 for the vacuum row) and every sector at k = 2, 3, 4
(two rational points each, 21 368 entries, all equal).

**Column groups** (states with one fixed column string; `col_lift.py`). On a column group the partial transpose does
nothing. The B-part is an ordinary Gram matrix ⟨MV′(b, y), D′⁻¹MV′(a, y)⟩, so the column group is the Schur complement
of 𝔹_col = [[A′_col, Wᵀ], [W, D′]] with W[z, a] = MV′(a, y)[z]. Exact checks: the identity and 𝔹_col ⪰ 0 at a
rational point, for all 14 column groups with at least two states at k = 3 and all 15 at k = 4. 𝔹_col has a
fixed-vector representation (q = 1, after an LP monomial rescaling where needed) for 12 of 14 at k = 3 and 13 of 15 at
k = 4. The four exceptions are exactly the columns with a letter 1 at the gauge copy g = 1 (010 and 210 at k = 3, 0100
and 2100 at k = 4), where the gauge has deleted the state (1@1, y). We have no proof for every k.

Taken together: inside row groups the structure is a Schur complement (Proposition R, proved). On column groups it is a
Gram-type Schur complement (exact checks). The obstruction to fixed-vector representations sits inside row groups and on
the crossed pairs, which are neither.

What this explains, and what it does not:
- It explains SF of every principal minor inside one row group, for every k (row 4 of the table at k = 4). It also
  explains the obstructed pairs of Lemma 2 that lie inside one row group (about 40–50 % of them, Section 3): a Schur
  complement of a fixed-vector matrix (a "Kron reduction") has unbounded entry ratios in general. The pairs between
  row groups are not explained.
- Minors that meet two row groups are not covered. Rows a ≠ b are coupled by the *crossed* term
  −⟨𝒯_bψ_a, E⁻¹𝒯_aψ_b⟩ (the partial transpose), plus the crossed swap term. Writing
  2⟨u, E⁻¹v⟩ = ⟨u + v, E⁻¹(u + v)⟩ − ⟨u, E⁻¹u⟩ − ⟨v, E⁻¹v⟩ turns the crossed sum into Schur-complement terms plus a
  block-diagonal term +Σ_{a≠b}⟨𝒯_bψ_a, E⁻¹𝒯_bψ_a⟩ with the wrong sign for a Schur complement. Polarising the other
  way moves the wrong sign to the antisymmetric part. So the lift of Proposition R does not extend directly to a whole
  block.
- Letter-2 rows: in the coordinates forced by the forms f_y, the hopping h_{t←r} has weight 1/μ_c² against a diagonal
  C_z/(μ_c² s_r) with no s_r in C_z (z_r = 1). That is the obstruction of Lemma 2 for C_a. So the proof does not carry
  over, although the row groups are SF-P at k ≤ 4 (rows 1, 4).

Related exact checks (`inv_test.py`, k = 2, 3):
- The inverse of every row group, in the coordinates η = μΨ, has every C_m PSD (q = 1). When C_a is diagonal on the
  sector this is what Woodbury gives, and C_a is always diagonal for the vacuum row (Lemma G). At k = 3 it also holds
  in the cases where C_a has hopping terms (letter-1 rows of (R∞)₂).
- The inverse of the whole block, and the Schur complement onto the rows with a letter, are obstructed for the main
  blocks. The whole-block inverse has width certificates for the letter-1 blocks of (R∞)₀ and (R∞)₂ and for (R∞)₁
  label (1,1,0,0).
- Exception: for the letter-1 block of (R∞)₀ at k = 3, the Schur complement onto the letter-1 rows has a fixed-vector
  representation.

## 6. A structural conjecture, and suggested routes

**Conjecture S.** For every k, j and every letter-1 or type I block, the gauge-fixed block is, up to a diagonal monomial
congruence, a Schur complement 𝔹/𝔹[Z] of a symmetric matrix 𝔹 with a fixed-vector representation with SF weights.

By Lemma 1 this gives SF-P of the block, hence every leading-minor certificate in every order and Conjecture P. With 𝔹
positive definite on the orthant it also gives the group conditions (Conjecture 12) for every k. Status: true for
every row group of a vacuum or letter-1 row (Proposition R); open for whole blocks. Consistent with all exact checks
here (rows 1–4). It is not implied by them.

Suggested routes:
- (i) Look for a lift of the crossed coupling. Both Lemma B terms between rows a ≠ b are of the form
  Σ⟨K_bψ_a, L K_aψ_b⟩ = ⟨V_sym, L V_sym⟩ − ⟨V_anti, L V_anti⟩. For the B-term (L = E_j⁻¹) the antisymmetric part
  enters Y with the wrong sign for a Schur complement. Woodbury with E = C + ½J⁺J⁻ gives, whenever C is invertible
  on the sector (always for j < k):

      ⟨u, E⁻¹u⟩ = ⟨u, C⁻¹u⟩ − ½⟨J⁻C⁻¹u, F⁻¹ J⁻C⁻¹u⟩,   F = I + ½J⁻C⁻¹J⁺.

  Here C⁻¹ is diagonal with SF entries 1/C_z, and F (on H_{j−1}) has the fixed-vector form of Lemma G. So the wrong-sign term
  splits into an SF-weighted diagonal term and a term of the right sign (a Schur complement over a new variable with
  matrix F). Write 𝕋ψ = (𝒯_bψ_a)_{(a,b)} for the vector indexed by ordered pairs of rows, and P_a for its
  antisymmetrisation in (a, b), so the wrong-sign term is ⟨P_a𝕋ψ, E⁻¹ P_a𝕋ψ⟩. The resulting candidate lift of the
  whole block has the top-left block A + Σ_z (P_a𝕋ψ)_z²/C_z. Here A = (M^Γ)^{⊗k} is a fixed-vector matrix with
  monomial weights in the rescaled basis: one copy of M^Γ is ½ on symmetric and 3/2 on antisymmetric letter pairs, and
  μ is the same on both letters of a pair.

  **We built this lift and it fails** (`lift_test.py`). It is built from ylib's own quantities. The vacuum row is
  treated as one more row, so the E_{j+1} and E_{j+2} terms are included; D′ = μDμ is block diagonal by charge. The
  pieces are:
  - x^T B′ x = Σ_{a,b}⟨V_ab, D′⁻¹V_ba⟩ with V_ab(x) = Σ_{y′} x_(a,y′) MV′(b, y′);
  - variables z_a (matrix D′) for the terms a = b;
  - z⁺_ab (matrix D′) for the symmetric parts;
  - for the antisymmetric parts, the Woodbury variables ζ_ab with
    ⟨Δ, D′⁻¹Δ⟩ = min_ζ [⟨Δ + Hζ, C′⁻¹(Δ + Hζ)⟩ + ⟨ζ, K⁻¹ζ⟩], where D′ = C′ + HKHᵀ, C′_z = 2^{−(k−c)} s^{N(z)} C_z,
    H[u + 1_t, u] = s_t and K = diag(2^{−(k−c)−1} s^{N(u)}).

  Checks, exact at a rational point: D′ = C′ + HKHᵀ, the Schur complement of the lift onto x equals the gauge-fixed
  block, and the lift is PSD. This holds for the blocks with j ≤ k − 2 (for j = k − 1 the string 1…1 has C = 0 and
  Woodbury does not apply). Results at k = 3:
  - letter-1 blocks of (R∞)₀ and (R∞)₁, and (R∞)₁ label (1,1,0,0): the symbolic lift (sizes 19, 31, 31) is obstructed
    against every monomial rescaling (exact Farkas certificates, objective −8/13, −40/69, −2/3);
  - (R∞)₁ label (0,1,1,0), which itself has a fixed-vector representation: after the LP rescaling the lift has no
    obstruction, but 12 of its 44 coefficient matrices are not PSD;
  - at k = 2, the letter-1 block of (R∞)₀: after the LP rescaling, 1 of 7 coefficient matrices is not PSD.

  So this natural lift is not the structure. Conjecture S, if true, needs a different lift. In Proposition R the
  coordinates are tuned so that c(y)/Δ(y) = 1; nothing like that happens for the pair variables.
- (ii) Test Conjecture S on the k = 3 blocks with an LP over candidate lifts: one copy of E_j per unordered pair of
  rows, unknown fixed vectors from a small dictionary, and coefficient matching as in `mono_psd.py`.
- (iii) Decide the two unobstructed blocks of row 14 (an SDP over multipliers q with den·q of positive coefficients).
- (iv) Run the all-subsets test on k = 4, j = 2, T = {0, 1} (2047 minors, about 2 h estimated from the 1023-minor run
  and the chain times).

## 7. Runs that failed or were stopped

- The first run of `inv_test.py 3` (01:57) grew to about 6 GB and was killed (by us and by the machine's watchdog).
  The cause was the floating-point LP of the monomial-rescaling search. For the inverse of the 8 × 8 letter-1 block
  of (R∞)₁ it built a constraint matrix with millions of columns. Since then `analyze.rescale_lp` refuses LPs with more
  than 60 000 variables ("LP skipped"). It is also not run when an exact width certificate already settles the
  question. Every run goes through `guard.sh`, which kills the process above 2.8 GB. The rerun stayed below 0.2 GB.
- A first foreground run of `inv_test.py 3` was cut by the 2-minute limit of the shell; it was rerun in the background.
- Nothing else was stopped. The longest exact run is the 1023-minor test (7 min, below 0.2 GB).

## 8. Files

| File | Content |
|---|---|
| `run_all.sh` | reruns every step, one Python process at a time under `guard.sh`; PASS/FAIL per step and ALL PASS; `LONG=1` adds the 1023-minor run and the random k = 4 sample |
| `guard.sh` | runs one Python process and kills it by its own PID above 2.8 GB |
| `ylib.py`, `flintfield.py` | copied from `relative-anchors/code/` |
| `blocks_sym_k4.pkl` | `blocks_sym.pkl` of `../middle_k4/` (the seven k = 4 blocks, s₀ = 1); `run_all.sh` rebuilds it with `../middle_k4/build_sym.py` if it is missing |
| `poslib.py`, `build_blocks.py` → `blocks_k2.pkl`, `blocks_k3.pkl` | blocks at k = 2, 3 for every j, gauge, checks (Section 1) |
| `psdtools.py`, `analyze.py` | loading, polynomial matrices, monomial coefficient matrices, exact PSD test with witness, obstruction tests |
| `allminors.py` | every principal minor (depth-first Bareiss over all subsets, factor stripping, point checks); random chains and random subsets |
| `rowgroup_minors.py` | every principal minor of every row group |
| `mono_psd.py` | monomial-PSD test (fixed-vector representations with weights s^m/den), vertex test, sign profiles |
| `newton_obstr.py` | Lemma 2 certificates and width certificates |
| `obstr_rows.py`, `obstr_cols.py` | split the obstructed pairs by row group, and by row group / column group / neither |
| `col_lift.py`, `lift_test_lib.py` | the Gram-type lift of the column groups (identity, PSD, obstruction, rescaling, monomial-PSD) |
| `rescale_farkas.py` | exact Farkas certificates against monomial rescalings |
| `mono_rescale.py` | LP search for monomial rescalings (floating point) followed by the exact monomial-PSD test |
| `inv_test.py` | inverse, vacuum-block inverse, Schur complement onto the rows, row-group inverses |
| `polya_mult.py`, `linear_mult.py`, `jet_obstr.py` | multipliers and jet obstructions for the two undecided blocks |
| `k2_rescaled.py` | the explicit non-monomial rescaling at k = 2 |
| `rowgroup_identity.py`, `rowgroup_lift.py` | exact checks of the algebra of Proposition R and of the explicit lifts |
| `lift_test.py` | the natural whole-block lift of Section 6: exact identity, PSD, obstruction, Farkas, rescaling, monomial-PSD |
| `logs/` | all logs (`run_*.log` from `run_all.sh`) |
