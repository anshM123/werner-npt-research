# Type I blocks of (R∞)_{k−1} with |U| ≥ 2: exact reduction, and the U-block for every u

Ansh Mishra and Aryan Senthilkumar, 4 November 2026. Nothing here is refereed.

Part A gives the closed form of every type I block of (R∞)_{k−1}, its exact structure, and the reduction

    block ⪰ 0 with kernel exactly span{p ⊗ f_h : h ∈ U}   ⟺   (C_u) and G ≻ 0.

Part B proves (C_u) for every u (Theorem C, computer-assisted through three exact Bernstein certificates). So these
blocks are now equivalent to G ≻ 0, a statement about the T-blocks; for |U| = 1 that is Proposition U1, and for
|U| ≥ 2 it is open.

## Part A. Closed form, exact structure and the reduction

Notation as in `relative-anchors/README.md` (Sections 1–3, Remark after Proposition U0, Proposition U1). T = copies carrying one letter 2 (m = |T| ≥ 1), U = [k] ∖ T (u = |U| ≥ 2), s_x = w_x²,
s_T = Σ_{t∈T} s_t, s_U = Σ_{h∈U} s_h, K = s_T + s_U/2. All forms are in units of 2Y (the form 2⟨Ψ, YΨ⟩ of Lemma B,
Z = 1).

### Summary of Part A

| Statement | Status |
|---|---|
| Closed form of the whole type I block for every k, T (Section A1) | **derived** from Lemmas B, N; **exact** agreement with the Lemma B builder (k = 3, 4, 5, eight (w, T)), float agreement with the direct construction (k = 4, 5, 6) |
| (S1) the X_T block is I − w_Tw_T^T ⊗ E_T^{−1}; (S2) the U-variables (X_U, q) see X_T only through ξ = X_T^T w_T | **proved** (read off the closed form), **exact** checks |
| (S3) after eliminating X_T, the U-block is the U-block of the *virtual* system with one T-copy of weight √s_T (any m) | **proved** (from S1, S2), **exact** checks (m = 2, 3; u = 2, 3) |
| (S4) the Schur complement G on the T-blocks (ρ_t, ζ_t) is pairwise: G_tt′ depends on (s_t, s_t′, s_T, w_U) only | **exact** checks (u = 2, 3, m = 4); structural reason as in u = 1 |
| Reduction: block ⪰ 0 with kernel exactly span{p ⊗ f_h} ⟺ (C_u) and G ≻ 0; G ≻ 0 ⟸ (A_u) + (B_u) with one constant H | **proved** (Section A3) |
| (C_u), (A_u), (B_u) | (C_u) is proved in Part B; (A_u), (B_u) are **open**; numerically true with the *same* extreme constants as u = 1 (sup A = 5/18, inf B = 3/4, dominance ratio ≤ 0.293), u = 2…5, m = 1…4, 9000 weight vectors; (C_u): now proved, Part B (Theorem C) |

So the u ≥ 2 case is **not** proved. What is new: an exact reduction of every type I block to (i) one positivity
statement (C_u) that depends only on (s_T, w_U), the same for every m (the T-copies act on the U-block as a single
copy of weight √s_T), and (ii) a pairwise block-dominance on the T-blocks that numerically has exactly the constants of
the proved u = 1 case.

### A1. Closed form (`cf.py`)

Variables: X_rh on P_{r,h} (r ∈ [k], h ∈ U), ρ_t on R_t, ζ_th on Z_{t,h}, q on Q. Column vectors: pattern T is R^U
(one-hole strings f_h); pattern T∖t is R^{U∪t}; pattern T∖{t,t′} is R^{U∪t∪t′}. On pattern P with free set F,
E_P = diag(s_x + s_P + s_F/2)_{x∈F} + ½ w_F w_F^T (from Lemma 1/Lemma B). Put ψ_r = X_r (row r of X, a vector on U),
y_t = ρ_t e_t + ζ_t (on U ∪ t), ξ = Σ_t w_t X_t. Lemma N gives, on one-hole strings with free set F,
- 𝒯_{(r,1)} v = ½(w·v) e_r + ½ w_r v + w_r v_r e_r if r ∈ F, and 𝒯_{(r,1)} v = w_r v if r carries a 2;
- 𝒯_{(t,2)} v = ½(w·v) e_t − ½ w_t v (the result lives on F ∪ t, pattern without t).

Then (Lemma B, quadratic-form version)

    2⟨Ψ,YΨ⟩ = ½Σ_{r≠h∈U} X_rh² + (3/2)Σ_{h∈U} X_hh² + ½(tr_U X + Σ_t ρ_t)² + (3/2)|ρ|² + Σ_{t,h}(X_th² − X_th ζ_th + ζ_th²)
             + 2q[ −tr_U X + (w_T^T X_T w_U + ½ w_U^T X_U w_U)/K − Σ_t (2Kρ_t + w_t w_U·ζ_t)/(2K − s_t) ] + 2q²
             − Σ_{r,r′∈[k]} ⟨𝒯_{(r′,1)}ψ_r, E_T^{−1} 𝒯_{(r,1)}ψ_{r′}⟩ − 2Σ_{r∈[k], t∈T} ⟨𝒯_{(t,2)}ψ_r, E_{T∖t}^{−1} 𝒯_{(r,1)} y_t⟩
             − Σ_{t≠t′} ⟨𝒯_{(t′,2)} y_t, E_{T∖tt′}^{−1} 𝒯_{(t,2)} y_{t′}⟩.

For u = 1 this is the form behind Proposition U1. Derivation: the A-part from the swap terms of Lemma B (the only
nonzero swaps are r = r′ ∈ U, the pairs (h, h′) of diagonal entries, (t, h) with X_th ζ_th, (h, t) with X_hh ρ_t, and
(t, t′) with ρ_tρ_t′); the q-row from σ_a⁻ and J⁻E_k^{−1}𝒯_a on the column 2_T 1_U (E_k diagonal there); the B-part
from Lemma N.

Check (`check_cf.py`, log `check_cf.log`, PASS): exact equality with 2 × `lemmab.LemmaB.typeI` for k = 3, 4, 5 and
eight (w, T) with u ≥ 2 and u = 1; relative error ≤ 2.2e-16 against `fastdirect.Direct` for k = 4, 5, 6.

### A2. Exact structure

**(S1), (S2).** X_T appears in the closed form only in |X_T|² − ⟨ξ, E_T^{−1}ξ⟩ (so the X_T block is
I − w_Tw_T^T ⊗ E_T^{−1}, the README block K in 2Y units), in −Σ X_th ζ_th, in the q-row through w_U·ξ, in the
U × T part of the first B-sum through Σ_{r∈U} ⟨X_r, E_T^{−1}𝒯_{(r,1)}ξ⟩, and in the cross B-sum (T-block terms).
So the U-variables (X_U, q) couple to X_T only through ξ.

**(S3) (virtual copy).** Write X_T = w_T η^T + X′_T with w_T^T X′_T = 0. Then ξ = s_T η, the X_T quadratic part is
|X′_T|² + s_T|η|² − s_T²⟨η, E_T^{−1}η⟩, and X′_T does not couple to (X_U, q). Hence the Schur complement of
{X_U, q, X_T} onto {X_U, q} is that of the system with T replaced by one copy τ of weight w_τ = √s_T (X_τ = w_τ η
gives the same ξ, the same |X_τ|² = s_T|η|², and E_T, K depend on T only through s_T). This Schur complement is the
U-block of the README Remark and of `typeIu2/explore.py`, where the smallest eigenvalue lives.

Checks (`check_struct.py`, log `check_struct.log`, PASS): (S1) and (S2) exactly for k = 4, 5 (four T); (S3) exactly
for (m, u) = (2, 2), (3, 2), (2, 3) at weights with s_T a rational square.

**(S4) (pairwise structure).** Eliminate X_T, the U-variables (X_U with P_{r0,h} deleted, and q), and call G the Schur
complement on the T-blocks (ρ_t, ζ_t). Then the blocks G_tt and G_tt′ depend only on (s_t, s_t′, s_T, w_U): every
term of the closed form is local in one t, a pair term (t, t′), or goes through the U-block, whose matrix depends on T
only through s_T (S3) and whose coupling to copy t depends on (w_t, s_T, w_U) only. Exact check (`check_struct.py`):
for u = 2 and u = 3, m = 4, changing (w_2, w_3) from (3/2, 1) to (3/5, 17/10) (same s_T, different multiset) leaves
G_00, G_01, G_11 unchanged in Q, while G_22 changes.

**Soft modes (numerical, `ex1.py`).** For s_T ≪ s_U the full block has exactly k − 1 small eigenvalues off the kernel,
each ≈ (0.8–0.9) s_T/s_U:
- u − 1 modes Ψ_a = Σ_{r,h∈U} a_r w_h P_{r,h}, a ⊥ w_U (the lift of Theorem T's second radical; before coupling,
  2⟨Ψ_a, YΨ_a⟩ = s_T ⟨a, E_0 E_T^{−1} a⟩ exactly, by 𝒯_{(r′,1)}(a_r w_U) = a_r E_0 e_{r′}; checked exactly in
  `check_soft.py`, k = 4, 5, 6);
- m modes Φ_t = Σ_h w_h (P_{t,h} + Z_{t,h}), one per t ∈ T (these are the soft modes of the u = 1 analysis).
The rescaling below removes the Φ_t; the Ψ_a stay inside (C_u).

### A3. The reduction

Fix r0 ∈ U. The kernel vectors p ⊗ f_h (Lemma N) have P_{r0,h}-coordinate w_{r0} ≠ 0, so the block is ⪰ 0 with kernel
exactly span{p ⊗ f_h} iff its restriction to {P_{r0,h} = 0 ∀h} is positive definite. Eliminate in the order X_T
(block I − w_Tw_T^T ⊗ E_T^{−1} ≻ 0, README), then the U-block. So the block is as claimed iff

- **(C_u)** the U-block of (S3), restricted to the complement of the kernel (equivalently: X_U ⊥ w_U c^T, or
  P_{r0,h} deleted), together with q, is positive definite. It depends only on (s_T, w_U), the same for all m; and
- **G ≻ 0**, G the Schur complement on the T-blocks.

For G use, per t, the coordinates (ρ_t, z_t, z′_t) with ζ_t = c z_t w_U/|w_U| + (z′_t ⊥ w_U), c = √(1 + s_U/s_T) (the
u = 1 rescaling ζ → cζ, applied to the w_U-component), and call the rescaled matrix Ĝ. By the weighted block
Gershgorin lemma of Proposition U1 (weights √s_t),

- **(A_u)** ‖Ĝ_tt′‖₂ ≤ H √(s_t s_t′)/s_T for all t ≠ t′, and
- **(B_u)** λ_min(Ĝ_tt) > H (1 − s_t/s_T) for all t (for m = 1: Ĝ_tt ≻ 0)

with one constant H imply G ≻ 0. For u = 1 these are exactly (A), (B) of Proposition U1 (proved there with H = 2/5).

### A4. Evidence for (C_u), (A_u), (B_u) (`sweep.py`, log `sweep.log`)

600 weight vectors per (u, m), u = 2…5, m = 1…4 (k ≤ 8): log-uniform with spreads up to e^{±6} in w, families with
s_T or s_U shrunk by up to e^{−12}, or one weight blown up by up to e^{12}. Coordinates for (C_u): X′ (X′w_U = 0,
w_U^TX′ = 0) and a ⊗ w_U/|w_U| (a ⊥ w_U), both orthonormal, and q. Results (`sweep.log`, PASS, 18 s, < 0.2 GB):

| quantity (worst over all (u, m) and samples) | value |
|---|---|
| \|Y(p ⊗ f_h)\| (relative) | ≤ 2e-11 |
| smallest eigenvalue of the full block off the kernel | > 0 down to roundoff (2e-17 at s_T/\|s\| ~ 1e-15) |
| λ_min(C_u) / (s_T/\|s\|), samples with s_T/\|s\| > 1e-8 | ≥ 0.500 |
| sup ‖Ĝ_tt′‖₂ s_T/√(s_t s_t′) | 0.2778 (u = 1: 5/18 = 0.27778) |
| inf λ_min(Ĝ_tt)/(1 − s_t/s_T), m ≥ 2 | 0.7500 (u = 1: 3/4) |
| weighted dominance ratio max_t Σ_t′ ‖Ĝ_tt′‖ √(s_t′/s_t)/λ_min(Ĝ_tt) | ≤ 0.2928 (u = 1: 0.2902) |
| negative control: the same ratio without the rescaling c (w_T shrunk by 10^{−3}) | 2·10^4 – 10^5 (fails, as it must) |

So numerically (A_u) and (B_u) hold with H = 2/5 for every tested u, with the same margins as u = 1, and (C_u)
suggests the clean bound

    (C_u′)  U-block ⪰ ½ (s_T/|s|) I on the complement of the kernel (orthonormal coordinates above).

The constant ½ is approached as s_T/|s| → 1 (`ex7.py`); as s_T → 0 the ratio tends to 8/9 (u = 2) and 15/16 (u = 3)
at equal U-weights and to 1 when one U-weight dominates, i.e. the soft modes Ψ_a have Rayleigh quotient
≈ (0.89–1) s_T/s_U. (Values below s_T/|s| ~ 1e-10 are at roundoff level and are not used.)

### A5. What remains

- (C_u) for every u (a statement about the virtual system k = u + 1, one T-copy of weight √s_T). Its soft modes are the
  Ψ_a; a proof needs a Theorem T-type argument with E_T = s_T I + E_0 (crossed-pair bound on the hard part) plus an
  exact treatment of the Ψ_a and their coupling to the column-sum direction μ = X_U^T w_U (in the gauge
  w_U^T X_U + w_τ X_τ = 0 the virtual row is X_τ = −μ/w_τ and contributes |μ|²/s_T, which makes μ stiff).
- (A_u) and (B_u) for u ≥ 2: blocks of size u + 1 depending on all of w_U; numerically they have the u = 1 extreme
  constants, so H = 2/5 should work, but no closed form uniform in u was obtained. For u = 2 the blocks depend on four
  normalised parameters, so the u = 1 route (closed forms + Bernstein certificates) is possible in principle but
  heavier (not attempted).

### Files of Part A

| File | Content |
|---|---|
| `cf.py` | closed form of the type I block (any number type) |
| `check_cf.py`, `check_cf.log` | closed form = Lemma B builder (exact), = direct construction (float); PASS |
| `check_struct.py`, `check_struct.log` | (S1)–(S4) exactly; PASS |
| `sweep.py`, `sweep.log` | numerical evidence for (C_u), (A_u), (B_u), full-block kernel, negative control; PASS |
| `check_soft.py`, `check_soft.log` | the exact Rayleigh value of the soft modes Ψ_a; PASS |
| `ex1.py`–`ex7.py` | exploration (soft modes, comparison matrices, bath test, location of the constant ½ in (C_u′)); no proof rests on them |
| `fastdirect.py`, `lemmab.py`, `bern.py` | shared with `../typeI_u1/` |

## Part B. The condition (C_u) for every u (Theorem C)

Notation as in `relative-anchors/README.md` (Lemmas B, N; Theorem T;
Remark after Proposition U0; Proposition U1) and Part A. U = copies without a letter 2, u = |U| ≥ 2; T = copies carrying one
letter 2, m = |T| ≥ 1; s_x = w_x²; S = Σ_{h∈U} s_h; σ = s_T = Σ_{t∈T} s_t > 0; K = σ + S/2. All forms are in units
of 2Y (as in Part A). Vectors a, c, w, … are indexed by U.

### Status of Part B

| Statement | Status |
|---|---|
| Lemma R: (C_u) ⟺ the principal block of the type I block on {P_{r,h}: r, h ∈ U} ∪ {Q} is positive definite | **proved** (Section B1; uses (S1), (S2) of Part A and Lemma N) |
| Closed form F(X) of that block with q eliminated; it depends on (σ, w_U) only | **derived**; exact agreement with the closed form `cf.py` of the whole block (7 cases, k ≤ 6) |
| Splitting X = X′ + a wᵀ + w cᵀ with the exact formulas (F1), (F2), (F3) | **derived**; exact checks u = 2…5 |
| Lemma H (hard part): F(X′) ≥ (5/12) Σ_{r≠h} X′_rh² + Σ_h (3/2 − ρ_h) X′_hh² | **proved** (Theorem T's crossed-pair bound) |
| Lemma S: F(X) ≥ σ Φ1(a, c), Φ1 an explicit quadratic form on R^{2u} | **proved**; exact checks (block formulas, the chain at random points) |
| Lemma A: the a-block of Φ1 is ⪰ (1 − t) I, t = σ/K | **proved** (by hand) |
| Conditions (I), (II), (III): three families of rational inequalities in three parameters | **proved** by exact Bernstein certificates (every coefficient positive on the root box; (III) with the corner rule at one vertex) |
| **Theorem C.** (C_u) holds for every u ≥ 1, every m ≥ 1 and all positive weights | **proved** (computer-assisted only through (I)–(III)) |
| Corollary Q: F(X) ≥ 9·10⁻⁵ (s_T/\|s\|) ‖X‖², uniformly in u (constant poor; numerically ≈ 0.42) | **proved** (same certificates with margins) |
| Hence, for every k and T: the type I block is ⪰ 0 with kernel exactly span{p ⊗ f_h : h ∈ U} ⟺ G ≻ 0 | **proved** (the reduction of Section A3, plus Theorem C) |
| G ≻ 0, i.e. (A_u) and (B_u), for u ≥ 2 | **open**; a candidate route (relaxing the trace-free hard part of the whole block) works numerically, Section B5 |

So one of the two halves of the u ≥ 2 problem is now closed for all u. The type I blocks with 1 ≤ |T| ≤ k − 2 are
**not** proved: the statement G ≻ 0 (the T-blocks) remains. Nothing here says anything about Conjecture 12 or the
endpoint problem beyond this block.

Inputs taken from earlier work and not re-derived here: the closed form of the whole type I block (`cf.py`, Section
A1, derived from Lemmas B and N; its exact checks are rerun), the structural facts (S1), (S2) (Part A,
rerun in `check_struct.py`), the null vectors p ⊗ f_h (Lemma N) and the reduction "block OK ⟺ (C_u) and G ≻ 0"
(Section A3). Exact point checks of Theorem C itself (`check_points.py`, 33 rational configurations,
u ≤ 5, m ≤ 3) and of the whole chain of Section B4 (`check_chain.py`) are sanity checks; the proof does not rest on them.

### B1. Lemma R: (C_u) is a statement about a principal block

Let M be the principal block of the type I block on the variables (X_U, q, X_T) (states P_{r,h}, r ∈ [k], h ∈ U, and
Q; ρ = ζ = 0). In Section A3, (C_u) says: after eliminating X_T (whose block is I − w_Tw_Tᵀ ⊗ E_T^{−1} ≻ 0),
the Schur complement on (X_U, q) is positive definite on the complement of {X_U = w_U cᵀ}. Equivalently (Schur
complement with a positive definite eliminated block): M ⪰ 0 with kernel exactly {κ_c : c ∈ R^U}, where κ_c = p ⊗
Σ_h c_h f_h has X_U = w_U cᵀ, X_T = w_T cᵀ, q = 0 (these are null vectors of the whole block, Lemma N, hence of M).

**Lemma R.** (C_u) holds iff the principal block of M on (X_U, q) (that is, X_T = 0) is positive definite.

*Proof.* Write X_T = w_T ηᵀ + X′_T with w_Tᵀ X′_T = 0. By (S1), (S2) of Part A (proved from the closed form, exact
checks rerun in `check_struct.py`), X′_T does not couple to (X_U, q), and the X_T block acts on X′_T as the identity and
does not couple X′_T to w_T ηᵀ (because w_Tᵀ X′_T = 0). Hence M(v) = M(X_U, q, w_T ηᵀ) + |X′_T|². Subtracting the null
vector κ_η: M(v) = M(X_U − w_U ηᵀ, q, 0) + |X′_T|². If the (X_U, q) block is positive definite this is ≥ 0 and vanishes
only if X_U = w_U ηᵀ, q = 0, X′_T = 0, i.e. v = κ_η. Conversely, if M ⪰ 0 with kernel {κ_c}, then M is positive definite
on {X_T = 0}, which meets the kernel only in 0 (w_T ≠ 0). ∎

So (C_u) depends only on (σ, w_U): the block on (X_U, q) involves T only through σ = s_T (below). `check_cu_form.py`
checks exactly that the formula below equals this principal block of `cf.TypeI` (k = 3…6, |T| = 1, 2, 3).

### B2. The form and its splitting

**Closed form** (the closed form of Section A1, at X_T = 0, ρ = ζ = 0). For X ∈ R^{U×U} with rows ψ_r,

    Φ(X, q) = ½‖X‖² + Σ_h X_hh² + ½(tr X)² + 2q(−tr X + wᵀXw/(2K)) + 2q² − 𝔅(X),
    𝔅(X) = Σ_{r,r′∈U} ⟨𝒯_{r′}ψ_r, E^{−1} 𝒯_rψ_{r′}⟩,   𝒯_{r′}ψ = ½(w·ψ) e_{r′} + w_{r′}(½ψ + ψ_{r′} e_{r′}),

with E = diag(d) + ½wwᵀ, d_h = s_h + S/2 + σ (this is E_T = σI + E_0). Since Φ_qq = 2, eliminating q gives
F(X) = Φ(X, 0) − ½(tr X − wᵀXw/(2K))². (C_u) ⟺ F ≻ 0 on R^{U×U}.

Put ρ_h = s_h/d_h, g = Σ_h ρ_h, η = w/d (entrywise), z = E^{−1}w, E_0 = E − σI.

**Splitting.** Every X is uniquely X = X′ + a wᵀ + w cᵀ with X′w = 0, wᵀX′ = 0, a ⊥ w, c ∈ R^U
(c = Xᵀw/S, a = Π(Xw)/S, X′ = ΠXΠ, Π = I − wwᵀ/S). With Ω = a wᵀ + w cᵀ and p = w·c:

- (F1) F(X′) = ½ Σ_{r≠h} X′_rh² + Σ_h (3/2 − ρ_h) X′_hh² + (1/(2+g)) Σ_{r,r′} ρ_r ρ_{r′} X′_{rr′} X′_{r′r}
  (the form of Theorem T, step 5, with d_h shifted by σ).
- (F2) F(Ω) = σ f(a, c), f = aᵀE_0E^{−1}a + cᵀE_0E^{−1}c − σp²/(2K²) + p (a·z) + 2 Σ_r a_r c_r w_r z_r.
  (For c = 0 this is the soft-mode value s_T⟨a, E_0E_T^{−1}a⟩ of Part A.)
- (F3) F(X′ + Ω) − F(X′) − F(Ω) = σ ⟨𝒢, X′⟩ (Frobenius), 𝒢 = 2E^{−1}diag(a∘w) + 2diag(c∘w)E^{−1} + p(E^{−1} − I/K).

Derivation. With v = Xw, μ = Xᵀw, β = wᵀXw and α_{rr′} = ½v_r + w_{r′}X_{rr′}, the B-term expands to
𝔅(X) = ¼μᵀE^{−1}μ + Σ_{r,h}(½β + w_rμ_r)E^{−1}_{rh}X_{rh} + Σ_{r,r′}E^{−1}_{rr′}α_{rr′}α_{r′r}, and
E^{−1} = D^{−1} − ηηᵀ/(2 + g), D = diag(d).
- (F1): for X′, v = μ = β = 0 and L = −tr X′; the last sum gives Σ_r ρ_r X′_rr² − (1/(2+g))Σρ_rρ_{r′}X′_{rr′}X′_{r′r}.
- (F2): by Lemma N, 𝒯_{r′}(a_r w + w_r c) = a_r E_0 e_{r′} + w_r P_{r′} with P_{r′} = 𝒯_{r′}c and Σ_{r′} w_{r′}P_{r′} = E_0c (on
  these strings Σ_{r∈U} w_r𝒯_r = E_0). Hence 𝔅(Ω) = aᵀE_0E^{−1}E_0a + 2Σ_r a_r⟨E^{−1}E_0w, P_r⟩ + cᵀE_0E^{−1}E_0c. The
  A-part is aᵀE_0a + cᵀE_0c + 2Σ_h s_h a_h c_h (a ⊥ w), and L(Ω) = −σp/K. Then use E_0 − E_0E^{−1}E_0 = σE_0E^{−1} and
  E^{−1}E_0w = w − σz.
- (F3): the A-part gives 2Σ_h X′_hh w_h(a_h + c_h) + p tr X′, the q-part −(σ/K) p tr X′, and the B-term
  2Σ_{r,r′}⟨𝒯_{r′}ψ^Ω_r, E^{−1}w_rN_rψ′_{r′}⟩ (𝒯_rψ′ = w_rN_rψ′ since X′w = 0; Σ_r s_rN_r = diag(d) − σI) equals
  2Σ_r a_rw_r(X′_rr − σ(E^{−1}X′)_rr) + p(tr X′ − σ⟨E^{−1}, X′⟩) + 2Σ_r w_rc_r(X′_rr − σ(X′E^{−1})_rr).
  Everything except the σ-terms cancels.

Exact checks (Fractions, random rational X′, a, c; u = 2…5): `check_cu_form.py`.

**Lemma H.** F(X′) ≥ (5/12) Σ_{r≠h} X′_rh² + Σ_h Λ_h X′_hh², Λ_h = 3/2 − ρ_h ≥ 5/6.

*Proof.* 2|X′_{rr′}X′_{r′r}| ≤ X′²_{rr′} + X′²_{r′r} gives F(X′) ≥ Σ_{r≠h}(½ − ρ_rρ_h/(2+g))X′_rh² + Σ_h Λ_h X′_hh².
Since g ≥ ρ_r + ρ_h, (r₁, r₂) ↦ r₁r₂/(2 + r₁ + r₂) is increasing, and ρ_h ≤ s_h/(s_h + S/2) = 2x/(1 + 2x) with
x = s_h/S, Theorem T's scalar lemma ((1+2x)(1+2y)(2 + ρ_a + ρ_b − 12ρ_aρ_b) = 2 + 6(x+y) − 32xy ≥ 0 for x + y ≤ 1) gives
ρ_rρ_h/(2+g) ≤ 1/12. Also ρ_h < 2/3. ∎ (`check_hardpart.py`.)

**Lemma S.** For a ⊥ w, put

    Φ1(a, c) = f(a, c) − (σ/4) [ (12/5) Σ_{r≠h} 𝒢_rh² + Σ_h 𝒢_hh² / Λ_h ].

Then F(X′ + awᵀ + wcᵀ) ≥ σ Φ1(a, c). If Φ1 ≻ 0 on {a ⊥ w} × R^U, then F ≻ 0.

*Proof.* By (F1)–(F3) and Lemma H, F(X) ≥ Σ_{r,h} λ_rh X′_rh² + σ Σ_{r,h} 𝒢_rh X′_rh + σf with λ_rh = 5/12 (r ≠ h),
λ_hh = Λ_h. Minimising each entry separately over all real X′_rh (a relaxation of X′w = 0, wᵀX′ = 0) gives
−(σ²/4) Σ 𝒢_rh²/λ_rh. If (a, c) ≠ 0 then F(X) ≥ σΦ1(a, c) > 0; if (a, c) = 0 then F(X) = F(X′) ≥ (5/12)‖X′‖². ∎

Numerically (`explore2.py`, `explore3.py`) this relaxation keeps a lot of room: the exact Schur complement onto the soft
space is ⪰ 0.75 f, the relaxed one (Φ1) ⪰ 0.556 f, the worst case being σ/S → ∞ with one dominant U-weight (5/9).

### B3. Normalised formulas for Φ1

Everything is homogeneous of degree 0 in (s_h, σ) (with a, c fixed), so normalise **K = 1**. Then t := σ ∈ (0, 1),
S = 2(1 − t), x_h := s_h ∈ (0, S], w_h = √x_h, d_h = 1 + x_h, ρ_h = x_h/d_h, g = Σρ_h ≤ S, G = 2 + g,
η_h² = x_h/d_h² = ρ_h/d_h, N = Σ η_h², Λ_h = 3/2 − ρ_h, e_h = 1/d_h − η_h²/G (= E^{−1}_hh > 0),
m_h = ρ_h + η_h²/G, β = 12/5, b_h = β η_h² (N − η_h²)/G², B = Σ b_h, Π_m = Σ m_h²/Λ_h, ν_h = ρ_h w_h.

On a ⊥ w (where α := Σρ_h w_h a_h = −K(a·η)):

    f = Σ_h [(1 − t + tρ_h)(a_h² + c_h²) + (4ρ_h/G) a_h c_h] + (t/G)[α² + γ² − 2pγ − (g/2)p²] − (2/G) p α,
    γ = Σ ρ_h w_h c_h,  p = Σ w_h c_h;
    𝒢_rh = −η_rη_h (2w_h a_h + 2w_r c_r + p)/G  (r ≠ h),   𝒢_hh = 2w_h e_h (a_h + c_h) − p m_h.

We use this expression of Φ1 as a quadratic form on all of R^{2u} (it agrees with Lemma S on a ⊥ w), and write
Φ1 = aᵀM_a a + cᵀM_c c + 2cᵀM_ca a with

    M_a = (1 − t)I + t diag(ρ) + (t/G) ννᵀ − t diag(β ρ_h² (N − η_h²)/G² + x_h e_h²/Λ_h),
    (M_c)_hh = 1 − t + tρ_h + tη_h²/G − t x_h/2 − (t/4) x_h [8b_h + B + 4e_h(e_h − m_h)/Λ_h + Π_m],
    (M_c)_hk = t w_h w_k π_hk,  π_hk = 1/(G d_h d_k) − 1/2 − ¼[B + 2b_h + 2b_k + Π_m − 2e_hm_h/Λ_h − 2e_km_k/Λ_k],
    (M_ca)_hh = (2 − x_h)ρ_h/G − t x_h e_h (e_h − m_h/2)/Λ_h − (tβ/(2G²)) x_h η_h² (N − η_h²),
    (M_ca)_hk = w_h w_k σ_hk,  σ_hk = −ρ_k/G + (t/2) e_k m_k/Λ_k − (tβ/G²) η_k² [η_h² + (N − η_k²)/2]   (h ≠ k).

`check_blocks.py` checks exactly (Fractions) that this block form equals Φ1 of Lemma S built from (F2), (F3) on
a ⊥ w (u = 2…5; the normalised products w_h w_k/K are rational for rational weights), and that F(X) − σΦ1(a, c) ≥ 0 at
random rational X. The entries were also compared in floating point with a direct polarisation of the form
(`explore10.py`, relative error 3e−15).

### B4. Theorem C: Φ1 ≻ 0, hence (C_u)

**Theorem C.** For every u ≥ 1, every σ > 0 and all positive U-weights, Φ1 ≻ 0 on R^{2u}. Hence F ≻ 0 (Lemma S), and
by Lemma R the condition (C_u) holds for every u, every T ≠ ∅ and all positive weights.

*Proof.* Three estimates and a per-coordinate reduction.

**(i) Lemma A: M_a ⪰ (1 − t) I.** M_a minus the positive term (t/G)ννᵀ is diagonal. Per coordinate,
x_h e_h²/Λ_h ≤ ρ_h/(d_hΛ_h) = ρ_h/(3/2 + x_h/2) ≤ (2/3)ρ_h (0 < e_h ≤ 1/d_h), and, with N − η_h² = Σ_{k≠h} ρ_k/d_k ≤ g − ρ_h,
β ρ_h² (N − η_h²)/G² ≤ (12/5) ρ_h · (g²/4)/G² ≤ (3/20) ρ_h (g < 2). So the diagonal is ≥ 1 − t + (11/60) tρ_h. ∎

**(ii) Signed weighted Gershgorin for M_c.** If π_hk ≤ 0 for all h ≠ k, then (2|c_hc_k| ≤ (w_k/w_h)c_h² + (w_h/w_k)c_k²)

    cᵀM_c c ≥ Σ_h c_h² [ (M_c)_hh + t Σ_{k≠h} x_k π_hk ] =: (1 − t) Σ_h LB_h c_h².

π_hk ≤ 0 holds: since B, b ≥ 0 and Π_m ≥ m_h²/Λ_h + m_k²/Λ_k,
π_hk ≤ Q(x_h, x_k, g) := 1/(G d_h d_k) − 1/2 + ¼[m_h(2e_h − m_h)/Λ_h + m_k(2e_k − m_k)/Λ_k], and

- **(III)** Q(x_h, x_k, g) ≤ 0 for x_h, x_k ≥ 0, x_h + x_k ≤ 2, ρ(x_h) + ρ(x_k) ≤ g ≤ 2 (ρ(y) = y/(1 + y)).

The signed row is explicit: with S′ = S − x_h, g′ = g − ρ_h = Σ_{k≠h} ρ_k, Σ_{k≠h} x_k/(G d_h d_k) = g′/(G d_h), so

    (1 − t) LB_h = (1 − t)² + tρ_h + t g/(G d_h)
                   − (t/4)[ S(B + Π_m) + (8x_h + 2S′) b_h + 4x_h e_h(e_h − m_h)/Λ_h − 2S′ e_hm_h/Λ_h
                            + Σ_{k≠h} x_k (2b_k − 2e_km_k/Λ_k) ].

**(iii) A Schur test for M_ca.** For any τ > 0 (2|c_ha_k| ≤ τ(w_k/w_h)c_h² + τ^{−1}(w_h/w_k)a_k²),

    2|cᵀM_ca a| ≤ (1 − t)[ τ Σ_h ROW_h c_h² + τ^{−1} Σ_k COL_k a_k² ],
    (1 − t) ROW_h = |(M_ca)_hh| + Σ_{k≠h} x_k |σ_hk|,   (1 − t) COL_h = |(M_ca)_hh| + Σ_{k≠h} x_k |σ_kh|.

Hence Φ1 ≥ (1 − t)[ Σ_k a_k² (1 − COL_k/τ) + Σ_h c_h² (LB_h − τ ROW_h) ], and Φ1 ≻ 0 follows from

- **(I)** LB_h − (4/5) ROW_h > 0 and **(II)** COL_h < 4/5, for every coordinate h of every configuration.

**(iv) Reduction to three parameters.** The global sums are bounded through quantities of the coordinate h and of the
rest R = U ∖ {h}, using only these facts (k ranges over R):

- M1. Σ_k x_kρ_k = S′ − g′ =: Δ′ (identity x − ρ(x) = xρ(x)).
- M2. ρ(S′) ≤ g′ ≤ S′ (ρ concave with ρ(0) = 0, hence subadditive; ρ(y) ≤ y).
- M3. Δ′/(1 + S′) ≤ Σ_k ρ_k² ≤ Δ′ (ρ_k = x_k/d_k, 1 ≤ d_k ≤ 1 + S′).
- M4. N − η_h² = Σ_k ρ_k/d_k = g′ − Σ_kρ_k² ∈ [0, N′], N′ := g′ − Δ′/(1 + S′).
- M5. B = (β/G²)(N² − Ση⁴) ≤ β(η_h² + N′)²/G²; Σ_k x_k b_k ≤ (β/G²)(η_h² + N′) Δ′ (x_kη_k² = ρ_k²).
- M6. Σ_k m_k²/Λ_k ≤ (2/3)(1 + 1/G)² Δ′ (m_k ≤ ρ_k(1 + 1/G), d_kΛ_k ≥ 3/2).
- M7. 0 ≤ (t/2)e_km_k/Λ_k ≤ ρ_k/G (uses g ≤ S = 2(1 − t): t(3 + g) ≤ t(5 − 2t) ≤ 3), so
  |σ_hk| ≤ ρ_k/G + (tβ/G²) η_k² [η_h² + (N − η_k²)/2].

With these: LB_h ≥ (the signed row with B, Π_m, b_h, Σ x_kb_k replaced by their bounds and −Σ x_ke_km_k/Λ_k dropped),
(1 − t)ROW_h ≤ |(M_ca)_hh| + Δ′/G + (tβ/G²)(η_h² + (η_h² + N′)/2)Δ′, and
(1 − t)COL_h ≤ |(M_ca)_hh| + S′ρ_h/G + (tβ/G²) η_h² (Δ′ + S′N′/2). The exact expressions are in `rowbounds2.py`
(function `bounds`). They depend only on the parameters

    S = 2s₀ ∈ [0, 2],  ξ = x_h/S ∈ [0, 1],  θ ∈ [0, 1] with g′ = ρ(S′) + θ(S′ − ρ(S′)),  ν ∈ [0, 1] with N − η_h² = νN′,

and are affine in ν, so ν ∈ {0, 1} suffices; |(M_ca)_hh| is handled by the two signs. `check_bounds_valid.py` checks
exactly, at 200 random rational configurations (u = 2…7, spread and one-dominant weights, S from 1/2000 to 1.99), that
every parameter lies in its range, that π_hk ≤ Q ≤ 0, M7, Lemma A's constant, and that every bound is valid.

**(v) Certificates (`build_conditions.py`, `prove_conditions.py`).** Each of (I) (two signs, ν = 0, 1), (II) (two signs,
ν = 0, 1) is a rational function of (s₀, ξ, θ) ∈ [0, 1]³; (III) is a rational function of (a, b, c) ∈ [0, 1]³ via
x_h = 2a, x_k = 2(1 − a)b, g = ρ_h + ρ_k + c(2 − ρ_h − ρ_k). SymPy forms numerator/denominator (numerators of degrees
(13, 12, 3) for (I), (7, 7, 2) for (II), (8, 4, 2) for (III); consistency with a direct Fraction evaluation is checked at
random rational points). Every denominator factor has all Bernstein coefficients positive on [0, 1]³; every numerator
has all tensor Bernstein coefficients positive on the root box (no subdivision needed), except that for (III) the
coefficient at the vertex (0, 0, 0) is 0, where −Q = 0 (x_h = x_k = g = 0, the limit t → 1), and the corner rule of
`bern.py` applies (Q < 0 elsewhere on the closed box). So (I), (II), (III) hold on the closed domains, which contain
every actual configuration (t ∈ (0, 1), M2, M4). Total time under 3 seconds.

Negative controls (`controls.py`): the same pipeline fails, with an explicit witness of the wrong sign, for (II) with
τ = 7/10 (sup COL ≈ 0.7116), for (I) with τ = 2, and for (III) with Q + 1/100. ∎

**Corollary Q (a uniform quantitative bound, poor constant).** For every u, all weights and all X ∈ R^{U×U},

    F(X) ≥ (25/267384) (σ/K) ‖X‖² ≥ 9·10⁻⁵ (s_T/|s|) ‖X‖²      (|s| = σ + S).

*Proof.* `prove_conditions.py` also certifies (I_q) LB − (4/5)ROW > 1/10 and (II_q) COL < 3/4 (same pipeline; (II_q)
needs one subdivision). Then Φ1 ≥ (1 − t)[(1/16)|a|² + (1/10)|c|²] ≥ ((1 − t)/16) ν², ν² = |a|² + |c|². Normalise K = 1
(F(X) is invariant under scaling all weights). For a ⊥ w, ‖Ω‖² = Sν² = 2(1 − t)ν² and ‖X‖² = ‖X′‖² + ‖Ω‖².
- (A) F(X) ≥ σΦ1 ≥ (t/32)‖Ω‖².
- (B) ‖𝒢‖_F ≤ √S(2|a| + 5|c|) ≤ √29 ‖Ω‖: ‖E^{−1}diag(a∘w)‖_F² = Σ x_r a_r² (E^{−2})_rr ≤ Σ ρ_r a_r² ≤ S|a|² (E ⪰ I), the
  same for diag(c∘w)E^{−1}, and |p| ‖E^{−1} − I‖_F ≤ √S|c| ‖diag(x) + ½wwᵀ‖_F ≤ (3/2)S^{3/2}|c| ≤ 3√S|c| (exact check in
  `check_blocks.py`: K²‖𝒢‖² ≤ 29 S ν² at random points). With Lemma H and f ≥ Φ1 ≥ 0,
  F(X) ≥ (5/12)‖X′‖² − t√29‖Ω‖‖X′‖ ≥ (5/24)‖X′‖² − (174/5)t²‖Ω‖².
- θ(A) + (1 − θ)(B) with θ = ½ if t ≤ 5/11136 and θ = 11136t/(5 + 11136t) otherwise gives a coefficient ≥ t/128 on ‖Ω‖²
  and ≥ min(5/48, 25/267384) on ‖X′‖². ∎

(Numerically the best constant is about 0.42 in the units σ/(σ + S), `sanity_cu.py`; the corollary only shows that a
bound of the form c·s_T/|s| holds uniformly in u.)

**Margins** (floating point, `explore16.py`, `sanity_cu.py`): min (LB − (4/5)ROW) ≈ 0.167, max COL ≈ 0.7116 (attained
for a single coordinate as t ≈ 0.93), min π-bound Q ≈ −5·10⁻⁵ near the corner t → 1. Over random weights
(u ≤ 8, σ/S over 14 decades): λ_min(Φ1) ≥ 0.667(1 − t) and λ_min of the (C_u) principal block ≥ 0.42 σ/(σ + S)
(the constant ½ of Part A's (C_u′) refers to a different gauge and normalisation and is not claimed here).

### B5. Consequence and what remains

With Theorem C, the reduction of Section A3 reads: for every k and every T with 1 ≤ m ≤ k − 1, the type I
block of (R∞)_{k−1} is ⪰ 0 with kernel exactly span{p ⊗ f_h : h ∈ U} **iff G ≻ 0**, where G is the Schur complement on
the T-blocks (ρ_t, ζ_t) after eliminating X_T, X_U (mod kernel) and q. For u = 1, G ≻ 0 is Proposition U1. For u ≥ 2,
G ≻ 0 is open; Part A's route (weighted block Gershgorin with (A_u), (B_u), H = 2/5) is numerically sound (u = 2…5,
same extreme constants as u = 1) but needs closed forms of the blocks Ĝ_tt, Ĝ_tt′, which involve the inverse of the
(C_u) block; no closed form uniform in u was obtained here (not attempted in the time).

A possible route for G ≻ 0 (numerical only, `explore17.py`, `explore18.py`): apply the relaxation of Lemma S to the
**whole** type I block, i.e. bound the hard U-part X′ by the Λ-diagonal form of Lemma H and eliminate it.
- Relaxing all of X′ fails for u ≥ 3 (eigenvalues of order −10³ s_T/|s|): X′ couples to ρ at order 1 through
  ½(tr X + Σρ_t)².
- Keeping the trace direction X′ = τ Π/(u − 1) exactly and relaxing only the trace-free part leaves a form on
  (τ, a, c, X_T, ρ, ζ, q) that is ⪰ 0 with kernel exactly span{p ⊗ f_h} in all samples (u = 2…5, m = 1…3, 75 weight
  vectors per u; smallest eigenvalue off the kernel ≥ 5·10⁻⁴ s_T/|s| in these non-orthonormal coordinates).
  That form has the same per-coordinate-plus-global structure as Φ1, so the method of Section B4 (with T-coordinates
  added) may extend; this was only tested with the subspace minimisation, not with an explicit coupling matrix.

Routes tried here that did not work:
- An exact decoupling of the soft modes in the virtual-copy gauge w̃ᵀX̃ = 0 (as in Theorem T's second radical): false
  (`explore4.py`, cross terms of order 0.1).
- The "three norms" condition λ_min(M_a)λ_min(M_c) ≥ ‖M_ca‖² with global λ_min(M_c) from the Gershgorin bound: fails
  narrowly after the moment bounds (0.7116² > 0.497); replaced by the per-row Schur test (iii).
- The first moment bounds (Σ_k ρ_k² ≤ ρ(S′)g′, the triangle bound for |π_hk|) made LB negative (−0.21) in the
  one-dominant-rest regime; replaced by the signed Gershgorin (π ≤ 0) and the identities M1, M3, M4.

### Files of Part B

| File | Content |
|---|---|
| `run_all.sh` | reruns every check below, one Python process at a time (each < 0.5 GB), prints PASS/FAIL per script and ALL PASS; logs in `logs/` (about 1 minute) |
| `cf.py`, `lemmab.py`, `fastdirect.py`, `check_cf.py`, `check_struct.py` | from Part A: closed form of the whole type I block and its exact checks, (S1)–(S4) |
| `bern.py`, `test_bern.py` | exact tensor Bernstein certificates (copy of the u = 1 module) and its unit tests |
| `cu_form.py` | the (C_u) form: F(X) and the principal block on (X_U, q) |
| `check_cu_form.py` | principal block of cf.TypeI = cu_form (exact, k ≤ 6), q-elimination, (F1)–(F3) exactly (u = 2…5) |
| `check_hardpart.py` | the scalar lemma of Lemma H and the identity d_hΛ_h = 3/2 + x_h/2 (SymPy) |
| `phi1.py`, `check_blocks.py` | Φ1 from (F2), (F3) and the normalised block formulas; exact equality, and F(X) ≥ σΦ1 at random points |
| `rowbounds2.py` | the per-coordinate bounds LB, ROW, COL of Section B4 and the pair function Q (generic number type) |
| `check_bounds_valid.py` | exact validity of every bound and auxiliary inequality at 200 configurations |
| `build_conditions.py`, `conditions.pkl` | conditions (I), (II), (III) as rational functions (regenerated by run_all.sh) |
| `prove_conditions.py` | the Bernstein certificates of (I), (II), (III), denominator signs, consistency checks |
| `controls.py` | negative controls (τ = 7/10, τ = 2, Q + 1/100) with witnesses |
| `check_chain.py` | exact: Φ1 − (1 − t)diag(1 − COL_k/τ, LB_h − τROW_h) ⪰ 0 (LDLᵀ) at 32 rational points, with the actual and with the bounded per-coordinate quantities |
| `check_points.py` | exact: the (C_u) principal block of cf.TypeI has positive LDLᵀ pivots at 33 rational configurations |
| `sanity_cu.py` | float sanity sweep (no claim rests on it) |
| `explore1.py` … `explore18.py`, `rowbounds.py` | exploration (numerical); no claim rests on them (`explore17.py`, `explore18.py`: the relaxation applied to the whole block, Section B5) |
| `logs/` | output of `run_all.sh` |
