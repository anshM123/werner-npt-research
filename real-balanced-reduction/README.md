# Real balanced projectors: an equivalent form of the endpoint

Ansh Mishra and Aryan Senthilkumar, 6 November 2026. Nothing here is refereed.

`support-and-projector-notes/` reduces the endpoint for every n to rank-two orthogonal projectors at odd n
(Proposition 4′) and to balanced planes, whose vectors all have maximally mixed one-copy marginals (Theorems 14–14″).
This note adds reality. Any negative rank-two witness becomes a real one of the same rank after an even number of added
copies, and the balancing can be done with real operators. So the endpoint for every n is equivalent to

    q_N(P) ≥ 0 for every odd N and every real rank-two orthogonal projector P whose range (and co-range) is balanced.

This is a reduction, not a positivity proof. It does not let one replace complex matrices by real ones at a fixed n.
The same reduction works at every fixed point t in (0, ½] of the Werner line, in particular at t = 2/5, where the
state is NPT (see "Below the endpoint", added 9 November 2026).

## Realification

Write b_n(X, Y) = ⟨X, M^{⊗n}(Y)⟩ and q_n(X) = b_n(X, X), with M(X) = X − Tr(X) I₃/2 on each copy. Let ω = e^{2πi/3},
f = (|00⟩ + ω|11⟩ + ω²|22⟩)/√3 and g = f̄: orthonormal maximally entangled two-qutrit vectors. In the basis
|f⟩⟨f|, |g⟩⟨g|, |f⟩⟨g|, |g⟩⟨f| their b₂-Gram matrix is

    (1/36) [[33, −3, 0, 0], [−3, 33, 0, 0], [0, 0, 24, 0], [0, 0, 0, 24]],

so a = 11/12 (diagonal), b = −1/12 (cross) and c = 2/3 (transitions). Tensoring r flag pairs gives a^r, b^r, c^r.

**Theorem R.** Let C = UΣV^† = Σ_j σ_j u_j v_j^† (rank k ≤ 2) on n copies, T = UΣV^T, F_r = f^{⊗r}, G_r = g^{⊗r},
R_r z = (z ⊗ F_r + z̄ ⊗ G_r)/√2 and C_r = Σ_j σ_j (R_r u_j)(R_r v_j)^T on n + 2r copies. Then C_r is real, has the same rank
and singular values as C, and

    q_{n+2r}(C_r) = ½ [ (11/12)^r q_n(C) + (−1/12)^r Re b_n(C, C̄) + (2/3)^r q_n(T) ].

*Proof.* R_r z has real coefficients and ⟨R_r x, R_r y⟩ = Re⟨x, y⟩, so the singular vectors stay orthonormal. Expand
C_r into its four flag blocks and use the Gram matrix above. ∎

**Corollary R.** Since ‖M^{⊗n}‖ = 1, |Re b_n(C, C̄)| and |q_n(T)| are at most ‖C‖². If q_n(C) = −d‖C‖² < 0, then
q_{n+2r}(C_r) ≤ ½(11/12)^r [−d + (1/11)^r + (8/11)^r] ‖C‖² < 0 as soon as (1/11)^r + (8/11)^r < d (for example
r > log(2/d)/log(11/8)). Parity is preserved. For positive semidefinite C (U = V) the output is positive semidefinite
with the same eigenvalues, so projectors stay projectors. A balanced input plane stays balanced for r ≥ 2: at the
original copies the orthogonality of F_r and G_r removes the cross branches, and at an added copy each diagonal branch
has marginal I₃/3.

## Real balancing with a uniform error

The twelve real orthogonal matrices S^t D_ε (S the cyclic shift of three coordinates, t = 0, 1, 2, and
ε ∈ {(1,1,1), (1,−1,−1), (−1,1,−1), (−1,−1,1)}) form a design: (1/12) Σ O^T H O = (Tr H/3) I₃. For N input copies use
the group 𝒢 = 𝒢₀^N (g₀ = 12^N elements), distinct binary labels of length m = ⌈log₂ g₀⌉, and the real orthonormal,
maximally entangled flags h₀ = (|00⟩ + |11⟩ + |22⟩)/√3 and h₁ = (|01⟩ + |12⟩ + |20⟩)/√3. For a group element O_j let H_j
be the product of the flags of its label, and V_ℓ x = g₀^{−1/2} Σ_j O_j x ⊗ H_j^{⊗ℓ} (ℓ ≥ 2), a real isometry whose whole
image is balanced.

**Theorem B.** For every real C on N copies, C̃ = V_ℓ C V_ℓ^T is real, has the same rank and singular values,
has balanced range and co-range, and

    | q_{N+2mℓ}(C̃) − (a_*^{mℓ}/g₀) q_N(C) | ≤ (a_*^{mℓ}/g₀)(g₀³ − 1)(8/11)^ℓ ‖C‖_F²,   a_* = 11/12.

*Proof.* At the endpoint the b₂-Gram matrix of the four flag operators of h₀, h₁ has diagonal entries a_* = 11/12, a
cross entry −1/12 and transition entries 2/3 and 0. Every entry other than the two all-equal diagonal ones has absolute
value at most (8/11)a_*. Exactly g₀ of the g₀⁴ quadruples of group indices are all equal. Each of them gives the
diagonal term, by the design and the covariance of M under real orthogonal conjugation. Every other quadruple differs
at some flag position and is suppressed by (8/11)^ℓ. ∎ So a negative real C stays negative after balancing once
(g₀³ − 1)(8/11)^ℓ is below its relative margin, and projectors stay projectors.

**Corollary (the equivalent target).** Start from a negative rank-two witness. Proposition 4′ of
`support-and-projector-notes/` gives a negative rank-two orthogonal projector with an odd number of copies. Theorem R
makes it real, and Theorem B balances it; both keep the projector property and add an even number of copies. So the
endpoint for every n is equivalent to the statement displayed at the top.

## Below the endpoint (added 9 November 2026)

The three steps above (dilation to a projector, realification, balancing) work at every fixed t with 0 < t ≤ ½. This
section is from an independent audit of this repository; we checked it (`code/below_endpoint_check.py`). Write
M_t(X) = X − t Tr(X) I₃, b_n(t; X, Y) = ⟨X, M_t^{⊗n}(Y)⟩, q_n(t, X) = b_n(t; X, X) and m = 1 − t. At t = ½ this is
the form above. M_t has eigenvalues 1 (traceless) and 1 − 3t (identity), so ‖M_t^{⊗n}‖ = 1 for 0 ≤ t ≤ 2/3, and
|b_n(t; X, Y)| ≤ ‖X‖_F‖Y‖_F.

*Dilation.* With the notation of Proposition 4′ of `support-and-projector-notes/` (C = Σ_i σ_i u_iv_i^*, λ = 2σ_1,
α_i + β_i = λ, α_iβ_i = σ_i², l_i = √α_i u_i, r_i = √β_i v_i, flags |0^k⟩ and |1^k⟩ on k added copies),

    q_{n+k}(t, P) = 2q_n(t, C) + m^k [q_n(t, D_0) + q_n(t, D_1)] + 2(−t)^k Re b_n(t; D_0, D_1),

because on one flag copy the off-diagonal units |0⟩⟨1|, |1⟩⟨0| have value 1, the diagonal units have value m, and the
two diagonal units pair to −t. Since D_0, D_1 ⪰ 0, t ≤ m and Tr D_0 + Tr D_1 = 2λ, the projector Π = P/λ satisfies
q_{n+k}(t, Π) ≤ 2q_n(t, C)/λ² + 4m^k. So a witness with q_n(t, C) = −δ‖C‖_F² gives a negative rank-two projector once
m^k < δ/8, with n + k of either parity. (At t = ½ the two factors m^k and (−t)^k combine into 2^{−k}q_n(D_0 ± D_1), as in
Proposition 4′; below the endpoint they must be kept apart.)

*Flags.* For the flag pairs (f, g) and (h₀, h₁) above, the b₂(t)-Gram matrix of the four flag operators is

    [[a, b, 0, 0], [b, a, 0, 0], [0, 0, c, 0], [0, 0, 0, c]],   a = 1 − 2t/3 + t²,   b = t² − 2t/3,   c = 1 − 2t/3,

checked exactly as polynomials in t. At t = ½ these are 11/12, −1/12 and 2/3. For 0 < t ≤ ½ we have a > c > |b|; the
gap a − c = t² is what drives the suppression below.

**Theorem R_t.** With C_r as in Theorem R,

    q_{n+2r}(t, C_r) = ½ [ a^r q_n(t, C) + b^r Re b_n(t; C, C̄) + c^r q_n(t, T) ].

So a negative witness stays negative as soon as (|b|/a)^r + (c/a)^r is below its relative margin. Rank, singular values
and the projector property are kept, and two copies are added at a time.

**Theorem B_t.** With the construction of Theorem B (the same design, labels and flags h₀, h₁; write L for the label
length, called m there), for every real C,

    | q_{N+2Lℓ}(t, C̃) − (a^{Lℓ}/g₀) q_N(t, C) | ≤ (a^{Lℓ}/g₀)(g₀³ − 1)(c/a)^ℓ ‖C‖_F².

The proof is the one of Theorem B: the g₀ quadruples with all four group indices equal give the main term, and every
other quadruple meets a label bit where the flag entry is b, c or 0, each at most (c/a)·a in absolute value.

**Corollary (the target at a fixed t).** For each fixed t in (0, ½], q_n(t, C) ≥ 0 for every n and every complex C of
rank at most two if and only if

    q_N(t, P) ≥ 0 for every odd N and every real rank-two orthogonal projector P whose range is balanced.

At t = 2/5, for example, a = 67/75, b = −8/75 and c = 11/15, and the statement for every odd N would make
(5I − 2F)/39, whose partial transpose has the eigenvalue −1/39, undistillable from any number of copies. That positivity
is not proved. Two remarks:
- On the target the two singular values are equal, so the sharp bound of Section 9 of `support-and-projector-notes/`
  reads q_N(t, P) ≥ 2(1 − 2t)(1 − t)^{N−1} there; at the endpoint it says only q_N(P) ≥ 0. In the two-copy sum-of-squares
  identity (the sum-of-squares lemma of `two-copy-thickness/`, after Fraser, Huber, Pozsgay and Vona) the singular-value variance term then
  vanishes in every gauge, and all of the slack sits in the contraction defects (the audit's observation).
- Realification also shows that the merging comparison q_d(t, C) ≥ ((1 − 2t)/2) q_{d′}(t, C) of Section 9 is only needed
  for real rank-two C: a complex negative witness at qutrit copies gives a real one, and merging keeps C real.

## Checks

- `code/below_endpoint_check.py` (`logs/below_endpoint_check.log`): the flag Gram matrices of (f, g) and (h₀, h₁) as
  exact polynomials in t (SymPy), and the design; the dilation identity, the projector property and the bound at
  t = 0.2, 0.35, 0.4, 0.5; Theorem R_t at the same t; the four-index expansion and the bound of Theorem B_t on a small
  instance (one input copy, a two-element group, one label bit, ℓ = 2, five copies in all); and the merging-telescope
  calibration of Section 9 of `support-and-projector-notes/` (floating point except where stated).
- `code/realification_exact.js` (`logs/realification_exact.log`), exact over Z[ω]: all sixteen flag Gram entries;
  direct contractions of the three-copy output for a complex non-Hermitian and a positive semidefinite rank-two example
  (rank and norm preserved, the identity exact); a balanced three-copy plane realified on seven copies, with every
  one-copy compression equal to I₂/3.
- `code/balancing_exact.js` (`logs/balancing_exact.log`), exact over Z[ω] and the integers: the flag Gram matrices of
  the balancing construction; the real orthogonal design; the erasure of cross labels. It also checks an interpolation
  family. For the maps M_z(X) = ((1 + z)/2) X + ((1 − z)/2) X^T − ½ Tr(X) I₃ (M₁ = M), 0 < z < 1, there are real rank-two
  projectors P_R on 1 + 2R copies with q < 0 (for example q_{11, 1/2}(P_5) = −330131/509607936), while at the endpoint
  z = 1 the same projector has q_{11}(P_5) = 49153/497664 > 0. So no argument that only uses the maps M_z with z < 1
  can reach the endpoint.
- `code/realify_check.py` (`logs/realify_check.log`), written separately: the identity of Theorem R from the
  definition of q_n (partial traces), for random complex rank-two C on n = 1, 2 copies with r = 1, 2 (floating point,
  agreement to 10⁻¹⁵), and C_r real, of rank two, with the same norm.

Run `bash code/run_all.sh` (Node, no packages; Python with NumPy).
