# Real balanced projectors: an equivalent form of the endpoint

Ansh Mishra and Aryan Senthilkumar, 6 November 2026. Nothing here is refereed.

`support-and-projector-notes/` reduces the endpoint for every n to rank-two orthogonal projectors at odd n
(Proposition 4′) and to balanced planes, whose vectors all have maximally mixed one-copy marginals (Theorems 14–14″).
This note adds reality. Any negative rank-two witness becomes a real one of the same rank after an even number of added
copies, and the balancing can be done with real operators. So the endpoint for every n is equivalent to

    q_N(P) ≥ 0 for every odd N and every real rank-two orthogonal projector P whose range (and co-range) is balanced.

This is a reduction, not a positivity proof. It does not let one replace complex matrices by real ones at a fixed n.

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

## Checks

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
