# Four copies at the qutrit Werner endpoint: Hermitian rank-two operators

Ansh Mishra and Aryan Senthilkumar

**What this proves.** For every Hermitian operator C of rank at most two on (C³)^⊗4,

    q_4(C) ≥ 0.        (Theorem 3)

In particular, q_4(P_u + P_v) ≥ 0 and q_4(P_u − P_v) ≥ 0 for every orthonormal pair u, v (Theorems 1 and 2).

**How.** Each proof is an exact rational dual certificate for a semidefinite relaxation. Each certificate is checked
in exact arithmetic by a short verification script that imports none of the code that produced it.

**What this does not prove.** It does not cover:
- non-Hermitian rank-two operators, so it is not the full four-copy endpoint;
- any number of copies other than four.

The proofs are computer-assisted. They have been checked independently by us, but they have not been refereed.

This file collects the statements, the definitions and the argument. The working notes are in `NOTES.md` (Section 3),
with the build logs in `code/exact4/BUILD_LOG.md`, `code/exact4m/BUILD_LOG.md` and `code/exact4L/BUILD_LOG.md`.

## 0. Setting

Let ρ = I − F/2 on C³ ⊗ C³. This is the qutrit Werner state at the endpoint α = 1/2 of its NPT, one-copy undistillable
range, unnormalised; F is the swap. For an operator C on H = (C³)^⊗n put

    q_n(C) = Σ_{J ⊆ {1..n}} (−1/2)^{|J|} ‖Tr_J C‖₂².

n-copy undistillability at the endpoint is equivalent to q_n(C) ≥ 0 for every C of rank at most two (Costa Rico;
Qi–Gupur–Wu–Guo; Song–Chen).

**What is known.**
- At two copies the inequality holds [FHPV, BGH].
- At three copies it is settled by Tomar's computer-assisted certificates [T1]. Wu and Zou prove it in every dimension
  for positive and for normal rank-two C [WZ].
- At four copies, Tomar's certificates reach α ≤ 0.4979, strictly below the endpoint [T1]. Clarisse reported numerical
  evidence at the endpoint [Cl].
- For every n the question is open.

A Hermitian C of rank two is a P_u + b P_v with u ⊥ v. Theorem 3 covers all real a, b. Theorems 1 and 2 are the
cases of equal modulus, a = b and a = −b; they were proved first, with separate certificates.

Writing W = ⊗_i (F_i − ½ I) on H ⊗ H, with F_i the swap of copy i, one has q_n(C) = Tr[(C ⊗ C†) W]. So

    q_n(a P_u + b P_v) = 2^{−n} (a² g(u,u) + b² g(v,v) + 2ab g(u,v)),   g(x,y) = ⟨x ⊗ y| ⊗_i (2F_i − I) |x ⊗ y⟩.

## 1. The relaxation

### Variables

Take four replicas of H. For a word c = c₁c₂c₃c₄ ∈ {u, v}⁴ put

    X_c = T(P_{c₁} ⊗ P_{c₂} ⊗ P_{c₃} ⊗ P_{c₄}),   T(Y) = ∫ (⊗_i U_i^⊗4) Y (⊗_i U_i^⊗4)† dU_1 … dU_n,

where U_i ∈ U(3) acts in the same way on the four replicas of copy i.
- By Schur–Weyl duality, T(Y) lies in the algebra spanned by the replica permutations, ⊗_i C[S_4]. Its coefficients
  are real when u and v are real.
- Since the replicas can be permuted, the five words uuuu, uuuv, uuvv, uvvv, vvvv suffice.
- The objective is linear in these variables:

      g(x, y) = Tr[X_{xy··} ((⊗_i (2F_i − I)) on replicas 1, 2)] =: f_xy.

  So q_4(aP_u + bP_v) = (a² f_uu + b² f_vv + 2ab f_uv)/16, and Theorems 1 and 2 concern
  c± = (f_uu + f_vv ± 2 f_uv)/16.

### Constraints

Every constraint below holds whenever u, v are orthonormal, i.e. at every *true point*.

- **(X) Positivity.** X_c ⪰ 0. A twirl of a positive operator is positive.
- **(M) Moments.** Let σ ∈ S_4 permute the replicas, on all copies at once. Then Tr[X_c σ] = Π over the cycles of σ of
  Tr(P_{c_{i₁}} P_{c_{i₂}} ⋯). This is 1 when each cycle carries a single letter, and 0 otherwise, because
  P_u P_v = 0.
- **(S) Pure letters.** X_c F_{ij} = X_c whenever c_i = c_j. Here F_{ij} swaps replicas i and j on every copy;
  F(x ⊗ x) = x ⊗ x, and T commutes with replica permutations.
- **(O) Orthogonality.** Tr_j[X_c F_{ij}] = 0 whenever c_i ≠ c_j, since Tr_2[(P_x ⊗ P_y) F] = P_x P_y = 0.
- **(G) Marginals.** Tr_4 X_{xyzu} = Tr_4 X_{xyzv}, because Tr P_u = Tr P_v = 1.
- **(C) Covariance.** For each x, y ∈ {u, v}:

      [[R(X_{xyuu}), R(X_{xyuv})], [R(X_{xyvu}), R(X_{xyvv})]] ⪰ 0,

  where R realigns replicas 3 and 4 on every copy, moving the bra of replica 3 to the ket side and the ket of replica 4
  to the bra side.
  - *Proof.* R(P_z ⊗ P_w) = |vec P_z⟩⟨vec P_w|. Hence the block matrix is the twirl of P_x ⊗ P_y ⊗ |η⟩⟨η|, with
    η = e₁ ⊗ vec P_u + e₂ ⊗ vec P_v.
  - After realignment, the twirl acts on replicas 3 and 4 as U ⊗ Ū, so it is still an average of unitary
    conjugations. The result is positive.

The certificates use only (X), (M), (S), (O), (G) and (C). The code also has partial-transpose blocks and
P_u + P_v ≤ I, but the certificates do not use them.

### Symmetry reduction and real parameters

- **Symmetries.** Averaging over permutations of the four copies maps true points to points with the same objective.
  So the variables can be taken symmetric under copy permutations. Together with (S), this leaves 1331 real parameters
  at n = 4.
- **Coordinates.** Each X_c is stored through its Fourier blocks. These use exact rational irreducible modules of
  C[S_4], with Gram matrices, built from integer Casimir eigenspaces.
- **Complex vectors.** For complex u, v the real part of the point x(u, v) is the average of x(u, v) and x(ū, v̄).
  Both are true points, and they have the same objective, since q_4 is invariant under complex conjugation. Every
  constraint map has rational coefficients. So the real part is feasible with the same objective, and real
  certificates cover complex pairs.

## 2. Certificates

### The form of a certificate

- **Data.** (S) is built into the parametrisation. The equalities (M), (O), (G) give 147 independent rows, E x = e.
  Each positivity constraint is W_k(x) ⪰ 0, with W_k linear.
- **Facial reduction.** Let N_k be a basis of the common kernel of W_k at the known zeros, i.e. true points with
  objective 0. Then W_k(x) ⪰ 0 implies F_k(x) := N_kᵀ W_k(x) N_k ⪰ 0, so the reduced constraints are still valid.
- **Certificate.** Rational y and rational S_k ⪰ 0 with

      c = Eᵀ y + Σ_k F_k*(S_k)   exactly, coordinate by coordinate,   and   e·y = 0.

- **Why it proves the theorem.** For every feasible x,

      c·x = y·E x + Σ_k Tr(S_k F_k(x)) = Σ_k Tr(S_k F_k(x)) ≥ 0.

  Every true point is feasible, so q_4 ≥ 0 there.
- **Why facial reduction is needed.** Without it no certificate is strictly feasible, because the zeros force every
  dual block to vanish on their ranges. With the full face the reduced dual is strictly feasible. A floating-point
  solution can then be rounded to a dyadic grid and corrected exactly, by solving the linear identity over Q.

### Theorem 1: P_u + P_v

Code: `code/exact4/`.
- Zeros: 11 exact integer zeros. These are products orthogonal on one copy or on all three remaining copies, rational
  rotations of those, and rotated interpolating zeros.
- Size after removing directions on which all data vanish: 256 blocks, the largest of size 33.
- Float margin: 0.248 (normalised).
- Exact correction: 1291 × 1291.
- Smallest LDLᵀ pivot: 9.4e-5.
- Entries: up to about 1060 bits.

### Theorem 2: P_u − P_v

Code: `code/exact4m/`.
- Zeros: 16 exact integer zeros: products orthogonal on two or on all four copies, the "inner" zeros of Section 3, and
  8 spin-flip zeros.
- Size: 272 blocks, the largest of size 32.
- Float margin: 0.09999.
- Exact correction: 1298 × 1298.
- Smallest pivot: 3.96e-4.
- Entries: up to 605 bits.
- The spin-flip zeros are essential: without them the face is strictly larger and the reduced dual is degenerate.

### Theorem 3: every Hermitian operator of rank at most two

**Theorem 3 (computer-assisted).** q_4(C) ≥ 0 for every Hermitian C of rank at most two on (C³)^⊗4. Equivalently,
q_4(a P_u + b P_v) ≥ 0 for all real a, b and every orthonormal pair u, v.

Code: `code/exact4L/` (build log `code/exact4L/BUILD_LOG.md`); logs: `logs/exact4L/`.

**One certificate for all ratios.** The certificate is matrix valued. It shows that the 2 × 2 matrix
[[f_uu·x, f_uv·x], [f_uv·x, f_vv·x]] is positive semidefinite on a lifted relaxation.
- Variables: Y₊, Y₋, Z ∈ R^1331, and ε₊, ε₋, ζ.
- Equalities: E Y₊ = e ε₊, E Y₋ = e ε₋, E Z = e ζ.
- Normalisation: [[ε₊, ζ], [ζ, ε₋]] ⪰ 0 and (ε₊ + ε₋)/2 = 1.
- For every block k:

      L_k = [[F₊,k(Y₊), G_k(Z)], [G_k(Z)ᵀ, F₋,k(Y₋)]] ⪰ 0.

  Here F±,k are the reduced constraints of Theorems 1 and 2, and G_k = N₊,kᵀ W_k N₋,k are the cross data between the two
  faces.
- Objective: ¼(c₊·Y₊ + c₋·Y₋ + 2d·Z), with d = (f_uu − f_vv)/16.

For an orthonormal pair, set x = x(u, v) and

    Y₊ = (a+b)² x,   Y₋ = (a−b)² x,   Z = (a² − b²) x,   ε₊ = (a+b)²,   ε₋ = (a−b)²,   ζ = a² − b².

This point is feasible, since L_k = MᵀW_k(x)M with M = [(a+b)N₊, (a−b)N₋]. Its objective is q_4(aP_u + bP_v).

**Certificate.** y₊, y₋, y_Z ∈ Q^147 and S_k = [[A_k, B_k], [B_kᵀ, C_k]] ⪰ 0 with, exactly,

    c₊/4 = Eᵀy₊ + Σ_k F₊,k*(A_k),   c₋/4 = Eᵀy₋ + Σ_k F₋,k*(C_k),   d/2 = Eᵀy_Z + Σ_k 2 G_k*(B_k),
    e·y₊ = e·y₋ = e·y_Z = 0.

Then the objective equals Σ_k Tr(S_k L_k) ≥ 0 on the lifted relaxation.

**Structure.** The cross data vanish at every zero of both kinds. Their ranges lie in the spans of the two
one-sign problems. So each lifted block lives on exactly the reduced spaces of Theorems 1 and 2 (1556 + 1564 = 3120
dimensions in all, 286 blocks, the largest of size 65).

**Numbers.**
- The float solve (Clarabel with faer) has smallest normalised lifted dual eigenvalue 9.6e-4.
- The exact correction is three square solves of size 1291, 1298 and 1301.
- Every S_k is PSD by exact LDLᵀ, with no zero pivots (smallest pivot 1.1e-5).
- Entries have up to about 1055 bits.

Theorem 3 contains Theorems 1 and 2 (the ratios (1, 1) and (1, −1)). The three certificates were built and checked
separately.

**What it does not cover.** Non-Hermitian rank-two C at four copies, which is the full four-copy endpoint, is still
open. Tomar's relaxation for general C gives −1.31e-3 there.

## 3. Zeros

**Contraction identity.** If u = a ⊗ φ and v = a′ ⊗ χ with a ⊥ a′ on one copy, then

    q_n(P_u + t P_v) = ½ q_{n−1}(P_φ − t P_χ).

Split the sum over J by whether that copy is traced. Iterating gives:
- **Products orthogonal on k copies, sharing any vector on the rest.** These are zeros for t = +1 when k is odd and for
  t = −1 when k is even.
- **Inner zeros of P_u − P_v.** u = a ⊗ φ, v = a′ ⊗ χ, with (φ, χ) any orthonormal basis of the range of a
  projector zero at n − 1 copies.

**Lemma (qubit subspaces).** Let Q_i ⊂ C³ be two-dimensional, and let C have range and co-range in Q_1 ⊗ ⋯ ⊗ Q_n.
Expand C = Σ_μ c_μ σ_μ in Pauli strings. Then

    q_n(C) = 2^n Σ_{μ ∈ {x,y,z}^n} |c_μ|².

*Proof.* Tr_J σ_μ = 0 unless μ_j = 0 for all j ∈ J, and then it is 2^{|J|} times the Pauli string on the remaining
copies. So ‖Tr_J C‖² = 2^n Σ_{μ: μ_J = 0} 2^{|J|} |c_μ|². In q_n the string μ therefore gets the weight
Σ_{J ⊆ Z(μ)} (−1)^{|J|}, where Z(μ) = {j : μ_j = 0}. This weight is 1 if Z(μ) = ∅ and 0 otherwise.

**Spin-flip zeros.** Let J = σ_y K be the spin flip on each Q_i and v = J^⊗n ū. The full-weight correlations of v are
(−1)^n times those of u. Hence:
- for even n, P_u − P_v is a zero whenever ⟨u, J^⊗n ū⟩ = 0;
- for odd n, P_u + P_v is a zero for every u.

At two copies, Tomar shows that every zero has range and co-range in a product of two-dimensional single-copy
subspaces [T2]. The lemma says exactly which operators on such products are zeros, for every n.

At four copies, 8 of 8 local searches for zeros of q_4(P_u − P_v) land on spin-flip pairs (`code/relaxations/tm1_*`).
Real integer points are easy to write down. With Y|x⟩ = (−1)^{w(x)} |x̄⟩ (σ_y^⊗4 in a real form), take u integer
with uᵀYu = 0 and v = Yu.

## 4. Checks

**Independent verification.** `verify4.py` and `verify4m.py` read only the saved data, the certificate and exact test
points. `verify4m.py` imports no pipeline code at all. They check, in exact rational arithmetic:
1. the identity c = Eᵀy + Σ_k F_k*(S_k) in every coordinate, and e·y = 0;
2. an exact LDLᵀ of every S_k;
3. at exact entangled true points:
   - q_4 is recomputed directly from the 81 × 81 matrix;
   - c·x(P) = q_4;
   - E x(P) = e;
   - the identity holds at the point with every term ≥ 0;
   - every reduced block F_k(x(P)) is PSD.

**Test points.**
- Theorem 1 passes at 8 points and Theorem 2 at 9. Of these, 6 and 7 respectively were generated after the
  certificates were built, from fresh seeds. They include near-zeros with q_4 = 2.4e-7 and 6.4e-5 (Theorem 1), and
  5.1e-6 and 1.8e-8 (Theorem 2).
- Theorem 3 is checked by `verifyL.py`. That script imports no pipeline code, and recomputes q_4(aP_u + bP_v) from the
  81 × 81 matrix for each ratio. It passes:
  - at 11 entangled points × 7 ratios (77 cases);
  - in a separate audit run, at 8 further points (6 of them near-zeros) × 11 ratios, including (1000, ±999)
    (88 cases; `logs/exact4L/audit_verifyL.log`).

**Other checks.**
- Negative controls (one multiplier changed; one dual diagonal entry negated) are rejected.
- Face data were recomputed with Python integers only: 5.2M, 7.0M and, for the cross data of Theorem 3, 15.2M
  entries.
- The engine agrees with the separately certified two- and three-copy code at random rational points: objective, all
  moments, all identity images, and the covariance blocks tested.
- A full rebuild reproduces the files byte for byte.

**To check:**

    cd code/exact4   && python verify4.py  p_ent p_ent2
    cd code/exact4m  && python verify4m.py p_ent p_ent2
    cd code/exact4L  && python verifyL.py exact4/p_ent exact4m/b_ent1

The first two take under a minute each, the third about a minute for these two points and 2–3 minutes for its default
11 points (Python 3, numpy, gmpy2; python-flint for `verify4.py`).

**What a referee has to trust.** Steps 1 and 2 are a complete proof that c·x ≥ 0 on the set {E x = e, F_k(x) ⪰ 0}.
What remains is that this set contains every true point, which is the validity of the constraint data.
- The saved data are produced by our code from the definitions in Section 1.
- They are supported by the cross-checks above and by step 3 at every test point.
- An independent re-derivation of E, e, c and the F_k from Section 1 would remove this last dependence. The saved
  pickles make a direct comparison possible.

## References

- [BGH] K. Bharti, R. Gajjala, T. Haug, arXiv:2607.24479.
- [Cl] L. Clarisse, The distillability problem revisited, Quantum Inf. Comput. 6, 539 (2006); arXiv:quant-ph/0510035.
- [FHPV] Fraser, Huber, Pozsgay, Vona, arXiv:2607.24309.
- [T1] T. S. Tomar, Undistillability of Werner states from three, four and five copies, Zenodo 23718035 (2026).
- [T2] T. S. Tomar, Equality cases and every-copy results for the Werner endpoint inequality, Zenodo 23718037 (2026).
- [WZ] T. Wu, Q. Zou, Sharp Plücker geometry for three-copy Werner distillation, arXiv:2608.02647.
