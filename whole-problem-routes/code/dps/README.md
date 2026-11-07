# Exact certificates for three, four and five copies of 𝒫_{3/4}

This folder proves Theorem 7.13 of [`whole-problem-routes/README.md`](../../README.md):

    ⟨x ⊗ y| C_{3/4}^{⊗n} |x ⊗ y⟩ ≥ 4^{−n} |x|² |y|²   for all x, y ∈ C^{3^n},   n = 3, 4, 5,

with equality at x = y = e₁^{⊗n}, and the same inequality with C_{81/125} and the constant 0 for n = 3. Here
C_a = aI − ¼(F + |Ω⟩⟨Ω|) is the Choi matrix of 𝒫_a(X) = a Tr(X) I − ¼(X + Xᵀ), so the left side is
⟨y|𝒫_a^{⊗n}(|x̄⟩⟨x̄|)|y⟩. The proof is computer-assisted: every certificate is a list of rational numbers, and its
validity is decided in exact arithmetic.

## What a certificate is

Let B = B₁…B_n be the outputs, B′ a second copy of B, S the swap of B and B′, P = (I + S)/2, and Γ_A, Γ_B the partial
transposes on all inputs and on all outputs. A certificate is a pair of rational operators Q_A, Q_B on A ⊗ B ⊗ B′ with

    Q_A′ = P Q_A P ⪰ 0,    Q_B ⪰ 0,    Z = P [ C_a^{⊗n} ⊗ I_{B′} − t I − Γ_A(Q_A′) − Γ_B(Q_B) ] P ⪰ 0.

For ψ = x ⊗ y ⊗ y one has Pψ = ψ, and

    ⟨xy|C_a^{⊗n}|xy⟩ |y|² − t |x|² |y|⁴ = ⟨ψ|Z|ψ⟩ + ⟨x̄yy|Q_A′|x̄yy⟩ + ⟨xȳy|Q_B|xȳy⟩ ≥ 0,

because Γ_A and Γ_B map the product |x⟩⟨x| ⊗ |y⟩⟨y| ⊗ |y⟩⟨y| to the same product with x, respectively the first y,
conjugated. This is the second level of the Doherty–Parrilo–Spedalieri hierarchy, with the extension on the output
side.

The operators are rational combinations of products D_{a₁} ⊗ … ⊗ D_{a_n} of the fifteen Brauer diagrams of one copy
(operators on the three strands A_s, B_s, B′_s of copy s). File format (`certs/`):

- lines `n`, `a` (numerator and denominator), `t`, and `L` (a common denominator);
- lines `D0` … `D14`: the diagrams, as perfect matchings of the six points 0, 1, 2 (the rows of the strands A_s, B_s,
  B′_s) and 3, 4, 5 (their columns);
- lines `QA` and `QB`: the integer numerators of the coefficients of D_{a₁} ⊗ … ⊗ D_{a_n} (15^n of them, in C order, so
  the first copy is the most significant index). The exact five-copy certificate uses the compressed lines `QAO` and
  `QBO`: one numerator for each orbit of index tuples under permutations of the copies and the transpose, the orbits
  ordered by the key min(sorted(a), sorted(τa)), where τ is the transpose acting on diagrams.

| file | n | a | t | remark |
|---|---|---|---|---|
| `cert_n3_a3_4_t1_64_exact.txt` | 3 | 3/4 | 1/64 | exact minimum; Q_B = 0 |
| `cert_n3_a3_4_t1_64_exactB.txt` | 3 | 3/4 | 1/64 | exact minimum; Q_A = 0 |
| `cert_n3_a3_4_t1953_125000.txt` | 3 | 3/4 | 1/64 − 10⁻⁶ | strictly feasible |
| `cert_n3_a3_4_t1_128.txt` | 3 | 3/4 | 1/128 | not symmetric under permutations of the copies |
| `cert_n3_a81_125_t0.txt` | 3 | 81/125 | 0 | |
| `cert_n4_a3_4_t1_256_exact.txt.xz` | 4 | 3/4 | 1/256 | exact minimum; Q_A = 0 |
| `cert_n4_a3_4_t1_512.txt.xz` | 4 | 3/4 | 1/512 | |
| `cert_n5_a3_4_t1_1024_exact.txt.xz` | 5 | 3/4 | 1/1024 | exact minimum; compressed lines |
| `cert_n5_a3_4_t1_2048.txt.xz` | 5 | 3/4 | 1/2048 | |

## How positivity is decided

Every product of diagrams commutes with SO(3) acting on the three strands of each copy. So an operator in their span
is positive semidefinite if and only if each of its compressions to the m = 0 vectors of a spin sector (j₁, …, j_n),
j ∈ {0, 1, 2, 3}, is positive semidefinite (Schur's lemma; per copy the multiplicities are 1, 3, 2, 1). These are
exact rational matrices of size at most 3^n.

At t = 4^{−n} every real product vector x = y = v₁ ⊗ … ⊗ v_n forces a zero of Q_A′, Q_B and Z. Per copy these vectors
span Sym³(C³) inside A_s ⊗ B_s ⊗ B′_s, which meets the sectors j = 1 and j = 3 in one m = 0 vector each. So each sector
with all j_s ∈ {1, 3} has one forced kernel vector, and the exact certificates have no other kernel.

At one copy there is a closed form. With Y = ¼(I − F_{AB}) ⊗ I_{B′} one has C_{3/4} ⊗ I − ¼I = Y + Γ_A(Y), and
PYP = ⅜(P − Π), where Π projects onto Sym³(C³). So Q_A′ = Z = ⅜(P − Π) and Q_B = 0 is a certificate with t = ¼.

## Verification

There are two sets of exact verifiers, written separately.

- `verify_independent.py` shares no code with the generating scripts. It rebuilds the diagram matrices from the
  matchings, the action of Γ_A, Γ_B, of S on either side and of the transpose (from explicit 27 × 27 matrices), the
  operator C_a ⊗ I from the definition of 𝒫_a, and the m = 0 spin bases (exact null spaces of the so(3) generators and
  Casimir). It forms Q_A′, Q_B and Z as exact integer tensors, checks that they are symmetric, and tests every sector
  in two ways. The first is fraction-free symmetric elimination with diagonal pivoting (used here on sectors of size
  at most 27). The second is a congruence test: the forced vector is checked to be an exact null vector and one
  coordinate on its support is deleted; for the rest an integer matrix M is built from floating-point eigenvectors,
  H = MᵀGM is computed exactly, and strict diagonal dominance of H with positive diagonal proves G ≻ 0. Logs
  `logs/independent_*.log`.
- `crosscheck3.py` checks the conventions at three copies in floating point. For random complex x, y it computes the
  left side from the definition of 𝒫_a (applied copy by copy, and through the 729 × 729 Choi matrix) and the right side
  from the certificate. They agree to 2·10⁻¹⁶ relative, every term is nonnegative, and both sides vanish at product
  minimisers (`logs/crosscheck3_*.log`).
- `exact_verify.py` (exact LDLᵀ over the rationals), `exact_verify_fast.py` (modular arithmetic, Sylvester's criterion,
  exact kernels) and `exact_verify_mod.py` (modular and memory-light, with an exact congruence and a Gershgorin test for
  the large sectors) come with the generating scripts. Logs `logs/verify_*`, `logs/fastverify_*`, `logs/modverify_*`.
  They read plain-text certificates, so run `python unpack_certs.py` first.

`bash run_checks.sh` reruns the independent verifier on all nine certificates, the three-copy crosscheck, two of the
other verifiers, and three negative controls (about ten minutes). Each control is an exact certificate with t raised just
above the product value, and each one must fail. Every check passes and every control fails
(`logs/run_checks.out`).

## How the certificates were found

Floating-point semidefinite programs, followed by exact rounding:

- `brauer.py` (the per-copy algebra: diagrams, partial transposes, products, spin bases), `blocks.py`, `level2_sdp.py`
  (the program in block coordinates, CVXPY with SCS), `make_cert.py`;
- `sym_sdp.py`, `sym_sdp5.py` (the program reduced by permutations of the copies and the transpose; Clarabel, and SCS
  for n = 5), `make_cert_sym.py`, `make_cert_orbit.py` (rounding at 2^{−40} to 2^{−48}, then exact symmetrisation);
- `sym_sdp_cones.py`, `sym_sdp5_cones.py` (one partial transpose only), `sym_sdp5_face.py` (the face with the forced
  zeros imposed);
- `exact_face.py`, `exact_face2.py`, `exact_face3.py`, `exact_face_ab.py`: the exact certificates at t = 4^{−n}, by
  the least correction (an exact linear system) that makes the forced vectors exact null vectors;
- `float_crosscheck*.py`: floating-point checks of the identity at n = 3, 4, 5. `analyze_face.py` looks at the
  optimal face.

Solver logs are in `logs/sdp/`. Reproduction (each program takes under two minutes at n ≤ 4; use
`bash guard25.sh LOG SECONDS script args` for the heavy ones):

    python level2_sdp.py 3 3 4 center 0.0078125 SCS
    python make_cert.py sol_n3_3_4_center.npz 3 3 4 1 128 40 certs/cert_n3_a3_4_t1_128.txt
    python sym_sdp.py 3 3 4 center 0.015624 1e-11 CLARABEL
    python make_cert_sym.py sol_sym_n3_3_4_center.npz 3 3 4 1953 125000 48 certs/cert_n3_a3_4_t1953_125000.txt
    python sym_sdp.py 3 81 125 center 0 1e-11 CLARABEL
    python make_cert_sym.py sol_sym_n3_81_125_center.npz 3 81 125 0 1 44 certs/cert_n3_a81_125_t0.txt
    python sym_sdp.py 4 3 4 center 0.001953125 1e-10 CLARABEL
    python make_cert_sym.py sol_sym_n4_3_4_center.npz 4 3 4 1 512 40 certs/cert_n4_a3_4_t1_512.txt
    python sym_sdp5.py 5 3 4 center 0.00048828125 1e-7 SCS          (200 s, 2.3 GB)
    python make_cert_orbit.py sol_sym_n5_3_4_center.npz 5 3 4 1 2048 40 certs/cert_n5_a3_4_t1_2048.txt
    python sym_sdp_cones.py 3 3 4 maxt 0 1e-10 CLARABEL A
    python exact_face2.py sol_sym_n3_3_4_maxt.npz 3 1 64 40 certs/cert_n3_a3_4_t1_64_exact.txt
    python sym_sdp5_face.py 3 3 4 center 0.015625 1e-10 CLARABEL B
    python exact_face3.py sol_sym_n3_3_4_center_B.npz 3 1 64 40 certs/cert_n3_a3_4_t1_64_exactB.txt
    python sym_sdp5_cones.py 4 3 4 maxt 0 1e-10 CLARABEL B
    python exact_face3.py sol_sym_n4_3_4_maxt_B.npz 4 1 256 44 certs/cert_n4_a3_4_t1_256_exact.txt
    python sym_sdp5.py 5 3 4 maxt 0 1e-6 SCS
    python exact_face_ab.py sol_sym_n5_3_4_maxt.npz 5 1 1024 40 certs/cert_n5_a3_4_t1_1024_exact.txt

The rounding steps depend on the floating-point solutions, so a rerun gives different but equally valid certificates.

## Numerical observations (floating point, not part of the proof)

The best t of the second level, with both partial transposes:

| n | a | best t | product value (a − ½)^n | value of the witness of Theorem 7.1 |
|---|---|---|---|---|
| 1, 2, 3, 4 | 3/4 | 4^{−n} (to 10⁻¹⁰) | 4^{−n} | |
| 5 | 3/4 | 0.00097657 (SCS, tolerance 10⁻⁶) | 0.0009765625 | |
| 3 | 0.70 | 0.0080000 | 0.008 | 0.0286 |
| 3 | 0.65 | 8.71·10⁻⁴ | 3.375·10⁻³ | 1.500·10⁻³ |
| 3 | 0.648 | 6.6·10⁻⁶ | | 5.45·10⁻⁴ |
| 3 | 0.647 | −4.26·10⁻⁴ | | 7.08·10⁻⁵ |
| 3 | 0.60 | −0.019625 | 0.001 | −0.019625 |

So at a = 3/4 the second level is exact for n ≤ 5, and at a = 0.6 it finds the witness value. The first level
(decomposability, by the linear program of Proposition 7.3 with the constant t added) is exact at n = 1, 2 but reaches
only t ≈ −0.0036 at n = 3, so three copies is where the extension is first needed. Near the threshold the second level
is not exact: at n = 3 it certifies positivity only from a ≈ 0.64798 on, above a₃ ≈ 0.646850. With one partial
transpose only: at n = 2, 3 either one reaches 4^{−n}; at n = 4 the transpose on the outputs reaches 1/256, the one on
the inputs only 0.0029257. With no partial transpose (the symmetric extension alone) the best t at n = 3 is −0.175.
Five copies is the largest case this implementation handles. At n = 6 there are 15⁶ ≈ 1.1·10⁷ coordinates and sectors
of size up to 729, which would need a reduction under the permutations of the copies.
