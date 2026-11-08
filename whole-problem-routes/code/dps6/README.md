# Six copies of 𝒫_{3/4}: the real second level, reduced by all symmetries

This folder proves Theorem 7.32 of [`whole-problem-routes/README.md`](../../README.md):

    ⟨y| 𝒫_{3/4}^{⊗6}(xx*) |y⟩ ≥ 4^{−7} |x|² |y|²   for all x, y ∈ C^{3^6},

with 𝒫_a(X) = a Tr(X) I − ¼(X + Xᵀ). Real vectors suffice, because 𝒫_a^{⊗n}(xx*) depends only on Re(xx*) and is real
symmetric. The proof is computer-assisted: the certificate is a list of rational numbers, and its validity is decided in
exact arithmetic by two separately written programs.

## What a certificate is

Each copy s has three strands: the input x_s and two copies y_s, y′_s of the output. Its fifteen Brauer diagrams D_0, …,
D_14 are the perfect matchings of the six points 0, 1, 2 (the rows of x_s, y_s, y′_s) and 3, 4, 5 (their columns). A
certificate gives one rational number z_μ for every multiset μ of n diagrams (38760 of them at n = 6). It defines

    Z = Σ_τ z_{sort(τ)} D_{τ₁} ⊗ ⋯ ⊗ D_{τ_n}        (τ runs over ordered n-tuples),

and claims

    (1)  ⟨x ⊗ y ⊗ y| Z |x ⊗ y ⊗ y⟩ = ⟨y|𝒫_a^{⊗n}(xxᵀ)|y⟩ |y|² − t |x|² |y|⁴   for all real x, y,
    (2)  Z ⪰ 0.

Together they give ⟨y|𝒫_a^{⊗n}(xxᵀ)|y⟩ ≥ t |x|²|y|². This is the second level of the Doherty–Parrilo–Spedalieri
hierarchy written for real vectors. The swap of the two x slots and the permutations of the four y slots do the work
that the partial transposes Γ_A, Γ_B and the symmetrizer P do in `../dps/`; numerically the two programs have the same
optimum.

**Why (1) is a finite check.** Both sides are quadratic forms in ψ = x ⊗ y ⊗ y with operators in the span of products of
diagrams. Two such operators give the same polynomial in (x, y) if and only if their averages over the group S_2 × S_4
agree, where S_2 permutes the two x slots and S_4 the four y slots, acting on every copy at once. For R³ the fifteen
diagrams of one copy are linearly independent, so the products of diagrams are a basis, and the averages agree if and
only if the coefficient sums agree on every orbit of tuples. The coefficients are symmetric under permutations of the
copies, so it is enough to compare, for each orbit of multisets, Σ_{μ ∈ orbit} (number of orderings of μ)·z_μ with the
same sum for C_a^{⊗n} ⊗ I − t I. At n = 6 there are 1042 orbits. The coefficients of C_a ⊗ I_{y′} in the diagram basis
are a for the identity and −¼ for the two diagrams that swap or cap x_s with y_s.

**Why (2) is a finite check.** Z commutes with O(3) acting on the three strands of each copy, and with permutations of
the copies. Per copy, the m = 0 vectors of the spin sectors j = 0, 1, 2, 3 have multiplicities 1, 3, 2, 1. Compressing to
these vectors (exact integer bases) splits Z into one block per spin type (k₀, k₁, k₂, k₃) (k_j copies of spin j), of
size 3^{k₁}2^{k₂}. Each type block commutes with S_{k₁} × S_{k₂}, and it is positive definite if and only if its
compressions to one copy of each irreducible V_{λ₁} ⊗ V_{λ₂} are (λ₁ ⊢ k₁ with at most three rows, λ₂ ⊢ k₂ with at most
two). At n = 6 these are 217 rational matrices of size at most 48.

## Files

| file | what it does |
|---|---|
| `sym_real_l2.py` | builds the reduced program: orbits, spin types, Gelfand–Tsetlin blocks (`gtbasis.py`); `save` writes the a-independent data |
| `sym_solve.py` | solves it (Clarabel): `maxt` maximises t, `center` fixes t and maximises the least block eigenvalue |
| `cert_l2.py` | `extract` recovers the diagram coefficients z from a solution, `round` makes them rational and exact (the residual of each orbit is put on one multiset of the orbit, then Z is replaced by (Z + Zᵀ)/2; `mix` combines the solution of largest t with the one of largest margin at t = 0), `export` writes the certificate, `check` verifies (1) and (2) in exact arithmetic (integer spin bases from `../tc3/exact_alg.py`, Young symmetrizer images, exact LDLᵀ) |
| `verify_l2.py` | a second, separately written exact verifier: it rebuilds the diagrams, the target and the orbits from the 27 × 27 matrices, recomputes the spin bases, reduces with Jucys–Murphy eigenvectors instead of Young symmetrizers, and tests (1) at random points of tensor rank two in floating point |
| `certs/` | the certificates (JSON, xz-compressed) |
| `run_checks.sh` | runs both verifiers on both certificates |
| `logs/` | the runs quoted below |

## Results

| file | n | a | t | orbits | blocks (largest) | checks |
|---|---|---|---|---|---|---|
| `certs/cert_n6_a3_4_t1_16384.json.xz` | 6 | 3/4 | 4^{−7} = 1/16384 | 1042 | 217 (48) | `logs/check_n6.log`, `logs/verify_n6.log`: PASS |
| `certs/cert_n4_a3_4_t1_512.json.xz` | 4 | 3/4 | 1/512 | 110 | 62 | test case (Theorem 7.13 has the sharp 1/256); PASS |

How the six-copy certificate was made (floating point, Clarabel; logs in `logs/`):

- the reduced program has 20 412 variables, 217 blocks and 1042 equations (`build_n6.log`, 17 minutes);
- maximising t gives t* = 1.374856·10⁻⁴ ≈ 0.563·4^{−6} (`solve_n6_maxt.log`). This is below the product value 4^{−6},
  so the second level is not exact at six copies (it is for n ≤ 5, Theorem 7.13);
- at t = 0 the least block eigenvalue can be pushed to 1.378·10⁻⁴ (`solve_n6_center_t0.log`);
- the certificate is θ·(first solution) + (1 − θ)·(second) with θ = 4^{−7}/t* = 0.4439. Its coefficients were rounded
  to 2^{−50}, which left orbit residuals of at most 1.2·10⁻⁸. Each residual was then put on one multiset of its orbit,
  and Z replaced by (Z + Zᵀ)/2 (`round_n6.log`). The orbit identity then holds exactly with t = 1/16384.

The smallest pivot in `check_n6.log` looks tiny relative to the diagonal (3·10⁻⁸). That is the scale of the integer
Young-symmetrizer basis, not a near-singularity. The pivots are exact rationals, and all of them are positive.

## Reproducing

    sh run_checks.sh

or one file at a time, for example

    python verify_l2.py certs/cert_n6_a3_4_t1_16384.json.xz
    python cert_l2.py 6 3 4 check certs/cert_n6_a3_4_t1_16384.json.xz

Both need numpy and python-flint. The second recomputes the integer spin bases with `../tc3/exact_alg.py` if
`exact_alg_k3.pkl` is not present. To rebuild the certificate from scratch:

    python sym_real_l2.py 6 3 4 CLARABEL save          # about 17 minutes, 1.2 GB
    python sym_solve.py 6 3 4 maxt CLARABEL             # about 13 minutes, 10 GB committed
    python sym_solve.py 6 3 4 center CLARABEL 0
    python cert_l2.py 6 3 4 extract maxt
    python cert_l2.py 6 3 4 extract center
    python cert_l2.py 6 3 4 round 50 1/16384 mix
    python cert_l2.py 6 3 4 export certs/cert_n6_a3_4_t1_16384.json.xz

The two solves need about 10 GB of memory each (most of it the factorisation of the dense equality rows).
