# Four-replica relaxations at the qutrit Werner endpoint

Ansh Mishra and Aryan Senthilkumar

**This folder does not prove the endpoint for every n, or for general rank-two C at four copies.** It proves one
four-copy case, with a computer-assisted exact certificate: q_4(P) ≥ 0 for every rank-two projector P (Section 3,
"Exact certificate at four copies"). Its results are:
- an exact rational certificate at four copies for rank-two projectors (not externally refereed);
- semidefinite relaxations of the endpoint inequality at four copies, with what they certify numerically;
- one exact identity about one of their constraints;
- a first-order reformulation on a slice;
- explicit failures of the natural four-copy versions of the three-copy proofs.

The main numerical findings:
- **Rank-two projectors (four copies).** A four-replica relaxation is exact at n ≤ 4. After facial reduction, the
  letter-resolved version of Section 3 has a strictly feasible reduced dual at n = 4, at both zero directions:
  - projectors: margin 2.7e-3 with the covariance block only, residual 1e-7;
  - opposite signs: margin 1.9e-3.

  For projectors we rounded this dual to an exact rational certificate, so q_4(P) ≥ 0 holds for every rank-two
  projector P (computer-assisted). The opposite-sign case is still numerical.
- **Hermitian rank-two operators.** A letter-resolved relaxation matches the conjectured sharp margin
  2^{−n}(σ₁ − σ₂)² at n = 4 on a grid of 13 ratios. A single lifted certificate covers all ratios at once (SCS, value
  ≈ 0 to 1e-7).

Throughout, d = 3, ρ = I − F/2 and q_n(C) = Σ_J (−1/2)^{|J|} ‖Tr_J C‖² for C of rank at most two on (C³)^⊗n. The
endpoint inequality at n copies (EP_n) is q_n(C) ≥ 0 for all such C.

## 1. Rank-two projectors, one letter

Let P be a rank-two projector. Twirl its fourth tensor power over local unitaries,
Z = E (U^⊗4) P^⊗4 (U^⊗4)*, with U = U_1 ⊗ … ⊗ U_n. Then Z lies in A^⊗n, where A is spanned by the 24 replica
permutations acting on (C³)^⊗4. Moreover Tr[Z (G ⊗ I ⊗ I)] = 4 q_n(P) with G = ⊗_i (F_i − I/2).

Constraints, each valid for every rank-two projector (`code/relaxations/rep4.py`):
- Z ≥ 0;
- the moments Tr[Z τ] = 2^{cycles(τ)};
- Z ≤ Tr_4 Z ⊗ I/2;
- one- and two-replica partial transposes;
- the identity from P² = P;
- the **covariance block** twirl(P ⊗ P ⊗ |vec P⟩⟨vec P|) ≥ 0. This is a realignment of Z, block-diagonalised in
  the walled Brauer algebra B_{3,1}(3) on each copy.

Results (Clarabel; residuals 1e-8 or better):

| n | Without the covariance block | Z ≥ 0, moments, covariance only | All constraints |
|---|---|---|---|
| 1 | | 2e-15 | |
| 2 | | 8e-10 | |
| 3 | −0.387 | −1.7e-8 | −3e-9 |
| 4 | −0.305 | −5.1e-8 | −7.1e-8 |
| 5 | | −7.47e-3 | −7.47e-3 (covariance and P² = P only) |

- **Controls.** Rank-three projectors reproduce the true negative minima: −1.5, −0.75, −1.375 (I ⊗ |Φ⟩⟨Φ|) and
  −0.9375 (I ⊗ |GHZ⟩⟨GHZ|).
- **Below the endpoint** (weight w < 1/2 in place of 1/2) the value equals the product value
  (2 − 4w)(1 − w)^{n−1} at n = 3, 4.
- So the covariance block makes this relaxation exact for projectors at n ≤ 4, and it is not exact at n = 5 in the
  variants that fit in memory.

## 2. The covariance block is already in Tomar's relaxation

Tomar's four-replica relaxation [T1] for general C lifts C = Tr_Q |a⟩⟨b| with a, b ∈ (C³)^⊗n ⊗ C². Then vec C is
(a ⊗ b̄) contracted with the maximally entangled vector on Q ⊗ Q. So the covariance blocks are compressions of
whole-replica partial transposes of his moment blocks:
- letter pair (a, a) comes from X31, transposed on replica 4;
- letter pair (a, b) comes from X22, transposed on replica 3;
- letter pair (b, b) comes from X13, transposed on replica 1.

For vector letters, realigning a replica pair is a partial transpose. A compression of a PSD block is PSD. So his
relaxation already contains every covariance constraint, and his value −1.31e-3 at four copies stands
(`lift_identity.py` checks the identity to 1.4e-17). The one-letter relaxation of Section 1 needs the block explicitly
only because it works with P directly, without the purifying qubit.

## 3. Hermitian rank-two operators, two letters

Write C = μ P_u + s ν P_v with u ⊥ v unit vectors, μ, ν ≥ 0 and s = ±1. The variables are the twirls X_l of
⊗_r P_{l_r} for letter strings l ∈ {u, v}^4; the canonical families are uuuu, uuuv, uuvv, uvvv, vvvv
(`herm2.py`, and `herm4.py` for the symmetry-reduced version). The objective is
μ² m_uu + ν² m_vv + 2 s μ ν m_uv, with m_xy = Tr[X (G ⊗ I ⊗ I)] on the matching family.

Constraints, each valid for every orthonormal pair:
- the moments;
- **pure letters:** X_l F_ij = X_l whenever l_i = l_j. In `herm4.py` this is built into the parametrisation: each
  Fourier block lives on the subspace fixed by that letter pattern's permutations;
- **orthogonality:** Tr_j[X_l F_ij] = 0 whenever l_i ≠ l_j (because P_u P_v = 0);
- shared three-replica marginals;
- one- and two-replica partial transposes;
- the coupled covariance block (the Gram matrix of vec P_u and vec P_v next to two further replicas);
- P_u + P_v ≤ I.

**Validation.** Random orthonormal pairs (complex and real) satisfy every constraint at n = 2, 3 and 4:
- reconstruction error ≤ 2e-16;
- smallest block eigenvalue ≥ −8e-17;
- moment and identity residuals ≤ 4e-15;
- the objective equals the directly computed q_n(C) to 12 digits.

We re-ran the n = 2 and n = 3 checks separately (`logs/audit_validate_herm4_n2.log`, `..._n3.log`). Exact zero
configurations also reconstruct at n = 4 to 1e-17.

**n = 2 and n = 3.** The relaxation is exact on every case tried, including both zeros. At n = 3 it returns
(|λ₁| − |λ₂|)²/8 for C = (0.8, −0.2), (0.6, 0.4) and (0.75, 0.25): 0.0450000, 0.0050001, 0.0312500. At n = 3 it stays
exact without the covariance block or the ≤ I bound; the identities do most of the work.

**n = 4.** The symmetry-reduced problem has 1331 parameters, 1902 PSD blocks (largest 72) and 80 210 rows. Building
takes about 260 s and each solve about 2 min, with a peak of 8.4 GB. The grid is C = P_u + t P_v
(`logs/herm4_n4_grid.log`, `logs/check_n4.log`). The corrected bound is −b·z − ‖Aᵀz + c‖_∞ ‖x‖₁, a lower bound for
every feasible point given exact problem data.

| t | Primal | Dual −b·z | Residual | Corrected bound | True minimum |
|---|---|---|---|---|---|
| −1 | −7.52e-8 | −7.74e-8 | 4.8e-6 | −1.1e-5 | 0 |
| −0.5 | 1.562563e-2 | 1.562562e-2 | 1.7e-5 | 1.56e-2 | 1/64 |
| −0.25 | 3.515750e-2 | 3.515748e-2 | 2.9e-5 | 3.51e-2 | 9/256 |
| 0 | 6.250001e-2 | 6.250001e-2 | 1.3e-6 | 6.25e-2 | 1/16 |
| 0.25 | 3.515890e-2 | 3.515886e-2 | 6.1e-5 | 3.50e-2 | 9/256 |
| 0.5 | 1.562519e-2 | 1.562519e-2 | 5.4e-6 | 1.56e-2 | 1/64 |
| 1 | −4.24e-7 | −4.36e-7 | 2.2e-5 | −4.0e-5 | 0 |

Six more ratios (primal values; `logs/herm4_n4_dense_grid.log`) also match (1 − |t|)²/16 to within 1.3e-6:
t = −0.9, −0.75, −0.1, 0.1, 0.75, 0.9 give 6.25108e-4, 3.90637e-3, 5.06256e-2, 5.06253e-2, 3.90750e-3 and
6.25126e-4. That makes 13 ratios in [−1, 1].

- **The minimum is (1 − |t|)²/16.**
  - For t > 0 it is attained by the products diag(1, t, 0) ⊗ |x⟩⟨x|^⊗3.
  - For t < 0 it is attained by e₀e₀* ⊗ xx* − |t| e₁e₁* ⊗ yy*, where x, y are product states orthogonal on an odd
    number of the remaining copies.
  - Local searches (12 restarts per point, `local_herm.py`) find the same values.
- **Margin conjecture.** (1 − |t|)²/16 is the sharp margin q_n(C) ≥ 2^{−n}(σ₁ − σ₂)² conjectured in
  `first-order-bound/`, with σ the singular values (here |eigenvalues|). That margin is a theorem at n = 2 [FHPV].
  For Hermitian C it is a theorem at n = 3: [WZ, Thm 1.3] for C ≥ 0 and [CR, Prop. 3] for opposite signs. At
  n = 4 the relaxation's lower bounds match it on the grid, to solver accuracy.
- **Control above the endpoint.** At w = 0.55 the relaxation returns −0.117 (t = 1) and −0.200 (t = −1), below the
  values −0.0182 and −0.0405 of explicit operators (`logs/audit_above_endpoint.log`). So it is not trivially
  nonnegative.
- **Zeros.** Suppose u = a ⊗ φ and v = a′ ⊗ χ, with a ⊥ a′ on one copy. Splitting the sum over J by whether that
  copy is traced gives

      q_n(P_u + t P_v) = ½ q_{n−1}(P_φ − t P_χ).

  This is the case X_01 = X_10 = 0 of the one-copy contraction formula in `contraction-notes/`. Iterating it
  gives every zero family used below, including interpolating ones and U(2) rotations of the pair (`zero_rot.py`,
  `zero_interp.py`). All of them have rational points.
- **Facial reduction at the zeros** (`fr_herm4.py`, `fr2_herm4.py`). Every dual certificate of value 0 must vanish on
  the ranges of all PSD blocks at the zeros. So we compress every block to its kernel at known zeros and re-solve.
  - t = 1 (rank-two projectors): 3 structured zeros (P_1 ⊗ ψ) and 2 found by local search. Of 2169 blocks, 1280
    are kept (dimension 10749 → 5116). The reduced problem is **Solved**: primal 1.1e-13, dual −9.4e-15, residual
    9.3e-10, **smallest eigenvalue of the reduced dual blocks 3.2e-4**. This dual is strictly feasible with a margin far
    above the residual, so an exact rational certificate should exist. It would prove EP_4 for every rank-two
    projector.
  - t = −1 (opposite signs): with the structured zeros alone the reduced dual is still degenerate (2.4e-11). Adding
    three zeros found by local search (12 zeros in all, including new interpolating families) gives **Solved**:
    primal 3.7e-12, dual −9.8e-15, residual 1.0e-7, **smallest reduced dual eigenvalue 1.9e-3**.
  - t = 1 again, with only the identities and the covariance block (no partial transposes, no ≤ I bound): Solved,
    residual 1.3e-7, smallest reduced dual eigenvalue 2.7e-3. So an exact projector certificate needs only the S_4
    blocks, the identities and the covariance blocks.
  - At n = 3 both signs give strictly feasible reduced duals (2.0e-3 and 5.6e-3).
- **One certificate for all ratios** (`lift2.py`). Lift every cone row to a 2 × 2 block, with variables
  X = (X11, X12, X22) and E ≥ 0, Tr E = 1, and minimise m_uu·X11 + 2m_uv·X12 + m_vv·X22. Since (x ⊗ vvᵀ, vvᵀ) is
  feasible, a value ≥ 0 certifies [[m_uu·x, m_uv·x], [m_uv·x, m_vv·x]] ≥ 0 on the whole relaxation, that is, every
  ratio and sign at once.
  - The lifted value is 4.4e-8 at n = 2 and 4.8e-7 at n = 3 (`logs/lift2_n2.log`, `logs/lift2_n3.log`).
  - At n = 4 the lifted problem has 1.4e8 nonzeros. The interior-point solver ran out of memory on its dense
    scaling blocks (total size 8e8). The first-order solver SCS, on a leaner assembly (`lift2s.py`), solved it:
    2325 iterations, primal −6.6e-8, dual +1.8e-9, residual 2.5e-6 (`logs/lift2s_n4.log`). The same assembly
    reproduces the Clarabel values at n = 2 and n = 3.
  - So at four copies one numerical certificate covers every Hermitian rank-two operator, for all ratios and both
    signs, to solver accuracy.

### Exact certificates at two and three copies (pipeline validation)

For rank-two projectors (t = 1) we turned the two-letter relaxation into exact rational certificates at n = 2 and
n = 3 (`exact_pilot_n2.py`, `exact_pilot_n3.py`, with `exact_site.py`, `exact_maps.py`, `exact_core.py`,
`exact_coreN.py`; python-flint and gmpy2). Both statements are known; the point is that the pipeline works end to end.
The steps:
1. **Exact variable space and data.** Per-site irreducible modules come from integer Casimir eigenspaces intersected
   with U(3) highest-weight coordinate subspaces, with rational Gram matrices. The swap identities are imposed by an
   exact basis of the fixed subspace, and site symmetry by exact linear conditions. The maps for the moments, the
   identities, the covariance blocks and the objective are exact.
2. **Exact zeros.** We use u = a ⊗ ψ, v = a′ ⊗ ψ, with a and a′ rows of rational orthogonal matrices (Cayley
   transforms). At n = 3 we add orthogonal products and rationally rotated interpolating zeros. Every zero is checked
   to have objective exactly 0, and the face is the exact common kernel.
3. **Float, then exact.** A floating-point solve of the face-reduced problem gives the dual, which we round to a 2^{−40}
   grid. An exact least-norm correction then makes c = Eᵀy + Σ_k M_k*(N_k S_k N_k^T) and e·y = 0 hold exactly.
4. **Exact positivity.** An exact LDLᵀ shows every reduced dual block S_k is PSD.
5. **Independent check.** On an exact entangled projector, c·x(P) = q_n(P) = y·E x(P) + Σ_k Tr(S_k M_k(x(P))), with
   y·E x(P) = 0 and every term ≥ 0.

| | n = 2 | n = 3 |
|---|---|---|
| parameters / identity rank | 55 / 32 | 255 / 64 |
| reduced blocks | 51 | 70 |
| smallest reduced-dual eigenvalue (float) | 0.093 | 0.242 |
| exact correction size | 4.6e-8 | 5.0e-8 |
| smallest exact LDLᵀ pivot | 0.121 | 1.22e-2 |
| largest rational (numerator / denominator bits) | 1198 / 1193 | 315 / 317 |
| run time | 4.5 s | 117 s |
| log | `logs/exact_pilot_n2.log` | `logs/exact_pilot_n3.log` |

Both runs end in CERTIFIED, and we re-ran both independently.

### Exact certificate at four copies: rank-two projectors

**Theorem (computer-assisted).** q_4(P) ≥ 0 for every rank-two orthogonal projector P on (C³)^⊗4.

To our knowledge this case was open. At three copies it follows from [WZ]. Numerical evidence for four copies of
Werner states goes back to [Cl]. The proof is an exact rational dual
certificate for the two-letter relaxation above, restricted to the face of its known zeros. It has been checked
independently but not externally refereed. Code: `code/exact4/` (build log `code/exact4/BUILD_LOG.md`); logs:
`logs/exact4/`.

**What is exact.**
- *Relaxation data at n = 4.* There are 1331 parameters. The 120 moments and the identities give 147 independent
  equalities E x = e; the identities have exact rank 146, equal to the floating-point rank. The objective is c. There
  are 256 face-reduced PSD blocks F_k: 98 replica blocks and 158 covariance blocks, the largest of size 33.
- *The face.* 11 exact rational zeros (the k1 and k3 product families, their rational rotations, and rotated
  interpolating zeros). Each has objective exactly 0, and q_4(P) = 0 computed directly from P. The face is the exact
  common kernel of the blocks at these zeros. Nine of the zeros already give the same face.
- *The certificate.* A vector y and blocks S_k with c = Eᵀy + Σ_k F_k*(S_k) exactly in all 1331 coordinates and
  e·y = 0. Every S_k is PSD by exact LDLᵀ (smallest pivot 9.4e-5). Entries have up to about 1060 bits.

So every feasible x has c·x = y·e + Σ_k Tr(S_k F_k(x)) ≥ 0.

Every real rank-two projector gives a feasible point x(P), averaged over permutations of the copies, with
c·x(P) = q_4(P). For complex P, the real part of x(P) is the average of the points for P and P̄. Both are feasible
with value q_4(P), because every constraint map has rational coefficients. Hence q_4(P) ≥ 0.

**How it was computed.** The pipeline is the one used at n = 3, with these changes:
- the data are computed modulo 27-bit primes and recovered by the Chinese remainder theorem with a-priori bounds;
- the covariance blocks are compressed by exact rational symmetry-adapted bases;
- blocks are cut to the span of their data by LLL-reduced integer bases;
- the float dual is rounded to a 2^−40 grid, then corrected by an exact 1291 × 1291 solve (largest change 6e-5).

The float face-reduced dual (Clarabel, 16 iterations) has smallest normalised eigenvalue 0.248. After the 13 exact
points (about 8 minutes in parallel), the build takes about 3 minutes, with a peak of 3.9 GB.

**Independent checks.**
- `verify4.py` reads only the saved data and certificate. It re-checks the identity and the LDLᵀ. At exact entangled
  true points it checks:
  - c·x(P) = q_4(P), with q_4(P) recomputed from the 81 × 81 matrix P;
  - E x(P) = e;
  - the identity at the point, with every term ≥ 0;
  - every face block is PSD there.

  It passes at eight points. Two were made during the build. Six were made afterwards, from fresh seeds, by a separate
  script (`points4_audit.py`): four generic points, and two near-zeros with q_4 = 2.4e-7 and 6.4e-5. Two further
  near-zeros exceeded that script's fixed prime list and were not completed (`logs/exact4/audit_pt_a_near2.log`,
  `..._near3.log`).
- `crosscheck4.py` recomputes 25 covariance face blocks (5.2M entries) with Python integers only, without primes.
- `test_engine.py` compares the n = 4 engine with the separately certified n = 2, 3 code (`exact_coreN.py`) at random
  rational points. It covers the objective, all 120 moments, all 210 identity images and 10 covariance blocks, with no
  mismatches.
- Negative controls: `verify4.py` rejects a certificate with one multiplier changed, and one with a dual diagonal entry
  negated.

**A pitfall.** The modular moment routine `exact_coreN.moments_modp` sums up to 243 products in one einsum at n = 4.
With 28-bit primes this can overflow int64 silently. `points4.py` uses 27-bit primes and asserts the bound on every
call. The n = 2 and n = 3 runs are not affected.

**To check it** (from `code/exact4/`; needs numpy, gmpy2 and python-flint):

    python verify4.py p_ent p_ent2 a_ent1 a_ent2 a_ent3 a_ent4 a_near1 a_near4

This takes about a minute and ends in VERIFIED. The saved files are Python pickles. The full rebuild in `BUILD_LOG.md`
regenerates all of them from scratch (it also needs scipy and Clarabel), in about 15 minutes on a laptop.

## 4. Positive operators with a fixed spectrum, one letter

These results use C = λ P_u + (1 − λ) P_v with the one-letter relaxation of Section 1, adapted to the spectrum (cubic
identity ω³ − ω² + λ(1 − λ) ω = 0; `rep4.py --spec`, `validate_spec.py`, `logs/spec_n4.log`).

| λ | n = 3 relaxation | n = 3 minimum | n = 4 relaxation | n = 4 product value |
|---|---|---|---|---|
| 0.5 | −2.8e-9 | 0 | −7.6e-9 | 0 |
| 0.6 | 4.83e-3 | 5.00e-3 | 2.29e-3 | 2.5e-3 |
| 0.75 | 3.09e-2 | 3.125e-2 | 1.53e-2 | 1.5625e-2 |
| 0.9 | 7.99e-2 | 8.00e-2 | 3.99e-2 | 4.0e-2 |

Away from λ = 1/2 this relaxation is positive but not exact. The two-letter relaxation of Section 3 is exact there.

## 5. Exclusive clicks and the pair-traceless slice

For the Schmidt-rank-two vector ψ = vec C, let K be the set of copies on which the measurement {Φ_i, I − Φ_i} gives
Φ_i (Φ_i is the maximally entangled projector of copy i). Then

    q_n(C) = ‖C‖² E[(−1/2)^{|K|}],    P(K) = W_{[n]∖K}(C)/‖C‖²,

where W_S(C) is the squared norm of the component of C that is traceless on S and proportional to the identity off S.

**The slice.** If Tr_{ij} C = 0 for all i < j, then every sector with two or more identity sites vanishes, so
q_n(C) = W_[n](C) − ½ Σ_i W_{[n]∖i}(C). On this slice EP_n is the single inequality Σ_i ‖Tr_i C‖² ≤ 2‖C‖²: a
Schmidt-rank-two state whose clicks are exclusive clicks with probability at most 2/3. Tomar's barrier profile
a = (1, 4/9, 0, …) lives here. If EP holds, that profile cannot come from a rank-two operator for n ≥ 5.

**FO1.** For every C of rank at most two, Σ_i W_{[n]∖i}(C) ≤ (2/3)‖C‖², i.e. P(exactly one click) ≤ 2/3. Equivalent
forms:
- the slope of t ↦ ⟨vec C|(I − 3tΦ)^⊗n|vec C⟩ at the PPT point t = 1/3 is at least −2‖C‖²;
- every Z = Σ_i I_i ⊗ X_i, with X_i fully traceless, has s₁(Z)² + s₂(Z)² ≤ (2/3)‖Z‖_F²;
- ⟨ψ|M_1|ψ⟩ ≤ 2/3 for Schmidt rank two, where M_1 = Σ_{|S|=1} M_S and M_S = ⊗_{i∈S} Φ_i ⊗_{i∉S} (I − Φ_i).

Facts:
- **EP_n ⟹ FO1_n at d = 3, exactly.** We have (2/3)I − M_1 = (2/3)(ρ^Γ)^⊗n + D with
  D = Σ_{|S|≥2} (2/3)(1 − (−1/2)^{|S|}) M_S ≥ 0. Checked to 6e-17 (`code/clicks/fo1_facts.py`).
- **FO1 ⟹ EP on the slice.**
- **Numerics.** The maximum is exactly 2/3 for n = 2, …, 6 at d = 3 (`logs/fo1_small.log`, `fo1_n5.log`,
  `fo1_n6.log`). It equals 1/2 = 2/d at (n, d) = (2, 4). It fails at d = 5 (0.44 > 0.40, attained by the product
  P_2 ⊗ |u⟩⟨u|), so FO1 is specific to d ≤ 4.
- **Tight families** (`logs/fo1_facts.log`).
  - P_2 ⊗ |0…0⟩⟨1…1|.
  - The W-type operator (|W⟩⟨1…1| + |0…0⟩⟨W′|)/√n. It lies on the slice and is a zero of q_n, but it is
    qubit-supported, so it belongs to a known zero family.
- **What fails.**
  - Rank-free and trace-norm versions fail already for a rank-one product at n = 2.
  - At n = 2, (2/3)I − M_1 is not a nonnegative combination of single-copy witnesses tensored with I − Φ, of
    (2/9)I − Φ_1Φ_2, and of the M_S. So two copies already need genuinely two-copy input.

FO1 is a conjecture for n ≥ 4; for n ≤ 3 it follows from the known three-copy endpoint.

## 6. The three-copy proofs do not extend in their simple forms

For orthonormal u, v put g(x, y) = ⟨x ⊗ y| ⊗_i (2F_i − I) |x ⊗ y⟩. This is 2^n times the bilinear form of q_n on
(P_x, P_y). Let p_k be the probability of exactly k local antisymmetric SWAP outcomes on u ⊗ v, and A_k(w) the same on
w ⊗ w. Then exactly

    −g(u,v) − 1 = Σ_{j≥1} (9^j − 1)(3 p_{2j+1} − p_{2j}),    g(w,w) − 1 = Σ_{j≥1} (9^j − 1) A_{2j}(w) ≥ 0,

and EP_n for Hermitian C is |g(u,v)| ≤ √(g(u,u) g(v,v)). The margin conjecture of Section 3 reads
|g(u,v)| ≤ 1 + √((g(u,u) − 1)(g(v,v) − 1)).

- **Wu and Zou's Plücker inequality** is the case j = 1 at n = 3: 3p₃ − p₂ ≤ √(a(u) a(v)).
  - Level by level it is false at n = 4: u = |0000⟩, v = |1111⟩ give 3/8 > 0, although the endpoint holds there.
  - Splitting by local pattern, by pattern size, by level with Cauchy–Schwarz, or cut by cut (one-site and 2|2
    linear entropies) is false too. The violations are 0.947, 0.952, 0.95 and 2–3 at n = 4.
  - The self terms are exact: g(w,w) − 1 − a Σ_k e_k − (2 − a) Σ_c e_c = 8(10 − a) A_4(w).
  - (`code/three_copy_mechanisms/`, `logs/levelwise_n3.log`, `jquick_n4.log`, `cutwise_I1_a2.log`,
    `plucker_n4_*.log`.)
- **Costa Rico's opposite-sign proof** amounts to the uniform bound g(u,v) ≤ 1. Since g(w,w) ≥ 1, that bound at n
  copies would give the opposite-sign case at n copies. It holds at n = 2, where g = 1 − 2(t₁ + t₂) with t_i = Tr
  ρ^u_i ρ^v_i ≥ 0, and at n = 3. But its maximum is exactly 5 at n = 4 and at least 7 at n = 5 (`logs/os_n4.log`,
  `os_n5.log`).

So from four copies on, both signs need the full coupling √(g(u,u) g(v,v)). The relaxation of Section 3 supplies it
numerically.

## 7. Status

- **Exact (computer-assisted).** The four-copy endpoint for rank-two projectors (Section 3), and the pipeline checks
  at two and three copies.
- **Numerical only** (SDP bounds at solver precision):
  - the opposite-sign case at four copies (strictly feasible after facial reduction, margin 1.9e-3; one of its zero
    families was found by local search and has no rational parametrisation yet);
  - the all-ratio lifted certificate at four copies;
  - everything at five copies.

Degenerate eigenvector bases differ between processes, so every exact step uses canonical exact bases.

The general rank-two case at four copies needs cross terms twirl(|α₁⟩⟨α₂| ⊗ |β₂⟩⟨β₁|). We estimate about 5–10k
parameters at n = 3 and about 50k at n = 4, beyond our machine. The Hermitian case at n = 5 would need about 6k
parameters and more than 20 GB. The endpoint for every n is open.

## Code and logs

| File | What it does |
|---|---|
| `code/relaxations/alg4.py`, `rep4.py` | The algebra A^⊗n, its blocks, and the one-letter relaxation of Sections 1 and 4 |
| `code/relaxations/validate.py`, `validate_spec.py` | True-point validation for the one-letter relaxation |
| `code/relaxations/dual_check.py`, `fr.py` | Dual checks and facial reduction for the one-letter relaxation |
| `code/relaxations/lift_identity.py` | Section 2: covariance blocks as partial transposes of Tomar's blocks |
| `code/relaxations/herm2.py`, `validate_herm2.py`, `ablate_herm2.py` | Two-letter relaxation (rows form), n ≤ 3 |
| `code/relaxations/herm4.py`, `validate_herm4.py`, `grid_herm4.py`, `check_n4.py` | Symmetry-reduced two-letter relaxation, the n = 4 grid, dual bounds and validation |
| `code/relaxations/fr_herm4.py`, `fr2_herm4.py`, `zero_rot.py`, `zero_interp.py` | Facial reduction at the zeros, and generators of zero families |
| `code/relaxations/exact_site.py`, `exact_maps.py` | Exact rational per-site irreps and maps (python-flint) |
| `code/relaxations/exact_core.py`, `exact_coreN.py`, `exact_pilot_n2.py`, `exact_pilot_n3.py` | Exact certificates at n = 2, 3 |
| `code/exact4/` | The exact certificate at n = 4 for rank-two projectors: engine (`core4.py`, `engine4.py`, `symm4.py`), data (`build_lin4.py`, `points4.py`, `face4.py`, `facedata4.py`, `span4.py`), certificate (`cert4.py`), checks (`verify4.py`, `crosscheck4.py`, `test_engine.py`, `points4_audit.py`); saved data and certificate as pickles |
| `code/relaxations/local_herm.py`, `local_min.py` | Local searches for the true minima |
| `code/relaxations/lift2.py`, `lift2s.py` | The all-ratio lifted certificate (Clarabel; lean assembly with SCS) |
| `code/relaxations/audit_above_endpoint.py`, `ablate_n4.py` | The control at w = 0.55, and constraint ablations |
| `code/clicks/` | Section 5 |
| `code/three_copy_mechanisms/` | Section 6 |

Main runs, from `code/relaxations/`:

    python rep4.py 4 --cons QTR1                                   # logs/n4_qtr1only.log
    python rep4.py 4 --cons LE,PT1,PT2,PT1LE,PT2LE,QTR1,QTR2,PROJ  # logs/n4_full.log
    python rep4.py 5 --cons QTR1,PROJ                              # logs/n5_qtr1_proj.log
    python grid_herm4.py 4 -1,1,-0.5,0.5,-0.25,0.25,0              # logs/herm4_n4_grid.log
    python check_n4.py 2                                           # logs/check_n4.log
    python lift2.py 3                                              # logs/lift2_n3.log
    python audit_above_endpoint.py 0.55                            # logs/audit_above_endpoint.log
    python ablate_n4.py ORTH,MARG,PT                               # logs/ablate_n4_pt.log
    python exact_pilot_n3.py                                       # logs/exact_pilot_n3.log (CERTIFIED)

and from `code/exact4/`:

    python verify4.py p_ent p_ent2                                 # logs/exact4/verify4.log (VERIFIED)

The scripts need Python 3 with NumPy, SciPy and Clarabel. The searches in Sections 5 and 6 also use PyTorch, and
the exact scripts use python-flint.

## References

- [T1] T. S. Tomar, Undistillability of Werner states from three, four and five copies, Zenodo 23718035 (2026).
- [T2] T. S. Tomar, Equality cases and every-copy results for the Werner endpoint inequality, Zenodo 23718037 (2026).
- [WZ] T. Wu, Q. Zou, Sharp Plücker geometry for three-copy Werner distillation, arXiv:2608.02647.
- [CR] P. Costa Rico, New partial trace inequalities and distillability of Werner states, Lett. Math. Phys. 115, 47
  (2025); arXiv:2310.05726.
- [BGH] K. Bharti, R. Gajjala, T. Haug, Two-copy nondistillability of Werner states: sharp partial-trace inequalities
  and finite-copy extensions, arXiv:2607.24479.
- [FHPV] Fraser, Huber, Pozsgay, Vona, arXiv:2607.24309.
- [Cl] L. Clarisse, The distillability problem revisited, Quantum Inf. Comput. 6, 539 (2006); arXiv:quant-ph/0510035.
- `first-order-bound/` in this repository: the margin conjecture q_n(C) ≥ 2^{−n}(σ₁ − σ₂)².
