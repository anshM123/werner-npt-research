# The bosonic sector of the Werner endpoint

Ansh Mishra and Aryan Senthilkumar, 30 October 2026

These notes study the all-copy endpoint form on operators whose range and co-range lie in the symmetric subspace. We
call this the bosonic sector. Results are labelled **proved** or **numerical**. The scripts in `code/` check every
formula and claim below, and their output is in `logs/`.

## Setting

For an operator C on (C^3)^{⊗n} and 0 ≤ t ≤ 1/2, let

    q_n(t, C) = Σ_J (−t)^{|J|} ‖Tr_J C‖_2²,

with J running over subsets of the n sites. The qutrit Werner state at α = 1/2 is n-copy undistillable if and only if
q_n(1/2, C) ≥ 0 for every C of rank at most two. The Schmidt-rank-two criterion is due to [DSSTT, DCLB], and this
partial-trace form to [CR23, Thm. 1] (see also [QGWG, Thm. V.1]; restated as [SC, Thm. 1]).

The bosonic sector is the set of C with P C P = C, where P projects onto Sym^n(C^3). It is the natural test class
for every n at once. The all-copy statement is still open even here.

## 1. Formulas (proved)

**(a) Annihilation form.** Let a_i : Sym^m → Sym^{m−1} be the annihilation operators in the occupation basis, and
A(C) = Σ_i a_i C a_i^*. For bosonic C every partial trace over j sites equals ((n−j)!/n!) A^j C, so

    q_n(t, C) = Σ_j (−t)^j (n−j)!/(j! n!) ‖A^j C‖².

**(b) SU(3) spectrum.**
- End(Sym^n C^3) is the direct sum of the irreducible components V_(k,k), k = 0, …, n, each of dimension (k+1)^3
  and multiplicity one.
- A^*A acts on the V_(k,k) component of End(Sym^m) as μ_{m,k} = (m−k)(m+k+2).
- Hence q_n(t, C) = Σ_k γ_{n,k}(t) ‖Π_k C‖², with

      γ_{n,k}(t) = Σ_j (−t)^j C(n−k, j) C(n+k+2, j) / C(n, j).

- At t = 1/2:
  - γ_{n,n} = 1;
  - γ_{n,n−1} = −1/(2n) and γ_{n,n−2} = 1/(2(n−1));
  - γ_{n,0} = 2^{−n−2} (1 + (−1)^n (2n+3)).
- The sign of γ_{n,k} is (−1)^{n−k}. We checked this exactly for n ≤ 39; we have not proved it for all n.
- The form is indefinite on End(Sym^n) for every n ≥ 1, since γ_{n,n−1} < 0 < γ_{n,n}. So a proof must use the rank.

**(c) Pair form.** Pair site i of the range with site i of the co-range, so that |C⟩⟩ ∈ (C^3 ⊗ C^3)^{⊗n}. Put
φ = Σ_i |ii⟩/√3 and B = I − 3t |φ⟩⟨φ|. For every C on (C^3)^{⊗n},

    q_n(t, C) = ⟨⟨C| B^{⊗n} |C⟩⟩,

because ‖Tr_J C‖² = 3^{|J|} ⟨⟨C| (|φ⟩⟨φ|)^{⊗J} ⊗ I |C⟩⟩. For bosonic C, |C⟩⟩ lies in Sym^n(C^9), on which B^{⊗n} acts as
the symmetric power of B.

Check: `code/check_formulas.py` covers (a) and (b), including the closed forms. The pair form (c) is checked in
`code/check_sector_theorem.py`, part (b), for a random C at n = 3.

## 2. Occupation-sector positivity (proved, every n)

Fix a unit vector z ∈ C^3 and put K = z^⊥. Let E_l be the projector onto the part of Sym^n(C^3) with exactly l
particles in mode z, that is, onto z^l · Sym^{n−l}(K).

**Theorem.** Let 0 ≤ t ≤ 1/2 and let C ∈ End(Sym^n C^3) have any rank.
1. If the range of C lies in some E_l, then q_n(t, C) ≥ 0.
2. The same holds if the co-range lies in some E_m.
3. More generally, q_n(t, C) ≥ 0 whenever, among the blocks E_l C E_m, at most one is nonzero on each diagonal
   l − m = δ.

In all three cases the bound is quantitative: q_n(t, C) ≥ (1 − 2t)^n ‖C‖². This margin is the sector analogue of the
endpoint-to-margin bounds in [T1, Lemma 2.3(b)], [T2, Lemma 6.2 and eq. (6.1)] and [MB, Remark 3.2(a)].

*Proof.*
- **Basis.** Let 0, 1 be an orthonormal basis of K. Take the orthonormal basis of C^3 ⊗ C^3

      f = (|00⟩ + |11⟩)/√2,  g = |zz⟩,  s_1 = (|00⟩ − |11⟩)/√2,  s_2 = |01⟩,  s_3 = |10⟩,
      s_4 = |0z⟩,  s_5 = |1z⟩,  s_6 = |z0⟩,  s_7 = |z1⟩.

  Then φ = √(2/3) f + √(1/3) g. So B is the identity on the span of the s_k, and on span{f, g} it is

      B_T = [[1 − 2t, −√2 t], [−√2 t, 1 − t]].

- **Action of B^{⊗n}.** Use the occupation basis |a, b, σ⟩ of Sym^n(C^9): a copies of f, b copies of g, and σ the
  occupations of s_1, …, s_7. Then B^{⊗n} maps |a, b, σ⟩ into the span of the |a′, b′, σ⟩ with a′ + b′ = a + b,
  through the symmetric power of B_T.
- **Statement 1.** The z-number on the range side of |a, b, σ⟩ is b + σ_6 + σ_7. If the range of C lies in E_l, then
  |C⟩⟩ is a combination of basis states with b + σ_6 + σ_7 = l. On these states b′ = b and a′ = a, so the
  compression of B^{⊗n} to them is diagonal. Its entries are the coefficients of f^a g^b in
  (B_ff f + B_gf g)^a (B_fg f + B_gg g)^b:

      d(a, b) = Σ_j C(a, j) C(b, j) (1 − 2t)^{a−j} (2t²)^j (1 − t)^{b−j}.

  For 0 ≤ t ≤ 1/2 every term is nonnegative. Keeping only j = 0 and using 0 ≤ 1 − 2t ≤ 1 − t ≤ 1 and a + b ≤ n gives
  d(a, b) ≥ (1 − 2t)^a (1 − t)^b ≥ (1 − 2t)^{a+b} ≥ (1 − 2t)^n. The pair vectorisation is an isometry, so
  q_n(t, C) = Σ d(a, b) |⟨⟨a, b, σ|C⟩⟩|² ≥ (1 − 2t)^n ‖C‖². The bound is far from sharp: numerically, at n = 2 and
  t = 0.4 the smallest sector minimum is 0.2, against (1 − 2t)² = 0.04 (`code/check_sector_theorem.py`, part (e)).
  It is empty at t = 1/2.
- **Statement 2** follows from q_n(t, C^*) = q_n(t, C).
- **Statement 3.** q_n is invariant under C ↦ U^{⊗n} C U^{*⊗n} with U = diag(1, 1, e^{iθ}) in the basis (0, 1, z).
  The block E_l C E_m picks up the phase e^{iθ(l−m)}. So blocks on different diagonals are orthogonal, both for the
  q-form and for the Hilbert–Schmidt inner product, and both are sums over the single blocks. Each block has range in
  one sector, so statement 1 applies to it. ∎

At t = 1/2 the entries are d(a, b) = C(b, a) 2^{−b}, and they vanish when a > b. So inside a fixed range-occupation
sector, q_n(1/2, C) = 0 exactly when every nonzero amplitude ⟨⟨a, b, σ|C⟩⟩ has a > b.

Remarks.
- **The case l = 0 needs no symmetry.** Tomar [T2, Lemma 2.6] proved that q_n(t, C) ≥ 0 for 0 ≤ t ≤ 1/2 whenever
  the range of C lies in a product W_1 ⊗ ⋯ ⊗ W_n of subspaces of dimension at most two, of any rank. The following is
  a strengthening.
  - Let S_z = {X ∈ End(C^3) : ⟨z|X|z⟩ = 0}. Then B(X) = ‖X‖² − t|Tr X|² satisfies B(X) ≥ (1 − 2t)‖X‖² on S_z, since
    Tr X = Tr X_KK and |Tr X_KK|² ≤ 2‖X_KK‖².
  - Hence q_n(t, C) ≥ (1 − 2t)^n ‖C‖² for every C in S_{z_1} ⊗ ⋯ ⊗ S_{z_n}, of any rank.
  - This includes every C whose range or co-range lies in a product of qubits K_1 ⊗ ⋯ ⊗ K_n.
- **For l ≥ 1 the sectors contain thick operators.** Their range has a one-site marginal of full rank. An example at
  n = 2 is the range span{Sym(z ⊗ e_0), Sym(z ⊗ e_1)}.
- **The argument is genuinely bosonic.** For distinguishable sites, the subspace "exactly l of the range sites equal
  z" is not positive. At n = 2, l = 1, t = 1/2, the exchange between (f at site 1, g at site 2) and (g at site 1,
  f at site 2) gives the block [[0, 1/2], [1/2, 0]]. Symmetrisation keeps only the combination with eigenvalue +1/2.
- **Blocks on one diagonal do interfere.** Two blocks on the same diagonal already give a negative eigenvalue at
  n = 2, namely −1/4. So the theorem does not give the bosonic sector by itself. Any proof for the full bosonic sector
  must use the rank to control this interference (Section 4).

Check: `code/check_sector_theorem.py`.
- It checks the diagonal structure directly for n ≤ 3: the largest off-diagonal entry is 6 × 10^{-17}, and the
  diagonal agrees with d(a, b) to 4 × 10^{-16}.
- It computes the minimum eigenvalue on every sector E_l for n ≤ 5.
- It checks all two- and three-block patterns for n ≤ 4.
- It computes the sector minima at t = 0.4 for n = 2, 3.

## 3. Numerical observations

These come from local optimisation, so they are evidence, not certificates.
- **The minimum is 0.** Over rank-two bosonic C, the minimum of q_n(1/2, C)/‖C‖² is 0, to within 10^{-15}, for
  n = 2, …, 7.
- **The zeros found are qubit-supported.**
  - At n = 3, every zero found has two equal singular values, and the one-site marginal of range plus co-range has
    spectrum (1/2, 1/2, 0).
  - Inside Sym^n(K) the form equals ‖P^{⊗n} C‖², where P projects onto the traceless operators on K, so it is
    positive semidefinite there.
  - So every rank-two C built only from SU(2) spins l with l ≢ n (mod 2) is a zero. These form a continuous family,
    which contains |0^n⟩⟨0^n| − (−1)^n |1^n⟩⟨1^n|.
  - Qubit-supported zeros are known for every n: every qubit-supported plane is a zero [T2, Lemma 2.7(a)]; see also
    the zero families I-III of [MB] and the anti-states of Wu and Zou [WZ, Prop. 8.2], whose n = 3 example is
    |000⟩⟨000| + |111⟩⟨111|. The characterisation through the SU(2) spins appears to be new.
- **The worst case is qubit-supported along the whole Werner family.**
  - For t ∈ {1/3, 0.4, 0.45, 0.48, 0.49, 1/2}, the minimum over all rank-two bosonic C equals the minimum over
    qubit-supported ones to six digits, for n ≤ 5.
  - The normalised value at |0^n⟩⟨0^n| − (−1)^n |1^n⟩⟨1^n| is (1 − t)^n − t^n = n 2^{1−n}(1 − 2t) + O((1 − 2t)²).
    For t ≥ 0.45 the minimum equals this value to six digits. For some smaller t the minimiser is a different
    qubit-supported operator, with a smaller value: n = 3, t = 1/3 (0.2349 against 0.2593); n = 4, t ≤ 0.4; n = 5,
    t ≤ 0.4.
- **The margin grows with thickness.**
  - Thickness here is the smallest eigenvalue of the one-site marginal of (P_ran + P_coran)/4. We impose
    thickness ≥ τ through a quadratic penalty, so we report the thickness each minimiser actually has.
  - At the smallest thickness sampled, the ratio of the minimum to the thickness is about 2 at n = 2 (thickness
    0.0053) and about 1 at n = 3 (thickness 0.0075). The ratio decreases at larger thickness.
  - In each minimiser the best qubit compression C_z = E_0 C E_0 has a smaller ratio q/‖·‖². This is consistent with
    the compression inequality min_z q(C_z)/‖C_z‖² ≤ q(C)/‖C‖², which would imply the bosonic sector.

What does not work:
- **Splitting Σ γ_k ‖Π_k C‖² into nonnegative adjacent pairs.** Only the top pair, ‖Π_n C‖² ≥ ‖Π_{n−1} C‖²/(2n),
  survives numerically for rank two.
- **Decomposable certificates.** One could try q_n written as a positive combination of the one-copy forms
  q_1(A^* C A), with A : C^m → Sym²(C^3), plus a positive semidefinite form. At n = 2, no such term with A ≠ 0
  vanishes on the whole zero set.
  - The zeros include |uu⟩⟨uu| − |vv⟩⟨vv| for every orthonormal pair u, v.
  - On these, A^* C A = |α⟩⟨α| − |β⟩⟨β| with α = A^*|uu⟩ and β = A^*|vv⟩, which has rank at most two.
  - q_1 vanishes on an operator of rank at most two only if it is 0 or a multiple of a rank-two projector. A
    difference of two rank-one positive operators is a multiple of such a projector only if it is 0.
  - So |α⟩⟨α| = |β⟩⟨β| for all orthonormal u, v. Hence A^*(u ⊗ u) = λ(u) w_0 with |λ| constant on unit vectors,
    where λ is a complex quadratic form on C^3.
  - A quadratic form on C^3 has an isotropic unit vector, so λ = 0. Then A^* vanishes on Sym²(C^3), and A = 0.
- **U(1) charge sectors** are indefinite.
- **Spaces of forms through m generic points** are not positive unless they are 2-dimensional.
  - Here the range is restricted to V_S = span{p^n : p ∈ S}^⊥ and the co-range is free.
  - We tried 30 random point sets S of each size m.
  - Every space of dimension at least 3 has a negative eigenvalue. The dimension-2 spaces (pencils) give small
    positive minima: +0.0019 at n = 2 and +0.0120 at n = 3.

Checks: `code/check_numerics.py`, `code/check_zeros_and_thickness.py` and `code/check_negative_results.py`.

## 4. What is open

The bosonic sector itself: is q_n(1/2, C) ≥ 0 for every rank-two C ∈ End(Sym^n C^3) and every n?

By Section 2, for each z the only possible source of negativity is interference between occupation blocks on the same
diagonal. For each spectator configuration σ, this interference is the quadratic form of the symmetric power
Sym^N(B_T) on the vector of amplitudes ⟨⟨a, N − a, σ|C⟩⟩. At t = 1/2 the eigenvalues of this power are (−1/2)^j. The
rank-two condition has to keep these vectors away from the eigenvectors with j odd.

## References

- [DSSTT] D. P. DiVincenzo, P. W. Shor, J. A. Smolin, B. M. Terhal, A. V. Thapliyal, Phys. Rev. A 61, 062312 (2000).
- [DCLB] W. Dür, J. I. Cirac, M. Lewenstein, D. Bruß, Phys. Rev. A 61, 062313 (2000).
- [CR23] P. Costa Rico, New partial trace inequalities and distillability of Werner states, Lett. Math. Phys. 115, 47
  (2025); arXiv:2310.05726.
- [QGWG] S.-Y. Qi, G. Gupur, Y.-C. Wu, G.-P. Guo, Unveiling NPT bound problem: from distillability sets to
  inequalities and multivariable insights, Phys. Rev. A 110, 012406 (2024); arXiv:2402.18037.
- [SC] Z. Song, L. Chen, A partial-trace matrix inequality and Werner-state distillability, arXiv:2607.23416.
- [T1] T. S. Tomar, Undistillability of Werner states from three, four and five copies, Zenodo 23718035 (2026).
- [T2] T. S. Tomar, Equality cases and every-copy results for the Werner endpoint inequality, Zenodo 23718037 (2026).
- [MB] I. MacCormack, N. Bao, Three copies of the Werner state: the sharp curve in dimensions three and four and a
  bound in every dimension (2026), github.com/physical-superintelligence/oqp-batch-2026.
- [WZ] T. Wu, Q. Zou, Sharp Plücker geometry for three-copy Werner distillation, arXiv:2608.02647.

## Code

The scripts need Python 3 with NumPy and SciPy.
- `code/bosonic.py`: shared helpers.
- `code/check_formulas.py`: Section 1.
- `code/check_sector_theorem.py`: Section 2.
- `code/check_numerics.py`, `code/check_zeros_and_thickness.py` and `code/check_negative_results.py`: Section 3.
