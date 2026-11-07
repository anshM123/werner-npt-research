# The permutation-symmetric sector

Checks for Proposition 7.16 of [`whole-problem-routes/README.md`](../../README.md): 𝒫_{3/4}^{⊗n} on inputs and outputs in
Sym^n(C³), the symmetric subspace of (C³)^{⊗n}.

## The certificate

Write C_{3/4} = I/4 + A/2 + B/2 on one pair (input, output), with A = (I − F)/2 and B = A^Γ = (I − |Ω⟩⟨Ω|)/2. Expanding
C_{3/4}^{⊗n} gives words in the letters I, A, B, and the compression of a word to Sym^n ⊗ Sym^n depends only on its letter
counts (k, s, t). The partial transpose Γ on the outputs exchanges s and t. The greedy operator Q_n is the compression of
the sum of the words with s > t, plus half of those with s = t ≥ 1 (with their weights 4^{−k}2^{−s−t}). Then

    Q_n + Γ(Q_n) = (Π ⊗ Π) C_{3/4}^{⊗n} (Π ⊗ Π) − 4^{−n} Π ⊗ Π,

so for real symmetric x, y one has ⟨y|𝒫_{3/4}^{⊗n}(xxᵀ)|y⟩ − 4^{−n}|x|²|y|² = 2⟨x ⊗ y|Q_n|x ⊗ y⟩. Positive semidefiniteness
of Q_n on Sym^n ⊗ Sym^n therefore proves the bound, and the kernel of Q_n gives the equality cases. In integers,
2·4^n·Q_n is the sum of the words in I, I − F, I − Ω with weight 2 for s > t and 1 for s = t ≥ 1.

Q_n commutes with rotations (the same rotation on all inputs and outputs). So it is positive semidefinite if and only if
its restriction to the vectors with J_z = 0 is, and the kernel is Sym^{2n}(R³) exactly when the kernel on that slice has
dimension n + 1.

## Two exact checks, written separately

- `verify_bosonic.py n`. The matrix of an operator in the basis f_c ⊗ f_d (f_c = the sum of the basis tensors with
  content c) is the coefficient array of its coherent-state symbol ⟨u^{⊗n} ⊗ w^{⊗n}|O|u′^{⊗n} ⊗ w′^{⊗n}⟩. For a word with
  counts (k, s, t) the symbol is X^k(α/2)^s(β/2)^t, with X = (u·u′)(w·w′), α = X − (u·w′)(w·u′) and β = X − (u·w)(u′·w′). The
  script expands these polynomials exactly (python-flint). It checks the identity above as a polynomial identity,
  (X + α + β)^n − X^n = (T + T|_{w↔w′})/2 with T the symbol of 2·4^n·Q_n. It computes the J_z = 0 slice as an exact null
  space and tests the restricted matrix by fraction-free elimination. Logs `logs/verify_bosonic_n*.log`.
- `exact_greedy.py n1 n2 …` (with `bos.py`, `exact_t2.py`). This builds the same matrices by an integer recursion over
  the copies, with an explicit basis of the rotation-invariant polynomials, and uses exact LDLᵀ on each parity block. Logs
  `logs/exact_greedy_a.log`, `logs/exact_greedy_b.log` and `logs/rerun_exact_greedy_2_7.log`. `exact_t2.py` checks the
  other closed-form split, by the first letter (logs `logs/exact_t2_*.log`).

Results (both checks agree):

| n | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|
| J_z = 0 slice | 8 | 16 | 29 | 47 | 72 | 104 | 145 | 195 |
| greedy Q_n | PSD, kernel 3 | PSD, 4 | PSD, 5 | PSD, 6 | PSD, 7 | PSD, 8 | PSD, 9 | **not** PSD |

So the bound holds with the stated equality cases for n ≤ 8. At n = 9 the greedy split fails, and so does the
first-letter split at n = 7 and 9.

## Up to 21 copies (`many_n/`)

From n = 9 on the greedy split fails, in a single spin block (L = 0 for odd n, L = 1 for even n). It is repaired by
changing the split only on the words with at most one identity letter (two at n = 21). The certificate is still in the
span of the compressed words with Q + Γ(Q) = K_n − 4^{−n}Π ⊗ Π, and every word not listed keeps its greedy weight. The
non-greedy split values are listed in `many_n/certificates.txt` (rationals with denominators at most 200). They were found
by a preconditioned semidefinite program (`cert_lean.py`, `restricted_sdp.py`) and are not canonical: the program has
many solutions, and at n = 21 the values range up to about ±50.

`many_n/verify_exact.py n CERT 200` decides positivity exactly in a coupled angular-momentum basis. Each spin L gives one
block, indexed by the harmonic degrees of the input and the output. The entries of a·b and b†·a are signed square roots of
rationals (Wigner–Eckart, oscillator reduced matrix elements, 6j symbols). After an exact rescaling of the basis the
blocks are rational, and exact LDLᵀ decides each one. The script also checks the identity Q + Γ(Q) = K′ − I as a
polynomial identity in the three contractions. `bash many_n/run_many_n.sh` reruns everything for n = 9, …, 21 (logs
`many_n/logs/`). For every n ≤ 21 the certificate is positive semidefinite with kernel exactly Sym^{2n}.

How far this is checked:
- n ≤ 9: two separate exact checks. The coupled-basis verifier reproduces the earlier exact results for the greedy split
  (positive semidefinite for n ≤ 8, not at n = 9). The n = 9 certificate also passes `verify_bosonic.py`, adapted to read
  a certificate (`many_n/verify_bosonic_cert.py`, `many_n/logs/verify_bosonic_cert_n9.log`).
- n = 10, …, 21: one exact check, by the coupled-basis verifier. Its formulas are uniform in n. `check_cpl.py` and
  `check_words.py` compare them with the monomial-basis code of this folder for n ≤ 5 and every word, and with the known
  spectra of the ladder operators.

Margins (floating point): the least eigenvalue of 4^n Q on the complement of Sym^{2n} is 1.09 at n = 9 and decreases
slowly to 0.97 at n = 21. No certificate of this kind can have it above 2. Indeed, for v a unit vector and c ⊥ v, the
vector ψ = Sym(v^{n−1}c) ⊗ vⁿ − vⁿ ⊗ Sym(v^{n−1}c) satisfies ⟨ψ|4^nQ|ψ⟩ = 2|ψ|² for every certificate. This follows from
expanding the product value to second order at real products, since Q kills Sym^{2n}. The best constant of these
decompositions equals the least value of Tr(K_n ρ) over states ρ on Sym^n ⊗ Sym^n with positive partial transpose.
Above n = 21 the certificate search ran out of memory; no rule valid for every n is known.

## Numerics (floating point)

- `num_min.py` (`logs/num_min.log`): alternating minimisation over symmetric x, y, n = 2…7, at a_n ± 0.002 and 3/4. The
  least value equals min{witness value, product value} to 10⁻¹⁶ in every case.
- `t3num.py` (`logs/t3num.log`): x symmetric, y arbitrary, n = 2…6, the same three values of a. Same result.
- `sdp_o3b.py`, `sosthr.py` (with `o3red.py`, `ybasis.py`, `zbasis.py`; logs `logs/sdp_o3b.log`, `logs/sosthr.log`): the
  best constant t in a decomposition K_n(a) − tI = Q + Γ(Q), Q ⪰ 0, on Sym^n ⊗ Sym^n (SCS, accurate to about 10⁻⁵). At
  a = a_n it is negative, about −0.011, −0.014, −0.0048 and −0.0019 for n = 2, …, 5. Decompositions certify positivity
  only from a ≈ 0.6232, 0.6775 and 0.6882 for n = 2, 3, 4 (against a_n = 0.6036, 0.6469 and 0.6704).
- `split9.py`, `split_scan.py` (logs `logs/split7.log`, `logs/split9scan.log`): a search over the weights given to the
  pairs of words (s, t), (t, s). At n = 9 it finds splits whose restricted matrix is positive definite once the forced
  kernel is removed (least eigenvalue about +5.9 in the scan's units), so this kind of certificate is not ruled out at
  n = 9. No split rule valid for every n is known.
