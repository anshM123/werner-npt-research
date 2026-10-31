# From m copies to n copies: an undistillability budget for Werner states

Ansh Mishra and Aryan Senthilkumar

This folder contains a short paper, with verification code, on multi-copy undistillability of Werner states. The
question behind it is whether there exist NPT states that cannot be distilled from any number of copies. This is
Problem 4 of P. Horodecki, Ł. Rudnicki and K. Życzkowski, *Five open problems in quantum information theory*, PRX
Quantum **3**, 010101 (2022).

## Background

The Werner state ρ_α ∝ I − αF on C^d ⊗ C^d is NPT for α > 1/d and one-copy distillable for α > 1/2.
Undistillability of the states in between is known in the following cases:
- **two copies at the endpoint α = 1/2, in every dimension:** settled in July 2026 by four independent groups (Fu–Gao–Park,
  Song–Chen, Fraser–Huber–Pozsgay–Vona, Bharti–Gajjala–Haug);
- **three copies at the endpoint for d = 3 and d = 4:** settled by exact computer certificates (an anonymous Zenodo
  deposit, and A.-X. Liu with C.-F. Qiao), replayed by MacCormack and Bao, who also formalised the d = 3 case in
  Lean and proved the sharp three-copy curve (1 − 2α)(1 − α)² for d ≤ 4;
- **three copies for α ≤ 37/80 in every dimension:** MacCormack and Bao;
- **three copies at the endpoint in every dimension, four and five copies close to it:** T. S. Tomar (Zenodo, October
  2026), by exact computer certificates.
  - Four copies: α ≤ 79/160 in every d, and α ≤ 0.4979 for qutrits.
  - Five copies: α ≤ 9/20 in every d, and α ≤ 0.47 for qutrits.
  - k copies: α ≤ 11/(5k) for k ≥ 6, and α ≤ 56/(25k) for k ≥ 9, in every d.

None of these methods says anything about more copies than the one it was built for, except Tomar's k-copy bounds.
Those are of order 2.2/k, so they give NPT windows only when k is below about 2.25d. The classical result of
DiVincenzo, Shor, Smolin, Terhal and Thapliyal (2000) gives, for each n, some non-explicit neighbourhood of the PPT
boundary that is n-copy undistillable. Its size was never determined. Qi, Gupur, Wu and Guo (arXiv:2402.18037)
and Bharti, Gajjala and Haug give explicit bounds that do not depend on d (about 0.48/n in the first case). Because
they are dimension-free, they only produce NPT windows when d is roughly larger than 2n.

## Result

**m-copy undistillability is a budget that can be spent across copies.**

Write κ = (1 − α)/(1 + α) and κ_PPT = (d − 1)/(d + 1). Suppose ρ_{α₀} is m-copy undistillable. Then a product
ρ_{α₁} ⊗ … ⊗ ρ_{αₙ} with all α_i ≤ α₀ cannot be distilled by any filter of Schmidt rank two whenever

    Σ_i (κ_PPT − κ_i)₊ ≤ m (κ_PPT − κ₀).

Separable copies cost nothing, and a copy of ρ_{α₀} costs exactly one unit. For identical copies, ρ_α is n-copy
undistillable for

    α ≤ [n(1 + α₀) + m(dα₀ − 1)] / [nd(1 + α₀) − m(dα₀ − 1)],

which at α₀ = 1/2 reads (3n + m(d − 2)) / (3nd − m(d − 2)).

Consequences:
- Qutrit Werner states with 1/3 < α ≤ 5/11 are NPT and four-copy undistillable, and α ≤ (n+1)/(3n−1) is n-copy
  undistillable.
- Ququart Werner states with 1/4 < α ≤ 3/7 are NPT and four-copy undistillable.
- In every dimension, the two-copy theorem alone gives four-copy undistillability for 1/d < α ≤ (d+4)/(5d+2).
- With Tomar's three-copy endpoint theorem as input, ρ_α is n-copy undistillable for α ≤ (n+d−2)/(nd−d+2), for
  every n ≥ 3 and every d.
- With Tomar's five-copy qutrit bound (α₀ = 47/100) as input and correlated separable factors, qutrit Werner states
  are six-copy undistillable for α ≤ 0.4469; the budget alone gives 1087/2441 ≈ 0.4453. Together with the correlated factors below, they are n-copy undistillable for
  α ≤ (851n + 1272)/(2553n − 1272) for every n ≥ 24, that is 1/3 + 0.664/n.
- For ququarts, his every-dimension five-copy bound (9/20) gives α ≤ (12913n + 39648)/(51652n − 39648) for n ≥ 24,
  that is 1/4 + 0.960/n.

All of these windows contain NPT states for every n at fixed d. Tomar's own k-copy bounds cannot give that.

The statement holds for arbitrary states, not only Werner states. If A is separable, B is m-copy undistillable and
A + vB is separable, then ⊗_i(A + u_iB) is one-copy undistillable whenever Σ_i(1 − v/u_i)₊ ≤ m.

A corollary holds for every number of copies at once. Call a site of a filter *thin* if its local supports overlap
in at most two dimensions (Tr(P_i Q_i) ≤ 2). A Werner factor with α ≤ 1/2 at a thin site becomes separable after
compression, so only thick sites spend budget. (This thin-site copy removal is Wu–Zou's Thm. 7.2 for three copies and
Tomar's Lemma 6.1, Zenodo 23718037, for every number of copies; the combination with the budget is new.) Hence, for every n, no filter with at most t thick sites distills
ρ_α^{⊗n} when α lies in the t-copy window; for qutrits this means α ≤ (t+1)/(3t−1).

More precisely, a site whose local supports overlap by τ costs (1 + α₀)(ατ − 1)₊ / ((α₀τ − 1)(1 + α)). That cost is zero
whenever ατ ≤ 1. Two things follow:
- If ρ_α is distillable at all, it needs more than m(κ_PPT − κ₀)/(κ_PPT − κ) copies.
- The distilling filter must have at least that many sites with overlap above 1/α.

For a qutrit Werner state with α = 0.34, for example, that is more than 67 copies.

The proof has two parts.
- A separable factor can be absorbed into the filter, so a product with at most m non-separable factors inherits
  the m-copy theorem.
- A combinatorial lemma: the polynomial ∏_i(s_i + z_i) is a nonnegative combination of products with at most m
  pure factors z_i iff Σ_i(1 − s_i)₊ ≤ m. Sufficiency is a pipage-rounding argument on the hypersimplex; necessity
  comes from comparing the two top coefficients.

The budget is therefore optimal for decompositions into products. Separable operations can do slightly better.
At α = 83/200, just above the budget value 7/17 ≈ 0.4118, two qutrit copies are an exact convex combination of
filtered copies of the nine-dimensional endpoint state. So four qutrit copies are undistillable for α ≤ 83/200 by
the two-copy theorem alone (applied in dimension 9). The certificate is in `code/pair_certificate_83_200.json`, with
the verifier `code/verify_pair_certificate.py`. Fed back into the budget with m = 4, it gives qutrit windows
(283n + 196)/(849n − 196) for n ≥ 4 from human-proved input alone, slightly above (3n + 2)/(9n − 2).

### Correlated separable factors

The budget is optimal only among products of single-site factors. A separable factor can also be correlated across
sites. In the factor Σ, Alice holds a maximally entangled pair on her two sites and Bob holds another on his. Σ is a
product across Alice and Bob, so it is separable. After twirling, though, it puts weight (d − 1)/(2d) on "both sites
antisymmetric", where products of single-site separable factors give at most 1/4.

Twirling every site and symmetrising over sites reduce the question to convolutions of small distributions, i.e. to a
linear program. Exact rational certificates (`code/correlated/`) then give:
- four qutrit copies up to α = 0.45776 from the three-copy theorem, and up to 0.41603 from the two-copy theorem
  alone (budget: 5/11 and 7/17);
- four ququart copies up to 0.43574 (budget: 3/7);
- certificates for n = 4, 6, 8, 10, 12, 16 and 24 copies, and for ququarts up to n = 24 (the excess over the
  PPT threshold is about 16% larger at large n).

Fed back into the budget, the n = 24 certificates enlarge the excess α_n − 1/3 of every window with n ≥ 24 by about
10%. For qutrits from
three copies this gives α ≤ (135403n + 149016)/(406209n − 149016), so a qutrit state at α = 0.34 would need more than
73 copies.

The gain grows with the dimension. Using only the two-copy theorem, four copies are undistillable for α ≤ 59/250 =
0.236 in every dimension. The previous dimension-free bounds were about 0.12. The budget alone gives (d+4)/(5d+2),
which tends to 1/5. The proof combines one exact certificate at d = 1000 (smaller d by compression) with a
certificate whose weights are explicit rational functions of d, positive for all d ≥ 1000
(`code/correlated/uniform_d_symbolic_weights_59_250_D1000.txt`).

Fed into the budget, this gives dimension-free n-copy windows α ≤ 236/(309n − 236) ≈ 0.76/n. An eight-copy uniform
certificate (α ≤ 0.1195 in every d) improves this to 1912/(2239n − 1912) ≈ 0.85/n for n ≥ 8, and a twelve-copy one
(α ≤ 0.0797 in every d) to 9564/(10797n − 9564) ≈ 0.89/n for n ≥ 12, against about 0.48/n previously. With MacCormack and Bao's every-dimension three-copy bound as input, the four-copy window is 0.335 in every
dimension, and the dimension-free n-copy windows become 268/(267n − 268) ≈ 1.0/n.

Every decomposition of this type, with at most m endpoint factors and any separable remainder, obeys
n(ρ − b)/(ρ + b) ≤ m(3 − b)/(3 + b), where ρ = (1+α)/(1−α) and b = (d+1)/(d−1). This ceiling also tends to 1/d.

The windows shrink to 1/d as the number of copies grows, so this does not settle the all-copy question.

**What Tomar's results supersede.** His computer-assisted every-dimension bounds are stronger than the following
every-dimension statements above:
- the four-copy window 0.236;
- the dimension-free windows ≈ 0.76/n, 0.85/n, 0.89/n and 1.0/n;
- the MacCormack–Bao-based values.

These statements remain correct. They are kept because of their weaker inputs: the first three use only the
human-proved two-copy theorem, and the last ones use the three-copy bound of MacCormack and Bao. At fixed d the
budget windows are still the only ones that stay above 1/d for every n.

## Status

- The proofs are complete and elementary.
  - The budget part has been checked independently twice, step by step and in exact rational arithmetic, the second
    time on about 31,000 randomized instances of the cone lemma (that check's code is not in this repository;
    `code/pipage_constructor.py` runs 400 instances per seed).
  - The correlated-factor part was checked independently once more. That check re-verified all certificates with
    separate code and re-derived the reduction, the ceiling and the uniform-in-d arguments.
  - The paper has not been refereed.
- The input theorems for m ≥ 3 are computer-assisted and also unrefereed. That includes Tomar's three-, four- and
  five-copy bounds. The m = 2 input has four independent human proofs.
- The revision that adds Tomar's inputs (30 October 2026) was checked independently:
  - all certificates with a non-endpoint core were re-verified with separate code;
  - the all-n formulas were re-derived;
  - every quoted number from Tomar's paper was compared with the source.
- The correlated-factor windows are exact rational certificates: computer-generated, but checkable by hand in
  principle, being identities between n + 1 rational numbers. In the certificate files the core factor is always
  called `E`, and an optional field `core_r0` gives its parameter (the default 3 is the endpoint).

## Contents

- `paper/`: LaTeX source and PDF.
- `code/verify_qutrit_certificate.py`: exact check of the explicit four-copy qutrit certificate at α = 5/11 on all 16
  parity sectors.
- `code/pipage_constructor.py`: runs the constructive proof of the cone lemma in exact arithmetic on random instances,
  and checks nonnegativity of the weights, the generator constraint and the polynomial identity.
- `code/hetero_cone_lp.py`: LP test of the budget boundary, feasible on the boundary and infeasible just beyond it.
- `code/moment_recursion_check.py`, `code/taylor_lemma_check.py`: checks of the alternative proof for identical copies.
- `code/overlap_budget_check.py`, `code/copies_needed_check.py`: checks of the local-overlap cost formula and of the
  copy lower bound.
- `code/copy_window_check.py`, `code/four_copy_cone_lp.py`: exact and LP computation of the windows from the
  product-decomposition cone.
- `code/pair_certificate_83_200.json`, `code/verify_pair_certificate.py`: the beyond-products certificate at
  α = 83/200 and its standalone verifier, in pure Python with exact fractions.
- `code/pair_sim_extract.py`: the search that produced that certificate, which also verifies it exactly.
- `code/correlated_core/`: certificates with a core other than the endpoint, `corr_cert_d{d}_m{m}_a{p}_{q}_n{n}.json`.
  The core is ρ at α₀ = p/q, for example Tomar's 47/100 and 9/20.
  - Generators: `corr_budget_cert_core.py` and `exact_simplex_cert_core.py`. Set the environment variable
    CORE_ALPHA to α₀.
  - Checker: `verify_core_certs.py`, an exact check of every file in the folder and of the all-n formulas. Run it from
    inside the folder.
- `code/correlated/`: correlated-factor certificates `corr_cert_d{d}_m{m}_n{n}.json`; the generators
  `corr_budget_cert.py` (LP plus exact solve) and `exact_simplex_cert.py` (exact rational simplex, used for n = 8, 16, 24); the standalone exact checker `verify_corr_cert.py` (pure Python, with an optional dense check
  of the factor weights); the threshold LP `named_factor_ncopy.py`; and `corr_sep_lp.py`, a search over arbitrary
  correlated separable factors (needs PyTorch). The every-dimension certificates are
  `uniform_d_symbolic_weights_*.txt` (weights as rational functions of d, produced by `uniform_d_cert.py`) together
  with the exact certificates at d = 1000, 10^5 and 10^6. `verify_uniform_symbolic.py`, written separately from the
  generator, rechecks each symbolic file: the identity in d, and positivity of every weight for d ≥ D0 by exact root
  counting. `dinf_limit_lp.py` is the d → ∞ limit of the same LP (numerical). With two-copy cores it levels off at
  nα ≈ 0.96, against the ceiling nα ≤ 1.

Everything except `pair_sim_extract.py` and `correlated/corr_sep_lp.py` needs only Python 3 with NumPy, SciPy
and SymPy; the plot also uses Matplotlib. Those two additionally need PyTorch.

Typical runs, from `code/`:

    python verify_qutrit_certificate.py
    python verify_pair_certificate.py
    python pipage_constructor.py
    python hetero_cone_lp.py 4 2 5
    python four_copy_cone_lp.py 4 3 3
    python correlated/verify_corr_cert.py correlated/corr_cert_d3_m3_n4.json --dense
    python correlated/verify_corr_cert.py correlated/corr_cert_d1000000_m2_n12.json
    python correlated/verify_uniform_symbolic.py correlated/uniform_d_symbolic_weights_797_10000_D1000000_n12.txt 12 2 3 1 797 10000 1000000
    cd correlated_core && python verify_core_certs.py

## License

MIT.
