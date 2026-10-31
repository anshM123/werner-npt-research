# A thickness margin for the two-copy Werner inequality

Ansh Mishra and Aryan Senthilkumar, October 2026.

## The result

The Werner state ρ ∝ I − F/2 is two-copy undistillable. Equivalently, every operator C of rank at most two on
H₁ ⊗ H₂ satisfies

    q₂(C) = ‖C‖₂² − ½(‖Tr₁ C‖₂² + ‖Tr₂ C‖₂²) + ¼|Tr C|² ≥ 0.

The equality cases have co-ranges whose one-copy reductions have rank at most two (T. S. Tomar). Tomar conjectured a
margin linear in the third eigenvalue λ₃ of those reductions, with constant 1/4. He noted that for two copies no
linear bound of this kind was known.

We prove

    q₂(C) ≥ min{τ₁, τ₂} ‖C‖₂² / (48 − 8√2),

where τ_l is the sum of all but the two largest eigenvalues of the reduction of the co-range projector on factor l.
Since τ_l ≥ λ₃, this gives Tomar's m₂(V) ≥ min_l λ₃⁽ˡ⁾(V)/37 for every plane V, in every dimension.

The main new ingredient is a sharp inequality. If A ≥ 0 on C² ⊗ K has qubit marginal I₂, then

    λ_max(A) ≤ 1 + ‖A₀₀ − A₁₁‖₂ / √2.

Its proof is an explicit dual certificate that turns out to be a rank-two projector. A refinement of the same
certificate brings in the tail τ of the reduced operator A₀₀ + A₁₁ on K beyond its two largest eigenvalues:

    ‖A₀₀ − A₁₁‖₂² ≥ 2(λ_max(A) − 1)(λ_max(A) − 1 + τ).

This gives the best possible constant 1/2 in the lower bound for 4‖Π(w₁ ⊗ w₂)‖² (Lemma 4).

If the range and co-range of C lie in one plane V, which includes every normal C, we prove Tomar's inequality with the
best possible constant, which beats his conjectured 1/4 (Theorem 5):

    q₂(C) ≥ (1 − 1/√2) min_l λ₃⁽ˡ⁾(V) ‖C‖₂²     (sharp in every dimension d ≥ 3),
    q₂(C) ≥ ((3 − √5)/4) max_l τ_l(V) ‖C‖₂²     (sharp as d → ∞),

and for C ≥ 0 also q₂(C) ≥ ¼ min_l τ_l(V) ‖C‖₂² (sharp as d → ∞). Tomar's extremal example is Hermitian and attains
the first bound. The key step is the identity q₂(C) = q₂(Re C) + q₂(Im C), which holds for every C. Wu and Zou use the
same splitting for three copies (arXiv:2608.02647, Cor. 6.2). In the first bound λ₃ cannot be replaced by τ: that
fails for indefinite Hermitian C from d = 4 on.

Tomar proved the same constant 1 − 1/√2 on a different class: every C whose range lies in a plane spanned by two
vectors with common Schmidt bases (his Theorem 7.2, Zenodo 23718037). Our Theorem 5 allows every plane, but needs both
the range and the co-range of C to lie in it.

## What this does not do

- For C whose range and co-range differ, the constant 1/37 is far from optimal. Numerically the best constant at
  d = 3 is 1 − 1/√2 ≈ 0.293, so Tomar's 1/4 appears to be correct for two copies, but we do not prove it in general.
- Numerically the margin is linear in max{τ₁, τ₂} for all C, with the constant (3 − √5)/4 ≈ 0.191 of Theorem 5,
  attained at the plane span{|00⟩, (|01⟩ + |12⟩)/√2} by an operator whose range is not in the plane. We prove this only
  for the plane class.
- The argument is specific to two copies. Undistillability of ρ from every number of copies remains open.

## Status

- Each step was checked by an independent re-derivation, with symbolic checks of the certificate. The scripts in
  `code/` run a grid check and random sampling; the adversarial searches (SDP, local optimisation) done during the
  independent check are not part of this folder. No counterexample was found.
- The note has not been refereed.

## Contents

- `paper/`: LaTeX source and PDF (8 pages).
- `code/check_thickness.py`: checks the certificate of Lemma 2 exactly on a grid, and tests Lemma 2, Lemma 4,
  Theorem 1 and Theorem 5 on random instances, including near-equality and near-degenerate ones. It needs Python 3
  with NumPy and runs in under a minute:

      python check_thickness.py

## License

MIT (see `../LICENSE`).
