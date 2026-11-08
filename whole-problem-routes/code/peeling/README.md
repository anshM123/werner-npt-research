# Product test vectors and peeling one copy (Proposition 7.30, Remark 7.31)

Everything here is floating point. Proposition 7.30 is proved in
[`whole-problem-routes/README.md`](../../README.md); the scripts only check it numerically.

- `product_probe.py`: for random mixed and pure states, states near products and a maximally entangled pair padded with
  products (n = 1, …, 5), and for random and locally optimised real product test vectors v, checks
  ⟨v|𝒫_{3/4}^{⊗n}(ρ)|v⟩ ≥ (¾)^n 3^{−Σ⟨v_s|ρ_s|v_s⟩} ≥ 4^{−n}3^{2Στ_s} ≥ ∏(¼ + ½τ_s). Log: `logs/product_probe.log`.
- `peel_search.py n rank mode starts seed`: adversarial search (L-BFGS over ρ = GGᵀ/Tr) for the statements of
  Remark 7.31: `tc` (Conjecture TC for mixed inputs), `peel1` and `peelbest` (peeling the least eigenvalue, at copy 1 and at
  the best copy), `vtc1` and `vtcbest` (the vertex bound of Remark 7.31(ii), at copy 1 and at the best copy). A printed
  ratio below 1 is a violation. Logs: `logs/peel_search_*.log`.
- `vertex_level2.py m a_num a_den rebit|qutrit`: the real second-level relaxation for positivity of
  Sym₂ ⊗ 𝒫_a^{⊗m} (copy 0 a plane, `rebit`) or 𝒫_{1/2} ⊗ 𝒫_a^{⊗m} (copy 0 a qutrit, `qutrit`), with O(2) blocks and the
  five Brauer relations at d = 2 for the plane copy. Log: `logs/vertex_level2.log`.

`run_checks.sh` reruns all of these (a few minutes). Python with NumPy, SciPy, PyTorch and CVXPY (Clarabel).
