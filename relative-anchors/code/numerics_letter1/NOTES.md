# Floating-point evidence for the letter-1 blocks

6 November 2026. Ansh Mishra and Aryan Senthilkumar. Exploration only: no statement of the README rests on these runs.

The letter-1 block of (R∞)_j involves only the letters 0 and 1, so it is also the qubit operator of the Remark after
Lemma G of the README. `l1float.py` builds it in floating point from Lemma B.

| File | What it does | Result |
|---|---|---|
| `validate.py` | `l1float.py` against the exact builder `../middle_j1/j1lib.py`, entry by entry, at general rational weights | k = 4, …, 7: largest difference 3.3·10⁻¹⁶; the known kernel p ⊗ H_j annihilated to 10⁻¹⁵ (`logs/validate.log`) |
| `randsearch.py` | For k = 9, 10, 11 and every 1 ≤ j ≤ k − 2: the smallest eigenvalue of 2^Z Y above the known kernel (p ⊗ H_j, plus the vector of Proposition K2 at j = k − 2), at random weights. The samples come from four families: log-normal; log-uniform over e^{±6}; equal weights with one weight scaled by e^{±8}; and two tiny weights | Positive in every sample (40, 24 and 8 samples per j); kernel eigenvalues below 10⁻¹³ in modulus (`logs/rand_k9.log`, `logs/rand_k10.log`, `logs/rand_k11.log`) |
| `gap.py` | Two weight values: the same smallest eigenvalue at q = 1, 1/2, 1/4, 1/10, 1/100, for k = 5, …, 9 and a few level sizes | No decrease with k at fixed q; roughly proportional to q for small q once the lighter level has at least two copies (`logs/gap.log`) |
| `bethe_fit.py` | Lemma BA: fits the Bethe states B(x_1)⋯B(x_c)\|vac⟩ to every eigenvector of E_c (k = 6 with c = 2, k = 7 with c = 3), then evaluates the Bethe equations and eigenvalues at the fitted roots | All eigenvectors fitted; Bethe equations to 3·10⁻⁷ or better; eigenvalues 2\|s\| − 2Σx_a to 10⁻¹³ (`logs/bethe_fit.log`) |
| `bethe_control.py` | The same fit after a random symmetric perturbation of E_c | Residuals of the size of the perturbation (`logs/bethe_control.log`) |

The small eigenvalues in `logs/rand_k*.log` occur at extreme weight vectors with some weights near 0. There the
gap is expected to go to 0, as described in Section 3e of the README.
