# Tensor powers of Pauli maps: decomposability and the GHZ condition (floating point)

Numerics behind the remark after Theorem 7.18 of [`whole-problem-routes/README.md`](../../README.md). A Pauli map Υ has
eigenvalues 1, λ₁, λ₂, λ₃ on I, X, Y, Z; its Choi matrix is Bell-diagonal, and decomposability of Υ^{⊗n} is a linear
program over the 4ⁿ strings of Bell labels (the partial transpose acts on Bell-diagonal operators by a fixed 4 × 4
matrix on each copy). The GHZ condition is eq. (7) of Filippov and Magadov, J. Phys. A 50, 055301 (2017), a necessary
condition for positivity of Υ^{⊗n}.

- `pauli_general.py n seed trials`: random maps that are neither completely positive nor completely copositive;
  counts of (GHZ condition, a single tensor-power certificate, decomposability) (`logs/pauli_general.log`, n = 2, 3).
- `fm_vs_lp.py n trials seed`: GHZ condition against decomposability (`logs/fm_vs_lp_n4.log`, n = 4).
- `fm_boundary.py n rays seed`: along random rays from the identity map, decomposability 0.1% inside and 0.1%
  outside the GHZ boundary (`logs/fm_boundary.log`, n = 3, 4).
- `ts_dopt.py n seed trials`: for eigenvalues (1, t, 0, s) at the GHZ boundary, the best single tensor-power
  correction D = dⁿ in P = (F + D − ΓD)/2 (`logs/ts_dopt.log`): enough at n = 2, not at n = 3.

In every sample the GHZ condition and decomposability agree. This is evidence, not a proof; the proved case is the
family (1, t, 0, t) of Theorem 7.18.
