# Linear-programming and replica checks behind two entries of `OBSTRUCTIONS.md`

Ansh Mishra and Aryan Senthilkumar, 4 November 2026. Python 3 with NumPy, SciPy, python-flint and, for `sdp3.py`,
CVXPY with the Clarabel solver.

| File | What it does |
|---|---|
| `exact_lp.py` | Exact rational simplex (python-flint) for the weight-enumerator relaxation of Section 7 of `../../NOTES.md`, with the range required to be t-uniform (A_j = B_j = 0 for 1 ≤ j ≤ t). Reproduces the certificates of Proposition 19, for example −5515759/22425216 at n = 9, t = 1. |
| `scan_threshold.py`, `logs/scan_threshold.log` | The least t with a positive optimum, for n = 9, …, 33. |
| `rep3.py`, `validate3.py` | Three-replica bookkeeping: the twirl of Π^{⊗3} in the irreducible blocks of C[S_3] at every copy. Checked against random rank-two projectors at n = 1, 2, 3: block positivity, Λ³ = 0, global replica symmetry, Tr T = 8, and the endpoint value. |
| `sdp3.py` | The three-replica relaxation (copy-symmetrised blocks, Λ³(ran Π) = 0, moments, A_j ≥ 0, B_j ≥ A_j, uniformity). `python sdp3.py 3 8 1`. |
| `cut_lp.py` | Feasibility of a Wu–Zou cut-inequality certificate at odd n: feasible at n = 3, infeasible at n = 5, 7. |
| `group_lp.py` | The same with the three-copy theorem added for every grouping into three blocks: infeasible at n = 5. |

The infeasibility statements are floating-point linear programs; the obstruction itself is proved in
`OBSTRUCTIONS.md` (tightness at the odd product pairs).
