# Squeeze orbits: an every-n reduction and a conjecture

Ansh Mishra and Aryan Senthilkumar, 1 November 2026

These notes give a short reduction of the endpoint statement EP_n to a local property of the form on orbits of local
filters, and the numerical evidence for that property. **Nothing here proves EP_n for any new n.** The reduction is
proved; the property it needs (Conjecture ORB below) is supported by numerics only.

## Setting

For an operator X on (C^3)^⊗n and 0 ≤ t ≤ 1 put

    q_t(X) = Σ_J (−t)^{|J|} ‖Tr_J X‖²,    R_t(X) = q_t(X) / ‖X‖²,

so EP_n is R_{1/2}(C) ≥ 0 for every C of rank at most two. A *squeeze* is G = G_1 ⊗ … ⊗ G_n with each G_s positive
definite on C^3, acting by X ↦ G X G. It preserves the rank. With P_s = G_s² the form on the orbit factorises copy by
copy:

    q_t(G C G) = ⟨C| ⊗_s S_t(P_s) |C⟩,   ‖G C G‖² = ⟨C| ⊗_s P_s ⊗ P_s^T |C⟩,
    S_t(P) = P ⊗ P^T − t |vec P⟩⟨vec P|,

so on the orbit R_t is a ratio of two forms that are quadratic in each P_s separately. Unitaries U_s act trivially on
R_t (the form is U ⊗ Ū invariant at each copy), so squeezes describe the whole orbit of local congruences
X ↦ g X g^*.

Call a rank-two C *thick* if, at every copy, the local supports (supports of the one-copy reductions) of both its range
and its co-range are all of C^3. Otherwise C is *thin* at some copy.

## The reduction

**Conjecture ORB_n (at t = 1/2).** No thick operator C of rank two is a local minimum of G ↦ R_{1/2}(G C G) at
G = I, where G ranges over squeezes.

**Theorem.** EP_{n−1} and ORB_n imply EP_n.

*Proof.* The set of operators of rank at most two and norm one is compact, so R_{1/2} attains its minimum on it at some
X_0. If X_0 has rank one, R_{1/2}(X_0) ≥ 2^{−n}, by the rank-one bound q_n(xy^*) ≥ 2^{−n}‖xy^*‖² (Song–Chen). If X_0
has rank two and is thin at some copy, then R_{1/2}(X_0) ≥ 0. This is Lemma 6.1 of Tomar (t = 1/2), which needs
EP_{n−1}. If X_0 has rank two and is thick, then X_0 is a global, hence local, minimum of R_{1/2} on its own squeeze
orbit (X_0 itself is the point G = I), which ORB_n excludes. So min R_{1/2} ≥ 0. ∎

Since EP_1 holds, ORB_n for every n would give EP_n for every n, and with it NPT bound entanglement.

ORB_n is stronger than what the proof uses: a local minimum of R on the whole rank-two variety is in particular a local
minimum on the squeeze orbit through it. It also contains Tomar's first-zero criterion: at t < 1/2, a first failure of
EP_n would produce a thick global minimum of R_t (Tomar, Theorem 6.3), which is a squeeze-orbit minimum.

## Numerical evidence

`code/orbit_descent.py` minimises R_t(G C G) over the squeezes, starting from G = I, by L-BFGS in the coordinates
H_s with G_s = exp(H_s). A run "stops at the boundary" if the end point becomes thin (third local eigenvalue → 0) or
degenerates to rank one (σ_2/σ_1 → 0); an interior stop would be a squeeze-orbit minimum.

| copies n | t | starting operators | runs | interior stops |
|---|---|---|---|---|
| 2 | 1/2 | random rank two / random Hermitian rank two | 8 / 8 | 0 / 0 |
| 3 | 1/2 | random / random projectors / near a thin zero | 6 / 6 / 8 | 0 / 0 / 0 |
| 4 | 1/2 | random | 4 | 0 |
| 3, 2 | 0.4, 0.7 | random | 4, 4 | 0, 0 |

(Logs in `logs/`. We also made about 190 development runs, not all scripted here, and none had an interior stop:
- n = 2 to 6;
- t from 0.1 to 0.7;
- random, Hermitian, projector, positive, bosonic, near-zero and saddle-point starts;
- a three-copy operator with an unusually small contraction ratio.

Three of those runs, at n = 2, were first flagged as interior. They turned out to be rank-one degenerations ending at
exactly 2^{−n}.)

The near-zero starts (a thin zero P_V ⊗ |x⟩⟨y| plus a small thick perturbation) are already almost thin (thickness
5·10⁻⁸ to 3·10⁻⁴), so they test little. Several of these runs do not move at all because the gradient is below
tolerance.

Two things stand out.

1. **The end values are the conjectured minima.** For t ≤ 1/2 every run ends, up to optimiser residuals of at most
   4·10⁻⁵, at the product value (1−2t)(1−t)^{n−1} or at the rank-one value (1−t)^n. The product value is the
   conjectured minimum of R_t over rank-two operators, attained by P_V ⊗ |c⟩⟨c|^{⊗(n−1)}. Examples: 0 at t = 1/2;
   0.072 at n = 3, t = 0.4. Development runs also gave 0.12 at n = 2, t = 0.4 and 0.01664 at n = 4, t = 0.45.
2. **The descent still reaches the boundary past the endpoint.** For t > 1/2 single copies are distillable and R_t is
   negative on thin operators. At n = 2, t = 0.7 the four runs end at −0.12 (twice; the product value), at −0.4 and
   near −0.38. Here −0.4 = 1 − 2t is the one-copy value, the conjectured minimum past the endpoint. All four end on the
   boundary. So ORB looks like a property of these forms for every t, not a consequence of EP.

A direct search for a counterexample at n = 2 (minimise the squeeze gradient and the negative part of the squeeze
Hessian over thick C, with thickness ≥ 0.05) never achieved both: the squared gradient stayed ≥ 3·10⁻², while the
smallest Hessian eigenvalue stayed between −0.015 and −0.14.

Seven copies: with C kept in factored form, both random starts at n = 7, t = 1/2 descend to the boundary at R → 0 (our
script `orbit_big.py`, not included here).

**Critical points are saddles.** Thick critical points of R on squeeze orbits do exist. We found them at n = 2 by
minimising the squeeze gradient, also with a bias towards a chosen value of R:
- a family near R ≈ 0.129;
- points at R = 0.2;
- bulk points with R between 0.5 and 0.9.

Every one has negative squeeze-Hessian directions: 3 to 10 of the 16, with the most negative curvature between −0.23
and −5.6. (An earlier version of these notes also listed points near the thin boundary with R ≈ 7·10⁻⁴. They were not
critical points. A thickness penalty in the search held them in place, and their squeeze gradient was about 10⁻³.
Polishing them with the gradient alone sends them to the thin boundary.) More on the low critical points below.

Past the endpoint (t = 0.7) we also found thick, nearly critical points with negative value (R ≈ −0.03). They too are
saddles. So "thick critical points have positive value" fails for t > 1/2, and a proof of ORB cannot be first-order
uniformly in t: it has to see the negative curvature.

## What ORB is not

* **It is not a property of every invariant form.** For q(X) = ‖X‖² − c|Tr X|² (c = 0.1, two copies), every
  positive semidefinite thick rank-two C has interior squeeze-orbit minima: there G C G is proportional to a rank-two
  projector and |Tr X|² = 2‖X‖² (`code/trace_form_check.py`, `logs/trace_form_check.log`; thickness 0.08–0.22 at the
  minimum). At odd n the Werner form contains a term of exactly this kind, (−t)^n |Tr X|², but descents started at
  thick projectors do not stop there. So a proof of ORB must use the structure of the Werner form.
* **It is not the compression conjecture.** That conjecture asks for a single finite compression to a plane at one copy
  that lowers R. That fails at some points. ORB only asks for an infinitesimal squeeze, possibly at second order, at
  thick points.
* **It does not hold for every t.** For large t and odd n the term (−t)^n |Tr X|² takes over and the trap above
  appears. Descents started at thick projectors and positive operators then stop inside, at thick near-projectors with
  very negative R. The threshold t_c(n) is:

  | n | t_c(n) | evidence |
  |---|---|---|
  | 3 | ≈ 1.48 | no stop up to 1.45; 2 of 6 starts stop at 1.5, all at 1.6 |
  | 5 | ≈ 1.045 | none at 1.03; 3 of 4 at 1.06 |
  | 7 | ≈ 0.97 | none at 0.95; both at 1.0 |
  | 9 | between 0.95 and 1.0 | none at 0.95; all at 1.0 (R = −0.0263, thickness 0.64) |
  | 2 | — | no stop up to t = 2, consistent with the good sign of the trace term at even n |

  The endpoint t = 1/2 lies well inside the good range for every n tested. A heuristic shows why the trap lives near
  t = 1 rather than near 1/2. For a random thick plane W, Page-type purities give
  R_t(P_W) ≈ (1 − t/3)^n + 2((1 − 3t)/3)^n. For odd n this is negative only for t above a value that tends to 1, the
  point where X ↦ Tr(X) I − tX stops being positive. At t = 1/2 it equals (5/6)^n − 2(1/6)^n > 0. At n = 9, t = 1 the heuristic gives −(2/3)^9 = −0.0260;
  the trapped value found is −0.0263. The thresholds level off near 1 (≈ 1.48, 1.045, 0.97, 0.97 for n = 3, 5, 7, 9).
  Whether t_c(n) > 1/2 for every n is part of the conjecture.

## A formula for the second variation

Fix a copy s and an orthonormal basis (e_a) there, and write X_ab = ⟨a|_s C |b⟩_s for the contractions, q′ for the
form on the other n − 1 copies and r = R_t(C). If C is a critical point of R_t on its squeeze orbit, the second
variation along the rebalancing squeeze exp(ε(|a⟩⟨a| − |b⟩⟨b|)) at copy s is proportional to

    V(a, b) = q_r(D_−) − q_r(O),    q_r = q_t − r‖·‖²,
    D_− = |a⟩⟨a| ⊗ X_aa − |b⟩⟨b| ⊗ X_bb,   O = |a⟩⟨b| ⊗ X_ab + |b⟩⟨a| ⊗ X_ba.

So a squeeze-orbit minimum needs the signed diagonal part D_− of every two-dimensional compression at every copy to beat
its off-diagonal part O. Averaged over all planes this reduces to EP_{n−1}, so a proof has to choose the plane. The first
variation at copy s is ⟨v|Γ_s|v⟩ along exp(ε|v⟩⟨v|) for a traceless Hermitian 3 × 3 matrix

    ⟨v|Γ_s|v⟩ = ½[Σ_b q_r(⟨v|C|b⟩) + Σ_a q_r(⟨a|C|v⟩)] − t Re B(Tr_s C, ⟨v|C|v⟩),

where B is the bilinear form of q′. At critical points Γ_s = 0 for every s. A natural explicit direction, squeezing along
the traceless part of the local support operator of C, lowers R at some copy for most operators but not all (197 of 200
random two-copy operators).

## The second variation at a critical point

Fix a copy s, put P at copy s and I at the other copies, and write F_s(P) = N(P) − r D(P) for the numerator and
denominator of R on the orbit (r = R_t(C)). This is a quadratic form in P, and its formula makes sense for every
Hermitian P, not only for positive ones:

    F_s(P) = Σ_{abcd} P_ca P_bd B′_r(X_cd, X_ab) − t q′(Σ_ab P_ba X_ab),    B′_r = B′ − r⟨·,·⟩.

Squeezes at copy s near I are P = exp(2εH), and F_s(exp(2εH)) = 4ε² F_s(H) + O(ε³) at a critical point. A squeeze-orbit
minimum at I therefore needs F_s ≥ 0 on all of Herm(3), indefinite H included, with F_s(I, ·) = 0. Equivalently, the
kernel K(u, v) = ½[q′_r(X_uv) + q′_r(X_vu)] − t Re B′(X_uu, X_vv) on unit vectors of C^3 must be positive
semidefinite, with zero row averages. Three exact facts about F_s (checked numerically to 10⁻¹⁰ or better by
`code/squeeze_hessian.py`):

1. **Trace.** For every C, with r = R_t(C), the trace over an orthonormal basis of Herm(3) is
   Tr F_s = q′_r(Tr_s C) − t Σ_ab q′(X_ab). When r ≤ 0, EP_{n−1} makes it nonnegative:
   Σ_ab q′(X_ab) ≥ 0, and then q_n(C) = r‖C‖² forces q′(Tr_s C) ≥ 0. So averaging over squeeze directions at one copy cannot find the descent, and the direction has to
   be chosen.
2. **Rank-one squeezes are monotone.** At a critical point, F_s(I − E) = F_s(E) for E = |e⟩⟨e|. Along the path
   G_p = (I − E) + √p E from C (p = 1) to its thin compression (p = 0), q_r(G_p C G_p) = F_s(E)(1 − p)² exactly,
   with F_s(E) = (1 − t) q′(X_ee) − r‖X_ee‖². The path descends iff R_{n−1}(X_ee) < r/(1 − t). At a hypothetical negative
   minimum (r < 0) this never happens, by EP_{n−1}: squeezing towards a thin operator one direction at a time only
   goes uphill there.
3. **Qubit planes.** For a plane K ⊂ C^3 with Pauli matrices σ_μ^K, Σ_μ S_t(σ_μ^K) = −(1 + 2t) S_{t′}(I_K) as operators
   on the qubit, with t′ = 1/2 + 3/(2(1 + 2t)) (t′ = 5/4 at t = 1/2). Hence

       Σ_μ F_s(σ_μ^K) = −(1 + 2t) q_{(t′ at s)}(C_K) − r(2‖Tr_K C‖² − ‖C_K‖²),

   where C_K is the compression of C to K at copy s, Tr_K C = Σ_{k∈K} ⟨k|C|k⟩, and q_{(t′ at s)} uses the parameter t′
   at copy s and t elsewhere. At t = 1/2 and r = 0 this is (3/2) q′(Tr_K C) − 2 q_n(C_K).

**What the low critical points show.** `code/squeeze_hessian.py n t target starts seed` searches for a thick
critical point with a bias towards a target value of R and a thickness floor of 0.02. It then removes both and polishes
with the squeeze gradient alone. A run counts as genuine if the squeeze gradient ends below 10⁻⁸ with thickness above
10⁻³. Results at t = 1/2 (logs in `logs/squeeze_hessian_n*_t*_R*.log`):

| n | targets that gave genuine critical points | targets where every run slid to the thin boundary |
|---|---|---|
| 2 | 0.15 (4/4), 0.2 (6/6), 0.4 (4/4) | 0.05 (4), 0.1 (4) |
| 3 | 0.08 (R = 0.080, 0.106), 0.15 (3/3), 0.3 (3/3) | 0.03 (3), 0.04 (3), 0.05 (3) |
| 4 | 0.05 (1/1), 0.1 (2/2) | 0.03 (2) |
| 5 | 0.05 (first run) | — |

So the lowest thick critical value lies between 0.10 and 0.13 at n = 2, between 0.05 and 0.08 at n = 3, and between
0.03 and 0.05 at n = 4. That is of order 2^{−n}. This is what a search found, not a bound.
(At t = 0.4, target 0.2 at n = 2 slid to the boundary in 4/4 runs.)

Every genuine critical point is a saddle, with most negative squeeze curvature between −0.13 and −2.5 (in the
coordinates of the script). At the low ones
the descent comes mostly from coupling between copies. Take the most negative direction of the full squeeze Hessian.
Along it, the single-copy blocks contribute a positive or small amount, and the cross-copy terms give most of the
curvature. Examples (single-copy part / cross-copy part):
- n = 2, R = 0.15: +0.50 / −1.32;
- n = 3, R = 0.08: +0.25 / −0.47;
- n = 4, R = 0.05: +0.003 / −0.134;
- n = 5, R = 0.05: −0.026 / −0.249.

Every single-copy block still has a negative eigenvalue, but a much weaker one. At n = 3, R = 0.08 they are between
−0.006 and −0.009, against −0.21 for the full Hessian.

We tested three explicit families of directions:
- (a) rank-one squeezes;
- (b) the trace over a qubit plane (formula 3 above);
- (c) the *top-two rebalance*. Here H_s = |b_s⟩⟨b_s| − |a_s⟩⟨a_s| at every copy s, where a_s and b_s are the largest and
  second-largest eigenvectors of the local support operator Tr_other(CC^* + C^*C), and the n weights are optimised
  jointly.

None of them works at every point. (a) and (b) see nothing at several low points. (c) is negative at all low points
(n = 2, R ≤ 0.2; n = 3, R = 0.08; n = 4 and 5). It is positive at 2 of the 4 n = 2 points at R = 0.4 and at 1 of the 3
n = 3 points at R = 0.15. So far every point is caught by at least one of the three. A proof of ORB that looks at one copy at a time
would have to find the weak single-copy directions. The cross-copy terms are where the curvature actually is.

## One copy at a time is not enough

Two natural strengthenings of ORB look only at one copy. Both are false at n = 3, t = 1/2. The operators are in `data/`
and `code/one_copy_counterexamples.py` checks them (output in `logs/one_copy_counterexamples.log`).

* **Single-copy version.** "If C is thick at copy s, then R over the squeeze orbit of copy s alone (other copies fixed)
  has no interior minimum." Since F_s is quadratic, an interior local minimum on that orbit is a global minimum of R over
  its closure. If this version held at one copy, EP_{n−1} would give EP_n at once: R(C) ≥ min over singular P ≥ 0. Random
  starts support it: single-copy descents reached the thin boundary in 26/26 runs at n = 2, 33/33 at n = 3, 8/8 at
  n = 4, and 24/24 at n = 2 for t = 1/3, 0.4 and 0.45. But an adversarial search finds a counterexample,
  `data/single_copy_min_n3.npy`, with R = 0.2968 and thickness 0.18–0.44 at every copy. Its copy-0 squeeze gradient is
  9·10⁻¹⁴, and F_0 has eigenvalues 0 (the direction I) and 5.6·10⁻⁴ to 0.022. So C is a strict minimum of R over the
  whole copy-0 orbit. (It is not critical at the other copies.)
* **Copywise version.** "No thick C is critical at every copy with every F_s positive semidefinite." Also false:
  `data/copywise_min_n3.npy` has R = 0.2623, thickness 0.22–0.42, squeeze gradient ≤ 1.1·10⁻¹³ at every copy, and each
  F_s has eigenvalues 0 (the direction I) and at least 2.9·10⁻⁴. No squeeze of any single copy lowers R, to second order
  or globally on that copy's orbit.

**ORB survives at both points.** The joint squeeze Hessian of the copywise example has 8 negative eigenvalues out of 24,
the smallest −0.57. Every pair of copies is already enough: the 16 × 16 blocks for the pairs (0,1), (0,2) and (1,2)
have smallest eigenvalues −0.33, −0.41 and −0.31. The same happens at four copies. `data/copywise_min_n4.npy` has
R = 0.2277, thickness 0.16–0.36, squeeze gradient ≤ 10⁻¹³ at every copy, and each F_s ≥ 1.1·10⁻⁴ on traceless
directions. Its joint Hessian has 12 negative eigenvalues out of 32 (smallest −0.39), and all six pair blocks are
negative (−0.11 to −0.22). So ORB is not a statement about single copies. A proof has to use the
coupling between copies. These points have R > 0. Restricted to negative values, the single-copy and copywise
statements are equivalent to EP_n itself, given EP_{n−1}. A thick negative global minimum would admit no descent at
all, and a non-minimum descends anyway. So the restricted versions do not give an easier target.

## The joint second variation

Let C be squeeze-critical with r = R_t(C), and take a 1-local direction H = Σ_s H_s (H_s acting at copy s). For
C(ε) = e^{εH} C e^{εH}, R(C(ε)) = r + ε² 𝓗(H)/‖C‖² + O(ε³), where

    𝓗(H) = 2 q_r({H, C}) − 2 Σ_s [ q_r(H_s C) + q_r(C H_s) ]
           + 8t² Σ_{s<s′} Re[ B″(Tr_{ss′} C, Tr_{ss′}(H_s H_{s′} C)) − B″(Tr_{ss′}(H_s C), Tr_{ss′}(H_{s′} C)) ].

Here {H, C} = HC + CH and B″ is the bilinear form on the other n − 2 copies (for n = 2, the product of the two
traces). `code/joint_second_variation.py` checks this against automatic differentiation to 10⁻⁸ or better. Expanding the
first term, the "self" parts cancel. So 𝓗(H) = 4 Σ_{k<l} Re B_r(V_k, V_l) + (the trace terms), summed over the 2n
rank-two pieces V ∈ {H_s C, C H_s}. The second variation is the total pairwise correlation of the left and right
multiplications. For one copy this gives F_s(H) = Re B_r(HC, CH).

At a thick global minimum with r < 0, every q_r(H_s C) and q_r(C H_s) is ≥ 0, since these operators have rank two. So
the middle sum pushes 𝓗 down, and everything depends on q_r of the rank-four tangent vector {H, C}.

At the copywise example the numbers are as follows. Along the best direction of each pair (s, s′), the single-copy
terms 4F_s(H_s) are only +0.005 to +0.013, because the anticommutator energy q_r({H_s, C}) ≈ 0.09–0.12 almost cancels
the commutator energy Σ_ab q′_r([H_s, C]_ab) ≈ 0.08–0.11. Meanwhile the tangent vectors {H_s, C} and {H_{s′}, C} at the
two copies have correlation Re B_r / √(q_r q_r) between −0.41 and −0.55, which gives a cross term of −0.15 to −0.23. The
trace terms are only ±0.02.

We also searched adversarially at n = 3 for a thick critical point where every pair block is positive semidefinite.
The search could push the pair blocks to −10⁻⁵, but only away from criticality (gradient ≥ 3·10⁻³) and at the
thickness floor 0.05 of the search. There the full Hessian was still −0.018. We did not reach a thick critical point
where all pairs are flat.

## Status

Reduction proved (above). ORB_n is a conjecture, tested numerically for n ≤ 7. A proof of ORB_n for all n would prove
NPT bound entanglement; we do not have one.

## Files

* `code/orbit_descent.py`: squeeze-orbit descent (`python code/orbit_descent.py n t starts seed ensemble`).
* `code/trace_form_check.py`: the control above.
* `code/squeeze_hessian.py`: thick squeeze-critical points, the split of their Hessian, the three direction tests and
  checks of the identities above.
* `code/one_copy_counterexamples.py` with `data/single_copy_min_n3.npy`, `data/copywise_min_n3.npy`,
  `data/copywise_min_n4.npy`: the operators above.
* `code/joint_second_variation.py` (with `data/critical_n2_R0.2.npy`): check of the joint second-variation formula.
* `logs/`: the runs in the tables and the control.

## References

* T. S. Tomar, *Equality cases and every-copy results for the Werner endpoint*, Zenodo (2026): Lemma 6.1 (thin copies)
  and Theorem 6.3 (first zeros).
* L. Song and L. Chen, arXiv:2607.23416 (rank-one bound).
