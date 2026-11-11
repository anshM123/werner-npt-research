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
| 5 | 0.05 (2/2) | — |

So the lowest thick critical value lies between 0.10 and 0.13 at n = 2, between 0.05 and 0.08 at n = 3, and between
0.03 and 0.05 at n = 4. That is of order 2^{−n}. This is what a search found, not a bound.
(At t = 0.4, target 0.2 at n = 2 slid to the boundary in 4/4 runs.)

On the whole rank-two variety, not only on squeeze orbits, the lowest thick critical points we found have closed-form
values. All of them are saddles (`code/exact_critical_points.py`, `logs/exact_critical_points.log`, data in `data/`).

| n | R | negative Hessian directions | structure |
|---|---|---|---|
| 2 | (9 − √41)/20 ≈ 0.12984 | 6 (smallest −0.58) | in the two-sector class of `../two-sector-class/` |
| 3 | (3 − √5)/8 ≈ 0.09549 | 30 (smallest −1/3) | local spectra (1/3, 1/3, 4/3) at every copy |
| 3 | (67 − √2329)/216 ≈ 0.08676 | 15 (smallest −0.28) | |C| in local eigenbases invariant under copy permutations |

The endpoint inequality is proved for every n on the two-sector class, so no operator of the two-copy shape gives a
counterexample at any number of copies. The three-copy points are not in that class. Adversarial searches for a thick
local minimum on the variety, at two and three copies, either stall at the thickness floor with many negative Hessian
directions or slide to thin global minima with R = 0.

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

## Natural families of directions are not enough

The descent at thick critical points is not confined to the obvious families of squeezes. Two of them fail at three
copies:
- the **torus**: traceless squeezes diagonal in the eigenbasis of the local support operator μ_s = Tr_other(CC^* + C^*C)
  at every copy (2n directions, an abelian family);
- the **weak planes**: squeezes in su(K_s), where K_s is the plane of the two smaller eigenvectors of μ_s
  (3n directions).

The two counterexamples:
- `data/torus_pd_n3.npy` is a thick squeeze-critical point: gradient 3·10⁻¹⁰, R = 0.2573, thickness 0.20. On the torus
  the Hessian is positive definite, with smallest eigenvalue +0.0029. The weak planes give −0.048 and the full squeeze
  Hessian −0.084.
- `data/suk_pd_n3.npy` is the other way round: gradient 9·10⁻¹¹, R = 0.3177, thickness 0.18. On the weak planes the
  Hessian is positive definite (+0.0025). The torus gives −0.063 and the full squeeze Hessian −0.151.

Both were found from `copywise_min_n3.npy` by minimising the squared squeeze gradient with a penalty on the restricted
Hessian, then polishing the gradient alone. Both are ordinary saddles. Every pair of copies descends there
(−0.046 to −0.129), and each point is caught by the families (a) and (b) of the previous section; the second is also
caught by (c).

At two copies the same searches stall far from criticality, with squared gradient about 0.025, so there these families
may be enough. At three copies a search of the same kind for a thick critical point with every pair block positive
semidefinite also stalls, at squared gradient ≥ 8·10⁻⁴. This agrees with the observation above that pairs of copies were
always enough. `code/restricted_families.py` checks both points by finite differences
(`logs/restricted_families.log`).

Two further observations (numerical only):
- **Pairs.** Pushed further (3000 steps, then L-BFGS on the penalised objective), the three-copy pair search gets stuck
  at squared gradient 6·10⁻⁶ > 0 with thickness 0.07. A four-copy run from `copywise_min_n4.npy` stalls at 1.1·10⁻³. So
  far some pair of copies always carries the descent.
- **Low values.** Both counterexamples above have fairly large values. If the torus search is also required to keep
  R ≤ 0.1 (the lowest thick critical values at three copies are about 0.074–0.087), every run pushes R up to the cap and
  gets stuck at squared gradient about 5·10⁻⁶. So we found no torus-flat thick critical point at low values. At the
  lowest points the descent is a rebalance inside the weak plane, with the same sign at every copy. Only low and
  nonpositive values matter for the reduction. There, however, the torus statement is equivalent to EP_n itself, given
  EP_{n−1}, so this locates the mechanism but does not prove it.

## Two-sided filters

The reduction also works for the larger group of local filters acting independently on range and co-range,
C ↦ (⊗_s g_s) C (⊗_s h_s)^T with g_s, h_s ∈ GL(3). It contains the squeezes and also the "similarities"
C ↦ G C G^{-1}. Its orbit closures again add only thin or rank-one operators. So the statement "no thick rank-two
operator is a local minimum of R on its two-sided orbit" (call it ORB2) also gives EP_{n−1} ⇒ EP_n. ORB2 is weaker
than ORB, so it is the easier target.

This group acts holomorphically. So at a local minimum, f = q − r‖·‖² has a positive semidefinite Levi form, which
here is just q_r restricted to the complex tangent space {XC + CY}, with X and Y 1-local. That is a statement about one
fixed linear space, with none of the second-order terms of the group action.

At a critical point of R on the whole rank-two variety, put L = Φ_r(C) = Φ(C) − rC, the residual of the Werner map
(L vanishes on the co-range of C, and its range is orthogonal to the range of C). For 1-local Hermitian K, the exact
second variations of q_r are:

    squeeze    e^{εK} C e^{εK}   :  ε² [ 2 q_r(KC) + 2 q_r(CK) + 4 b(K) ],
    similarity e^{εK} C e^{−εK}  :  ε² [ −4 b(K) ],          b(K) = Re⟨L, K C K⟩.

For K at a single copy, in its eigenbasis (this needs only criticality for squeezes at that copy),

    b(K_s) = −½ Σ_ab (k_a − k_b)² q′_r(X_ab).

`code/two_sided_identities.py` checks all three at a critical point (finite differences, and exactly for the last).
At a thick global minimum with r < 0, every q_r(KC), q_r(CK) and q′_r(X_ab) is ≥ 0, since these operators have rank
two. So every single-copy similarity and every one-sided term points the wrong way. Descent there has to come from
the cross-copy part of b,

    Re⟨L, K_s C K_{s′}⟩ = Σ_i Re⟨Π_A^⊥ K_s A_i, L Π_B^⊥ K_{s′} B_i⟩     (C = Σ_i A_i B_i^*),

which pairs the leakage of the range at copy s with the leakage of the co-range at copy s′ through L. Numerically
neither family works alone. An adversarial search at n = 2 finds thick critical points with no similarity descent
(largest b about +5·10⁻⁷). There the squeeze Hessian is −2.9 to −4.0 and the two-sided Hessian about −0.8 to −1.0.

**One copy is not enough here either.** Consider two-sided filters at a single copy,
C ↦ (g ⊗ 1)C(h ⊗ 1)^T with g, h ∈ GL(3). They reach the thin boundary from every random start we tried (16/16 at
n = 2 and 3), including from both squeeze examples above. But an adversarial search finds a strict local minimum at
n = 3: `data/two_sided_single_copy_min_n3.npy`, with R = 0.3111 and thickness 0.26–0.33 at every copy. The gradient over
(X, Y) at copy 0 is 1.5·10⁻¹³. On the 24 directions where R is not constant, the Hessian has smallest eigenvalue
0.010. The excluded directions are the unitary conjugations and scalings, whose span has rank 12; R is constant along
them to 10⁻¹⁶. Along 300 random directions R goes up, and perturbed descents return to the same value.
`code/two_sided_single_copy_min.py` checks all of this. (A first version of this check removed one direction too
many, because a QR factorisation of the rank-deficient set of invariant directions was used; a rank-revealing
projection is needed.) At copies 1 and 2 the gradient is 0.10–0.12, so the full two-sided orbit still descends. Both
the squeeze version and the two-sided version of the orbit property are therefore genuinely multi-copy.

**The Levi form alone is not enough.** At a point C critical for the two-sided filters, C itself lies in V_C with
q_r(C) = 0, so "q_r is not positive semidefinite on V_C" would be a pure eigenvalue statement, and with it ORB2 would
follow. It is false at two copies. `data/levi_positive_n2.npy` is a thick rank-two operator with
- R = 0.143505, thickness 0.17 and σ₂/σ₁ = 0.44;
- two-sided gradient 1.4·10⁻¹⁵;
- q_r positive on the 28-dimensional complement of C in V_C, with smallest eigenvalue 0.0037.

It is still a saddle. Along e^{εX} C e^{εY} the second variation of q_r is

    q_r(XC + CY) + Re⟨L, X²C + 2XCY + CY²⟩,    L = Φ_r(C),

and replacing (X, Y) by (iX, iY) flips the sign of the second term. So at a point where the Levi form is positive, any
descent comes from the curvature term ⟨L, ·⟩. Here the real second variation has 8 negative directions, and the
smallest is −3.64 (per unit tangent norm). Along the best of them the Levi part is +0.85 and the curvature part is
−4.49. `code/levi_positive_point.py` checks these numbers (`logs/levi_positive_point.log`). So a proof of ORB2 has to
control the pairing of L with the second-order motion, not only q_r on the tangent space.

## The co-range optimised first (added 4 November 2026)

ORB moves the co-range together with the range. Optimising the co-range first gives a weaker property that is still
enough, and it removes the one-copy traps found above. For a two-dimensional range plane U ⊂ (C^3)^{⊗n} put

    m(U) = min over all co-ranges of R_{1/2}(C) = λ_min [Ψ^{⊗n}(w_r w_r'^*)]_{r,r'},

with (w_1, w_2) an orthonormal basis of U and Ψ(Y) = Tr(Y) I − Y/2 on each copy. Local unitaries leave m unchanged, so
on the orbit {(⊗_s G_s)U : G_s ∈ GL(3)} it is a function of the positive parts G_s^*G_s.

**Conjecture ORB_m (at t = 1/2).** No range plane U that is thick at every copy (range local supports all of C^3) is a
local minimum of m on its orbit at a value m(U) < 2^{−n}.

**Theorem.** EP_{n−1} and ORB_m for n copies imply EP_n.

*Proof.* Rank-one operators satisfy q_n ≥ 2^{−n}‖·‖² (Song–Chen). The Grassmannian of planes is compact and m is
continuous, so m attains its minimum at some U*. If U* is thin at some copy, every C with range U* is thin there, and
Lemma 6.1 of Tomar with EP_{n−1} gives q_n(C) ≥ 0. Otherwise U* is thick. It is a minimum of m on its own orbit, so
ORB_m forces m(U*) ≥ 2^{−n} > 0. ∎

The threshold 2^{−n} is there because m can have flat pieces. On single-excitation planes, m equals 2^{1−n} on an
open set, where the best co-range has rank one (`povm-minimum-principle/`). Values coming from a rank-one co-range
are at least 2^{−n} by the rank-one bound, so the threshold excludes them. A local minimum of m on the orbit of U, at a thick U whose optimal co-range is thick,
would give a local minimum of R on the squeeze orbit of the corresponding C, since R(GCG) ≥ m(GU) ≥ m(U) = R(C). So
ORB implies ORB_m at such points. ORB_m is the weaker statement.

On single-excitation and tangent planes the orbit of a thick plane is the interior of a set of qubit POVMs, and ORB_m
becomes the minimum principle of `povm-minimum-principle/`.

Numerics (`code/orbm_descent.py`, `code/orbm_hard_starts.py`, `code/orbm_perturbed.py`):
- **Random planes.** L-BFGS on the filter parameters, starting at G = I. Every run ends at a thin plane (range
  local-support ratio 10⁻⁷ or smaller) with m → 0: 10 of 10 at n = 2, 12 of 12 at n = 3, and 3 of 3 at n = 4
  (`logs/orbm_descent.log`).
- **Hard starts.** We started from the range planes of every stored operator above, and of the tangent-class and
  shared-shift operators of `tangent-class/` (17 operators at n = 2, 3, 4). From all of them m descends to 0 at a thin
  plane (`logs/orbm_hard_starts.log`). That includes the copywise and single-copy squeeze minima, which are traps for
  ORB-type one-copy arguments but not for m.
  - The three exact critical points on the rank-two variety do not move, since the gradient vanishes there. After a
    perturbation of size 0.02 in the filter parameters, m drops below the critical value at once and descends to 0
    (`logs/orbm_perturbed.log`). So they are saddles for m too.

These are descents, not an adversarial search for local minima of m; see `povm-minimum-principle/` for such a search
on single-excitation planes.

## Weak planes and the residual sandwich (added 10 November 2026)

We looked for the descent direction itself at every thick squeeze-critical point we could produce: the ten points saved
above, 35 new ones at t = 1/2 (n = 2, 3, 4, 5; R from 0.05 to 0.6) and 14 at t = 0.7 with negative values (n = 2, 3;
R from −0.03 to −0.15). That is 59 points in all: 28 at n = 2, 24 at n = 3, 6 at n = 4 and one at n = 5. Six of the new
ones are in `data/` (`critical_n4_R0.05.npy`, `critical_n5_R0.05.npy`, `critical_t0.7_n2_neg_a.npy`,
`critical_t0.7_n2_neg_b.npy`, `critical_t0.7_n3_neg_a.npy`, `critical_t0.7_n3_neg_b.npy`), and
`code/weak_plane_descent.py` reproduces the numbers below on the saved points (`logs/weak_plane_descent.log`).

**The cross term, exactly.** Let C be squeeze-critical with ‖C‖ = 1, r = R_t(C) and Y = Φ_t(C) − rC, where
Φ_t = (id − tD)^{⊗n}. For H = A + B with A at copy s and B at copy s′ (G = e^{εH}), the cross entry of the squeeze
Hessian is X(A, B)/2 with

    X(A, B) = 4 Re B_r({A, C}, {B, C}) + 4 Re⟨Y, {AB, C}⟩ + 4 Re⟨Y, ACB + BCA⟩,    B_r(U, V) = ⟨U, Φ_t(V) − rV⟩,

and the diagonal entry is S(A) = 2 q_r({A, C}) + 2 Re⟨Y, A²C + 2ACA + CA²⟩. The script checks the cross entry against
finite differences of R along e^{εH} C e^{εH} (−0.129519 against −0.129519 at t = 1/2, +0.843847 against +0.843847 at
t = 0.7). The middle term involves the two-copy operator AB. At a critical point of the whole rank-two variety it
vanishes, because Y is then orthogonal to every XC + CX′. The last term, the *sandwich*, pairs Y with C multiplied on
the left at one copy and on the right at the other. It is the cross-copy part of b in "Two-sided filters" above.

At a critical point of the whole variety the two-sided second variation also simplifies. There Φ_t(C) = rC + L, with
L mapping the orthogonal complement of the co-range into the orthogonal complement of the range. So ⟨L, X²C⟩ =
⟨L, CY²⟩ = 0, and along e^{εX} C e^{εY} the second variation of q_r is exactly q_r(XC + CY) + 2 Re⟨L, XCY⟩. At a global
minimum, XC and CY have rank two, so q_r(XC) ≥ 0 and q_r(CY) ≥ 0. Rescaling X and Y and choosing a phase shows that
the absence of two-sided descent is equivalent to

    |B_r(XC, CY) + ⟨L, XCY⟩|² ≤ q_r(XC) q_r(CY)    for all 1-local X (acting on the left) and Y (on the right).

A thick global minimum with r < 0 is therefore excluded exactly when some pair (X, Y) violates this inequality.

**Where the descent lives.** Let K_s be the plane of the two smallest eigenvectors of the local support operator
μ_s = Tr_other(CC^* + C^*C) (the *weak plane* at copy s), and σ_μ (μ = 1, 2, 3) its Pauli matrices in that
eigenbasis. Pair two copies,

    H^μ = σ_μ^{(s)} + (Oσ)_μ^{(s′)},    O ∈ diag(±1, ±1, ±1) · {1, (σ_1 ↔ σ_2)}    (16 identifications),

and let PT = Σ_μ (second variation along H^μ) = Σ_μ [S(σ_μ^{(s)}) + S((Oσ)_μ^{(s′)}) + X(σ_μ, (Oσ)_μ)].
- **The paired weak-plane squeeze descends at 58 of the 59 points**, for some pair of copies and some O, including all
  14 points with negative value. The exception is `data/suk_pd_n3.npy`, which was constructed so that every
  weak-plane squeeze block is positive semidefinite (R = 0.318). Choosing O instead by the polar decomposition of the
  sandwich matrix N_{μν} = 4 Re⟨Y, σ_μ C σ′_ν + σ′_ν C σ_μ⟩ also works at 57 of them.
- **The sandwich drives it.** Splitting PT (polar identification) into the single-copy part, the B_r part, the
  two-copy Y part and the sandwich, the sandwich is negative at every point and is the largest negative contribution at most of them. The
  single-copy part is negative at the exact critical points of the whole variety but positive at several points near
  the thin boundary.
- **The weak planes matter at the hard points.** With random planes in place of the weak ones (20 choices per point),
  the same construction descends for at most 15% of the choices at the exact critical points
  `variety_critical_n2_low`, `variety_critical_n3_golden`, `variety_critical_n3_low` and at the points with negative
  value. At some other points random planes do about as well (up to 90%).
- **Averaged and simpler versions fail.**
  - No linear combination of 14 natural local features descends at every point. The features were range and co-range
    marginals, their products and inverses, partial traces and residual marginals, with the same coefficients at
    every copy.
  - Summing over all nine Weyl operators at two copies gives a closed form for the single-copy factors,
    Σ_W q_r(W^{(s)} C) = (15/2) Σ_ab q_{n−1}(X_ab) − 9r at t = 1/2. But the resulting criterion fails at every low
    point tested, with ratio 0.27 to 0.78 against the threshold 1.
  - Multiplying C on either side by a swap of two copies, plain or twisted between the weak planes, never lowers R at a
    point with R ≤ 0.15.
  - "Some direction inside the block-diagonal family su(K_s) ⊕ ℝ(P_{K_s} − 2P_{d_s}) descends" holds at every point,
    but so does the same statement for random four-dimensional families. Only the explicit paired direction says
    something.

So the descent at thick critical points is concrete. Squeeze the weak planes of two copies along matched Pauli
directions, and the residual L, paired with the range leakage at one copy and the co-range leakage at the other, does
the rest. A proof that this paired squeeze (or the polar version) always descends at thick critical points with
r ≤ 0 would prove ORB_n at the points that matter, and with it EP_n for every n. We have no such proof. The weak planes
are spectral subspaces of the marginals, so any argument has to control them at an unknown point.

## Product-qubit compressions (added 10 November 2026)

At copy s let Q_s ⊂ C^3 be a plane with unit normal ℓ_s and put P = P_{Q_1} ⊗ ⋯ ⊗ P_{Q_n}. On an operator X whose
rows and columns lie in Q_1 ⊗ ⋯ ⊗ Q_n the endpoint form sees, at each copy, only the compression of I − ½|Ω⟩⟨Ω|
(Ω = Σ_a |aa⟩) to Q_s ⊗ Q̄_s. That compression is I − ½|Ω_Q⟩⟨Ω_Q|, with eigenvalues 1, 1, 1, 0, so q_n(1/2, X) ≥ 0
for every such X, of any rank. The same holds when the rows lie in a plane A_s and the columns in a plane B_s: the
compression is then I − ½|vec(P_A P_B)⟩⟨vec(P_A P_B)| and ‖P_A P_B‖_F² ≤ 2. This suggests a strengthening of EP_n
that needs no induction,

    (CC_n)   every rank-two C has planes Q_1, …, Q_n with PCP ≠ 0 and R(PCP) ≤ R(C).

Here PCP is the limit of G_λ C G_λ with G_λ = ⊗_s (P_{Q_s} + λ P_{ℓ_s}) as λ → 0, so CC_n asks for a point of the
squeeze-orbit closure of C that lives on a product of qubits and is no worse than C.

- **One copy: true.** Take Q to be the range of C. Then PCP = C P_Q, and the rest C − C P_Q = C P_ℓ has its rows in Q
  and its column along ℓ. On those entries the form is the identity and does not couple to C P_Q, and ⟨ℓ|C = 0. So
  q(C) = q(C P_Q) + ‖C P_ℓ‖² and ‖C‖² = ‖C P_Q‖² + ‖C P_ℓ‖², and R(C) ≥ R(C P_Q) because R ≤ 1. (C P_Q ≠ 0 because
  two planes in C^3 meet.)
- **Three copies: false numerically.** `data/compression_counterexample_n3.npy` has R = 0.0945311529, while the
  smallest value of R over product-qubit compressions that we find is 0.0977396233. That value is reached repeatedly
  by 1581 L-BFGS runs, including all 81 starts at eigenvectors of the one-copy marginals of CC^*, C^*C and their sum,
  and again by 400 BFGS runs in an angle parametrisation with q computed from partial traces
  (`logs/qubit_compression_verify_n3.log`, `logs/qubit_compression_check_n3.log`). These are local searches in
  floating point; we have no certified global lower bound over the twelve-dimensional set of compressions. The point was found by maximising min_Q R(PCP) − R(C) from
  `data/variety_critical_n3_low.npy`, the critical point of the rank-two variety with R = (67 − √2329)/216, where the
  margin is already down to −0.013 (`logs/qubit_compression_points.log`).
- **Searches from random starts miss it.** Gradient ascents of the same gap reach −0.0013 at n = 2, −0.025 at n = 3
  and −0.20 at n = 4, and near the zeros at n = 2 the best compression stays 0.013 to 0.195 times R(C) as C → zero
  (`logs/qubit_compression_adversary*.log`, `logs/qubit_compression_nearzero_n2.log`). Every saved critical point,
  n = 2 to 5, has a negative gap.
- **One copy at a time fails already at two copies.** Ascents find rank-two C for which every compression of copy 0
  raises R (by 0.021), every compression of the rows of copy 0 raises R (0.096), every compression of either single copy
  raises R (0.006), and, with copy 0 already a qubit, every compression of copy 1 raises R (0.030)
  (`logs/qubit_compression_onecopy_*.log`).
- **Independent planes for rows and columns** survive at the counterexample, where the best such compression has
  R = 0.0155. But this version is weak: at `data/variety_critical_n3_low.npy` the best one already has
  R = 3.6·10⁻⁵ (`logs/qubit_compression_independent.log`), so its minimum sits near zero and it says little beyond
  EP_n. Arbitrary rank-two filters are weaker still: along the squeeze orbit of a random rank-two C at n = 2, descent
  runs to operators thin at both copies, with R below 10⁻⁸ in five of eight trials
  (`logs/qubit_compression_orbit_n2.log`; of the others, two stalled at points thin at both copies and one overflowed).

So the same-plane compression is not a route to EP_n, and the independent-plane version is not a lever.

## Status

Reduction proved (above). ORB_n is a conjecture, tested numerically for n ≤ 7. The weaker ORB_m, with the co-range
optimised first, also suffices; it is a conjecture, tested by descents for n ≤ 4. A proof of either for all n would
prove NPT bound entanglement; we do not have one.
The product-qubit compression inequality CC_n (section above) holds at one copy and fails numerically at three
(local searches, not certified).

## Files

* `code/orbit_descent.py`: squeeze-orbit descent (`python code/orbit_descent.py n t starts seed ensemble`).
* `code/trace_form_check.py`: the control above.
* `code/squeeze_hessian.py`: thick squeeze-critical points, the split of their Hessian, the three direction tests and
  checks of the identities above.
* `code/one_copy_counterexamples.py` with `data/single_copy_min_n3.npy`, `data/copywise_min_n3.npy`,
  `data/copywise_min_n4.npy`: the operators above.
* `code/joint_second_variation.py` (with `data/critical_n2_R0.2.npy`): check of the joint second-variation formula.
* `code/two_sided_identities.py` (with `data/variety_critical_n2.npy`): the squeeze and similarity second variations at a
  critical point on the rank-two variety.
* `code/two_sided_single_copy_min.py` (with `data/two_sided_single_copy_min_n3.npy`): a strict local minimum for two-sided
  filters at one copy (n = 3).
* `code/levi_positive_point.py` (with `data/levi_positive_n2.npy`): a two-copy two-sided critical point with positive
  Levi form that still descends through the curvature term.
* `code/exact_critical_points.py` (with `data/variety_critical_n2_low.npy`, `data/variety_critical_n3_golden.npy`,
  `data/variety_critical_n3_low.npy`): the exact low thick critical points on the rank-two variety.
* `code/orbm_descent.py`, `code/orbm_hard_starts.py`, `code/orbm_perturbed.py`: the co-range-optimised orbit property
  ORB_m (logs `logs/orbm_*.log`).
* `code/weak_plane_descent.py` (with the saved critical points, including `data/critical_n4_R0.05.npy`,
  `data/critical_n5_R0.05.npy` and `data/critical_t0.7_n*_neg_*.npy`): the cross second variation and the paired
  weak-plane squeeze of the section above (`logs/weak_plane_descent.log`).
* `code/qubit_compression.py`, `code/qubit_compression_check.py` (with `data/compression_counterexample_n3.npy`): the
  product-qubit compressions of the section above (`logs/qubit_compression_*.log`).
* `code/restricted_families.py` (with `data/torus_pd_n3.npy`, `data/suk_pd_n3.npy`): three-copy critical points at which
  the torus, resp. the weak-plane family, has no descent.
* `logs/`: the runs in the tables and the control.

## References

* T. S. Tomar, *Equality cases and every-copy results for the Werner endpoint*, Zenodo (2026): Lemma 6.1 (thin copies)
  and Theorem 6.3 (first zeros).
* L. Song and L. Chen, arXiv:2607.23416 (rank-one bound).
