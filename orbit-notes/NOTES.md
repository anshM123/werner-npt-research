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
minimising the squeeze gradient, also with a bias towards small R:
- a family near R ≈ 0.129;
- points close to the thin boundary with R ≈ 7·10⁻⁴;
- bulk points with R between 0.5 and 0.8.

Every one has negative squeeze-Hessian directions: 3 to 10 of the 16, with the most negative curvature between −0.34
and −5.6 in the bulk and about −0.003 near the boundary.

Past the endpoint (t = 0.7) we also found thick, nearly critical points with negative value (R ≈ −0.03). They too are
saddles. So "thick critical points have positive value" fails for t > 1/2, and a proof of ORB cannot be first-order
uniformly in t: it has to see the negative curvature. At most critical points the descent is already visible in
squeezes that are diagonal in the eigenbasis of the local support operator Tr_other(CC^* + C^*C), one copy at a time.
This fails close to the thin boundary, where the descent uses other directions.

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

## Status

Reduction proved (above). ORB_n is a conjecture, tested numerically for n ≤ 7. A proof of ORB_n for all n would prove
NPT bound entanglement; we do not have one.

## Files

* `code/orbit_descent.py`: squeeze-orbit descent (`python code/orbit_descent.py n t starts seed ensemble`).
* `code/trace_form_check.py`: the control above.
* `logs/`: the runs in the table and the control.

## References

* T. S. Tomar, *Equality cases and every-copy results for the Werner endpoint*, Zenodo (2026): Lemma 6.1 (thin copies)
  and Theorem 6.3 (first zeros).
* L. Song and L. Chen, arXiv:2607.23416 (rank-one bound).
