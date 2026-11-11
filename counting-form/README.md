# Counting forms of the endpoint

Ansh Mishra and Aryan Senthilkumar, 11 November 2026. Nothing here is refereed.

The three reductions below come from an independent audit (10–11 November 2026). We checked every quantitative step
we could isolate numerically (`code/check_counting_form.py`, `logs/check_counting_form.log`). The limiting arguments
were read and agree with these checks, but we have not re-derived every estimate line by line.

Notation: H_N = (I − Ω/2)^{⊗N} with Ω = Σ_a |aa⟩⟨bb|, q_N(C) = ⟨vec C, H_N vec C⟩ = Σ_J (−1/2)^{|J|} ‖Tr_J C‖², and
d(a, b) = q_m(ab^T) for real vectors a, b on m copies.

## 1. Uniform superpositions on two disjoint sets

The endpoint holds for every N if and only if, for every N and all disjoint U, V ⊆ {0, 1, 2}^N,

    Σ_{S ⊆ [N]} (−1)^{|S|} 2^{N−|S|} Σ_{i,j} ( N_S^{UV}(i, j) − N_S^{VU}(i, j) )² ≥ 0,

where N_S^{UV}(i, j) counts the assignments s on the deleted copies S with (i, s) ∈ U and (j, s) ∈ V. The left side is
exactly 2^N q_N(1_U 1_V^T − 1_V 1_U^T), so the inequality is an integer statement about pairs of sets.

How a negative witness becomes one of this kind:
- The real skew reduction (`real-balanced-reduction/`, `whole-problem-routes/`, Proposition S) gives real u, v with
  q_n(uv^T − vu^T) < 0.
- Two binary flags a = (00 − 11)/√2, b = (01 + 10)/√2 have ⟨ab, H_2 ab⟩ = ⟨ab, H_2 ba⟩ = 1/2. Replacing (u, v) by
  (u ⊗ a, v ⊗ b) makes the supports disjoint and halves the score.
- Rational approximation, then integer vectors p, q.
- On a binary code of even length m put x = |+⟩^{⊗m} and y(z) = 2^{−m/2}(−1)^{Σ_{i<j} z_i z_j}. Then
  d(x, x) = 2^{−m}, d(x, y) = d(y, x) = (3/4)^m and d(y, y) = 1/2 + 2^{−m}. The vectors (x ± y)/2 have disjoint
  supports and equal amplitudes. Sums of products of them over sign tuples turn the integer coefficients into
  uniform superpositions.
- Every term containing an x flag is O((3/4)^{m/2}), by |⟨ab, H cd⟩| ≤ √(d(a, d) d(c, b)) (partial transpose and
  Cauchy–Schwarz for the positive H^Γ). The surviving term is 2^{−4t} d(y, y)^t q_n(pq^T − qp^T) < 0.

## 2. Cross-distances near N/2

It is even enough to prove the inequality, for each fixed B and all large N, for sets whose cross-distances
d(u, v) (u ∈ U, v ∈ V) all lie in [N/2 − B, N/2 + B]. The flags a₊ = (00 + 11)/√2 and b have direct value 1/2 and
crossed value −1/2. Two of them in a row have direct and crossed values 1/4 and 1/4, so appending them multiplies
the wedge score by exactly 1/4. They also add exactly 2 to every cross-distance on 4 added copies.

## 3. Near a product of qubits

A negative witness also gives, at length N = n + L, orthonormal real u_L, v_L at squared distance at most C 2^{−N/2}
from one product-binary subspace, with 2^N q_N(u_L ∧ v_L) ≤ −c < 0. With cat tags y₁ = (0^L + 1^L)/√2 and
y₂ = (w + w̄)/√2 (|w| = L/2, L ≡ 0 mod 4), and binary p, q with zero wedge score,

    u_L = √(1 − t) p ⊗ y₁ + √t u ⊗ y₂,   v_L = √(1 − t) q ⊗ y₁ + √t v ⊗ y₂.

The Gram of the tag dyads y_i y_j^T under H_L, in the order (11, 12, 21, 22), is

    [[d, 0, 0, b²], [0, b, b, 0], [0, b, b, 0], [b², 0, 0, d]],   b = 2^{−L/2},   d = 1/2 + b².

The score is (1 − t)² d q(p ∧ q) + t² d q(u ∧ v) + O(bt), and t = K b with K large makes it negative at the scale
2^{−N}.

So a uniform o(1) lower bound for 2^N q_N on such neighborhoods would prove the endpoint. The construction also says
something about method. A counterexample can sit next to a binary zero as a small, almost decoupled piece. An
argument that controls only a fixed number of orders in the distance to the binary subspace therefore cannot settle
the problem.

## What this does not do

These are reformulations. None of them proves positivity, and the integer inequality of Section 1 is open.

## Checks

`python code/check_counting_form.py` (floating point; the identities are exact algebra):
1. the counting identity against q from partial traces (random U, V at N = 2, 3, 4, and two examples at N = 2);
2. the flag constants at m = 2, 4, 6;
3. the Cauchy–Schwarz bound on random real vectors;
4. the disjoint flags and the halving of the score;
5. the window flags and the factor 1/4;
6. the tag Gram at L = 4 and the decomposition of the score, with its error below 20bt.
