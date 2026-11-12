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

## 4. Pairs separated by a ternary hyperplane (every N)

Added 12 November 2026, from an independent derivation; checked here (`code/check_coset_separation.py`).

**Proposition.** Let a ∈ F₃^N and α ≠ β in F₃. If u is supported on {x : a·x = α} and v on {x : a·x = β}, with
arbitrary real or complex entries, then

    2^N q_N(uv^T − vu^T) = 2 Σ_{S ⊆ [N]} 2^{|S|} ‖P₋^{(S)}(u ⊗ v)‖² ≥ 2 ‖u‖² ‖v‖²,

where P₋^{(S)} antisymmetrises the two factors at the copies in S. In particular the counting inequality of Section 1
holds, with margin 2|U||V|, at every N, for any two sets U, V lying in two distinct parallel affine hyperplanes of
F₃^N. Equivalently, U and V lie in distinct cosets of a ternary linear code L (take a ∈ L^⊥ with a·(p − q) ≠ 0), and
U, V may be arbitrary nonempty subsets of their cosets.

*Proof.* Let Z|x⟩ = ω^x|x⟩ with ω = e^{2πi/3}, and Z_a = Z^{a_1} ⊗ ⋯ ⊗ Z^{a_N}, a product unitary. The form q_N is
invariant under C ↦ Z_a C Z_a^†. Since Z_a u = ω^α u and Z_a v = ω^β v, the two terms transform as
Z_a(uv^T)Z_a^† = ω^{α−β} uv^T and Z_a(vu^T)Z_a^† = ω^{β−α} vu^T. These characters differ because 2(α − β) ≠ 0 in F₃,
so the cross term ⟨vec(uv^T), H_N vec(vu^T)⟩ vanishes and q_N(uv^T − vu^T) = q_N(uv^T) + q_N(vu^T). For rank one,
‖Tr_J(xy^T)‖² = ⟨x ⊗ y, F_J (x ⊗ y)⟩ with F_J the swap of the two factors at the copies in J, so
q_N(xy^T) = ⟨x ⊗ y, (I − F/2)^{⊗N} x ⊗ y⟩, and 2(I − F/2) = I + 2P₋. The term S = ∅ gives the bound. ∎

The argument needs the third root of unity. Over F₂ the two characters are both −1, and indeed U = {00, 11},
V = {01, 10} (separated by the binary parity x₁ + x₂, but by no ternary hyperplane) has a nonzero cross term
(−2, against the value 32 for 4q₂). The binary wedge zeros of the endpoint live exactly where this cancellation is
unavailable.

Scope: this is an every-N result for a structured class. Generic disjoint pairs lie in no pair of parallel
hyperplanes (a set with more than 3^{N−1} points lies in no hyperplane at all), and the integer inequality for them
stays open. In the
checks the ratio 2^N q_N / (2‖u‖²‖v‖²) was never below 2, so the stated margin is probably not sharp.

## What this does not do

Sections 1–3 are reformulations and Section 4 covers a special class. None of them proves positivity in general,
and the integer inequality of Section 1 is open.

## Checks

`python code/check_counting_form.py` (floating point; the identities are exact algebra):
1. the counting identity against q from partial traces (random U, V at N = 2, 3, 4, and two examples at N = 2);
2. the flag constants at m = 2, 4, 6;
3. the Cauchy–Schwarz bound on random real vectors;
4. the disjoint flags and the halving of the score;
5. the window flags and the factor 1/4;
6. the tag Gram at L = 4 and the decomposition of the score, with its error below 20bt.

`python code/check_coset_separation.py` (Section 4; `logs/check_coset_separation.log`): 326 random coset pairs at
N = 2, …, 5 (indicator vectors and random real vectors), cross term zero, identity to 7·10⁻¹⁶, and the control pair
above.
