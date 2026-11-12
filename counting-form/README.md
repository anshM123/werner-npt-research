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

**Any separating product symmetry works.** Only the characters were used. Let 𝒰 = U₁ ⊗ ⋯ ⊗ U_N be a product
unitary with 𝒰u = λu and 𝒰v = μv. If λμ̄ ∉ {1, −1}, the same computation gives the identity and the bound. With
𝒰 = ⊗_s diag(e^{iφ_s(0)}, e^{iφ_s(1)}, e^{iφ_s(2)}) this covers any U, V lying in two level sets of an additive phase
Φ(x) = Σ_s φ_s(x_s) mod 2π whose values θ_U, θ_V satisfy θ_U − θ_V ∉ πZ. Two instances:
- ternary hyperplanes, with φ_s = 2πa_s x_s/3;
- two different level sets of any real additive function F(x) = Σ_s f_s(x_s), with φ_s = t f_s and t small. For
  example U ⊆ {Σ_s x_s = c₁} and V ⊆ {Σ_s x_s = c₂} with integers c₁ ≠ c₂, sets of up to about 3^N/√N points.

**The general form: balanced rectangles** (added 12 November 2026, from an independent derivation; checked here).
Using all diagonal product unitaries at once gives the sharpest version. Per copy,
⟨ab|F − I/2|cd⟩ = [a = d][b = c] − ½[a = c][b = d], so the cross term is

    ⟨vec(uvᵀ), H_N vec(vuᵀ)⟩ = Σ u_x u_{x′} v_y v_{y′} ∏_s c_s,

summed over x, x′ ∈ U and y, y′ ∈ V with {y_s, y′_s} = {x_s, x′_s} as multisets at every copy s (a balanced
rectangle). Here c_s = 1 where (y_s, y′_s) = (x′_s, x_s) ≠ (x_s, x′_s), c_s = −½ where (y_s, y′_s) = (x_s, x′_s) with
x_s ≠ x′_s, and c_s = ½ where all four letters agree. Hence:
- if U and V admit no balanced rectangle, then q_N(u∧v) = 2q_N(uvᵀ) ≥ 2^{1−N}‖u‖²‖v‖², for arbitrary real
  entries;
- in general q_N(u∧v) ≥ 2(2^{−N}‖u‖²‖v‖² − 𝓔), with 𝓔 the same sum taken with absolute values.

A separating phase excludes balanced rectangles, because summing Φ over a rectangle gives 2θ_U ≡ 2θ_V. The converse
fails. U = {200, 222, 102, 221} and V = {021, 202, 120, 211} have no balanced rectangle and no real additive
separator. Every one-copy marginal of their span has rank three, the span contains no vector that is a product across
any single copy (so Section 5 does not apply either), and q₃(1_U ∧ 1_V) = 33/2. The full endpoint is exactly the
statement that the signed rectangle sum never exceeds q_N(uvᵀ). This section settles the pairs for which the sum is
empty.

The ratio must avoid −1. A binary parity gives θ_U − θ_V = π. Indeed U = {00, 11}, V = {01, 10}, which are separated
by x₁ + x₂ mod 2 but by no ternary hyperplane and no additive real function, have a nonzero cross term (−2, against
the value 32 for 4q₂). This matches the binary wedge zeros of the endpoint, where no such cancellation is available.

Scope: an every-N result for a structured class. Generic large disjoint pairs lie on no such pair of level sets,
and the integer inequality for them stays open. In the checks the ratio 2^N q_N/(2‖u‖²‖v‖²) was never below 2, so
the stated margin is probably not sharp.

## 5. Peeling a copy at which the plane contains a product vector

Added 12 November 2026; checked in `code/check_coset_separation.py` (part 4).

**Lemma.** Let u, v be real vectors on N copies, and suppose span{u, v} contains a vector w ⊗ ũ, with w ∈ R³ at one
copy s and ũ arbitrary on the other N − 1 copies. Take u = w ⊗ ũ with |w| = 1, and rotate copy s by a real orthogonal
matrix that sends w to e₀. Write v = Σ_k e_k ⊗ ṽ_k at that copy, and x∧y = xy^T − yx^T. Then

    q_N(u∧v) = ½ q_{N−1}(ũ∧ṽ₀) + 2 Σ_{k=1,2} q_{N−1}(ũ ṽ_k^T),

and the rank-one scores are at least 2^{1−N}‖ũ‖²‖ṽ_k‖². So the real skew endpoint at N − 1 copies implies it at N
copies on every such plane. A first failure needs a range plane with no vector of Schmidt rank one across any single
copy.

*Proof.* q_N is invariant under real orthogonal conjugation at one copy. Split C = u∧v into copy-s blocks C_{jk}.
The nonzero ones are C₀₀ = ũ∧ṽ₀, C_{0k} = ũṽ_k^T and C_{k0} = −ṽ_kũ^T (k = 1, 2). Since H₁ = I − ½|Ω⟩⟨Ω|,

    q_N(C) = Σ_{j≠k} q_{N−1}(C_{jk}) + Σ_j q_{N−1}(C_{jj}) − ½ q_{N−1}(Σ_j C_{jj}).

Only C₀₀ survives on the diagonal, which gives the factor ½, and transposed rank-one blocks have equal scores. ∎

Iterating gives an unconditional every-N class: planes that keep containing such a vector after each peel. It
includes every plane through a fully product vector, which is the skew case of the Wu–Zou product-ray theorem
(`support-and-projector-notes/`). The step does not extend to planes that are only thin at copy s, that is, planes
containing a vector of Schmidt rank two there. The diagonal part is then ½ q_{N−1}(C₀₀ − C₁₁), a skew operator of
rank up to four that the endpoint does not control; compare the dephasing counterexample in `contraction-notes/`.

## What this does not do

Sections 1–3 are reformulations, and Sections 4 and 5 cover special classes. None of them proves positivity in
general, and the integer inequality of Section 1 is open.

## Checks

`python code/check_counting_form.py` (floating point; the identities are exact algebra):
1. the counting identity against q from partial traces (random U, V at N = 2, 3, 4, and two examples at N = 2);
2. the flag constants at m = 2, 4, 6;
3. the Cauchy–Schwarz bound on random real vectors;
4. the disjoint flags and the halving of the score;
5. the window flags and the factor 1/4;
6. the tag Gram at L = 4 and the decomposition of the score, with its error below 20bt.

`python code/check_coset_separation.py` (Sections 4 and 5; `logs/check_coset_separation.log`):
1. 326 random coset pairs at N = 2, …, 5 (indicator vectors and random real vectors): cross term zero, identity to
   7·10⁻¹⁶;
2. the binary-parity control above;
3. 100 pairs on level sets of random integer additive functions, N = 2, 3, 4: cross term zero;
4. the peel identity of Section 5 on 60 random planes, N = 2, 3, 4, to 9·10⁻¹⁶;
5. the balanced-rectangle formula against the cross term on 80 random pairs (N = 2, 3, to 10⁻¹⁶), zero cross term on
   all 69 rectangle-free pairs, and the 3-copy example (no rectangle, q₃ = 33/2).
