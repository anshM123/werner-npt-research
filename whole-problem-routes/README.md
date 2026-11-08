# Whole-problem routes: exact reformulations of the endpoint

Ansh Mishra and Aryan Senthilkumar, 6 November 2026. Nothing here is refereed.

By `support-and-projector-notes/` and `real-balanced-reduction/`, the endpoint for every number of copies is equivalent to

    q_N(P) ≥ 0 for every odd N and every real rank-two orthogonal projector P with balanced range,

where q_N(C) = Σ_J (−1/2)^{|J|} ‖Tr_J C‖² and b_N(X, Y) = Σ_J (−1/2)^{|J|} ⟨Tr_J X, Tr_J Y⟩. This note records five
exact reformulations of that whole statement. Four restate the Werner endpoint. The fifth is a sufficient condition for
NPT bound entanglement that goes through a different state. **None of them is proved here.** Each reduces the problem to
one inequality, stated below, that would settle it.

Section 7 adds what we can prove about the fifth route. It gives explicit witnesses that bound the n-copy positivity
threshold of the map from below at every n, and the exact threshold at two copies. It also states Conjecture QT: the
threshold at n copies is (1 + 2^{1−1/n})/4. That conjecture would settle the route, and so the existence of NPT bound
entanglement. It is open from three copies on, and proved for every n on a real plane, where the threshold is exactly
(1 + 2^{1−1/n})/4 (Theorem 7.18). At the parameter that matters, a = 3/4, inputs supported on a plane at a copy are
entanglement broken there, so every copy where the input is thin can be removed (Proposition 7.8). The map is positive
at six copies (Theorem 7.32, computer-assisted), so only inputs with full support at every copy, on at least seven
copies, remain (Corollary 7.14). With a different parameter at each copy, the
plane threshold becomes the rule ∏_{s∈S}(2a_s − ½) ≥ ½ for every set S of copies (Theorem 7.19). We conjecture the same
rule for the qutrit map (Conjecture QT_m, proved for two copies), and more: that the least value is always attained on
planes (Conjecture PM). The extreme points of this rule are exactly 2-positivity statements (Proposition 7.20), and the
plane witnesses are local minima in the full qutrit space (Proposition 7.22). Finally, from any critical point of the
value, squeezing one copy toward a plane moves the value monotonically, with the sign set by one contraction
inequality (Proposition 7.23). So a failure at a = 3/4 would have to be a strict minimum along every one-copy squeeze
(Corollary 7.24). Such points exist away from the minimum (a balanced pair of spin singlets at four copies,
Observation 7.26), so squeezes alone cannot finish the argument. In every search, inputs with full support at every
copy sit far above the product value (Observation 7.27), which suggests a multiplicative bonus for thickness
(Conjecture TC). The bonus holds at two copies (Proposition 7.28) and, by an exact computer-assisted
certificate with a variable for the least eigenvector of each marginal, at three copies (Theorem 7.29). It is open for
every larger n, and Theorem 7.29 does not extend to more copies by itself. For every n the bonus holds against product
test vectors, for every input (Proposition 7.30). Peeling one copy at a time reduces its version for mixed inputs to a
bound on the vertex map 𝒫_{1/2} ⊗ 𝒫_{3/4}^{⊗(n−1)} at a suitable copy. Positivity of that map is real 2-positivity of
𝒫_{3/4}^{⊗(n−1)}, and that is where the whole difficulty of this route sits (Remark 7.31).

## 1. The cross ratio, and fusion of planes

For a real orthonormal pair (u, v) on n copies put z = (u + iv)/√2 (so zᵀz = 0), C_z = zz*, P = uuᵀ + vvᵀ = C_z + C_z̄ and

    Δ_n = 4 q_n(C_z),   ε_n = 4 b_n(C_z, C_z̄),   so that   q_n(P) = (Δ_n + ε_n)/2.

Both numbers depend only on the plane, and Δ_n = 4⟨z ⊗ z|(I − F/2)^{⊗n}|z ⊗ z⟩ ≥ 2^{2−n}. In real terms, with
𝓛_n = (Id − T/2)^{⊗n} (T the transpose on one copy; ⟨X, 𝓛_n X⟩ is a positive definite norm on real matrices),
H₁ = uuᵀ − vvᵀ, H₂ = uvᵀ + vuᵀ and J = uvᵀ − vuᵀ: Δ_n = ‖H₁‖²_𝓛 + ‖H₂‖²_𝓛 and ε_n = ‖P‖²_𝓛 − ‖J‖²_𝓛.

**Proposition F (fusion).** For real orthonormal pairs (u₁, v₁) on n₁ copies and (u₂, v₂) on n₂ copies let
u* = (u₁ ⊗ u₂ − v₁ ⊗ v₂)/√2 and v* = (u₁ ⊗ v₂ + v₁ ⊗ u₂)/√2. They are real and orthonormal,
P* = (P₁ ⊗ P₂ − J₁ ⊗ J₂)/2, and

    Δ*_{n₁+n₂} = Δ_{n₁} Δ_{n₂} / 4,   ε*_{n₁+n₂} = ε_{n₁} ε_{n₂} / 4.

If both planes are balanced, so is the fused plane.

*Proof.* z* = z₁ ⊗ z₂, and q and b are multiplicative on tensor products because M^{⊗(n₁+n₂)} = M^{⊗n₁} ⊗ M^{⊗n₂}. Each
copy belongs to one factor, so the one-copy compressions on the four-dimensional product of the two planes are scalar,
and they stay scalar on the plane of (u*, v*). ∎

**Corollary F.** Let Γ_n be the supremum of |ε_n|/Δ_n over real balanced planes on n copies. If log Γ_n = o(n) along
the odd n, then Γ_n ≤ 1 for every odd n, and the endpoint holds for every number of copies. Conversely the endpoint
gives |ε_n| ≤ Δ_n, since C_z + C_z̄ and C_z − C_z̄ both have rank two.

*Proof.* A plane with ratio r > 1 on m copies (m odd), fused with itself k times (k odd), is a balanced plane on km
copies with ratio r^k. ∎

So a polynomial bound, or any subexponential bound, on the cross ratio would settle the endpoint. The bound that comes
for free is Γ_n ≤ 2^n.

## 2. Monotone minima and padding

Let μ_n be the minimum of q_n(C) over C of rank at most two with ‖C‖_F = 1.

**Proposition P.** μ_{n+1} ≤ μ_n and −1/2 ≤ μ_n ≤ 0. Hence the endpoint holds for every n if and only if μ_N ≥ −ε_N
with ε_N → 0 along some unbounded set of N.

*Proof.* q_{n+1}(C ⊗ |0⟩⟨1|) = q_n(C), and C ⊗ |0⟩⟨1| has the same rank and norm. The lower bound is the least
eigenvalue of (I − ½|Ω⟩⟨Ω|)^{⊗n}; the upper bound is the one-copy zero |0⟩⟨0| + |1⟩⟨1|, padded. If μ_n = −δ < 0, then
μ_N ≤ −δ for all N ≥ n. ∎

The same works inside the real balanced class, with the real GHZ projectors D_r = |g_r⟩⟨g_r|,
g_r = (|0^r⟩ + |1^r⟩ + |2^r⟩)/√3 (balanced for r ≥ 2), as padding. At a general coefficient t,

    q_{r,t}(D_r) = ((1 − t)^r + 2(1 + (−t)^r))/3,   so at t = 1/2:   2/3 + 2^{−r} (r even),  2/3 − 2^{−r}/3 (r odd).

Let b_N be the minimum of q_N(P)/2 over real balanced rank-two projectors at odd N. If b_m = −δ < 0, then
b_{m+r} ≤ −δ q_r(D_r) ≤ −2δ/3 for every even r ≥ 2. The balanced three-copy code span{Σ_a |aaa⟩, |012⟩ + |120⟩ + |201⟩}
has q₃(P)/2 = 3/4, and padding it gives b_N ≤ (3/4)(11/12)^{(N−3)/2}. So the endpoint holds if and only if
lim sup_{N odd} b_N = 0, and a lower bound b_N ≥ −ε_N with ε_N → 0 along one unbounded set of odd N would suffice.

## 3. The real skew form

**Proposition S.** Let C = σ₁u₁v₁ᵀ + σ₂u₂v₂ᵀ be real of rank two on n copies (singular value form, s = σ₁ + σ₂), and on
r further copies put

    a = (√σ₁ u₁ ⊗ |0^r⟩ + √σ₂ v₂ ⊗ |1^r⟩)/√s,   b = (√σ₂ u₂ ⊗ |0^r⟩ − √σ₁ v₁ ⊗ |1^r⟩)/√s,   A_r = abᵀ − baᵀ.

Then a, b are real and orthonormal, and with J_U = u₁u₂ᵀ − u₂u₁ᵀ, J_V = v₁v₂ᵀ − v₂v₁ᵀ,

    q_{n+r}(A_r) = s^{−2} [ 2 q_n(C) + σ₁σ₂ 2^{−r} ( q_n(J_U) + q_n(J_V) + 2(−1)^r b_n(J_U, J_V) ) ].

The bracketed error is at most 8σ₁σ₂2^{−r} in absolute value. So a negative real C gives a negative real skew rank-two
matrix after r added copies, for r large and of either parity. With the realification and balancing of
`real-balanced-reduction/` (which keep skew symmetry), the endpoint for every n is equivalent to
**q_N(J) ≥ 0 for every odd N and every real skew rank-two J with balanced support.** Since iJ is Hermitian of rank two
and q_N(iJ) = q_N(J), this is a Hermitian statement too.

## 4. A spectral form on spherical harmonics

Put V = (R³)^{⊗N}. For S ⊆ [N] let L_{ij,S} = Σ_a (x_{ai}∂_{aj} − x_{aj}∂_{ai}) rotate the coordinates of V_{Sᶜ}
identically for every index a of V_S, let 𝖢_{Sᶜ} = −Σ_{i<j} L_{ij,S}² and 𝒜_N = Σ_S (−1/2)^{|S|} 𝖢_{Sᶜ}.

**Proposition H.** (i) For real x, y,
Σ_S (−1/2)^{|S|} Σ_{i<j} [Σ_a (x_{ai}y_{aj} − x_{aj}y_{ai})]² = ½ q_N(xyᵀ − yxᵀ).
(ii) q_N ≥ 0 on all real skew rank-two matrices if and only if 𝒜_N is positive semidefinite on every space ℋ_m(V) of
harmonic polynomials of degree m ≥ 1.
(iii) On m tensor replicas,
𝒜_N = m[(5/2)^N − (1/2)^N] I + 2 Σ_{a<b} [ ⊗_s (F_{ab,s} − I/2) − ⊗_s (Ω_{ab,s} − I/2) ],
with F_{ab,s} the swap and Ω_{ab,s} the unnormalised contraction of replicas a, b at copy s.
(iv) On ℋ_m(V), 𝒜_N ⪰ (m/2^N)[(26 − m)5^{N−2} + m − 2] I.

*Proof.* (i) The inner sum is ½‖Tr_S(xyᵀ − yxᵀ)‖². (ii) Rotations preserve degree, and
⟨f, 𝒜_N f⟩ = ∫ Σ_S (−1/2)^{|S|} Σ_{i<j} (L_{ij,S} f)², whose integrand is the form of (i) at (x, ∇f(x)). So
nonnegativity of the form gives 𝒜_N ⪰ 0. Conversely (sketch), a negative value at a tangent covector gives, at high
frequency, an oscillating test function with negative energy; approximating it by polynomials and decomposing into
harmonic degrees gives a negative block. (iii) For so(r) on
m replicas the Casimir is m(r − 1)I + 2Σ_{a<b}(F_ab − Ω_ab), because Σ_{i<j} L_ij ⊗ L_ij = Ω − F; summing over S with
r = 3^{N−|S|} gives the formula. (iv) On ℋ_m the global swap of two replicas is +1, so ⊗_s(F_{ab,s} − I/2) has
eigenvalues 2^{−N}3^r with r even, all ≥ 2^{−N}; the global contraction annihilates harmonic tensors, so
⊗_s(Ω_{ab,s} − I/2) ⪯ ¼(5/2)^{N−2} there. Insert both bounds into (iii). ∎

The bound (iv) is positive for every N exactly when m ≤ 26, and the four-copy theorem of `four-replica-notes/` (applied
to iJ, padded by |0⟩⟨0|) gives every block for N ≤ 4. What is left of this route is N ≥ 5 and m ≥ 27. Bounding each
replica pair separately loses a factor m², and nothing here removes that loss.

## 5. A tensor-stable map

For a real parameter a let 𝒫_a(X) = a Tr(X) I₃ − ¼(X + Xᵀ) on M₃(C). Its Choi matrix C_a = aI₉ − ¼(F + |Ω⟩⟨Ω|) is
invariant under partial transpose, with eigenvalues a − 1 (on Ω), a − ¼ (five times) and a + ¼ (three times).

**Proposition M.**
1. 𝒫_a is positive if and only if a ≥ 1/2, 2-positive if and only if a ≥ 3/4, and completely positive (equivalently
   completely copositive) if and only if a ≥ 1. 𝒫₁ is entanglement breaking:
   𝒫₁(X) = (15/4) ∫_{S²} (I − uuᵀ)(uᵀXu) dμ(u).
2. The set of a for which 𝒫_a is tensor-stable positive (𝒫_a^{⊗n} positive for every n) is an interval [a*, ∞) with
   3/4 ≤ a* ≤ 1.
3. **If a* < 1, NPT bound entanglement exists.** Precisely: if 𝒫_a is tensor-stable positive for some a < 1, the Werner
   state proportional to (a − 1/16)I − (5/16)F is NPT and undistillable from every number of copies. For a = 9/10 this
   is ρ = (67I − 25F)/528, with ρ^Γ having the eigenvalue −1/66.
4. Put m = a − 1/2 and μ_n(a) = min ⟨y, 𝒫_a^{⊗n}(xx*) y⟩ over unit x, y. Then 𝒫_a is tensor-stable positive if and
   only if μ_n(a) = m^n for every n, and if and only if s_N(z) = zᵀ𝒫_a^{⊗N}(zzᵀ)z ≥ 0 for all real z at every N in some
   unbounded set (or ≥ −ε_N a^N‖z‖⁴ with ε_N → 0 along it).

*Proof.* 1. On a product vector the Choi expectation is a − ¼(|⟨x, y⟩|² + |xᵀy|²) ≥ a − 1/2. A Schmidt-rank-two unit
vector has ⟨F⟩ ≤ 1 and |⟨Ω|ψ⟩|² ≤ 2, with both attained at (|11⟩ + |22⟩)/√2. The integral follows from
E[u_iu_j] = δ_ij/3 and E[u_iu_ju_ru_s] = (δ_ijδ_rs + δ_irδ_js + δ_isδ_jr)/15.
2. Adding tD (D(X) = Tr(X)I) keeps tensor-stable positivity: expand and trace out the D-positions. For a < 3/4,
compression to a real qubit gives a positive qubit map that is neither completely positive nor copositive, which cannot
be tensor-stable: a tensor-stable positive map with a two-dimensional input or output space is completely positive or
completely copositive [MRW, Thm. 2]. Here [MRW] is Müller-Hermes, Reeb, Wolf, J. Math. Phys. 57, 015202 (2016),
arXiv:1502.05630.
3. By the same theorem every compression of 𝒫_a^{⊗n} to a qubit input is completely positive or copositive, and
T ∘ 𝒫_a^{⊗n} = 𝒫_a^{⊗n} makes it completely positive in either case; so 𝒫_a^{⊗n} is 2-positive. Its Choi matrix C_a^{⊗n} is
then positive on Schmidt-rank-two vectors. The U ⊗ U twirl of C_a is (a − 1/16)I − (5/16)F (the twirl of |Ω⟩⟨Ω| is
(I + F)/4), and since C_a^Γ = C_a, the partial transpose of the twirled state is an average of local unitary conjugates of
C_a. So (W_a^Γ)^{⊗n} is positive on Schmidt-rank-two vectors for every n, which is n-copy undistillability; W_a^Γ has the
eigenvalue a − 1 < 0 on Ω. ([MRW, Thm. 4] already gives NPT bound entanglement from any tensor-stable positive map
that is neither completely positive nor completely copositive; the argument here names the state.)
4. Write 𝒫_a = mD + R + T ∘ R with R(X) = ¼(Tr(X)I − Xᵀ), completely positive with the Kraus operators
(|i⟩⟨j| − |j⟩⟨i|)/2, each of rank two. If Φ = 𝒫_a^{⊗n} is 2-positive (as in 3) then R ⊗ Φ and (T ∘ R) ⊗ Φ are positive,
so (𝒫_a ⊗ Φ)(X) ⪰ m I ⊗ Φ(Tr₁X) ⪰ m μ_n(a) I; product vectors attain m^n. For the scalar criterion, with
z_k = x ⊗ |0⟩^{⊗k} + y ⊗ |1⟩^{⊗k},

    s_{n+k}(z_k) = m^k [s_n(x) + s_n(y)] + 2a^k b_n(x, y) + 4(−1/4)^k κ_n(x, y),

b_n(x, y) = yᵀΦ_n(xxᵀ)y, κ_n(x, y) = xᵀΦ_n(xyᵀ)y. Dividing by 2a^k and letting k → ∞ gives b_n(x, y) ≥ 0, i.e.
positivity on real rank-one inputs, which extends to all positive inputs because Φ_n commutes with the transpose. ∎

The open step of this route is a lower bound μ_n(a) ≥ m^n for every n at one a < 1 (at a = 9/10: μ_n ≥ (2/5)^n).
At a = 3/4 it is proved for n ≤ 5 (Theorem 7.13).

[MRW, Thm. 3] also shows that a non-trivial tensor-stable positive map exists if and only if one exists in the
one-parameter family 𝒲_p ⊗ (ϑ ∘ 𝒲_p) on M_d ⊗ M_d, p ∈ [−1, 0), where
𝒲_p(X) = ((d − p) Tr(X) I − (1 − dp) Xᵀ)/(d² − 1) and ϑ is the transpose. The map 𝒫_a is a different candidate: it acts
on one qutrit, it is invariant under the transpose on either side, and its tensor powers have the extra O(3)
symmetry used throughout Section 7.

For context: tensor-stable positive maps that are neither completely positive nor completely copositive are known to
exist over the hypercomplex numbers, and tensor-stable positivity is undecidable for the matrix multiplication tensor
(van der Eyden, Netzer, De las Cuevas, J. Phys. A 55, 264006 (2022), arXiv:2110.02113). Over the complex numbers no
such map is known; one would give NPT bound entanglement [MRW, Thm. 4].
It does not need the Werner endpoint itself. A related qutrit family was studied by Filippov and Magadov, who
considered spin-1 maps that rescale the antisymmetric part of X and remove its traceless symmetric part, and showed
that positivity of the tensor square forces the square of the map to be completely positive [S. N. Filippov,
K. Yu. Magadov, Lobachevskii J. Math. 39, 65 (2018), arXiv:1907.11878]. The map 𝒫_a does the opposite: it removes the
antisymmetric part and multiplies the traceless symmetric part by −½.

## 6. Numerical observations on the map route

Everything in this section is floating point, and nothing above rests on it. Write m₂(n) for the minimum of
⟨χ|C_a^{⊗n}|χ⟩ over unit vectors χ of Schmidt rank at most two across inputs | outputs, so that 𝒫_a^{⊗n} is
2-positive exactly when m₂(n) ≥ 0.

**A structural identity.** C_a = (a − ½)I + N with N = −½P_Ω + ¼P_sym0 + ¾P_anti, where P_Ω, P_sym0 and P_anti are
the projectors onto Ω, the traceless symmetric and the antisymmetric part of C³ ⊗ C³. N does not depend on a. For one
copy the product minimum is a − ½ and m₂(1) = a − ¾.

**Observation 1 (the threshold does not move).** At a = 3/4 our searches found no Schmidt-rank-two vector with a
negative value for n = 2, …, 5. At a = 0.74 they find negative values at once (n = 1, 2).

**Observation 2 (a multiplicative margin).** For a = 0.8 and a = 0.9 and n = 1, …, 5, the minimum found equals
(a − ¾)(a − ½)^{n−1} to eight digits. This is the product value: one copy at its Schmidt-rank-two minimum, the others at
real product vectors. Entangled vectors never did better.

**Conjecture MM.** For every a ≥ 3/4 and every n, m₂(n) = (a − ¾)(a − ½)^{n−1}.

(Proposition 7.11 below proves m₂(2) ≥ 0, that is 2-positivity of 𝒫_a^{⊗2}, for a ≥ (23 + 6√35)/68 ≈ 0.8602, and
m₂(n) ≥ 0 for n = 3, 4, 5 when a ≥ 0.91, 0.93, 0.95.)

Conjecture MM at one a < 1, or even only m₂(n) ≥ 0 for every n at one a < 1, makes every tensor power of 𝒫_a
2-positive, hence 𝒫_a tensor-stable positive (Proposition M, item 4), and so gives NPT bound entanglement (item 3), with
the Werner state proportional to (a − 1/16)I − (5/16)F.

**Two copies.** For unit χ, ⟨χ|C_a^{⊗2}|χ⟩ − (a − ¾)(a − ½) = (a − ½)(⟨N ⊗ I⟩ + ⟨I ⊗ N⟩ + ¼) + ⟨N ⊗ N⟩ is affine in a.
So Conjecture MM at n = 2, for all a ∈ [¾, 1], is equivalent to two statements: C_{3/4}^{⊗2} is 2-block-positive, and
C_1^{⊗2} − ⅛I is 2-block-positive. At a = 3/4, 2C_{3/4} = M + P_anti with M = I − ½ΩΩ* the operator of the Werner
endpoint: the 3/4 map is the endpoint form with every antisymmetric copy weighted twice.

**Calibration.** The same search finds the failure of 2-positivity of 𝒫_{9/10} ⊗ 𝒮, where 𝒮 has Choi matrix
I₉ − Σ_{i=1}^{5} v_iv_iᵀ/|v_i|² for the product vectors v₁ = e₀ ⊗ (e₀ − e₁), v₂ = (e₀ − e₁) ⊗ e₂, v₃ = e₂ ⊗ (e₁ − e₂),
v₄ = (e₁ − e₂) ⊗ e₀, v₅ = (e₀ + e₁ + e₂)^{⊗2}. 𝒮 is completely positive and copositive, yet the tensor product is not
2-positive (value −1.6·10⁻³). So 2-positivity together with transpose invariance is not an invariant that an induction
can carry; any proof of Conjecture MM must use more of the structure of 𝒫_a.

`code/map_numerics.py` (`logs/map_numerics.log`; `python map_numerics.py big` gives `logs/map_numerics_big.log`)
reproduces these numbers: alternating minimisation over the two Schmidt pairs, random complex starts.

## 7. The positivity threshold of the map at each number of copies

| # | Statement | Status |
|---|---|---|
| 7.1 | 𝒫_a^{⊗n} is not positive for a < a_n = (1 + 2^{1−1/n})/4, every n (explicit witness) | **proved** |
| 7.2 | n = 2: positive iff a ≥ (1 + √2)/4; not decomposable below (7 + 2√15)/22 | **proved** |
| 7.3 | decomposability of 𝒫_a^{⊗n} is a linear program; its thresholds tend to 1 | reduction **proved**; thresholds numerical (LP) |
| 7.4 | the threshold is a_n for n ≤ 6 (searches) and n = 9, 15, 21 (matrix product states); general matrix product states for n ≤ 12 find nothing below the product value at a_n + 0.002 and at 3/4 | numerical |
| QT | for a ≥ 0, 𝒫_a^{⊗n} is positive iff a ≥ a_n, every n (would give NPT bound entanglement) | **open** from n = 3 (at a = 3/4 proved for n ≤ 5, see 7.13; on real planes proved for every n, see 7.18) |
| 7.7–7.10 | at a ≥ 3/4 copies where the input or the output vector is thin are free; only inputs thick at every copy, n ≥ 3, remain | **proved** |
| 7.11 | 𝒫_a^{⊗n} is 2-positive for a ≥ 0.8602 (n = 2) and a ≥ 0.91, 0.93, 0.95 (n = 3, 4, 5) | **proved** (n = 2 by hand; n = 3–5 exact certificates) |
| 7.12 | the two-copy Werner endpoint and the two-copy case of the map at 3/4 are the ends of one family | reduction **proved**; the range [1, 2] numerical |
| 7.13 | 𝒫_{3/4}^{⊗n} is positive for n = 3, 4, 5, with least output eigenvalue exactly 4^{−n}; 𝒫_a^{⊗3} is positive for a ≥ 81/125 = 0.648 | **proved** (computer-assisted: exact rational certificates, two separate exact verifiers) |
| 7.14 | a failure of tensor-stable positivity of 𝒫_{3/4} needs n ≥ 7 copies, with the input and the output vector thick at seven or more common copies | **proved** (the step from six to seven uses Theorem 7.32) |
| 7.15 | the second level of the symmetric-extension hierarchy is exact at a = 3/4 for n ≤ 5, but not near a_n; its threshold ℓ_n increases with n and is ≈ 0.720 at n = 5; at n = 6 it still certifies a = 3/4 (Theorem 7.32), and in floating point it is no longer exact there | numerical (ℓ_{n+1} ≥ ℓ_n **proved**) |
| 7.16 | on permutation-symmetric inputs and outputs: positivity fails below a_n for every n; at a = 3/4 the least value is exactly 4^{−n} for n ≤ 21, attained only at real products | **proved** (the first part every n; the second computer-assisted: n ≤ 9 by two separate exact checks, 10 ≤ n ≤ 21 by one); the threshold a_n on this sector numerical (n ≤ 7) |
| 7.17 | for a ≥ 3/4 and every n, inputs that become thin one copy at a time (for example every open matrix product state of bond dimension two, such as the AKLT chain) have least output eigenvalue at least the product value (a − ½)^n | **proved** (by hand; with Theorem 7.13, up to five copies may stay thick) |
| 7.18 | on a real plane (the Pauli family (1, t, 0, t) of Filippov and Magadov) the n-copy threshold is exactly a_n, every n, and positivity equals decomposability there | **proved** (explicit decomposition) |
| 7.19 | on real planes with different parameters: 𝒬_{a₁} ⊗ ⋯ ⊗ 𝒬_{a_n} is positive iff decomposable iff ∏_{s∈S}(2a_s − ½) ≥ ½ for every S, every n | **proved** (explicit decomposition) |
| QT_m | the same condition decides positivity of 𝒫_{a₁} ⊗ ⋯ ⊗ 𝒫_{a_n} | **open** for n ≥ 3 (necessity and n = 2 proved; n = 3, 4 numerical) |
| PM | the least value of 𝒫_{a₁} ⊗ ⋯ ⊗ 𝒫_{a_n} (a_s ≥ ½) is attained on planes, with an explicit formula m(a); implies QT_m | **open** for n ≥ 3 (n = 2 proved; n = 3, 4 numerical; at a = ¾ proved for n ≤ 5) |
| 7.20 | the vertices of QT_m are 2-positivity: 𝒫_{1/2} ⊗ Φ is positive iff Φ is 2-positive on real inputs | **proved** |
| 7.21 | the complex (half-plane, Lee–Yang) version of QT_m is false from n = 5 (with every Re a_s > 3/4, from n = 7) | **proved** (explicit zeros) |
| 7.22 | at every point of the facet ∏(2a_s − ½) = ½, the plane witness is a second-order local minimum in the full qutrit space, every n | **proved** |
| 7.23 | from a joint critical point of the value, squeezing one copy toward a plane moves the normalised value monotonically all the way to the compression; one contraction inequality gives the sign; the trace of the one-copy squeeze Hessian has a closed form (every n, every a) | **proved** |
| 7.24 | a smallest failure at a ≥ 3/4 would be a strict minimum along every one-copy squeeze (thick critical points of that kind exist, Observation 7.26, so this is not a reduction) | **proved** |
| 7.25 | the thick critical points found at n = 3, 4 are saddles along a one-copy squeeze, and the squeeze Hessian trace is negative there; some one-copy compression lowers the least output eigenvalue in random and adversarial tests (n = 2, 3, 4), but not at the four-copy partner vector of Observation 7.26, where every one-copy compression raises it | numerical (critical points found by search; a failure would need n ≥ 7) |
| 7.26 | a balanced critical pair at n = 4 (spin singlets, value (151 − √6805)/1536) rises along every one-copy and every joint squeeze, yet is a saddle | identities **exact**, Hessians numerical |
| 7.27 | inputs thick at every copy sit far above the product value; the per-copy factor under a thickness floor hardly depends on n; singlet minima are 7.7 to 109 times 4^{−n} for n ≤ 7 | numerical |
| 7.28 | at two copies, λ_min(𝒫_{3/4}^{⊗2}(yyᵀ)) ≥ (¼ + ½τ)² for every input y of thickness τ | **proved** |
| 7.29 | at three copies, λ_min(𝒫_{3/4}^{⊗3}(yyᵀ)) ≥ (1/64)(1 + (98/27)Σ_sτ_s) ≥ ∏_s(¼ + ½τ_s): Conjecture TC at n = 3 | **proved** (computer-assisted, exact rational certificate, two verifiers) |
| 7.30 | for every input, pure or mixed, and every real product test vector, the value is at least 4^{−n}3^{2Σ_sτ_s} ≥ ∏_s(¼ + ½τ_s), every n; so TC holds whenever the least output eigenvalue is attained at a product vector (for example for inputs diagonal in a product basis) | **proved** |
| 7.31 | peeling one copy (𝒫_{3/4} = ¼D + 𝒫_{1/2}): TC for mixed inputs at n follows from n − 1 and a bound on the vertex map 𝒫_{1/2} ⊗ 𝒫_{3/4}^{⊗(n−1)} at one copy; the operator form of the peeling is equivalent to real 2-positivity | identity **proved**; the bounds numerical (n ≤ 4); level 2 does not certify the first open vertex case |
| 7.32 | 𝒫_{3/4}^{⊗6}(xx*) ⪰ 4^{−7}‖x‖² I for every x: 𝒫_{3/4}^{⊗6} is positive; in floating point the second level is not exact at six copies (optimum ≈ 0.563·4^{−6}) | the bound **proved** (computer-assisted: exact rational certificate, two separate exact verifiers); the non-exactness numerical |
| 7.33 | Conjecture TC for every separable input, pure or mixed, every n: 𝒫_{3/4}^{⊗n}(σ) ⪰ 4^{−n}3^{2Σ_sτ_s(σ)} I; convexity cannot carry this to entangled inputs | **proved** |
| TC | λ_min(𝒫_{3/4}^{⊗n}(yyᵀ)) ≥ ∏_s(¼ + ½τ_s(y)), τ_s(y) the least eigenvalue of y's one-copy marginal; it contains QT at 3/4 | **open** (n = 2, 3 proved: 7.28, 7.29; every n against product test vectors: 7.30; n = 4 adversarial searches) |

Section 6 was about 2-positivity. Here we look at positivity itself, copy by copy. Write

    a_n = (1 + 2^{1 − 1/n})/4,   so a_1 = 1/2, a_2 = (1 + √2)/4 ≈ 0.6036, a_3 ≈ 0.6468, a_4 ≈ 0.6704, …,

an increasing sequence with limit 3/4 (and 3/4 − a_n ≈ (ln 2)/(2n)). Real vectors suffice throughout: 𝒫_a^{⊗n}(xx*)
depends only on Re(xx*) and is real symmetric.

**Theorem 7.1 (witnesses at every n).** For every n ≥ 1 and every a with 0 ≤ a < a_n, 𝒫_a^{⊗n} is not positive.
Explicitly, let
p = (e₁ + ie₂)/√2, let x = (p^{⊗n} + p̄^{⊗n})/√2 (a real unit vector, the GHZ state in the eigenbasis of
Y = e₁e₂ᵀ − e₂e₁ᵀ), and let y = x for odd n and y = (Y ⊗ I^{⊗(n−1)})x for even n. Then

    ⟨y| 𝒫_a^{⊗n}(xxᵀ) |y⟩ = 2^{−n} [ (2a − 1/2)^n − 1/2 ].

*Proof.* Everything lies in span{e₁, e₂} at every copy, where 𝒫_a compresses to 𝒬_a(X) = a Tr(X) I₂ − ¼(X + Xᵀ). In the
real basis I, Z, X, Y of M₂(R), 𝒬_a is diagonal with eigenvalues 2a − ½, −½, −½, 0. Expanding xxᵀ = 2^{−n} Σ_P c_P P over
real Pauli strings (c_P = ⟨x|P|x⟩, and likewise d_P = ⟨y|P|y⟩) gives

    ⟨y|𝒬_a^{⊗n}(xxᵀ)|y⟩ = 2^{−n} Σ_{P ∈ {I,Z,X}^n} (2a − 1/2)^{n − |P|} (−1/2)^{|P|} c_P d_P.

Since p*p̄ = ⟨p|Z|p⟩ = ⟨p|X|p⟩ = 0, ⟨p|Z|p̄⟩ = 1 and ⟨p|X|p̄⟩ = −i, one gets c_I = 1, c_P = cos(π #X(P)/2) for
P ∈ {Z, X}^n, and c_P = 0 for every other non-identity string in {I, Z, X}^n. So exactly 2^{n−1} strings, all of full
weight, have c_P = ±1. For odd n, d = c. For even n, conjugation by Y on the first copy flips the sign of Z and X there,
so d_P = −c_P on all full-weight strings. In both cases the sum is (2a − ½)^n − 2^{n−1}·2^{−n} = (2a − ½)^n − ½, which
is negative for 0 ≤ a < a_n. ∎

This makes the bound a* ≥ 3/4 of Proposition M explicit. The witnesses live on one plane per copy, where 𝒫_a is a
qubit map that is completely positive exactly when a ≥ 3/4.

On a real plane, 𝒬_a/(2a − ½) is the unital Pauli map with eigenvalues λ₁ = λ₃ = −1/(4a − 1) and λ₂ = 0 (on X, Z and
Y). Filippov and Magadov studied n-tensor-stable positivity of Pauli maps [FM, J. Phys. A 50, 055301 (2017),
arXiv:1604.01716]. Their necessary condition from n-qubit GHZ states (Proposition 7 there) gives, for this family,
|λ₁| ≤ 2^{−1+1/n}, which is exactly a ≥ a_n (their Example 3 treats this family). So the witnesses of Theorem 7.1 are
their GHZ obstruction, carried from the plane to the qutrit map. They show that the condition is also sufficient for
this family at n = 2 (with an explicit decomposition, their Example 1), and give numerical evidence at n = 3.
Theorem 7.18 below proves sufficiency for every n.

**Theorem 7.2 (two copies).** 𝒫_a^{⊗2} is positive if and only if a ≥ (1 + √2)/4. For
(1 + √2)/4 ≤ a < (7 + 2√15)/22 ≈ 0.6703 it is positive but not decomposable.

*Proof.* By the covariance 𝒫_a(OXOᵀ) = O𝒫_a(X)Oᵀ (O real orthogonal) and the singular value decomposition, it is enough
to take x = Σ_i s_i e_i ⊗ e_i with Σ s_i² = 1. With E_ij = e_ie_jᵀ and S_ij = E_ij + E_ji,

    𝒫_a^{⊗2}(xxᵀ) = Σ_i s_i² (aI − ½E_ii)^{⊗2} + ⅛ Σ_{i<j} s_i s_j S_ij ⊗ S_ij,

which is block diagonal: one block on span{|kk⟩} and one on each span{|ij⟩, |ji⟩}, i < j.
- The first block is M_A = diag(a² − (a − ⅛)s_k²) + ⅛ssᵀ, and with b = a − ¼,
  zᵀ(M_A − (a − ½)²)z = Σ_{k<l} [ b(s_l²z_k² + s_k²z_l²) + ¼ s_ks_l z_kz_l ] ≥ (2a − ¾) Σ_{k<l} |s_ks_lz_kz_l| ≥ 0.
- The block on {|ij⟩, |ji⟩} has eigenvalues a² − (a/2)(s_i² + s_j²) ± ⅛ s_is_j ≥ a² − a/2 − 1/16 = (a − ¼)² − ⅛, which
  is ≥ 0 exactly when a ≥ (1 + √2)/4. At s = (e₁ + e₂)/√2 this value is attained, which is the case n = 2 of Theorem 7.1.

For non-decomposability, let P_Ω, P_s0, P_an be the projectors onto Ω, the traceless symmetric and the antisymmetric part
of C³ ⊗ C³. The operator ρ = (5 P_Ω⊗P_Ω + 10 P_Ω⊗P_an + 4 P_s0⊗P_s0 + 10 P_an⊗P_Ω)/29 is positive with positive partial
transpose, and Tr(C_a^{⊗2} ρ) = 15(44a² − 28a − 1)/116, which is negative for a < (7 + 2√15)/22. ∎

**Proposition 7.3 (decomposability is a linear program).** 𝒫_a^{⊗n} is decomposable (its Choi matrix is C₁ + C₂^Γ with
C₁, C₂ ⪰ 0) if and only if there is c₂ ≥ 0 in R^{3^n} with c_C − G^{⊗n} c₂ ≥ 0, where
c_C = (a − 1, a − ¼, a + ¼)^{⊗n} and

    G = [[1/3, 5/3, −1], [1/3, 1/6, 1/2], [−1/3, 5/6, 1/2]]

is the partial transpose written in the basis P_Ω, P_s0, P_an (G² = I). By symmetry c₂ may be taken to depend only on how
many copies carry each projector, which leaves about n²/2 variables.

*Proof.* Twirling over O(3)^n commutes with the partial transpose and fixes C_a^{⊗n}, so C₁ and C₂ may be taken
invariant. The commutant of O ⊗ O on C³ ⊗ C³ is spanned by the three projectors, so invariant operators are diagonal in
their products, and positivity means nonnegative coefficients. ∎

Solving this program (floating point, HiGHS) gives the decomposability threshold of 𝒫_a^{⊗n}:

| n | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|---|---|---|
| decomposable from a = | 0.5 | 0.6703 | 0.7562 | 0.8069 | 0.8401 | 0.8637 | 0.8947 | 0.9142 | 0.9276 |
| a_n | 0.5 | 0.6036 | 0.6468 | 0.6704 | 0.6853 | 0.6954 | 0.7085 | 0.7165 | 0.7219 |

The thresholds increase toward 1, and they have to. A map all of whose tensor powers are decomposable is completely
positive or completely copositive [A. Müller-Hermes, J. Math. Phys. 59, 102203 (2018), arXiv:1805.11570, Theorem 1.1],
and 𝒫_a is neither for a < 1. (Müller-Hermes's Section 4 also reduces decomposability of tensor products of
Werner–Holevo maps to a linear program, by the same twirling argument as here.) So from n = 3 on, 𝒫_{3/4}^{⊗n} is not
decomposable, and a proof of tensor-stable
positivity at any a < 1 cannot go through decomposability of the tensor powers. Nor can it go through the step of
Proposition M(4), which gives positivity of 𝒫_a^{⊗n} from 2-positivity of 𝒫_a^{⊗(n−1)}. The program with a qubit
ancilla shows that 𝒫_a^{⊗k} ⊗ id₂ is decomposable only at a = 1, for every k ≤ 6. So that 2-positivity is never
decomposable either.

The same thresholds rule out a natural induction. Suppose 𝒫_a^{⊗(k+1)} were a sum of maps Φ ∘ (𝒫_a^{⊗j} ⊗ ℰ) ∘ Ψ with
j ≤ k, Φ and Ψ completely positive and ℰ entanglement breaking. Then every tensor power would be such a sum (by
induction on the number of copies), and so 𝒫_a would be tensor-stable positive. For k ≤ 2 this cannot happen at any
a < 1. The building blocks 𝒫_a^{⊗j}, j ≤ 2, are decomposable for a ≥ 0.6703, these operations keep decomposability,
and so every tensor power would be decomposable, against the table. (For a < 0.6703 the map is not even
tensor-stable positive, by Theorem 7.1.) A certificate of this inductive kind would need a building block with three
or more copies, where the tensor powers stop being decomposable.

On a plane the picture is different. The same program for 𝒬_a, with G₂ = [[½, 1, −½], [½, 0, ½], [−½, 1, ½]] and
c_C = (a − ¾, a − ¼, a + ¼)^{⊗n}, gives decomposability thresholds equal to a_n to eight digits for every n ≤ 12.
Together with Theorem 7.1 this says that, for n ≤ 12, 𝒬_a^{⊗n} is positive exactly when a ≥ a_n (computer-assisted,
floating point). Theorem 7.18 proves it for every n, with an explicit decomposition.

**Observation 7.4 (numerical).** For n = 2, …, 6 the positivity threshold of 𝒫_a^{⊗n} found by search is a_n. At
a_n − 0.001 the search finds negative values, and at a_n + 0.001 it finds none. The searches use alternating minimisation
seeded at the vector of Theorem 7.1, at random rotations and perturbations of it in all three directions of every copy,
at random starts, and by continuation in a. Every negative minimiser found has one-copy rank 2, so it lives on one plane
per copy. Bisections without these seeds agree to five digits for n = 2, 3. We also tried inputs with full support at
every copy: the qutrit GHZ state, ε_{ijk}, the AME(4,3) state, the GHZ state in the Fourier basis with its real and
imaginary parts, and tensor products of these with each other and with the witness of Theorem 7.1. For n = 3, …, 6 their
least output eigenvalues stay well above zero at a_n + 0.001 and at 3/4. Descents started from them (n ≤ 5) end at the
product value or on the plane branch. For larger n we searched over translation-invariant real matrix product
states x and y of bond dimension 2 and 3 (n = 9, 15, 21). At a_n − 0.003 the least value found is exactly that of
the witness of Theorem 7.1. At a_n + 0.003 and at 3/4 nothing negative appears (the least values are the product value
at n = 9, 15, and positive values of order 10⁻¹¹ at n = 21, where the optimiser stops short of it). We also searched over
general (not translation-invariant) open-boundary real matrix product states x and y of bond dimension 4 and 6, by
alternating one-site sweeps: for fixed y, 𝒫_a^{⊗n}(yyᵀ) is a matrix product operator, x is updated site by site as
its lowest eigenvector, and then the roles of x and y are exchanged. For n = 6, 8, 10, 12, at a_n + 0.002 and at
3/4, with 19 starts each (random, random with full-support marginals, and the witness), the least value found in
every run is the product value (a − ½)^n. We also tested n → ∞ for translation-invariant inputs. For periodic matrix
product states x and y the value is Tr(𝔼ⁿ) for a fixed transfer operator 𝔼. If 𝒫_a^{⊗n} is positive for every n, the
series Σ_n Tr(𝔼ⁿ)zⁿ has nonnegative coefficients, so by Pringsheim's theorem the spectral radius of 𝔼 is an eigenvalue.
Maximising, over bond dimensions 2, 3 and 4 at a = ¾, the ratio of the largest modulus among the other eigenvalues to
the largest positive eigenvalue gives at most 1 (reached by the plane witness and by states of period two), never more.
At a = 0.74 the plane witness gives 1/(2a − ½) ≈ 1.020, and the search finds it (`code/ti_transfer_search.py`,
`logs/ti_transfer_*.log`).

**Conjecture QT.** For every n and every a ≥ 0, 𝒫_a^{⊗n} is positive if and only if a ≥ a_n. (For negative a the
restriction matters: when a ≤ −1/4, −𝒫_a is completely positive, so 𝒫_a^{⊗n} is positive for every even n.)

Since a_n < 3/4 for every n, Conjecture QT makes 𝒫_{3/4} tensor-stable positive. By Proposition M(3), the Werner state
proportional to 11I − 5F would then be NPT and undistillable from every number of copies. In words, the conjecture says
that the only obstruction is the one on planes, where 𝒫_{3/4} is completely positive. Its case n = 2 is Theorem 7.2.
The threshold is open for every n ≥ 3 (at n = 3 it lies in [a₃, 0.648], Theorem 7.13). At a = 3/4, where by
Proposition 7.3 the positivity it asserts is not decomposable from n = 3 on, Theorem 7.13 proves it for n ≤ 5 and
Theorem 7.32 for n = 6, and Corollaries 7.9 and 7.14 reduce the rest to seven or more copies at which both the input and
the output vector have full support.

**Remark 7.5 (an exact sum over planes, and why it is not enough).** As maps on M₃(C),

    𝒫_a = Σ_π ι_π ∘ 𝒬′_a ∘ c_π,   𝒬′_a(X) = a Tr(X) I₂ − ¼(X + Xᵀ) − ½(a − ½) diag(X₁₁, X₂₂),

where π runs over the three coordinate planes, c_π compresses to π and ι_π embeds back. (Both sides vanish on
antisymmetric matrices; on real rank-one inputs and outputs the identity reads
a|x|²|y|² − ½(x·y)² = Σ_π [a|x_π|²|y_π|² − ½(x_π·y_π)²] − (a − ½)Σ_k x_k²y_k².) The identity holds copy by copy, so 𝒫_a^{⊗n} is the sum, over the 3^n assignments of a
plane to each copy, of compressions of 𝒬′_a^{⊗n}. It cannot prove Conjecture QT term by term. If 𝒬′_a were completely
positive, so would be 𝒫_a; indeed the Choi matrix of 𝒬′_a has least eigenvalue −(1 − a)/2. So for a < 1 the qubit map
𝒬′_a is neither completely positive nor completely copositive, and its tensor powers are eventually not positive
[MRW, Thm. 2]. A proof has to use cancellation between different plane assignments.

**Remark 7.6 (a form in terms of collisions).** The six lines through opposite vertices of an icosahedron form a
projective 4-design, with Σ_j v_jv_jᵀ = 2I. Hence 𝒫_a(X) = Σ_{j=1}^6 Tr(E_jX) Q_j with E_j = v_jv_jᵀ/2 and
Q_j = (a + ¼)I − (5/4)v_jv_jᵀ, and for any states ρ, σ on n copies

    ⟨σ, 𝒫_a^{⊗n}(ρ)⟩ = (a + ¼)^n E[(−r)^N],   r = (9/4 − a)/(a + ¼),

where N is the number of copies at which two independent outcome strings agree, one from measuring ρ and one from
measuring σ with the product measurement {E_j}^{⊗n}. Positivity of 𝒫_a^{⊗n} says that this generating function of the
collision count is nonnegative at −r. The value a = 1 is r = 1, a = 3/4 is r = 3/2, and a = 1/2 is r = 7/3.

**Lemma 7.7 (inputs on a plane are entanglement breaking for a ≥ 3/4).** Let V: C² → C³ be an isometry onto any
two-dimensional subspace, real or not, and Φ_V(Y) = 𝒫_a(VYV*). Its Choi matrix is

    C_V = aI − ¼(|ω_V⟩⟨ω_V| + S_V),   ω_V = Σ_i e_i ⊗ Ve_i,   S_V = Σ_{ij} |e_i⟩⟨e_j| ⊗ |V̄e_j⟩⟨V̄e_i|.

The partial transpose of the output exchanges the two terms, so C_V^Γ = C_V. Since ‖ω_V‖² = 2 and S_V is a partial
isometry, C_V ⪰ (a − ¾)I. So for a ≥ 3/4 the matrix C_V is positive with positive partial transpose, hence separable
(it lives on C² ⊗ C³), and Φ_V is entanglement breaking. For a < 3/4 and a real plane, C_V has the eigenvalue a − ¾. ∎

**Proposition 7.8 (thin copies are free).** Let a ≥ 3/4, let ρ be a state on (C³)^{⊗n}, and let K be the set of copies
at which the one-copy reduced state of ρ has rank 3. Write μ_k(a) for the least eigenvalue of 𝒫_a^{⊗k}(σ) over states σ
on k copies (μ₀ = 1). If μ_{|K|}(a) ≥ 0, then

    λ_min(𝒫_a^{⊗n}(ρ)) ≥ (a − ½)^{n−|K|} μ_{|K|}(a) ≥ 0.

By Theorem 7.2 (whose proof gives μ₂(a) = min{(a − ½)², (a − ¼)² − ⅛} for a ≥ (1 + √2)/4) this holds for every n
whenever |K| ≤ 2. At a = 3/4 it gives λ_min ≥ 4^{−n} for every input thin at all but two copies, which is the value of
real product inputs: on this class the least output eigenvalue is exactly multiplicative.

*Proof.* For s ∉ K let V_s be an isometry onto a plane that contains the support of ρ_s. Then ρ = Wρ′W* with
W = ⊗_{s∉K} V_s ⊗ I_K and ρ′ ⪰ 0, so 𝒫_a^{⊗n}(ρ) = (T ⊗ 𝒫_a^{⊗|K|})(ρ′) with T = ⊗_{s∉K} Φ_{V_s}. T is entanglement breaking
(Lemma 7.7), so (T ⊗ id)(ρ′) = Σ_i σ_i ⊗ ρ_i with σ_i, ρ_i ⪰ 0, and with Ψ = 𝒫_a^{⊗|K|},

    (T ⊗ Ψ)(ρ′) = Σ_i σ_i ⊗ Ψ(ρ_i) ⪰ μ_{|K|}(a) Σ_i Tr(ρ_i) σ_i ⊗ I = μ_{|K|}(a) T(Tr_K ρ′) ⊗ I ⪰ μ_{|K|}(a) λ_min[T] I.

This is the argument of [MRW, Thm. 9], which only needs Ψ(σ) ⪰ μ Tr(σ) I with μ ≥ 0. Finally
λ_min[T] = (a − ½)^{n−|K|}: each Φ_{V_s} has least output eigenvalue a − ½ over states (¼(VYV* + (VYV*)ᵀ) has norm at most
½ for a state Y, with equality on real ones), and the least output eigenvalue over states is multiplicative for
entanglement breaking maps [MRW, Thm. 9]. ∎

**Corollary 7.9.** 𝒫_{3/4} is tensor-stable positive, and then NPT bound entanglement exists, if and only if for every n
the operator 𝒫_{3/4}^{⊗n}(ρ) is positive for every state ρ whose one-copy reduced states all have rank 3. A smallest
failure, if there is one, needs n ≥ 3 copies and an input that is thick (rank 3) at every copy (n ≥ 7 by
Corollary 7.14).

*Proof.* Induction on n. An input with a thin copy is covered by Proposition 7.8 and the case of fewer copies; n ≤ 2 is
Theorem 7.2. ∎

The output vector can be thinned as well. 𝒫_a is self-adjoint for the Hilbert–Schmidt inner product, so
⟨y|𝒫_a^{⊗n}(xx*)|y⟩ = ⟨x|𝒫_a^{⊗n}(yy*)|x⟩. Let x be thin outside K_x and y thin outside K_y, with a ≥ 3/4. Then
⟨y|𝒫_a^{⊗n}(xx*)|y⟩ ≥ 0 as soon as 𝒫_a^{⊗|K_x ∩ K_y|} is positive. To see this, decompose as in the proof of Proposition 7.8
on the copies outside K_x: the value becomes Σ_i Tr[τ_i 𝒫_a^{⊗|K_x|}(ρ_i)] = Σ_i Tr[𝒫_a^{⊗|K_x|}(τ_i) ρ_i], where
τ_i = Tr_{∉K_x}[(σ_i ⊗ I) yy*] ⪰ 0. Each τ_i is thin wherever y is, so Proposition 7.8 makes 𝒫_a^{⊗|K_x|}(τ_i) ⪰ 0. So a
smallest failure needs both the input and the output vector thick at every copy.

**Corollary 7.10 (a statement about minimisers).** Let f_n(x) = λ_min(𝒫_{3/4}^{⊗n}(xxᵀ)) on real unit vectors x. Then
𝒫_{3/4} is tensor-stable positive if and only if, for every n ≥ 3, f_n attains its minimum at some input that is thin
(one-copy marginal of rank at most 2) at one copy or more.

*Proof.* If the minimum is attained at an input thin at some copy, Proposition 7.8 and positivity of fewer copies (by
induction) make that minimum nonnegative. Conversely, tensor-stable positivity gives min f_n = (1/4)^n (Proposition M(4)),
which the product vector e₁^{⊗n} attains. ∎

So it would be enough to show that no input with full support at every copy is a global minimiser. Numerically this is
what happens. At n = 3, inputs whose one-copy marginals all have least eigenvalue at least τ give values no lower than
0.116, 0.087, 0.058 and 0.045 for τ = 0.25, 0.15, 0.08 and 0.04 (projected gradient search; the product value is
1/64 ≈ 0.016), and the constraint is active at every optimum found.

**Proposition 7.11 (2-positivity of tensor powers above 3/4).**
(i) For every a ≥ a* = (23 + 6√35)/68 ≈ 0.8602 (the larger root of 272a² − 184a − 43), 𝒫_a ⊗ 𝒫_a is 2-positive. For
a ∈ [a*, 1] there is the identity

    C_a ⊗ C_a = λ(a) [M ⊗ S + S ⊗ M] + Σ_{j,k} w_jk(a) P_j ⊗ P_k,    λ(a) = 3(1 − a)(4a + 1),

with S = P_s0/10 + P_an/6, the O(3) twirl of |e₁e₂⟩⟨e₁e₂| and so a separable operator, and

    w_ΩΩ = (1 − a)²,  w_Ωs = w_sΩ = (2/5)(1 − a)²,  w_Ωa = w_aΩ = 0,  w_ss = (272a² − 184a − 43)/80,
    w_sa = w_as = 3(4a + 1)(28a − 23)/80,  w_aa = 5(4a + 1)(4a − 3)/16,

all nonnegative on [a*, 1].
(ii) 𝒫_a^{⊗n} is 2-positive for a ≥ 91/100 (n = 3), 93/100 (n = 4) and 19/20 (n = 5). At those values, C_a^{⊗n} minus
a nonnegative rational combination of operators M ⊗ S₁ ⊗ … ⊗ S_{n−1} (M at one copy, twirls of product states at the
others) is a nonnegative combination of the products P_J (exact certificates with 10, 28 and 94 terms).

*Proof.* M is positive on Schmidt-rank-two vectors (the endpoint at one copy: |Tr Z|² ≤ 2‖Z‖² for rank Z ≤ 2). So is
M ⊗ σ for separable σ = Σ|α⟩⟨α| ⊗ |β⟩⟨β|, because
contracting a Schmidt-rank-two vector with the product vector α ⊗ β leaves Schmidt rank at most two. The remainders
are positive semidefinite. Above the stated values write C_a = C_{a₀} + (a − a₀)I and expand: each term is
C_{a₀}^{⊗K} ⊗ I, which is covered by the certificates for |K| copies (their thresholds increase with the number of
copies) and the separability of I. ∎

As far as we know these are the first proofs of 2-positivity of tensor powers of 𝒫_a for some a < 1. By
Proposition 7.3, 𝒫_a^{⊗k} ⊗ id₂ is not decomposable for any a < 1, so decomposability cannot give them. They do not reach
a = 3/4. The best thresholds of this kind of certificate are 0.8602, 0.9031, 0.9259 and 0.9400 for n = 2, …, 5
(floating-point linear programming), and they increase with n. Through Proposition M(4) they give positivity of
𝒫_a^{⊗(n+1)}, but less than decomposability already gives there.

**Remark 7.12 (the two-copy endpoint and the map are two ends of one family).** Put X_w = M + (w − 1)P_an, so that
X₁ = M and X₂ = 2C_{3/4}. For a fixed vector, ⟨χ|X_{w₁} ⊗ X_{w₂}|χ⟩ is affine in w₁ and in w₂. So positivity on
Schmidt-rank-two vectors for all (w₁, w₂) ∈ [1, 2]² reduces to the corners (1, 1), (1, 2) and (2, 2). The corner (1, 1)
is the two-copy Werner endpoint, settled in July 2026 (references in `copy-amplification/README.md`). The corner (2, 2)
is 2-positivity of 𝒫_{3/4}^{⊗2}, which would give positivity of 𝒫_{3/4}^{⊗3} by Proposition M(4). Numerically, X_w ⊗ X_w
is positive on Schmidt-rank-two vectors exactly for w ∈ [1, 2]: the minimum is −1.2·10⁻² at w = 0.95, zero to rounding
at w = 1, 1.5 and 2, and −1.3·10⁻³ at w = 2.01. At w = 2 the zeros found form two families, and every vector in them has
a thin copy. Thin copies are covered by a proved statement. If one of the four one-copy marginals (A₁, B₁, A₂ or B₂) of
a Schmidt-rank-two χ has rank at most 2, then ⟨χ|X_{w₁} ⊗ X_{w₂}|χ⟩ ≥ 0 for all w₁, w₂ ∈ [1, 2]. Indeed, the
restriction of X_w to a plane on the input side is positive with positive partial transpose, by Lemma 7.7 for X₂ and
directly for X₁, and these properties pass to convex combinations. So the restriction is separable on C² ⊗ C³. The
output side follows from FX_wF = X_w, which makes the restriction to C³ ⊗ (plane) the swap of the restriction to
(plane) ⊗ C³. One then contracts χ with the product vectors of the separable decomposition and uses one-copy
positivity of X_{w₂} on Schmidt-rank-two vectors.

**Theorem 7.13 (three, four and five copies at a = 3/4; computer-assisted).** For n = 3, 4, 5,

    𝒫_{3/4}^{⊗n}(ρ) ⪰ 4^{−n} Tr(ρ) I   for every positive ρ on (C³)^{⊗n},

and the least eigenvalue equals 4^{−n} at ρ = (e₁e₁ᵀ)^{⊗n}. So 𝒫_{3/4}^{⊗n} is positive and μ_n(3/4) = 4^{−n} for every
n ≤ 5. Moreover, 𝒫_a^{⊗3} is positive for every a ≥ 81/125 = 0.648, so with Theorem 7.1 the three-copy threshold
lies in [a₃, 0.648] = [0.6468503…, 0.648].

*Proof.* Let B′ be a second copy of the outputs B = B₁…B_n, S the swap of B and B′, P = (I + S)/2, and Γ_A, Γ_B the
partial transposes on all inputs and on all outputs. For each case, `code/dps/certs/` contains rational operators
Q_A, Q_B ⪰ 0 with PQ_AP = Q_A such that

    Z = P [ C_a^{⊗n} ⊗ I_{B′} − t I − Γ_A(Q_A) − Γ_B(Q_B) ] P ⪰ 0,

with t = 4^{−n} at a = 3/4 and t = 0 at a = 81/125, n = 3. For ψ = x ⊗ y ⊗ y we have Pψ = ψ, so

    ⟨xy|C_a^{⊗n}|xy⟩ |y|² − t |x|² |y|⁴ = ⟨ψ|Z|ψ⟩ + ⟨x̄yy|Q_A|x̄yy⟩ + ⟨xȳy|Q_B|xȳy⟩ ≥ 0,

and ⟨xy|C_a^{⊗n}|xy⟩ = ⟨y|𝒫_a^{⊗n}(|x̄⟩⟨x̄|)|y⟩. This is the second level of the Doherty–Parrilo–Spedalieri hierarchy,
with the extension on the output side. The operators are rational combinations of products of the fifteen Brauer
diagrams of one copy (on its three strands A_s, B_s, B′_s), so they commute with SO(3) acting on each copy, and
positivity reduces to exact rational matrices of size at most 3^n, one for each assignment of spins j_s ∈ {0, 1, 2, 3}.
At t = 4^{−n} the real product vectors force one kernel vector in each block with all j_s ∈ {1, 3}, and the certificates
have no other kernel. Positivity of every block was decided in exact arithmetic by two separately written verifiers
(`code/dps/README.md`).
For a > 81/125 write 𝒫_a = 𝒫_{81/125} + (a − 81/125)D with D(X) = Tr(X)I and expand 𝒫_a^{⊗3}. Up to the order of the
copies, each term is 𝒫_{81/125}^{⊗k} ⊗ D^{⊗(3−k)}. It is positive because 𝒫_{81/125}^{⊗k} is positive for k = 3 (the
certificate), k = 2 (Theorem 7.2, as 81/125 > (1 + √2)/4) and k = 1, and (Φ ⊗ D)(X) = Φ(Tr₂X) ⊗ I. ∎

**Corollary 7.14.** For every n, 𝒫_{3/4}^{⊗n}(ρ) ⪰ 4^{−n} Tr(ρ) I for every state ρ whose one-copy reduced states have
rank 3 at five copies or fewer, and 𝒫_{3/4}^{⊗n}(ρ) ⪰ 4^{−n−1} Tr(ρ) I if they have rank 3 at exactly six copies. Also
⟨y|𝒫_{3/4}^{⊗n}(xx*)|y⟩ ≥ 0 whenever x and y are thick (one-copy marginal of rank 3) at six or fewer common copies. A
smallest failure of tensor-stable positivity of 𝒫_{3/4}, if there is one, needs n ≥ 7 copies and an input and an
output vector that are both thick at seven or more common copies.

*Proof.* Proposition 7.8 with μ_k(3/4) = 4^{−k} for k ≤ 5 (Theorem 7.2 for k ≤ 2, Theorem 7.13 for k = 3, 4, 5) and
μ₆(3/4) ≥ 4^{−7} (Theorem 7.32), and the thinning of the output vector after Corollary 7.9. ∎

**Remark 7.15 (what the certificates suggest).** At a = 3/4 the second level of the hierarchy is exact for n ≤ 5: its
optimal t is the true minimum 4^{−n}. The first level (decomposability, the linear program of Proposition 7.3 with the
constant t added) is exact at n = 1, 2 (t = 1/4, 1/16) but reaches only t ≈ −0.0036 at n = 3. The optimal dual point
of that linear program is an explicit state. Let ρ = Σ_J p_J P_J / Tr P_J, where J runs over the assignments of P_Ω, P_s0,
P_an to the three pairs (input, output) of copies. The weights are 1/27, 10/117, 4/13 and 200/351 on the letter
multisets {Ω, Ω, Ω}, {Ω, Ω, s0}, {Ω, an, an} and {s0, s0, s0}, spread evenly over arrangements. Then ρ has positive
partial transpose, and Tr[C_{3/4}^{⊗3}ρ] = −3/832 (exact, `code/ppt_state_n3.py`). By Theorem 7.13, ρ is therefore an
explicit PPT entangled state on (C³)^{⊗3} ⊗ (C³)^{⊗3} that C_{3/4}^{⊗3} detects. At one copy the second level has a
closed form: with Y = ¼(I − F_{AB}) ⊗ I_{B′} one has
C_{3/4} ⊗ I − ¼I = Y + Γ_A(Y) and PYP = ⅜(P − Π), where Π projects onto Sym³(C³) ⊂ A ⊗ B ⊗ B′. We have no closed form
for n ≥ 2, and an argument for every n would need one. Two things block a direct induction. The cone
PSD + Γ_A(PSD) + Γ_B(PSD) is not closed under tensor products: a product such as Z ⊗ Γ_A(Q′) is in none of the three
parts. And the symmetric subspace of (B₁B₂) ⊗ (B′₁B′₂) is not the product of the two one-copy symmetric subspaces: it
also contains antisymmetric ⊗ antisymmetric. At n = 2 the second level is exact for every a ≥ a₂ (its value equals
μ₂(a) = min{(a − ½)², (a − ¼)² − ⅛} at every a tested). At n = 3 it is not exact near the threshold: it certifies
positivity only from a ≈ 0.64798 on, above a₃ ≈ 0.646850, and it is exact again from a ≈ 0.66 on. At n = 4 it needs
the partial transpose on the outputs: with Γ_A alone the best t is 0.0029257 < 1/256. At n = 3 either partial transpose
alone suffices, and with neither (the symmetric extension alone) the best t is −0.175.

Let ℓ_n be the least a at which the second level certifies positivity of 𝒫_a^{⊗n}. Then ℓ_{n+1} ≥ ℓ_n. Indeed, a
feasible point of the dual program at n copies (a state ρ on A ⊗ B ⊗ B′ with PρP = ρ, Γ_A(ρ) ⪰ 0 and Γ_B(ρ) ⪰ 0),
tensored with the state of the real product vector e₁ ⊗ e₁ ⊗ e₁ at one more copy, is feasible at n + 1 copies, and
its value is multiplied by a − ½ > 0. The primal program is strictly feasible, so a negative optimum at n copies has a
dual point of negative value, and this point stays negative at n + 1. Numerically ℓ₂ = a₂, ℓ₃ ≈ 0.6480, ℓ₄ ≈ 0.686 and
ℓ₅ ≈ 0.720, while a₄ ≈ 0.6704 and a₅ ≈ 0.6853 (at n = 5 and a = 0.72 the best t is −5.8·10⁻⁵). The increments 0.044,
0.038 and 0.034 shrink slowly. At six copies the second level still certifies positivity at a = 3/4, so ℓ₆ < 3/4
(Theorem 7.32), but in floating point it is no longer exact there: the solver's optimum is t ≈ 1.3749·10⁻⁴ ≈
0.563·4^{−6}, below the value 4^{−6} of real products. So at a = 3/4 the method appears to be close to its limit. These numbers are floating point
(`code/dps/logs/sdp/threshold*.log`, `code/dps6/logs/`).

**Proposition 7.16 (permutation-symmetric inputs and outputs).** Let Sym^n ⊂ (C³)^{⊗n} be the symmetric subspace.
(i) For every n and every a < a_n there are x, y ∈ Sym^n with ⟨y|𝒫_a^{⊗n}(xx*)|y⟩ < 0.
(ii) For n ≤ 21 and all x, y ∈ Sym^n,

    ⟨y| 𝒫_{3/4}^{⊗n}(xx*) |y⟩ ≥ 4^{−n} |x|² |y|²   (computer-assisted).

For real x, y, equality holds only when x and y are both multiples of v^{⊗n} for one real unit vector v. Hence, for
n ≤ 21, the least value on Sym^n is (a − ½)^n for every a ≥ 3/4.

*Proof.* (i) The vectors of Theorem 7.1 are symmetric. Since Yp = ip and Yp̄ = −ip̄, (Y ⊗ I)x = i(p^{⊗n} − p̄^{⊗n})/√2.
(ii) Real vectors suffice, as before. Write C_{3/4} = I/4 + A/2 + B/2 on one pair, with A = (I − F)/2 and
B = A^Γ = (I − |Ω⟩⟨Ω|)/2. Expanding C_{3/4}^{⊗n} gives words in I, A, B. Compressed to Sym^n ⊗ Sym^n, a word depends only on
its letter counts (k, s, t), and the partial transpose Γ on the outputs exchanges s and t. Let Q_n be the compression of
the words with s > t, plus half of those with s = t ≥ 1, with their weights. Then
Q_n + Γ(Q_n) = (Π ⊗ Π)C_{3/4}^{⊗n}(Π ⊗ Π) − 4^{−n}Π ⊗ Π, so for real symmetric x, y the difference
⟨y|𝒫_{3/4}^{⊗n}(xxᵀ)|y⟩ − 4^{−n}|x|²|y|² equals 2⟨x ⊗ y|Q_n|x ⊗ y⟩. For n ≤ 8, Q_n is positive semidefinite with kernel
exactly Sym^{2n}(R³). For 9 ≤ n ≤ 21 the same holds after the split is changed on the words with at most two letters I
(explicit rational values in `code/symmetric_sector/many_n/certificates.txt`); any split with weights y and 1 − y on
the classes (k, s, t) and (k, t, s) keeps the identity. Positivity and the kernel were decided in exact arithmetic,
block by block: for n ≤ 9 by two separately written checks, and for 10 ≤ n ≤ 21 by one check in a coupled
angular-momentum basis (`code/symmetric_sector/`). If x ⊗ y lies in Sym^{2n}, then x and y are multiples of one v^{⊗n}. For a > 3/4,
expand 𝒫_a = 𝒫_{3/4} + (a − 3/4)D as in Theorem 7.13; the partial traces of symmetric vectors are supported on symmetric
subspaces. ∎

The inequality in (ii) is new on this sector for 6 ≤ n ≤ 21 (for n ≤ 5 it follows from Theorem 7.13). On the full
space no decomposition of this kind exists from n = 3 on (Proposition 7.3); on Sym^n one exists for every n ≤ 21. The
greedy split fails from n = 9 on, in a single spin block. The repaired splits are not canonical, and no rule valid for
every n is known. The least eigenvalue of 4^nQ_n off Sym^{2n} is about 1 for these repairs. If every split is left
free and this margin is maximised, it increases with n (floating point): 1.55 at n = 4, 1.805 at n = 12, 1.825 at
n = 20, with shrinking increments, apparently toward a limit near 1.83–1.85. No certificate of this kind can push it
above 2: for a unit vector v and c ⊥ v, the
vector ψ = Sym(v^{n−1}c) ⊗ vⁿ − vⁿ ⊗ Sym(v^{n−1}c) has ⟨ψ|4^nQ_n|ψ⟩ = 2|ψ|² for every certificate. This follows by
expanding the product value to second order at real products. The best constant of these decompositions equals the
least value of Tr(K_nρ), K_n = (Π ⊗ Π)C_{3/4}^{⊗n}(Π ⊗ Π), over states ρ on Sym^n ⊗ Sym^n with positive partial
transpose.

Numerically, on Sym^n the least value is min{witness value, product value} for n ≤ 7 at a_n ± 0.002 and at 3/4. The same
holds with x symmetric and y arbitrary for n ≤ 6. So Conjecture QT appears to hold on this sector. Decompositions of
this kind cannot prove it near the threshold, though. On Sym^n ⊗ Sym^n at a = a_n the best constant is negative, about
−0.011, −0.014, −0.005 and −0.002 for n = 2, …, 5. At n = 2 the threshold is nevertheless a₂, by the singular value
argument of Theorem 7.2.

**Proposition 7.17 (inputs that thin out one copy at a time, every n).** Let a ≥ 3/4 and m ≤ 5. Call a pair (x, y) of
vectors on n copies *reducible* if n ≤ m, or if one of the two vectors, say x, has a one-copy reduced state of rank at
most 2 at some copy s and the pair (⟨e|_s x, ⟨f|_s y) on the other n − 1 copies is reducible for all e, f ∈ C³ (and the
same with the roles of x and y exchanged). Here ⟨e|_s x is the contraction of x with e at copy s. For every reducible
pair,

    ⟨y| 𝒫_a^{⊗n}(xx*) |y⟩ ≥ (a − ½)^n |x|² |y|².

With m ≤ 2 the proof uses only Theorem 7.2, and with m ≤ 5 it uses Theorem 7.13. In particular, if for some ordering of
the copies x has Schmidt rank at most 2 across each cut {1, …, k} | {k + 1, …, n} with k ≤ n − m, then
𝒫_a^{⊗n}(xx*) ⪰ (a − ½)^n |x|² I. This holds for every open matrix product state of bond dimension two, at every n. An
example is the open AKLT chain. Its one-copy reduced states have full rank at every copy except the two ends, so for
n ≥ 8 copies (more than five thick copies) Corollary 7.14 does not reach it.

*Proof.* Induction on n. Let x be thin at copy s, let V be an isometry from C² onto a plane that contains the support of
its reduced state there, and write x = (V ⊗ I)x′. By Lemma 7.7, Φ_V(Y) = 𝒫_a(VYV*) is entanglement breaking. So its
Choi matrix is a finite sum of rank-one products, and Φ_V(Y) = Σ_i p_i ⟨α_i|Y|α_i⟩ |β_i⟩⟨β_i| with p_i > 0. Applying
this at copy s,

    ⟨y| 𝒫_a^{⊗n}(xx*) |y⟩ = Σ_i p_i ⟨y_i| 𝒫_a^{⊗(n−1)}(x_i x_i*) |y_i⟩,   x_i = ⟨Vα_i|_s x,   y_i = ⟨β_i|_s y.

Each pair (x_i, y_i) is reducible, so by induction the right side is at least (a − ½)^{n−1} Σ_i p_i |x_i|² |y_i|². This
sum equals Tr[σ Φ_V(ρ′)], where ρ′ and σ are the reduced states of x′ and y at copy s. It is at least (a − ½)|x|²|y|²,
because Φ_V(Y) ⪰ (a − ½) Tr(Y) I for Y ⪰ 0 (proof of Proposition 7.8). If y is the thin vector, use
⟨y|𝒫_a^{⊗n}(xx*)|y⟩ = ⟨x|𝒫_a^{⊗n}(yy*)|x⟩. For n ≤ m the bound is Theorem 7.2 (whose proof gives μ₂(a) = (a − ½)² for
a ≥ 3/4) or Theorem 7.13, extended to a > 3/4 by the expansion 𝒫_a = 𝒫_{3/4} + (a − 3/4)D used there. For the
statement about cuts: the Schmidt rank across {1} | {2, …, n} bounds the rank of the reduced state at copy 1, and
contracting copy 1 with a vector does not raise the Schmidt rank across the later cuts. ∎

This is the argument of [MRW, Thm. 9] (used in Proposition 7.8) applied once per copy, with the conditional states after
each step pure, so that thinness only has to hold after the earlier copies have been contracted.
A failure of tensor-stable positivity at a = 3/4 needs a pair that is not
reducible. In particular, in every ordering of the copies, both x and y must have Schmidt rank at least three across
some cut {1, …, k} | {k + 1, …, n} with k ≤ n − 5. `code/bond_two_check.py` (floating point, `logs/bond_two_check.log`)
illustrates the bound for the open AKLT chain and random open matrix product states of bond dimension two, n = 4, …, 8.
For the AKLT chain the least output eigenvalue is 20 to 800 times the product value.

**Theorem 7.18 (the threshold on a real plane, every n).** Let 𝒬_a(Y) = a Tr(Y) I₂ − ¼(Y + Yᵀ) be the compression of 𝒫_a
to a real plane. For every n ≥ 1 and every a ≥ 0 the following are equivalent: (i) 𝒬_a^{⊗n} is positive;
(ii) 𝒬_a^{⊗n} is decomposable; (iii) a ≥ a_n. In the language of [FM]: for every real t, the tensor power Υ_t^{⊗n} of
the Pauli map with eigenvalues 1, t, 0, t on I, X, Y, Z is positive if and only if |t| ≤ 2^{−1+1/n}, and then it is
decomposable. So Conjecture QT holds for inputs and outputs that lie in real planes at every copy (the planes of the
input and of the output may differ; see after the proof).

*Proof.* On the plane 𝒬_a has the Pauli eigenvalues 2a − ½, −½, 0, −½, so 𝒬_a = (2a − ½) Ad_Y ∘ Υ_t with t = 1/(4a − 1),
and a ≥ a_n means u := 2t ≤ 2^{1/n}. If a < a_n, Theorem 7.1 applies, and (ii) implies (i). For (iii) ⇒ (ii), first
let u = 2^{1/n}. The Choi matrix of Υ_t is Bell-diagonal with eigenvalues (1 + u, 1, 1, 1 − u)/4 on Φ⁺, Φ⁻, Ψ⁺, Ψ⁻.
The partial transpose maps a Bell-diagonal operator to a Bell-diagonal one, and acts on its eigenvalues by
γ = ½[[1, 1, 1, −1], [1, 1, −1, 1], [1, −1, 1, 1], [−1, 1, 1, 1]]. For a string J of Bell labels let F(J), F̃(J) and
R(J) be the products over the copies of the entries of (1 + u, 1, 1, 1 − u), (1 + u, 1, 1, u − 1) and (2/u, 1, 1, 0).
Directly, γ fixes (1 + u, 1, 1, 1 − u), maps (1 + u, 1, 1, u − 1) to u·(2/u, 1, 1, 0), and maps (2/u, 1, 1, 0) to
(1 + u, 1, 1, u − 1)/u. Let P be the Bell-diagonal operator with eigenvalues 4^{−n}p, where p = ½(F + F̃) − R. The
partial transpose acts copy by copy, so P^Γ has eigenvalues 4^{−n}[½(F + uⁿR) − u^{−n}F̃], and

    P + P^Γ has eigenvalues 4^{−n}[F + (½uⁿ − 1)R + (½ − u^{−n})F̃] = 4^{−n}F,

which are those of the Choi matrix of Υ_t^{⊗n}. Moreover p ≥ 0. With k₀ and k₃ the numbers of Φ⁺ and Ψ⁻ in J,
p = (1 + u)^{k₀}(1 − u)^{k₃} when k₃ ≥ 2 is even, p = 0 when k₃ is odd, and p = (1 + u)^{k₀} − (2/u)^{k₀} ≥ 0 when
k₃ = 0, because u ≥ 1. So Υ_t^{⊗n} is decomposable at u = 2^{1/n}. For smaller t ≥ 0, Υ_t is a convex combination of
Υ_{t_n} (t_n = 2^{−1+1/n}) and the completely depolarizing map Υ₀, which is entanglement breaking. Expanding the tensor
power gives terms Υ_{t_n}^{⊗k} ⊗ (entanglement breaking), k ≤ n, and Υ_{t_n}^{⊗k} is decomposable by induction on n
(as t_n ≤ t_k). Decomposability survives tensoring with entanglement-breaking maps and composing with Ad_Y. (The
restriction a ≥ 0 is needed: for even n and a ≤ (1 − 2^{1−1/n})/4 < 0 the tensor power is positive although a < a_n.)
If the input lies in a real plane V_x and the output in a real plane V_y at some copy, the compressed map there is
Y ↦ V_yᵀ𝒫_a(V_xYV_xᵀ)V_y = K𝒬_a(Y)Kᵀ + a Tr(Y)(I − KKᵀ) with K = V_yᵀV_x, a contraction. Expanding the tensor power,
every term is a completely positive map applied after 𝒬_a^{⊗k} (k ≤ n), tensored with entanglement-breaking maps,
so it is positive for a ≥ a_n as well. ∎

For n = 2 the threshold is [FM, Example 3, via their Proposition 4], and the decomposition of Υ ⊗ Υ at
t = 1/√2 is already [FM, Example 1, eq. (18)]. At n = 3 [FM] check sufficiency only numerically (their Section VIII),
although their abstract calls that case fully characterized. At n = 2 the equivalence of (i) and (ii) is also a case
of Müller-Hermes's theorem that every positive tensor square of a qubit map is decomposable [A. Müller-Hermes,
J. Math. Phys. 62, 092204 (2021), arXiv:2006.14543]. He leaves open (his Section 9) whether a positive qubit map can
have a positive but non-decomposable tensor cube, and reduces that question to Pauli maps. Theorem 7.18 answers it
for this family, for every number of copies; it does not settle the general question.
The certificate is one case of a general identity. If C^Γ = C, then for every one-copy operator C̃ the operator
Q = ½[C^{⊗n} + C̃^{⊗n} − (Γ C̃)^{⊗n}] satisfies Q + Q^Γ = C^{⊗n}, because the partial transpose acts copy by copy.
Theorem 7.18 takes C̃ = |C|, for which ½(C^{⊗n} + |C|^{⊗n}) ⪰ 0 always, and the content is that the subtracted term
is dominated exactly when uⁿ ≤ 2. A positive semidefinite Q of this kind is a decomposition of C^{⊗n}, so for the
qutrit map the identity cannot help from three copies on (Proposition 7.3).
The theorem also proves the decomposability
thresholds of 𝒬_a computed after Proposition 7.3 for every n. The certificate was found by computing the
minimum-norm solution of that linear program in the form P + P^Γ for n = 2, 3, 4 and reading off its pattern.
`code/plane_threshold_check.py` (`logs/plane_threshold_check.log`) checks the matrix γ from explicit partial transposes,
the three identities exactly, the identity P + P^Γ = C^{⊗n} exactly in Q(2^{1/n}) for n ≤ 4 and with explicit
4ⁿ × 4ⁿ matrices for n ≤ 5, and that it fails 1% above the threshold. On a plane the threshold question is thus one about
decomposable maps. For qutrits it is not, from three copies on (Proposition 7.3), and Conjecture QT remains open there.

The same seems to hold for every Pauli map, which would settle the question that [FM] leave open (they give necessary
and, separately, sufficient conditions for every n). For random Pauli maps (eigenvalues 1, λ₁, λ₂, λ₃, neither
completely positive nor completely copositive), the exact decomposability of Υ^{⊗n} (a linear program over Bell
strings) agrees with their GHZ condition ([FM, Proposition 7, eq. (38) of arXiv v3]) in every sample for n = 2, 3, 4 (184, 184 and 140 maps). Along
65 random rays from the identity map, n = 3 and 4, the decomposability boundary and the GHZ boundary agree to 0.1%.
So, numerically, Υ^{⊗n} is positive exactly when it is decomposable, and exactly when the GHZ condition holds; this is
evidence for a positive answer to Müller-Hermes's question above. For n = 2 the first equivalence is his theorem. Beyond the family of Theorem 7.18 we have no proof. The
certificate of that theorem is not enough in general: for eigenvalues (1, t, 0, s) with t ≠ s, a single tensor-power
correction works at n = 2 but not at n = 3 (floating point; `code/pauli_maps/`).

**Theorem 7.19 (different parameters on planes, every n).** For a₁, …, a_n > 0 the following are equivalent:
(i) 𝒬_{a₁} ⊗ ⋯ ⊗ 𝒬_{a_n} is positive; (ii) it is decomposable; (iii) for every nonempty S ⊆ {1, …, n},

    ∏_{s∈S} (2a_s − ½) ≥ ½.

Only the copies with a_s < 3/4 matter in (iii): it says that a_s ≥ ½ for every s and ∏_{s: a_s<3/4} (2a_s − ½) ≥ ½. In
the language of [FM], with τ_s = 2|t_s|: a product Υ_{t₁} ⊗ ⋯ ⊗ Υ_{t_n} of maps from the family (1, t, 0, t) is positive
iff ∏_{s∈S} τ_s ≤ 2 for every S, and then it is decomposable. Theorem 7.18 is the case of equal parameters. Products
with maps outside the family can behave differently: Müller-Hermes [J. Math. Phys. 62, 092204 (2021), Theorem 8.4]
shows that the product of the map (1, t, t, t) with a map (1, b, 0, b) of this family can be positive without being
decomposable.

*Proof.* If (iii) fails for S, take the vectors of Theorem 7.1 on the copies in S (for |S| = 1, x_s = y_s = e₁) and real
unit vectors x_s ⊥ y_s at the other copies. The value is 2^{−|S|}[∏_{s∈S}(2a_s − ½) − ½] ∏_{s∉S} a_s < 0. For
(iii) ⇒ (ii) put τ_s = 1/(2a_s − ½) ∈ (0, 2], so that 𝒬_{a_s} = (2a_s − ½) Ad_Y ∘ Υ_{τ_s/2}, and let L = {s : τ_s > 1};
(iii) says ∏_{s∈L} τ_s ≤ 2. For s ∉ L, Υ_{τ_s/2} is completely positive and equal to its composition with the transpose,
so it is entanglement breaking, and tensoring with it keeps decomposability. The one-pair identities in the proof of
Theorem 7.18 hold at each copy with its own u = τ_s. So if ∏_{s∈L} τ_s = 2, the operator P built over the copies in L from
F = ⊗(1 + τ_s, 1, 1, 1 − τ_s), F̃ = ⊗(1 + τ_s, 1, 1, τ_s − 1) and R = ⊗(2/τ_s, 1, 1, 0) satisfies
P + P^Γ = F + (½∏τ_s − 1)R + (½ − 1/∏τ_s)F̃ = F, and p ≥ 0 as before because every τ_s ≥ 1. If ∏_{s∈L} τ_s < 2, raise
τ_{s₀} for one s₀ ∈ L to 2/∏_{s∈L∖s₀} τ_s ≤ 2. Then Υ_{τ_{s₀}/2} is a convex combination of Υ at the raised value and the
depolarizing map Υ₀. The second term is an entanglement-breaking map tensored with a product over L ∖ s₀, which is
decomposable by induction on |L|. ∎

`code/plane_products_check.py` (`logs/plane_products_check.log`, PASS) checks the identity and p ≥ 0 exactly at eight
rational parameter vectors with product 2 (n ≤ 4), with explicit matrices for n ≤ 3, a control with product 21/10, the
witness values, and the sign of the least value at 40 random parameter vectors (n = 3).

The same condition seems to govern the qutrit map.

**Conjecture QT_m.** For a₁, …, a_n > 0, 𝒫_{a₁} ⊗ ⋯ ⊗ 𝒫_{a_n} is positive iff condition (iii) of Theorem 7.19 holds.

Conjecture QT is the case of equal parameters. The witnesses of Theorem 7.19 live on planes, so (iii) is necessary. It
is also sufficient for n = 2. The reduction of Theorem 7.2 to x = Σ s_i e_i ⊗ e_i works with different parameters. The
first block satisfies zᵀM_Az − (a₁ − ½)(a₂ − ½)|z|² = Σ_{k<l}[c(s_l²z_k² + s_k²z_l²) + ¼ s_ks_lz_kz_l] with
c = ½(a₁ + a₂) − ¼ ≥ ⅛. On span{|ij⟩, |ji⟩} it is enough to take s_i² + s_j² = 1 (scaling s_i, s_j up lowers the diagonal
and raises the off-diagonal entry). With s_i² = ½ + δ the block is [[K − δΔ, ⅛s_is_j], [⅛s_is_j, K + δΔ]], where
K = a₁a₂ − ¼(a₁ + a₂) and Δ = ½(a₂ − a₁), and its determinant is K² − 1/256 − δ²(Δ² − 1/64). At δ² = ¼ this is
a₁a₂(a₁ − ½)(a₂ − ½) ≥ 0, and at δ = 0 it is ≥ 0 exactly when K ≥ 1/16, which is the case S = {1, 2} of (iii) (the
cases |S| = 1 say a_s ≥ ½). If Δ² > 1/64 the determinant
decreases in δ², so it is ≥ 0 throughout; otherwise its minimum is at δ = 0. For n = 3 and n = 4, alternating
minimisation at 21 and 17 parameter vectors, some within 10⁻³ of the boundary on either side, always finds a least value
of the predicted sign (`code/multiparam_numerics.py`, `logs/multiparam_real_n3.log` and `logs/multiparam_real_n4.log`;
3000 random pairs for n = 2 in `logs/multiparam_two.log`).

The searches suggest more: the least value itself seems to be attained on planes.

**Conjecture PM (the least value is planar).** For a₁, …, a_n ≥ ½, the least value of ⟨y|𝒫_{a₁} ⊗ ⋯ ⊗ 𝒫_{a_n}(xxᵀ)|y⟩
over unit x, y is

    m(a) = min{ ∏_s (a_s − ½),  min_{|S|≥2} 2^{−|S|}[∏_{s∈S}(2a_s − ½) − ½] ∏_{s∉S} a_s },

the value of aligned real products, or of the vectors of Theorem 7.1 on the copies in S with orthogonal real products at
the others. In particular it is attained by vectors that lie on one real plane at every copy.

Conjecture PM implies Conjecture QT_m, because m(a) ≥ 0 exactly when (iii) holds. At a_s = ¾ it says that the least
value is 4^{−n}, which Theorem 7.13 proves for n ≤ 5. For n = 2 it follows from the computation above: the first block
has least eigenvalue (a₁ − ½)(a₂ − ½), and the least eigenvalue K − √(1/256 + δ²(Δ² − 1/64)) of a 2 × 2 block is
monotone in δ², so it is attained at δ = 0, where it is ¼[(2a₁ − ½)(2a₂ − ½) − ½], or at δ² = ¼, where it is
min{a₂(a₁ − ½), a₁(a₂ − ½)} ≥ (a₁ − ½)(a₂ − ½). Alternating minimisation from random and thick starts (no planar
seeds) never goes below m(a): at 11 parameter vectors for n = 3 and 8 for n = 4, from a_s = ½ (deep inside the
non-positive region) to a_s = ¾, the least value found equals m(a) to 3·10⁻⁸. The minimisers found are thin at every
copy, except at the vertex (½, ¾, ¾), where the least value 0 is attained on a large set (`code/planar_minimum.py`,
`logs/planar_minimum_n3.log`, `logs/planar_minimum_n4.log`). At n = 5 the same search does not go below m(a) either, but
it does not always reach it.

**Proposition 7.20 (the vertices are 2-positivity).** Let Φ = 𝒫_{a₂} ⊗ ⋯ ⊗ 𝒫_{a_n} with real a_s. Then 𝒫_{1/2} ⊗ Φ is
positive iff id₂ ⊗ Φ is positive on real inputs. In particular the point (½, ¾, …, ¾) satisfies Conjecture QT_m iff
𝒫_{3/4}^{⊗(n−1)} is 2-positive on real inputs, and these statements hold for every n iff 𝒫_{3/4} is tensor-stable
positive.

*Proof.* For 2 × 2 matrices Tr(Y)I₂ − Y = σ_yYᵀσ_y and Tr(Y)I₂ − Yᵀ = σ_yYσ_y, so 𝒬_{1/2}(Y) = ¼σ_y(Y + Yᵀ)σ_y. The
output of Φ is invariant under the transpose, and Φ maps real symmetric matrices to real symmetric ones. So for real
X ⪰ 0, (𝒬_{1/2} ⊗ Φ)(X) = ½(σ_y ⊗ I)Z(σ_y ⊗ I) with Z = (id₂ ⊗ Φ)(X), and positivity of 𝒫_{1/2} ⊗ Φ, compressed to a real
plane at the first copy, gives Z ⪰ 0. Conversely 𝒫_{1/2} = R + T ∘ R, where R(X) = ¼(Tr(X)I − Xᵀ) is completely positive
with real Kraus operators of rank two (Proposition M, item 4). If id₂ ⊗ Φ is positive on real inputs, then so is R ⊗ Φ,
and (T ∘ R) ⊗ Φ = T ∘ (R ⊗ Φ). Real inputs suffice for 𝒫_{1/2} ⊗ Φ, as for 𝒫_a^{⊗n}. For the last sentence:
tensor-stable positivity gives 2-positivity of every power (Proposition M, item 3), and real 2-positivity of
𝒫_{3/4}^{⊗(n−1)} gives positivity of 𝒫_{3/4}^{⊗n} as in item 4. ∎

So Conjecture QT_m says that copies with a_s ≥ 3/4 cost nothing and that the others share a budget of log 2:
Σ_s max(0, −log(2a_s − ½)) ≤ log 2. In the coordinates t_s = log(2a_s − ½) the conjectured region is a polyhedron, its
vertices are the n points of Proposition 7.20 (one t_s = −log 2, the others 0), and the point of Conjecture QT,
t_s = −(log 2)/n, is their barycentre. So if the set of parameters at which 𝒫_{a₁} ⊗ ⋯ ⊗ 𝒫_{a_n} is positive is convex
in these coordinates, real 2-positivity of 𝒫_{3/4}^{⊗(n−1)} gives Conjecture QT at n copies with its sharp threshold a_n.
We do not know whether it is convex.

Numerically the vertices hold for n ≤ 6: the searches of Section 6 (Observation 1) find no Schmidt-rank-two vector with a
negative value of C_{3/4}^{⊗m} for m ≤ 5. They are harder for the method of Theorem 7.13, though. At the three-copy
vertex (½, ¾, ¾) the second level of the hierarchy reaches only t ≈ −4.7·10⁻³ (floating point,
`code/dps/level2_multi.py`, `code/dps/logs/sdp/vertex_n3_scs.log`; at the two-copy vertex it reaches t = 0), while at
a = ¾ it certifies five copies. So 2-positivity of 𝒫_{3/4} ⊗ 𝒫_{3/4} is not proved here.

**Remark 7.21 (there is no complex version).** The value f(a₁, …, a_n) = ⟨y|𝒫_{a₁} ⊗ ⋯ ⊗ 𝒫_{a_n}(xxᵀ)|y⟩ is multi-affine in
the parameters. At n = 3 (and in shorter runs at n = 4), searches over x, y and complex parameters find no zero with Re a₁ > ½ and Re a_s > 3/4 for
s ≥ 2: the largest Re a₁ of such zeros is ½, attained by products and approached by the plane witness. A zero-free region
of this kind would give the real statements by continuity, but it fails from n = 5. For
x = √p e₂^{⊗n} + √(1 − p) e₁^{⊗n} and y = e₁^{⊗n} the cross terms vanish and f = p∏a_s + (1 − p)∏(a_s − ½). At n = 5
this is zero for a₁ = 0.5300 + 0.0995i, a₂ = ⋯ = a₅ = 0.7600 + 0.4445i and p ≈ 0.0220, whose real parts satisfy (iii).
The reason is that arg a_s − arg(a_s − ½) can be close to −90° when Re a_s is close to ½ and close to −30° when
Re a_s > 3/4, and the copies reach −180° together. With every Re a_s > 3/4 the same family needs seven copies:
a_s = 0.76 + 0.2510i for every s and p ≈ 0.0038 give a zero at n = 7
(`multiparam_numerics.py complex`). So there is no half-plane (Lee–Yang) version of Conjecture QT_m; it is a statement
about real parameters.

Given Conjecture QT_m for fewer copies, the value is nondecreasing in every a_s on the conjectured region, so at n copies
it is enough to look at the facet ∏_s(2a_s − ½) = ½ with every a_s ∈ [½, ¾] (in the coordinates t_s this is a simplex,
spanned by the vertices of Proposition 7.20). The vectors of Theorem 7.1 give the value 0 at every point of this facet.
The next proposition shows that they cannot be pushed below 0 by leaving the plane, to second order.

**Proposition 7.22 (the plane witness is locally stable, every n).** Let a_s ∈ [½, ¾] with ∏_s(2a_s − ½) = ½, let x₀, y₀
be the vectors of Theorem 7.1 on the plane spanned by e₁, e₂, and f(x, y) = ⟨y|𝒫_{a₁} ⊗ ⋯ ⊗ 𝒫_{a_n}(xxᵀ)|y⟩. Then
f(x₀, y₀) = 0, the first derivatives of f vanish there, and the Hessian of f at (x₀, y₀) is positive semidefinite on all
of (R³)^{⊗n} × (R³)^{⊗n}. If every a_s < ¾, it is positive on perturbations that leave the plane at some copy, except
along the 2n directions that rotate the plane at a single copy (copies with a_s = ¾ add further zero directions).

*Proof.* With M_x = ⊗_s𝒫_{a_s}(x₀x₀ᵀ) and M_y = ⊗_s𝒫_{a_s}(y₀y₀ᵀ), the second-order part of
f(x₀ + δx, y₀ + δy) is Q = ⟨δx|M_y|δx⟩ + ⟨δy|M_x|δy⟩ + 2⟨δy|⊗_s𝒫_{a_s}(x₀δxᵀ + δx x₀ᵀ)|y₀⟩. Split every vector by
the set T of copies at which it points along e₃ (in the plane at the other copies). 𝒫_a maps a matrix supported on the
plane to a matrix that is block diagonal for plane ⊕ e₃, with entry a·Tr at e₃, and it maps e_ie₃ᵀ (i ≤ 2) to
−¼(e_ie₃ᵀ + e₃e_iᵀ). So M_x, M_y and the cross term are block diagonal in T. For T = ∅ the block is the Hessian of the
plane problem at a zero, positive semidefinite by Theorem 7.19; the same theorem gives M_x, M_y ⪰ 0, hence M_x y₀ = 0
and M_y x₀ = 0. For T ≠ ∅, the reduced states of x₀ and y₀ on the other copies are ½(P^{⊗} + P̄^{⊗}) with P = pp*, and
𝒬_a(pp*) = (a − ¼)I₂, so the diagonal part of the block is A_T(|δx_T|² + |δy_T|²) with
A_T = ∏_{s∈T}a_s∏_{s∉T}(a_s − ¼). Expand δx_T = Σ_σ α_σ σ over σ ∈ {p, p̄}^{T^c} (α_σ̄ = ᾱ_σ), and δy_T likewise with
β_σ. The cross term becomes 2(−¼)^{|T|} Σ_σ w_σ β_σ ᾱ_σ, where w_σ = ∏_s c_s(σ_s) + ∏_s c̄_s(σ_s) for odd n and
w_σ = −i(∏_s c_s(σ_s) − ∏_s c̄_s(σ_s)) for even n, with c_s(p) = c̄_s(p̄) = a_s − ¼ and c_s(p̄) = c̄_s(p) = −½. So the
block splits into 2 × 2 forms with least eigenvalue
A_T − 4^{−|T|}|w_σ|. Put r_s = 1/(2a_s − ½) ≥ 1 and R = ∏_{s∉T} r_s, which equals 2∏_{s∈T}(2a_s − ½) on the facet.
Then |w_σ| ≤ (x + R/x)∏_{s∉T}(a_s − ¼) for some x ∈ [1, R], and x + R/x ≤ 1 + R. So it suffices that
∏_{s∈T}4a_s ≥ 1 + 2∏_{s∈T}(2a_s − ½). For |T| = 1 this is an identity, and a 2 × 2 form is singular only when
x ∈ {1, R}. If every a_s < ¾ (all r_s > 1), that means σ = p^{⊗} or σ = p̄^{⊗}, the two rotations at that copy; copies with
a_s = ¾ (r_s = 1) let further σ qualify. For |T| = t ≥ 2 the left side is at least 2^t ≥ 4 and the right side at most 3. ∎

So, to second order, no direction out of the plane lowers the value near the witnesses, at any point of the facet and
for every n; in particular at a_s = a_n, the threshold of Conjecture QT. This matches Observation 7.4, where every
negative minimiser found below a_n lives on one plane per copy. `code/plane_witness_hessian.py`
(`logs/plane_witness_hessian.log`, PASS) computes the full Hessian for n = 2, …, 6 at five facet points each and compares
the least positive eigenvalue off the plane with the closed form A_T − 4^{−|T|}|w_σ| (agreement to 10⁻⁹). At the points
with every a_s < ¾ (the centre and the random points) the kernel has dimension 2n + 3 (rotations and the two scalings);
at the edge points and vertices, where some a_s = ¾, it is larger.

Corollary 7.10 asks whether a thick input can be a minimiser. The next proposition says what a single copy can do about
it. At a critical point, squeezing one copy toward a plane moves the value monotonically, all the way to the
compression.

**Proposition 7.23 (squeeze paths through critical points, every n).** Let a be real and let x, y be real unit vectors
on n copies with 𝒫_a^{⊗n}(yyᵀ)x = λx and 𝒫_a^{⊗n}(xxᵀ)y = λy, where λ = ⟨y|𝒫_a^{⊗n}(xxᵀ)|y⟩. (Every local minimum
of f(x, y) = ⟨y|𝒫_a^{⊗n}(xxᵀ)|y⟩ on the product of the unit spheres is such a joint critical point.) Fix a copy s and a
unit w ∈ R³, let ξ = ⟨w|_s x and η = ⟨w|_s y be the contractions at copy s, α = |ξ|², β = |η|², and for 0 ≤ v ≤ 1 let
G_v = I − (1 − √(1 − v))wwᵀ act at copy s (G_1 projects copy s onto w^⊥). Then

    f(G_v x, G_v y) / (|G_v x|² |G_v y|²) = λ + v² [ (a − ½) f_{n−1}(ξ, η) − λαβ ] / [ (1 − vα)(1 − vβ) ],

with f_{n−1}(ξ, η) = ⟨η|𝒫_a^{⊗(n−1)}(ξξᵀ)|η⟩. So the normalised value is monotone along the whole squeeze, and the sign
of the single number (a − ½) f_{n−1}(ξ, η) − λαβ decides whether it goes up or down.

(ii) Let Q_s be the second derivative at t = 0 of f(e^{tH}x, e^{tH}y)/(|e^{tH}x|²|e^{tH}y|²), a quadratic form in the
real symmetric 3 × 3 matrix H acting at copy s. Its trace over an orthonormal basis of the traceless H is

    −4 [ λ (4a − 1 + 2 Tr ρ^x_s ρ^y_s) − (a − ½)(1 + 4a) τ_s ],   τ_s = ⟨ρ^x_{−s}, 𝒫_a^{⊗(n−1)}(ρ^y_{−s})⟩,

where ρ^x_s, ρ^y_s are the one-copy reduced states at copy s and ρ^x_{−s}, ρ^y_{−s} those of the other n − 1 copies. In
particular, if λ (4a − 1 + 2 Tr ρ^x_s ρ^y_s) > (a − ½)(1 + 4a) τ_s at some copy, then (x, y) is not a local minimum.

*Proof.* Write x = x′ + w ⊗ ξ and y = y′ + w ⊗ η, with x′ and y′ supported on w^⊥ at copy s, and put
x_u = x′ + u w ⊗ ξ, y_u = y′ + u w ⊗ η, so that G_v x = x_u and G_v y = y_u for u = √(1 − v). The value f(x_u, y_u)
is a polynomial of degree four in u. It is even, because the reflection R = I − 2wwᵀ at copy s maps (x_u, y_u) to
(x_{−u}, y_{−u}) and 𝒫_a(RXRᵀ) = R𝒫_a(X)Rᵀ. So f(x_u, y_u) = A + Bu² + Cu⁴. The coefficient of u⁴ comes from w ⊗ ξ and
w ⊗ η alone, so C = ⟨ww|C_a|ww⟩ f_{n−1}(ξ, η) = (a − ½) f_{n−1}(ξ, η). At u = 1, A + B + C = λ. The derivative at u = 1
is 2⟨w ⊗ ξ|𝒫_a^{⊗n}(yyᵀ)|x⟩ + 2⟨w ⊗ η|𝒫_a^{⊗n}(xxᵀ)|y⟩ = 2λ(α + β), by the critical equations and ⟨w ⊗ ξ|x⟩ = α.
So 2B + 4C = 2λ(α + β), and with t = 1 − u² = v,

    f(x_u, y_u) = λ [1 − (α + β)t] + C t²,   |x_u|² |y_u|² = (1 − αt)(1 − βt) = 1 − (α + β)t + αβ t².

Subtracting λ times the second line from the first gives the formula.
(ii) Scaling H ↦ H + cI does not change the normalised value, so Q_s(H) = Q_s(H − (Tr H/3)I) and the trace over the
traceless H is (15/2) E_w[Q_s(wwᵀ)] for a uniformly random unit w (E[(wwᵀ − I/3) ⊗ (wwᵀ − I/3)] is 2/15 times the
projection onto traceless symmetric matrices). Expanding (i) to second order along e^{−t wwᵀ} (u = e^{−t}, v ≈ 2t) gives
Q_s(wwᵀ) = −8 [λαβ − (a − ½) f_{n−1}(ξ, η)]. Finally E[αβ] = (1 + 2 Tr ρ^x_s ρ^y_s)/15, and with Z_{(ik),(jl)} =
⟨x_i ⊗ y_k|C_a^{⊗(n−1)}|x_j ⊗ y_l⟩ (x_i = ⟨e_i|_s x and so on), E[f_{n−1}(ξ, η)] = (Tr Z + Tr FZ + ⟨Ω|Z|Ω⟩)/15, where
Tr Z = τ_s and λ = a Tr Z − ¼(Tr FZ + ⟨Ω|Z|Ω⟩). Putting these together gives the trace. ∎

**Corollary 7.24.** Let a ≥ ¾, and suppose μ_m(a) = (a − ½)^m for every m < n (μ_m as in Proposition 7.8; at a = ¾
Theorem 7.13 gives this for m ≤ 5). If μ_n(a) < (a − ½)^n, then every pair of unit vectors at which f attains μ_n is a
joint critical point, both vectors are thick at every copy, and for every copy s and every unit w

    (a − ½) f_{n−1}(⟨w|_s x, ⟨w|_s y) > μ_n |⟨w|_s x|² |⟨w|_s y|²,

so that every one-copy squeeze strictly raises the normalised value along its whole path. (Given μ_{n−1}, the display
holds automatically at any critical point with value below (a − ½)^n, so this part only says where a minimiser would
have to sit.) Hence 𝒫_a is tensor-stable positive (and μ_n(a) = (a − ½)^n for every n) if, for every n, no joint
critical point of f_n at which both vectors are thick at every copy, at any value, is a strict minimum along every
one-copy squeeze; that is, if each such point has a copy s and a direction w with
(a − ½) f_{n−1}(⟨w|_s x, ⟨w|_s y) ≤ f(x, y) |⟨w|_s x|² |⟨w|_s y|².

*Proof.* A minimiser is a joint critical point (Lagrange conditions on each sphere). It is thick at every copy on both
sides: if x has a thin copy, Proposition 7.8 gives λ_min(𝒫_a^{⊗n}(xxᵀ)) ≥ (a − ½)μ_{n−1} = (a − ½)^n, and if y has one,
the same applies to ⟨y|𝒫_a^{⊗n}(xxᵀ)|y⟩ = ⟨x|𝒫_a^{⊗n}(yyᵀ)|x⟩. Thickness also gives α, β < 1. The compression G_1x,
G_1y is thin at copy s, so the left side of the formula at v = 1 is at least (a − ½)^n > μ_n = λ. The last sentence
follows by induction on n, starting from μ₁ = a − ½. ∎

This is the map's version of the squeeze-orbit picture of `orbit-notes/`, where rank-one squeezes through critical points
of the endpoint form are monotone in the same way (item 2 there). For the endpoint form, the statements that look at
one copy at a time are false at n = 3 (thick operators that are minima over every single-copy orbit; see "One copy at a
time is not enough" there), while the joint statement survives. For the map, such points exist even for joint
squeezes: Observation 7.26 gives a balanced critical pair at four copies that rises along every one-copy squeeze and
along every joint squeeze of all four copies. So the last sentence of the corollary cannot be used. What a proof needs
is only the version for thick critical points below (a − ½)^n, and that version is equivalent to the problem itself.

**Observation 7.25 (numerical).** At n = 3 and a = ¾, minimising the residual of the critical equations from random
starts gives thick joint critical points: five of them, with values 0.178 to 0.210 (far above 1/64) and least one-copy
marginal eigenvalues at least 0.08. A second search first pulls the value toward a target (0.05, 0.1 or 0.15) with a
floor on the thickness and then solves the critical equations alone. It finds thick critical points only with values
0.122 to 0.174, none lower (one run aimed at 0.05 ended at a thin point; `code/thick_critical_points.py`,
`logs/thick_critical_points.log`). At every thick critical point found, some copy and direction have
λ − (a − ½) f₂(ξ̂, η̂) > 0 (λ the critical value, ξ̂, η̂ the normalised contractions), with margins 0.076 to 0.12 in the second search. So each
is a saddle along a one-copy squeeze, which then lowers the value all the way to a thin pair. A plain census from
random starts (`code/critical_census.py`) finds the same picture at n = 3 (ten thick critical points, values 0.151 to
0.239, margins at least 0.11; `logs/critical_census_n3.log`) and at n = 4 (six thick critical points, values 0.113 to
0.118, against the product value 1/256, margins at least 0.07; `logs/critical_census_n4.log`). The trace of
Proposition 7.23(ii) gives the same conclusion without choosing a direction. Write r_s for the ratio
(a − ½)(1 + 4a)τ_s / [λ(4a − 1 + 2 Tr ρ^x_s ρ^y_s)]; the trace is negative exactly when r_s < 1. At every copy of every
thick critical point in `logs/squeeze_trace.log`, r_s ≤ 0.70 (eleven points at n = 3) and r_s ≤ 0.66 (two points at
n = 4). At aligned products r_s = 1, and at the plane GHZ pairs of Theorem 7.1, which are thin, r_s = 4/3 at a = ¾ for
every n (there the trace is positive, and these pairs descend in other directions). None of this shows that thick
local minima cannot exist. The inequality r_s < 1 is checked only at the critical points that the searches happened
to find, at n ≤ 4, while a failure would need n ≥ 7. It is not a general property of critical points, since it fails
at the plane GHZ pairs. And at the minimiser of a failure it would fail automatically, since there λ < (a − ½)^n while
τ_s ≥ (a − ½)^{n−1}. Proving r_s < 1 at every thick critical point, for every n, is as hard as the problem itself.
Descents from twelve random starts all end at thin points (eleven at 1/64). A stronger statement also held in
every test. Let g(y) = λ_min(𝒫_{3/4}^{⊗n}(yyᵀ))/|y|². Then some one-copy compression G_1 has g(G_1y) ≤ g(y). We
checked this for all Schmidt vectors at n = 2, where y can be taken as Σ s_i e_i ⊗ e_i (the largest gap
min g(G_1y) − g(y) on a 25 × 25 grid of s is −9·10⁻⁵, approached only near products). We also checked random y at n = 3
and 4, and adversarial hill climbing, which gives largest gaps −1.3·10⁻² at n = 3 and −1.6·10⁻² at n = 4. At n = 2 the
best direction was the middle Schmidt direction in every test. If this compression statement held at every n, then with
Proposition 7.8 (g(y) ≥ (a − ½)μ_{n−1} when y has a thin copy) it would give μ_n = 4^{−n} for every n by induction.
Keeping the optimal input and compressing it together with y also worked, but adversarial searches drive that margin to
zero as y becomes thin. The direction cannot be found by averaging over an orthonormal basis {w_k} (this fails already
at n = 2), although for every x, y and basis Σ_k f(G^{(k)}x, G^{(k)}y) = f(x, y) + (a − ½) Σ_k f_{n−1}(⟨w_k|_s x, ⟨w_k|_s y),
with G^{(k)} the compression of copy s onto w_k^⊥ (this is Remark 7.5).
The compression statement fails at four copies (added 8 November 2026). For the partner vector y = p₂ − c(p₁ + p₃) of Observation 7.26,
g(y) = (151 − √6805)/6·4^{−4} ≈ 11.418·4^{−4}, and every one-copy compression, at any copy and in any direction (the
input being re-optimised), gives g ≈ 11.733·4^{−4}. All compressions give the same value, because y is invariant under
rotations of all copies at once and the copy permutations that fix it act transitively. A nearby vector whose one-copy
marginals are not degenerate also has a positive gap at every copy (0.23 to 0.33 times 4^{−4}), so the failure is not
an artefact of the symmetry. The failure reaches the range n ≥ 7 that an induction would need: for products,
g(y ⊗ y′) = g(y) g(y′) (both outputs are positive semidefinite), and a one-copy compression acts inside one factor, so
on y ⊗ y (eight copies, thick at every copy) every one-copy compression raises g as well. At three copies,
compressing along a fixed eigenvector of the one-copy marginal already
fails: along the middle eigenvector at an input with smallest marginal eigenvalues 0.06, 0.25 and 0.20 (gap +0.28·4^{−3}
at every copy), and along the smallest eigenvector at a nearly thin input (+0.0017·4^{−3}). In adversarial searches the
best of the three eigenvectors never failed at three copies. So the compression route, like the squeezes of
Observation 7.26, ends at four copies.
Compressing two copies at once survives these tests. At the partner vector the best compression of a pair of copies,
each onto a plane, lowers g by 3.80·4^{−4}, and by 4.30·4^{−4} for the two copies of an Ω-pair of p₂. A gradient ascent
on y that tries to make every pair compression raise g (four runs of 60 steps, from near the partner and from random
inputs) kept the best pair compression at least 1.9·4^{−4} below g(y); from random inputs it only moves y toward
thinner vectors. At six copies, at the pair that minimises the value over spin singlets (48.28·4^{−6}, which is also
the least value of g at that test vector over all inputs), a one-copy compression already lowers g, by 5.76·4^{−6},
and an adversarial search over spin-singlet test vectors at six copies (four runs) found none at which every one-copy
compression raises g (the largest gap found is −6.4·4^{−6}). So far the four-copy partner and its products are the only
inputs we know that defeat one-copy compressions.
If for every y some pair of copies could be compressed onto planes without raising g, induction on n would give
μ_n = 4^{−n} for every n (Proposition 7.8 with two thin copies). But at the minimiser of a hypothetical failure every
such compression would have to raise g, so this statement is no easier than the problem itself. We record it as a
pattern, not as a route.
`code/compression_rules/compression_rules_check.py` (`logs/` there, PASS) checks the four-copy counterexample (all copies,
200 directions each, and one value in 30-digit arithmetic) and the two three-copy inputs (floating point);
`compression_pairs_check.py`, `adversary_pairs_n4.py` and `singlet_search_n6.py` there give the pair compressions and the
six-copy singlet search (floating point).
`code/squeeze_paths.py` (`logs/squeeze_paths.log`, PASS) checks the formula of Proposition 7.23 at minima and at thick
saddles for a = ¾ and 0.7 (agreement to 10⁻⁸, limited by the residual of the critical points), and the compression
gaps at n = 2, 3, 4. `code/squeeze_trace.py` (`logs/squeeze_trace.log`, PASS) checks the trace of 7.23(ii) against
the closed form of the Hessian and against finite differences, at products, plane GHZ pairs and thick saddles
(a = ¾ and 0.7), and computes the ratios r_s.

**Observation 7.26 (a balanced critical pair that defeats squeezes).** At n = 4 let p₁ = Ω₀₁Ω₂₃, p₂ = Ω₀₂Ω₁₃ and
p₃ = Ω₀₃Ω₁₂, where Ω_{st} = Σ_i e_i ⊗ e_i on copies s, t (unnormalised). Let x ∝ p₃ − p₁ and y ∝ p₂ − c(p₁ + p₃), with
c = (√6805 − 71)/36 ≈ 0.3192, both normalised. Then, exactly,

    𝒫_{3/4}^{⊗4}(yyᵀ) x = λx,   𝒫_{3/4}^{⊗4}(xxᵀ) y = λy,   λ = (151 − √6805)/1536 ≈ 0.04460,

and every one-copy reduced state of x and of y is I/3. The rest is numerical. λ is the least eigenvalue of
𝒫_{3/4}^{⊗4}(yyᵀ). Every one-copy squeeze raises the normalised value: the margin λ − ¼ f₃(ξ̂, η̂) is at most −0.0168 at
every copy. So does every joint squeeze: the squeeze Hessian over all traceless directions at all four copies is
positive definite, with least eigenvalue 0.0118. The full Hessian on the product of spheres has five negative eigenvalues
(−0.0028), so the pair is a saddle of f, but only through directions that are not local filters. Both vectors are spin
singlets (invariant under O(3) acting on all copies at once), and λ is the least value of f over pairs of singlets at
n = 4. It is 11.4 times the product value 1/256. Consequently, for the map, no squeeze-orbit argument (one copy, or all
copies jointly) can rule out thick minimisers. For the endpoint form the joint version survives at the one-copy
counterexamples of `orbit-notes/`; here it does not. `code/thick_sector.py` (`logs/thick_sector.log`, PASS) checks the
identities in exact arithmetic (modulo 36c² + 142c − 49) and the rest in floating point.

**Observation 7.27 (thick inputs sit far above the product value).** All of this is floating point.
(a) *A thickness floor.* Minimise f_n over pairs whose one-copy marginals all have least eigenvalue at least τ (a penalty
method from a few starts, so the values are upper bounds for the constrained minima). The per-copy factor
(least value)^{1/n} hardly depends on n:

| τ | n = 2 | n = 3 | n = 4 | n = 5 |
|---|---|---|---|---|
| 0.02 | 0.315 | 0.299 | 0.298 | 0.300 |
| 0.05 | 0.364 | 0.349 | 0.344 | 0.347 |
| 0.1 | 0.441 | 0.398 | 0.393 | 0.397 |
| 0.2 | 0.510 | 0.454 | 0.434 | 0.463 |
| 0.3 | 0.545 | 0.497 | 0.459 | 0.488 |

The product value corresponds to the factor ¼ (`code/thickness_floor.py`, `logs/thickness_floor.log`).
(b) *Singlets.* The least value of f over pairs of spin singlets is 23/48, 7/32, 0.04460, 0.02962, 0.01179 and 0.006630
for n = 2, …, 7. That is 7.7, 14, 11.4, 30, 48 and 109 times 4^{−n}. The per-copy factor is at least 0.4596, the value
of Observation 7.26, which tensor powers of that pair attain at every multiple of four (`logs/thick_sector.log`). Every
balanced pair we found (all one-copy marginals I/3) has a per-copy factor of at least 0.4596.
(c) *A multiplicative bonus.* With τ_s the smaller of the two thicknesses at copy s, the two-sided bound
⟨y|𝒫_{3/4}^{⊗n}(xxᵀ)|y⟩ ≥ ∏_s(¼ + ½τ_s) held in adversarial searches from random, near-thin and structured starts. With
a larger coefficient it fails: ¼ + 0.7τ_s at the pair of Observation 7.26, and ¼ + τ_s at the least balanced pair at
n = 2 (value 11/36). Near thin pairs the gain is linear in τ: a bonus √τ_s/3 fails there.
(d) *Exact values on balanced inputs.* At three copies, with ε the antisymmetric tensor, 𝒫_{3/4}^{⊗3}(εεᵀ/6) has
spectrum 7/64 (three times), 12/64 (seven times) and 14/64 (seventeen times), from the exact characteristic polynomial;
the per-copy factor is 0.4782. At six copies, group the copies into three pairs and put p = e₀ ⊗ e₀ and
e = (e₁ ⊗ e₁ + e₂ ⊗ e₂)/√2 on a pair. The input y = (pee + epe + eep)/√3 is balanced, and its least output eigenvalue is
9/1024 = 36·4^{−6} (floating point, fifteen digits), a per-copy factor of 0.4543. That is the lowest known for a
balanced input, still about 9% above the factor 5/12 that Conjecture TC asks for there. Its bottom test vector is
planar at every copy, so the pair is not balanced on both sides and (b) is unchanged. On two pairs the same
construction gives (63 − √2029)/512, a root of x² − 126x + 1940 for 512W (exact). (`code/exact_extremal_values.py`,
`logs/exact_extremal_values.log`.)

**Proposition 7.28 (the bonus at two copies).** Let y be a real unit vector on two copies and τ the least eigenvalue of
its one-copy reduced state. Then

    λ_min(𝒫_{3/4}^{⊗2}(yyᵀ)) ≥ min{1/16 + ¾τ, 1/8 + 7τ/16} ≥ (¼ + ½τ)².

*Proof.* By local rotations y = Σ_i s_i e_i ⊗ e_i, and τ = min_i s_i². By the proof of Theorem 7.2, 𝒫_{3/4}^{⊗2}(yyᵀ) is
block diagonal. On the first block, zᵀ(M_A − 1/16)z = Σ_{k<l}[½(s_l²z_k² + s_k²z_l²) + ¼s_ks_lz_kz_l]. Since
¼|s_ks_lz_kz_l| ≤ ⅛(s_l²z_k² + s_k²z_l²), each term is at least (3/8)(s_l²z_k² + s_k²z_l²) ≥ (3/8)τ(z_k² + z_l²), and the
sum is at least ¾τ|z|². On the block for i < j, with k the third index, the least eigenvalue is
9/16 − (3/8)(s_i² + s_j²) − ⅛s_is_j ≥ 9/16 − (7/16)(1 − s_k²) ≥ 1/8 + 7τ/16, using s_is_j ≤ ½(s_i² + s_j²). Finally
1/16 + ¾τ − (¼ + ½τ)² = τ(2 − τ)/4 ≥ 0, and 1/8 + 7τ/16 − (¼ + ½τ)² = 1/16 + 3τ/16 − τ²/4 > 0 for τ ≤ ⅓. ∎

So at two copies a thick input raises the least output eigenvalue above the product value 1/16, whatever the output
vector. Numerically the coefficient ½ can be raised to 3(√(13/48) − ¼) ≈ 0.811 at two copies, with the maximally
entangled y as the worst case (`code/thickness_bonus.py`, `logs/thickness_bonus.log`).

**Conjecture TC (a bonus for thickness).** For every n and every real unit y on n copies,

    λ_min(𝒫_{3/4}^{⊗n}(yyᵀ)) ≥ ∏_s (¼ + ½ τ_s(y)),

where τ_s(y) is the least eigenvalue of the one-copy reduced state of y at copy s. Equivalently,
⟨y|𝒫_{3/4}^{⊗n}(xxᵀ)|y⟩ ≥ ∏_s (¼ + ½ τ_s(y)) for all unit x; by symmetry the same holds with the thickness of x.

It contains Conjecture QT at a = ¾, with the least value 4^{−n} (take τ_s ≥ 0), and it is tight at products. It says
that a copy where the input has full support pays strictly more than the factor ¼ of a thin copy. Proposition 7.28 is
its case n = 2. Unlike thinness, the bonus does not pass through the peeling of Proposition 7.8, because the conditional
states there can be thinner than their average. So it is not inductive in an obvious way. Evidence: adversarial
searches at n = 3, 4 (`code/thickness_bonus.py`), and for the weaker two-sided form at n = 2, 3 (`code/thick_sector.py`).
The data of Observation 7.27 are consistent with it. With coefficient 0.7 instead of ½ it fails at the singlet pair of
Observation 7.26.

**Theorem 7.29 (the bonus at three copies; computer-assisted).** For every real unit vector y on three copies,

    λ_min(𝒫_{3/4}^{⊗3}(yyᵀ)) ≥ (1/64)(1 + (98/27)(τ₁ + τ₂ + τ₃)) ≥ ∏_s (¼ + ½τ_s),

with τ_s = τ_s(y) as in Conjecture TC. So Conjecture TC holds at n = 3.

The second inequality is elementary. ∏_s(1 + 2τ_s) − 1 − (98/27)Σ_sτ_s is affine in each τ_s, so its maximum on
[0, ⅓]³ is at a vertex. With k coordinates equal to ⅓ it is (5/3)^k − 1 − 98k/81, which is negative for k = 1, 2 and
zero for k = 3.

*Proof (computer-assisted).* Write t = 1/64 and κ = 98/27. For real x, y on three copies and real w₀, w₁, w₂ ∈ R³ let
W = |w₀|²|w₁|²|w₂|², let ρ^y_s = Tr_{≠s} yyᵀ, and let

    P(x, y, w) = f(x, y)|y|²W − t|x|²|y|⁴W − tκ|x|²|y|² Σ_s ⟨w_s|ρ^y_s|w_s⟩ W/|w_s|²,    f(x, y) = ⟨y|𝒫_{3/4}^{⊗3}(xxᵀ)|y⟩.

The certificate is an exact polynomial identity

    P = Avg_σ [ ⟨Ψ|Z₀ ⊗ I_w + Z₁|Ψ⟩ + Tr(Λ(x, y) G(y, w₀)) |w₁|²|w₂|² ] ∘ σ,

where σ runs over the permutations of the copies, Ψ = x ⊗ y ⊗ y ⊗ w with the strands x, y, y′, w at each copy,
G(y, w₀) = |w₀|²ρ^y_0 − ⟨w₀|ρ^y_0|w₀⟩I and Λ(x, y)_{ij} = ⟨e_i ⊗ x ⊗ y|L|e_j ⊗ x ⊗ y⟩. Here Z₀ acts on x, y, y′ at every
copy, Z₁ on x, y, y′, w at copy 0 and on x, y, y′ at copies 1 and 2, and L on (e_i, x, y) at copy 0 and on (x, y) at
copies 1 and 2. All three are positive semidefinite. If each w_s is a least eigenvector of ρ^y_s, then G ⪰ 0, so the
right side is nonnegative. At unit x, y, w_s this is f(x, y) ≥ t(1 + κΣ_sτ_s), the first inequality.

Z₀, Z₁ and L are rational combinations of Brauer diagrams, so they commute with O(3) at each copy. Positivity
therefore reduces to exact Gram blocks on the m = 0 vectors of each spin sector. Real product vectors x = y force
zeros: Sym³(R³) at a copy with strands x, y, y′; span{e ⊗ e ⊗ e ⊗ w : w ⊥ e} (spins 1 to 4) at the copy of Z₁ that
carries w; and e ⊗ e at the two-strand copies of L. The certificate has exactly these kernels and is positive definite
on the rest; the least LDLᵀ pivot of the reduced blocks is 0.0069. The identity is checked on orbit sums. On real
vectors, the quadratic form of a product of diagrams depends only on its orbit under exchanging the two x legs,
permuting the four y legs, exchanging the two legs of w_s at a copy, and (because of the average) permuting the copies.
At d = 3 the 105 diagrams on four strands span only 91 dimensions, and the identity uses 19 of the operator relations
among them. All 190 orbit equations hold exactly. There are two verifiers, written separately. The second rebuilds the
diagrams, the target, the group, the spin bases and the forced kernels from scratch, and takes only the operators from
the certificate. Both sides also agree at random real points, computed from the definition of 𝒫_{3/4}, to 3.4·10⁻¹⁵
relative (`code/tc3/`). ∎

The variable w is what makes this work. Simple stand-ins for τ_s in terms of x and y fail. For example
½(1 − Tr ρ^x_s ρ^y_s), which is at least τ_s(y) for unit vectors, gives an operator that is not block positive: take a
maximally entangled y on a plane and the singlet x in the same plane. With w but without the second output copy y′,
the same scheme reaches only κ ≈ 0.33 at two copies. With y′ it reaches the exact best constant 5 at two copies,
attained at the maximally entangled vector. At three copies it reaches κ ≈ 4.125, while the best constant is 5 to
search accuracy, attained at e ⊗ Ω̂ (`code/tc3/kappa_scan.py`, `code/tc3/tc3_numerics.py`).

Theorem 7.29 says nothing about more copies, and positivity at three copies was already known (Theorem 7.13). What is
new is the quantitative bonus for thick inputs. The bonus does not pass through conditioning on a copy. The vector
y = ½[e₁ ⊗ (e₁e₁ + e₂e₃) + e₂ ⊗ (e₁e₂ + e₃e₃)] is thin at the first copy and has τ = ¼ at the other two, yet every
conditional vector (fᵀ ⊗ I)y with f in the plane of the first copy has Schmidt coefficients (1/√2, 1/√2, 0). So peeling
as in Proposition 7.8, with Proposition 7.28 applied to the conditional vectors, gives only 1/32, below the TC value
9/256; the true value is 0.0774 (`code/tc3/tc3_numerics.py`). We have not tried four copies. There Z₁ would have about
3·10⁵ block coordinates.

**Proposition 7.30 (the bonus against product test vectors, every n).** Let ρ be any state on (C³)^{⊗n}, pure or mixed,
with one-copy reduced states ρ_s and τ_s = λ_min(ρ_s). For every real unit product vector v = v₁ ⊗ ⋯ ⊗ v_n,

    ⟨v|𝒫_{3/4}^{⊗n}(ρ)|v⟩ ≥ (¾)^n 3^{−Σ_s⟨v_s|ρ_s|v_s⟩} ≥ 4^{−n} 3^{2Σ_sτ_s} ≥ ∏_s (¼ + ½τ_s).

In particular Conjecture TC holds, for mixed inputs too, whenever the least eigenvalue of 𝒫_{3/4}^{⊗n}(ρ) is attained
at a real product vector. This covers every ρ that is diagonal in a real product basis, because the output is then
diagonal in the same basis.

*Proof.* 𝒫_{3/4} is self-adjoint for the Hilbert–Schmidt inner product, so the value is Tr ρ ⊗_s 𝒫_{3/4}(v_s v_sᵀ),
and 𝒫_{3/4}(v_s v_sᵀ) = ¾(I − ⅔ v_s v_sᵀ) = ¾·3^{−v_s v_sᵀ}. The projectors v_s v_sᵀ act on different copies and
commute, so the product is (¾)^n 3^{−N} with N = Σ_s v_s v_sᵀ. Jensen's inequality for the distribution of N in the
state ρ gives Tr ρ 3^{−N} ≥ 3^{−Tr ρN}, and Tr ρN = Σ_s ⟨v_s|ρ_s|v_s⟩ ≤ Σ_s λ_max(ρ_s) ≤ Σ_s (1 − 2τ_s). Finally
3^{2τ} ≥ 1 + 2τ. ∎

The exponential form does not extend to entangled test vectors. At the singlet pair of Observation 7.26, where all
one-copy marginals are I/3, the value (151 − √6805)/1536 ≈ 0.0446 is below 4^{−4}3^{8/3} ≈ 0.073. Padding that pair with
real products shows that a bound 4^{−n} exp(κ Σ_s τ_s), for all inputs and test vectors, needs
κ ≤ ¾ log((151 − √6805)/6) ≈ 1.83, while near thin inputs κ ≥ 2 would be needed to give Conjecture TC. So a failure of
Conjecture TC needs an entangled test vector and an input with coherences in every product basis, and the coherent
terms have to be controlled with their signs.

**Remark 7.31 (peeling one copy).** Write D(X) = Tr(X)I, so that 𝒫_{3/4} = ¼D + 𝒫_{1/2}. For every state ρ and every
copy s,

    𝒫_{3/4}^{⊗n}(ρ) = ¼ I_s ⊗ 𝒫_{3/4}^{⊗(n−1)}(Tr_s ρ) + (𝒫_{1/2}^{(s)} ⊗ 𝒫_{3/4}^{⊗(n−1)})(ρ),

where 𝒫_{1/2}^{(s)} acts on copy s and 𝒫_{3/4}^{⊗(n−1)} on the others.
(i) The operator inequality 𝒫_{3/4}^{⊗n}(ρ) ⪰ ¼ I_s ⊗ 𝒫_{3/4}^{⊗(n−1)}(Tr_s ρ) for all ρ is positivity of
𝒫_{1/2} ⊗ 𝒫_{3/4}^{⊗(n−1)}. By Proposition 7.20 it is equivalent to real 2-positivity of 𝒫_{3/4}^{⊗(n−1)}, and
tensor-stable positivity is equivalent to it holding for every n. Iterated, it gives
λ_min(𝒫_{3/4}^{⊗n}(ρ)) ≥ 4^{−k} λ_min(𝒫_{3/4}^{⊗(n−k)}(ρ_K)) for every set K of n − k copies.
(ii) Conjecture TC for mixed inputs at n copies follows from Conjecture TC for mixed inputs at n − 1 copies, applied to
Tr_s ρ, together with the bound λ_min((𝒫_{1/2}^{(s)} ⊗ 𝒫_{3/4}^{⊗(n−1)})(ρ)) ≥ ½τ_s ∏_{t≠s}(¼ + ½τ_t) at one suitable
copy s. No thickness is lost, because the one-copy marginals of Tr_s ρ are those of ρ. The two factors add up to
¼ + ½τ_s, which is where the coefficient ½ in Conjecture TC comes from.
(iii) Numerically (adversarial searches over mixed states, floating point, `code/peeling/`): no violation of
Conjecture TC for mixed inputs was found at n = 2, 3, 4. The smallest ratios found are 1, 1 and 1.22; at n = 2, 3 the
searches end at pure products, where the ratio is exactly 1, and at n = 4 the search is weaker. The bound in (ii) holds
at n = 2, 3 when s is chosen well (smallest ratios found 1.090 and 1.239), but fails for a fixed copy (0.947 at n = 2).
The stronger peeling λ_min(𝒫_{3/4}^{⊗n}(ρ)) ≥ (¼ + ½τ_s) λ_min(𝒫_{3/4}^{⊗(n−1)}(Tr_s ρ)) fails (0.907 at n = 2, and
0.954 even for the best copy). So an induction may carry the thickness of the marginals, but not the least eigenvalue
of the marginal's output.
(iv) With every τ_s = 0, the bound in (ii) contains the vertex positivity of (i). Its first open case, 2-positivity of
𝒫_{3/4}^{⊗2}, is not certified by the second level of the hierarchy even when the vertex copy is restricted to a plane.
For inputs and outputs in a plane, 𝒫_{1/2} acts as ½J Sym₂(·) Jᵀ, with J the rotation by a right angle and
Sym₂(X) = (X + Xᵀ)/2. The real
level-2 relaxation of Sym₂ ⊗ 𝒫_{3/4}^{⊗2} has best value −8.3·10⁻³, and that of 𝒫_{1/2} ⊗ 𝒫_{3/4}^{⊗2} on the full qutrit
has −4.6·10⁻³ (the complex value quoted after Proposition 7.20 is −4.7·10⁻³); with one qutrit copy both are 0 to solver
accuracy (`code/peeling/vertex_level2.py`). So the reduction does not make Conjecture TC easier than tensor-stable
positivity. It moves the whole difficulty into the vertex map.

**Theorem 7.32 (six copies at a = 3/4; computer-assisted).** For all x, y ∈ C^{3^6},

    ⟨y|𝒫_{3/4}^{⊗6}(xx*)|y⟩ ≥ 4^{−7} |x|² |y|².

So 𝒫_{3/4}^{⊗6} is positive, and μ₆(3/4) lies in [4^{−7}, 4^{−6}].

*Proof.* Real vectors suffice. For real x, y put ψ = x ⊗ y ⊗ y, with a second copy B′ of the outputs as in
Theorem 7.13. The file `code/dps6/certs/cert_n6_a3_4_t1_16384.json.xz` gives a rational number z_μ for each of the
38 760 multisets μ of six of the fifteen Brauer diagrams of one copy. Let Z = Σ_τ z_{sort(τ)} D_{τ₁} ⊗ ⋯ ⊗ D_{τ₆}, with τ
running over ordered six-tuples. Then

    ⟨ψ|Z|ψ⟩ = ⟨y|𝒫_{3/4}^{⊗6}(xxᵀ)|y⟩ |y|² − 4^{−7} |x|² |y|⁴   for all real x, y,   and   Z ⪰ 0.

The identity is one between polynomials in (x, y). Two operators in the span of products of diagrams give the same
polynomial exactly when their averages over S₂ × S₄ agree, where S₂ permutes the two x slots and S₄ the four y slots,
at every copy at once. The fifteen diagrams are linearly independent on R³, so this is one linear equation for each
orbit of multisets, 1042 in all, and they hold exactly. Z commutes with O(3) at every copy and with permutations of the
copies. Compressing to the m = 0 vectors of each spin type (k₀, k₁, k₂, k₃) and reducing by S_{k₁} × S_{k₂} leaves 217
rational matrices of size at most 48, and each has an exact LDLᵀ factorisation with positive pivots. Two separately
written programs check all of this in exact arithmetic (`code/dps6/README.md`). ∎

This is the second level used for Theorem 7.13, written for real vectors: the swap of the two x slots and the
permutations of the four y slots do the work of Γ_A, Γ_B and P, and numerically the two versions have the same optimum.
At six copies the floating-point optimum is t ≈ 1.3749·10⁻⁴ ≈ 0.563·4^{−6}, below the product value, so unlike n ≤ 5 the
level appears not to be exact. The solver reports this optimum as "optimal_inaccurate", and we have not proved it; the
theorem uses only the exact certificate at t = 4^{−7}. The certificate combines a solution at that optimum with one of largest margin at t = 0 (least block
eigenvalue 1.378·10⁻⁴), rounded to 2^{−50} and then corrected exactly. A seventh copy is out of reach this way: the
same reduction has 2889 orbits there, and the threshold ℓ_n of Remark 7.15 is probably above 3/4 by then.

**Proposition 7.33 (separable inputs, every n).** Let σ = Σ_k p_k v^k v^{k*} be a separable state on (C³)^{⊗n}, with
product vectors v^k = ⊗_s v^k_s, and let τ_s(σ) be the least eigenvalue of its one-copy marginal σ_s. Then

    𝒫_{3/4}^{⊗n}(σ) ⪰ 4^{−n} 3^{2Σ_s τ_s(σ)} I ⪰ ∏_s(¼ + ½τ_s(σ)) I.

So Conjecture TC holds for every separable input, pure or mixed.

*Proof.* For a unit vector v_s, 𝒫_{3/4}(v_sv_s*) = ¾(I − ⅔R_s) with R_s = Re(v_sv_s*), so 0 ⪯ R_s ⪯ I and Tr R_s = 1. On
[0, 1] the convex function 3^{−r} lies below its chord 1 − ⅔r, so I − ⅔R_s ⪰ 3^{−R_s}, and 𝒫_{3/4}^{⊗n}(vv*) ⪰ (¾)^n 3^{−N}
with N = Σ_s R_s (one commuting term per copy). For a unit test vector x with one-copy marginals ρ^x_s, Jensen's
inequality for the spectral measure of N in the state x gives ⟨x|𝒫_{3/4}^{⊗n}(vv*)|x⟩ ≥ (¾)^n 3^{−Σ_s Tr(ρ^x_s R_s)}.
Averaging over k and using Jensen once more, ⟨x|𝒫_{3/4}^{⊗n}(σ)|x⟩ ≥ (¾)^n 3^{−Σ_s Tr(ρ^x_s Re σ_s)}. Finally
Tr(ρ^x_s Re σ_s) ≤ λ_max(Re σ_s) ≤ λ_max(σ_s) ≤ 1 − 2τ_s(σ), and 3^{2τ} ≥ 1 + 2τ. ∎

By self-adjointness this is Proposition 7.30 averaged over a separable decomposition. It does not reach entangled
inputs, and no convexity argument can bring it there. A sufficient condition for TC at an input x would be
𝒫_{3/4}^{⊗n}(xxᵀ) ⪰ 𝒫_{3/4}^{⊗n}(σ) for some separable σ with marginals at least as thick. This fails at the singlet
pair of Observation 7.26: the bottom eigenvector y of 𝒫_{3/4}^{⊗4}(xxᵀ) is balanced, so every separable σ gives
⟨y|𝒫_{3/4}^{⊗4}(σ)|y⟩ ≥ (¾)⁴3^{−4/3} ≈ 0.0731, while ⟨y|𝒫_{3/4}^{⊗4}(xxᵀ)|y⟩ = (151 − √6805)/1536 ≈ 0.0446. The
coherent terms have to be handled with their signs (`code/exact_extremal_values.py`).

## Checks

`code/route_checks.py` (`logs/route_checks.log`), written from the definitions (floating point to 10⁻¹⁰ unless marked
exact):
- the real formulas for Δ and ε and the identity q(P) = (Δ + ε)/2 at n = 1, 2, Δ ≥ 2^{2−n}; fusion at (n₁, n₂) = (1, 1),
  (2, 1), (2, 2); the balanced three-copy code and its six-copy fusion (one-copy compressions zero for all eight
  traceless directions);
- q₃(C ⊗ |0⟩⟨1|) = q₂(C); the GHZ values h_r(t) for r = 2, …, 5 at t = 1/2, 3/5; the exact values h_r(1/2) for
  r ≤ 11; q₅(P ⊗ D₂) = (3/2)(11/12);
- Proposition S at (n, r) = (1, 1), (1, 2), (2, 1), (1, 3);
- Proposition H(i) at N = 2, 3, and the algebra of H(iv) exactly (positive for m ≤ 26 at N ≤ 40, and negative at some
  N ≤ 40 for each m = 27, …, 60);
- for 𝒫_a: the Choi matrix and its spectrum at four values of a; the two-positivity threshold (and 300 random
  Schmidt-rank-two vectors); the Kraus form of R; the measure-and-prepare formula at a = 1 exactly; the twirl and the
  spectrum of ρ = (67I − 25F)/528 and of ρ^Γ; the recursion
  B_{n+1}(x, y) = m Σ_{ij} B_n(x_j, y_i) + ½ Σ_{i<j} [B_n(x_j, y_i) + B_n(x_i, y_j) − 2y_iᵀΦ(x_jx_iᵀ)y_j] at n = 1, 2;
  the flag identity of 4 at n, k ∈ {1, 2} and a ∈ {3/5, 3/4, 9/10, 1}; the product values (a − 1/2)^n.
- Control: for x = |00⟩ + |11⟩, y = |01⟩ − |10⟩ and z = x ⊗ |0⟩ + y ⊗ |1⟩, s₃(z) = ((4a − 1)³ − 4)/4, which is
  −157/500 at a = 3/5 (outside the range of item 2, where 𝒫_a is not 2-positive), so the scalar criterion does detect a
  failing parameter.

Run `bash code/run_all.sh` (Python with NumPy). The numerics of Section 6 are not part of it (about 30 minutes).

`code/map_thresholds.py` checks Section 7, one mode at a time, each ending with PASS (logs in `code/logs/`):
- `witness`: Theorem 7.1 in exact integer arithmetic on the planes for n = 1, …, 8 at six rational values of a, after
  checking that the compression of 𝒫_a^{⊗n} to the planes is 𝒬_a^{⊗n} (n ≤ 4); and, in floating point, that y is a
  bottom eigenvector of the qutrit output for n ≤ 6 (`map_thresholds_witness.log`);
- `two`: Theorem 7.2 exactly (SymPy): the block decomposition, the identity for the first block, the bound for the
  others and its root, the matrix G, and the PPT operator ρ with Tr(C_a^{⊗2}ρ) and its roots (`map_thresholds_two.log`);
- `lp`: the table of Proposition 7.3 for both maps, n ≤ 12 (`map_thresholds_lp.log`);
- `search`: Observation 7.4 for n = 2, …, 5 (`map_thresholds_search.log`, about half an hour); n = 6 (dimension 729)
  was run once with `qutrit_search.py 6 4`, the same procedure (`qutrit_search_n6.log`, about an hour); the inputs with
  full support are in `qutrit_structured.py` (`qutrit_structured.log`);
  the matrix product state searches are `mps_search.py` (`mps_<n>_<D>.log`, arguments in their first lines);
  the general matrix product state searches are `dmrg_search.py n D 12 30 a…` (`dmrg_n<n>_D<D>.log`);
- `planes`: the identity of Remark 7.5 at n = 1, 2 and the least eigenvalue of the Choi matrix of 𝒬′_a
  (`map_thresholds_planes.log`);
- `thin`: Lemma 7.7 at 5 000 random complex planes (the formula for C_V, C_V^Γ = C_V, positivity at a = 3/4, and the
  negative control at a = 0.74 on a real plane), and Proposition 7.8 at n = 3, 4 on random inputs thin at all but two
  copies and by descent inside that class, which reaches 4^{−n} and not below, and the thinning of the output vector
  (x thin only at copy 1, y thin only at copy 2, n = 3: the least value found is 4^{−3}) (`map_thresholds_thin.log`); the numbers after Corollary 7.10 come from `thick_inputs.py 3` (`thick_inputs_n3.log`);
- `collision`: Remark 7.6 at n = 1, 2 with random entangled states (`map_thresholds_collision.log`).

`code/two_positivity/` checks Proposition 7.11 and Remark 7.12 (logs in `code/two_positivity/logs/`, all PASS):
- `two_copies.py`: the identity of Proposition 7.11(i) exactly (SymPy, coefficient by coefficient in a, and as an 81 × 81
  rational matrix identity), the closed forms of w_jk and their signs on [a*, 1];
- `twirls.py`: exactly, the twirls of the product states |e₁e₁⟩, |ww⟩, |ww̄⟩, |e₁e₂⟩ (w = (e₁ + ie₂)/√2) in the basis
  P_Ω, P_s0, P_an;
- `certificates.py n p q`: the certificates of Proposition 7.11(ii) at a = p/q, found by a linear program, rounded to
  rationals and checked exactly (files in `certs/`; runs for (2, 87, 100), (3, 91, 100), (4, 93, 100), (5, 19, 20));
- `matrix_check.py`: a second check, which rebuilds the four separable factors as averages over the 60 rotations of the
  icosahedral group (whose twirl equals the O(3) twirl here, and which gives explicit separable decompositions) and
  verifies the n = 2 and n = 3 certificates as explicit 81 × 81 and 729 × 729 matrices;
- `numerics.py` (with `common.py`): the floating-point observations of Remark 7.12.

`code/dps/` contains the certificates of Theorem 7.13 and their verifiers (details in `code/dps/README.md`):
- `verify_independent.py CERT`: exact verification from the definition of 𝒫_a and the matchings in the file, by
  exact elimination and by an exact congruence test, sector by sector (`logs/independent_*.log`, all PASS);
- `crosscheck3.py CERT`: at three copies, the identity behind the certificate on random entangled complex vectors,
  with the left side computed from the definition of 𝒫_a (agreement to 2·10⁻¹⁶);
- `exact_verify.py`, `exact_verify_fast.py`, `exact_verify_mod.py`: the verifiers that came with the generating
  scripts (exact rational elimination; modular arithmetic with Sylvester's criterion; a memory-light modular version);
- `run_checks.sh` runs all of these, plus negative controls (t raised above 4^{−n}, which must fail), in about ten
  minutes (`logs/run_checks.out`).

`code/dps6/` proves Theorem 7.32 (details in `code/dps6/README.md`). `cert_l2.py 6 3 4 check CERT` checks the orbit
identity exactly and positivity with integer spin bases, Young-symmetrizer blocks and exact LDLᵀ. `verify_l2.py CERT`
is a separately written exact verifier: it rebuilds the target and the orbits from the 27 × 27 matrices, reduces with
Jucys–Murphy eigenbases, and tests the identity at random points of tensor rank two. Both PASS (`code/dps6/logs/`), and
`run_checks.sh` reruns them in about four minutes. `code/exact_extremal_values.py` checks Observation 7.27(d),
Proposition 7.33 at three copies on random separable inputs, and the remark after it
(`logs/exact_extremal_values.log`, PASS).

`code/symmetric_sector/` checks Proposition 7.16 (details in its README): `verify_bosonic.py n` (exact, through
coherent-state symbols and an exact null space for the rotation-invariant slice) and `exact_greedy.py` (exact, through an
integer recursion over the copies) both find the greedy operator positive semidefinite with kernel Sym^{2n} for
n = 2, …, 8 and not positive semidefinite at n = 9. `many_n/run_many_n.sh` verifies the certificates for
n = 9, …, 21 exactly in a coupled angular-momentum basis (`many_n/verify_exact.py`, logs `many_n/logs/`), and
`many_n/verify_bosonic_cert.py` checks the n = 9 certificate a second way. The numerics of the proposition are
`num_min.py`, `t3num.py`, `sdp_o3b.py` and `sosthr.py`.

`code/plane_threshold_check.py` checks Theorem 7.18 (exact in Q(2^{1/n}) for n ≤ 4, explicit matrices for n ≤ 5, a
control above the threshold; `logs/plane_threshold_check.log`).

`code/bond_two_check.py n samples seed` illustrates Proposition 7.17 in floating point (the AKLT chain and random open
matrix product states of bond dimension two and, for comparison, three; `logs/bond_two_check.log`, n = 4, …, 8).

`code/plane_products_check.py` checks Theorem 7.19 (exact identities and p ≥ 0 for rational parameters with product 2,
n ≤ 4; explicit matrices for n ≤ 3; a control; the witness values; signs at random parameters;
`logs/plane_products_check.log`). `code/multiparam_numerics.py` gives the numerics of Conjecture QT_m (modes `two`,
`real 3 5`, `real 4 11`) and Remark 7.21 (mode `complex`), logs `logs/multiparam_*.log`.
`code/plane_witness_hessian.py` checks Proposition 7.22 (full Hessians for n ≤ 6 against the closed form;
`logs/plane_witness_hessian.log`). `code/squeeze_paths.py` checks Proposition 7.23 at minima and at thick saddles for
n = 3, and the compression numerics of Observation 7.25 (`logs/squeeze_paths.log`); `code/thick_critical_points.py`
gives the biased search for thick critical points there (`logs/thick_critical_points.log`). `code/thick_sector.py`
checks Observation 7.26 (exact identities and Hessians), the singlet minima and the tests of Conjecture TC
(`logs/thick_sector.log`); `code/thickness_floor.py` gives the table of Observation 7.27(a)
(`logs/thickness_floor.log`); `code/thickness_bonus.py` checks Proposition 7.28 on a grid of Schmidt
coefficients and searches for violations of Conjecture TC at three and four copies (`logs/thickness_bonus.log`). `code/tc3/` proves Theorem 7.29: the certificate is
`code/tc3/cert/tc3_certificate.json.xz`, `verify_tc3.py` is the separately written exact verifier, `check_float.py` the check
at random points, and `run_checks.sh` reruns everything (`code/tc3/README.md`, logs in `code/tc3/logs/`). `code/planar_minimum.py n starts seed` gives the numerics of Conjecture PM
(`logs/planar_minimum_n3.log`, `logs/planar_minimum_n4.log`). `code/peeling/` checks Proposition 7.30 at random and
locally optimised product test vectors (`product_probe.py`), runs the searches of Remark 7.31 (`peel_search.py`) and
the level-2 relaxations of its part (iv) (`vertex_level2.py`); logs in `code/peeling/logs/`, and `run_checks.sh` reruns
all of them in a few minutes.

Python with NumPy, SciPy (HiGHS) and SymPy; `code/dps/` and `code/symmetric_sector/` also use python-flint, and CVXPY
with Clarabel and SCS for the semidefinite programs.
