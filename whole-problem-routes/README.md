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
entanglement. It is open from three copies on. At the parameter that matters, a = 3/4, inputs supported on a plane at a
copy are entanglement broken there, so every copy where the input is thin can be removed (Proposition 7.8). Only inputs
with full support at every copy, on at least three copies, remain.

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
It does not need the Werner endpoint itself.

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

The thresholds increase toward 1. So from n = 3 on, 𝒫_{3/4}^{⊗n} is not decomposable, and a proof of tensor-stable
positivity at any a < 1 cannot go through decomposability of the tensor powers. Nor can it go through the step of
Proposition M(4), which gives positivity of 𝒫_a^{⊗n} from 2-positivity of 𝒫_a^{⊗(n−1)}. The program with a qubit
ancilla shows that 𝒫_a^{⊗k} ⊗ id₂ is decomposable only at a = 1, for every k ≤ 6. So that 2-positivity is never
decomposable either.

On a plane the picture is different. The same program for 𝒬_a, with G₂ = [[½, 1, −½], [½, 0, ½], [−½, 1, ½]] and
c_C = (a − ¾, a − ¼, a + ¼)^{⊗n}, gives decomposability thresholds equal to a_n to eight digits for every n ≤ 12.
Together with Theorem 7.1 this says that, for n ≤ 12, 𝒬_a^{⊗n} is positive exactly when a ≥ a_n (computer-assisted,
floating point; exact certificates are not written out here).

**Observation 7.4 (numerical).** For n = 2, …, 6 the positivity threshold of 𝒫_a^{⊗n} found by search is a_n. At
a_n − 0.001 the search finds negative values, and at a_n + 0.001 it finds none. The searches use alternating minimisation
seeded at the vector of Theorem 7.1, at random rotations and perturbations of it in all three directions of every copy,
at random starts, and by continuation in a. Every negative minimiser found has one-copy rank 2, so it lives on one plane
per copy. Bisections without these seeds agree to five digits for n = 2, 3. We also tried inputs with full support at
every copy: the qutrit GHZ state, ε_{ijk}, the AME(4,3) state, the GHZ state in the Fourier basis with its real and
imaginary parts, and tensor products of these with each other and with the witness of Theorem 7.1. For n = 3, …, 6 their
least output eigenvalues stay well above zero at a_n + 0.001 and at 3/4. Descents started from them (n ≤ 5) end at the
product value or on the plane branch.

**Conjecture QT.** For every n, 𝒫_a^{⊗n} is positive if and only if a ≥ a_n.

Since a_n < 3/4 for every n, Conjecture QT makes 𝒫_{3/4} tensor-stable positive. By Proposition M(3), the Werner state
proportional to 11I − 5F would then be NPT and undistillable from every number of copies. In words, the conjecture says
that the only obstruction is the one on planes, where 𝒫_{3/4} is completely positive. Its case n = 2 is Theorem 7.2. It
is open for every n ≥ 3, where, by Proposition 7.3, the positivity it asserts at a = 3/4 is not decomposable. At
a = 3/4, Corollary 7.9 below reduces it to inputs that have full support at every copy.

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
at which the one-copy reduced state of ρ has rank 3. If 𝒫_a^{⊗|K|} is positive, then 𝒫_a^{⊗n}(ρ) ⪰ 0. By Theorem 7.2 this
holds for every n whenever |K| ≤ 2.

*Proof.* For s ∉ K let V_s be an isometry onto a plane that contains the support of ρ_s. Then ρ = Wρ′W* with
W = ⊗_{s∉K} V_s ⊗ I_K and ρ′ ⪰ 0, so 𝒫_a^{⊗n}(ρ) = (⊗_{s∉K} Φ_{V_s} ⊗ 𝒫_a^{⊗|K|})(ρ′). A tensor product of entanglement
breaking maps is entanglement breaking, so it has the form X ↦ Σ_k Tr(F_kX) R_k with F_k, R_k ⪰ 0, and then
(Σ_k Tr(F_k ·) R_k ⊗ Ψ)(ρ′) = Σ_k R_k ⊗ Ψ(Tr₁[(F_k ⊗ I)ρ′]) ⪰ 0 for every positive map Ψ. ∎

**Corollary 7.9.** 𝒫_{3/4} is tensor-stable positive, and then NPT bound entanglement exists, if and only if for every n
the operator 𝒫_{3/4}^{⊗n}(ρ) is positive for every state ρ whose one-copy reduced states all have rank 3. A smallest
failure, if there is one, needs n ≥ 3 copies and an input that is thick (rank 3) at every copy.

*Proof.* Induction on n. An input with a thin copy is covered by Proposition 7.8 and the case of fewer copies; n ≤ 2 is
Theorem 7.2. ∎

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
- `planes`: the identity of Remark 7.5 at n = 1, 2 and the least eigenvalue of the Choi matrix of 𝒬′_a
  (`map_thresholds_planes.log`);
- `thin`: Lemma 7.7 at 5 000 random complex planes (the formula for C_V, C_V^Γ = C_V, positivity at a = 3/4, and the
  negative control at a = 0.74 on a real plane), and Proposition 7.8 on random inputs thin at all but two copies, n = 3, 4
  (`map_thresholds_thin.log`); the numbers after Corollary 7.10 come from `thick_inputs.py 3` (`thick_inputs_n3.log`);
- `collision`: Remark 7.6 at n = 1, 2 with random entangled states (`map_thresholds_collision.log`).

Python with NumPy, SciPy (HiGHS) and SymPy.
