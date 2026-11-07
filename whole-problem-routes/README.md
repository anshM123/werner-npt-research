# Whole-problem routes: exact reformulations of the endpoint

Ansh Mishra and Aryan Senthilkumar, 6 November 2026. Nothing here is refereed.

By `support-and-projector-notes/` and `real-balanced-reduction/`, the endpoint for every number of copies is equivalent to

    q_N(P) ≥ 0 for every odd N and every real rank-two orthogonal projector P with balanced range,

where q_N(C) = Σ_J (−1/2)^{|J|} ‖Tr_J C‖² and b_N(X, Y) = Σ_J (−1/2)^{|J|} ⟨Tr_J X, Tr_J Y⟩. This note records five
exact reformulations of that whole statement. Four restate the Werner endpoint. The fifth is a sufficient condition for
NPT bound entanglement that goes through a different state. **None of them is proved here.** Each reduces the problem to
one inequality, stated below, that would settle it.

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
