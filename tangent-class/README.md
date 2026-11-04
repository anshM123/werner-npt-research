# Tangent planes: one-copy reductions of the endpoint, and where they fail

Ansh Mishra and Aryan Senthilkumar, November 2026.

This folder studies the endpoint inequality, for every n, on a natural class of operators: rank-two operators whose
range lies in the tangent space of a product vector. The class contains every operator whose range lies in the single
excitations over a product vector (`single-excitation-range/`), with any co-range, so it contains the sectors c ≥ 3
that are still open there. We prove two reductions of EP on this class to conditions on 9 × 9 matrix pencils at single
copies (Theorems 3 and 3′).

**Correction (4 November 2026).** The per-copy condition we first conjectured (Conjecture SC, and its strong form
SC∞) is **false**. We give thick operators in the class whose pencil at one copy is positive definite at λ = −½: exactly
at n = 3, and in floating point with a wide margin at n = 4, 5 and 6 (`code/check_sc_counterexamples.py`). Earlier
versions of this folder listed numerical support for SC. That support came from random and balanced sampling, and it
missed these operators. The reductions only need the condition at some copy of each operator (Conditions SC′ and SC′∞
below). At every counterexample another copy satisfies it. Adversarial searches have not found an operator at which
every copy fails, but the margins they find go to zero at thin operators. SC′ and SC′∞ are open. Nothing here settles
the every-n endpoint.

Throughout, q_n(C) = Σ_J (−½)^{|J|} ‖Tr_J C‖² and R(C) = q_n(C)/‖C‖² on (C^3)^{⊗n}, and "rank two" means rank at most two.

## The class

Fix a product vector e = e_1 ⊗ ⋯ ⊗ e_n. Its tangent space is

    T_e = span{ e_1 ⊗ ⋯ ⊗ e_{s−1} ⊗ y ⊗ e_{s+1} ⊗ ⋯ ⊗ e_n : s = 1, …, n, y ∈ C^3 },

of dimension 2n + 1: the vector e and all single excitations over it. Let 𝒯_n be the set of rank-two operators whose
range lies in T_e for some product vector e; the co-range is arbitrary. In coordinates with e = |0…0⟩, the range vectors
are u_r = a_r e + Σ_s (U_s e_r)@s, where U_s is a 2 × 2 block mapping a two-dimensional reference space into
E_s = span{|1⟩, |2⟩} at copy s and a_r are vacuum amplitudes. The single-excitation class is a = 0.

Three properties make the class suitable for induction:
- It is invariant under local filters on either side, C ↦ (⊗G_s) C (⊗K_s), since G maps T_e onto T_{Ge}.
- Its unit sphere is compact, since e ↦ T_e is continuous and rank ≤ 2 is a closed condition.
- It is closed under removing a thin copy (Lemma 2).

We call C *thin* at copy s if the local support of its range at s (the support of the copy-s marginal of the range
projector) has dimension at most two. For C ∈ 𝒯_n this means rank U_s ≤ 1, or U_t = 0 for every t ≠ s.

## Proposition 1 (Fock form)

Let D = ⊗_s diag(½, 1, 1), F = Σ_r |r⟩ ⊗ D^{1/2} v_r for X = Σ_r u_r v_r^*, and Ā_s = (U_s^* U_s)^T. Let Z_s act on
F by contracting the reference index r and the letter of copy s with U_s, leaving the vacuum at copy s:
Z_s F = Σ_{r,z} conj(U_s[z, r]) ⟨z|_s f_r. Then, for arbitrary U, a and co-range,

    q_n(X) = ⟨F| Σ_u Ā_u ⊗ (2|0⟩⟨0|_u + Π_{E_u}) |F⟩ − Σ_s ‖Z_s F‖² + ½ ‖Σ_s Z_s F − √2 (⟨a| ⊗ I) F‖².

*Proof.* Write the 2-positivity form as ⟨φ|(id ⊗ Φ^{⊗n})(ψψ^*)|φ⟩ with Φ(Y) = Tr(Y)I − Y/2, ψ = Σ_r |r⟩u_r and
φ = Σ_r |r⟩v_r, and expand Φ^{⊗n}(ψψ^*) = Σ_K (−½)^{|K|} ρ_{RK} ⊗ I_{K^c} over the copies K that keep their letters.
For ψ in R ⊗ T_e, the marginal is ρ_{RK} = |χ_K⟩⟨χ_K| + Ā_{K^c} ⊗ |0_K⟩⟨0_K|, with
χ_K = |a⟩|0_K⟩ + Σ_{s∈K} ψ_s|0_{K∖s}⟩. Summing over K gives D-weighted terms: ⊗_u (I − ½|0⟩⟨0|) = D, and similar
products with one or two factors removed. Conjugating by D^{−1/2} turns them into the one-body term, the star term
−Σ_s‖Z_s F‖², and the square. Its cross term comes from |a⟩|0⟩ against ψ_s. ∎

Two special cases are worth recording.
- **a = 0** (single excitations). The square is ½‖𝒵F‖² with 𝒵 = Σ_s Z_s. With the collective annihilators
  𝔠_r = Σ_s Σ_z conj(U_s[z, r]) |0⟩⟨z|_s, one has [𝔠_r^*, 𝔠_{r'}] = Σ_s |U_s e_r⟩⟨U_s e_{r'}|_s − Σ_s A_s[r', r] |0⟩⟨0|_s.
  This gives the anti-normal-ordered form q_n = ‖F‖² + X(F) − ½‖𝒵F‖² for Σ_s U_s^*U_s = I, where X(F) is the
  swap pairing ⟨(I ⊗ 𝐂)F| SWAP |(I ⊗ 𝐂)F⟩ with 𝐂 = Σ_r |r⟩ ⊗ 𝔠_r^*.
- **a ≠ 0.** The vacuum amplitude only displaces the collective annihilator: 𝔠_r → 𝔠_r − √2 conj(a_r).

`code/check_fock_identity.py` checks the identity against explicit partial traces at n = 2, 3, 4 with random
unnormalised U, a and co-range. The largest relative error is 5·10⁻¹⁶.

## Lemma 2 (thin-copy removal inside the class)

Let C be any rank-two operator whose range has local support at copy s inside a plane V. Then

    q_n(C) = ½ Σ_{k=1}^{6} q_{n−1}(⟨ψ_k|_s C |ψ_k^⊥⟩_s) + Σ_{x ∈ {v_1, v_2}} q_{n−1}(⟨x|_s C |v_3⟩_s).

Here v_1, v_2 is an orthonormal basis of V and v_3 ⊥ V, ψ_k are the six Pauli eigenstates of V in that basis, and ψ_k^⊥
is the state orthogonal to ψ_k inside V. If C ∈ 𝒯_n, every contraction has rank ≤ 2 and range in T_{e'}, where e' is
e with copy s removed. So EP on 𝒯_{n−1} gives q_n(C) ≥ 0 whenever C ∈ 𝒯_n is thin at some copy.

*Proof.* The copy-s kernel of q_n, restricted to range letters in V, is I_{V⊗C^3} − ½|Φ_V⟩⟨Φ_V| with
Φ_V = Σ_j v_j ⊗ v̄_j. On V ⊗ V̄ this is the partial transpose of the two-qubit operator I − F/2. That operator equals
3∫dψ |ψ⟩⟨ψ| ⊗ |ψ^⊥⟩⟨ψ^⊥|, and the six Pauli states form a 3-design, so they reproduce the integral exactly. The rest of
the kernel is I_V ⊗ |v̄_3⟩⟨v̄_3|. Expanding the kernel in these product vectors expresses q_n(C) as a positive
combination of (n−1)-copy values of contractions. Contracting a vector of T_e at copy s gives a vector of T_{e'}. ∎

`code/check_thin_removal.py` checks the identity and the range and rank of every contraction at n = 2, 3, 4. The
largest relative error is 4·10⁻¹⁶.

## Theorem 3 (reduction)

Let C ∈ 𝒯_n and let s be a copy. For M ∈ PD_3 acting on copy s, put G = M^{1/2}. Then M ↦ 𝒬_s(M) := q_n(GCG) and
M ↦ 𝒩_s(M) := ‖GCG‖² are quadratic forms in M on the nine-dimensional space Herm_3:

    𝒬_s(M) = Σ_{a,a',b,b'} M_{aa'} M_{b'b} q_{n−1}(C_{ab}, C_{a'b'}) − ½ q_{n−1}(Tr_s(C M)),    𝒩_s(M) = Tr(C^* M C M),

where C_{ab} = ⟨a|_s C |b⟩_s and q_{n−1}(·, ·) is the polarised form.

**Conjecture SC (false).** We first conjectured that if C ∈ 𝒯_n (n ≥ 4) is thick at copy s, then 𝒬_s − λ𝒩_s is not
positive semidefinite on Herm_3 for any λ ∈ [−½, 0). This fails at n = 4, 5 and 6 (see *Counterexamples* below).

**Condition SC′(n).** No C ∈ 𝒯_n that is thick at every copy has 𝒬_s − λ𝒩_s ⪰ 0 at every copy s for one common
λ ∈ [−½, 0).

SC′ is implied by the false Conjecture SC, and by the weaker statement that every such C has one copy whose pencil is
PSD for no λ ∈ [−½, 0).

**Theorem 3.** If SC′(n) holds, then EP on 𝒯_{n−1} implies EP on 𝒯_n. EP holds on 𝒯_n for n ≤ 3 (Tomar's three-copy
theorem). Consequently, SC′(n) for every n ≥ 4 gives q_n(C) ≥ 0 for every C ∈ 𝒯_n and every n. In particular it would
give the endpoint for every single-excitation range with any co-range.

*Proof.*
1. Let m = min R over the unit sphere of 𝒯_n. The minimum exists by compactness. Note m ≤ 0, since
   q_n(P_{E_s} ⊗ |0…0⟩⟨0…0|) = q_1(P_{E_s}) · q_{n−1}(·) = 0. Also m ≥ −½, the smallest eigenvalue of Λ^{⊗n} with
   Λ(Y) = Y − ½Tr(Y)I.
2. Suppose m < 0, and let C* attain it.
3. If C* is thin at some copy, Lemma 2 and EP on 𝒯_{n−1} give q_n(C*) ≥ 0, a contradiction.
4. Otherwise C* is thick at every copy. Fix s. Since GC*G ∈ 𝒯_n for every M ∈ PD_3, M = I minimises
   𝒬_s(M)/𝒩_s(M) over PD_3.
5. Put f(ε) = 𝒬_s(I + εH) − m𝒩_s(I + εH) for H ∈ Herm_3. Then f ≥ 0 near 0 and f(0) = 0, so f''(0) ≥ 0. Hence
   𝒬_s − m𝒩_s ⪰ 0 at every copy s, with the same m ∈ [−½, 0). This contradicts SC′(n).

No boundary analysis is needed. ∎

The reduction uses the condition only at the value λ = m of a hypothetical negative minimum. Restricting λ to [−½, 0)
uses only m ≤ 0 (an explicit zero) and the spectral bound. Any proof of SC′ must be sign-free: it may use EP on
𝒯_{n−1}, but not q_n ≥ 0 on 𝒯_n. The common value of λ matters: over all rank-two operators, not just 𝒯, the
analogue of SC′ is false (*What does not work*).

### A second route, without induction

A stronger condition needs no induction and no input from three copies.

**Conjecture SC∞ (false).** We first conjectured that if C ∈ 𝒯_n is thick at copy s, there is H ∈ Herm_3 with
𝒩_s(H) = 0 and 𝒬_s(H) < 0. This fails at n = 3, 4, 5 and 6 (see *Counterexamples* below).

**Condition SC′∞(n).** Every C ∈ 𝒯_n that is thick at some copy has a thick copy s and an H ∈ Herm_3 with
𝒩_s(H) = 0 and 𝒬_s(H) < 0.

Such an H exists at copy s exactly when 𝒬_s − λ𝒩_s is positive semidefinite for no real λ (ν_s < 0 below), provided
𝒩_s is indefinite. This follows from convexity of the joint numerical range of two real quadratic forms in nine
variables (Brickman).

**Theorem 3′.** If SC′∞(n) holds for every n, then q_n ≥ 0 on 𝒯_n for every n.

*Proof.* Take a thick copy s with such an H.
1. H is not semidefinite. Otherwise q(GCG) would be a concave quadratic in ε along M = I + εH, while ‖GCG‖² stays
   positive and affine, so R → −∞, against R ≥ −½.
2. So M = I + εH is positive definite exactly on an interval (ε_−, ε_+) and singular at both ends.
3. On that interval, ‖GCG‖² is affine and positive and q(GCG) is concave, so R is quasi-concave. Hence R(C) is at
   least the smaller of its limits at the two ends.
4. The end operators are thin at s and still lie in 𝒯. An end where GCG = 0 is a proper algebraic condition, and it
   can be avoided by perturbing H inside {𝒩_s = 0}, keeping 𝒬_s(H) < 0.
5. Squeezes at other copies keep a thin copy thin. Repeating at a good thick copy each time, every branch ends, after
   at most n steps, at an operator that is thin at every copy.
6. Such an operator has range in a product of planes, so q ≥ 0 by separability of 2 × 3 PPT states. ∎

`code/check_line_iteration.py` runs this construction on random tangent operators (three at n = 3, one at n = 4). It
takes the direction H from a small semidefinite program at each thick copy and follows both ends of every line. Along
all 36 lines, R at the start is at least the smaller end value, with a gap of at least 0.045. All 40 branches end at
operators that are thin at every copy, and the smallest value there is +0.088.

## Numerics

The scripts report ν_s = max_λ λ_min(𝒬_s − λ𝒩_s)/‖𝒩_s‖. ν_s < 0 means the pencil at copy s is PSD for no λ, and
ν_s > 0 means it is positive definite for some λ. Scripts are in `code/`, outputs in `logs/`.

### Counterexamples to SC and SC∞

Take the W-type range u_z = Σ_t |z at copy t⟩ (z = 1, 2), so U_t = I at every copy and a = 0, and choose the co-range
adversarially. `code/check_sc_counterexamples.py` checks four operators C = u_1 g_1^T + u_2 g_2^T. All four have rank
two and are thick at every copy, for range and co-range alike (smallest local-support ratio 0.07–0.22). At copy 0,
𝒬_0 − λ𝒩_0 is positive definite at λ = −½.
- **n = 3.** The co-range is integer:

      g_1 = −|000⟩ − 3|001⟩ − 3|010⟩ + |020⟩ + 3|100⟩ − |101⟩ − |112⟩ − |121⟩ + 2|211⟩,
      g_2 = |000⟩ − 3|002⟩ − 3|020⟩ − 2|122⟩ + 3|200⟩ + |202⟩ + |212⟩ + |221⟩.

  The pencil is positive definite exactly on about (−2.31, −0.36). In exact rational arithmetic, all nine leading
  principal minors of 𝒬_0 − λ𝒩_0 are positive at λ = −1 and λ = −½. R(C) = 0.32.
- **n = 4, 5, 6.** The co-ranges are Gaussian-integer vectors (`data/sc_counterexample_corange_n*.npy`). The values of
  λ_min(𝒬_0 + ½𝒩_0)/‖𝒩_0‖ are +4.2·10⁻³, +5.6·10⁻⁴ and +4.9·10⁻⁴, far above rounding error. The forms agree with
  direct evaluation of q(GCG) and ‖GCG‖² to 10⁻¹⁴. Since n ≥ 4, these refute Conjecture SC as stated.

Such points are common. With the W-type range fixed, maximising ν_0 over the co-range alone reaches ν_0 > 0, at
thickness 0.16–0.20, from 5 of 8 random starts at n = 3. It reaches it from none of 6 at n = 2 (`code/search_corange.py`,
`logs/search_corange.log`). At n = 3 and n = 4, the other copies of the stored operators have no PSD λ
(ν_s ≤ −0.09). So the counterexamples satisfy SC′ and SC′∞.

*Normal form.* At a thick copy s, a filter A = [[α, β^T], [0, B]] (vacuum first) maps 𝒯 to itself. It changes the
pencil at s only by H ↦ A^{−*}HA^{−1}. Together with a change of basis of the range, this allows U_s = I and a = 0.

### The surviving conditions

- **Short adversarial searches** inside 𝒯 maximise min_s ν_s (for SC′∞) or the smallest margin at one common
  λ ∈ [−½, −0.01] (for SC′). Thickness of range and co-range at every copy is held at τ or above. The runs use L-BFGS
  from random, W-type and counterexample starts: 12 runs at n = 3 and 2 at n = 4 (`code/search_surviving_conditions.py`,
  `logs/search_surviving_conditions.log`). Every run ended at the thickness floor, with no violation. At n = 3 the best
  values scale with the floor: about −0.004 at τ = 0.02, −0.01 to −0.06 at τ = 0.05, and −0.05 to −0.10 at τ = 0.1. The
  n = 4 runs ended at −0.099 (τ = 0.05) and −0.154 (τ = 0.1). So, if SC′∞ holds, it holds with no uniform margin.
  These searches are short, and the method is the same kind that missed the counterexamples before, so this is weak
  evidence.
- **Random tangent operators** (worst co-range at n = 3, 4; random co-range at n = 5). ν_s < 0 at every copy. These
  operators do not come close to the counterexamples above.
- **Scale invariance.** For A ∈ GL_3 at copy s, 𝒬_{ACA^*}(H) = 𝒬_C(A^*HA), and likewise for 𝒩
  (`code/check_exact_facts.py`). Only the sign of ν_s is invariant, so sampling can use the compact slice of balanced
  operators, where Tr_{−s}(CC^* + C^*C) ∝ I.
  - Four regimes (generic, nearly thin U_0, a = 0, range = co-range) give largest balanced values −0.27, −0.16 and
    −0.13 at n = 2, 3, 4 (`code/check_balanced.py`).
  - Three samples per regime at n = 5 give −0.12 (generic), −0.12 (a = 0) and −0.064 (range = co-range).
  - In the range = co-range regime the value halves at each step, −0.52, −0.25, −0.13, −0.064 for n = 2, …, 5.
  - Random sampling of this kind did not find the counterexamples above.
- **Targeted test.** The general three-copy operator with a single-copy squeeze minimum at every copy
  (`orbit-notes/data/copywise_min_n3.npy`) has only 67% of its range inside its best tangent space. Projecting it into
  that tangent space and maximising min_s ν_s inside the class runs to thin points, with local-support ratios ≤ 3·10⁻⁴
  and ν → 0 from below. The same holds at n = 4, starting from `copywise_min_n4.npy` (41% inside).
- **Exact facts at thin points** (`code/check_exact_facts.py`). At a thin copy, H = |h⟩⟨f| + |f⟩⟨h| (h in the support
  V, f the missing direction) has 𝒩_s(H) = 0 and 𝒬_s(H) = −½ q_{n−1}(⟨h|C|f⟩). The phase average over e^{iθ}|f⟩⟨w| + h.c. gives
  E𝒬 = 2Re q'(C_ff, C_ww) − ½[q'(C_fw) + q'(C_wf)] and E𝒩 = 2Re⟨C_ff, C_ww⟩.

### A larger class

Theorem 3′ uses only two properties of the class: invariance under single-copy squeezes, and that operators thin at
every copy are products of planes. So any squeeze-invariant class on which SC′∞ holds would do.

One candidate is the bond-two class. Its range planes lie in
span{a_1 ⊗ ⋯ ⊗ a_{s−1} ⊗ y ⊗ b_{s+1} ⊗ ⋯ ⊗ b_n : s, y} for product vectors a and b.
- a = b gives the tangent class.
- a = |0…0⟩, b = |1…1⟩ contains the GHZ plane span{|0…0⟩, |1…1⟩}, a known extremal case.
- At each copy, the range component orthogonal to span{a_s, b_s} has rank one with range a_{<s} ⊗ b_{>s}.

In `code/check_bond_two.py`, balanced margins at random points are at most −0.21 at n = 3 and −0.15 at n = 4. They stay
below −0.17 next to the GHZ plane. Hill-climbs on the balanced margin, not included here, plateau at −0.16 to −0.18 at
n = 3 and −0.14 at n = 4. The class contains 𝒯, so the counterexamples to SC∞ above lie in it as well. Only SC′∞ is
open there.

## What does not work

Four things fail.
- **Larger classes.** In every natural enlargement of 𝒯 we tried, even SC′∞ fails: there are thick operators whose
  pencil is positive definite for some λ at every copy. For ranges with at most two excitations over a product it
  fails as follows. Projecting the operator above into T^(2) and re-optimising inside that class gives a thick operator
  (local-support ratios 0.50, 0.37, 0.35, stored in `data/two_excitation_pencil_psd_n3.npy`) whose pencil is positive
  definite for some λ at every copy (ν_s = +0.037). For general rank-two C, single-copy squeeze minima exist at three
  and four copies (`orbit-notes/`). The known examples have 98–99% of their range inside some T^(n−1) but only
  41–69% inside a tangent space. Planes that merely meet a tangent space also fail: they contain Wu and Zou's
  planes through a product vector and are closed under squeezes and thin-copy removal. Starting from the same
  three-copy operator, re-optimising inside that class reaches thickness ≥ 0.2 with ν_s ≈ +0.02 at every copy.
  Product anchoring alone is not enough either. Take ranges in ⊗_t B_t + Σ_s C φ_s ⊗ f^(s), with planes B_t,
  φ_t ⊥ B_t, and product anchors f^(s) inside the planes. This contains the bond-two class, and it has a thick
  operator with ν_s > 0 at every copy and R = 0.30 (`data/anchored_pencil_psd_n3.npy`). The tangent and bond-two
  ranges span only 2n + 1 dimensions, against 2^n + n here.
- **One shared shift for every copy.** Run over all rank-two operators, the minimiser argument of Theorem 3 has more
  to work with: at a thick negative minimiser the pencils at all copies are PSD at the same λ = m. The sign-free
  version, "no thick operator has all its pencils PSD at one shared λ ∈ [−½, 0)", is false.
  `data/common_shift_pd_n3.npy` is a thick three-copy operator (local-support ratios of range and co-range ≥ 0.30)
  whose three pencils are all positive definite at λ = −½. They are PSD on [−2.39, −0.038], [−5.11, +0.050] and
  [−7.00, −0.026]. Its value is R = 0.265, so what the argument really uses is the coupling λ = R(C), and for
  negative values that is equivalent to EP itself (`code/check_common_shift.py`). The same failure occurs away from
  the endpoint. At Werner parameter t, the argument would need the shift to avoid [1 − 3t, 0). At t = 0.45, 0.40,
  0.35 and 0.335 there are thick three-copy operators whose three pencils are all positive definite at λ = 1 − 3t
  (`code/check_common_shift_t.py`, `data/common_shift_t*_n3.npy`). So, for general rank-two operators, the one-copy
  minimiser argument fails at every t in (1/3, ½] that we tried.
- **Structured certificates.** Several natural certificates for the one-copy condition fail adversarially as
  standalone arguments, at points where the full 9 × 9 pencil still has a negative direction:
  - tilts of the vacuum alone;
  - the vacuum squeeze diag(t, 1, 1) alone;
  - compressing along the weakest direction of U_s;
  - the U(1)-twirled certificate built from the traceless excitation block and phase-averaged tilts.
- **Global certificates for the single-excitation sectors.** None of these closes the sectors c ≥ 3 of
  `single-excitation-range/`:
  - a decomposable (PSD plus partial transpose) certificate for the biquadratic form, which fails already at n = 2;
  - block-norm Perron bounds;
  - dual star (monogamy) bounds;
  - a scalar anisotropic central-spin bound.

## Files

- `code/common.py`: q_n, tangent-class operators, and the single-copy squeeze forms.
- `code/check_fock_identity.py`: Proposition 1. Log: `logs/check_fock_identity.log`.
- `code/check_thin_removal.py`: Lemma 2. Log: `logs/check_thin_removal.log`.
- `code/check_squeeze_condition.py`: random sampling of ν_s and the two-excitation contrast. Log:
  `logs/check_squeeze_condition.log`.
- `code/check_sc_counterexamples.py` with `data/sc_counterexample_corange_n{4,5,6}.npy`: the counterexamples to SC and
  SC∞ (n = 3 exact, needs sympy). Log: `logs/check_sc_counterexamples.log`.
- `code/search_corange.py`, `code/search_surviving_conditions.py`: the adversarial searches above (W-type range with free
  co-range; SC′ and SC′∞ inside 𝒯). Logs: `logs/search_corange.log`, `logs/search_surviving_conditions.log`.
- `code/check_exact_facts.py`: scale invariance, the thin-point tilt identity and the phase-averaged tilt formula. Log:
  `logs/check_exact_facts.log`.
- `code/check_balanced.py`: balanced-slice margins. Log: `logs/check_balanced.log`.
- `data/two_excitation_pencil_psd_n3.npy`: the thick two-excitation operator with ν_s > 0 at every copy, written in
  the basis where e = |000⟩.
- `code/check_bond_two.py`: the bond-two class and the anchored contrast. Log: `logs/check_bond_two.log`.
- `data/anchored_pencil_psd_n3.npy`: the thick anchored operator with ν_s > 0 at every copy.
- `code/check_line_iteration.py`: the construction of Theorem 3′ followed on examples. Log:
  `logs/check_line_iteration.log`.
- `code/check_common_shift.py` with `data/common_shift_pd_n3.npy`: a general thick operator whose pencils are all
  positive definite at λ = −½. Log: `logs/check_common_shift.log`.
- `code/check_common_shift_t.py` with `data/common_shift_t*_n3.npy`: the same at Werner parameters t = 0.45, 0.40,
  0.35 and 0.335. Log: `logs/check_common_shift_t.log`.

## Status

- **Proved:** Proposition 1 (identity, checked numerically), Lemma 2 (identity, checked numerically), Theorem 3 (the
  reduction; re-derived by an independent check) and Theorem 3′.
- **False:** Conjectures SC and SC∞ (explicit counterexamples, exact at n = 3).
- **Open:** Conditions SC′ and SC′∞, and with them the endpoint on 𝒯_n for n ≥ 4.
- **Not journal peer-reviewed.**
