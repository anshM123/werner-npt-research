# Tangent planes as qubit POVMs: a minimum principle

Ansh Mishra and Aryan Senthilkumar, 4 November 2026.

This note rewrites the endpoint inequality on single-excitation and tangent range planes, with any co-range, as a
statement about one function on qubit POVMs. Its extreme points are products of planes, where the inequality is
known. So every n would follow from a minimum principle: the function has no local minimum inside a face of the POVM
set. The reformulation and the sector identities below are proved (and checked numerically). The minimum principle is
a **conjecture**. Descents and adversarial searches at n = 2–5 found no violation. Nothing here settles the every-n
endpoint.

Throughout, q_n(C) = Σ_J (−½)^{|J|} ‖Tr_J C‖², R(C) = q_n(C)/‖C‖², "rank two" means rank at most two, and
Ψ(Y) = Tr(Y) I − Y/2 on one copy. Then q_n(C) = ⟨V|(id ⊗ Ψ^{⊗n})(|𝐮⟩⟨𝐮|)|V⟩ for C = Σ_r u_r v_r^*, with
|𝐮⟩ = Σ_r |r⟩u_r and |V⟩ = Σ_r |r⟩v_r.

## 1. The POVM picture

Fix e = |0…0⟩ and E_s = span{|1⟩, |2⟩} at copy s. A tangent range plane is spanned by

    u_r = a_r e + Σ_s (U_s e_r)@s,      r = 1, 2,

with U_s : C² → E_s and vacuum amplitudes a ∈ C² (`tangent-class/`). The single-excitation planes of
`single-excitation-range/` are the case a = 0. A change of basis of the range plane lets us assume the gauge

    |ā⟩⟨ā| + Σ_s A_s = I,      A_s = U_s^* U_s.

Then ⟨u_r|u_r'⟩ = δ_rr', so ‖C‖² = ‖V‖², and the minimum over all co-ranges is an eigenvalue:

    m := min over co-ranges of R(C) = λ_min K,      K = [Ψ^{⊗n}(u_r u_r'^*)]_{r,r'}.

**Proposition 1.** m depends only on the tuple (A_0 = |ā⟩⟨ā|, A_1, …, A_n), up to simultaneous unitary conjugation.
This tuple is a qubit POVM whose vacuum element A_0 has rank at most one. For single-excitation planes it is an
arbitrary n-outcome qubit POVM.

*Proof.* Write U_s = W_s A_s^{1/2} with W_s unitary on E_s. The local unitary ⊗_s(|0⟩⟨0| + W_s) commutes with Ψ^{⊗n}
and is absorbed by the co-range, so m does not see W_s. The phase of a is a local unitary of the same kind. A
remaining unitary change of range basis conjugates every A_s simultaneously. ∎

**Proposition 2 (extreme points).** If every A_s (s ≥ 1) has rank at most one, or the POVM is (I, 0, …, 0), then m ≥ 0.
Every extreme point of the set of POVMs (for a fixed vacuum direction) is of this kind.

*Proof.* Then the range has local support span{|0⟩, x_s} at each copy, so it lies in a product of planes, and
q_n ≥ 0 by separability of 2 × 3 PPT states. A qubit POVM with an element of rank two is extreme only if that element
is I. ∎

So the endpoint on single-excitation planes, for one n, holds if m attains its minimum over the (compact, convex) POVM
set at an extreme point. The same holds on tangent planes, over the convex set P_ψ of POVMs whose vacuum element is a
multiple of a fixed rank-one projector |ψ⟩⟨ψ|.

## 2. Sectors

For single-excitation planes, K preserves the number c of excitations of co-range strings. On sector c,
2^{n−c} K equals

    K̃_c = I + Σ_{s vacuum} Ā_s − ½ Σ_{s excited} |Ũ_s⟩⟨Ũ_s| + ½ Σ_{t excited, s vacuum} Z_s^* Z_t,

with |Ũ_s⟩ = Σ_{r,z} U_s[z, r] |r⟩|z⟩_s and Z_t = Σ_{r,z} conj(U_t[z, r]) ⟨r| ⊗ |0⟩⟨z|_t. This is the sector formula of
`single-excitation-range/`, written in the POVM gauge.

**Proposition 3.**
1. K̃_0 = 2I, so sector 0 contributes the constant 2^{1−n} and m ≤ 2^{1−n} everywhere.
2. K̃_n ⪰ 0, since ‖Σ_s |Ũ_s⟩⟨Ũ_s|‖ ≤ Σ_s Tr A_s = 2.
3. K̃_c ⪰ 2I − 3 Σ_{t excited} Ā_t. So sector c is nonnegative whenever Σ_{t∈σ} A_t ≤ (2/3) I for every set σ of c
   copies.

*Proof of 3.* Each Ā_t ⊗ I − ½|Ũ_t⟩⟨Ũ_t| is positive semidefinite, by Cauchy–Schwarz with |Ũ_t⟩ = (Ā_t^{1/2} ⊗ I)|Ω⟩.
The hopping term is ½Σ_ρ [|Σ_t τ_t(ρ)|² − Σ_t |τ_t(ρ)|²], where τ_t(ρ) is Z_t F evaluated at the string ρ. Drop the
square and use Σ_{s vacuum} Ā_s = I − Σ_{t excited} Ā_t. ∎

Sectors c ≤ 2 are Theorem 1 of `single-excitation-range/`. So on single-excitation planes, m ≥ 0 reduces to
μ_c := λ_min K̃_c ≥ 0 for 3 ≤ c ≤ n − 1. The first open case is n = 4, c = 3.

`code/check_reformulation.py` checks, at n = 2, 3, 4:
- m = λ_min K against the direct minimum over co-ranges, for single-excitation and tangent planes;
- invariance under local unitaries and under simultaneous conjugation;
- the explicit sector formula (to 3·10⁻¹⁵), and that blocks between different sectors vanish;
- items 1–3 of Proposition 3, and m ≥ 0 at extreme POVMs.

Log: `logs/check_reformulation.log`.

## 3. The conjecture

**Conjecture MP.** For every n and every 1 ≤ c ≤ n, the sector function μ_c has no local minimum in the relative
interior of a face of positive dimension of the POVM set, at a value below its infimum over the relative boundary of
that face.

By Propositions 2 and 3 and induction over faces, MP gives μ_c ≥ 0 for every c, hence the endpoint for every
single-excitation range plane with any co-range, for every n. Only the open sectors 3 ≤ c ≤ n − 1 are needed; the
statement for the other sectors is included because the numerics support it there too. That includes every sector c ≥ 3 left open in
`single-excitation-range/`. The analogue on the convex sets P_ψ would give it for every tangent plane.

MP is stated sector by sector because m itself has flat pieces. Where the constant sector 0 is the smallest, m equals
its maximum value 2^{1−n} on an open set. At n = 5, two of five random descents started on such a plateau and stayed
there (`logs/descents.log`). These plateaus sit at the top value of m and play no role in its minimum.

MP is a minimum principle on a single group orbit. On thick single-excitation planes, the range-side filters
⊗_s(|0⟩⟨0| + B_s), B_s ∈ GL(2), act transitively. Their balanced (Kempf–Ness) points are the weighted W-planes. There
Theorem 2 of `single-excitation-range/` proves the endpoint, for weights within a factor two of each other. MP is in
the spirit of the orbit property of `orbit-notes/` (no local minima on filter orbits), but with the co-range optimised
first and one sector at a time. We do not claim that either statement implies the other.

## 4. What does not hold

Natural stronger statements fail, so a proof of MP has to use more than convexity:
- **Concavity of m** on the POVM set fails. At n = 2 the function is φ(a, b), with (a, b) the eigenvalues of A_1. φ has
  positive Hessian eigenvalues at a quarter of interior grid points, near a corner.
- **Quasi-concavity** fails. At n = 2 the segment (0.1313, 0.1234) → (0.2617, 0.0501) dips from 0.0909 and 0.0917 to
  0.0902 in between (`code/check_n2_landscape.py`).
- **Quasi-concavity along one-parameter squeeze subgroups** U_s ↦ e^{τK_s} U_s fails at n = 2 (2 of 80 random
  curves). It held on all tested curves at n = 3 (40) and n = 4 (12). Sector by sector
  (`code/check_sector_geodesics.py`), the top sector μ_n dips on every tested curve, by up to 0.36, and μ_2 dips once
  in 10 curves at n = 4. μ_1 at n = 3, 4 and μ_3 at n = 4 never dip. So no convexity along curves is available in
  general, even though, as below, no sector shows an interior local minimum.
- **Fixed co-range.** For a fixed co-range, A ↦ q_n(C(A, V)) is not concave, and is not quasi-concave along squeeze
  subgroups (about half of random curves dip). So the co-range has to be optimised first.
- **One-copy squeeze conditions.** On tangent planes, every one-copy squeeze condition we formulated is false
  (`tangent-class/`).

## 5. Evidence for MP

- **n = 2.** φ vanishes at the four corners and equals 1/4 at the centre, the W-plane. The centre is the only interior
  critical point found, and it is a maximum. All 40 local minimisations from interior starts end at a corner
  (`logs/check_n2_landscape.log`).
- **Descents** (`code/search_minimum_principle.py`, L-BFGS from random interior points, `logs/descents.log`):
  - single-excitation planes: 30 of 30 runs at n = 3 and 16 of 16 at n = 4 end at an extreme POVM with m = 0. At
    n = 5, three of five end on the boundary with m = 0, and two stay on the plateau m = 2^{1−n};
  - tangent planes with a vacuum amplitude: 16 of 16 at n = 3 and 8 of 8 at n = 4 end on the boundary (some A_s
    singular) with m = 0;
  - the open sector c = 3 at n = 4: 12 of 12 runs end at an extreme POVM with μ_3 = 0; c = 2 at n = 4: 6 of 6; c = 3 at
    n = 5: 3 of 3, and 3 of 3 more in `logs/descent_n5_c3.log`; c = 4 at n = 5: 3 of 3 (`logs/descent_n5_c4.log`).
    The other sectors (c = 1 at n = 3, 4; c = 2 at n = 4; the top sector c = n at n = 3, 4) behave the same: 24 of 24
    runs end on the boundary with value 0 (`logs/descent_sectors_n3_n4.log`);
  - on a face: with A_1 = w|ψ⟩⟨ψ| for a fixed ψ, the descents of μ_3 (n = 4, 5) and μ_4 (n = 5) all end at an
    extreme point of that face, every other element of rank one, with value 0 (14 of 14; `code/search_face.py`,
    `logs/face_descents.log`).
- **Starts at hard points.** The tangent-class operators that are squeeze minima at every copy
  (`tangent-class/data/copywise_squeeze_min_n*.npy`) have m = 0.097 (n = 3) and 0.098 (n = 4) once the co-range is
  optimised. All ten descents from them end at m = 0.
- **Adversary for strict local minima** (`logs/adversary.log`). It maximises
  φ_h(A) = min_Δ [f(A + hΔ) − f(A)]/h over interior POVMs, where Δ runs over 2(4n − 4) + 40 tangent directions, and
  φ_h > 0 would mean a strict local minimum. Results:
  - n = 2: −3.9·10⁻³, at the W-plane;
  - n = 3: −0.13 to −0.18 over 12 starts and two scales h;
  - n = 4: 0 for m, which is the plateau m = 2^{1−n}, and −0.127 to −0.157 over 3 starts for the sector function μ_3.

  The search uses Nelder–Mead in 8n dimensions, which is weak, so these numbers are evidence only.

## Files

- `code/povm_form.py`: the form, the POVM gauge, K, m, the sector functions μ_c and parametrisations.
- `code/check_reformulation.py`: Propositions 1–3 numerically. Log: `logs/check_reformulation.log`.
- `code/check_n2_landscape.py`: the n = 2 landscape (no concavity, no quasi-concavity, no interior minimum). Log:
  `logs/check_n2_landscape.log`.
- `code/search_minimum_principle.py`: descents and the φ_h adversary. Logs: `logs/descents.log`, `logs/adversary.log`,
  `logs/descent_n5_c3.log`, `logs/descent_n5_c4.log`, `logs/descent_sectors_n3_n4.log`.
- `code/search_face.py`: descents on the face A_1 = w|ψ⟩⟨ψ|. Log: `logs/face_descents.log`.
- `code/check_sector_geodesics.py`: sector functions along squeeze subgroups. Log: `logs/check_sector_geodesics.log`.

## Status

- **Proved:** Propositions 1–3 (identities, also checked numerically).
- **Conjecture:** MP. It would give the endpoint on all single-excitation range planes with any co-range for every n,
  and on tangent planes with its analogue on P_ψ.
- **Not proved, not journal peer-reviewed.**
