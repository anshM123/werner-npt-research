# Tangent planes: an induction that closes with one copy at a time

Ansh Mishra and Aryan Senthilkumar, November 2026.

This folder reduces the endpoint inequality, for every n, on a natural class of operators to a 9 × 9 matrix condition at
a single copy. The class contains every operator whose range lies in the single excitations over a product vector
(`single-excitation-range/`), with any co-range, so it contains the sectors c ≥ 3 that are still open there. The
reduction and the identities below are proved. The 9 × 9 condition (Conjecture SC) is supported by numerics at three,
four and five copies but is **not proved**. Nothing here settles the every-n endpoint.

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

**Conjecture SC.** If C ∈ 𝒯_n (n ≥ 4) is thick (not thin) at copy s, then 𝒬_s − λ𝒩_s is not positive semidefinite on Herm_3
for any λ ∈ [−½, 0).

**Theorem 3.** If Conjecture SC holds for n copies, then EP on 𝒯_{n−1} implies EP on 𝒯_n. EP holds on 𝒯_n for n ≤ 3
(Tomar's three-copy theorem). Consequently, Conjecture SC for every n ≥ 4 gives q_n(C) ≥ 0 for every C ∈ 𝒯_n and every
n. In particular it would give the endpoint for every single-excitation range with any co-range.

*Proof.*
1. Let m = min R over the unit sphere of 𝒯_n. The minimum exists by compactness. Note m ≤ 0, since
   q_n(P_{E_s} ⊗ |0…0⟩⟨0…0|) = q_1(P_{E_s}) · q_{n−1}(·) = 0. Also m ≥ −½, the smallest eigenvalue of Λ^{⊗n} with
   Λ(Y) = Y − ½Tr(Y)I.
2. Suppose m < 0, and let C* attain it.
3. If C* is thin at some copy, Lemma 2 and EP on 𝒯_{n−1} give q_n(C*) ≥ 0, a contradiction.
4. Otherwise C* is thick at every copy. Fix s. Since GC*G ∈ 𝒯_n for every M ∈ PD_3, M = I minimises
   𝒬_s(M)/𝒩_s(M) over PD_3.
5. Put f(ε) = 𝒬_s(I + εH) − m𝒩_s(I + εH) for H ∈ Herm_3. Then f ≥ 0 near 0 and f(0) = 0, so f''(0) ≥ 0. Hence
   𝒬_s − m𝒩_s ⪰ 0 with m ∈ [−½, 0), which contradicts Conjecture SC.

No boundary analysis is needed. ∎

The reduction uses the conjecture only at the value λ = m of a hypothetical negative minimum. Restricting λ to [−½, 0)
uses only m ≤ 0 (an explicit zero) and the spectral bound. Any proof of SC must be sign-free: it may use EP on
𝒯_{n−1}, but not q_n ≥ 0 on 𝒯_n.

### A second route, without induction

The strong form of the conjecture needs no induction and no input from three copies.

**Conjecture SC∞.** If C ∈ 𝒯_n is thick at copy s, there is H ∈ Herm_3 with 𝒩_s(H) = 0 and 𝒬_s(H) < 0.

What the numerics test is that 𝒬_s − λ𝒩_s is positive semidefinite for no real λ (ν_s < 0 below). When 𝒩_s is
indefinite, that gives such an H, by convexity of the joint numerical range of two real quadratic forms in nine
variables (Brickman).

**Theorem 3′.** If SC∞ holds for every n, then q_n ≥ 0 on 𝒯_n for every n.

*Proof.* Take a thick copy s and such an H.
1. H is not semidefinite. Otherwise q(GCG) would be a concave quadratic in ε along M = I + εH, while ‖GCG‖² stays
   positive and affine, so R → −∞, against R ≥ −½.
2. So M = I + εH is positive definite exactly on an interval (ε_−, ε_+) and singular at both ends.
3. On that interval, ‖GCG‖² is affine and positive and q(GCG) is concave, so R is quasi-concave. Hence R(C) is at
   least the smaller of its limits at the two ends.
4. The end operators are thin at s and still lie in 𝒯. An end where GCG = 0 is a proper algebraic condition, and it
   can be avoided by perturbing H inside {𝒩_s = 0}, keeping 𝒬_s(H) < 0.
5. Squeezes at other copies keep a thin copy thin. Repeating at each thick copy, every branch ends, after at most n
   steps, at an operator that is thin at every copy.
6. Such an operator has range in a product of planes, so q ≥ 0 by separability of 2 × 3 PPT states. ∎

## Evidence for Conjecture SC

The numerics test the stronger statement with every λ ∈ [−5, 5]. They report
ν_s = max_λ λ_min(𝒬_s − λ𝒩_s)/‖𝒩_s‖, and ν_s < 0 means the pencil is PSD for no λ. Scripts are in `code/`, outputs in
`logs/`.
- **Random tangent operators** (worst co-range at n = 3, 4; random co-range at n = 5). ν_s < 0 at every copy. Values
  are about −0.1 at n = 5; they approach 0 only where the raw local support is nearly thin.
- **Scale invariance.** For A ∈ GL_3 at copy s, 𝒬_{ACA^*}(H) = 𝒬_C(A^*HA), and likewise for 𝒩
  (`code/check_exact_facts.py`). Only the sign of ν_s is invariant, so the conjecture can be checked on the compact
  slice of balanced operators, where Tr_{−s}(CC^* + C^*C) ∝ I. There the margins do not close near thin points.
  Sampling four regimes (generic, nearly thin U_0, a = 0, range = co-range) gives largest balanced values −0.27, −0.16
  and −0.13 at n = 2, 3, 4 (`code/check_balanced.py`). Independent adversarial searches on the balanced slice, not
  included here, found the same sign.
- **Targeted test.** The general three-copy operator with a single-copy squeeze minimum at every copy
  (`orbit-notes/data/copywise_min_n3.npy`) has only 67% of its range inside its best tangent space. Projecting it into
  that tangent space and maximising min_s ν_s inside the class runs to thin points, with local-support ratios ≤ 3·10⁻⁴
  and ν → 0 from below. The same holds at n = 4, starting from `copywise_min_n4.npy` (41% inside).
- **Exact facts at thin points** (`code/check_exact_facts.py`). At a thin copy, H = |h⟩⟨f| + |f⟩⟨h| (h in the support
  V, f the missing direction) has 𝒩_s(H) = 0 and 𝒬_s(H) = −½ q_{n−1}(⟨h|C|f⟩). The phase average over e^{iθ}|f⟩⟨w| + h.c. gives
  E𝒬 = 2Re q'(C_ff, C_ww) − ½[q'(C_fw) + q'(C_wf)] and E𝒩 = 2Re⟨C_ff, C_ww⟩.

## What does not work

Three things fail, and they show that the conjecture uses the tangent structure.
- **Larger classes.** The one-copy statement fails for every natural enlargement we tried. For ranges with at
  most two excitations over a product it fails as follows. Projecting the operator above into T^(2) and re-optimising inside that class gives a thick operator
  (local-support ratios 0.50, 0.37, 0.35, stored in `data/two_excitation_pencil_psd_n3.npy`) whose pencil is positive
  definite for some λ at every copy (ν_s = +0.037). For general rank-two C, single-copy squeeze minima exist at three
  and four copies (`orbit-notes/`). The known examples have 98–99% of their range inside some T^(n−1) but only
  41–69% inside a tangent space. Planes that merely meet a tangent space also fail: they contain Wu and Zou's
  planes through a product vector and are closed under squeezes and thin-copy removal. Starting from the same
  three-copy operator, re-optimising inside that class reaches thickness ≥ 0.2 with ν_s ≈ +0.02 at every copy.
- **Structured certificates.** Several natural certificates for the conjecture fail adversarially as standalone
  arguments, though the full 9 × 9 pencil still has a negative direction there:
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
- `code/check_squeeze_condition.py`: the evidence for Conjecture SC and the two-excitation contrast. Log:
  `logs/check_squeeze_condition.log`.
- `code/check_exact_facts.py`: scale invariance, the thin-point tilt identity and the phase-averaged tilt formula. Log:
  `logs/check_exact_facts.log`.
- `code/check_balanced.py`: balanced-slice margins. Log: `logs/check_balanced.log`.
- `data/two_excitation_pencil_psd_n3.npy`: the thick two-excitation operator with ν_s > 0 at every copy, written in
  the basis where e = |000⟩.

## Status

- **Proved:** Proposition 1 (identity, checked numerically), Lemma 2 (identity, checked numerically), Theorem 3 (the
  reduction; re-derived by an independent check) and Theorem 3′.
- **Not proved:** Conjectures SC and SC∞, and with them the endpoint on 𝒯_n for n ≥ 4.
- **Not journal peer-reviewed.**
