# The letter-1 block of (R∞)_{k−2}: exact reduction, soft modes, sharpness

Ansh Mishra and Aryan Senthilkumar, 4 November 2026. Nothing here is refereed. Notation as in
`relative-anchors/README.md` (Sections 1–3): s_t = w_t², |s| = Σ_t s_t, x_t = s_t/|s|, all w_t > 0, k ≥ 4.

## Status

| Statement | Status |
|---|---|
| Goal (Proposition K2⁺): the letter-1 block of (R∞)_{k−2} is ⪰ 0 with kernel exactly p ⊗ (letter-1 H_{k−2}) ⊕ span(vac ⊗ J⁻\|1…1⟩), every k ≥ 4, all weights | **not proved** (also not for k = 4 at general weights) |
| Lemma R1 (exact elimination of the vacuum row): Goal ⟺ F ≻ 0 on the complement of p ⊗ H_{k−2}, F = Y11 − ½Y10Y01 explicit | **proved**, all k; exact check k = 4, 5, 6 |
| Lemma J (commutation): 𝒯_r J⁻ = E σ_r⁻ − σ_r⁻ C + 2w_r J⁻ \|0⟩⟨0\|_r on letter-1 strings | **proved**, all k; exact check (as matrices H_{k−1} → H_{k−2}) k = 4, 5, 6 |
| Lemma W: E_{k−2}^{−1} = 𝒩^{−1} − ½𝒩^{−1}J⁻M^{−1}J⁺𝒩^{−1} with the k × k Cauchy-type matrix M; E_{k−2}^{−1}J⁻ = 𝒩^{−1}J⁻M^{−1} | **proved**, all k; exact check k = 4, 5, 6 |
| Lemma BJ: the B-term on the J⁻ part is Σ_{h,g} ⟨G_h X_{·g}, E^{−1} G_g X_{·h}⟩ | **proved**, all k; exact check k = 4, 5, 6 |
| **Proposition S** (soft modes, exact Rayleigh form): F(Ψ(a,h)) = 2s_h\|a\|² − 4s_h²⟨u, E^{−1}u⟩ − ½\|c\|² with explicit c = O(s_h) | **proved**, all k; exact check k = 4, 5, 6 |
| **Corollary S** (sharpness): on the complement of its kernel the block has smallest eigenvalue ≤ min_h s_h/(2(\|s\| − s_h)) | **proved**, all k |
| The block is ⪰ 0 with exactly the claimed kernel | numerically true: k = 4, 5, 6, 7, 610 weight vectors (spreads up to 10^{±3} in w), `sweep.log` |
| Hard part: F ⪰ 0.45 I on R = (S ⊕ p ⊗ H_{k−2})^⊥ | numerical only (inf ≈ 0.4506 over the sweep, units of 4Y) |

So what we have is an exact reduction and an exact description of the soft directions, not a proof.
What it pins down: the gap of the block is at most min_h s_h/(2(|s| − s_h)) (a proof), and numerically at least about
0.24 min_h x_h. So every proof has to be sharp to first order in the smallest weight. The soft directions are known
exactly: they are the J⁻-lifts Ψ(a, h) of the Theorem T vectors of the (k − 1)-copy subgroups, one family of dimension
k − 2 for each copy h. Their Rayleigh form is explicit, and everything nonlocal goes through one k × k Cauchy-type
matrix M.

## 1. The block and its reduced form

Letter-1 column strings are labelled by their sets of holes (zeros). Level k − 2 consists of the pairs P = {h, h′},
level k − 1 of the single holes h. The block has the states P_{r,P} = (1@r, holes P), r ∈ [k], k·C(k,2) of them, and
Q_h = (vac, hole h), k of them. Write Ψ = Σ_r e_r ⊗ ψ_r + vac ⊗ φ, with ψ_r ∈ R^{pairs}, φ ∈ R^k, and X[r, P] = ψ_r(P).
All forms below are in units of 2^Z Y = 4Y (Lemma B, Z = 2).

On letter-1 strings (Lemma N):
- 𝒯_r|P⟩ = ½w_r|P⟩ + ½Σ_{t∈P} w_t|P − t + r⟩ if r ∉ P, and 𝒯_r|P⟩ = 2w_r|P⟩ if r ∈ P;
- 𝒯_r|h⟩ = ½w_h|r⟩ + ½w_r|h⟩ if h ≠ r, and 𝒯_r|r⟩ = 2w_r|r⟩;
- E := E_{k−2} = 𝒩 + ½J⁻J⁺ with 𝒩|P⟩ = (½|s| + s_P)|P⟩;
- E_1 := E_{k−1} = diag(½|s| + s_h) + ½ww^T;
- J⁻|h⟩ = Σ_{t≠h} w_t|{h,t}⟩, σ_r⁻|h⟩ = |{h,r}⟩ (0 if r = h), C = 2s_h on |h⟩ and 2s_P on |P⟩.

**The A-part.** It is ½Σ_{r, P∌r} X[r,P]² + Σ_x [(3/2)|Y_x|² + ½(Σ_a Y_x[a])²], where Y_x[a] := X[a, {a,x}] for
a ≠ x. This is the swap term of Lemma B plus 2 − n_r.

**The B-part.** It is −Σ_{a,b} ⟨𝒯_bψ_a, E^{−1}𝒯_aψ_b⟩.

**The coupling.** It is −2⟨c(X), φ⟩ with c(X) = Σ_r (σ_r⁺ − 𝒯_r^T E_1^{−1} J⁺) ψ_r ∈ R^k.

**The vacuum-row block.** It is 2(I − ww^T/|s|) (Lemma G, E_k = ½|s| on 1…1).

`mid.py` builds this block for any number type. `check_closed_form.py` and `sweep.py` compare it with the direct
construction from M, V_p, D (`single-excitation-range/code/check_bath_reduction.group_ops`), with deviation at most
2.2 × 10^{−16} for k = 4, 5. Both scripts print PASS.

**Lemma R1 (vacuum row eliminated).** Let F(X) := ⟨X, Y11 X⟩ − ½|c(X)|², that is, F = Y11 − ½Y10Y01. Then
c(X) ⊥ w for every X, and for every Ψ = (X, φ), with φ_⊥ the part of φ orthogonal to w,

    ⟨Ψ, YΨ⟩ = F(X) + 2|φ_⊥ − ½c(X)|².

Hence the block is ⪰ 0 with kernel exactly p ⊗ H_{k−2} ⊕ span(vac ⊗ J⁻|1…1⟩) iff F ⪰ 0 with kernel exactly
p ⊗ H_{k−2} = {X[r,P] = w_r g(P)}, that is, iff F ≻ 0 on {Σ_r w_r ψ_r = 0}.

*Proof.* Proposition K2 says (0, w) is in the kernel, so Y10 w = 0 and Y00 w = 0. Hence c(X) = −Y01X ⊥ w, and
Y00 = 2Π_{w⊥}. Complete the square in φ_⊥. If (X, φ) is in the kernel, then F(X) = 0 and φ_⊥ = ½c(X). For
X = p ⊗ g we have c(X) = 0, since (p ⊗ g, 0) is in the kernel by Lemma N. ∎

## 2. Exact structure

**Lemma J (commutation).** On letter-1 strings, for every r,
𝒯_r J⁻ = E σ_r⁻ − σ_r⁻ C + 2w_r J⁻ |0⟩⟨0|_r. Here E and C act at the lower level. On H_{k−1}:
𝒯_r J⁻|h⟩ = E|{h,r}⟩ − 2s_h|{h,r}⟩ + 2w_r δ_{rh} J⁻|h⟩.

*Proof.* We have 𝒯_r = ½J⁺σ_r⁻ + 2w_r|0⟩⟨0|_r, and σ_r⁻ commutes with J⁻. Next, ½J⁺J⁻ = E − C on the lower level, and
Cσ_r⁻ = σ_r⁻(C + 2s_r). Finally |0⟩⟨0|_r J⁻ = J⁻|0⟩⟨0|_r + w_rσ_r⁻. The two terms ±2s_rσ_r⁻ cancel. ∎

**Lemma W (the only nonlocal object).** Put M := I + ½J⁺𝒩^{−1}J⁻ on H_{k−1}. Then
- M_hg = ½w_hw_g/(½|s| + s_h + s_g) for h ≠ g, and M_hh = 1 + ½Σ_{t≠h} s_t/(½|s| + s_h + s_t);
- E^{−1} = 𝒩^{−1} − ½𝒩^{−1}J⁻M^{−1}J⁺𝒩^{−1} (Woodbury);
- E^{−1}J⁻ = 𝒩^{−1}J⁻M^{−1}.

With a_h = s_h + |s|/4 the off-diagonal part of M is ½W[1/(a_h + a_g)]W, a Cauchy matrix. Its diagonal part is not of
the Cauchy form, so M^{−1} has no evident closed form.

*Proof.* The last identity: E^{−1}J⁻ = 𝒩^{−1}J⁻[I − ½M^{−1}J⁺𝒩^{−1}J⁻] = 𝒩^{−1}J⁻[I − M^{−1}(M − I)]. ∎

Consequently 2^Z⟨Ψ, YΨ⟩ splits into local terms (diagonal 𝒩^{−1}; E_1^{−1} by Sherman–Morrison) plus the crossed term
+½Σ_{a,b} ⟨β_ab, M^{−1}β_ba⟩ with β_ab = J⁺𝒩^{−1}𝒯_bψ_a ∈ R^k. This term is not definite: it is the symmetric part
minus the antisymmetric part. Bounding the two parts separately loses the cancellation described after Proposition S,
so it does not work.

**Lemma BJ (the J⁻ part).** Let ψ_r = Σ_h X_rh J⁻|h⟩ with X ∈ R^{k×k} (the J⁻ part of G1; dimension k², and it
contains the k kernel vectors X = wc^T). Then the B-term is

    Σ_{a,b} ⟨𝒯_bψ_a, E^{−1}𝒯_aψ_b⟩ = Σ_{h,g} ⟨G_h X_{·g}, E^{−1} G_g X_{·h}⟩,   G_h b := (E − 2s_h) σ_b⁻|h⟩ + 2w_h b_h J⁻|h⟩,

where X_{·g} is column g and σ_b⁻|h⟩ = Σ_{t≠h} b_t|{h,t}⟩. *Proof.* Lemma J gives 𝒯_tJ⁻|h⟩ =: g_t^h, so
𝒯_tψ_r = Σ_h X_rh g_t^h, and the crossed sum regroups. ∎

**Proposition S (soft modes).** Fix h and let a ∈ R^k with a_h = 0 and a ⊥ w. Put Ψ(a,h) := Σ_r a_r e_r ⊗ J⁻|h⟩ (so
X = a e_h^T), and u := σ_a⁻|h⟩. Then Ψ(a,h) ⊥ p ⊗ H_{k−2}, |Ψ(a,h)|² = (|s| − s_h)|a|², and

    F(Ψ(a,h)) = 2s_h|a|² − 4s_h² ⟨u, E^{−1}u⟩ − ½|c|²,   c = 4s_h [ ½⟨E_1^{−1}e_h, a⟩ w + (E_1^{−1}e_h) ∘ a ∘ w ].

*Proof.*
- **The A-part.** X[r, {h,t}] = a_r w_t and X vanishes on pairs without h. The P ∌ r sum gives
  ½Σ_r a_r²(|s| − s_h − s_r). Only Y_h is nonzero, with Y_h[r] = a_r w_r and Σ_r a_r w_r = 0. Total:
  ½(|s| − s_h)|a|² + Σ_r s_r a_r².
- **The B-part.** By Lemma BJ with X = a e_h^T it is ⟨G_h a, E^{−1}G_h a⟩. Since a_h = 0, G_h a = (E − 2s_h)u. So it equals
  ⟨u, Eu⟩ − 4s_h|u|² + 4s_h²⟨u, E^{−1}u⟩. Here |u|² = |a|², J⁺u = w_h a + (a·w)e_h = w_h a, and
  ⟨u, 𝒩u⟩ = Σ_t a_t²(½|s| + s_h + s_t). So ⟨u, Eu⟩ = (½|s| + (3/2)s_h)|a|² + Σ_t s_t a_t².
- **Subtract.** A − B = 2s_h|a|² − 4s_h²⟨u, E^{−1}u⟩.
- **The coupling.** c_g = Σ_r a_r[⟨J⁻h, σ_r⁻g⟩ − ⟨J⁺J⁻h, E_1^{−1}𝒯_r g⟩], with J⁺J⁻ = 2(E_1 − C) on H_{k−1}. The terms
  without E_1^{−1} vanish because a ⊥ w and a_h = 0. Further, Σ_r a_r𝒯_r = ½aw^T + diag(a ∘ w). This gives the formula
  for c. Moreover ⟨c, w⟩ = 4s_h⟨E_1^{−1}e_h, E_1 a⟩ = 4s_h a_h = 0. ∎

Note what happens to the dangerous term 2w_h b_h J⁻|h⟩ of G_h. It is O(w_h), not O(s_h), and in the crossed sum it is
paired with a zero row coefficient (a_h = 0). It cancels between the symmetric and the antisymmetric parts. This is why
the crossed structure must be kept intact.

**Corollary S (sharpness).** For every k ≥ 3 and all weights, the smallest eigenvalue of the letter-1 block of
(R∞)_{k−2} on the orthogonal complement of p ⊗ H_{k−2} ⊕ span(vac ⊗ J⁻|1…1⟩) is at most min_h s_h/(2(|s| − s_h))
(in units of Y).

*Proof.* Take v = (Ψ(a,h), ½c(Ψ(a,h))). It is orthogonal to the kernel space, because c ⊥ w. By Lemma R1,
⟨v, Yv⟩ = F(Ψ(a,h))/4 ≤ s_h|a|²/2, and |v|² ≥ (|s| − s_h)|a|². ∎

So no bound of the Theorem T type that is uniform in the weights can work. A proof has to be exact to first order in
min_h s_h, in the directions Ψ(a, h).

## 3. The exact reduction that is left, with numerics

Split the complement of p ⊗ H_{k−2} in G1 orthogonally as S ⊕ R, where S = span{Ψ(a,h) : h ∈ [k], a_h = 0, a ⊥ w}.
S has dimension k(k − 2): the copies h give independent column patterns X = a e_h^T. R has dimension
k(k² − 4k + 5)/2 (10 for k = 4, 25 for k = 5). By Lemma R1 the goal is equivalent to:

    F_RR ≻ 0   and   Σ := F_SS − F_SR F_RR^{−1} F_RS ≻ 0 on S.

The entries of F_SS are explicit by Proposition S and Lemma BJ. The diagonal blocks are as in Proposition S. The
off-diagonal blocks are 2F(Ψ(a,h), Ψ(b,g)) = (A-cross) − 2⟨G_h b, E^{−1}G_g a⟩ − ⟨c(Ψ(a,h)), c(Ψ(b,g))⟩, and here
b_h ≠ 0, so the J⁻|h⟩ term contributes at order w_h w_g.

Numerics (`sweep.py`, `ex2.py`, `ex3.py`, `ex5.py`; float, no claim rests on them; units 4Y unless stated):

| quantity | value |
|---|---|
| the full block ⪰ 0, kernel dimension exactly C(k,2) + 1 | yes: k = 4…7, 610 weight vectors, min eigenvalue ≥ −5 × 10^{−16} relative, \|Yκ\| ≤ 1.2 × 10^{−16} |
| gap / (min_h s_h/(2(\|s\| − s_h))) (Y units; ≤ 1 by Corollary S) | in [0.369, 1.000] (k = 4), [0.501, 1.000] (k = 5), [0.746, 1.000] (k = 6, 7) over the sweep; 0.344, 0.368, 0.386, 0.400 at equal weights, k = 4…7 (`ex6.py`) |
| gap / x_min (Y units) | ≥ 0.241 (k = 4), 0.291 (k = 5), 0.40 (k = 6, 7) over the sweep; 0.2295, 0.2300, 0.2317, 0.2335 at equal weights, k = 4…7 |
| λ_min(F on R) | ≥ 0.4506 over the whole sweep (k = 4…7); tends to ½ when one or two weights dominate |
| F on S_h, divided by x_h (orthonormal basis of S_h) | eigenvalues in [0.50, 2.00]; → 2 as x_h → 0; ≈ ½–⅔ for a dominant copy |
| one weight ε → 0 (`ex5.py`) | the k − 2 lowest eigenvalues are (s_h/(2\|s\|))·κ_i + o(s_h), with κ_i ∈ [0.90, 0.97] depending on the other weights; the eigenvectors lie in S_h up to 10^{−4} |
| cross terms | ‖F(S_h, S_g)‖ ≤ 0.67 √(x_h x_g); ‖F(S_h, R)‖ ≤ 0.57 √x_h |

The last row explains why naive block dominance fails. The rescaled soft blocks (eigenvalues ≥ ½) are not larger than
the rescaled couplings (up to 0.67 between two large copies, 0.57 to R). Both couplings are genuinely O(w_h): the
Schur correction is of the same order s_h as F_SS. This is the κ_i < 1 in the ε → 0 row. So the Schur complement Σ must
be computed, not estimated by norms. The boundary problem is a clean intermediate target. At the face s_h = 0, the first-order coefficient matrix
K_h = lim_{s_h→0} (Σ restricted to S_h)/x_h (orthonormal basis of S_h, units 4Y) is a (k − 2) × (k − 2) matrix
determined by the (k − 1)-copy subgroup [k] ∖ h. Its eigenvalues are 2κ_i, with κ_i ∈ [0.90, 0.97] in the examples of
`ex5.py`. Proving K_h ≻ 0 for every subgroup and every weight vector on it is a necessary condition for the goal, and
it is open.

## 4. What remains

- **The goal itself (Proposition K2⁺), for every k, including k = 4 at general weights.**
- A route suggested by the structure above:
  1. A Theorem T-type lower bound F ⪰ c·I on R. Numerically c = 0.45 works for k ≤ 7; the R-directions contain the
     whole ker J⁺ part, where every 𝒯 is local (Lemma N).
  2. The first-order (face) problem K_h ≻ 0, which lives on the (k − 1)-subgroup.
  3. An exact treatment of Σ on S, rescaled by 1/√s_h, for instance by a pairwise structure of Σ in the copies
     (h, g), as in the proof of Proposition U1. This was not tested.
- For k = 4 alone: Σ is 8 × 8 and F_RR is 10 × 10, with entries rational in w through det E (degree 6) and det M. This
  is smaller than the 18 × 18 reduced block (whose determinant has degree about 70–80 per variable), but no certificate
  was attempted.
- Dead ends met here:
  - bounding the symmetric and antisymmetric parts of the crossed E^{−1} (or M^{−1}) term separately: it loses an
    O(s_h) amount on Ψ(a, h), the same order as the gap;
  - plain or weighted (√x_h) block Gershgorin on S ⊕ R or on Σ: the ratio reaches 1.0.

## Files (run `bash run_all.sh`; one Python process at a time, each < 0.3 GB)

| File | Content |
|---|---|
| `mid.py` | closed form of the block (Lemma B, j = k − 2), any number type |
| `direct.py` | direct construction from M, V_p, D (via `check_bath_reduction.group_ops`) |
| `check_closed_form.py`, `.log` | closed form = direct construction (float, k = 4, 5); the kernel vectors exactly (k = 4, 5, 6); PASS |
| `check_exact.py`, `.log` | exact: Lemma R1, Lemma J, Lemma W, Lemma BJ, Proposition S (k = 4, 5, 6, two weight vectors each, one with w_0 = 1/40); PASS |
| `sweep.py`, `.log` | float sweep k = 4…7 (PSD, kernel dimension, gap against Corollary S, hard part), direct spot checks; PASS |
| `ex1.py` … `ex6.py` | exploration (J⁻-part spectra, S/R split, Schur complement on S, dominant weights, ε → 0 asymptotics, equal weights); no claim rests on them |
