# Squared partial traces and the nuclear norm

Ansh Mishra and Aryan Senthilkumar

This repository contains a short paper, with verification code, on a new inequality for partial traces and what it
says about multi-copy distillability of Werner states.

## The inequality

For every operator Z on a finite-dimensional tensor product H₁ ⊗ … ⊗ Hₙ (not necessarily Hermitian, of any rank),

    Σᵢ ‖Trᵢ Z‖₂²  ≤  (n − 1) ‖Z‖₂² + ‖Z‖₁²,

where Trᵢ is the partial trace over the i-th factor, ‖·‖₂ is the Hilbert–Schmidt norm and ‖·‖₁ the trace norm. It is
sharp: equality holds for Z = P ⊗ |c⟩⟨c| with P ≥ 0 on one factor and c a product vector.
- For n = 2 the inequality follows from a recent bound of Costa Rico and Wolf (arXiv:2507.18278, Prop. 9).
- For n ≥ 3 it is new as far as we know.

The proof is an induction on n. It reduces everything to one statement about block operators: for
Z = Σ |a⟩⟨b| ⊗ Z_ab,

    Σ_ab ‖Z_ab‖₁² + ‖Σ_a Z_aa‖₂²  ≤  Σ_ab ‖Z_ab‖₂² + ‖Z‖₁².

We prove it with a balanced factorisation Z = PQ*. One step uses the second-compound estimate of Wu and Zou
(arXiv:2608.02647, Lemma 5.1).

## Consequences

- **The first-order Werner bound.** For every operator C of rank at most two,
  Σᵢ ‖Trᵢ C‖₂² ≤ (n + 1) ‖C‖₂², for every n and every local dimension.
  - This is Assumption (i) of Chattopadhyay and Sarkar (quant-ph/0609050, 2006).
  - It is the first-order term of the n-copy distillability problem for Werner states at the endpoint α = 1/2.
  - For rank r the constant is n + r − 1, and this is sharp.
- **All of Chattopadhyay and Sarkar's assumed bounds.** For every m, the sum of ‖Tr_J C‖₂² over |J| = m is at most
  [C(n−1, m) + 2C(n−1, m−1)] ‖C‖₂² for rank-two C. The case m = n − 1, which involves only one-site marginals, follows
  from Costa Rico and Wolf.
- **A Han inequality for the linear entropy.** For every quantum state, Σᵢ S_L(ρ with factor i traced out) ≥
  (n − 1) S_L(ρ), where S_L = 1 − Tr ρ². This holds even though strong subadditivity fails for S_L.
- **Schmidt-number witnesses.** On (Cᵈ ⊗ Cᵈ)^⊗n, the operator (n + r − 1) I − Σᵢ P⁺ on pair i is nonnegative on all
  vectors of Schmidt rank ≤ r. Here P⁺ = Σ_jk |jj⟩⟨kk|. For r < d the operator is negative on some vector of Schmidt
  rank r + 1.
- **A sharp rank-two bound when C and C* share a two-dimensional range.** For every C of rank two with
  ran C + ran C* of dimension two (C = VMV* with V an isometry from C²; in particular every normal C), Σᵢ ‖Trᵢ C‖₂² ≤ n ‖C‖₂² + ½ |Tr C|². For n = 2 this is the known two-copy
  inequality. In our searches the bound holds for every rank-two C, but we have no proof when ran C + ran C* is larger.
- **A sharp two-factor bound with a flatness term.** For C of rank at most r on H₁ ⊗ H₂,
  ‖Tr₁ C‖₂² + ‖Tr₂ C‖₂² ≤ ‖C‖₂² + (1 − 1/r) ‖C‖₁² + (1/r) |Tr C|².
  - This follows from the proof of Fraser, Huber, Pozsgay and Vona (arXiv:2607.24309). For r = 2 the paper gives an
    explicit sum-of-squares identity.
  - It is an equality for ρ ⊗ |c⟩⟨c| with ρ ≥ 0 of rank r and any spectrum.
  - For r = 2 it reads q₂(C) ≥ ¼(σ₁ − σ₂)², which is stronger than the two-copy inequality.
  - The n-factor analogues are stated as conjectures, with the numerical evidence. The analogue of the partial-trace
    form, Σᵢ ‖Trᵢ C‖₂² ≤ (n − 1)‖C‖₂² + ½(‖C‖₁² + |Tr C|²), is proved for normal C.
  - Its rank-r version, with (1 − 1/r)‖C‖₁² + (1/r)|Tr C|², is proved for r = 1 and every n. In that case
    Σᵢ ‖Trᵢ(xy*)‖₂² ≤ (n − 1)‖x‖²‖y‖² + |⟨y,x⟩|², by induction from the two-factor case.
- **A thickness margin for two copies.** For C of rank at most two on H₁ ⊗ H₂,
  q₂(C) ≥ min{τ₁, τ₂} ‖C‖₂² / (48 − 8√2),
  where τ_l is the sum of all but the two largest eigenvalues of the reduction, on factor l, of the projector onto the
  co-range of C.
  - In T. S. Tomar's notation this is m₂(V) ≥ min_l λ₃⁽ˡ⁾(V)/37: the case k = 2 of his Conjecture 7.1, with a smaller
    constant than his 1/4. He noted that no linear bound of this kind was known for two copies.
  - The proof uses a new sharp inequality: λ_max(A) ≤ 1 + ‖A₀₀ − A₁₁‖₂/√2 for A ≥ 0 on C² ⊗ K with Tr_K A = I₂,
    and a refinement, ‖A₀₀ − A₁₁‖₂² ≥ 2(λ_max − 1)(λ_max − 1 + τ) with τ the tail of A₀₀ + A₁₁ beyond its two
    largest eigenvalues. The refinement gives the best possible constant in the pair lemma.
  - If the range and co-range of C lie in one plane V (for instance C normal), the paper proves Tomar's inequality with
    the best possible constant: q₂(C) ≥ (1 − 1/√2) min_l λ₃⁽ˡ⁾(V) ‖C‖₂², which beats his conjectured 1/4. Tomar's
    extremal example is Hermitian and attains it. It also proves q₂(C) ≥ ((3 − √5)/4) max_l τ_l(V) ‖C‖₂² there.
    Tomar proved the constant 1 − 1/√2 on a different class of planes (his Theorem 7.2), with no condition on the
    co-range.
  - The paper now also derives the two-factor sum-of-squares identity that this rests on.

The paper also records a graded form of the conjecture behind the endpoint. The first-order bound is the first of n
inequalities which together give the "no collective gain" conjecture on the whole parameter range. The last of them is
the endpoint itself.

For three copies the middle inequality is an explicit nonnegative combination of three things:
- the three-copy endpoint form;
- the two-copy theorem applied to the three ways of splitting the copies into one and two;
- |Tr C|² ≤ 2‖C‖₂².

So, given T. S. Tomar's computer-assisted three-copy theorem, the "no collective gain" conjecture holds for three
copies on the whole range. For four or more copies a linear program finds no such combination.

## What this does not do

It does not settle whether Werner states with 1/d < α ≤ 1/2 are undistillable from every number of copies. The
endpoint inequality needs lower bounds on the even-order terms, and the paper shows numerically that natural rank-free
versions of it fail. The rank-two constraint has to enter in an essential way there, while at first order it does not.

## Status

- The proof is short and elementary.
- It has been checked independently, step by step, including degenerate cases, with separate code and adversarial
  numerical searches.
- The paper has not been refereed.

## Contents

- `paper/`: LaTeX source and PDF.
- `code/check_proof_R.py`: checks every step of the proof of the block lemma, and the lemma itself, on 3000 random
  block operators, including near-equality ones; also checks the main inequality on random operators.
- `code/check_R_degenerate.py`: the block lemma on 4000 larger or degenerate cases (zero blocks, rank one, real,
  nilpotent, block permutations, indefinite Hermitian).
- `code/check_corC.py`: the bound for partial traces over m sites, for all m.
- `code/check_witness.py`: the Schmidt-number witnesses, including the rank-(r + 1) vectors that violate them.
- `code/rung_lp.py`: the linear program for the graded levels (derives the middle level at n = 3; finds no combination
  for n = 4, 5). Run as `python rung_lp.py n`.
- `code/check_remarks.py`: the statements of the Remarks section:
  - the graded identity, and the three-copy identity behind the middle inequality;
  - the sharp rank-two bound for normal operators and for operators C = VMV* whose range equals that of C*;
  - the sharp two-factor bound for r = 2, 3, and its r = 2 sum-of-squares identity;
  - the n-factor bound with a flatness term for normal C, and its rank-one case for every n.
- `code/check_thickness.py`: the two lemmas and the proposition behind the two-copy thickness margin (the explicit
  certificate, the inequalities on random instances, and the equality family), and the bound for normal operators.

Everything needs only Python 3 with NumPy; `rung_lp.py` also needs SciPy. Typical run, from `code/`:

    python check_proof_R.py
    python check_R_degenerate.py
    python check_corC.py
    python check_witness.py
    python check_remarks.py
    python check_thickness.py
    python rung_lp.py 4

## License

MIT.
