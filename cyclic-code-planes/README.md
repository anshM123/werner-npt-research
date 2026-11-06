# Cyclic qutrit code spaces

Ansh Mishra and Aryan Senthilkumar, 6 November 2026. Nothing here is refereed.

For every number of copies n ≥ 4 there is an explicit three-dimensional subspace Q_n of (C³)^{⊗n} on which the endpoint
form is positive definite. It holds for operators supported on Q_n on both sides at every n ≥ 4, and for arbitrary
co-range at n = 4, …, 13, n = 18, …, 26 and every n ≥ 30. Every vector of Q_n has maximally mixed one-copy marginals, so
these spaces are thick at every copy. The positivity is uniform with an explicit margin, so it extends to an explicit
neighbourhood of every plane inside Q_n in the whole Grassmannian. The statements hold for every rank, not only rank
two. They do not settle the endpoint for general planes.

## The code

Let ω = e^{2πi/3}, and for n ≥ 4 let

    |g_n⟩ = 3^{−n/2} Σ_{x ∈ F₃ⁿ} ω^{Σ_i x_i x_{i+1}} |x⟩     (indices mod n),

the graph state of the n-cycle. Let Q_n = span{Z^{t}⊗…⊗Z^{t} |g_n⟩ : t ∈ F₃} (Z|x⟩ = ω^x|x⟩) and P_n its projector. The
three vectors are orthonormal. The stabiliser of Q_n has label space S = {(u, Gu) : Σ_i u_i = 0}, with G the adjacency
matrix of the cycle.

**Balanced marginals.** A Pauli label (x, z) normalises S exactly when z − Gx is a constant vector b·(1, …, 1). A
nonidentity Pauli on one vertex cannot do this: a vertex away from it and its neighbours forces b = 0, then x = 0, then
z = 0. So P_n E P_n = 0 for every nonidentity one-copy Pauli E, and every vector of Q_n has one-copy marginals I₃/3.

## The form as a transfer trace

With Φ(X) = X − ½ Tr(X) I on one qutrit, q_n(C) = ⟨C, Φ^{⊗n}(C)⟩_HS; Φ is 1 on nonidentity Paulis and −½ on I, so
‖Φ^{⊗n}‖ = 1. The operators P_n E, one Pauli E from each coset of S, are an orthogonal basis of {C : P_n C = C}, each of
squared norm 3, and different cosets have disjoint Pauli expansions. So the compressed form is diagonal in this basis.
A coset is labelled by a = Σ_i x_i and d = z − Gx ∈ F₃ⁿ, and its eigenvalue is

    λ_{a,d} = Z_{a,d} / (3^{n−1} 2ⁿ),   Z_{a,d} = Σ_{u : Σu_i = a} Π_i f(u_i, u_{i−1} + u_{i+1} + d_i),

with f(y, v) = −1 if y = v = 0 and f(y, v) = 2 otherwise. Over R = Z[t]/(t³ − 1) let T_b (b ∈ F₃) be the 9 × 9 matrix
on states (x, y) ∈ F₃² with entries (T_b)_{(x,y),(y,z)} = f(y, x + z + b) t^y. Expanding closed walks gives

    Z_{a,d} = [t^a] Tr(T_{d₁} ⋯ T_{d_n}).

Matrices over R with nonnegative coefficients are closed under products; if every coefficient of A is ≥ m > 0 and every
coefficient of B is ≥ m′ > 0, every coefficient of AB is ≥ 27mm′.

## Results

**Theorem A (both sides in Q_n, every n ≥ 4).** For every nonzero C = P_n C P_n, q_n(C) > 0, with no rank restriction.

*Proof.* Only normaliser cosets occur, so d = b(1, …, 1) is constant. Every coefficient of T_b^m is positive for
b = 0, 1, 2 and m = 8, …, 15 (least coefficient 432). Every n ≥ 8 is 8q + r with 8 ≤ r ≤ 15, so T_b^n is coefficientwise
positive; n = 4, …, 7 are checked directly. ∎ The least eigenvalues at n = 4, …, 7 are 144/432, 639/2592, 3888/15552 and
15903/93312.

**Theorem B (range in Q_n, arbitrary co-range).** For every nonzero C with P_n C = C, q_n(C) > 0, at n = 4, …, 13,
n = 18, …, 26 and every n ≥ 30. Equivalently, the compression of the partially transposed Werner operator to
Q_n ⊗ (C³)^{⊗n} is positive definite.

*Proof.* For all 3¹² words of length 12 and all 3¹³ of length 13, every coefficient of T_{d₁}⋯T_{d_ℓ} is positive (least
coefficients 365 472 and 2 773 116). For lengths 18, …, 23 write the product as PQ with P of length 12 and Q of length
r = 6, …, 11. Let m_i and M_i be the least and largest coefficients in row i of P. Over all P, m_i ≥ 365 472 and
M_i/m_i ≤ 4 764 744/365 472. Then every coefficient of PQ is at least m_i S⁺_j − M_i S⁻_j ≥ 365 472 S⁺_j − 4 764 744 S⁻_j,
where S^±_j are the sums of the positive and of the absolute negative coefficients in column j of Q. An exhaustive check
over all Q of each length r = 6, …, 11 shows this bound is strictly positive. Every n ≥ 30 is a sum of lengths from
{12, 13, 18, …, 23}; for example 30 = 18 + 12, then add 12s. The trace is cyclic, so every coset trace is positive. The
sizes n = 4, …, 13 are checked by exhaustive enumeration of all cosets, and 18, …, 26 are lengths or sums of the positive
blocks. ∎ The sizes 14, …, 17 and 27, …, 29 are not covered by these certificates.

**Explicit margins.** For a decomposition n = ℓ₁ + ⋯ + ℓ_m into the block lengths, with g_ℓ the certified least
coefficient of length ℓ, q_n(C) ≥ μ_n ‖C‖_F² with μ_n = 9 · 27^{m−1} Π g_{ℓ_i} / (3^{n−1} 2ⁿ).

**Corollary (a neighbourhood of every code plane).** Let L₀ ⊂ Q_n be a plane at a size covered by Theorem B, with
orthonormal frame W₀. If W is the orthonormal frame of a plane L with ‖W − W₀‖ ≤ μ_n/4, then q_n(C) ≥ (μ_n/2)‖C‖_F² for
every C with range in L and any co-range. *Proof.* Write C = WX, C₀ = W₀X; since ‖Φ^{⊗n}‖ = 1,
|q_n(C) − q_n(C₀)| ≤ ‖C − C₀‖(‖C‖ + ‖C₀‖) ≤ 2‖W − W₀‖‖X‖². ∎

## Checks

- `code/transfer_certificates.js` (`logs/transfer_certificates.log`) checks:
  - the homogeneous powers T_b^m for m = 8, …, 15;
  - every word of length 12 and 13 exhaustively (531 441 and 1 594 323 words, no product with a nonpositive
    coefficient);
  - the coset minima for n = 4, …, 13;
  - the transfer formula against a direct construction of the code projector and of all partial traces, at 31 random
    cosets on four and five copies.
- `code/tail_certificates.js` (`logs/tail_certificates.log`) checks the row ratio of the length-12 products, the six
  column-sum bounds for r = 6, …, 11, and the length decompositions of n = 30, …, 41.
- `code/cyclic_check.py` (`logs/cyclic_check.log`), written separately from the transfer formula, builds the code
  vectors and the operator M^{⊗n} from their definitions and computes the least eigenvalue of the compressed form. At
  n = 4 and n = 5, one-sided and two-sided, it gives 1/3 = 144/432 and 639/2592, the values of the transfer formula
  (floating point, 12 digits).

All certificates use exact integer arithmetic. Run `bash code/run_all.sh`; it needs Node (no packages) and Python with
NumPy and SciPy.

## Relation to other folders

`compression-notes/` has an exact formula on the code spaces of pure qutrit codes. It gives positivity on the
five-qutrit code space and on code spaces whose effective parameter is small. The cyclic code here is a family over all
n, with certificates uniform in n.
