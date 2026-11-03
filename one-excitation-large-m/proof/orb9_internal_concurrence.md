# An internal two-qubit concurrence bound for actual star

This is a rigorous finite-dimensional reformulation for arbitrary
rank-two-Y data, including simultaneous outside-plane fields. It is
not a proof of universal star. The concurrence estimate needed after
division by the actual positive denominator remains open. No numerical
search is used here.

## Actual factors and replica convention

Put H=(C^2)^(tensor m), and identify its two replicas with the row and
column sector spaces. Write

    Y=UMV*,    A=M^(-1),    X=R A L*,

where U,V:C^2->C^3 are isometries, M is invertible, and

    R_a=c_a V+E_a,      L_b=d_b U+F_b.

Here R,L:C^2->C^3 tensor H stack the sector coordinates. Choose the
scalar components by free Frobenius projection, so
Tr(V*E_a)=Tr(U*F_b)=0. The residual costs are exactly

    p=qtilde(P)=sum_a ||E_a||_F^2,
    q=qtilde(Q)=sum_b ||F_b||_F^2.

Nothing restricts the residuals to their two-dimensional free planes.
The letters R,L are full rectangular factors; they are not density
operators. All internal vectorizations use the order (00,01,10,11).

Let F_0 swap the two free qutrit replicas, let F_s swap sector qubit s,
and put F_A=product_s F_s. On physical replicas define

    W_0=I-F_0/2,
    D_m=product_s(2I-F_s)-2I+F_A.

Both are positive semidefinite. For D_m, the common swap eigenspace
with k negative swaps has eigenvalue 3^k-2+(-1)^k, which is zero for
k=0,1 and positive for k>=2.

For a rank-one physical operator xy*, the exact replica identities in
Lemma R of the two-sector notes (a working note not included here) give

    Sigma_full(xy*)-Gamma(xy*)
       =<x tensor y|W_0 tensor D_m|x tensor y>.

The canonical permutation regrouping free and sector replicas is
understood. Polarizing the product projectors yields the following
four-by-four internal matrix:

    B=(R* tensor L*)(W_0 tensor D_m)(R tensor L)>=0,

    H(A):=Sigma_full(X(A))-Gamma(X(A))
        =<vec A|B^(Gamma_2)|vec A>
        =Tr[B (|vec A><vec A|)^(Gamma_2)].             (1)

Gamma_2 transposes the second internal qubit. In particular there is
no conjugate on L in the formula for B. The usual conjugate on the
right factor of vec(R A L*) disappears under this partial transpose.

Equivalently, form the actual four-by-four Hermitian Gram matrices
G,S of Gamma,Sigma on X_ij=R_i L_j*. Then

    S-G=B^(Gamma_2).                                 (2)

Both G and S are positive semidefinite: every X(A) has rank at most
two, so the existing two-copy theorem for Gamma and thin-copy theorem
for Sigma apply. Matrix B is independent of A, but depends on all
the actual factor fields and their sector correlations.

## The pure determinant inequality and concurrence

For any two-by-two coefficient matrices a,b,

    <vec a|(|vec b><vec b|)^(Gamma_2)|vec a>
          >=-2|det a| |det b|.                       (3)

To prove this, use local Schmidt unitaries to put b=diag(s,t), s,t>=0.
They preserve both determinant magnitudes. The left side is

    s^2|a_00|^2+t^2|a_11|^2
       +2st Re(conjugate(a_01)a_10)
    >=2st(|a_00 a_11|-|a_01 a_10|)
    >=-2st|a_00 a_11-a_01 a_10|.

There is no division, so zero Schmidt values and singular a are
included. For real a the exact square identity is

    left side+2st det(a)=(s a_00+t a_11)^2.

Define unnormalized two-qubit concurrence by

    C(B)=inf_(B=sum_j |vec b_j><vec b_j|)
                    sum_j 2|det b_j|.

Its closed form is max(0,lambda_1-lambda_2-lambda_3-lambda_4), where
the lambda_i are the decreasing nonnegative square roots of the
eigenvalues of B (sigma_y tensor sigma_y) conjugate(B)
(sigma_y tensor sigma_y). This is the homogeneous form of
[Wootters' two-qubit concurrence formula](https://arxiv.org/abs/quant-ph/9709029).

Apply (3) to every ensemble and then take the infimum. Equation (1)
gives

    H(A)>=-|det A| C(B).                              (4)

Thus whenever Gamma(A)+Sigma(A)>0,

    sqrt(Gamma(A))-sqrt(Sigma(A))
       <= |det A| C(B)/(sqrt(Gamma(A))+sqrt(Sigma(A))). (5)

If both forms vanish the star gap is zero and no division is needed.
The positive axis contributions have not been dropped in (4)-(5).
This is why (4) does not assert a determinant factor for the isolated
mixed polarization refuted in
the actual mixed-term obstruction (a working note not included here).

## The axes are separable

If E=0, the range of X(A) lies in C^3 tensor c. Lemma R' therefore
gives H(A)>=0 for every A. Hence B and B^(Gamma_2) are positive, so
the internal two-qubit B is separable. The same holds if F=0.

A constructive version also explains the dimension dependence.
Contract D_m on the left sector vector c to obtain a positive
right-sector operator D_c. Compress the left free replica by V. The
resulting two-by-three operator (V* tensor I)W_0(V tensor I) has
partial transpose I_6-|vec V><vec V|/2>=0 and is separable. Its
product decomposition, tensor D_c, followed by the local right map
Z->L* Z L, proves that B is separable. The adjoint argument covers
F=0.

This does not give a product decomposition of uncompressed W_0 on
three-by-three space. In fact

    W_0^(Gamma)=I_9-|vec I_3><vec I_3|/2

has eigenvalue -1/2 on vec I_3. Any general both-outside argument
which assumes a product frame of W_0 would be invalid.

## A four-by-four sufficient matrix condition

Put K=V*U, and q_K(Z)=||Z||_F^2-|Tr(K Z)|^2/2. In the standard
principal-plane convention K=diag(1,k). Let J_K be the
positive four-by-four matrix satisfying

    q_K(adj A)=<vec A|J_K|vec A>.

The adjugate is linear for two-by-two matrices. Since M=A^(-1),

    q_1(Y)|det A|^2=q_K(adj A).                       (6)

For fixed actual R,L, a sufficient condition for star for every
invertible A is the finite Loewner inequality

    C(B)^2 J_K <= (2^m-2)^2 p q (G+S).               (7)

Indeed (6)-(7), together with
(sqrt(Gamma)+sqrt(Sigma))^2>=Gamma+Sigma, imply star through (5).
This is only a sufficient condition. No claim is made that all
actual factors satisfy (7), or that (7) is necessary for star.
It is a four-by-four target retaining the actual sector correlations,
rather than an unrestricted four-by-four model substituting for them.

There is also a less restrictive exact one-parameter formulation of
this concurrence-based sufficient condition. Put

    T=C(B)^2 J_K/[(2^m-2)^2 p q],

when p q>0. The uniform inequality

    <a|T|a><=(sqrt(<a|G|a>)+sqrt(<a|S|a>))^2

for all internal a is equivalent to

    z^2 G+z(G+S-T)+S>=0 for every real z>0.           (8)

Indeed for nonnegative g,s,
(sqrt(g)+sqrt(s))^2=inf_(z>0)[(1+z)g+(1+1/z)s].
The infimum includes the limits z->0 or infinity when one form
vanishes. Quantifying over all a and all z can be interchanged,
which gives (8) after multiplication by z. This is a degree-two
matrix polynomial of order four, with all sector correlations
already contained in its coefficients. It remains only a
sufficient route to star because (4) itself can be strict.

### There is no exact common-kernel obstruction when concurrence is positive

If C(B)>0, then G+S is positive definite. This conclusion is stronger
than merely observing it in the example, but does not give a uniform
lower eigenvalue as the factors vary.

For completeness, a non-positive-semidefinite block-positive Hermitian
operator W on C^2 tensor C^2 is nonsingular. Suppose instead that
0!=vec K lies in ker W, and choose vec V with strictly negative
W-expectation. If K has rank two, det(V+zK) has nonzero quadratic
coefficient. If K has rank one, choose V within the open set of
negative vectors so that its linear coefficient
Tr(adj(V)K) is nonzero; that linear polynomial in the entries of V
is not identically zero for nonzero K. In either case det(V+zK)
is nonconstant and has a complex root z_0. Since W vec K=0, the
vector vec(V+z_0K) has the same negative expectation. It is nonzero
and has rank at most one, contradicting block positivity. This is
the d=2 case of
the maximal-block-positivity nonsingularity lemma (a working note not included here).

Now W=B^(Gamma_2) is block positive because B>=0. Positive
concurrence makes B entangled, so the two-qubit PPT criterion makes
W non-positive-semidefinite. It is therefore nonsingular. But

    ker(G+S)=ker G intersection ker S
            is contained in ker(S-G)=ker W.

Hence ker(G+S)=0. In particular one cannot refute (7) by finding an
actual factor field with C(B)>0, a common null vector of G and S,
and positive q_K(adj A) on that vector. Any quantitative failure
of (7), if one exists, would have to occur with G+S positive
definite. Near-separable B can still make its coercivity arbitrarily
small, so the observation does not establish a uniform estimate.

## Exact actual example and scalar-baseline dependence

Use three sector qubits and set

    c=|000>,    w=|110>+|101>,    U=V=[|0>,|1>],
    R=L=[h|0,c>+|2,w>, h|1,c>],    h>0.

The costs are p=q=2 for every h. The actual 27-Pauli evaluation gives

    B(h)=[ 4h^2+4    0       0       0 ]
         [   0      4h^2   -2h^2    0 ]
         [   0     -2h^2    4h^2    0 ]
         [   0       0       0       0 ].

It is positive. Its partial transpose has a two-by-two principal
minor with determinant -4h^4, hence it is entangled for every h>0.
The eigenvalues of B conjugate-spin-flip(B) are
36h^4,4h^4,0,0, so

    C(B(h))=4h^2.

Consequently no bound on C(B) using only p and q can be uniform in
the scalar baselines. This does not obstruct (5), which retains
the denominator. It also does not refute star.

At h=1 the exact actual matrices are

    G=[25/2  0  0   1/2]     S=[41/2  0  0  -3/2]
      [  0   5  0    0 ]       [  0   9  0    0 ]
      [  0   0  5    0 ]       [  0   0  9    0 ]
      [ 1/2  0  0   1/2],      [-3/2  0  0   1/2].

Removing the left residual gives B=diag(4,0,4,0), and removing the
right residual gives B=diag(4,4,0,0), both separable as predicted.
For h=1, (7) is positive semidefinite at m=3, p=q=2 and C(B)=4;
this is an exact check on these fixed factor fields, not an
all-field proof.

The checker `orb9_internal_concurrence.py` verifies all these actual
Pauli matrices, positivity by every principal minor, both axis
partial transposes, a genuinely complex internal vectorization test,
the scaled-baseline spin-flip characteristic polynomial, and the
stated example of (7). The remaining task is a uniform bound after
division by the denominator in (5), or another use of the actual
structure of B. No such bound is proved here.
