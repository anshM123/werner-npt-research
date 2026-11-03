# A sharp pure-state bound for the free-weighted merged kernel

This extends the pure-state estimate behind
[induction16](induction16_qubit_merged_gap.md) to arbitrary entanglement
between the free copy and both sector factors. It proves a sharp
rank-two curvature bound and improves constants in the rank-one-Y
bridge argument. It does **not** prove the conjectured universal
square-root gap or MB_4.

Let the physical Hilbert space be C^3_free tensor C^2_new tensor H.
On two replicas put

    W=I-F_free/2,
    K=2W(I-F_new)(I-F_H).

For an operator X let

    H(X)=Gamma_bar(X)-Gamma(X),

using the actual merged forms from
[star19](star19_energy_aware_recursion.md). The replica operator K is
positive, since its commuting factors are positive. The induced
quadratic form H is nonnegative on rank-one operators; it can be
negative on rank-two operators because the right singular indices
are crossed.

## 1. Sharp bound on every pure diagonal contraction

For all physical vectors u,v, with arbitrary free/sector entanglement,

    0<=H(uv*)=<u tensor v|K|u tensor v>
       <=2||u||^2||v||^2.                               (1)

The coefficient two is independent of dim H.

To prove it, first work on C^3_free tensor H and set

    A=W(I-F_H)>=0.

For vectors a,b in that space, define

    alpha=Tr(rho_(a,free) rho_(b,free)),
    beta =Tr(rho_(a,H) rho_(b,H)),
    theta=|<a,b>|^2,

where reduced operators are not normalized when the vectors are not
unit. The swap identity gives

    <a tensor b|A|a tensor b>
       =||a||^2||b||^2-beta-alpha/2+theta/2.             (2)

Both alpha and beta are squared Frobenius norms of partial traces of
the rank-one operator ab*. Both of those partial traces have rank at
most three: one is a 3-by-3 matrix and the other is a sum of three
rank-one matrices. Their trace is <b,a>. Thus Schatten Cauchy gives

    theta<=3alpha, theta<=3beta,
    beta+alpha/2>=theta/2.

Consequently

    ||sqrt(A)(a tensor b)||<=||a|| ||b||.                (3)

The dimension-three free factor is used at exactly this trace bound.

Now split the new qubit, reordering tensor factors for notation:

    u=|0> tensor x+|1> tensor y,
    v=|0> tensor z+|1> tensor w,

where x,y,z,w lie in C^3_free tensor H. Directly applying I-F_new gives

    H(uv*)=2||sqrt(A)(x tensor w-y tensor z)||^2.         (4)

The ordinary triangle inequality, (3), and scalar Cauchy imply

    H(uv*)
      <=2(||x|| ||w||+||y|| ||z||)^2
      <=2(||x||^2+||y||^2)(||z||^2+||w||^2),

which is (1). This does not discard the free swap term: inequality
(2) retains it through both reduced overlaps.

The bound is sharp with genuinely free-entangled vectors. For dim H=3
let omega=(|00>+|11>+|22>)/sqrt(3) on free tensor H, and take

    u=|0>_new tensor omega,
    v=|1>_new tensor omega.

Here alpha=beta=1/3 and theta=1 in (2), so the expectation in (2)
is one and H(uv*)=2.

## 2. Crossed polarizations and rank-two curvature

Let B_H be the sesquilinear polarization of the quadratic H. For
rank-one X=ab* and Y=cd*, the exact crossed-index identity is

    B_H(X,Y)=<a tensor d|K|c tensor b>.

Cauchy for the positive physical K, followed by (1), gives

    |B_H(X,Y)|<=2||X||_F||Y||_F.                        (5)

For a rank-two SVD X=sigma_1 u_1v_1*+sigma_2 u_2v_2*, the two
diagonal terms H(u_iv_i*) are nonnegative. Applying (5) only to the
off-diagonal terms proves the dimension-independent bound

    Gamma_bar(X)-Gamma(X)>=-4sigma_1sigma_2.             (6)

In the other direction, (5) yields

    |H(X)|<=2||X||_*^2<=4||X||_F^2                     (rank X<=2).

The first inequality also holds in arbitrary rank. More generally the
nonnegative diagonal terms give

    H(X)>=-2(||X||_*^2-||X||_F^2)

in every rank, by an SVD and (5).

Equation (6) settles one branch of the conjectured auxiliary gap. If

    Gamma(X)>=(sigma_1+sigma_2)^2,

then Gamma_bar>=Gamma-4sigma_1sigma_2, and monotonicity of
sqrt(t)-sqrt(t-4sigma_1sigma_2) proves

    sqrt(Gamma)-sqrt(Gamma_bar)<=2sigma_2.                (7)

The remaining branch Gamma<(sigma_1+sigma_2)^2 is not proved here.
The sector-only flat margin in induction16 does not directly apply to
these free-weighted forms. In particular Gamma_bar can be less than
(sigma_1-sigma_2)^2: a fully product rank-one projector has
Gamma_bar=1/2. Thus (6) alone cannot justify repeating that proof's
second branch.

## 3. A sharp free-entangled rank-two family

With omega as above, put u_0=|0>_new tensor omega and
u_1=|1>_new tensor omega. They are orthonormal. For sigma_1>=sigma_2>=0
take

    X=sigma_1 u_0u_0*+sigma_2 u_1u_1*.

The actual partial traces give exactly

    Gamma=(11/6)(sigma_1^2+sigma_2^2)+(1/3)sigma_1sigma_2,
    Gamma_bar=(11/6)(sigma_1-sigma_2)^2,
    H=-4sigma_1sigma_2.                                 (8)

Thus (6) is sharp for every pair of singular values, even with
nonproduct singular vectors. This family lies in the resolved branch:
Gamma-(sigma_1+sigma_2)^2=(5/6)(sigma_1-sigma_2)^2>=0.
At sigma_1=sigma_2=1 it has Gamma=4, Gamma_bar=0 and gap two, so the
conjectured global auxiliary coefficient two cannot be improved.

## 4. Uniform energy upper bounds used with this lemma

For rank X<=2 the existing two-copy trace-rank argument gives

    0<=Gamma(X)<=2||X||_F^2.

For completeness, writing A0=Tr_free X and B0=Tr_sector X, their
ranks are at most six and three. Hence
|Tr X|^2<=6||A0||^2 and <=3||B0||^2, which proves the upper bound.

The qubit product-frame identity of star19 additionally gives

    Gamma_bar(X)=sum_j w_j Gamma_prev(X_(u_j,v_j))
      <=2sum_j w_j ||X_(u_j,v_j)||^2
      =4||X||^2-2||Tr_new X||^2<=4||X||^2.              (9)

Every contracted X_(u_j,v_j) still has rank at most two. These are
actual forms; positivity uses only the two-copy theorem. Combining
(5), (9), and the fixed-axis decomposition in the independent
[rank-one-Y argument](induction17_rank_one_Y_bridge.md) improves its
finite bridge constant to 50. That downstream argument
must also control the scalar-axis interactions; (1) or (6) alone is not
a proof of the full bridge.

`star20_weighted_qubit_kernel.py` checks the physical new-qubit
antisymmetry identity, the complex crossed polarization, the sharp
free-entangled examples, and the displayed actual partial traces with
exact arithmetic. No optimization or floating-point evidence enters
the proofs.

## 5. A finite-amplitude concurrence bound at a null scalar baseline

Here is a separate consequence useful for the remaining full-Y problem.
It bounds the internal concurrence at every finite amplitude on a stated
scalar-baseline stratum; it is not a finite MB_B estimate.

For a purification matrix Z with four internal columns, write
G(Z)=Z*Z for its two-qubit Gram state and C(G(Z)) for homogeneous
two-qubit concurrence. If G(Z_1) and G(Z_2) are separable, then

    C(G(Z_1+Z_2))
      <=min(||Z_1||_F^2,||Z_2||_F^2)+2||Z_1||_F||Z_2||_F
      <=3||Z_1||_F||Z_2||_F.                           (10)

Indeed pad the common environment and apply an HJW isometry making
the rows z_(1,j) a product-state ensemble of G(Z_1). Apply the same
isometry to Z_2. For a pure two-qubit coefficient matrix, concurrence
is 2|det z|. Determinant polarization and scalar Cauchy give

    C(z_(1,j)+z_(2,j))
       <=C(z_(2,j))+2||z_(1,j)|| ||z_(2,j)||,

because det z_(1,j)=0. Also C(z_(2,j))<=||z_(2,j)||^2.
Summing and using the convex roof gives the first bound with
||Z_2||^2 in place of the minimum; exchange the roles to obtain
the minimum. The final inequality uses min(a,b)<=sqrt(ab).

For arbitrary Z,W, the same argument, using an optimal (or arbitrarily
near-optimal) concurrence ensemble of G(Z), gives

    C(G(Z+W))<=C(G(Z))+2||Z||_F||W||_F+||W||_F^2.        (11)

The only concurrence facts used here are its pure-state determinant
formula, convex-roof definition, and the common-environment HJW
property; see the primary derivation of
[two-qubit concurrence](https://arxiv.org/abs/quant-ph/9709029).

Apply this to arbitrary two-column actual fields

    R_a=c_a V+E_a, L_b=d_b U+F_b,
    U*U=V*V=I_2,
    Tr(V*E_a)=Tr(U*F_b)=0,
    p=||E||_F^2, q=||F||_F^2,
    C_0=||c||, D_0=||d||.

The free planes U,V can differ. E,F can contain arbitrary simultaneous
in-plane and outside components. Suppose only that the sector scalar
baseline is null:

    <c tensor d|(I-F_new)(I-F_H)|c tensor d>=0.           (12)

Let B be the actual internal positive two-qubit matrix

    B=(R* tensor L*)K(R tensor L).

With Z=K^(1/2)(R tensor L), split its purification into scalar-scalar,
E-scalar, scalar-F and E-F components. Positivity of K and (12) kill
the scalar-scalar component exactly. Denote the other three by
Z_E,Z_F,Z_EF. The first two Gram states are separable by the actual
scalar-axis comparison: a fixed sector vector on one range makes
B and its partial transpose positive, and the two-qubit PPT criterion
applies. Equation (1), summed over internal columns, gives

    ||Z_E||_F^2<=4D_0^2 p,
    ||Z_F||_F^2<=4C_0^2 q,
    ||Z_EF||_F^2<=2pq.                                  (13)

Use (10) for Z_E+Z_F and then (11) to add Z_EF. The result is the
uniform finite-amplitude bound

    C(B)<=12 C_0 D_0 sqrt(pq)
            +4sqrt(2)(D_0 sqrt(p)+C_0 sqrt(q))sqrt(pq)
            +2pq.                                      (14)

No common-product assumption on c,d, no coincidence of free planes,
and no condition number bound is imposed. Condition (12) is essential
to this proof. For an arbitrary Schur matrix A, local-filter covariance
of pure determinants and the convex roof give

    C(B_(RA,L))=|det A| C(B_(R,L)).

Thus (14) gives a determinant-compensated bound on the negative part
of Gamma_bar(RAL*)-Gamma(RAL*) on this null-baseline stratum.
Converting (14) alone to a uniform square-root bridge requires a lower
energy estimate through the mixed scalar/residual degenerations.
The subsequent [rectangle theorem](star21_concurrence_rectangle.md)
handles positive scalar-baseline curvature, and the
[uniform bridge](induction24_free_completion_uniform_bridge.md) supplies
the energy conversion with coefficient 2129/10. Formula (14) remains a sharper
concurrence bound on its null-baseline stratum; it is not by itself a
proof of the coefficient-four bridge. The unresolved universal auxiliary
gap branch identified above is also unchanged.
Actual star at m=3 through 16, ORB, and unrestricted endpoint induction
remain open.
