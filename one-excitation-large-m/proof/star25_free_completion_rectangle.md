# A uniform free rectangle from coherent product completion

For the actual free qutrit Werner operator, the universal free rectangle
has coefficient 13/2. Through the exact channel reduction, this proves
coefficient 13 for the full first-three-component physical concurrence
rectangle, and coefficient 21 for the scalar-baseline term of the full
physical rectangle. All free planes, complex scalar amplitudes, and
simultaneous outside residuals are included.

The proof uses the sharp coherent-addition theorem in
[orb24](orb24_coherent_separable_sum.md). It does not prove the proposed
sharp free coefficient two, MB_4, ORB, or the unrestricted all-copy
endpoint. Any merger or actual-star consequence needs the separate
energy argument. That argument is now supplied in
[induction24](induction24_free_completion_uniform_bridge.md): it proves
the uniform merger coefficient 2129/10 and actual star for every m>=17,
hence the specified one-free-copy endpoint class for n>=18. General
star at m=3 through 16 remains open, as do ORB and unrestricted EP.

## 1. Product compressions are separable

Put W=I-F_free/2 on C3 tensor C3. For arbitrary 3-by-2 matrices R,L,
define

    B(R,L)=(R tensor L)*W(R tensor L).

This two-qubit positive matrix is separable, including singular R,L.
Indeed, for row-major vectorization and every complex 2-by-2 T,

    <vec T|B(R,L)^Gamma|vec T>
        =||RTL*||^2-|Tr(RTL*)|^2/2 >=0.                (1)

The last inequality holds because rank(RTL*)<=2 and
|Tr A|<=||A||_*<=sqrt(2)||A||_F on rank-at-most-two matrices.
Thus B(R,L) is PPT, hence separable in two-qubit dimensions. Its trace
also satisfies the exact identity

    Tr B(R,L)=||R||^2||L||^2-||R*L||^2/2
              <=||R||^2||L||^2.                      (2)

The Frobenius norms in (1)-(2) do not require isometric R,L. These
facts apply to every product compression used below.

## 2. The free coefficient 13/2

Let V,U:C2 to C3 be arbitrary isometries, and let E,F be arbitrary
complex 3-by-2 matrices with Tr(V*E)=Tr(U*F)=0. For any complex alpha put

    Z=alpha V tensor U+E tensor U+V tensor F,
    e=||E||, f=||F||, h=max(e,f), a=|alpha|.

Then homogeneous concurrence obeys

    C(Z*WZ)<=(13/2)ef.                                (3)

If e=0 or f=0, Z is itself a product operator and (1) gives concurrence
zero. Hence assume e,f>0. Use the same positive square root of W for
all purifications in each coherent sum. The theorem of orb24 says

    C(G(X+Y))<=2||X||||Y||                            (4)

whenever both individual Grams are separable, with arbitrary relative
purification environment.

First suppose a<=9h/10. If h=e, split

    Z=(alpha V+E) tensor U+V tensor F.

Both terms give separable Grams by (1). Equations (2),(4) and the
trace orthogonality give

    C(Z*WZ)<=2sqrt[(4a^2+2e^2)(2f^2)]
             =4ef sqrt(1+2a^2/e^2).

If h=f, use instead E tensor U+V tensor (alpha U+F), with the
symmetric estimate. In either case

    C(Z*WZ)<=4ef sqrt(1+2a^2/h^2)
             <=(2/5)sqrt(262)ef <(13/2)ef.            (5)

The last comparison follows from 4192<4225.

Now suppose a>=9h/10. Here alpha is nonzero. Complete the product:

    Z_full=(alpha V+E) tensor (U+F/alpha),
    Z=Z_full-(E tensor F)/alpha.                      (6)

Both terms in the second line again have separable Grams by (1).
Their relative minus sign and the phase of alpha are allowed in (4).
Using (2) and the two trace-orthogonality relations,

    C(Z*WZ)
      <=2ef sqrt[(2+e^2/a^2)(2+f^2/a^2)]
      <=2ef(2+h^2/a^2)
      <=(524/81)ef <(13/2)ef.                         (7)

The final comparison is 1048<1053. This proves (3) for every alpha,
including alpha=0 by the first case. There is no condition-number or
positive-baseline assumption.

## 3. The actual full-field physical rectangle

Use K_H=2W_free tensor D, where
D=(I-F_new)(I-F_H), and the actual two-column fields

    R=V tensor c+E, L=U tensor d+F,
    Tr(V*E_z)=Tr(U*F_z)=0 for every sector coordinate z.

Write C0=||c||, D0=||d||, p=||E||^2, q=||F||^2. Denote the four
purification components from their products by Z0,Z1,Z2,Z3, where
Z0 is scalar-scalar and Z3 is residual-residual. Put

    B012=G(Z0+Z1+Z2), B=G(Z0+Z1+Z2+Z3).

Then (3) and the exact complex channel reduction of
[star24](star24_free_rectangle_reduction.md), Sections 2-3, give

    C(B012)<=13 C0 D0 sqrt(pq).                       (8)

For clarity, the channel operators Omega=J_new tensor J_ij are real
symmetric, D=sum |vec Omega><vec Omega|, and each channel contributes
twice a free Gram. Its fields are

    alpha_Omega=c^T Omega d,
    E_Omega=sum_z (Omega d)_z E_z,
    F_Omega=sum_z (Omega^T c)_z F_z.

The transposes are essential for complex scalar vectors. The exact
moments obey

    sum ||E_Omega||^2<=D0^2 p,
    sum ||F_Omega||^2<=C0^2 q.

Concurrence subadditivity, its homogeneity, and Cauchy--Schwarz apply
(3) with the overall factor two to yield (8). This includes H of
arbitrary finite dimension and all zero scalar amplitudes; no divisions
by a channel amplitude are used.

To include the fourth component use the general purification bound

    C(G(X+Y))<=C(G(X))+2||X||||Y||+||Y||^2.

The actual component trace bounds, proved in
[star21](star21_concurrence_rectangle.md), are

    ||Z0||<=2sqrt(2)C0 D0, ||Z1||<=2D0 sqrt(p),
    ||Z2||<=2C0 sqrt(q),   ||Z3||<=sqrt(2pq).

Consequently the universal full-field bound is

    C(B)<=sqrt(pq)[21 C0 D0
                  +4sqrt(2)(D0 sqrt(p)+C0 sqrt(q))
                  +2sqrt(pq)].                       (9)

Here 21=13+8: the extra eight is the scalar/double-residual cross term
and must be retained. Both-outside fields, unequal planes, entangled
scalar sector vectors and nonzero residual zero-coordinates are all
included. A Schur matrix A contributes the usual exact factor |det A|
by local-filter covariance of two-qubit concurrence.

## Exact audit and remaining target

The companion checker verifies the product-compression partial
transpose and trace identities on complex rational nonprincipal data,
the exact complex product-completion identity, PSD/PPT of all product
Grams, both cutoff comparisons, and actual full-field channel and
purification identities. It also checks all four physical component
trace estimates on a datum where all four components are nonzero.
The general conclusions follow from (1)-(9), rather than from finite
controls. The checker performs no search or optimization.

An independent read of this proof and combined rerun of the bridge
checker are recorded in
the orb review (a working note not included here).

The sharp free coefficient two, equivalently physical coefficient
four for B012, remains the eleven-dimensional pencil target in
star24. The improved finite coefficient above does not settle it.
