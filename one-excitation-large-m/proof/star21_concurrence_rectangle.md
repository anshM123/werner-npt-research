# A uniform concurrence rectangle estimate for both scalar axes

This proves the concurrence interpolation estimate needed between the
positive and null scalar-baseline regimes. It applies to arbitrary
scalar sector vectors, unequal free planes, and simultaneous outside
components. It is a bound on the actual positive internal two-qubit
matrix; conversion to a square-root merger bridge is supplied by the
separate [energy argument](induction24_free_completion_uniform_bridge.md),
now with coefficient 2129/10 and actual-star consequence for every m>=17.
The sharper physical rectangle used there is proved by
[free completion](star25_free_completion_rectangle.md). This
note makes no unrestricted all-copy endpoint claim.

Use homogeneous two-qubit concurrence: on a pure coefficient matrix x
it is 2|det x|, and on positive matrices it is the convex roof of that
quantity. The standard common-environment HJW property lets one apply
the same padded environment isometry to several purifications. These
facts are given in the primary derivation of
[two-qubit concurrence](https://arxiv.org/abs/quant-ph/9709029).
All Frobenius norms and state traces below are unnormalized.

## 1. Two coherent-sum estimates

Write G(Z)=Z*Z for a purification with four internal columns. If
G(Z_1),G(Z_2) are separable, then

    C(G(Z_1+Z_2))
      <=min(z_1^2,z_2^2)+2z_1z_2<=3z_1z_2,
    z_i=||Z_i||_F.                                      (1)

Choose a common HJW environment making the rows of Z_1 a product
ensemble of its Gram state. For a product row e and any row f,
determinant polarization gives

    C(e+f)<=C(f)+2||e|| ||f||<=||f||^2+2||e|| ||f||.

Summing and applying scalar Cauchy proves the bound with z_2^2;
interchange the roles for the minimum. The final inequality uses
min(z_1^2,z_2^2)<=z_1z_2. Padding permits the common isometry even when
the required product ensemble has more terms than the original rows.

This common-environment step is elementary: if two factors Z,E have
Z*Z=E*E, the map taking each column of Z to the corresponding column
of E is a well-defined isometry of their column spans. After padding
the environments it extends to a unitary; applying that same unitary
to all other purifications preserves their norms, Grams and cross terms.

For arbitrary Z,W the same proof, starting with an optimal or
arbitrarily near-optimal concurrence ensemble of G(Z), gives

    C(G(Z+W))<=C(G(Z))+2||Z||_F||W||_F+||W||_F^2.        (2)

No assumption that the rows of W form an optimal ensemble is made.

## 2. Transport below a separable positive matrix

**Transport lemma.** If S is separable and positive definite,
0<=B<=S, D=S-B, and kappa=lambda_max(S)/lambda_min(S), then

    C(B)<=(4sqrt(kappa)+1) Tr D.                        (3)

Set

    E=S^(-1/2) D S^(-1/2), 0<=E<=I,
    Q=(I-E)^(1/2), T=S^(1/2) Q S^(-1/2).

Then T S T*=B. Choose a product ensemble S=sum_j e_j e_j*, and put
f_j=(T-I)e_j. The vectors e_j+f_j form a valid ensemble of B,
although T need not be a local filter. Let F=sum_j ||f_j||^2.
Since e_j are product vectors, the pure-state estimate in Section 1
gives

    C(B)<=2sqrt((Tr S)F)+F.

The exact transport error satisfies, with f(E)=I-sqrt(I-E),

    F=Tr[S f(E)^2]<=Tr[S E^2]
      <=||E||_op Tr(SE)=||E||_op Tr D.

Here f(E)^2<=E^2 and E^2<=||E||_op E are spectral inequalities;
their traces against positive S remain valid without commutativity.
Since ||E||_op<=min(1,Tr D/lambda_min(S)), both

    F<=Tr D,
    F<=(Tr D)^2/lambda_min(S)

hold. Using the second bound in the square-root term and the first in
the last term, together with Tr S<=4lambda_max(S), proves (3).
D=0 is included without division by Tr D.

The later sharp separable-majorant theorem (a working note not included here)
improves (3) to C(B)<=c(kappa)Tr D, with optimal c(kappa)=kappa/2
for 1<=kappa<=2 and min(kappa-1,2sqrt(kappa)) for kappa>=2.
The transport estimates below remain valid and are the ones used
in the existing coefficient-400 bridge; the new coefficient alone
does not establish the sharp rectangle or MB_4.

The proof also gives a sharper estimate when ||E||_op<=ell^2<1.
Put a=1/(1+sqrt(1-ell^2)); then f(E)<=a E and hence

    F<=a^2 (Tr D)^2/lambda_min(S),
    F<=a^2 ell^2 Tr D.                                 (3a)

These two inequalities, rather than a condition-number factor on F,
will be used in the rectangle estimate.

## 3. Abstract two-axis rectangle theorem

Suppose Z_0,Z_1,Z_2 share a common environment and have four internal
columns. Assume

    G(Z_0+tZ_1) and G(Z_0+tZ_2) are separable for every complex t.

Suppose additionally that either Z_0=0, or G(Z_0) is positive definite
with condition number at most three. Then

    C(G(Z_0+Z_1+Z_2))<=36||Z_1||_F||Z_2||_F.             (4)

If z_1=0 or z_2=0, this follows from the assumed separable axis.
If Z_0=0, it follows from (1). Otherwise put

    z_i=||Z_i||_F, r=z_1+z_2,
    s_0=sigma_min(Z_0)>0.

The limiting Gram states G(Z_1),G(Z_2) are separable because the
separable cone is closed and the two entire axes are separable.

Put eta=21/50. First suppose r<=eta s_0. Set theta=z_1/r,
and take the separable state

    S=theta G(Z_0+Z_1/theta)
        +(1-theta)G(Z_0+Z_2/(1-theta)).

With B=G(Z_0+Z_1+Z_2), direct expansion gives

    S-B=G(sqrt(z_2/z_1)Z_1-sqrt(z_1/z_2)Z_2)>=0,
    Tr(S-B)<=4z_1z_2.                                  (5)

Both perturbations Z_1/theta and Z_2/(1-theta) have Frobenius norm r.
The singular-value triangle inequalities therefore put each axis
Gram, and hence S, between

    (1-eta)^2 s_0^2 I and (sqrt(3)+eta)^2 s_0^2 I.

The Gram factor on the right of (5) has Frobenius norm at most
2sqrt(z_1z_2)<=r. Thus ||D||_op<=r^2 and

    ||E||_op<=ell^2, ell=eta/(1-eta)=21/29,
    sqrt(kappa(S))<=(sqrt(3)+eta)/(1-eta),
    a=1/(1+sqrt(1-ell^2))=29/49.

Apply (3a), Tr S<=4lambda_max(S), and Tr D<=4z_1z_2 to get

    C(B)<= [16a(sqrt(3)+eta)/(1-eta)+4a^2ell^2]z_1z_2
         =[(39200sqrt(3)+18228)/2401]z_1z_2
         <36z_1z_2.

If instead r>eta s_0, then z_0<=2sqrt(3)s_0
<2sqrt(3)r/eta. Apply the first, unsymmetrized determinant estimate
in Section 1 to the separable pairs (Z_0+Z_1,Z_2) and
(Z_0+Z_2,Z_1). It gives

    C(B)<=min[z_2^2+2(z_0+z_1)z_2,
              z_1^2+2(z_0+z_2)z_1]
        <=3z_1z_2+2z_0 min(z_1,z_2)
        <=[3+8sqrt(3)/eta]z_1z_2
         =[3+(400/21)sqrt(3)]z_1z_2
         <36z_1z_2,

where (z_1+z_2)min(z_1,z_2)<=2z_1z_2 was used. This completes the uniform transition
through every relative size of the baseline and the two perturbations.

## 4. Application to the actual free-weighted merged kernel

Let R,L be two-column fields on C^3_free tensor C^2_new tensor H,
written

    R_a=c_a V+E_a, L_b=d_b U+F_b,
    U*U=V*V=I_2,
    Tr(V*E_a)=Tr(U*F_b)=0,
    p=||E||_F^2, q=||F||_F^2, C_0=||c||, D_0=||d||.

The isometries U,V may have different free planes. No support,
commutation, or outside-component restriction is imposed on E,F.
Use the actual positive physical operator

    K=2(I-F_free/2)(I-F_new)(I-F_H),
    B=(R* tensor L*)K(R tensor L).

The corresponding actual forms satisfy
Gamma_bar(RAL*)-Gamma(RAL*)>=-|det A|C(B), as in the internal
concurrence reduction. The following theorem concerns C(B) itself.

**Full-field rectangle bound.** At every finite scalar amplitude,

    C(B)<=152 C_0 D_0 sqrt(pq)
          +4sqrt(2)(D_0 sqrt(p)+C_0 sqrt(q))sqrt(pq)
          +2pq.                                        (6)

To prove it, split the purification K^(1/2)(R tensor L) into

    Z_0 =K^(1/2)((cV) tensor (dU)),
    Z_1 =K^(1/2)(E tensor (dU)),
    Z_2 =K^(1/2)((cV) tensor F),
    Z_3 =K^(1/2)(E tensor F).

Both entire axes G(Z_0+tZ_1),G(Z_0+tZ_2) are separable. Indeed the
unchanged factor has a fixed sector vector. The actual scalar-axis
comparison makes the two-qubit Gram and its partial transpose
positive; the two-qubit PPT criterion gives separability.

One can see this directly before the remaining local compression.
Conditioning the sector operator (I-F_new)(I-F_H)>=0 on its fixed
sector vector gives a positive operator on the other sector. Compressing
one free leg of I-F_free/2 by a two-column isometry gives a 2-by-3
positive operator whose partial transpose is I-(1/2)|v><v|, with
||v||^2=2, hence positive. Tensoring the conditioned sector operator
and applying the remaining two-column local filter preserves this PPT
property. The final internal state is two qubits, so it is separable.

For the scalar baseline put

    h=2<c tensor d|(I-F_new)(I-F_H)|c tensor d> >=0.

Then G(Z_0) is h times the compression of I-F_free/2 by V tensor U.
Its eigenvalues lie in [h/2,3h/2]. When h>0 it is full rank with
condition number at most three. It is also separable: its partial
transpose has the form h[I-(1/2)|k><k|], with
||k||^2=||V*U||_F^2<=2. When h=0, positivity gives Z_0=0. Thus all
h, including zero, satisfy the hypotheses of (4).

The pure weighted-kernel estimate from
[star20](star20_weighted_qubit_kernel.md), summed over internal
columns, gives

    z_0<=2sqrt(2) C_0 D_0,
    z_1<=2D_0 sqrt(p), z_2<=2C_0 sqrt(q),
    z_3<=sqrt(2pq).                                     (7)

Equation (4) bounds the concurrence of G(Z_0+Z_1+Z_2) by
144 C_0D_0 sqrt(pq). Add Z_3 using (2) and
||Z_0+Z_1+Z_2||<=z_0+z_1+z_2. The extra cross term involving Z_0
is at most 8C_0D_0 sqrt(pq); the other terms are exactly bounded
by the remaining terms in (6). This explains the coefficient 152.

If p=0 or q=0, the actual state lies on a separable scalar axis and
both sides of (6) vanish. If C_0=0 or D_0=0, the proof still applies
with Z_0 and the corresponding first-order component zero. No division
by any of these quantities occurs in the final statement.

Finally, a local filter on the first internal qubit gives

    C(B_(RA,L))=|det A|C(B_(R,L)).                        (8)

For invertible A this follows directly from the pure determinant
formula and the bijection between decompositions before and after the
filter. Singular A follows by continuity. Thus (6) is uniform in the
Schur matrix condition number after the exact determinant factor in
(8) is retained. No energy denominator has been discarded.

## Exact audit and scope

`star21_concurrence_rectangle.py` checks the complete symbolic identity
(5), a noncommuting rational transport example with an explicit
product ensemble, and actual positive- and null-baseline fields with
unequal free planes and outside components. It verifies the physical
purification Gram, the scalar-baseline eigenvalue bounds, PPT axes,
and all component trace bounds with exact arithmetic.

The finite checks audit the identities used in this dimension-independent
proof; they are not used to infer (3), (4), or (6). A square-root merger
bridge is proved in star22 and improved in induction24 using the separate
energy argument. Actual star at m=3 through 16, ORB, unrestricted endpoint induction, and the
general all-copy NPT endpoint remain open.
