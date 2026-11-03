# Free completion gives a uniform bridge coefficient 2129/10

This proof is retained as an intermediate result. The subsequent
[combined bridge theorem](orb28_combined_uniform_bridges.md) improves
the coefficient to 164, with the same general actual-star threshold m>=17,
and improves the entire rank-at-most-one-Y subclass to m>=12.

For every actual rank-at-most-two block operator C=[X,P;Q,Y] in the
one-free-copy merger setting, with arbitrary finite previous space H,

    sqrt(q_1(Y))[sqrt(Gamma(X))-sqrt(Gamma_bar(X))]
       <=(2129/10)sqrt(qtilde(P)qtilde(Q)).               (1)

The exact energy-aware recursion consequently proves actual star for
every m>=17, and the one-free-copy, one-excitation endpoint class for
n=m+1>=18. This is a proper class. General star at m=3 through 16,
the sharp merger coefficient four, ORB, and unrestricted endpoint
positivity and induction remain open.

This improves induction23 (a working note not included here)
without changing its canonical gauge, adjugate comparison, zero cases
or rank-preserving energy path. The new input is the sharp generic
coherent-sum theorem of [orb24](orb24_coherent_separable_sum.md):

    G(Z1),G(Z2) separable
       => C(G(Z1+Z2))<=2||Z1||||Z2||.                  (2)

It permits arbitrary relative purification environments and does not
require simultaneous product rows. The free-completion argument below
is derived in [star25](star25_free_completion_rectangle.md).

## 1. A free-only rectangle bound with coefficient 13/2

Let W=I-F_free/2 on C3 tensor C3. For arbitrary complex 3-by-2 R,L,

    B_(R,L)=(R* tensor L*)W(R tensor L)

is a separable two-qubit matrix. Indeed B is positive since W>0.
For every 2-by-2 test matrix T, its partial-transpose expectation is

    <vec T|B_(R,L)^Gamma|vec T>
       =||R T L*||^2-|Tr(R T L*)|^2/2>=0.

The last inequality follows from rank(R T L*)<=2 and
|Tr Z|<=||Z||_*<=sqrt(2)||Z||. Thus B is PPT and hence separable
in two-qubit dimensions, including all singular factors. Its trace is

    Tr B_(R,L)=||R||^2||L||^2-||R*L||^2/2
                  <=||R||^2||L||^2.                     (3)

Let V,U:C2 to C3 be arbitrary isometries. Let E,F be arbitrary complex
3-by-2 matrices with Tr(V*E)=Tr(U*F)=0, and let alpha be any complex
scalar. Write e=||E||, f=||F|| and

    Q=alpha V tensor U+E tensor U+V tensor F.

Then

    C(Q*WQ)<=(13/2)ef.                                 (4)

If e=0 or f=0, Q is a product operator, so (4) follows from (3).
Otherwise h=max(e,f)>0. First suppose |alpha|<=9h/10. If h=e,
split Q=(alpha V+E) tensor U+V tensor F. The two W-weighted
purifications have separable Grams by (3). Trace orthogonality gives
||alpha V+E||^2=2|alpha|^2+e^2. Equations (2)-(3) yield

    C(Q*WQ)<=4ef sqrt(1+2|alpha|^2/e^2)
             <=(2/5)sqrt(262)ef<(13/2)ef.              (5)

For h=f, group Q instead as E tensor U+V tensor (alpha U+F).
The same estimate holds with f in the denominator. No equality of
the planes or restriction on outside components is used.

If |alpha|>=9h/10, then alpha is nonzero. Complete the product:

    Q=(alpha V+E) tensor (U+F/alpha)-E tensor F/alpha.

Both terms again have separable W-weighted Grams, including the minus
sign on the second purification. Apply (2), (3) and trace orthogonality:

    C(Q*WQ)
      <=2ef sqrt[(2+e^2/|alpha|^2)(2+f^2/|alpha|^2)]
      <=(524/81)ef<(13/2)ef.                           (6)

This proves (4) for every complex alpha, including zero. The two strict
scalar comparisons are 16*262<65^2 and 2*524<13*81.

## 2. The actual full-field concurrence coefficient is 21

Keep the actual physical kernel and scalar projections from induction23:

    K_H=2W_free(I-F_new)(I-F_H),
    R=V tensor a+r, L=U tensor c+delta,
    alpha_r=||r||, beta=||delta||, A=||a||, C0=||c||.

Let Z0,Z1,Z2,Z3 be the four purification components with zero, left,
right and double residuals, respectively. The real antisymmetric
channel decomposition of [star24](star24_free_rectangle_reduction.md)
transfers any free coefficient b into a first-three-component physical
coefficient 2b. Specifically, for each real symmetric channel
Omega=J_new tensor J_ij, the free operator is

    Q_Omega=(a^T Omega c)V tensor U
                  +r_Omega tensor U+V tensor delta_Omega,
    r_Omega=sum_z(Omega c)_z r_z,
    delta_Omega=sum_z(Omega^T a)_z delta_z.

The exact Gram identity and channel moments are

    G(Z0+Z1+Z2)=2sum_Omega Q_Omega*W Q_Omega,
    sum_Omega||r_Omega||^2<=C0^2 alpha_r^2,
    sum_Omega||delta_Omega||^2<=A^2 beta^2.

These statements retain the complex transposes in the bilinear channel
coefficients. The moment bounds follow from the actual qubit source
0<=D_x<=||x||^2I, not from a generic positive-rest replacement.
Concurrence subadditivity, (4), and scalar Cauchy therefore give

    C(G(Z0+Z1+Z2))<=13AC0 alpha_r beta.                 (7)

For dim H<2 the source is zero and the same statement holds directly.
All zero scalar and residual cases are included without division.

The unchanged physical component bounds are

    ||Z0||<=2sqrt(2)AC0, ||Z1||<=2C0 alpha_r,
    ||Z2||<=2A beta,    ||Z3||<=sqrt(2)alpha_r beta.

Adding Z3 by the general HJW perturbation estimate proves

    C(B)<=alpha_r beta[21AC0
                 +4sqrt(2)(C0 alpha_r+A beta)+2alpha_r beta],
    B=(R* tensor L*)K_H(R tensor L).                    (8)

The leading coefficient is 13+8=21, including the scalar/double-residual
cross term. As before, H(R M L*)>=-|det M|C(B). The proof uses actual
free and qubit-source structure at both stages.

## 3. Tightened functional coefficient 146/5

The unscaled functional lemma of induction23 remains valid with 30
replaced by 146/5:

    sqrt(Gamma(F L*))-sqrt(Gamma_bar(F L*))
       <=(146/5)||f1||||delta||,
    F=[f1,f2], L=U tensor c+delta.                       (9)

Use exactly its compensated-curvature estimate and safe-axis cases,
but take tau0=151/20. The small-case coefficient becomes

    (171/20)(sqrt(2)+2)<146/5.

The sum of the two energy lower bounds in the large case has coefficient
(16361sqrt(2)-12880)/22801, and the curvature coefficient is 1972/151.
Their ratio is

    297772/(16361sqrt(2)-12880)<146/5.

The strict inequalities follow respectively from
242^2-2*171^2=82 and
2*2388706^2-3369340^2=59380673272. Thus (9) holds on the same full
rank-at-most-two range space, with no Y assumption or bound on f2.

## 4. Uniform bridge conversion and all degeneracies

For rank-two Y retain the exact canonical identities of induction23:

    Y=U diag(1,t)V*, M=diag(1,s), s=1/t,
    X=R M L*, P=R V*, Q=U L*,
    p=alpha_r^2, q=beta^2, y=q_1(Y)<=2,
    q_1(Y^+)>=s^2 y/2.

Take tau=47/8. If A sqrt(y)<=tau alpha_r or C0 sqrt(y)<=tau beta,
apply (9) or its adjoint. The scaled gap is at most

    (146/5)(tau+sqrt(2))alpha_r beta
       =(3431/20+(146/5)sqrt(2))alpha_r beta
       <(2129/10)alpha_r beta.                         (10)

The bound on the first column is independent of s, as in induction23.
The strict comparison is 827^2-2*584^2=1817>0.

Otherwise both scalar amplitudes are large. For
X0=(V tensor a)M(U tensor c)*, both energies are at least
s^2 y A^2 C0^2/2. Retain both energies in the denominator and use
the same two-step path X -> R M(U tensor c)* -> X0. Every step has
a common rank-at-most-two range or co-range. Hence

    sqrt(Gamma(X))+sqrt(Gamma_bar(X))
       >=s sqrt(y)AC0 (1329sqrt(2)-1632)/2209.           (11)

The denominator is strictly positive. For H(X)<0, use the determinant
bound, (8), and alpha_r/A,beta/C0<=sqrt(y)/tau. The scaled gap is
at most

    alpha_r beta [21+16/tau+4/tau^2]
          /[sqrt(2)-(sqrt(2)+2)(2/tau+sqrt(2)/tau^2)]
       =alpha_r beta 52661/(1329sqrt(2)-1632)
       <(2129/10)alpha_r beta.                         (12)

The last strict comparison is
2*2829441^2-4001138^2=2367449918>0. For H(X)>=0 the gap is
nonpositive. This proves (1) when rank Y=2 and y,p,q are nonzero.

All remaining cases are unchanged and exhaustive. If y=0, the left
side of (1) is zero. If p=0 or q=0, a scalar safe axis makes the gap
nonpositive. Rank Y<=1 is covered by the already proved coefficient
50, including singular B-sector factors by continuity. Scaling the
entire original C restores its pre-normalization singular value with
the same quadratic homogeneity on both sides. No varying-field limit
or ill-conditioned full-rank Y is excluded.

## 5. The m=17 base comparison and all larger m

The exact recurrence is now

    K_m<=sqrt(2)K_(m-1)+2129/10, K_1=0.

At m=17 it gives

    K_17<=(108579/2)(1+sqrt(2))<131070=2^17-2.

Indeed 153561^2-2*108579^2=2182239>0. For every m>=17,

    (3/2)(2^m-2)+2129/10
       =3*2^(m-1)+2099/10<=2^(m+1)-2,

because 2^(m-1)>=2119/10. This proves actual star for every m>=17
and only the specified endpoint class for n>=18. It does not settle
the user's all-n NPT target.

The companion checker verifies all new scalar comparisons, product-Gram
partial-transpose and trace identities on complex unequal-plane data,
both free-completion branches, the sharp coherent-sum identities and
actual physical channel controls, the rank-preserving energy path,
and the all-larger-m comparison. It performs no numerical search.
