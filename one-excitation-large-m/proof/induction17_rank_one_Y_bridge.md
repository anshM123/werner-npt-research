# A uniform merged bridge for every rank-one Y

This coefficient-50 proof and its exact gauge are retained as intermediate
results. The [combined bridge theorem](orb28_combined_uniform_bridges.md)
uses the new functional estimate to give coefficient 22sqrt(2)<32 and
actual star for every rank-at-most-one-Y datum at every m>=12.

For every actual rank-at-most-two block operator C=[X,P;Q,Y] with
rank(Y)<=1, the merged forms of [star19](star19_energy_aware_recursion.md)
satisfy

    sqrt(q_1(Y))[sqrt(Gamma(X))-sqrt(Gamma_bar(X))]
       <=50 sqrt(qtilde(P)qtilde(Q)).                     (1)

The constant is independent of the previous dimension, free planes,
scalar parts, and simultaneous outside components. It is not claimed
sharp. The energy-aware recursion then proves actual star for this
entire rank-one-Y class at every m>=13. The later
[universal merger theorem](star22_uniform_merger_bridge.md) covers either
rank of Y at every m>=19, using the unscaled estimate proved here.
General smaller m, ORB and the all-copy endpoint remain open.

The proof uses qualitative two-copy endpoint positivity and the actual
qubit merger kernels. It uses no numerical search or generic rest
replacement. An independent review is in
orb19 (a working note not included here).

## 1. Actual forms, energy bounds and two safe axes

Let K=C2_new tensor H, with H arbitrary finite-dimensional, and put

    H(X)=Gamma_bar(X)-Gamma(X).

Here H(X) denotes a quadratic form; H without an argument denotes the
previous Hilbert space. Its positive physical replica kernel is

    K_H=2W_free(I-F_new)(I-F_H),  W_free=I-F_free/2.

For rank-one E=ab* and F=cd*, the crossed-index identity is

    H(E,F)=<a tensor d,K_H(c tensor b)>.

The sharp pure-state estimate from
[star20](star20_weighted_qubit_kernel.md) is

    0<=H(E)<=2||E||_F^2,
    |H(E,F)|<=2||E||_F||F||_F.                           (2)

For a self-contained proof, set A0=W_free(I-F_H)>=0 on two replicas
of C3_free tensor H. For vectors x,w on that space, expansion gives

    <x tensor w,A0(x tensor w)>
      =||x||^2||w||^2-Tr(rho_(x,H)rho_(w,H))
          -Tr(rho_(x,free)rho_(w,free))/2+|<x,w>|^2/2.

Both reduced overlaps are squared Frobenius norms of partial traces
of xw*. Each such partial trace has rank at most three and trace
<w,x>. Thus each overlap is at least |<x,w>|^2/3, proving

    ||A0^(1/2)(x tensor w)||<=||x||||w||.

Write u=|0>x+|1>y and v=|0>z+|1>w in the new-qubit basis. Its
antisymmetric subspace is one-dimensional, so exactly

    H(uv*)=2||A0^(1/2)(x tensor w-y tensor z)||^2
       <=2(||x||||w||+||y||||z||)^2<=2||u||^2||v||^2.

Cauchy--Schwarz for the positive physical K_H, applied to the crossed
pairs, proves the second bound in (2). The operator norm of K_H itself
is not asserted to be two.

On every rank-at-most-two input Z the actual energies satisfy

    0<=Gamma(Z)<=2||Z||_F^2,
    0<=Gamma_bar(Z)<=4||Z||_F^2.                          (3)

Positivity of Gamma is the two-copy endpoint theorem; that of
Gamma_bar follows from its positive qubit product-frame sum of
rank-two previous merged Gamma forms. For the first upper bound,

    Gamma(Z)=2||Z||^2-||Tr_free Z||^2-||Tr_sector Z||^2
                 +|Tr Z|^2/2.

The partial-trace ranks are at most six and three respectively.
Their squared norms are therefore at least |Tr Z|^2/6 and
|Tr Z|^2/3. For the second, insert the product frame for 2I-F_new
from star19 and use the first bound and the exact norm moment

    sum_j w_j||Z_(u_j,v_j)||^2
       =2||Z||^2-||Tr_new Z||^2<=2||Z||^2.

For rank-one E there is also

    Gamma(E)>=||E||_F^2/2,                               (4)

since its physical kernel W_free(2I-F_new F_H) is at least I/2.
On any fixed two-dimensional range space, Gamma and Gamma_bar are
nonnegative quadratic forms, so their square roots are seminorms.
All triangle inequalities below stay within one such space.

We need one further safe family:

    H(Z)>=0 if ran Z* lies in U tensor c,
             dim U<=2, c in K fixed.                     (5)

Condition the positive sector kernel (I-F_new)(I-F_H) on c in the
right replica. Diagonalizing that positive contraction expresses
H(Z) as a positive sum of q_1-values of free matrices whose
co-ranges lie in U. Each has rank at most two, so each value is
nonnegative. This is the merger version of Lemma R-prime. The
vector c may be entangled across new and H. Positivity holds on
the entire linear space in (5).

## 2. A gauge retaining the exact scalar cancellation

Factor C=[R_A;R_B][L_A;L_B]* with two columns. Y=0 makes (1)
immediate. If rank Y=1, at least one B-sector factor has rank one.
Take the adjoint if necessary so that rank L_B=1.

First suppose R_B has column rank two. An invertible internal change
of factor basis makes R_B an isometry. A subsequent internal unitary
puts L_B=[v,0], preserving that isometry. A common physical free
unitary then gives

    R_B=[e_0,e_1], L_B=[v,0],
    R_A=[x_1,x_2], L_A=[y_1,y_2],
    Y=e_0v*, P=x_1v*, Q=e_0y_1*+e_1y_2*,
    X=x_1y_1*+x_2y_2*.                                  (6)

Internal changes preserve every block. The common physical unitary
preserves all forms and costs. With y_(i,a) denoting free component
a of y_i, define

    c=(y_(1,0)+y_(2,1))/2,
    y_1=e_0 tensor c+delta_1,
    y_2=e_1 tensor c+delta_2,
    delta_(1,0)+delta_(2,1)=0.

Then the potentially large scalar part of Q cancels exactly:

    q:=qtilde(Q)
      =||y_1||^2+||y_2||^2-||y_(1,0)+y_(2,1)||^2/2
      =||delta_1||^2+||delta_2||^2=:beta^2.                (7)

Put alpha=||x_1||, T=||x_2|| and S=||c||. Every free block of P
has rank at most one, giving

    p:=qtilde(P)>=||v||^2 alpha^2/2,
    q_1(Y)<=||v||^2.                                    (8)

No lower bound by half the full norm is assumed for qtilde(Q).

## 3. Compensated curvature and the square-root bound

Decompose X=A+B+C0+D0, where

    A=x_2(e_1 tensor c)*, B=x_1(e_0 tensor c)*,
    C0=x_2 delta_2*,      D0=x_1 delta_1*.

Their ranges all lie in span{x_1,x_2}. Also

    ||A||=TS, ||B||=alpha S,
    ||C0||<=T beta, ||D0||<=alpha beta.

The form H is nonnegative on span{A,B} by (5), and on span{A,C0}
because every combination there has rank at most one. Cauchy--Schwarz
on those two restrictions, followed by (2), gives

    H(A+B+C0)>=-2sqrt(H(B)H(C0))-2|H(B,C0)|
                  >=-8alpha beta TS.

For clarity, the first inequality follows by writing
a=sqrt(H(A)), b=sqrt(H(B)), c0=sqrt(H(C0)) and using

    a^2+b^2+c0^2-2a(b+c0)
        =(a-b-c0)^2-2bc0>=-2bc0.

Add D0, retain H(D0)>=0, and apply (2) to its three crossed terms:

    H(X)>=-12alpha beta TS-4alpha beta(alpha S+beta T).   (9)

If alpha=0, X has rank at most one and H(X)>=0. If beta=0,
condition (5) applies. Suppose both are positive. Three exhaustive
cases now bound the square-root gap.

If T<=9alpha, compare X with the safe axis X0=A+B. The error
E=C0+D0 has norm at most 10alpha beta. The seminorm triangles
on the common range space and (3) give

    sqrt(Gamma(X))-sqrt(Gamma_bar(X))
       <=sqrt(Gamma(E))+sqrt(Gamma_bar(E))
       <=10(sqrt(2)+2)alpha beta.                        (10)

If S<=9beta, use the rank-one axis X0=A+C0 and error E=B+D0.
It has the same norm bound, proving (10) again.

Otherwise T>=9alpha and S>=9beta. Equations (3)-(4) and the reverse
triangle on the same range space give

    sqrt(Gamma(X))
       >=TS/sqrt(2)-sqrt(2)(alpha S+T beta+alpha beta)
       >=43TS/(81sqrt(2)).                              (11)

Equation (9) gives -H(X)<=(116/9)alpha beta TS. If H(X)>=0 there
is nothing to prove. Otherwise divide by
sqrt(Gamma(X))+sqrt(Gamma_bar(X)), which is at least (11), to get

    sqrt(Gamma(X))-sqrt(Gamma_bar(X))
       <=(1044sqrt(2)/43)alpha beta.                     (12)

Finally multiply (10) or (12) by sqrt(q_1(Y)) and use (7)-(8).
The respective coefficients of sqrt(pq) are

    20+20sqrt(2) and 2088/43,

both strictly less than 50. This proves (1) in the nonsingular-factor
case. The relative sizes of T and S have never been bounded.

## 4. Singular factors, recurrence and scope

If R_B is singular, first put L_B=[v,0] by an internal basis change.
Its first R_B column is nonzero since Y is nonzero. Perturb only the
second R_B column to get column rank two. This leaves Y, X and P
unchanged, keeps rank C at most two, and makes Q converge to its
original value. Apply (1) and pass to the limit. The constants are
independent of the factor gauge and its condition number. This also
includes rank-one full C. Taking adjoints covers the other initial
choice of rank-one B-sector factor.

The contractions in star19 preserve Y. Restricting its energy-aware
recursion to rank(Y)<=1 and applying (1) gives

    K_m<=sqrt(2)K_(m-1)+50, K_1=0,
    K_m<=50(1+sqrt(2))[2^((m-1)/2)-1].

At m=13 this is

    K_13<=3150(1+sqrt(2))<8190=2^13-2,

using sqrt(2)<8/5. The comparison continues for all larger integers:
if K_m<=2^m-2, then sqrt(2)<3/2 gives
K_(m+1)<=3*2^(m-1)+47<=2^(m+1)-2 whenever 2^(m-1)>=49.
Thus actual star holds for every rank-one-Y datum at every m>=13.

`induction17_rank_one_Y_bridge.py` audits the exact gauge and block
identities, Q-cost cancellation, both safe-axis Gram matrices,
physical crossed contractions, sharp pure-kernel estimates, and all
constants and recurrence comparisons. It uses exact arithmetic and
no optimization. Finite checks audit the displayed proof; they do
not replace its arbitrary-dimension argument.

The rank-two-Y case has no rank-one B-sector factor in general, so
(6)-(8) do not apply to it. Its finite merger bound is instead proved
by the separate scalar-rectangle and adjugate argument in star22.
This result leaves m=3 through 12 open for general rank-one Y;
neither theorem proves all-n EP.
