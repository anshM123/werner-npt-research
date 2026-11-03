# A full-safe-axis functional estimate with coefficient 22

The later [two-flat-margin theorem](orb37_flat_margin_functional_bridge.md)
improves this coefficient to 163/10 and proves merger 125. The proof
below remains a valid source of its full-purification curvature bound.

Use the actual merger forms Gamma,Gamma_bar on operators on
C3_free tensor C2_new tensor H, with arbitrary finite H, as in
[star19](star19_energy_aware_recursion.md). Write H(X)=Gamma_bar(X)-Gamma(X).
For arbitrary two-column F=[f1,f2], an isometry U:C2 to C3, a sector
vector c, and an arbitrary two-column delta, put

    L0=U tensor c, L=L0+delta, X=FL*,
    alpha=||f1||, T=||f2||, S=||c||, beta=||delta||.

Then

    sqrt(Gamma(X))-sqrt(Gamma_bar(X))<=22alpha beta.   (1)

No bound on f2, no scalar projection condition on delta, and no
orthogonality of the columns of F is assumed. This improves the
functional coefficient 146/5 in induction24. It does not by itself
prove MB_4, a new actual-star copy threshold, ORB, or general EP.

## 1. Inputs and the valid full-purification curvature bound

The established actual estimates are

    0<=Gamma(Z)<=2||Z||^2,
    0<=Gamma_bar(Z)<=4||Z||^2,                  rank Z<=2,
    H(ab*)>=0,  H(ab*)<=2||a||^2||b||^2.             (2)

On a fixed range of dimension at most two, each square root is a
seminorm. The scalar co-range axis has H(F L0*)>=0. More strongly,
its internal two-qubit physical Gram is separable; the same scalar
co-range property holds after every internal determinant-one filter.
These are the actual-kernel results used in induction23 and star20.

For completeness, retain the entire physical purification. Split
F=R0+R1, with R0=[0,f2], R1=[f1,0], and let

    ZA=Z(R0,L0), ZB=Z(R1,L0),
    ZC=Z(R0,delta), ZD=Z(R1,delta),
    Z(R,L)=sqrt(K)(R tensor L),
    K=2W_free(I-F_new)(I-F_H).

The first two components have a separable joint Gram. ZA,ZC have
first internal-qubit support |1>; ZB,ZD have support |0>. The
row-supported concurrence lemma (a working note not included here)
therefore gives, in two **full-purification** steps,

    C(G(ZA+ZB+ZC+ZD))
      <=2||ZB||||ZC||+2||ZA+ZC||||ZD||
      <=2zB zC+2(zA+zC)zD.                           (3)

The row lemma applies also when its starting Gram is entangled, so
the second step is legitimate. The rank-one bound in (2), summed
over the two internal columns on the other side, gives

    zA<=2TS, zB<=2alpha S,
    zC<=sqrt(2)T beta, zD<=sqrt(2)alpha beta.

The identity-matrix concurrence witness is the actual value H(X).
Consequently (3) implies

    -H(X)<=8sqrt(2)alpha beta TS+4alpha T beta^2.       (4)

Nothing is discarded when an operator product such as R1delta*
vanishes: ZD remains in the physical purification. Nor does reversing
the addition order justify replacing the last term by an unproved
minimum of asymmetric expressions. The exact crossed-index example
in the row-lemma note is included in the checker below.

## 2. Compare with the entire safe operator

Set Xsafe=F L0*, E=Fdelta*, so X=Xsafe+E. All three operators have
range inside ran F, of dimension at most two. Thus both seminorm
triangles are valid on this **single common range**, and

    ||E||<=||F||_op beta<=sqrt(T^2+alpha^2)beta.        (5)

The proved two-factor flat margin from
[induction16, Section 2](induction16_qubit_merged_gap.md) is

    Gamma(Z)>= (sigma1(Z)-sigma2(Z))^2/2,
                                                   rank Z<=2.

Since L0*L0=S^2 I2, the singular values of Xsafe are S times those
of F. Moreover

    (sigma1(F)-sigma2(F))^2
      =T^2+alpha^2-2sqrt(det(F*F))>=(T-alpha)^2.

For T>=alpha, the safe-axis property therefore gives

    sqrt(Gamma(Xsafe)),sqrt(Gamma_bar(Xsafe))
                              >=S(T-alpha)/sqrt(2).

Apply the reverse triangle to both energies and use (2),(5):

    sqrt(Gamma(X))+sqrt(Gamma_bar(X))
      >=sqrt(2)S(T-alpha)
                      -(sqrt(2)+2)sqrt(T^2+alpha^2)beta. (6)

This is a comparison with the full rank-two safe operator, not the
single rank-one summand f2(u1 tensor c)*. It avoids the extra
alpha*S error of that earlier comparison. No potentially rank-four
error is passed through a seminorm triangle.

## 3. Exhaustive cases and exact coefficient comparison

If alpha=0, X has rank at most one, and if beta=0 it lies on the
scalar safe axis. Both give a nonpositive gap directly. Otherwise
take the following three cases.

**Small T:** T<=(63/10)alpha. Compare X to Xsafe. Its gap is
nonpositive, so the two triangles and (5) give

    gap<= (sqrt(2)+2)sqrt(T^2+alpha^2)beta
       <=(sqrt(2)+2)sqrt(4069)alpha beta/10
       <22alpha beta.                               (7)

**Small S:** S<=(27/5)beta. Compare X to the rank-one safe operator

    Xrank1=f2(u1 tensor c+delta2)*.

Its error is f1(u0 tensor c+delta1)*, whose norm is at most
alpha(S+beta). Both operators again have range inside ran F. Hence

    gap<=(sqrt(2)+2)alpha(S+beta)
        <=(32/5)(sqrt(2)+2)alpha beta<22alpha beta.    (8)

**Both large:** T>(63/10)alpha and S>(27/5)beta. Put
u=alpha/T<=10/63 and v=beta/S<=5/27. Equation (6) bounds the
energy sum below by TS d, where

    d=sqrt(2)(1-10/63)
           -(sqrt(2)+2)(5/27)sqrt(1+(10/63)^2)
     =[1431sqrt(2)-5(sqrt(2)+2)sqrt(4069)]/1701.

If H(X)>=0 the gap is nonpositive. Otherwise divide (4) by the
positive denominator from (6), obtaining

    gap<=alpha beta [8sqrt(2)+20/27]/d<22alpha beta.  (9)

Here all strict numerical comparisons have short rational certificates.
The bounds

    140/99<sqrt(2)<99/70, sqrt(4069)<319/5

follow by squaring positive integers. They imply the two small-case
coefficients are less than 76241/3500 and 3824/175, respectively,
both below 22. For the large case they give

    1701[22d-(8sqrt(2)+20/27)]
      =17874sqrt(2)-1260-110(sqrt(2)+2)sqrt(4069)
      >17874(140/99)-1260-110(239/70)(319/5)
       =21139/385>0.

This also proves d>0. Thus no unsupported numerical spectral bound
or scalar optimization is used. The cases exhaust all data, including
T=0,S=0, and prove (1).

## Verification and scope

The companion exact checker audits the row-supported determinant
identity, the valid full purification and all four component norms,
the crossed-index control with R1delta*=0 but ZD nonzero, the
full-safe singular-value and common-range energy inequalities, and
the three rational coefficient certificates. There is no search.

The theorem is an unconditional functional estimate for the actual
forms, in arbitrary finite H and without scalar-residual projection
assumptions. Inserting 22 into a uniform bridge still requires its
separate canonical two-case conversion; no new bridge constant or
copy threshold is claimed here. In particular the unresolved
negative interval in the general EP/ORB problem is not excluded.
