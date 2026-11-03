# Uniform bridge 164, and bridge 32 when rank Y is at most one

The later [flat-margin conversion](orb37_flat_margin_functional_bridge.md)
improves the general bridge to 125 and the rank-one-Y bridge to
(163/10)sqrt(2)<24, proving actual star for m>=16 and m>=11,
respectively. The valid earlier proof below supplies its canonical
gauges, field estimate and degeneracy treatment.

For every actual rank-at-most-two block operator C=[X,P;Q,Y] in the
one-free-copy merger setting, with arbitrary finite previous space H,

    sqrt(q1(Y))[sqrt(Gamma(X))-sqrt(Gamma_bar(X))]
          <=164sqrt(qtilde(P)qtilde(Q)).               (1)

For the full subclass rank(Y)<=1, a stronger bound holds:

    sqrt(q1(Y))[sqrt(Gamma(X))-sqrt(Gamma_bar(X))]
          <=22sqrt(2)sqrt(qtilde(P)qtilde(Q))
          <=32sqrt(qtilde(P)qtilde(Q)).                (2)

Consequently actual star holds for every m>=17 in the general class,
and for every m>=12 in the entire rank-at-most-one-Y class. The
corresponding one-free-copy, one-excitation endpoint classes have
n=m+1>=18 and n>=13, respectively. These are proper classes; neither
statement proves unrestricted EP, ORB, general induction or MB_4.
The first copy threshold is unchanged from induction24.

The new inputs are the actual full-safe-axis functional estimate 22
from [orb27](orb27_full_safe_functional.md) and the fourth-component
concurrence bound C(B3)<=pq from
[orb26](orb26_arbitrary_gram_coherent_triangle.md). The canonical
factor gauges and actual energy recursion are the ones already
proved in induction17, induction24 and star19.

## 1. Every rank-at-most-one Y: exact gauge and all degeneracies

The proof in [induction17, Section 2](induction17_rank_one_Y_bridge.md)
applies without any change of hypotheses. Factor the full block with
two columns, take an adjoint if necessary, and suppose the B-sector
right factor has rank one. If the other B-sector factor has column
rank two, its whitening, an internal unitary, and a common physical
free unitary put the factors exactly in the gauge

    Y=e0 v*, P=x1 v*, Q=e0 y1*+e1 y2*,
    X=x1y1*+x2y2*,
    y1=e0 tensor c+delta1, y2=e1 tensor c+delta2.

The scalar c is the orthogonal projection used there, so, with
alpha=||x1|| and beta^2=||delta1||^2+||delta2||^2,

    q:=qtilde(Q)=beta^2,
    p:=qtilde(P)>=||v||^2 alpha^2/2,
    q1(Y)<=||v||^2.                                  (3)

Only P is given the rank-one trace coercivity bound. No half-norm
bound is assumed for the generally rank-two blocks of Q. The
unscaled functional22 theorem applies to this exact F=[x1,x2] and
L=U tensor c+delta. Multiplying by sqrt(q1(Y)) and using (3) gives

    scaled gap<=22||v||alpha beta<=22sqrt(2)sqrt(pq),

which proves the first inequality in (2). The second follows, for
example, from sqrt(2)<99/70, giving 22sqrt(2)<1089/35<32.

The remaining factor cases are the same as in induction17, Section 4.
Y=0 gives zero on the left. If the other B-sector factor is singular,
first put the rank-one factor in the form [v,0]. Its partner's first
column is nonzero when Y is nonzero. Perturb only that partner's
second column to obtain column rank two. X,P,Y remain fixed; Q and
its cost converge, and the full operator continues to have rank at
most two. Apply the uniform bound and pass to the limit. No bounded
inverse or bounded whitening gauge is assumed. Adjoint covariance
covers the other initial factor orientation. Zero costs are included
directly in the functional theorem and its safe axes.

## 2. Full-rank Y: the unchanged exact canonical reduction

For rank Y=2, scale the entire original block by the reciprocal of
the largest singular value of Y. Both sides of (1) have quadratic
homogeneity. The exact Schur representation of induction24 is

    Y=U diag(1,t)V*, 0<t<=1, M=diag(1,s), s=1/t,
    X=RML*, P=RV*, Q=UL*,
    R=V tensor a+r, L=U tensor c+delta,
    A=||a||, D=||c||, alpha=||r||, beta=||delta||,
    p=alpha^2, q=beta^2, y=q1(Y)<=2,
    q1(Y+)>=s^2 y/2.                                 (4)

The last comparison uses the unrestricted-complex principal-plane
adjugate identity; the planes need not coincide and their singular
bases need not be principal bases. No lower bound on t is used.

The exact full-field estimate from orb26 is

    C(B)<=alpha beta[21AD
                  +4sqrt(2)(Dalpha+Abeta)+alpha beta]. (5)

In particular its last term is pq, not 2pq. Its proof retains the
full fourth purification component and the actual kernel K.

If y=0, (1) is immediate. If alpha=0 or beta=0, the scalar safe axis
makes the gap nonpositive. Remove only these cases before divisions.

## 3. The two full-rank cases give 164

Set tau=301/50. If A sqrt(y)<=tau alpha, apply functional22 to
F=RM. Its first column is R1, with ||R1||<=A+alpha independent
of s. Hence the scaled gap is at most

    22(tau+sqrt(2))alpha beta<164alpha beta.           (6)

If D sqrt(y)<=tau beta, use the adjoint argument; the same estimate
holds. This does not bound either second column.

Otherwise both amplitudes are large. Retain both actual energies
and the same path as induction24,

    X -> RML0* -> R0ML0*,
    R0=V tensor a, L0=U tensor c.

Each step has a common rank-at-most-two range or co-range. The
whole error may have rank four, so no direct triangle is applied
to it. The energy denominator is bounded below by

    sqrt(Gamma(X))+sqrt(Gamma_bar(X))
      >=s sqrt(y)AD d_tau,
    d_tau=sqrt(2)-(sqrt(2)+2)(2/tau+sqrt(2)/tau^2)
         =(55501sqrt(2)-65200)/90601>0.               (7)

For H(X)>=0 the desired gap is nonpositive. Otherwise use
H(X)>=-s C(B), (5), alpha/A,beta/D<=sqrt(y)/tau, and y<=2.
This gives

    scaled gap<=alpha beta[21+16/tau+2/tau^2]/d_tau
      =alpha beta 2148421/(55501sqrt(2)-65200)
      <164alpha beta.                               (8)

The exact comparisons require only 140/99<sqrt(2)<99/70. They imply

    22(tau+sqrt(2))<28622/175<164,
    164(55501sqrt(2)-65200)-2148421
       =9102164sqrt(2)-12841221>3022081/99>0.

The latter proves positivity of (7) as well. These bounds include
arbitrary complex fields, entangled scalar sector vectors, both
outside components, and arbitrarily ill-conditioned full-rank Y.
Together with Section 1 and the zero cases they prove (1).

## 4. The two exact recurrences and their scopes

The [energy-aware recursion](star19_energy_aware_recursion.md) has
unchanged base K1=0. With (1),

    Km<=sqrt(2)K_(m-1)+164,
    K17<=41820(1+sqrt(2))<104550<131070=2^17-2.

If Km<=2^m-2, the comparison at m+1 follows from

    (3/2)(2^m-2)+164=3*2^(m-1)+161<=2^(m+1)-2

when 2^(m-1)>=163, automatic for m>=17. Thus (1) proves actual
star at all m>=17. This improves the bridge coefficient but not the
previous general m-threshold.

For (2), the recursion contracts only the A sector and **preserves
Y**, hence preserves rank(Y)<=1 at every step. Using its rounded
coefficient32 gives

    Km<=sqrt(2)K_(m-1)+32,
    K12<=2016+992sqrt(2)<3504<4094=2^12-2.

The comparison at the next integer follows whenever 2^(m-1)>=31,
again automatic for m>=12. Therefore actual star holds at every
m>=12 for this entire rank-at-most-one-Y class. No assumption is
made on X,P,Q beyond the original rank-at-most-two full block.

The arguments here do not settle general star at m=3..16 or the
rank-at-most-one-Y cases at m=3..11. They also do not transfer these
large-m class statements back to smaller m, or to unrestricted
rank-two endpoint operators. In particular no all-n NPT conclusion
is inferred.

## Exact audit

The companion checker executes the new functional and concurrence
audits, checks all new rational comparisons and both recurrences,
rechecks the exact rank-one-Y factor gauge and cost cancellation,
and verifies the bridge on the established full-rank actual controls
with unequal complex planes and a rank-four whole error. It performs
no optimization or numerical search. The all-dimension proof is the
argument above, not an extrapolation of the controls.
