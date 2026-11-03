# A uniform merged bridge for all rank-two block data

The coefficient-400 proof below is retained. The later
[coherent-completion theorem](induction24_free_completion_uniform_bridge.md)
improves its coefficient to 2129/10 and proves actual star for every
m>=17. The current remaining general cases are m=3 through 16;
ORB and unrestricted endpoint induction are still open.

For every actual rank-at-most-two block operator C=[X,P;Q,Y] in the
one-free-copy setting, on C3_free tensor ((C2_new tensor H) direct-sum
an excitation line), the merged forms satisfy

    sqrt(q_1(Y))[sqrt(Gamma(X))-sqrt(Gamma_bar(X))]
       <=400 sqrt(qtilde(P)qtilde(Q)).                     (1)

H has arbitrary finite dimension. The bound permits either rank of Y,
different free planes, arbitrary scalar fields and simultaneous outside
components. Its coefficient is coarse. Combining it with the exact
energy-aware recursion of [star19](star19_energy_aware_recursion.md)
proves actual star for **every m>=19**. It does not prove the required
coefficient at m=3 through 18, ORB, or the general all-n endpoint.

The proof uses the all-field concurrence rectangle from
[star21](star21_concurrence_rectangle.md), the unscaled functional
estimate in [induction17](induction17_rank_one_Y_bridge.md), and the
actual two-copy positivity inputs. No numerical search is used.

## 1. Two inputs with their exact scope

Write arbitrary two-column physical fields as

    R=V tensor a+r, L=U tensor c+delta,
    U*U=V*V=I_2,
    Tr(V*r_z)=Tr(U*delta_z)=0,
    alpha=||r||, beta=||delta||, A=||a||, C0=||c||.

All norms here are Frobenius norms unless explicitly labeled operator
norms. The scalar maps have operator norms A and C0, respectively.
For the actual positive internal two-qubit matrix

    B=(R* tensor L*)K_H(R tensor L),
    K_H=2W_free(I-F_new)(I-F_H), W_free=I-F_free/2,

star21 proves, without a scalar-baseline assumption,

    C(B)<=alpha beta[152 A C0
             +4sqrt(2)(C0 alpha+A beta)+2alpha beta].     (2)

C(B) is homogeneous two-qubit concurrence. Its determinant bound is

    H(R M L*)>=-|det M| C(B),
    H(Z)=Gamma_bar(Z)-Gamma(Z).                          (3)

This follows by applying the pure determinant inequality to any ensemble
of B and taking its convex roof, as in
[orb9](orb9_internal_concurrence.md). No finite dual optimizer is assumed.

The second input is a functional inequality already proved before the
cost conversion in induction17. For **any** factors F=[f1,f2] and
L=U tensor c+delta, not requiring any Y or B-sector assumption,

    sqrt(Gamma(F L*))-sqrt(Gamma_bar(F L*))
       <=35||f1|| ||delta||.                            (4)

Indeed the proof there gives the maximum of 10(sqrt(2)+2) and
1044sqrt(2)/43, both below 35. Every compared matrix has range in
span{f1,f2}. Thus (4) is valid independently of the rank-one-Y gauge
used later in that note. Applying it to an adjoint is also allowed.

We retain the energy estimates

    0<=Gamma(Z)<=2||Z||^2, 0<=Gamma_bar(Z)<=4||Z||^2
                                                    (rank Z<=2),
    Gamma(uv*)>=||u||^2||v||^2/2.                       (5)

Their proofs, including partial-trace rank bounds and the qubit-frame
norm moment, are in star20 and induction17.

## 2. Full-rank Y and a uniform adjugate comparison

Y of rank at most one is already covered by coefficient 50. For rank-two
Y, homogeneity permits normalizing its largest singular value to one.
One may scale the entire C by that positive singular value's reciprocal;
both sides of (1) then have the same quadratic scaling. Write

    Y=U D V*, D=diag(1,t), 0<t<=1,
    M=D^(-1)=diag(1,s), s=1/t.

The exact rank-two block identities give

    P=R V*, Q=U L*, X=R M L*,
    R=V tensor a+r, L=U tensor c+delta,
    p=qtilde(P)=alpha^2, q=qtilde(Q)=beta^2,
    y=q_1(Y)<=1+t^2<=2.                                (6)

The fields R and L are unrestricted; their decompositions are orthogonal
free scalar projections. In particular ||R||_op<=A+alpha and
||L||_op<=C0+beta.

We need the comparison

    q_1(Y^+)>=s^2 y/2.                                 (7)

Here Y^+=V M U*. To prove it without aligning the singular vectors
with the principal planes, choose new internal bases for U,V in which
their overlap is K=diag(1,k), 0<=k<=1. In those bases Y=U A0 V* for
an arbitrary invertible 2-by-2 A0, with |det A0|=t. Put

    q_K(Z)=||Z||^2-|Tr(K Z)|^2/2.

For Z=[[a,b],[c,d]], direct expansion gives

    2q_K(Z)-q_K(adj Z)
       =|ka-d|^2/2+(1-k^2)|d|^2+|b|^2+|c|^2>=0.

Apply this also to adj Z, whose adjugate is Z. Thus
q_K(adj A0)>=q_K(A0)/2. Since
q_1(Y^+)=q_K(adj A0)/t^2 and q_1(Y)=q_K(A0), (7) follows.
This retains arbitrary principal-plane and singular-vector orientations.

If y=0, (1) is immediate. If alpha=0 or beta=0, one field is exactly
scalar, so its fixed sector range or co-range is a safe axis and H(X)>=0.
We may therefore assume y,alpha,beta>0.

## 3. Small scalar amplitudes require no condition-number estimate

If A sqrt(y)<=10alpha, apply (4) with F=R M. Its first column is
R1=v0 tensor a+r1, of norm at most A+alpha. Therefore

    sqrt(y)[sqrt(Gamma(X))-sqrt(Gamma_bar(X))]
       <=35sqrt(y)(A+alpha)beta
       <=35(10+sqrt(2))alpha beta<400alpha beta.         (8)

The second column of F can grow like 1/t; (4) imposes no bound on it.
This is what avoids a singular-value-ratio loss.

If C0 sqrt(y)<=10beta, apply (4) to

    X*=(L M)R*.

The right factor R has scalar plane V and residual norm alpha; the
first column of L M is L1, of norm at most C0+beta. The adjoint
invariance of Gamma and Gamma_bar gives (8) again. No inverse estimate
on t was used in either case.

## 4. Both large scalar amplitudes retain the actual energy denominator

It remains to consider

    A sqrt(y)>=10alpha, C0 sqrt(y)>=10beta.              (9)

Let R0=V tensor a, L0=U tensor c, X0=R0 M L0*. Its product structure
and (7) imply

    Gamma(X0)>=q_1(Y^+) A^2 C0^2>=s^2 y A^2 C0^2/2.

Use **two** energy triangles, never one on the potentially rank-four
difference X-X0. Set X1=R M L0*. X and X1 have a common range of
dimension at most two. X1 and X0 have a common co-range of dimension
at most two. Consequently (5) gives

    sqrt(Gamma(X))
      >=sqrt(Gamma(X0))-sqrt(2)||r M L0*||
                              -sqrt(2)||R M delta*||
      >=s[sqrt(y)A C0/sqrt(2)
                   -sqrt(2)(C0 alpha+A beta+alpha beta)]
      >=s sqrt(y)A C0[1/sqrt(2)-2sqrt(2)/10-2/100]
       >(2/5)s sqrt(y)A C0.                            (10)

The operator norms of R0,L0, rather than their Frobenius norms, justify
the middle line. Every triangle stays on a rank-at-most-two linear
space. In the last line we used y<=2, (9), and sqrt(2)>7/5.

If H(X)>=0 the bridge is immediate. Otherwise (3), |det M|=s, and
the denominator in (10) give

    sqrt(y)[sqrt(Gamma(X))-sqrt(Gamma_bar(X))]
       <=(5/2) C(B)/(A C0).

Insert (2) and use alpha/A,beta/C0<=sqrt(y)/10:

    (5/2) C(B)/(A C0)
      <=(5/2)alpha beta[152+4sqrt(2)(alpha/A+beta/C0)
                                      +2alpha beta/(A C0)]
      <=(5/2)(152+16/10+4/100)alpha beta
       =(3841/10)alpha beta<400alpha beta.              (11)

Together with (8), the zero cases and the rank-one-Y theorem, this
proves (1) for every actual rank-at-most-two C. No limit with varying
fields or ill-conditioned Y is omitted.

## 5. Consequence for actual star, and the remaining requirements

The energy-aware recurrence now holds without a field restriction:

    K_m<=sqrt(2)K_(m-1)+400, K_1=0,
    K_m<=400(1+sqrt(2))[2^((m-1)/2)-1].

At m=19,

    K_19<=204,400(1+sqrt(2))
         <511,000<524,286=2^19-2,

using sqrt(2)<3/2. If K_m<=2^m-2, then

    K_(m+1)<=(3/2)(2^m-2)+400
           =3*2^(m-1)+397<=2^(m+1)-2

whenever 2^(m-1)>=399, which holds for every m>=19. Actual star is
therefore proved for all m>=19, hence the one-free-copy, one-excitation
endpoint class for n=m+1>=20. This is a proper operator class.

This coefficient-400 proof does not establish the required star
coefficients at m=3 through 18; the later theorem linked above covers
m=17,18 as well. The constant-four merger candidate is still unproved; this
argument gives a finite constant rather than its sharp value. General
ORB and the unrestricted induction EP_(n-1)=>EP_n are separate open
requirements. No all-n NPT solution follows from this theorem.

`star22_uniform_merger_bridge.py` checks the adjugate identity on
unrestricted complex entries, canonical blocks with unequal free planes,
the rank-two triangle spaces (including rank-four whole errors), actual
merged forms and positive-gap controls, and every scalar coefficient
and all-larger-m comparison. The all-dimension proof is the displayed
argument; the exact finite controls are its implementation audit.
