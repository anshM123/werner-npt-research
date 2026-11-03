# Exact reduction of the full first-order rectangle to a free-only lemma

This reduces the universal sharp first-order concurrence target to a
fixed qutrit-only problem. It does not prove the target, MB_4, the
remaining small-m actual-star inequalities, ORB, or the unrestricted
endpoint. The established universal merger coefficient is now 2129/10
in [induction24](induction24_free_completion_uniform_bridge.md).
The weaker free coefficient 13/2, physical first-three coefficient 13,
and full physical scalar-baseline coefficient 21 are proved in
[star25](star25_free_completion_rectangle.md); the sharp target below
remains open.

## The fixed free-only target

Let V,U:C2 to C3 be isometries, let E,F be arbitrary complex 3-by-2
matrices with Tr(V*E)=Tr(U*F)=0, and let alpha be any complex number.
On the two free qutrits put

    W=I-F_free/2,
    Z=alpha V tensor U+E tensor U+V tensor F.

The proposed finite-dimensional inequality is

    C(Z*WZ)<=2||E||_F||F||_F.                         (FR_2) [OPEN]

There is no bound on alpha, no equality of the free planes, and no
support restriction on the outside components of E,F. As throughout,
C denotes homogeneous two-qubit concurrence, with pure cost 2|det|.
The coefficient two is necessary: for U=V the standard plane,
E=N=|0><1|, F=N*, alpha=0, the sharp matrix in star23 divided by
two has concurrence two and both residual norms are one. The
coincident-plane proof (a working note not included here) establishes
(FR_2) when the two planes coincide.

## 1. The full physical target

For arbitrary finite H, let

    D=(I-F_new)(I-F_H), K=2W tensor D,
    R_a=c_a V+E_a, L_b=d_b U+F_b,
    Tr(V*E_a)=Tr(U*F_b)=0,
    p=sum_a||E_a||_F^2, q=sum_b||F_b||_F^2,
    C_0=||c||, D_0=||d||.

The first three components of the actual purification are

    Z_012=K^(1/2)[(cV) tensor (dU)
                         +E tensor (dU)+(cV) tensor F],
    B_012=Z_012*Z_012.

The sharp full-field conjecture is

    C(B_012)<=4C_0 D_0 sqrt(pq).                      (P_4) [OPEN]

The common-plane theorem (a working note not included here) proves
this when U,V have the same free plane. The precise equivalence is:

    (FR_2) for all free-only data
       <=> (P_4) for all physical data and arbitrary H
       <=> (P_4) for all physical data with H=C2.

## 2. Real channel conventions and the forward reduction

Choose an orthonormal basis of H, and set

    J_new=[[0,1],[-1,0]],
    J_ij=|i><j|-|j><i| for i<j,
    Omega_ij=J_new tensor J_ij.

These coefficient matrices on C2_new tensor H are real and symmetric.
With vectorization in the order (new index,H index),

    D=sum_(i<j)|vec Omega_ij><vec Omega_ij|.           (1)

No normalization factors are suppressed: each J has squared Frobenius
norm two, so each Omega has squared Frobenius norm four. Fix one
channel Omega. Projecting the physical sector row against vec Omega
produces the free 9-by-4 operator

    Z_Omega=alpha_Omega V tensor U
                         +E_Omega tensor U+V tensor F_Omega,
    alpha_Omega=c^T Omega d,
    E_Omega=sum_a (Omega d)_a E_a,
    F_Omega=sum_b (Omega^T c)_b F_b.                   (2)

These are transposes, not adjoints, in the scalar bilinear contraction.
Indeed the projected physical row is
sum_(a,b)conjugate(Omega_ab) R_a tensor L_b, and Omega is real.
The matrices and scalar vectors c,d,E_a,F_b may all be complex.
Trace projection is preserved by the linear combinations in (2).
The physical first-order Gram is exactly

    B_012=2sum_Omega Z_Omega* W Z_Omega.                (3)

For any sector vector d, the channel moment matrix satisfies

    M_d=sum_Omega (Omega d)(Omega d)*
       =D_(conjugate(d))<=||d||^2 I.                  (4)

Here D_x is the contraction of the positive source D on its second
replica x. As proved in star23,

    D_x=||x||^2 I-rho_new(x)-rho_H(x)+|x><x|,
    0<=D_x<=||x||^2 I.

The conjugation in (4) follows directly from the projected row
convention and must be retained for complex fields. The upper bound
is valid for every x, so that conjugation changes no norm estimate.
If G_E[a,a']=Tr(E_a*E_a'), then G_E>=0, and (4) gives

    sum_Omega||E_Omega||_F^2=Tr(G_E M_d)<=D_0^2 p,
    sum_Omega||F_Omega||_F^2<=C_0^2 q.                 (5)

Assuming (FR_2), concurrence homogeneity, subadditivity and scalar
Cauchy in (3),(5) yield

    C(B_012)<=4sum_Omega||E_Omega||||F_Omega||
             <=4C_0D_0sqrt(pq),

which is (P_4). All zero-cost and zero-scalar cases are included; no
division is used. When dim H<2 the source D is zero.

## 3. Exact converse embedding, including zero alpha

Take H=C2. There is a single real channel
Omega=J_new tensor J_H; in particular Omega_11,00=Omega_00,11=1.
Fix any data V,U,E,F,alpha in the free-only problem.

For alpha!=0 choose

    c=alpha |11>, d=|00>,
    E_11=E, F_00=F/alpha,
    every other residual block zero.

Then (2) gives exactly alpha_Omega=alpha, E_Omega=E,
F_Omega=F. Furthermore

    C_0=|alpha|, D_0=1,
    p=||E||^2, q=||F||^2/|alpha|^2,
    C_0D_0sqrt(pq)=||E||||F||.

For alpha=0, instead choose

    c=d=|00>, E_11=E, F_11=F,
    every other residual block zero.

Again (2) gives alpha_Omega=0, E_Omega=E, F_Omega=F, now with
C_0=D_0=1 and p=||E||^2,q=||F||^2. Both embeddings preserve the
trace-free residual conditions and every outside component.

Thus in either case B_012=2Z*WZ and the right side of (P_4) is
4||E||||F||. Dividing the homogeneous concurrence inequality by
two proves (FR_2). This establishes the claimed equivalence, not the
inequalities themselves. A violation of (FR_2) would therefore be an
actual counterexample to (P_4), while remaining separate from actual
star because the full double-residual component and energy denominator
have not been included.

## 4. Eliminate the scalar baseline in the concurrence dual

There is also an exact way to remove alpha from the remaining free
target. Put

    Q_0=V tensor U, Q_1=E tensor U+V tensor F,
    A_0=Q_0*WQ_0, C_01=Q_0*WQ_1, A_1=Q_1*WQ_1.

For a complex 2-by-2 matrix T with det T=1, vectorized in the order
(00,01,10,11), define

    a(T)=<vec T|A_0^Gamma|vec T>,
    b(T)=<vec T|C_01^Gamma|vec T>,
    c(T)=<vec T|A_1^Gamma|vec T>.

The a,c values are real. The determinant-normalized concurrence dual
of star18 (a working note not included here) shows that (FR_2) for all
alpha is equivalent to the following scalar Hermitian pencil:

    [[a(T), b(T)],
     [conjugate(b(T)), c(T)+2||E||||F||]] >=0
                    for every det T=1.              (6) [OPEN]

Indeed the dual expectation for alpha Q_0+Q_1 is exactly

    a|alpha|^2+2Re(conjugate(alpha)b)+c.

When a>0, its minimum over alpha is c-|b|^2/a, at alpha=-b/a.
If a=0, b=0: separate scalar-axis PPT for alpha V tensor U+E tensor U
and alpha V tensor U+V tensor F forces both baseline cross terms to
vanish. Thus the zero case of (6) requires c+2||E||||F||>=0 and does
not require a division by a. No optimizer over noncompact SL(2,C)
is assumed to exist.

The baseline part is explicit. With K=V*U,

    A_0^Gamma=I-|vec K><vec K|/2,
    a(T)=||T||^2-|Tr(K* T)|^2/2
          >=||T||^2/2-1,                              (7)

because K is a contraction, |Tr(K*T)|<=||T||_*, and
||T||_*^2=||T||_F^2+2 when det T=1. Thus (6) retains the positive
baseline contribution that an isolated mixed-term bound would lose.
One may homogenize (6) to all T by replacing its bottom-right
tolerance with 2||E||||F|| |det T|. The singular-T endpoint follows
by continuity; alternatively its dual vector is product, so the
original positive free Gram already gives nonnegative expectations
for every alpha there. The resulting determinant condition is an
explicit quartic inequality in the four complex entries of T, with
the determinant modulus retained. It remains unproved for general
unequal-plane data.

## 5. Remove the residual fields: a fixed eleven-dimensional pencil

One can also quantify out E,F exactly. Fix orthonormal Frobenius bases
E_1,...,E_5 and F_1,...,F_5 of their trace-free 3-by-2 spaces relative
to V,U. The spaces each have complex dimension five. Define eleven
free operators

    Q_0=V tensor U,
    Q_i=E_i tensor U for 1<=i<=5,
    Q_(5+j)=V tensor F_j for 1<=j<=5.

For det T=1, form the Hermitian eleven-by-eleven matrix

    H_ij(T)=<vec T|(Q_i*WQ_j)^Gamma|vec T>.

Then (FR_2) is exactly equivalent to

    H(T)+diag(0,lambda I_5,lambda^(-1)I_5)>=0
      for every lambda>0, det T=1, and free planes.    (8) [OPEN]

To verify the quantifiers, write E=sum_i x_i E_i,F=sum_j y_j F_j.
The concurrence-dual expectation is the quadratic form of H on
(alpha,x,y). Adding the displayed diagonal adds
lambda||E||^2+lambda^(-1)||F||^2. Its infimum over lambda>0 is
2||E||||F||, including limiting endpoints when one norm vanishes.
Therefore requiring this quadratic form nonnegative for every
alpha,x,y and every lambda is equivalent to the original homogeneous
concurrence bound for every E,F,alpha. It is neither a relaxation nor
an exchange of a maximum with a minimum.

The two six-by-six principal blocks indexed by {0,1,...,5} and
{0,6,...,10} are positive by the two actual scalar-axis arguments.
They do not imply positivity of the full pencil (8). Principal-angle
coordinates and independent internal unitary basis choices reduce
the remaining free-plane parameter to k in [0,1]. Thus the entire
open target can be posed on this fixed eleven-dimensional pencil
depending only on k,T,lambda, rather than on arbitrary residual
fields or a sector dimension. Singular-T endpoints may again be
included by multiplying the cost diagonal by |det T|.

## 6. All-orientation coefficient structure and exact axis kernels

For further analytic work one may diagonalize T by independent
internal unitary changes, but must then retain arbitrary relative
orientations of U,V. Set T=diag(t,t^(-1)), t>0, and write

    A=V*E, B=U*F, Tr A=Tr B=0,
    k=Tr(TV*U), e=Tr(TU*E), r=Tr(TV*F), m=Tr(TE*F).

No simultaneous principal-plane diagonalization is imposed here.
The exact dual expectation splits as h_W=h_I-h_swap/2, where

    h_I=|alpha|^2||T||^2
          +2Re[conjugate(alpha)(Tr(T^2 A)+Tr(T^2 B))]
          +||ET||^2+||FT||^2+2ReTr(T A* T B),

    h_swap=|alpha|^2|k|^2+|e|^2+|r|^2
          +2Re[conjugate(alpha)(k e+conjugate(k)r)]
          +2Re[conjugate(k)m].                        (9)

These follow by expanding the two product-operator terms I and F_free.
For example, the swapped coefficient between A_1 tensor B_1 and
A_2 tensor B_2 has partial transpose
|vec(A_1*B_2)><vec(A_2*B_1)|. Formula (9) therefore retains the
complex conjugations needed for arbitrary free-plane orientations.

The two-by-two Cayley-Hamilton identity, with det T=1, gives

    T A* T=Tr(T A*)T+A*,
    Tr(T A* T B)=Tr(A*B)+conjugate(Tr(TA))Tr(TB),
    Tr(T^2 A)=Tr(T)Tr(TA).                            (10)

The same identities hold for B. They isolate several unbounded-t
coefficients into low-rank trace terms, but do not by themselves
prove the full pencil positive.

Each scalar-axis block has an additional exact rank property. For a
variable 3-by-2 matrix R,

    <vec T|[(R* tensor U*)W(R tensor U)]^Gamma|vec T>
        =||RT||^2-|Tr(U*RT)|^2/2.

This positive quadratic form on the six-dimensional R space has
rank five and kernel R=U T^(-1). In the coordinates R=alpha V+E,
with Tr(V*E)=0, its null vector has residual part

    E_0=U T^(-1)-[Tr(V*UT^(-1))/2]V.

Similarly the opposite scalar axis has kernel L=V T^(-1), with
trace-free residual

    F_0=V T^(-1)-[Tr(U*VT^(-1))/2]U.

For a(T)>0, eliminating the common alpha leaves two positive
five-by-five axis Schur complements, each of rank four, with these
explicit null residuals. Their mixed coupling is not controlled by
their separate positivity. Establishing positivity of the whole
cost-compensated pencil (8), including that coupling and arbitrary
orientation, remains the missing argument.

The companion exact checker verifies the real channel decomposition,
the complex moment convention, and both converse embeddings against
the actual physical purification. It also checks the explicit baseline
partial transpose and the complex alpha minimization identity. It
verifies the eleven-dimensional pencil evaluation and its exact cost
diagonal on complex rational fields, as well as (9) in a nonprincipal
complex orientation and both axis-kernel identities. It performs no optimization and
does not infer (FR_2) or (8) from finite tests.
