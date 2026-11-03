# A coherent concurrence triangle and its actual fourth-component bound

For arbitrary matrices X,Y with four internal columns and a common
purification environment, write G(X)=X*X and use homogeneous two-qubit
concurrence. Then

    C(G(X+Y))<=C(G(X))+C(G(Y))+2||X||_F||Y||_F.       (1)

No separability hypothesis is needed. The sharp coherent-two theorem
in [orb24](orb24_coherent_separable_sum.md) is the special case where
both input concurrences vanish. No priority claim is made for (1).

For the actual physical merger kernel, a second consequence is

    C((E* tensor F*)K(E tensor F))
       <=4sqrt(det(E*E)det(F*F))<=||E||_F^2||F||_F^2. (2)

This keeps the exact two-column wedge geometry. Together, (1)-(2)
replace a fourth-component trace error by its actual concurrence and
retain an explicit witness-dependent nonnegative square. They do not
prove the sharp free rectangle, MB_4, ORB, or unrestricted endpoint
positivity or induction.

## 1. The shifted root-fidelity inequality

Let S,T be arbitrary positive two-qubit matrices and let
W=(|vec A><vec A|)^Gamma with det A=1. Put

    c_S=C(S), c_T=C(T), a=Tr(WS), b=Tr(WT).

The following strengthening of orb24 holds:

    ||sqrt(S) W sqrt(T)||_*
      <=sqrt(Tr S Tr T)+sqrt((a+c_S)(b+c_T)).          (3)

Both shifted expectations are nonnegative by the concurrence dual.
All singular and zero matrices are included.

Normalize both traces to one; homogeneity restores the general case.
Apply the same local Schmidt reduction, commuting four-phase twirl,
negative-modulus middle-coherence replacement, and common SWAP average
as in orb24. The root fidelity increases throughout. The original
witness expectations a,b can only decrease in the negative-modulus
step. Keep c_S,c_T fixed at their **original** values during these
reductions.

The only step requiring additional verification is the concurrence
control. The commuting local twirl cannot increase concurrence, by
convexity and local unitary invariance. For each twirled state, its
single middle coherence zeta can be made negative real by a local
diagonal unitary, so this replacement preserves that state's
concurrence. Let x,y be its outer populations and z=|zeta|. The
concurrence dual on the witnesses W_r gives

    c_S>=max(0,sup_(r>0) [2z-rx-y/r])
        =2max(0,z-sqrt(xy)).                          (4)

This also covers x=0 or y=0 by taking an infimum over positive r;
no finite minimizing r is required. SWAP averaging preserves x,y,z,
so the same lower bound remains available after balancing the middle
populations. An exact X-state concurrence formula is unnecessary.

The final scalar quantities are those of orb24:

    c=r sqrt(xX)+r^(-1)sqrt(yY)+sqrt(uU)+sqrt(vV),
    a_0=rx+r^(-1)y+u-v,
    b_0=rX+r^(-1)Y+U-V,
    x+y+u+v=X+Y+U+V=1,
    v-u=2z, V-U=2w.

Put A0=sqrt(rx)-sqrt(y/r), B0=sqrt(rX)-sqrt(Y/r).
Equation (4) now gives

    a_0+c_S-A0^2=2sqrt(xy)-2z+c_S>=0,
    b_0+c_T-B0^2>=0.                                 (5)

The crossed-coordinate cancellation and ordinary Cauchy--Schwarz give

    c-A0B0=sqrt(xY)+sqrt(yX)+sqrt(uU)+sqrt(vV)<=1.

Hence c<=1+sqrt((a_0+c_S)(b_0+c_T)). Since the original fidelity
is no larger, and a>=a_0,b>=b_0 with all shifted factors nonnegative,
this proves (3). This argument does not assume S or T is PPT.

## 2. The triangle and its square remainder

Let S=G(X), T=G(Y), z_X=||X||, z_Y=||Y||. Arbitrary polar
isometries of the two purifications imply

    |Tr(W X*Y)|<=||sqrt(S) W sqrt(T)||_*.

Thus for every determinant-normalized dual witness W,

    -Tr(W G(X+Y))
      <=c_S+c_T+2z_Xz_Y
        -(sqrt(a+c_S)-sqrt(b+c_T))^2.                 (6)

Taking the positive part of the dual supremum proves (1). The
supremum need not be attained. The coefficient two is necessary
already for product rows X=<00|,Y=<11|. The bound with both input
concurrences nonzero is also attained: X=Y=<Phi+| gives 4=1+1+2.

Any known upper bounds eta_S>=c_S, eta_T>=c_T may replace c_S,c_T
**inside both shifted terms** of (6). This follows directly from
(3), whose right side is monotone in the two nonnegative shifts.
In particular the valid witness-sensitive version is

    -Tr(W G(X+Y))
      <=eta_S+eta_T+2z_Xz_Y
        -(sqrt(a+eta_S)-sqrt(b+eta_T))^2.              (7)

The square in (6) or (7) is allowed to vanish. There is no claimed
uniform positive lower bound on it.

## 3. The actual double-residual concurrence costs at most pq

Use the actual kernel of [star20](star20_weighted_qubit_kernel.md):

    K=2W_free(I-F_new)(I-F_H), W_free=I-F_free/2,
    H(Z)=Gamma_bar(Z)-Gamma(Z).

For arbitrary finite H and arbitrary physical two-column factors E,F,
put

    B3=(E* tensor F*)K(E tensor F),
    p=||E||^2, q=||F||^2,
    eta_3=4sqrt(det(E*E)det(F*F)).

These may in particular be the trace-orthogonal residual factors in
the full physical rectangle; the proof needs no residual condition.
The exact physical crossed-index identity gives for every det-one A

    Tr(W_A B3)=H(E A F*).

Star20 proves H(Z)>=-4sigma_1(Z)sigma_2(Z) for every rank-at-most-two
Z, retaining the actual free and qubit-source kernel. The nonzero
singular values of EAF* satisfy the exact determinant identity

    sigma_1(EAF*)sigma_2(EAF*)
       =|det A|sqrt(det(E*E)det(F*F)).                 (8)

Consequently -Tr(W_A B3)<=eta_3 for every dual witness, proving
C(B3)<=eta_3. The two-by-two arithmetic-geometric mean inequality
det(E*E)<=p^2/4, and its F analogue, yield eta_3<=pq. Singular
factors have eta_3=0; in that case B3 is separable. There is no
inverse, condition-number bound, or limiting witness assumption.

The last coefficient is sharp for unrestricted two-column E,F.
For H=C3 take omega=(|00>+|11>+|22>)/sqrt(3) on free tensor H,
u_i=|i>_new tensor omega, and E=F=[u_0,u_1]. The exact compression is

    B3=2(I_4-F_internal)=4|Psi-><Psi-|,
    p=q=2, C(B3)=4=eta_3=pq.

This sharp control is not asserted to satisfy the extra trace-
orthogonality conditions relative to any specified scalar free plane.
The inequality itself applies to all such residual fields.

## 4. What this changes for the full physical rectangle

Write B012=G(Z0+Z1+Z2), B3=G(Z3), B=G(Z0+Z1+Z2+Z3).
Equation (1), rather than the older trace-only perturbation estimate,
gives

    C(B)<=C(B012)+eta_3+2||Z0+Z1+Z2|| ||Z3||.          (9)

The more informative version is (7), applied before discarding its
witness-dependent square. If a=Tr(W B012), b=Tr(W B3), then for
any eta_012>=C(B012),

    -Tr(W B)<=eta_012+eta_3+2||Z012||||Z3||
                 -(sqrt(a+eta_012)-sqrt(b+eta_3))^2. (10)

For example, the proved eta_012=13CDsqrt(pq) from star25 and the
unchanged component trace bounds imply

    C(B)<=sqrt(pq)[21CD
                   +4sqrt(2)(Dsqrt(p)+Csqrt(q))
                   +sqrt(pq)].                      (11)

Here C,D are scalar-sector norms. The final term is pq, rather than
the previous 2pq. The determinant version eta_3 in (9)-(10) can be
strictly smaller still. More importantly, (10) keeps a compensation
that a scalar concurrence upper bound discards. No relation strong
enough to control that square by the required actual merger energies
is proved here. Thus neither (10) nor the improved scalar bound (11)
settles MB_4 or unrestricted induction.

The exact square can vanish on actual nonzero physical data, even
with the exact input concurrences. Take H=C2, U=V the standard free
plane, N=|0><1|, and sector vectors

    c=d=|00>, e=|11>+|00>/2,
    E_z=e_z VN, F_z=e_z VN*.

All residual trace conditions hold. The single real channel has
`c^T Omega d=0`, `e^T Omega d=c^T Omega e=1`, and
`e^T Omega e=1`. Therefore the actual component Grams are

    B012=2[|00+11><00+11|+|10><10|], C(B012)=4,
    B3=2|10><10|,                        C(B3)=0.

For A_t=[[0,t],[-1/t,0]], t>0, det A_t=1, the two witness
expectations are `a=-4+2/t^2`, `b=2/t^2`. Thus
`a+C(B012)=b+C(B3)` and the square in (6) vanishes for every t.
This rules out a positive remainder bound valid for all actual
fields. These data are not asserted to be a negative critical point;
they do not exclude a stronger bound using that additional structure.

The companion exact checker verifies the new shifted scalar
cancellation, complex phase controls and X-block characteristic
identity, both sharp coherent controls, the two-column determinant
identity, the sharp actual K-compression, and the actual vanishing-
square control above. It performs no
search. The proofs remain valid for arbitrary finite H; the fixed
controls only verify the algebra and conventions used in them.
