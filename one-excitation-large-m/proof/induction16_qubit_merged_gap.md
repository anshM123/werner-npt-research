# A sharp merged gap on a qubit times an arbitrary dimension

This proves a sharp auxiliary inequality for rank-two operators on
C2 tensor H. It supplies MB_4 from
[star19](star19_energy_aware_recursion.md) when the full A-sector
operator X has a common rank-one free factor. It does **not** prove
MB_4 for arbitrary X, or even for all rank-one-Y block data. The
sharp general bridge MB_4 and all-copy endpoint remain open. A later
[uniform merger theorem](star22_uniform_merger_bridge.md) proves the
coarser general bridge MB_400 and actual star for every m>=19.

The only quantitative two-copy input needed below is proved explicitly
in Section 2. The endpoint statement EP_2 alone would not supply that
margin. The same balanced-factor sum of squares is already derived in
the first-order paper, in
"Sharp forms with a flatness term"; its two-factor statement is
separate from the subsequent n-factor conjectures.

## 1. Sharp auxiliary gap

For rank(Z)<=2 on C2 tensor H, put

    gamma(Z)=2||Z||_F^2-|Tr Z|^2,
    s(Z)=4q_2(Z)
        =4||Z||_F^2-2||Tr_1 Z||_F^2-2||Tr_2 Z||_F^2+|Tr Z|^2.

If sigma_1>=sigma_2>=0 are its singular values, then

    sqrt(gamma(Z))-sqrt(s(Z))<=2 sigma_2.                 (1)

Both square roots are real. Gamma>=0 follows from
|Tr Z|<=sigma_1+sigma_2; positivity of s also follows from Section 2.
The constant two is sharp, even on C2 tensor C2.

Define

    h(Z)=||Z||_F^2-||Tr_1 Z||_F^2-||Tr_2 Z||_F^2+|Tr Z|^2,
    D=(I-F_1)(I-F_2),

where the swaps act between two replicas of C2 tensor H. For unit
u,v, the swap trick gives

    h(uv*)=<u tensor v,D(u tensor v)>.

In fact

    0<=h(uv*)<=1.                                       (2)

The lower bound follows from D>=0. For the upper bound, condition D
on u in its first replica. The resulting operator on the second is

    K_u=I-rho_(u,1) tensor I-I tensor rho_(u,2)+|u><u|.

The Schmidt rank of u is at most two. In a Schmidt basis
u=sqrt(lambda_0)|00>+sqrt(lambda_1)|11>, lambda_0+lambda_1=1,
the only nonzero part within its 2-by-2 Schmidt support is the
rank-one projector onto
sqrt(lambda_1)|00>+sqrt(lambda_0)|11>. On H directions outside
that support, K_u has eigenvalues lambda_0 or lambda_1. The two
crossed Schmidt vectors have eigenvalue zero. Thus 0<=K_u<=I.
The Schmidt-rank-one case follows by setting one lambda to zero;
H of dimension one gives D=0 directly. This proves (2), in arbitrary
complex local bases, without a dimension-dependent constant.

Expand an SVD Z=sigma_1 u_1v_1*+sigma_2 u_2v_2*. The full physical
replica identity, including the crossed right indices, is

    h(Z)=sum_(i,j) sigma_i sigma_j
               <u_i tensor v_j,D(u_j tensor v_i)>.

Its diagonal terms are nonnegative. Cauchy--Schwarz for positive D,
followed by (2), bounds the modulus of the off-diagonal term by one.
Consequently

    h(Z)>=-2 sigma_1 sigma_2,
    s(Z)=gamma(Z)+2h(Z)>=gamma(Z)-4 sigma_1 sigma_2.       (3)

Now split into two exact cases.

* If |Tr Z|<=sigma_1-sigma_2, then
  gamma>=(sigma_1+sigma_2)^2. The function
  sqrt(x)-sqrt(x-4 sigma_1 sigma_2) decreases on this range.
  Equation (3) therefore gives a gap at most
  (sigma_1+sigma_2)-(sigma_1-sigma_2)=2 sigma_2.
* If |Tr Z|>=sigma_1-sigma_2, then
  gamma<=(sigma_1+sigma_2)^2. The margin proved next gives
  s>=(sigma_1-sigma_2)^2, again proving (1).

The zero and rank-one cases are included, with sigma_2=0.

## 2. Self-contained two-factor flat margin

For an operator Z of rank at most two on **any** finite tensor product
H_1 tensor H_2, one has

    q_2(Z)>= (sigma_1-sigma_2)^2/4.                       (4)

Here the two factors may have different dimensions. To prove this,
first choose a balanced factorization

    Z=x_1y_1*+x_2y_2*,
    <x_i,x_j>=<y_i,y_j>=g_ij,
    <y_i,x_i>=Tr Z/2, i=1,2.

For existence, start with X0=U sqrt(Sigma), Y0=V sqrt(Sigma)
from an SVD, padding a zero column if necessary. For the 2-by-2
matrix M=Y0*X0, write M=Tr(M)I/2+m dot sigma in the Pauli basis.
Choose a unit real Bloch vector perpendicular to both Re(m) and Im(m).
Its orthonormal eigenbasis gives a unitary Q with both diagonal
entries of Q*MQ equal to Tr(M)/2. The columns of X0Q,Y0Q have the
required properties and common Gram matrix G=Q*Sigma Q.

Let S(W)=||Tr_1 W||^2+||Tr_2 W||^2. The exact identity is

    ||Z||_F^2+||Z||_*^2/2+|Tr Z|^2/2-S(Z)
       =(g_11-g_22)^2/2
          +h(x_1y_1*)+h(x_2y_2*)+h(x_2x_1*)+h(y_2y_1*)
          +S(x_2x_1*-y_2y_1*).                          (5)

All h terms here are nonnegative in any two dimensions, since they
are expectations of D>=0; no qubit upper bound is used in (5).
To verify the identity, use

    S(ab*)=||a||^2||b||^2+|<a,b>|^2-h(ab*),

and the swap identity that transforms the mixed term in S(Z) into

    S(x_2x_1*)+S(y_2y_1*)-S(x_2x_1*-y_2y_1*).

Expanding and using the common Gram matrix and equal diagonal
overlaps gives (5). Hence its left side is nonnegative. Finally
||Z||_*^2=||Z||_F^2+2 sigma_1 sigma_2 shows that its left side
equals 2q_2(Z)-(sigma_1-sigma_2)^2/2. This proves (4).

Thus the margin in the second case of Section 1 is established by a
sum of squares, rather than inferred from qualitative EP_2.

## 3. MB_4 for a common rank-one free factor

Let C=[X,P;Q,Y] have rank at most two in the setting of star19, with
X acting on C3_free tensor C2_new tensor H. Suppose

    X=(a b*) tensor Z                                   (6)

for free vectors a,b and a rank-at-most-two Z on C2 tensor H.
No factorization of Y, P or Q is assumed. Then the actual bridge is

    sqrt(q_1(Y)) [sqrt(Gamma(X))-sqrt(Gamma_bar(X))]
       <=4 sqrt(qtilde(P)qtilde(Q)).                    (7)

Indeed write f=a b*. The definitions give

    Gamma(X)=q_1(f) gamma(Z),
    Gamma_bar(X)=q_1(f) s(Z).

Since q_1(f)<=||f||_F^2, (1) gives

    sqrt(Gamma(X))-sqrt(Gamma_bar(X))<=2 sigma_2(X).      (8)

If rank X<=1, the left gap is nonpositive by the positive rank-one
physical merger kernel, so (7) is immediate. Suppose rank X=2.
Then rank C=rank X, and the exact block-rank identities give

    ran P subset ran X,   ran Q* subset ran X*,
    Y=Q X^+ P.

In particular each free 3-by-3 block of P has range in the line
Ca, and each free block of Q has co-range in Cb. They all have rank
at most one. The rank-one trace inequality therefore gives

    p:=qtilde(P)>=||P||_F^2/2,
    q:=qtilde(Q)>=||Q||_F^2/2.                           (9)

Schatten Holder and the pseudoinverse bound also give

    ||Y||_*<=||Q||_F ||X^+P||_F
           <=||P||_F ||Q||_F/sigma_2(X).

Finally sqrt(q_1(Y))<=||Y||_* since rank Y<=2, so multiplying
(8) by sqrt(q_1(Y)) and applying (9) proves (7). This includes
rank-one Y, zero Y and singular cross blocks. It does not assert
the coercivity in (9) for arbitrary X: the common free factor (6)
is what forces every free cross block to have rank at most one.

The coefficient four is sharp in this class. The exact witness from
star19 has X=|0><1| tensor (|00><00|-|11><11|),
q_1(Y)=4, p=q=1, Gamma=4 and Gamma_bar=0, attaining (7).

This is a **merged two-sector** bound. It does not replace Gamma_bar
by the full m-qubit Pauli form Sigma_m. In fact the counterexample in
star_simple_counterexample.py already
has a rank-one common free factor and refutes that full-Pauli auxiliary
coefficient two at m=3. There is no conflict: merging the previous
sector changes the comparison form. Also, the free-factor class above
already has p/||P||^2,q/||Q||^2>=1/2, so the existing
trace-margin theorem (a working note not included here) covers its actual
star inequality for m>=4. No new all-m star threshold is asserted here.

As a direct corollary, (7) holds for **every product operator**
X=F_free tensor Z_sector of rank at most two. The remaining nonzero
case has rank F<=2 and rank Z=1. Here q_1(F)>=0 and h(Z)>=0, so
Gamma_bar(X)-Gamma(X)=2q_1(F)h(Z)>=0 and the bridge is immediate.
This extension imposes no additional conditions on the full blocks.

## Exact checks and remaining gap

`induction16_qubit_merged_gap.py` checks the physical crossed-index
formula for h on fixed rational-complex SVD data on C2 tensor C3,
both pure-state-inversion bounds, the complete balanced-factor SOS
(5), and the sharp auxiliary and block witnesses. No search or
floating-point eigenvalue calculation is used. These checks audit
the displayed dimension-independent proofs; finite samples are not
used to infer the inequalities.

For a general rank-one Y, one two-column B-sector factor has rank
one, but the other may have rank two. Then free blocks of one cross
field can have rank two and qtilde can vanish through its scalar
rank-two-projection direction. Thus (9), and the common-free-factor
reduction itself, need not hold. Controlling that case, and general
free-sector entanglement in X, remains the missing part of MB_4.
