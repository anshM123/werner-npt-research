# A universal flat margin for the actual merged three-factor form

For every rank-at-most-two operator X on C3_free tensor C2_new tensor H,
where H is any finite-dimensional Hilbert space, let

    Gamma(X)=2q_2(X),
    Gamma_bar(X)=4q_3(X).

Here the first q_2 merges new tensor H into its second factor, whereas
q_3 retains all three factors, with each local endpoint map Id-E/2.
These are exactly the physical merger forms in
[star19](star19_energy_aware_recursion.md). If s>=t>=0 are the two
singular values, padded by a zero when needed, then

    Gamma_bar(X)>=(s-t)^2/2.                         (1)

In particular an actual Gamma_bar-zero of rank two has equal singular
values. This excludes the unequal-singular-value null boundary present
in some abstract comparison counterexamples. The constant one half is
sharp. A separate consequence below is the uniform auxiliary bound
2sqrt(2) sigma2; its constant is not claimed optimal. Neither result
proves the desired coefficient-two auxiliary gap or MB_4.

## 1. Two inputs already proved in the repository

For rank-at-most-two Z on any two finite tensor factors, the explicit
balanced-factor sum of squares in
[induction16, Section 2](induction16_qubit_merged_gap.md) proves

    q_2(Z)>=(sigma1(Z)-sigma2(Z))^2/4.               (2)

This quantitative fact is stronger than positivity alone. Therefore
the previous merged form satisfies

    Gamma_rest(Z)=2q_2(Z)>=(sigma1-sigma2)^2/2.

The exact qubit product-frame identity from star19 is

    Gamma_bar(X)=sum_j w_j Gamma_rest(<u_j|X|v_j>), (3)
    2I-F_new=sum_j w_j |u_j><u_j| tensor |v_j><v_j|.

Every contracted operator has rank at most two, so all summands in (3)
are nonnegative. The two-copy theorem suffices; no general three-copy
flat-margin assertion is assumed.

## 2. Two qubit axes balance the opposite factor traces

Choose balanced SVD factors

    X=A B*, A=U sqrt(Sigma), B=V sqrt(Sigma),
    A*A=B*B=Sigma=diag(s,t).

The reduced new-qubit density matrices of AA* and BB* have the same
trace S=s+t. Write their real Bloch vectors as a,b in R3. Choose TWO
orthogonal real unit axes n1,n2 perpendicular to a+b. This is possible
also when a+b=0, in which case any orthonormal pair works. Complete them
to an orthonormal triple n1,n2,n3. For either of the first two axes let
u,v be its antipodal orthonormal qubit vectors. Then

    ||A_u||^2=||B_v||^2=m_plus,
    ||A_v||^2=||B_u||^2=m_minus,
    m_plus+m_minus=S,
    A_u=<u|A, B_u=<u|B.                            (4)

Indeed the difference of the first two cross traces is
(a+b) dot n/2. The complementary traces then agree because both whole
factors have trace S. The individual left/right traces need not all be
S/2; retaining the possibly different m_plus,m_minus permits two axes.

The six antipodal atoms belonging to any orthonormal triple of Bloch
axes are a positive frame in (3), each with weight one. Explicitly, if
u,v give its third axis, the other four atoms can be written

    x_j=(u+i^j v)/sqrt(2),
    y_j=(u-i^j v)/sqrt(2), j=0,1,2,3,

also each with weight one. Direct expansion gives

    2I-F=P_u tensor P_v+P_v tensor P_u
                          +sum_j P_(x_j) tensor P_(y_j). (5)

Choosing the harmless phases of u,v gives any desired orientation of
the other two axes. Thus the frame has total weight six. The two axes
n1,n2 provide FOUR cross-balanced atoms. For one antipodal pair, (2)
bounds its two energy terms below by half the sum of its squared
singular-value gaps; the third pair's energies remain nonnegative.

    Gamma_bar(X)>=sum_(first two pairs) delta_atom^2/2. (6)

The next identity lower bounds the gap sum for EACH balanced pair.

## 3. The equal-trace two-by-two Gram identity

For positive-semidefinite two-by-two K,L with Tr K=Tr L=m, the
singular-value gap of an operator with left/right Grams K,L obeys

    delta(K,L)^2=Tr(KL)-2sqrt(det K det L)
       =||K+L-mI||^2/2+(sqrt(det K)-sqrt(det L))^2. (7)

The identities remain valid for singular Grams. The two singular
values squared are the eigenvalues of the Gram product, so their
squared sum is Tr(KL) and product is sqrt(det K det L).
Expanding the Frobenius norm and using the two-by-two characteristic
identity proves the second equality. Hence

    delta(K,L)^2>=||K+L-mI||^2/2.                   (8)

Apply this to

    K=A_u*A_u, L=B_v*B_v, Tr K=Tr L=m_plus,
    K'=A_v*A_v=Sigma-K, L'=B_u*B_u=Sigma-L,
    Tr K'=Tr L'=m_minus.

Put D1=K+L-m_plus I and D2=K'+L'-m_minus I. By (4),

    D1+D2=2Sigma-SI,
    ||D1+D2||^2=2(s-t)^2.

The parallelogram identity now gives

    delta1^2+delta2^2
       >=(||D1||^2+||D2||^2)/2
       >=||D1+D2||^2/4=(s-t)^2/2.                 (9)

The two balanced pairs each satisfy (9). Combining them with (6)
proves (1).

More explicitly, its nonnegative remainder is the sum of the two
third-axis frame energies, the four remainders
Gamma_rest(X_atom)-delta_atom^2/2 from (2), and, for each balanced pair,

    ||D1-D2||^2/8
       +(sqrt(det K)-sqrt(det L))^2/2
       +(sqrt(det K')-sqrt(det L'))^2/2.

This retains all complex orientations. The zero case is immediate;
rank one is included by a zero singular value and padded factors.
If H has dimension one, the proof still applies, with the corresponding
two-factor statements and qubit contractions unchanged.

For sharpness take any fixed product sector vector and
X=(s|0><0|+t|1><1|)_free tensor P_new tensor P_H. Direct substitution
gives Gamma_bar(X)=(s-t)^2/2, so (1) cannot have a larger coefficient.

## 4. A universal auxiliary corollary and its limits

The actual weighted pure-kernel estimate in
[star20](star20_weighted_qubit_kernel.md) gives

    Gamma_bar(X)-Gamma(X)>=-4st.

For a positive gap, write b=Gamma_bar(X). Since
sqrt(b+4st)-sqrt(b) decreases with b, (1) gives

    sqrt(Gamma)-sqrt(Gamma_bar)
       <=sqrt((s-t)^2/2+4st)-(s-t)/sqrt(2)
       <=2sqrt(2)t.                               (10)

For t>0 the last inequality follows by squaring
sqrt(r^2+6r+1)<=r+3, r=s/t>=1, whose squared slack is eight.
The case t=0 is pure comparison. This improves the immediate
sqrt(2)+2 constant obtained by rank-one truncation and the separate
pure-energy norm bounds. Optimality of 2sqrt(2) is not claimed.

In particular (1) does not imply the desired coefficient-two auxiliary
gap. It also does not supply the separate P,Q cross-cost coercivity
needed to turn an auxiliary gap into MB_4. Universal MB_4, unrestricted
ORB, the general copy-induction step and the all-copy NPT endpoint
remain open.

## 5. Sharp flatness transfers through any additional qubit factors

The factor argument above is independent of the particular rest energy:
for every rank-at-most-two X on C2 tensor K there is an octahedral
frame, chosen from its balanced SVD factors, such that

    sum_(six atoms) [sigma1(X_uv)-sigma2(X_uv)]^2
                              >=[sigma1(X)-sigma2(X)]^2. (11)

Let F be any Hermitian quadratic form on operators on K satisfying
F(Z)>=alpha[sigma1(Z)-sigma2(Z)]^2 for all rank-at-most-two Z, with
alpha>=0. Tensor its superoperator with the actual qubit endpoint
Id-E/2. The same exact frame identity, with its scalar normalization,
is

    F_qubit(X)=sum_ab F_bil(X_ab,X_ab)-F(Tr_new X)/2
              =(1/2)sum_(six atoms) F(X_uv).

Since all X_uv have rank at most two, (11) proves

    F_qubit(X)>=alpha[sigma1(X)-sigma2(X)]^2/2.       (12)

The new qubit need not be a product factor of X. If F attains its
bound, tensoring its equality example with a rank-one qubit projector
attains (12), because the one-qubit form of that projector is 1/2.

Starting from q_1(Z)>= (s-t)^2/2, or the arbitrary-dimension two-factor
bound q_2(Z)>=(s-t)^2/4, and iterating (12) gives the exact all-n result

    q_n(X)>=2^(-n)[sigma1(X)-sigma2(X)]^2            (13)

whenever AT MOST TWO local factors have dimension greater than two.
Factors of dimension one only multiply the form by 1/2 and present
no difficulty. Tensor-factor order is immaterial. The coefficient is
sharp when at least one local factor has dimension at least two:
take a positive rank-two operator on that factor and rank-one
projectors on all remaining factors.

In particular the actual full Pauli form on one free qutrit and m
sector qubits has the normalization Sigma_m=2^m q_(m+1). Therefore

    Sigma_m(X)>=[sigma1(X)-sigma2(X)]^2/2            (14)

for every m and every rank-at-most-two X, with sharp coefficient.
This establishes the dimension-restricted instance of the flat-margin
candidate recorded in the first-order paper; it does not extend that
candidate to three or more qutrit factors. Positivity for these
qubit-containing classes was already available by tensoring separable
qubit factors. The conclusion here is the quantitative sharp margin.

Neither (13) nor (14) supplies the P,Q cross-cost estimate or settles
the all-qutrit copy-induction problem.

## 6. The nullspace on any fixed left plane is a maximally entangled code

Fix an isometry R:C2 -> physical space and identify X=R W with its
two-by-d coefficient matrix W, d=6 dim H. The compressed actual form
Gamma_bar is positive semidefinite on this whole coefficient space.
Let K be its complex nullspace. By (1), every W in K satisfies

    W W*=||W||^2 I2/2.

Complex polarization inside the linear space K therefore gives

    W Z*=Tr(W Z*) I2/2 for all W,Z in K.             (15)

For example expand the preceding identity for W+zZ for all complex z
and equate its z and conjugate(z) coefficients. Choose a basis W_j
of K with Frobenius inner products <W_i,W_j>=2delta_ij. Equation (15)
becomes W_i W_j*=delta_ij I2. Thus the two rows of all W_j form an
orthonormal family of 2dim(K) vectors in C^d. Consequently

    dim(K)<=floor(d/2).

After a single unitary on the right coefficient space, these basis
matrices have rows e_j^T and e_(k+j)^T, where k=dim(K). Equivalently
K consists of two mutually orthogonal isometric copies of one common
row-vector space. This classifies the compressed nullspace, without
classifying which physical left planes admit it.

The bound is attained: for the orb32 plane R=C2_new tensor phi_free,H,
the nullspace consists exactly of I_new tensor |phi><chi| for arbitrary
chi in C3_free tensor H. Its dimension is 3dim(H)=d/2.

In particular, for every unit input vector a, the row penalty
Pi_a=|a><a| tensor I_d has compression I_K/2. The exact auxiliary pencil
from induction32 (a working note not included here) consequently obeys

    P_K[delta B-delta(1-delta)G+4(1-delta)Pi_a]P_K
       =(1-delta)[2I_K-delta G_K]
       >=2(1-delta)^2 I_K, 0<delta<1.              (16)

Here G and B are the fixed-support Gamma/Gamma_bar matrices. The last
step uses the established actual rank-two bound Gamma(X)<=2||X||^2,
so G_K<=2I_K. Equation (16) does not control the pencil cross block
between K and its orthogonal complement. It proves no full-pencil
positivity or sharp auxiliary bound.

## Exact audit

The companion checker verifies the six-atom frame as an exact complex
matrix identity, the general equal-cross-trace Gram and complementary
parallelogram identities with m_plus possibly different from m_minus,
actual full physical energy aggregation on a complex, nonorthogonal
internal-Gram control, and its balanced antipodal trace pairs. It
checks zero and rank-one controls, the symbolic sharp family for all
s,t, the scalar bound (10), and the qubit-transfer/Pauli normalizations.
It also checks a maximal-dimensional actual nullspace and its complex
row-penalty compression. It uses exact rational/algebraic arithmetic
and no optimization.
