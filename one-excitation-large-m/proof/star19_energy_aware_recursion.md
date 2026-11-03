# Energy-aware recursion and the merged-sector bridge

This is an exact conditional reduction for **all** rank-at-most-two block
data, including rank-one Y and simultaneous outside components.
The subsequent [uniform merged-bridge theorem](induction24_free_completion_uniform_bridge.md)
proves a finite bridge coefficient 2129/10 and, through this recursion,
actual star for every m>=17. The coefficient-four bridge and general
actual star at m=3 through 16 remain open. The conditional discussion
below records the sharper target and its exact obstructions; coefficient
two for the bridge is false. ORB and the unrestricted all-copy endpoint
and induction remain separate open requirements.

## Actual forms and the merged bridge

Write the sector Hilbert space as C^2 tensor H, where H has arbitrary
finite dimension D. For X on C^3 tensor C^2 tensor H put

    T(Z)=sum_(a,b) q_1(Z_ab)
        =||Z||_F^2-(1/2)||Tr_free Z||_F^2,
    Gamma(X)=2T(X)-q_1(Tr_(C^2 tensor H) X),
    Gamma_bar(X)=4T(X)-2T(Tr_(C^2) X)-2T(Tr_H X)
                   +q_1(Tr_(C^2 tensor H) X).                 (1)

The notation T uses whatever sector coordinates remain after a partial
trace. In replica language these are the actual contractions with

    W_free tensor (2I-F_new F_H),
    W_free tensor (2I-F_new)(2I-F_H),
    W_free=I-F_free/2,

respectively. Hence

    Gamma_bar-Gamma
      =contraction[2W_free(I-F_new)(I-F_H)].                 (2)

For rank-two X, Gamma is nonnegative by the two-copy endpoint theorem.
Gamma_bar is four times the three-copy endpoint form on dimensions
3,2,D, but its nonnegativity also needs only the two-copy theorem:
insert a positive product frame for the separable qubit operator
2I-F_new, as in (4). Every contracted X_(u,v) has rank at most two,
and its previous merged Gamma is nonnegative by the two-copy theorem
on dimensions 3,D. This proves Gamma_bar>=0. The general three-copy
theorem [T1] is not needed. Positivity of the replica operator in (2)
does **not** make its contraction with a rank-two X nonnegative.

The needed dimension-independent bridge, with an unspecified constant B,
is

    sqrt(q_1(Y)) [sqrt(Gamma(X))-sqrt(Gamma_bar(X))]
      <=B sqrt(qtilde(P) qtilde(Q))                         (MB_B)

for every actual rank-at-most-two C=[X,P;Q,Y]. Y may have rank one
or two. No inverse or Schur identity is required in this formulation.

## Exact all-field aggregation

Factor C=[R_A;R_B][L_A;L_B]* with two-column factors. For unit vectors
u,v in the new C^2, contract that index in R_A and L_A only, leaving
R_B,L_B unchanged. The resulting C_(u,v) still has rank at most two,
has the same Y, and has

    X_(u,v)=<u|X|v>.

Let p_u=qtilde(P_u), q_v=qtilde(Q_v). There are positive semidefinite
2-by-2 matrices P_cost,Q_cost such that

    p_u=u* P_cost u, q_v=v* Q_cost v,
    Tr P_cost=p, Tr Q_cost=q,
    p=qtilde(P), q=qtilde(Q).

These Gram matrices remain positive even though q_1 is indefinite on all
3-by-3 matrices: every linear combination of the relevant blocks has
rank at most two, where q_1 is nonnegative. Thus this step includes
arbitrary scalar parts, noncommuting fields, outside components, and
rank-one Y.

Take any normalized positive product frame

    2I-F_new=sum_j w_j |u_j><u_j| tensor |v_j><v_j|.

Its total weight is six. The exact cost moment is

    sum_j w_j p_(u_j) q_(v_j)
      =2pq-Tr(P_cost Q_cost)<=2pq.                         (3)

The actual energy identities are

    sum_j w_j Gamma_(m-1)(X_j)=Gamma_bar(X),
    sum_j w_j Sigma_(m-1)(X_j)=Sigma_m(X).                 (4)

Here H=(C^2)^(m-1) in (1), and Sigma is the full nonidentity Pauli
form on the indicated qubits. The previous-copy Gamma is the two-copy
merged form, not the full Pauli form. Identity (4) follows directly by
inserting the product frame into the commuting replica operators.

Suppose the previous-copy star inequality holds uniformly with some
coefficient K_(m-1), not necessarily its desired coefficient. For each
frame atom it gives

    sqrt(q_1(Y) Gamma_j)
      <=sqrt(q_1(Y) Sigma_j)+K_(m-1) sqrt(p_j q_j).

The ordinary Euclidean triangle inequality, with frame weights absorbed
into the vectors, and (3)-(4) prove

    sqrt(q_1(Y)) [sqrt(Gamma_bar)-sqrt(Sigma_m)]
      <=K_(m-1) sqrt[2pq-Tr(P_cost Q_cost)]
      <=sqrt(2) K_(m-1) sqrt(pq).                        (5)

This also handles q_1(Y)=0 or pq=0 without dividing by them. The factor
sqrt(2) in the last **cost-moment bound** is sharp: orthogonally supported
rank-one P_cost,Q_cost have Tr(P_cost Q_cost)=0. This does not assert
sharpness of the complete energy-aware inequality (5).

Combining (MB_B) with (5) would give the genuine all-field recurrence

    K_m<=sqrt(2) K_(m-1)+B,
    K_1=0,
    K_m<=B [2^((m-1)/2)-1]/[sqrt(2)-1].                  (6)

Thus a universal bridge with B=4 would imply actual star for every
m>=5: (6) gives K_5=12(1+sqrt(2))<30=2^5-2, and the comparison
continues to hold at every larger integer m. It would still leave the
required m=3,4 inequalities unresolved. More generally, **any** universal
finite B would imply actual star for all sufficiently large m. Such a
constant is now proved separately in star22; this note itself supplies
the recursion and the lower-bound/obstruction calculations.

For the continuation comparison, if K_m<=2^m-2 at m>=5, then
sqrt(2)K_m+4<=2(2^m-2)+2=2^(m+1)-2. Thus the single exact
comparison at m=5 is sufficient for the stated conditional conclusion.

## The bridge must have coefficient at least four

Let c=|00>, w=|11>, and let e denote the excitation sector. In two
free-copy basis vectors |0>,|1>, take the rank-two factors

    R_1=|0>c+|1>e,     L_1=|1>c+|0>e,
    R_2=-|0>w+|1>e,    L_2=|1>w+|0>e.

Then

    X=|0><1| tensor (|c><c|-|w><w|),
    Y=2|1><0|,
    q_1(Y)=4, p=q=1, Gamma=4, Gamma_bar=0.

The bridge ratio is exactly four. This is the repository's explicit
constant-four witness, now evaluated as a merged bridge. It directly
refutes (MB_2).

This lower bound also applies when the bridge is restricted to rank-two
Y. Replace the second B-sector columns by

    R_(B,2)=|1>+epsilon|2>,
    L_(B,2)=|0>+epsilon|2>, epsilon>0,

and leave X unchanged. The full C retains rank at most two; Y now has
rank two, and

    q_1(Y)=4+2epsilon^2+epsilon^4/2,
    p=q=1+epsilon^2,
    bridge_ratio^2=4(4+2epsilon^2+epsilon^4/2)/(1+epsilon^2)^2 ->16.

So any rank-two-Y universal bridge constant is at least four as well.
This limiting example does not refute coefficient four.

## The arbitrary previous dimension reduces to at most 24

For the bridge alone, expand each of the two columns of R_A,L_A across
the free dimension three and the new qubit. Their H-coefficients are
at most 24 vectors total. A single orthonormal compression onto their
joint span preserves X, its rank, every norm and partial trace in (1),
the costs p,q, and Y. R_B,L_B are untouched. Consequently a universal
bridge is equivalent to its restriction to previous dimension at most
24 (or exactly 24 after harmless padding). This reduction is independent
of m and does not assume rank-two Y. It is a fixed-dimensional target,
not a verification of the inequality on that target.

## A sharp bridge case that is proved

The [qubit merged-gap theorem](induction16_qubit_merged_gap.md) proves
(MB_4) whenever X is a product F_free tensor Z_sector and rank X<=2,
for arbitrary H and either rank of Y. Its main estimate is

    sqrt(2||Z||^2-|Tr Z|^2)-sqrt(4q_2(Z))<=2sigma_2(Z)

on C^2 tensor H. The proof retains the actual two-sector form and
uses an explicit two-factor flatness sum of squares. The coefficient
four is sharp already in this product class. General X may entangle
the free and sector operator factors, and is not covered. This result
does not replace Gamma_bar by the full Pauli Sigma_m or prove the
universal bridge required in (6).

There is also a [uniform all-field rank-one-Y bridge](induction17_rank_one_Y_bridge.md):
coefficient 50 works for every actual datum with rank(Y)<=1 and arbitrary
H, including coherent scalar cancellation and simultaneous outside
fields. New-qubit contraction leaves Y unchanged, so (6) closes on this
entire class and proves actual star for every m>=13. The general
bridge is proved with coefficient 2129/10 in
[induction24](induction24_free_completion_uniform_bridge.md);
coefficient four is still open. The earlier coefficient-400 proof is
retained in [star22](star22_uniform_merger_bridge.md).
The needed pure weighted-kernel
estimate, valid with arbitrary free/sector entanglement, is proved in
[star20](star20_weighted_qubit_kernel.md).

## A limitation of a different energy comparison

It is not possible to bound Gamma_bar by a multiple of Gamma+Sigma_m.
For m=3, let Phi=(|00>+|11>)/sqrt(2) on the new qubit and the first
previous qubit, and set

    X=|0><0|_free tensor |Phi><Phi| tensor I_2.

This has rank two and exactly

    Gamma=0, Sigma_3=0, Gamma_bar=3.

The remaining previous factor I_2 kills every full three-qubit Pauli
coefficient. This does not refute the bridge: its left side is negative.
It only rules out an intermediate energy comparison that would discard
the compensation between the two terms of the D recursion.

`star19_energy_aware_recursion.py` checks the actual Pauli aggregation,
the complete cost-moment identity, the merged formulas, the coefficient
four witness and its rank-two-Y perturbation, and this zero-energy
example using exact SymPy arithmetic. None of these checks proves
(MB_4) or universal star. The separate induction24 proof makes (6)
unconditional with B=2129/10 and establishes actual star for every m>=17.

## The stronger squared-difference shortcut is false

Dropping the positive energy denominator would suggest the stronger
inequality

    q_1(Y) [Gamma-Gamma_bar]<=16pq.                       (7)

This already fails for a specialization of the mixed large-baseline
family in star9 (a working note not included here). Let U=V be the
common isometric embedding of C^2 into C^3, set c=|00>, w=|11>, and
write N=|0><1| and A=[[0,-1],[1,0]]. Choose

    R_(A,c)=hV, R_(A,w)=VN,
    L_(A,c)=hV, L_(A,w)=VN*,
    Y=V A^(-1) V*, h>0,

with every other A-sector block zero. The right-hand blocks are
L_(A,c)*=hV* and L_(A,w)*=N V*, so the Schur X=R_A A L_A* has

    X_cc=h^2 VAV*, X_cw=h V|1><1|V*,
    X_wc=h V|0><0|V*, X_ww=VNV*.

The full C is an actual rank-at-most-two operator; for example its
two-column factors are [R_A A;V] and [L_A;V(A^(-1))*]. Direct evaluation
gives

    p=q=1, q_1(Y)=2,
    Gamma=2h^4+4h^2+1,
    Gamma_bar=2h^4+2h^2+1.

Consequently the left side of (7) is 4h^2 and is unbounded. At h=3
it is exactly 36>16. The actual merged-bridge ratio, however, is

    sqrt(2) [sqrt(2h^4+4h^2+1)-sqrt(2h^4+2h^2+1)] ->1.

Thus this is not a counterexample to (MB_4). It shows why the
denominator sqrt(Gamma)+sqrt(Gamma_bar), retained in the original
square-root bridge, cannot be discarded. The example has no outside
components and uses only previous dimension two.

## One bounded adversary of the remaining bridge

`star19_merged_bridge_search.py` ran one authorized batch of 32
trajectories for 1,200 Adam steps, with previous dimension four,
unrestricted rank-two-Y Schur fields, unequal free planes, arbitrary
scalar sector vectors, simultaneous outside residuals, and M singular
value ratios down to approximately 10^(-3). Residual costs were
normalized before scalar baselines were added. CPU partial-trace versus
new-qubit Pauli evaluation and three directional gradient audits passed
before CUDA. All 32 saved trajectory-best points also passed an
independent CPU evaluation after the run.

The maximum merged-bridge ratio was 3.9880219661, below four, in
11.52 seconds. Its residuals were nearly in-plane: outside squared
fractions 3.41e-6 and 3.98e-6. Twelve saved trajectory-best points had
both outside fractions at least .01; the maximum among those was
2.3285669334. Thus the global maximum chiefly approaches the already
known rank-one-Y boundary, rather than providing strong evidence about
the full outside geometry. The search is finite evidence only.

Original two-column factors and all trajectory-best parameters are
archived in
`results/star19_merged_bridge_b32_s20261119.npz`, with independent audit
records in the companion JSON. There is no numerical violation to
certify and no additional job is running. The remaining analytic target
is the sharp candidate (MB_4) on the fixed 3-by-2-by-24 geometry.
The subsequent star22 argument proves the coarser (MB_400), improved
to (MB_(2129/10)) in induction24, by adding
a concurrence rectangle and an energy argument. Neither the cost moment
nor the two-copy positivity inputs supply that bridge by themselves.

Reproduce the exact algebra audit from the repository root with

    .\.venv-cuda\Scripts\python.exe research\star19_energy_aware_recursion.py

The single numerical run used

    .\.venv-cuda\Scripts\python.exe research\star19_merged_bridge_search.py

The search intentionally refuses to overwrite its existing archive;
pass `--output` with a fresh prefix to reproduce it. Its batch size
32 and 1,200 optimization steps are fixed in the script. The exact
coefficient-four lower bound and the finite near-four result should not
be confused with a proof of (MB_4).
