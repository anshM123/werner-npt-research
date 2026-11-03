# Sharp coherent addition of two separable two-qubit Grams

Let `G(Z)=Z*Z`, let `||Z||` denote the Frobenius norm, and use
homogeneous two-qubit concurrence, normalized by
`C(|x><x|)=2|det x|`. The two matrices below may have any common
purification environment; no simultaneous product-row decomposition
is assumed.

**Theorem.** If `G(Z1)` and `G(Z2)` are separable, then

    C(G(Z1+Z2)) <= 2 ||Z1|| ||Z2||.                    (1)

The coefficient two is optimal. This proves the previously open generic
coherent-sum target in star23 (a working note not included here), Section 2.
It does not by itself prove the sharp physical rectangle, MB_4, ORB,
or the unrestricted all-copy endpoint statement.

The reduction below gives a scalar fidelity inequality, which has a
short crossed-coordinate proof.

## 1. A separable root-fidelity inequality

For positive matrices define root fidelity by

    F(S,T)=||sqrt(S)sqrt(T)||_*.

For every determinant-one concurrence witness

    W=(|vec A><vec A|)^Gamma, det A=1,

and every pair of separable positive two-qubit matrices `S,T`,

    ||sqrt(S) W sqrt(T)||_*
      <= sqrt(Tr S Tr T)+sqrt(Tr(WS) Tr(WT)).          (2)

The witness expectations in the second term are nonnegative, since
`W` is nonnegative on product vectors. Singular states are allowed.

It suffices by homogeneity to prove (2) for `Tr S=Tr T=1`. Local
Schmidt unitaries put the witness in the form

    W_r = [[r,0,0,0],
           [0,0,1,0],
           [0,1,0,0],
           [0,0,0,1/r]],                 r>0.        (3)

These unitaries preserve separability and all quantities in (2).
Write

    c(S,T)=||sqrt(S) W_r sqrt(T)||_*
          =F(S,W_r T W_r),
    a=Tr(W_r S), b=Tr(W_r T).

Root fidelity is jointly concave. One self-contained way to see the
direction used here is its block-positive characterization:

    F(P,Q)=max Re Tr X,
           [[P,X],[X*,Q]] >= 0.

Indeed the block condition is equivalent to
`X=sqrt(P) K sqrt(Q)` for a contraction `K`, with singular cases
obtained by continuity; maximizing the trace gives the nuclear norm.
Convex combinations of feasible blocks prove joint concavity.
Consequently `c` is jointly concave in `S,T`.

### First reduction: four-phase twirling

Let `U=diag(1,i)` and average both states over conjugation by
`U^j tensor U^j`, `j=0,1,2,3`. Every such unitary commutes with `W_r`.
Thus this common twirl increases `c`, preserves `a,b`, and preserves
separability. The resulting states have the forms

    S=[[x,0,0,0],          T=[[X,0,0,0],
       [0,s1,zeta,0],         [0,t1,xi,0],
       [0,bar(zeta),s2,0],    [0,bar(xi),t2,0],
       [0,0,0,y]],            [0,0,0,Y]].

Their traces are one. Positivity and PPT give, respectively,

    |zeta|^2<=s1 s2,       |zeta|^2<=xy,
    |xi|^2<=t1 t2,         |xi|^2<=XY.                (4)

The formula for the root fidelity of two positive `2 by 2` matrices,
`F(P,Q)^2=Tr(PQ)+2sqrt(det P det Q)`, gives exactly

    c=r sqrt(xX)+r^(-1)sqrt(yY)
      +sqrt(s1 t2+s2 t1+2Re(zeta xi)
             +2sqrt((s1s2-|zeta|^2)(t1t2-|xi|^2))).  (5)

The product is `zeta xi`, without a conjugate: the middle block of
`W_r T W_r` is obtained by swapping the two middle coordinates.

### Second reduction: negative real coherences

Replace `zeta` by `-|zeta|` and `xi` by `-|xi|`. This increases (5),
because `Re(zeta xi)<=|zeta||xi|`. It decreases both witness
expectations, since

    a=rx+r^(-1)y+2Re(zeta),
    b=rX+r^(-1)Y+2Re(xi).

The new states remain positive and PPT by (4), hence separable in
two-qubit dimensions. Their new `a,b` remain nonnegative; for example

    rx+r^(-1)y-2|zeta|
      >= (sqrt(rx)-sqrt(y/r))^2 >=0.

It is therefore enough to prove the desired upper bound for these
new states: the left side has increased and the right side has
decreased.

### Third reduction: common SWAP averaging

Average both new states with their conjugates by SWAP. SWAP commutes
with `W_r`, and it preserves separability. Joint concavity again
increases `c` while preserving the new witness expectations. The
middle diagonal entries become equal; the negative real coherences
remain unchanged. Put

    u=(s1+s2)/2-|zeta|, v=(s1+s2)/2+|zeta|,
    U=(t1+t2)/2-|xi|,   V=(t1+t2)/2+|xi|.

These numbers are nonnegative, and

    x+y+u+v=X+Y+U+V=1,
    v-u<=2sqrt(xy), V-U<=2sqrt(XY).                  (6)

The final middle blocks commute, so the three quantities are

    c=r sqrt(xX)+r^(-1)sqrt(yY)+sqrt(uU)+sqrt(vV),
    a=rx+r^(-1)y+u-v,
    b=rX+r^(-1)Y+U-V.                               (7)

### The scalar inequality

Set the possibly signed real numbers

    A0=sqrt(rx)-sqrt(y/r),
    B0=sqrt(rX)-sqrt(Y/r).

By (6), `a>=A0^2` and `b>=B0^2`. The exact cancellation

    c-A0 B0=sqrt(xY)+sqrt(yX)+sqrt(uU)+sqrt(vV)<=1   (8)

is Cauchy--Schwarz applied to the two unit vectors
`(sqrt x,sqrt y,sqrt u,sqrt v)` and
`(sqrt Y,sqrt X,sqrt U,sqrt V)`. Hence

    c<=1+A0B0<=1+|A0B0|<=1+sqrt(ab).

No common sign assumption on `A0,B0` is needed. This proves the
normalized version of (2). Homogeneity proves the general case.

## 2. Arbitrary relative purification and concurrence

Put `S=G(Z1)`, `T=G(Z2)`, and `zj=||Zj||`. For each fixed witness
let `a=Tr(WS)>=0`, `b=Tr(WT)>=0`. Polar factorization of the two
purifications gives

    |Tr(W Z1*Z2)|<=||sqrt(S) W sqrt(T)||_*.

Here the intervening environment overlap is a contraction. There
is no requirement that either given purification have product rows,
and no requirement of a common change of rows. By (2),

    -Tr[W G(Z1+Z2)]
      <=-a-b+2(z1z2+sqrt(ab))
       =2z1z2-(sqrt(a)-sqrt(b))^2
      <=2z1z2.                                      (9)

The homogeneous two-qubit concurrence dual is

    C(B)=max(0,sup_(det A=1) -Tr[W_A B]).

See the repository dual note (a working note not included here), or
Eq. 66 of [Eisert, Brandao and Audenaert](https://arxiv.org/html/quant-ph/0607167).
Taking this supremum in (9) proves (1); no dual maximizer is needed.
Zero purifications cause no exception.

For sharpness take the one-row purifications `Z1=<00|`, `Z2=<11|`.
Their individual Grams are product projectors, their norms are one,
and their sum has homogeneous concurrence two.

## Verification and scope

The companion checker verifies the finite twirl, commutations, the
complex middle-block fidelity trace and determinant identities, the
SWAP reduction, all scalar cancellations and the four-dimensional
Cauchy--Schwarz sum of squares, the final square in (9), a boundary
control, and the sharp example. It performs no optimization or search.
Concavity and the stated PSD/PPT implications are proved above rather
than inferred from sample calculations.

This theorem is generic in the two separable Gram matrices and their
relative environment. It improves an ingredient of the concurrence
rectangle argument. Any improvement to the final merger constants
must still be derived through that argument. The unresolved sharp
physical rectangle, sharp merger bound, and general EP/ORB problems
are not conclusions of (1).
