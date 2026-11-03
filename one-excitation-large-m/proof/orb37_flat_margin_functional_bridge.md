# Flat-margin functional bound 163/10 and uniform merger 125

The sharp merged flat margin of [orb33](orb33_universal_merged_flat_margin.md)
improves the unrestricted functional coefficient 22 of
[orb27](orb27_full_safe_functional.md) to 163/10. Its conversion through
the exact canonical gauges of [orb28](orb28_combined_uniform_bridges.md)
gives the dimension-independent merger bounds

    sqrt(q1(Y))[sqrt(Gamma(X))-sqrt(Gamma_bar(X))]
       <=125 sqrt(qtilde(P)qtilde(Q)),                    (1)

and, when rank(Y)<=1, the stronger coefficient

    (163/10)sqrt(2)<24.                                  (2)

The energy-aware recurrence therefore proves actual original star
for every m>=16, and for every m>=11 when rank(Y)<=1. The proper
one-free-copy, one-excitation endpoint classes have n=m+1>=17 and
n>=12, respectively. General m=3 through 15, rank-one-Y m=3 through
10, ORB and unrestricted induction remain open. No all-n NPT theorem
is asserted.

The written proof and exact checker were independently reviewed in
induction36 (a working note not included here).

## 1. Both actual forms have the same sharp flat margin

Use the functional variables and hypotheses of orb27:

    F=[f1,f2], L0=U tensor c, L=L0+delta, X=FL*,
    alpha=||f1||, T=||f2||, S=||c||, beta=||delta||.

All fields are unrestricted and the preceding space is any finite
dimension. On rank-at-most-two operators, the two-factor theorem and
the new merged theorem give

    Gamma(Z), Gamma_bar(Z)
       >=[sigma1(Z)-sigma2(Z)]²/2.                        (3)

If S>beta, then L*L>=(S-beta)²I, while ||L||op<=S+beta.
Thus the ordinary product singular-value bounds yield

    sigma1(FL*)>=sigma1(F)(S-beta)>=T(S-beta),
    sigma2(FL*)<=sigma2(F)(S+beta)<=alpha(S+beta).          (4)

The last inequality also follows by removing the first column of F
and using the best rank-one approximation inequality. When the
following right-hand side is positive, (3)-(4) imply

    sqrt(Gamma(X))+sqrt(Gamma_bar(X))
       >=sqrt(2)[TS-alpha S-beta(T+alpha)].               (5)

This estimate uses the actual singular values of X. It does not
apply a seminorm triangle to a rank-four error.

## 2. Functional coefficient 163/10

The unchanged valid full-purification curvature bound from orb27 is

    -H(X)<=8sqrt(2)alpha beta TS+4alpha T beta²,
    H=Gamma_bar-Gamma.                                   (6)

When alpha or beta vanishes, the established rank-one or scalar
safe-axis theorem gives a nonpositive gap. Otherwise split into the
following exhaustive cases.

If T<=(14/3)alpha, the existing comparison with the full safe
operator F L0* gives

    gap<=(sqrt(2)+2)sqrt(T²+alpha²)beta
        <=(sqrt(2)+2)sqrt(205)alpha beta/3
        <(163/10)alpha beta.                             (7)

If S<=(19/4)beta, the universal auxiliary coefficient 2sqrt(2)
from orb33 and the product bound in (4) give, without requiring
S>beta,

    gap<=2sqrt(2)sigma2(X)
        <=2sqrt(2)alpha(S+beta)
        <=(23/2)sqrt(2)alpha beta
        <(163/10)alpha beta.                             (8)

In the remaining case, set u=alpha/T<3/14 and v=beta/S<4/19.
The bracket in (5), divided by TS, is at least

    1-3/14-(4/19)(1+3/14)=141/266>0.

For H>=0 there is no positive gap. Otherwise divide (6) by (5):

    gap<=alpha beta [8sqrt(2)+16/19]/[141sqrt(2)/266]
        <(2288/141)alpha beta
        <(163/10)alpha beta.                             (9)

The first strict comparison uses 1/sqrt(2)<5/7. For (7), the
rational bounds sqrt(2)<99/70 and sqrt(205)<7159/500 give
an upper coefficient 1711001/105000<163/10. For (8), use
2277/140<163/10. All radical bounds follow by squaring positive
integers. Hence (7)-(9) prove, in arbitrary previous dimension,

    sqrt(Gamma(F L*))-sqrt(Gamma_bar(F L*))
       <=(163/10)alpha beta.                             (10)

## 3. Convert the functional bound to the general merger

The factor gauges, degeneracies and exact full-purification estimate
are precisely those already proved in orb28; only constants change.
For rank(Y)<=1, its Section 1 applies (10) in place of coefficient
22 and gives (2). In particular (163/10)sqrt(2)<24.

For rank(Y)=2, use the unchanged normalized variables of orb28:

    y=q1(Y)<=2, alpha=sqrt(p), beta=sqrt(q),
    A=||a||, D=||c||, s>=1,
    X=(V tensor a+r)diag(1,s)(U tensor c+delta)*.

Take tau=25/4. If A sqrt(y)<=tau alpha, or the adjoint case
D sqrt(y)<=tau beta, (10) gives the scaled gap at most

    (163/10)(tau+sqrt(2))alpha beta
       <(174899/1400)alpha beta<125alpha beta.           (11)

In the complementary case the old denominator comparison and old
concurrence estimate are retained in full. They give

    scaled gap<=alpha beta [21+16/tau+2/tau²]/d_tau,
    d_tau=sqrt(2)-(sqrt(2)+2)(2/tau+sqrt(2)/tau²)
          =(393sqrt(2)-432)/625,
    21+16/tau+2/tau²=14757/625.                          (12)

The exact rational lower bound sqrt(2)>7/5 gives

    125d_tau-14757/625>18/625>0.                         (13)

This proves d_tau>0 and the coefficient 125 in (1).
The zero-cost and singular-factor limits require no new argument:
they are included in the uniform functional estimate and the same
canonical perturbation/adjoint proof of orb28. No condition-number
bound on Y is imposed.

## 4. Recurrences and precise scope

Starting from K1=0, the unchanged energy-aware recurrence is

    Km<=sqrt(2)K_(m-1)+125.

At m=16 its geometric sum gives

    K16<=125(255+127sqrt(2))
        <111375/2<65534=2^16-2.                          (14)

For m>=16, the next-step comparison follows from

    (3/2)(2^m-2)+125=3*2^(m-1)+122<=2^(m+1)-2,

since 2^(m-1)>=124. Thus actual star holds for every m>=16.

The recurrence preserves Y, so the entire rank(Y)<=1 subclass
uses the rounded coefficient 24 at every step. At m=11,

    K11<=744(1+sqrt(2))<1860<2046=2^11-2.                (15)

Propagation requires only 2^(m-1)>=23, automatic thereafter.
These are the large-m conclusions stated at the beginning; they do
not imply smaller-copy results or unrestricted endpoint positivity.

## Verification

The companion checker reruns the exact sharp-flat-margin audit and
the valid functional-curvature audit, checks the Gram comparisons
underlying (4) on a complex nonorthogonal-column control with entangled sector c, and
checks all new constants, polynomial expressions and recurrences in
exact arithmetic. The dimension-independent singular-value argument
and exhaustive cases above establish the theorem; the finite control
is a normalization check. There is no optimization or numerical sweep.
