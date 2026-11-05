"""Helpers for the coefficient-sign certificates (exact arithmetic, python-flint)."""
import flint


def pure_positive(P, pos):
    """True if P has a monomial with positive coefficient whose exponents vanish outside the index set pos."""
    return any(cf > 0 and all(ex[k] == 0 for k in range(len(ex)) if k not in pos) for ex, cf in P.to_dict().items())


def nonneg_report(P, pos):
    """(number of terms, number of negative coefficients, has a positive monomial in the variables pos only)."""
    d = P.to_dict()
    neg = sum(1 for cf in d.values() if cf < 0)
    return len(d), neg, pure_positive(P, pos)


def positive_denominator(D, pos):
    """D is a product of a positive constant and polynomials with only positive coefficients, each containing a
    monomial in the variables pos only; then D > 0 whenever the variables pos are > 0 and the others are >= 0."""
    cst, fac = D.factor()
    ok = cst > 0 and all(all(cf > 0 for cf in f.coeffs()) and pure_positive(f, pos) for f, _ in fac)
    return ok, [(len(f), f.total_degree(), m) for f, m in fac]


def amgm_repair(P, ia, ib):
    """For every negative coefficient -c at a monomial m * x^k y^k (x, y the variables ia, ib, k >= 1), subtract
    (c/2) m x^(k-1) y^(k-1) (x - y)^2, which is >= 0 for x, y >= 0.  Returns (remainder, list of repairs).
    P = remainder + sum of the subtracted terms, so P >= 0 if the remainder has only nonnegative coefficients."""
    ctx = P.context()
    gens = ctx.gens()
    x, y = gens[ia], gens[ib]
    rem = P
    reps = []
    for ex, cf in P.to_dict().items():
        if cf < 0 and ex[ia] == ex[ib] and ex[ia] >= 1:
            k = ex[ia]
            ex2 = list(ex); ex2[ia] = k - 1; ex2[ib] = k - 1
            mono = ctx.from_dict({tuple(ex2): 1})
            term = (-cf / 2) * mono * (x - y) * (x - y)
            rem = rem - term
            reps.append((ex, cf))
    return rem, reps
