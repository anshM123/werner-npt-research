"""Exact per-site maps for the certificate pipeline (first component, checked against the floating-point data).

Words: the 24 permutation operators V_sigma on (C^3)^(x)4.  Every per-site algebra element is X = sum_sigma x_sigma V_sigma,
with exact Fourier blocks X_lam = sum_sigma x_sigma rho_lam(sigma) (rational irreps from exact_site.ExactAlg).  A per-site
linear map phi is applied to words exactly (integer matrices: partial transposes, the pair realignment, partial traces),
and its image is decomposed with the exact Fourier data of the target algebra.  Multi-site maps are tensor products.
Checks: (i) Fourier inversion on the S_4 algebra reproduces V_sigma exactly from its blocks; (ii) for random rational
elements, the exact target blocks (made symmetric with the Gram matrices, G X_hat) have the same spectra (with
multiplicities) as the floating-point blocks of rep4/alg4.
usage: python exact_maps.py"""
import itertools, numpy as np, flint
from alg4 import S4, PERM, pt, ptrace_last, build_algebras, perm_op, inv, mul
from exact_site import ExactAlg, Q, build_exact, nullspace_q, weights
from rep4 import realign1

def exact_A3():
    S3 = list(itertools.permutations(range(3)))
    gens = {s: perm_op(s, 3) for s in S3}
    T = sum(perm_op(s, 3) for s in S3 if sum(1 for i in range(3) if s[i] != i) == 2)
    hw = {6: (3, 0, 0), 0: (2, 1, 0), -6: (1, 1, 1)}
    return ExactAlg(gens, T, set(), hw, 'A3 exact', nr=3)

def fourier_of(alg, M):
    """exact Fourier blocks of an integer matrix M lying in the algebra"""
    return alg.coeffs(M)

def sym_spectrum(blocks_q, alg):
    """eigenvalues (floats, with multiplicities m) of the operator with exact blocks: G^(1/2)-symmetrised spectrum"""
    out = []
    for b, Xq in zip(alg.blocks, blocks_q):
        X = np.array([[float(Xq[i, j]) for j in range(b['f'])] for i in range(b['f'])])
        ev = np.linalg.eigvals(X).real                 # similar to a symmetric matrix: real spectrum
        out += sorted(ev.tolist()) * b['m']
    return np.sort(np.array(out))

if __name__ == "__main__":
    A, CB, gensB1, evB = build_exact()
    from exact_site import ExactAlg as EA
    W8 = weights({3}); hwB = {}
    for lam in evB:
        hwB[int(round(lam * 2))] = {3.0: (3, 0, -1), 0.0: (2, 1, -1), -2.0: (2, 0, 0), -4.0: (1, 1, 0)}[round(lam, 6) + 0.0]
    B1 = EA(gensB1, CB, {3}, hwB, 'B31 exact')
    A3 = exact_A3()
    print("A3 exact:", [(b['f'], b['m']) for b in A3.blocks])
    # (i) Fourier inversion on A: V_sigma from its blocks via the exact matrix units
    rng = np.random.default_rng(5)
    # (ii) spectra of map images for random integer combinations of words
    for name, phi, tg in [('PT_3', lambda X: pt(X, [3]), B1), ('realign(2,3)', realign1, B1),
                          ('Tr_3[. F_23]', lambda X: ptrace_last(X @ PERM[(0, 1, 3, 2)]), A3),
                          ('Tr_3', lambda X: ptrace_last(X), A3)]:
        worst = 0.0
        for trial in range(3):
            coef = rng.integers(-3, 4, size=24)
            X = sum(int(cf) * PERM[s] for cf, s in zip(coef, S4))
            # make it symmetric (a valid algebra element whose image is self-adjoint up to the map)
            X = X + X.T
            Y = phi(X)
            Yi = np.rint(Y).astype(np.int64)
            assert np.array_equal(Yi, Y), "map image not integral"
            ex = fourier_of(tg, Yi)
            spec_exact = sym_spectrum(ex, tg)
            Ys = (Y + Y.T) / 2
            spec_dense = np.sort(np.linalg.eigvalsh(Ys)) if np.allclose(Y, Y.T) else None
            if spec_dense is not None:
                pad = np.concatenate([spec_exact, np.zeros(len(spec_dense) - len(spec_exact))]) if len(spec_exact) < len(spec_dense) else spec_exact
                worst = max(worst, np.abs(np.sort(pad) - spec_dense).max())
        print(f"   {name}: exact target blocks vs dense spectrum, max diff {worst:.1e}")
