"""End-to-end exact check of the certificates: for every cert_*.pkl given on the command line, at 4 random rational
weight vectors (s0 = 1) the block is built directly over Q (ylib point builder, independent of the symbolic pipeline),
the gauge coordinates are deleted, the states are put in the certificate order, and every leading principal minor is
computed by exact elimination over Q.  It must equal core_r * prod_i f_i^{x_ir - r E_i} / cst^r evaluated at the
point (the formula whose signs the certificate checks), and be positive.
Usage: python check_minors_points.py [j T cert.pkl ...]   (default: the seven certificates of NOTES.md)"""
import sys, pickle, random, time
import flint
from k4lib import block_states, gauge_keep, pt_builder, label

t0 = time.time()
args = sys.argv[1:] or ["1", "01", "cert_j1_T01_nat.pkl", "2", "012", "cert_j2_T012_nat.pkl", "1", "0", "cert_j1_T0_nat.pkl",
                        "2", "01", "cert_j2_T01_nat.pkl", "1", "-", "cert_j1_T-_nat.pkl", "2", "-", "cert_j2_T-_nat.pkl",
                        "2", "0", "cert_j2_T0_nat.pkl"]
ok = True
random.seed(21)
for q in range(0, len(args), 3):
    j = int(args[q])
    T = () if args[q + 1] == '-' else tuple(int(c) for c in args[q + 1])
    C = pickle.load(open(args[q + 2], 'rb'))
    ctx = flint.fmpz_mpoly_ctx.get(tuple(C['names']), 'lex')
    ctxq = flint.fmpq_mpoly_ctx.get(tuple(C['names']), 'lex')
    facs = [(ctx.from_dict(d), e) for d, e in C['den_factors']]
    cores = [(ctx.from_dict(d), ex) for d, ex in C['cores']]
    cst = C['den_const']
    S = block_states(j, T)
    keep = gauge_keep(S, j, T)
    labs = [label(S[a]) for a in keep]
    order = C['order']
    assert [labs[i] for i in order] == C['labels']
    n = len(order)
    good_all = True
    for trial in range(4):
        sv = [flint.fmpq(1)] + [flint.fmpq(random.randint(1, 60), random.randint(1, 60)) for _ in range(3)]
        Bd = pt_builder(sv)
        Mfull = Bd.matrix(S)
        A = [[Mfull[keep[order[a]]][keep[order[b]]] for b in range(n)] for a in range(n)]
        # exact leading minors by elimination
        mins = []
        Dr = flint.fmpq(1)
        W = [row[:] for row in A]
        for p in range(n):
            pv = W[p][p]
            Dr = Dr * pv
            mins.append(Dr)
            if pv == 0:
                break
            for i in range(p + 1, n):
                if W[i][p] == 0:
                    continue
                f = W[i][p] / pv
                for l in range(p + 1, n):
                    W[i][l] -= f * W[p][l]
        pt = sv[1:]
        # evaluate the certified formula with rational arithmetic
        good = len(mins) == n
        fvals = [(ctxq.from_dict({m: flint.fmpq(int(v)) for m, v in f.to_dict().items()})(*pt), E) for f, E in facs]
        for r, (c, ex) in enumerate(cores, 1):
            val = ctxq.from_dict({m: flint.fmpq(int(v)) for m, v in c.to_dict().items()})(*pt)
            for (fv, E), x in zip(fvals, ex):
                val = val * fv ** (x - r * E) if x - r * E >= 0 else val / fv ** (r * E - x)
            val = val / flint.fmpq(cst) ** r
            good &= r <= len(mins) and val == mins[r - 1] and val > 0
        good_all &= good
    ok &= good_all
    print(f'j = {j}, T = {T}, {args[q + 2]}: all {n} leading minors = certified formula and > 0 at 4 rational points: '
          f'{good_all}  [{time.time() - t0:.1f}s]', flush=True)
print('PASS' if ok else 'FAIL')
