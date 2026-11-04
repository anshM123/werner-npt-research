"""Check the charge-twirl splitting and the sector formula against direct partial traces.
(a) For X with range in H_1 and co-range in H_c: 2^{n-c+1} q_n(X) = sector form, n = 2..5, c = 1..n, random rank-two and
    random full-rank X.
(b) For X with range in H_1 and co-range spread over several sectors: q_n(X) = sum_c q_n(X P_c)."""
import numpy as np
from common import q_direct, strings, weight, single_excitations, sector_form, embed

rng = np.random.default_rng(2026)
worst_a = worst_b = 0.0
for n in range(2, 6):
    rows = single_excitations(n)
    for c in range(1, n + 1):
        cols = [st for st in strings(n) if weight(st) == c]
        for rank2 in (True, False):
            for trial in range(2 if n == 5 else 3):
                if rank2:
                    U = rng.normal(size=(len(rows), 2)) + 1j * rng.normal(size=(len(rows), 2))
                    V = rng.normal(size=(len(cols), 2)) + 1j * rng.normal(size=(len(cols), 2))
                    Xs = U @ V.conj().T
                else:
                    Xs = rng.normal(size=(len(rows), len(cols))) + 1j * rng.normal(size=(len(rows), len(cols)))
                lhs = 2 ** (n - c + 1) * q_direct(embed(Xs, n, rows, cols), n)
                rhs = sector_form(Xs, n, c, rows, cols)
                err = abs(lhs - rhs) / max(1.0, abs(lhs))
                worst_a = max(worst_a, err)
        print(f"(a) n={n} c={c}: max relative error so far {worst_a:.1e}", flush=True)
for n in range(2, 6):
    rows = single_excitations(n)
    cols = [st for st in strings(n)]
    U = rng.normal(size=(len(rows), 2)) + 1j * rng.normal(size=(len(rows), 2))
    V = rng.normal(size=(len(cols), 2)) + 1j * rng.normal(size=(len(cols), 2))
    Xs = U @ V.conj().T
    total = q_direct(embed(Xs, n, rows, cols), n)
    parts = 0.0
    for c in range(0, n + 1):
        mask = np.array([weight(st) == c for st in cols])
        parts += q_direct(embed(Xs * mask[None, :], n, rows, cols), n)
    err = abs(total - parts) / abs(total)
    worst_b = max(worst_b, err)
    print(f"(b) n={n}: q(X) = {total:+.10f}, sum over sectors = {parts:+.10f}, relative error {err:.1e}", flush=True)
print(f"RESULT: sector formula max relative error {worst_a:.1e}; charge splitting max relative error {worst_b:.1e}")
assert worst_a < 1e-9 and worst_b < 1e-9
print("PASS")
