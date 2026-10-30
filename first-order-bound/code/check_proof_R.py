"""Step-by-step numerical check of the block lemma (Lemma 6 of the paper, eq. (4)) and of the main inequality.
Mapping to the paper: s0 = the balanced factorisation facts; s1, s2 = the two inequalities in eq. (5); s3 = a
cross-block form of eq. (6) (the paper uses only its diagonal case a = b); R = eq. (4); spade = Theorem 1.
Original description:
For random block operators Z (d blocks of size m), with the balanced SVD factorisation Z = P Q^*:
  (s0) Z_ab = P_a Q_b^*, sum_a R_a = sum_b S_b = Sigma, ||Z||_1 = Tr Sigma
  (s1) delta(Z_ab) := ||Z_ab||_1^2 - ||Z_ab||_2^2 <= sqrt((a^2-p)(b^2-q)),  a = Tr R_a, p = Tr R_a^2, b = Tr S_b, q = Tr S_b^2
  (s2) sqrt((a^2-p)(b^2-q)) <= a b - sqrt(p q)
  (s3) |<Z_aa, Z_bb>| <= (p_a p_b q_a q_b)^{1/4}
  (R)  sum_ab ||Z_ab||_1^2 + ||sum_a Z_aa||_2^2 <= sum_ab ||Z_ab||_2^2 + ||Z||_1^2
and (spade) for random Z on n sites."""
import itertools
import numpy as np
rng = np.random.default_rng(2026)
tn = lambda M: np.linalg.svd(M, compute_uv=False).sum()
fn2 = lambda M: np.real(np.vdot(M, M))
worst = {k: -np.inf for k in ('s0', 's1', 's2', 's3', 'R', 'spade')}

for trial in range(3000):
    d = int(rng.integers(2, 5)); m = int(rng.integers(1, 6))
    D = d * m
    kind = trial % 4
    if kind == 0:
        Z = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))
    elif kind == 1:
        r = int(rng.integers(1, 4))
        Z = (rng.normal(size=(D, r)) + 1j * rng.normal(size=(D, r))) @ (rng.normal(size=(r, D)) + 1j * rng.normal(size=(r, D)))
    elif kind == 2:
        al = rng.normal(size=d) + 1j * rng.normal(size=d); al /= np.linalg.norm(al)
        W = rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m))
        Z = np.kron(np.outer(al, al.conj()), W) + 10 ** rng.uniform(-4, -1) * (rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D)))
    else:
        A = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D)); Z = A @ A.conj().T
    U, s, Vh = np.linalg.svd(Z)
    keep = s > 1e-12 * s[0]
    U, s, Vh = U[:, keep], s[keep], Vh[keep]
    P = U * np.sqrt(s); Q = Vh.conj().T * np.sqrt(s)
    Pa = [P[a * m:(a + 1) * m] for a in range(d)]; Qa = [Q[a * m:(a + 1) * m] for a in range(d)]
    Zb = lambda a, b: Z[a * m:(a + 1) * m, b * m:(b + 1) * m]
    R = [p.conj().T @ p for p in Pa]; S = [q.conj().T @ q for q in Qa]
    Sig = np.diag(s)
    e0 = max(np.abs(Zb(a, b) - Pa[a] @ Qa[b].conj().T).max() for a in range(d) for b in range(d))
    e0 = max(e0, np.abs(sum(R) - Sig).max(), np.abs(sum(S) - Sig).max(), abs(tn(Z) - s.sum()))
    worst['s0'] = max(worst['s0'], e0 / s.sum())
    al_ = [np.real(np.trace(x)) for x in R]; p_ = [np.real(np.trace(x @ x)) for x in R]
    be_ = [np.real(np.trace(x)) for x in S]; q_ = [np.real(np.trace(x @ x)) for x in S]
    scale = s.sum() ** 2
    for a in range(d):
        for b in range(d):
            delta = tn(Zb(a, b)) ** 2 - fn2(Zb(a, b))
            bound1 = np.sqrt(max(al_[a] ** 2 - p_[a], 0) * max(be_[b] ** 2 - q_[b], 0))
            worst['s1'] = max(worst['s1'], (delta - bound1) / scale)
            worst['s2'] = max(worst['s2'], (bound1 - (al_[a] * be_[b] - np.sqrt(p_[a] * q_[b]))) / scale)
            ip = abs(np.vdot(Zb(a, a), Zb(b, b)))
            worst['s3'] = max(worst['s3'], (ip - (p_[a] * p_[b] * q_[a] * q_[b]) ** 0.25) / scale)
    lhs = sum(tn(Zb(a, b)) ** 2 for a in range(d) for b in range(d)) + fn2(sum(Zb(a, a) for a in range(d)))
    rhs = sum(fn2(Zb(a, b)) for a in range(d) for b in range(d)) + tn(Z) ** 2
    worst['R'] = max(worst['R'], (lhs - rhs) / scale)

# (spade) on random operators, several (d, n)
for trial in range(400):
    d = int(rng.integers(2, 4)); n = int(rng.integers(2, 5 if d == 2 else 4))
    D = d ** n
    Z = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))
    if trial % 2:
        r = int(rng.integers(1, 3)); Z = (rng.normal(size=(D, r)) + 1j * rng.normal(size=(D, r))) @ (rng.normal(size=(r, D)) + 1j * rng.normal(size=(r, D)))
    T = Z.reshape([d] * (2 * n))
    lhs = 0.0
    for i in range(n):
        Ti = np.trace(T, axis1=i, axis2=n + i)
        lhs += fn2(Ti)
    rhs = (n - 1) * fn2(Z) + tn(Z) ** 2
    worst['spade'] = max(worst['spade'], (lhs - rhs) / tn(Z) ** 2)

for k, v in worst.items():
    print(f'{k:6s}: max violation (normalised) = {v:+.3e}')
