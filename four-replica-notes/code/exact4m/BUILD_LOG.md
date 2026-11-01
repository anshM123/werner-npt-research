EXACT4M: exact rational certificate for q_4(P_u - P_v) >= 0, u, v orthonormal in (C^3)^(x)4 (t = -1), two-letter relaxation
==============================================================================================================================

Scope: the t = +1 pipeline of ../exact4 (EXACT4_NOTES.md) adapted to opposite signs.  Same relaxation: X_c =
twirl((x)_r |c_r><c_r|), c in {uuuu, uuuv, uuvv, uvvv, vvvv}; X >= 0, moments, swap identities (built into the variable
space), ORTH, MARG, covariance (COV) blocks; no partial-transpose blocks.  Only the objective changes:
  q_4(P_u + t P_v) = (fu + t^2 fv + 2 t fuv) / 2^4,   here t = -1:   c = (fu + fv - 2 fuv) / 2^4.
True points use P = P_u - P_v everywhere (points, cert, verify).  Nothing in ../exact4 or ../relaxations was
modified; the per-site pieces are imported from ../relaxations exactly as in ../exact4.

Log (times from the shell clock, EST)
-------------------------------------
S0 (2026-11-01 08:15) started; read EXACT4_NOTES.md, the exact4 pipeline and the floating-point t = -1 run
  (../relaxations/fr_tm1_qubit.py, ../../logs/fr_tm1_qubit_n4.log: X + identities + COV with the zeros k2 (2), k4 (2), inner (4) and
  6 real spin-flip zeros: Solved, smallest reduced-dual eigenvalue 3.9e-2; without the spin-flip zeros 5.4e-10).

E1 (2026-11-01 08:19) folder and exact linear data.
  core4.py, engine4.py, symm4.py, peakmem.py copied unchanged from ../exact4 (md5 identical).  The engine, the
  variable space, the symmetry-adapted COV compressions and all constraint data are therefore those validated in
  ../exact4 (test_engine.py, crosscheck4.py there).
  build_lin4m.py (= build_lin4.py with the t = -1 objective; 11 s, 0.49 GB private): 1331 parameters, 150 groups,
  6186 raw identity rows of rank 146 mod 1048573 (exact rank 146), 147 independent equalities E x = e (entries <= 37
  bits).  c has 689 nonzero coordinates, 599 of them differ from the t = +1 objective.  The t = +1 objective is
  saved as cplus (checks only).  Exact comparison with ../exact4/lin4.pkl: E, e and the 120 moment rows are identical,
  cplus equals the certified t = +1 objective, c differs.

E2 (2026-11-01 08:32) exact zeros and entangled test points (points4m.py; launched 08:22, 8 processes at a time,
  each <= 0.24 GB private, 1.5-8 min each, all 18 done in 10 min).
  points4m.py = points4.py (same moment code: 27-bit primes with the asserted int64 bound, CRT, exact Fourier
  transform, exact projection) with the t = -1 families; all vectors real and integer:
    k2     u = a (x) b (x) psi, v = a2 (x) b2 (x) psi: (a, a2), (b, b2) random ordered pairs of rows of integer Cayley
           orthogonal matrices (cayley_int), psi in {-3..3}^9 redrawn until det(psi as 3 x 3) != 0 (entangled).
           Zero because q_n(P_{a (x) phi} + t P_{a2 (x) chi}) = q_{n-1}(P_phi - t P_chi) / 2 for a orthogonal to a2:
           q_4(P_u - P_v) = q_3(P_{b psi} + P_{b2 psi}) / 2 = q_2(P_psi - P_psi) / 4 = 0.
    k4     orthogonal products on all four sites: q_4 = q_1(P_d + P_d2) / 8 = 0.
    inner  u = a (x) phi, v = a2 (x) chi, (phi, chi) = (r1 p + r2 q, -r2 p + r1 q), p = b (x) c (x) d,
           q = b2 (x) c2 (x) d2, (r1, r2) from PYTH with random order and sign: q_4 = q_3(P_phi + P_chi) / 2 and
           P_phi + P_chi = P_p + P_q, a t = +1 zero at three copies.
    sf     u supported on {0,1}^4, u_x in {-1, 0, 1} (rng.integers(-1, 2)) on the 7 pairs {x, xbar} with x_1 = 0,
           u_0000 = 1, u_1111 = -S, S = sum over those 7 pairs of (-1)^w(x) u_x u_xbar; v = Y u with
           Y|x> = (-1)^w(x)|xbar> (sigma_y^(x)4, a real symmetric involution).  Then u^T Y u = 0, so <u, v> = 0 and
           |u| = |v|; P_v = Y P_u Y has the same full-weight Pauli correlations as P_u (even n), so q_4(P_u - P_v) = 0.
    ent    u = p, v = (p.p) q - (p.q) p, p, q in {-2..2}^81 (as in exact4).
  Seeds and names: k2 101 z_k2_a, 102 z_k2_b; k4 111 z_k4_a, 112 z_k4_b; inner 121 z_in_a (pair 4,-3), 123 z_in_b
  (8,15), 124 z_in_c (20,-21), 126 z_in_d (24,7); sf 131..138 z_sf_a..z_sf_h; ent 197 p_ent, 198 p_ent2.
  The inner seeds were chosen so that the four Pythagorean pairs are distinct: x is invariant under local
  orthogonal maps, and an inner zero depends (up to them) only on the unordered pair {|r1|, |r2|} (seeds 121 and 122
  both give (4,-3), so 122 was skipped).  For the same reason all k4 zeros give one point: x(z_k4_a) = x(z_k4_b)
  exactly (an exact check of the local invariance of the pipeline).
  Every zero (16): in the variable space (exact reconstruction), all 120 moments exact, E x = e exact, objective
  c . x exactly 0, and q_4(P_u - P_v) computed directly from the 81 x 81 matrix exactly 0.  Moment primes: sf 1-2,
  k2 3-4, k4 5, inner 5-7 (largest entries: inner 43200, k4 5544, sf 2); x up to 65 bits.
  Objective cross-check at the entangled points (exact): c . x(P) = q_4(P_u - P_v) = 0.920386902566... (p_ent) and
  0.997269114416... (p_ent2), q_4 recomputed directly from P = P_u - P_v; and cplus . x(P) = q_4(P_u + P_v)
  (0.973144905689..., 0.946926709096...), so the three pieces fu, fv, fuv enter with the right weights.  The cplus
  check also holds at all 16 zeros.  x of p_ent2 up to 112 bits.

E3 (2026-11-01 08:36) exact face at the 16 zeros, saturation, face-reduced data.
  face4m.py (19 s, 0.65 GB private): x* = sum_z lambda_z x_z integer, up to 59 bits.  X blocks: dimension 709 -> 493.
  COV class blocks: |Yhat(x*)| <= 2^105, |What| <= 2^129, 5 primes; dimension 5028 -> 4119 (182 class blocks with
  W(x*) = 0); exact integer kernel bases (entries <= 3 bits), W N = 0 checked exactly.
  Saturation (face4m.py --nosave; a subset gives a face containing the full one, so equal dimension means equal face;
  X / COV dimensions):
    without z_sf_g, z_sf_h (6 sf)          493 / 4119   (identical)
    only z_sf_a..d (4 sf)                  493 / 4119   (identical)
    without z_in_c, z_in_d (2 inner)       493 / 4119   (identical)
    structured only (no sf)                552 / 4285   (larger face: the spin-flip zeros are needed, as in floats)
    sf only                                507 / 4205
  facedata4m.py (98 s, 3.34 GB private): 379 COV face blocks (max r 64, svec 45495), |Fhat| <= 2^65, 3 primes;
  168 COV and 104 X face blocks with nonzero data; exact face check: F_k(x_z) = 0 for all 272 blocks and all 16
  zeros (0 failures); numerators <= 46 bits.

E4 (2026-11-01 08:37) vacuous directions removed exactly.
  span4m.py (26 s, 1.64 GB private): 165 of the 272 blocks have directions on which all data vanish; total block
  dimension 2784 -> 1564, largest block 32, svec rows 9592; B entries <= 3 bits, numerators <= 48 bits; row-space
  completeness checked exactly.

E5 (2026-11-01 08:39) n = 4, t = -1: CERTIFIED (cert4m.py on facedata4mr.pkl, 40 s, 2.49 GB private; independent
  re-check verify4m.py).
  Float problem: 147 equality rows + 9592 svec rows (61 scalar + 211 PSD blocks, largest 32) x 1331 parameters;
  33 null directions removed by whitening (objective component 2.0e-6 against ||c|| = 4.3e5, relative 5e-12).
  Clarabel (chordal decomposition off): Solved in 17 iterations, primal -1.0e-12, dual -2.4e-15, ||A^T z + c||
  1.1e-9, smallest reduced-dual eigenvalue 0.09999 (Jacobi-normalised blocks; attained on Xuuuv(3,3,3,3), size 14;
  next 0.149 on Xuuvv(3,3,3,3)); solve 23 s.
  Rounding to dyadic grids (40 bits below each block's largest entry): exact residual 4.5e-5 in the identity,
  1.2e-12 in e . y.  Column subset by pivoted QR on the 1332 x 9663 float equation matrix: rank 1298 of 1332 (diag(R)
  ratio 1.3e-3 at the cut, 2.2e-15 next), 57 multipliers y and 1241 dual entries; exact 1298 x 1298 solve (flint,
  Dixon) in 1 s; max |delta| 2.2e-5.  Afterwards ALL 1332 equations hold exactly: c = E^T y + sum_k F_k^*(S_k) and
  e . y = 0.
  Exact LDL^T: all 272 dual blocks PSD (smallest pivot 3.96e-4).  Certificate entries up to 605/605 bits (y up to
  604/605 bits); data entries: E <= 37 bits, face-block numerators <= 48 bits.
  verify4m.py (reads only lin4m.pkl, facedata4mr.pkl, cert4m.pkl, pts/*.pkl; imports no pipeline code; q_4 computed
  by its own partial-trace routine on the integer matrix |u|^2 |v|^2 (P_u - P_v), which agrees with
  exact_coreN.exact_qn on 5 test pairs): (a) e . y = 0 and the identity in all 1331 coordinates exactly; (b) all 272
  dual blocks PSD by exact LDL^T; (c) at p_ent and p_ent2: u . v = 0, q_4(P_u - P_v) recomputed from P equals the
  stored value, c . x(P) = q_4(P_u - P_v), E x(P) = e, c . x(P) = y . E x(P) + sum_k Tr(S_k F_k(x(P))) with all 272
  terms >= 0, every face block F_k(x(P)) PSD (exact LDL^T): VERIFIED (9 s, 0.58 GB private).
  Negative controls (negctl4m.py, then verify4m.py --cert ...): one multiplier changed (y_0 + 1/1024, e_0 = 1/50625):
  e . y = 0 fails, identity fails in 1154 coordinates, identity at p_ent fails -> VERIFICATION FAILED; one dual
  diagonal entry negated (block COVuv(2,2,3,3)2|2, entry 2.97): identity fails in 288 coordinates, block not PSD ->
  VERIFICATION FAILED.
  Observation from the second control: the identity at p_ent still held there, because for that block the
  functional x -> F_k(x)_00 lies in the row space of [E | e] with right-hand side 0 (exact rank test), i.e. this
  diagonal entry (and with it the whole row, by PSD-ness at true points) vanishes on the affine space E x = e; three
  diagonal entries of this block vanish at both entangled points.  So the face could be reduced further using the
  equalities.  This does not affect validity (the identity is exact for all x) and the reduced dual is strictly
  feasible as it stands.
  crosscheck4m.py 25 (51 s, 5.75 GB private): 25 random covariance face blocks, 7,030,009 numerator entries recomputed
  with Python integers only: all equal to the CRT values; the reduced data equal B Fhat B^T.
  Statement certified: q_4(P_u - P_v) >= 0 for every orthonormal pair u, v in (C^3)^(x)4.  Real pairs: x(u, v) is
  feasible and c . x(u, v) = q_4(P_u - P_v), so q_4 = y . e + sum_k Tr(S_k F_k(x)) >= 0.  Complex pairs by the
  realification argument of the n = 2, 3 and t = +1 runs: the real part of the twirled two-letter moments of (u, v)
  is the average of the true points of (u, v) and (conj u, conj v), both orthonormal pairs, hence feasible, and its
  objective value is (q_4(P_u - P_v) + q_4(conj(P_u - P_v))) / 2 = q_4(P_u - P_v).

E6 (2026-11-01 08:44) reproducibility.
  logs/rerun4m.sh (build_lin4m -> face4m -> facedata4m -> span4m -> cert4m -> verify4m on the saved points;
  rerun4m.log, 3.3 min): same numbers and CERTIFIED / VERIFIED again; lin4m.pkl, face4m.pkl, facedata4m.pkl,
  facedata4mr.pkl, float_dual4m.pkl and cert4m.pkl are byte-identical to the first run (md5 in logs/md5_first_run.txt
  and at the end of rerun4m.log).
  Peak memory (private bytes): build 0.49 GB, points <= 0.24 GB each (8 at a time), face 0.65 GB, facedata 3.34 GB,
  span 1.64 GB, cert 2.49 GB, verify 0.58 GB, crosscheck 5.75 GB.  Wall: build 11 s, points 1.5-8 min each (18 in
  10 min), face 19 s, facedata 98 s, span 26 s, cert 40 s, verify 9 s, crosscheck 51 s.
  Not rerun here: test_engine.py (engine, variable space and constraint maps are the byte-identical ../exact4 modules,
  validated there against exact_coreN); the only new linear datum, the objective, satisfies c . x = q_4(P_u - P_v)
  exactly at all 18 points (nonzero at the two entangled ones), and cplus . x = q_4(P_u + P_v) at all 18 as well.
  Files: cert4m.pkl (y and S_k as (num, den)); cert4m_neg_y.pkl, cert4m_neg_S.pkl are the deliberately corrupted
  copies of the negative controls (regenerated by negctl4m.py); logs/launch_points.log is the scheduler log of the 18
  points4m.py runs (at most 8 at once), logs/points_pids.txt their wrapper PIDs.

E7 (2026-11-01 09:01) independent audit of the saved certificate.
  verify4m.py re-run at p_ent and p_ent2 (VERIFIED), and at five further exact true points made afterwards with
  points4m_audit.py from fresh seeds (../../logs/exact4m/audit_*.log): b_ent1-3 (generic entangled pairs,
  q_4(P_u - P_v) = 0.939, 0.946, 1.000) and b_nsf1, b_nsf2 (a spin-flip zero scaled by 4 plus a small integer
  perturbation, then v orthogonalised against u exactly; q_4 = 0.299, 0.384).  VERIFIED at all seven points.
  Then (09:11) b_nk21 (perturbed k2 zero, q_4 = 5.12e-6) and b_nin1 (perturbed inner zero, q_4 = 1.82e-8; vectors with
  entries up to 1.7e14, 11-13 primes per family): VERIFIED (../../logs/exact4m/audit_verify4m_more.log).  Nine points.
  The repository copy of verify4m.py, reading facedata4mr.pkl.gz, gives the same output.

Reproduce (from this folder; ../relaxations must be present):
  python build_lin4m.py                                                                  (lin4m.pkl)
  python points4m.py k2 101 z_k2_a ; k2 102 z_k2_b ; k4 111 z_k4_a ; k4 112 z_k4_b ; inner 121 z_in_a ;
         inner 123 z_in_b ; inner 124 z_in_c ; inner 126 z_in_d ; sf 131 z_sf_a ; sf 132 z_sf_b ; sf 133 z_sf_c ;
         sf 134 z_sf_d ; sf 135 z_sf_e ; sf 136 z_sf_f ; sf 137 z_sf_g ; sf 138 z_sf_h ; ent 197 p_ent ;
         ent 198 p_ent2                                                                  (pts/*.pkl, one process each)
  python face4m.py z_k2_a z_k2_b z_k4_a z_k4_b z_in_a z_in_b z_in_c z_in_d z_sf_a z_sf_b z_sf_c z_sf_d z_sf_e z_sf_f
         z_sf_g z_sf_h                                                                   (face4m.pkl)
  python facedata4m.py ; python span4m.py ; python cert4m.py                            (cert4m.pkl; CERTIFIED line)
  python verify4m.py p_ent p_ent2 ; python crosscheck4m.py 25                           (independent re-checks)
  (the repository ships lin4m.pkl, face4m.pkl, facedata4mr.pkl.gz, cert4m.pkl and pts/*.pkl, so verify4m.py runs
  without the rebuild; crosscheck4m.py needs facedata4m.pkl from facedata4m.py)
  python negctl4m.py ; python verify4m.py --cert cert4m_neg_y.pkl p_ent ;
         python verify4m.py --cert cert4m_neg_S.pkl p_ent                               (both must FAIL)
  python face4m.py --nosave <subset of the zeros>                                       (saturation checks)

STATUS: n = 4, t = -1: CERTIFIED. Exact identity c = E^T y + sum_k F_k^*(S_k) with e . y = 0 (all 1331 coordinates and
the dual objective), c = (fu + fv - 2 fuv)/16 the objective of P = P_u - P_v, all 272 dual blocks PSD by exact LDL^T,
independent exact checks (verify4m.py) at two entangled true points passed, negative controls rejected, byte-identical
rerun.  Statement: q_4(P_u - P_v) >= 0 for all orthonormal u, v in (C^3)^(x)4 (complex pairs via realification).
Computer-assisted; not externally refereed.  Validity rests on the exact constraint data (same code as ../exact4,
validated there; face data cross-checked here).
