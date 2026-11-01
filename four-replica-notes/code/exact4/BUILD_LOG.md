EXACT4: exact rational certificate for q_4(P) >= 0, P a rank-2 projector on (C^3)^(x)4 (t = +1), two-letter relaxation
=====================================================================================================================

Scope: scale the certified n = 3 pipeline (../relaxations/exact_pilot_n3.py; ../../NOTES.md, Section 3) to n = 4.
Relaxation (as at n = 2, 3): X_c = twirl((x)_r |c_r><c_r|), c in {uuuu, uuuv, uuvv, uvvv, vvvv}; X >= 0, moments,
swap identities (built into the variable space), ORTH, MARG, covariance blocks; objective m_uu + m_vv + 2 m_uv.
All code for this step lives in this folder (imports the exact per-site pieces from ../relaxations).

Log (times from the shell clock, EST)
-------------------------------------
S0 (2026-11-01 06:13) started from the n = 3 pipeline.

E1 (2026-11-01 06:28) exact engine for n sites written and validated.
  core4.py / engine4.py: exact variable space (RREF basis of the site-invariant symmetric matrices, so a coordinate is
  the entry at a pivot), per-site maps as integer matrices (denominator lcm 8, numerators <= 3 bits), and a shared
  tree over the sites that produces the images in ALL unsorted target patterns at once; the image of the site orbit in
  a sorted target is the sum over orbit representatives with the site axes permuted back.  Tree arithmetic is float64 on
  integers with an a-priori bound < 2^53 (asserted per group), hence exact.
  n = 4 variable space: 1331 parameters, 175 groups (largest K = 96 parameters, t <= 21, D <= 81), integer tensors
  <= 13 bits, sum K D^2 = 4.8M entries; builds in 4 s.
  test_engine.py: at a random rational point of the common variable space, our objective, all 120 moment functionals,
  the raw images of all ORTH/MARG maps (every sorted target) and every raw covariance block agree EXACTLY with the
  certified exact_coreN data: n = 2 (84 identity images, 30 covariance blocks) and n = 3 (140, 60), 0 mismatches.
  symm4.py: exact symmetry-adapted compressions of the covariance blocks: group T x Stab(mu) (T = {id, slot swap 01}
  for uu/vv, trivial for uv; same groups as the float code), rational primitive idempotents ((1 +- tau)/2 times Young
  symmetrisers), U = LLL-reduced integer basis of range(R(e)); 385 class blocks, largest 72 (as in the float run),
  isotypic dimension count checked for every block.  Any rational compression is a valid constraint; invariance only
  makes it lossless.

E2 (2026-11-01 06:44) exact linear data and exact points at n = 4.
  build_lin4.py (9 s): objective c and the 120 moment functionals (tree with scalar per-site maps), raw ORTH/MARG
  images for all sorted A3 targets (6186 integer rows after column scaling by dX); a maximal independent subset modulo
  p = 1048573 has 146 rows (= the floating-point identity rank 146), exact rank 146; exact RREF of [moments; identities]
  gives 147 independent equalities E x = e (entries <= 37 bits).  Saved lin4.pkl.
  PITFALL FOUND: the inherited modular moment code (exact_coreN.moments_modp) sums up to 3^5 = 243 products in one
  pairwise einsum at n = 4 (96 of the 13824 heads), so with 28-bit primes int64 can overflow silently (243 p^2 > 2^63).
  points4.py uses 27-bit primes (243 (2^27)^2 < 2^63) and asserts the bound on every einsum call.
  points4.py: exact moments (CRT, bound prod_r ||c_r||_1^2), exact integer Fourier transform, site symmetrisation,
  exact projection; about 4-6 min per point, run in parallel.  Zeros (all real, rational): k1 (seeds 11, 12), k3 (13, 14),
  rotated k1 (21, 22), rotated k3 (23, 24), rotated interpolating (31, 32, 33); entangled test point p_ent (seed 97:
  u = p, v = (p.p) q - (p.q) p, p, q in {-2..2}^81).  Each zero: in the variable space (exact reconstruction), moments
  exact, E x = e exact, objective exactly 0, and q_4(P) computed directly from P is exactly 0.  p_ent: in the variable
  space, E x = e exact, c . x(P) = q_4(P) = 0.998715701733... exactly (q_4 computed directly from the 81 x 81 matrix P).

E3 (2026-11-01 06:49) exact face and exact face-reduced constraint data.
  face4.py (25 s): x* = sum_z lambda_z x_z (lambda_z = denominator of x_z, so x* is an integer vector, <= 88 bits) over
  the 11 zeros; ker W(x*) = intersection of ker W(x_z) because every W(x_z) is PSD (true points).
  X blocks: Y_g(x*) exact (flint), dimension 709 -> 537.  Covariance class blocks: the numerators of W(x*) computed
  modulo 6 primes (27-bit) and recovered by CRT with an a-priori bound (|What| <= 2^158 < prod p / 2); dimension
  5028 -> 4067 (168 class blocks with W(x*) = 0); exact integer kernel bases (flint nullspace + LLL), entries <= 3 bits,
  W N = 0 checked exactly.  With 9 zeros (without z_rint_b, z_rint_c) the face is identical: saturated.
  facedata4.py (97 s): F_k(e_i) = N^T W_k(e_i) N for every parameter, computed as (Gt U N_r)^T Ytil_rs,i (U N_s)
  modulo 3 primes from the float64-exact tree images, CRT with the a-priori bound |Fhat| <= 2^64 (3 primes = 81 bits);
  numerators <= 45 bits.  Blocks whose face data vanish identically are dropped (vacuous): 158 covariance face blocks
  and 98 X face blocks remain (largest 65).  Exact face check: F_k(x_z) = 0 for all 256 blocks and all 11 zeros.

E4 (2026-11-01 06:53) first float solve, and exact removal of vacuous directions.
  First run of cert4.py on facedata4.pkl: float problem 31427 rows x 1331 (7 scalar + 249 PSD blocks), 40 null
  directions (as in the float run; objective component 1.3e-6 against ||c|| = 4.3e5, i.e. relative 3e-12), Clarabel
  Solved in 18 iterations, dual 1.2e-15, ||A^T z + c|| 3.3e-9, smallest reduced-dual eigenvalue 0.234.  The run then
  stopped in the rounding: 155 of the 256 face blocks have directions on which every parameter's data vanish (623 zero
  columns plus non-coordinate ones), and the Jacobi scaling blew up there (inf in the back-transformed dual).
  span4.py (26 s): for each block, B = LLL-reduced integer basis of the row space of all its data matrices (rows chosen
  independent mod 1048573, exact RREF); completeness checked EXACTLY (F_i Kb = 0 for an exact basis Kb of ker B, every
  data matrix); reduced data B F_i B^T (a valid compression; a dual S' of the reduced block is S = B^T S' B, and strict
  dual feasibility is preserved).  Total block dimension 2764 -> 1556, largest block 33, svec rows 31280 -> 9322,
  B entries <= 3 bits, numerators <= 46 bits.  Saved facedata4r.pkl.

E5 (2026-11-01 07:02) n = 4, t = +1: CERTIFIED (cert4.py on facedata4r.pkl, 40 s; independent re-check verify4.py).
  Float problem: 147 equality rows + 9322 svec rows (28 scalar + 228 PSD blocks, largest 33) x 1331 parameters; 40 null
  directions removed by whitening (objective component 1.2e-6 against ||c|| = 4.3e5).  Clarabel (chordal decomposition
  off): Solved in 16 iterations, primal -8.9e-12, dual -4.7e-16, ||A^T z + c|| 3.9e-9, smallest reduced-dual eigenvalue
  0.248 (Jacobi-normalised blocks); solve 21 s.
  Rounding: y and every S_k to dyadic grids (40 bits below each block's largest entry); exact residual 1.2e-4 in the
  identity (relative 3e-10), 2.7e-12 in e.y.  Column subset by pivoted QR on the 1332 x 9347 float equation matrix:
  rank 1291 of 1332 (diag(R) ratio 1.3e-3 at the cut, 2.2e-15 next), 55 multipliers y and 1236 dual entries;
  rows chosen by pivoted QR; exact 1291 x 1291 solve (flint, Dixon) in 2 s; max |delta| 6.0e-5.
  Afterwards ALL 1332 equations hold exactly: c = E^T y + sum_k F_k^*(S_k) and e . y = 0.
  Exact LDL^T: all 256 dual blocks PSD (smallest pivot 9.4e-5).  Certificate entries up to 1061/1061 bits (y up to
  1053/1055 bits); data entries: E <= 37 bits, face-block numerators <= 46 bits.
  Entangled true points (exact): p_ent (seed 97) q_4(P) = 0.998715701733..., p_ent2 (seed 98) q_4(P) = 0.938967775077...
  (q_4 recomputed in verify4.py directly from the 81 x 81 matrix P): in the variable space, E x(P) = e,
  c . x(P) = q_4(P), c . x(P) = y . E x(P) + sum_k Tr(S_k F_k(x(P))) with y . E x(P) = 0 and all 256 terms >= 0, and every
  face block F_k(x(P)) PSD by exact LDL^T.
  verify4.py (reads only lin4.pkl, facedata4r.pkl, cert4.pkl, pts/*.pkl): VERIFIED for both points.
  Cross-checks of the n = 4 data: crosscheck4.py recomputes the face numerators of 25 random covariance blocks with Python
  integers only (5.17M entries): all equal to the CRT values, and the reduced data equal B Fhat B^T.  test_engine.py 4:
  at a random rational point of the common variable space, the objective and the raw covariance blocks COVuu(0,0,1,2),
  COVuv(0,1,2,3) agree exactly with the certified exact_coreN code; the extended run (test_engine.py 4 more) also
  confirms all 120 moment functionals (0 mismatches; identity images and 8 more covariance targets running).
  Reproducibility: full rerun build_lin4 -> face4 -> facedata4 -> span4 -> cert4 -> verify4 gives the same numbers.
  Peak memory (psutil): build_lin4 1.1 GB private, face4 1.3 GB, facedata4 3.9 GB, span4 2.3 GB, cert4 3.7 GB,
  verify4 1.2 GB; points4 well below 1 GB each.  Wall: lin 10 s, points 4-8 min each (11 in parallel: about 8 min),
  face 24 s, facedata 95 s, span 26 s, cert 40 s, verify 8 s.
  Statement certified: q_4(P) >= 0 for every rank-2 orthogonal projector P on (C^3)^(x)4 (complex P included by the
  realification argument of the n = 2, 3 runs: the real part of the twirled two-letter moments of P is the average of
  the true points of P and conj(P), both with value q_4(P)).

E6 (2026-11-01 07:16) cross-checks completed.
  test_engine.py 4 more (1096 s): at a random rational point of the common variable space, our n = 4 data agree EXACTLY
  with the certified exact_coreN code for the objective, all 120 moment functionals, all 210 (ORTH/MARG map, sorted A3
  target) images, and 10 raw covariance blocks COVuu(0,0,1,2), (1,1,1,1), (0,1,3,3), (2,2,3,3), COVuv(0,1,2,3),
  (2,3,3,3), (1,1,2,2), (0,0,0,0), COVvv(3,3,3,3), (1,2,3,3) (targets with stabilisers 1, Z2, Z2 x Z2, S3, S4): 0 failures.
  Negative controls: verify4.py rejects a certificate with one multiplier changed (identity and e.y fail) and one with a
  dual diagonal entry negated (identity fails, block not PSD).

E7 (2026-11-01 07:36) independent audit of the saved certificate.
  verify4.py re-run at p_ent and p_ent2 (../../logs/exact4/audit_verify4.log), and at six further exact true points made
  afterwards with points4_audit.py from fresh seeds (../../logs/exact4/audit_pt_*.log, audit_verify4_more.log):
  a_ent1-4 (generic entangled pairs, q_4 = 0.940, 0.908, 1.006, 1.001) and two near-zeros, a_near1 (q_4 = 2.4e-7) and
  a_near4 (q_4 = 6.4e-5), each a zero scaled by 40 plus a small integer perturbation.  VERIFIED at all eight points.
  Two more near-zeros (a_near2, a_near3) stopped inside the point generator with "not enough primes" (its fixed prime
  list is too short for their integer sizes), so they were never checked; this is a capacity limit of the generator.
  The repository copy of verify4.py, reading facedata4r.pkl.gz, gives the same output.

Reproduce (from this folder; ../relaxations must be present):
  python test_engine.py 2 ; python test_engine.py 3 ; python test_engine.py 4 more     (engine vs exact_coreN)
  python build_lin4.py                                                                   (lin4.pkl)
  python points4.py k1 11 z_k1_a ; k1 12 z_k1_b ; k3 13 z_k3_a ; k3 14 z_k3_b ; rk1 21 z_rk1_a ; rk1 22 z_rk1_b ;
         rk3 23 z_rk3_a ; rk3 24 z_rk3_b ; rint 31 z_rint_a ; rint 32 z_rint_b ; rint 33 z_rint_c ; ent 97 p_ent ;
         ent 98 p_ent2                                                                   (pts/*.pkl, one process each)
  python face4.py z_k1_a z_k1_b z_k3_a z_k3_b z_rk1_a z_rk1_b z_rk3_a z_rk3_b z_rint_a z_rint_b z_rint_c   (face4.pkl)
  python facedata4.py ; python span4.py ; python cert4.py                               (cert4.pkl; CERTIFIED line)
  python verify4.py p_ent p_ent2 ; python crosscheck4.py 25                             (independent re-checks)
  (the repository ships lin4.pkl, face4.pkl, facedata4r.pkl.gz, cert4.pkl and pts/*.pkl, so verify4.py runs without
  the rebuild; crosscheck4.py needs facedata4.pkl from facedata4.py)

STATUS: n = 4, t = +1: CERTIFIED. Exact identity c = E^T y + sum_k F_k^*(S_k) with e . y = 0 (all 1331 coordinates and
the dual objective), all 256 dual blocks PSD by exact LDL^T, independent exact checks at two entangled true points passed.
Computer-assisted; not externally refereed.  Validity rests on the exact constraint data (validated as above).
