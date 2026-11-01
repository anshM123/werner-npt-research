EXACT4L: exact rational certificate for q_4(a P_u + b P_v) >= 0 for all real a, b and orthonormal u, v in (C^3)^(x)4
=====================================================================================================================

Scope: combine the two certified n = 4 pipelines (../exact4, t = +1; ../exact4m, t = -1) into one certificate for
every ratio (a, b), i.e. for every Hermitian operator of rank <= 2 on (C^3)^(x)4.  Same relaxation, same variable space
(1331 parameters), same 147 equalities E x = e, same blocks (X groups, COV class blocks) and the same engine
(core4.py, engine4.py, symm4.py copied unchanged, md5 identical).  Nothing in ../exact4, ../exact4m or ../relaxations was
modified; their saved files are only read.

Lifted problem in (+, -) coordinates (float reference: ../relaxations/lift_fr3.py, ../../logs/lift_fr3_n4_cov.log):
  variables Y+, Y-, Z in R^1331, eps+, eps-, zeta;
  E Y+ = e eps+,  E Y- = e eps-,  E Z = e zeta,  [[eps+, zeta], [zeta, eps-]] >= 0,  (eps+ + eps-)/2 = 1;
  block k:  L_k = [[F+_k(Y+), G_k(Z)], [G_k(Z)^T, F-_k(Y-)]] >= 0,  F+-_k(x) = N+-_k^T W_k(x) N+-_k,  G_k(x) = N+_k^T W_k(x) N-_k;
  objective (1/4) (c+ . Y+ + c- . Y- + 2 d . Z),  c+- = (fu + fv +- 2 fuv)/16,  d = (fu - fv)/16.
Lifted true point for a pair (u, v) and a ratio (a, b): Y+ = (a+b)^2 x, Y- = (a-b)^2 x, Z = (a^2 - b^2) x,
eps+ = (a+b)^2, eps- = (a-b)^2, zeta = a^2 - b^2 (x = x(u, v)); then L_k = M^T W_k(x) M with M = [(a+b) N+_k, (a-b) N-_k],
so L_k >= 0, and the objective is
  (1/4) [(a+b)^2 (fu + fv + 2 fuv) + (a-b)^2 (fu + fv - 2 fuv) + 2 (a^2 - b^2)(fu - fv)] . x / 16
  = (1/4) [4 a^2 fu + 4 b^2 fv + 8 a b fuv] . x / 16 = (a^2 fu + b^2 fv + 2 a b fuv) . x / 16 = q_4(a P_u + b P_v).
Certificate: y+, y-, yZ in Q^147 and S_k = [[A_k, B_k], [B_k^T, C_k]] >= 0 with, exactly,
  c+/4 = E^T y+ + sum_k F+_k^*(A_k),  c-/4 = E^T y- + sum_k F-_k^*(C_k),  d/2 = E^T yZ + sum_k 2 G_k^*(B_k),
  e . y+ = e . y- = e . yZ = 0      (G_k^*(B)_i = Tr(B^T G_k(e_i))).
Then objective = y+ . E Y+ + y- . E Y- + yZ . E Z + sum_k Tr(S_k L_k) = sum_k Tr(S_k L_k) >= 0 on the lifted relaxation
(y . E Y = (e . y) eps = 0), in particular q_4(a P_u + b P_v) >= 0 at every lifted true point.

Log (times from the shell clock, EST)
-------------------------------------
S0 (2026-11-01 09:09) started; read EXACT4_NOTES.md, EXACT4M_NOTES.md, both pipelines and ../relaxations/lift_fr3.py with its
  log (float lifted problem with the float code's blocks: Solved, smallest lifted dual eigenvalue 9.24e-4).
  Sizes of the two exact faces (face4.pkl / face4m.pkl): COV class blocks 385, of which 376 have both faces nonempty,
  3 only N+, 3 only N-, 3 none; sum of d+ + d- over COV blocks 8186 (largest 129); X groups 150, sum d+ 537, sum d- 493
  (largest d+ + d- = 42).  The class compressions U are identical in both face files.

E1 (2026-11-01 09:18) exact linear data (build_linL.py, 15 s, 1.14 GB private).
  fu, fv, fuv by functional_row exactly as in ../exact4/build_lin4.py (fu = Tr[X_uuuu W01], fv = Tr[X_vvvv W01],
  fuv = Tr[X_uuvv W12]); c+, c-, d; E, e recomputed (147 independent equalities, entries <= 37 bits).
  Exact comparison: c+ equals the certified objective of ../exact4/lin4.pkl; c- equals the certified objective of
  ../exact4m/lin4m.pkl (and its saved cplus equals c+); E, e and the 120 moment rows are identical in both files;
  c+ - c- = fuv/4.  d = (fu - fv)/16 has 90 nonzero coordinates.  Saved linL.pkl.

E2 (2026-11-01 09:28) exact cross data (facedataL.py, 486 s, 7.02 GB private).
  For every block the full compression N^T W N with N = [N+ | N-] is computed, so the F+ and F- parts double as a
  check of the new code against the certified face data.  COV blocks: numerator Hhat = H0 + H0^T,
  H0 = sum_{r,s} (Gt U N_r)^T Ytil_rs (U N_s) (W = W0 + W0^T as in facedata4.py, so the cross part is
  N+^T W0 N- + (N-^T W0 N+)^T), modulo 3 primes of 27 bits with the a-priori CRT bound |Hhat| <= 2^65 (bound recomputed
  for the joint bases); three passes xy = uu (152 blocks), uv (78 blocks, largest d+ + d- = 129), vv (152) to bound the
  memory.  X blocks: N^T Y_k N exactly (flint).
  Exact comparison: the F+ parts equal ../exact4/facedata4.pkl in all 7745 block-groups and the F- parts equal
  ../exact4m/facedata4m.pkl in all 8036 block-groups (0 mismatches; blocks / groups absent there vanish here).
  Cross data: 490 blocks with d+ + d- > 0 (382 COV, 108 X), 242 of them with G not identically zero; 57,952,701
  numerator entries (COV numerators <= 46 bits).
  Which zeros kill what (worked out, then tested exactly).  At a t = +1 zero, W(x_z) >= 0 and N+ lies in its kernel,
  so W(x_z) N+ = 0; at a t = -1 zero, W(x_z) N- = 0.  Hence G(x_z) = N+^T W(x_z) N- = 0 at BOTH kinds of zeros, while
  F+(x_z) = 0 only at t = +1 zeros and F-(x_z) = 0 only at t = -1 zeros.  Exact test: G_k(x_z) = 0 for all 490
  blocks at all 27 zeros, F+_k(x_z) = 0 at the 11 t = +1 zeros, F-_k(x_z) = 0 at the 16 t = -1 zeros (0 failures);
  and, as expected, F+_k(x_z) != 0 at some t = -1 zero for 139 of the 256 t = +1 face blocks, F-_k(x_z) != 0 at some
  t = +1 zero for 170 of the 272 t = -1 face blocks.  (So the lifted point of a t = +1 zero at ratio (1, 1), and of a
  t = -1 zero at ratio (1, -1), has every lifted block equal to 0.)  Saved facedataL.pkl (169 MB).

E3 (2026-11-01 09:33) exact lifted span reduction (spanL.py, 246 s, 4.97 GB private).
  V+_k = row space of the stack [F+_k(e_i); G_k(e_i)^T] (all i), V-_k = row space of [F-_k(e_i); G_k(e_i)] (rows chosen
  independent mod 1048573, exact RREF, LLL); completeness checked EXACTLY (every stacked matrix times an exact basis of
  ker B+- is 0); reduced data B+ F+ B+^T, B- F- B-^T, B+ G B-^T (common denominator per block and group).
  Result: total dimension 9216 -> 3120; 286 lifted blocks kept (largest 65, svec 35595, 242 with both parts);
  B entries <= 3 bits, numerators <= 48 bits.  The separate spans of F+ alone and F- alone, recomputed here, have total
  dimension 1556 and 1564 (= span4.py and span4m.py), and 1556 + 1564 = 3120: in every block the column spaces of the
  cross data G_k(e_i) already lie in the F+ span and their row spaces in the F- span, so the lifted blocks live on exactly
  the reduced spaces of the two certified problems.  Saved facedataLr.pkl.

E4 (2026-11-01 09:34 - 09:37) first float solves FAILED numerically (diagnosed, logs kept in logs/).
  Float lifted problem (certL.py): 36040 rows x 3996 columns (Y+, Y-, Z, eps+, eps-, zeta), nnz 14.4M; zero cone 442
  rows (3 x 147 equalities + normalisation), the 2 x 2 eps-block, 22 scalar and 264 PSD lifted blocks (largest 65,
  sum svec^2 20.7M); blocks Jacobi-normalised.  Gram matrix of the columns: rank 3893 of 3996 (103 null directions, as
  in the float reference), objective component on the null space 3.8e-15 (relative).  Column norms 1.6e-3 .. 19,
  scaled Gram eigenvalues (non-null) 1.2e-6 .. 1, |c_j| / column norm up to 1.4e5: the RREF coordinates are badly
  scaled for an interior-point method.
    (i) independent columns as they are (as in lift_fr3.py), faer, chordal off: NumericalError at iteration 4
        (logs/certL_plain_try1.log);
    (ii) same with columns scaled to unit norm: NumericalError at iteration 1 (logs/certL_colscale_try2.log).
  Fix: whitening, as ../exact4/cert4.py does, but per column group (Y+, Y-, Z separately, ranks 1291 / 1298 / 1301;
  eps columns only scaled), which keeps the block rows sparse outside their own group (46.2M nonzeros).

E5 (2026-11-01 09:45) n = 4, lifted (all ratios): CERTIFIED (certL.py, default --mode whiten; 452 s, 10.16 GB private).
  Float solve: Clarabel (faer, 24 threads, chordal decomposition off, tolerances 1e-10): Solved in 19 iterations
  (289 s), primal +4.9e-10, dual +5.6e-10, ||A^T z + c|| 8.0e-8; eps-block dual [0, 0, 0], normalisation multiplier
  -5.6e-10.  Polish (minimum-norm correction of the dual residual over all columns, pseudo-inverse of the scaled Gram
  matrix, twice): ||A^T z + c|| 3.5e-10 (against |c| up to ~1e5).  Smallest lifted dual eigenvalue (Jacobi-normalised
  blocks) 9.607e-4 at Xuuvv(3,3,3,3) (r+ 21, r- 21); no block below 1e-6; smallest ++ part 6.6e-2, -- part 4.3e-2.
  (Float reference lift_fr3.py with its own blocks: 9.24e-4.)
  Rounding to dyadic grids (40 bits below each block's / each y-vector's largest entry): exact residuals 2.1e-5 (+),
  1.4e-5 (-), 5.5e-6 (Z) in the identities, |e . y| <= 8.8e-13.
  Exact correction, separately for the three decoupled identities (1332 equations each; columns chosen by pivoted QR
  on the float equation matrix, rows by pivoted QR, exact square solve with flint Dixon):
    +: 9347 columns, rank 1291 (diag R ratio 1.3e-3 at the cut, 2.3e-15 next), 55 y's + 1236 dual entries, 2 s, max |delta| 1.0e-5;
    -: 9663 columns, rank 1298 (1.3e-3 / 2.1e-15), 57 y's + 1241 dual entries, 1 s, max |delta| 6.9e-6;
    Z: 16633 columns, rank 1301 (1.6e-3 / 1.6e-15), 54 y's + 1247 cross entries B, 2 s, max |delta| 3.3e-6.
  (The equation ranks of + and - are those of exact4 / exact4m, 1291 and 1298.)  Afterwards ALL 3 x 1332 equations
  hold exactly: c+/4 = E^T y+ + sum_k F+_k^*(A_k), c-/4 = E^T y- + sum_k F-_k^*(C_k), d/2 = E^T yZ + sum_k 2 G_k^*(B_k),
  e . y+ = e . y- = e . yZ = 0.
  Exact LDL^T: all 286 lifted dual blocks S_k = [[A_k, B_k], [B_k^T, C_k]] PSD, no zero pivot, smallest pivot 1.11e-5
  (LDL 98 s).  Certificate entries up to 1054 / 1056 bits (numerator / denominator).
  In-pipeline check at ../exact4/pts/p_ent (exact): E x = e, c+ . x = stored q_4(P_u + P_v); for the ratios (1,0),
  (1,1), (1,-1), (2,1) the lifted objective equals y-terms (= 0) + sum_k Tr(S_k L_k) with all 286 terms >= 0.
  Saved certL.pkl (y+, y-, yZ and S_k as (num, den)) and float_dualL.pkl; md5 in logs/md5_first_run.txt.
