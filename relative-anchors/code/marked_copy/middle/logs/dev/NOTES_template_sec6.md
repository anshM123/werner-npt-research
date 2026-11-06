## 6. The families [q0 ≥ 4, t ≥ 4, i = b] for b = 2{{B3LIST}} (`run_extra.sh`{{B3SCRIPT}})

**Theorem N_b.** Let b ∈ {2{{B3SET}}}. For every k, j with q0 = j − b ≥ 4 and t = k − 1 − j − b ≥ 4, and every σ > 0, the
reduced form Q_b is ⪰ 0 with kernel exactly span(κ_0, κ_1) = K_j ∩ V_b.

*Proof (computer-assisted).* The certificate of Section 3 on the family with q0 = 4 + x, t = 4 + y symbolic and i = b,
exactly as for Theorem N1: all sixteen vectors present, all four blocks split, parts 12 + 4, gauge U_{1,0}, U_{1,1};
now all five Johnson nodes of the window enter (b + e ≥ 0 for e ≥ −2). The pivots, in the order taken, with the terms
of numerator / denominator:
- b = 2 (`logs/run_mid_i2.log`, 1 032 s): S_{10} 359/280, S_{00} 976/797, S_{01} 976/797, S_{11} 976/886; G_1 8/12,
  G_0 4/14, O_0 175/182, O_1 141/208, W_{11} 7 449/7 563, W_{10} 12 129/12 185, W_{01} 18 416/17 788,
  W_{00} 19 070/18 546, U_{00} 38 513/40 737, U_{01} 39 873/40 986.
{{B3PIVOTS}}
Every pivot has numerator and denominator with positive coefficients only, and a nonzero part free of x, y. ∎

**Corollary P_b.** For every k ≥ 2, every j = −1, …, k and every σ > 0, on the isotypic parts of all S_{k−1}-irreps
(k − 1 − i, i) and (k − 2 − i, i, 1) with i ≤ {{BMAX}}, the letter-1 block of (R∞)_j with squared weights (σ, 1, …, 1) is
⪰ 0, with kernel exactly the part of K_j there (and the same for every weight vector with all squared weights but one
equal).

*Proof.* As for Corollary P1, it suffices to treat Q_i for 2 ≤ i ≤ {{BMAX}}. If q0 ≥ 4 and t ≥ 4, Theorem N_i. If
q0 = a ≤ 3 and t ≥ 4, the member lies in the band family [q0 = a, t ≥ 4, i] of Theorem B (proved there for a + i ≤ 8,
and here a + i ≤ {{AIMAX}}); if q0 ≥ 4 and t = c ≤ 3, in the band family [q0 ≥ 4, t = c, i] (proved for c + i ≤ 7, and
here c + i ≤ {{AIMAX}}); if both are ≤ 3, k = q0 + t + 2i + 1 ≤ {{KSMALL}} and Theorem M applies.
`coverage_mid.py 80 32 band:{{BNEXT}}` checks the bookkeeping for k ≤ 80: every form not covered has 9 ≤ j ≤ k − 9,
i ≥ {{BNEXT}} and k ≥ 41. ∎

The band families stop covering the edges at i = 5 ([q0 ≥ 4, t = 3, i = 5] has c + i = 8), so beyond i = 4 the edge
families with i fixed would have to be added as well.
