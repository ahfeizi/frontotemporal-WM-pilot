# E. Global null: count of cells with within-cell FWER hits

N = 22 subjects; 48 cells x 30 pairs; seed = 1; N_perm = 5000.
Per-subject contrasts computed from EEG files directly (not battery_persubject.csv).
Within-cell threshold: 95th pct of max|t| over 30 pairs under H0
(from same 5000 sign-flip patterns; compare exact ~3.67 in §4.2).

Observed: **3** cells pass FWER (cf. manuscript: 2 cells).
Cell-specific thresholds: min = 2.972, max = 3.274, median = 3.187.

## Distribution of count under joint null
  P(count = 0) = 0.1612
  P(count = 1) = 0.2246
  P(count = 2) = 0.1984
  P(count = 3) = 0.1700
  P(count = 4) = 0.1034
  P(count = 5) = 0.0696
**P(count >= 2) = 0.6142** (replaces "expectation of 2.4" in §4.2)

Outputs: E_global_null.csv (per-cell), E_global_null_dist.csv (distribution).
