# F. Pair-level hits table

Threshold: within-cell FWER (max-t sign-flip) p < 0.05.
Source: unified_contrasts_ds004117.csv, level==pair (1440 rows).
Pairs passing FWER in any cell: **2**

- beta_low | PLV | wavelet | csd | pair FCz-T7 | t = 3.9516 | p_FWER = 0.0153
- theta | dwPLI | hilbert | csd | pair F3-FT8 | t = 3.8619 | p_FWER = 0.0132

Columns in F_pair_hits.csv: hit_band, hit_measure, hit_decomp, hit_ref, hit_pair,
  hit_t, hit_p_fwer, other_cell (band__measure__decomp__ref), other_t, other_p.
One row for the hit cell itself (other_cell == hit cell) plus 47 rows (other cells),
per hit pair.
