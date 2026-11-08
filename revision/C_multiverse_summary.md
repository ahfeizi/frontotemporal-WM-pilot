# C. Multiverse summary (group-level Python portion)

N = 22; 48 cells (12 pipelines × 4 bands); BF₀₁ ceiling = 4.4854.
Sampling SE of dz (approx, under H0) = 1/√N = 0.2132.

## Holm correction across 48 cells
Significant cells (Holm p < 0.05): 0

## SD of dz across 12 pipelines, by band
- theta     : SD = 0.0961  (SD/SE = 0.45)
- beta      : SD = 0.1185  (SD/SE = 0.56)
- beta_low  : SD = 0.2201  (SD/SE = 1.03)
- beta_high : SD = 0.0621  (SD/SE = 0.29)

## Main effects: mean dz by factor level, beta band
  measure: ImC=-0.1095, PLV=0.0020, dwPLI=0.0630
  decomposition: hilbert=-0.0643, wavelet=0.0345
  reference: average=0.0372, csd=-0.0670

Bootstrap dz_range (beta, 2000 resamples, seed 1): median = 0.8393, 95% CI [0.3855, 1.4216]
Bootstrap bf01_range (beta, 2000 resamples, seed 1): median = 3.8413, 95% CI [1.9093, 4.4498]

Outputs: C_forest_data.csv, C_variance_decomp.csv