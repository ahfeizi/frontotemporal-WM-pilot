# B. Prior sensitivity — JZS Cauchy scale r = 0.5, √2/2, 1.0

N = 22; cells = 48 (12 pipelines × 4 bands).

## Attainable ceiling BF₀₁ at t = 0, N = 22
- r = 0.5000: ceiling BF₀₁ = 3.3617
- r = 0.7071: ceiling BF₀₁ = 4.4854
- r = 1.0000: ceiling BF₀₁ = 6.1257

## Beta band (13–30 Hz) BF₀₁ range across 12 pipelines
- r = 0.5000: BF₀₁ min = 1.7309, max = 3.3612  (ceiling = 3.3617)
- r = 0.7071: BF₀₁ min = 2.1676, max = 4.4847  (ceiling = 4.4854)
- r = 1.0000: BF₀₁ min = 2.8304, max = 6.1247  (ceiling = 6.1257)

## Lowest-BF₀₁ pipeline in beta for each r
- r = 0.5000: dwPLI wavelet average
- r = 0.7071: dwPLI wavelet average
- r = 1.0000: dwPLI wavelet average

## Qualitative conclusion
  r = 0.5000: min BF₀₁ 1.7309 < 3 → at least one cell remains anecdotal.
  r = 0.7071: min BF₀₁ 2.1676 < 3 → at least one cell remains anecdotal.
  r = 1.0000: min BF₀₁ 2.8304 < 3 → at least one cell remains anecdotal.

Pipeline ordering preserved; qualitative conclusion stable: one outlier cell remains anecdotal (BF₀₁ < 3) regardless of prior scale.

Output: B_prior_sensitivity.csv