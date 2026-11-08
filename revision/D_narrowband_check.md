# D. Narrowband check

Fractional bandwidth = (f_high − f_low) / f_center.

## Fractional bandwidth by band
- theta       : 4–7 Hz, f_c = 5.5 Hz, fractional BW = 0.545
- beta        : 13–30 Hz, f_c = 21.5 Hz, fractional BW = 0.791
- beta_low    : 13–20 Hz, f_c = 16.5 Hz, fractional BW = 0.424
- beta_high   : 21–30 Hz, f_c = 25.5 Hz, fractional BW = 0.353

## |dz(Hilbert) − dz(wavelet)| per band (mean ± max across 6 measure×ref cells)
- theta       : mean = 0.0959,  max = 0.2005
- beta        : mean = 0.1143,  max = 0.2573
- beta_low    : mean = 0.1458,  max = 0.3566
- beta_high   : mean = 0.0228,  max = 0.0435

## Full table (n=4 bands; see D_narrowband_check.csv)
band         frac_BW  measure ref      dz_H     dz_W     |diff| 
theta           0.545 PLV     average  +0.4174  +0.3168   0.1006
theta           0.545 PLV     csd      +0.2982  +0.4987   0.2005
theta           0.545 dwPLI   average  +0.2406  +0.3344   0.0937
theta           0.545 dwPLI   csd      +0.3762  +0.4437   0.0674
theta           0.545 ImC     average  +0.2123  +0.1999   0.0125
theta           0.545 ImC     csd      +0.3267  +0.2261   0.1006
beta            0.791 PLV     average  -0.0039  +0.0935   0.0974
beta            0.791 PLV     csd      -0.0638  -0.0179   0.0459
beta            0.791 dwPLI   average  +0.0171  +0.2744   0.2573
beta            0.791 dwPLI   csd      -0.0895  +0.0498   0.1393
beta            0.791 ImC     average  -0.1286  -0.0291   0.0995
beta            0.791 ImC     csd      -0.1170  -0.1635   0.0465
beta_low        0.424 PLV     average  +0.0026  +0.1670   0.1644
beta_low        0.424 PLV     csd      -0.3352  -0.0378   0.2974
beta_low        0.424 dwPLI   average  +0.3980  +0.3935   0.0045
beta_low        0.424 dwPLI   csd      -0.1562  +0.2004   0.3566
beta_low        0.424 ImC     average  +0.0850  +0.0822   0.0028
beta_low        0.424 ImC     csd      -0.0887  -0.1381   0.0494
beta_high       0.353 PLV     average  +0.0196  +0.0051   0.0145
beta_high       0.353 PLV     csd      +0.0078  +0.0032   0.0045
beta_high       0.353 dwPLI   average  -0.0406  -0.0381   0.0026
beta_high       0.353 dwPLI   csd      -0.1349  -0.1027   0.0322
beta_high       0.353 ImC     average  -0.0979  -0.1414   0.0435
beta_high       0.353 ImC     csd      -0.0957  -0.1352   0.0395

Output: D_narrowband_check.csv