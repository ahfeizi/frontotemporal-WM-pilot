# Sensitivity analysis at the fixed N (two-sided paired t, d_z)

| Dataset | N | alpha | MDE at 80% / 90% / 95% power | Power at d_z 0.2 / 0.4 / 0.5 / 0.8 |
|---|---|---|---|---|
| D1 ds003702 | 47 | 0.05 | 0.417 / 0.483 / 0.537 | 0.269 / 0.766 / 0.919 / 1.0 |
| D1 ds003702 | 47 | 0.05/24 (Holm step 1) | 0.603 / 0.671 / 0.727 | 0.038 / 0.317 / 0.569 / 0.983 |
| D2 ds004117 | 22 | 0.05 | 0.626 / 0.725 / 0.807 | 0.146 / 0.433 / 0.609 / 0.947 |
| D2 ds004117 | 22 | 0.05/48 (Holm step 1) | 1.009 / 1.118 / 1.208 | 0.007 / 0.051 / 0.111 / 0.5 |
- D1 ds003702: TOST (+/-0.4) declares equivalence with probability 0.707 when true d_z = 0.
- D2 ds004117: TOST (+/-0.4) declares equivalence with probability 0.122 when true d_z = 0.
- Analytic (scipy noncentral t); no random numbers, so no seed.
