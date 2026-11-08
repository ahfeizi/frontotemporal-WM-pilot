# G. FCz beta power: exact vs parametric p discrepancy

## Connectivity cells with largest |exact p − parametric p| (top 5)
- beta_high ImC wavelet average: t=-0.6632, p_param=0.5144, p_exact=0.5569, divergence=0.0425
- beta_low PLV wavelet csd: t=-0.1774, p_param=0.8609, p_exact=0.8933, divergence=0.0324
- beta_high ImC hilbert average: t=-0.4591, p_param=0.6509, p_exact=0.6709, divergence=0.0200
- theta dwPLI hilbert average: t=1.1287, p_param=0.2718, p_exact=0.2877, divergence=0.0159
- theta ImC hilbert average: t=0.9960, p_param=0.3306, p_exact=0.3453, divergence=0.0147

## FCz beta power (manuscript §4.1: t=-2.23, exact p=0.017, param p=0.037)
N included = 22
Computed: dz = -0.4651, t = -2.1814, p_param = 0.0407
Skewness = -0.2502
Shapiro-Wilk W = 0.9498, p = 0.3130 (reject normality: 0)
LOO dz range = [-0.5531, -0.4042]
Without sub-024: dz = -0.4174, t = -1.9126

Output: G_cell_discrepancy.csv (22 rows: participant, contrast, loo_dz)
