#!/usr/bin/env python3
"""Sensitivity (minimum detectable effect) analysis at the fixed sample sizes.

Both N's were fixed by the original public datasets (D1 ds003702, N = 47;
D2 ds004117, N = 22), so this reports what the paired design could detect,
not a sample-size calculation. Analytic, two-sided one-sample t test on the
within-participant differences (the primary test of Section 2.5), effect in d_z.

    python brm_sensitivity_power.py      -> brm_sensitivity_power.csv / .md

Alpha levels: 0.05 (single cell) and the first Holm step of each dataset's
region-level family (0.05/24 for D1, 0.05/48 for D2), which is the most
stringent threshold any cell in that family has to pass.
TOST: probability of declaring equivalence at the manuscript's bounds
(d_z = +/-0.4, alpha 0.05 per one-sided test) when the true d_z is 0.
"""
from pathlib import Path
import numpy as np
from scipy import stats, optimize

OUT = Path(__file__).parent
DATASETS = {"D1 ds003702": (47, 24), "D2 ds004117": (22, 48)}
POWERS = (0.80, 0.90, 0.95)
EFFECTS = (0.2, 0.4, 0.5, 0.8)
BOUND = 0.4


def power(dz, n, alpha):
    df, nc = n - 1, dz * np.sqrt(n)
    tc = stats.t.ppf(1 - alpha / 2, df)
    # lower tail via symmetry, P(T < -tc | nc) = P(T > tc | -nc); scipy's nct.cdf
    # returns NaN deep in the lower tail at large nc, where the true value is ~0
    lower = np.nan_to_num(stats.nct.sf(tc, df, -nc), nan=0.0)
    return stats.nct.sf(tc, df, nc) + lower


def mde(n, alpha, target):
    return optimize.brentq(lambda d: power(d, n, alpha) - target, 1e-6, 2.5, xtol=1e-10)


def tost_power_at_zero(n, bound=BOUND, alpha=0.05):
    # Bounds are +/-bound x sample SD (Section 2.5: one-sided t tests on the mean
    # difference against the scaled bounds). Both one-sided tests reject iff
    # |t_obs| < sqrt(n)*bound - t_crit, and under d_z = 0 t_obs is central t(n-1),
    # so the rate is exact. Sample-SD scaling matches util_contrast_stats.tost_p
    # ("bounds expressed in sample-SD units"); seeded Monte Carlo agrees (0.122 at N = 22).
    df = n - 1
    tc = stats.t.ppf(1 - alpha, df)
    lim = np.sqrt(n) * bound - tc
    if lim <= 0:
        return 0.0
    return stats.t.cdf(lim, df) - stats.t.cdf(-lim, df)


rows = []
for name, (n, fam) in DATASETS.items():
    for alpha, lab in ((0.05, "0.05"), (0.05 / fam, f"0.05/{fam} (Holm step 1)")):
        for p in POWERS:
            rows.append((name, n, lab, "MDE d_z", f"power {p:.2f}", round(mde(n, alpha, p), 3)))
        for e in EFFECTS:
            rows.append((name, n, lab, "power", f"d_z {e}", round(power(e, n, alpha), 3)))
    rows.append((name, n, "0.05 per side", "TOST equivalence rate at true d_z = 0",
                 f"bounds +/-{BOUND}", round(tost_power_at_zero(n), 3)))

hdr = "dataset,N,alpha,quantity,condition,value"
(OUT / "brm_sensitivity_power.csv").write_text(
    hdr + "\n" + "\n".join(",".join(map(str, r)) for r in rows) + "\n", encoding="utf-8")

md = ["# Sensitivity analysis at the fixed N (two-sided paired t, d_z)", "",
      "| Dataset | N | alpha | MDE at 80% / 90% / 95% power | Power at d_z 0.2 / 0.4 / 0.5 / 0.8 |",
      "|---|---|---|---|---|"]
for name, (n, fam) in DATASETS.items():
    for lab in ("0.05", f"0.05/{fam} (Holm step 1)"):
        m = [r[5] for r in rows if r[0] == name and r[2] == lab and r[3] == "MDE d_z"]
        pw = [r[5] for r in rows if r[0] == name and r[2] == lab and r[3] == "power"]
        md.append(f"| {name} | {n} | {lab} | {' / '.join(map(str, m))} | {' / '.join(map(str, pw))} |")
for name, (n, _) in DATASETS.items():
    v = [r[5] for r in rows if r[0] == name and r[3].startswith("TOST")][0]
    md.append(f"- {name}: TOST (+/-{BOUND}) declares equivalence with probability {v} when true d_z = 0.")
md.append("- Analytic (scipy noncentral t); no random numbers, so no seed.")
(OUT / "brm_sensitivity_power.md").write_text("\n".join(md) + "\n", encoding="utf-8")
print("\n".join(md))
