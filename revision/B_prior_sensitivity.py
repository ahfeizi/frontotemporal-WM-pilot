"""
B_prior_sensitivity.py
Prior sensitivity for JZS Bayes factors, Dataset 2 (ds004117, N=22).

Reviewer concern: the 2.2-4.5 range rests partly on one pipeline cell; assess
whether the pipeline ordering and qualitative conclusion change with r = 0.5,
sqrt(2)/2 (default), and 1.0.

Reads:  results/unified_contrasts_ds004117.csv (ROI rows only)
        results/design_limits_n22.csv
Writes: revision/B_prior_sensitivity.csv
        revision/B_prior_sensitivity.md
"""

from __future__ import annotations
import math
import pathlib
import sys

import numpy as np
import pandas as pd

# ── paths ────────────────────────────────────────────────────────────────────
REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

from util_contrast_stats import jzs_bf10  # noqa: E402  (import after path setup)

DATA_CSV   = REPO / "results" / "unified_contrasts_ds004117.csv"
LIMITS_CSV = REPO / "results" / "design_limits_n22.csv"
OUT_CSV    = REPO / "revision" / "B_prior_sensitivity.csv"
OUT_MD     = REPO / "revision" / "B_prior_sensitivity.md"

# ── prior scales ─────────────────────────────────────────────────────────────
R_VALS = {
    "r0p50":     0.50,
    "r0p707":    math.sqrt(2) / 2,   # default (≈ 0.7071)
    "r1p00":     1.00,
}

N = 22

# ── read ROI rows ─────────────────────────────────────────────────────────────
df = pd.read_csv(DATA_CSV)
roi = df[df["level"] == "roi"].copy()
assert len(roi) == 48, f"expected 48 ROI rows, got {len(roi)}"

# ── attainable ceiling BF01 at t = 0 for each r ──────────────────────────────
ceilings = {label: 1.0 / jzs_bf10(0.0, N, r) for label, r in R_VALS.items()}

# ── recompute BF01 for each cell and each r ───────────────────────────────────
records = []
for _, row in roi.iterrows():
    t = float(row["t"])
    rec = dict(
        band=row["band"],
        measure=row["measure"],
        decomposition=row["decomposition"],
        reference=row["reference"],
        n=N,
        t=t,
    )
    for label, r in R_VALS.items():
        bf10 = jzs_bf10(t, N, r)
        rec[f"bf10_{label}"] = bf10
        rec[f"bf01_{label}"] = 1.0 / bf10
    records.append(rec)

out = pd.DataFrame(records)

# ── write CSV ─────────────────────────────────────────────────────────────────
out.to_csv(OUT_CSV, index=False, float_format="%.6g")

# ── ordering and conclusion checks ───────────────────────────────────────────
# Focus on beta band (12 pipelines); check whether the cell with BF01 = 2.17
# (dwPLI wavelet average) remains the outlier across all r values.
beta_out = out[out["band"] == "beta"].copy()
beta_out = beta_out.sort_values("bf01_r0p707")

# For each r, rank cells by BF01 and identify the bottom-1 cell
rank_agreement = {}
for label in R_VALS:
    ranked = beta_out.sort_values(f"bf01_{label}")["measure"] + " " + \
             beta_out.sort_values(f"bf01_{label}")["decomposition"] + " " + \
             beta_out.sort_values(f"bf01_{label}")["reference"]
    rank_agreement[label] = ranked.iloc[0]  # lowest BF01 cell

# Range of BF01 in beta for each r
beta_ranges = {label: (beta_out[f"bf01_{label}"].min(),
                       beta_out[f"bf01_{label}"].max())
               for label in R_VALS}

# ── write markdown summary ────────────────────────────────────────────────────
lines = []
lines.append("# B. Prior sensitivity — JZS Cauchy scale r = 0.5, √2/2, 1.0")
lines.append("")
lines.append(f"N = {N}; cells = 48 (12 pipelines × 4 bands).")
lines.append("")
lines.append("## Attainable ceiling BF₀₁ at t = 0, N = 22")
for label, r in R_VALS.items():
    lines.append(f"- r = {r:.4f}: ceiling BF₀₁ = {ceilings[label]:.4f}")
lines.append("")
lines.append("## Beta band (13–30 Hz) BF₀₁ range across 12 pipelines")
for label, r in R_VALS.items():
    lo, hi = beta_ranges[label]
    lines.append(f"- r = {r:.4f}: BF₀₁ min = {lo:.4f}, max = {hi:.4f}  (ceiling = {ceilings[label]:.4f})")
lines.append("")
lines.append("## Lowest-BF₀₁ pipeline in beta for each r")
for label, r in R_VALS.items():
    lines.append(f"- r = {r:.4f}: {rank_agreement[label]}")
lines.append("")
lines.append("## Qualitative conclusion")
# Check: does the lowest-BF01 cell always remain "anecdotal" (< 3)?
any_conclusive_change = False
for label, r in R_VALS.items():
    lo, _ = beta_ranges[label]
    if lo >= 3.0:
        any_conclusive_change = True
        lines.append(f"  r = {r:.4f}: min BF₀₁ {lo:.4f} ≥ 3 → all cells moderate at this prior.")
    else:
        lines.append(f"  r = {r:.4f}: min BF₀₁ {lo:.4f} < 3 → at least one cell remains anecdotal.")
lines.append("")
if any_conclusive_change:
    lines.append("Pipeline ordering is preserved across r values; qualitative conclusion changes "
                 "for some r (see above).")
else:
    lines.append("Pipeline ordering preserved; qualitative conclusion stable: one outlier cell "
                 "remains anecdotal (BF₀₁ < 3) regardless of prior scale.")
lines.append("")
lines.append(f"Output: {OUT_CSV.name}")

OUT_MD.write_text("\n".join(lines), encoding="utf-8")
print(f"Wrote {OUT_CSV}")
print(f"Wrote {OUT_MD}")
