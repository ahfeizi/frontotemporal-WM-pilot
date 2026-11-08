"""
C_multiverse_summary.py  — group-level portion (Python on existing CSVs)
Dataset 2 (ds004117, N = 22).

Produces:
  (1) Forest-plot data: per-cell dz, 95% CI, BF01, BF01 ceiling, Holm-corrected p, Holm sig.
  (2) Variance of dz across the 12 pipelines within each band, compared with sampling SE.
      Main effects: mean dz by measure, decomposition, and reference.

Note: participant bootstrap of dz range and BF01 range in beta requires per-participant
ROI values.  If revision/roi_per_participant.csv does not exist, this script prints a
warning and skips the bootstrap (the companion MATLAB script saves that file).

Reads:  results/unified_contrasts_ds004117.csv (ROI rows only)
        results/design_limits_n22.csv
Writes: revision/C_forest_data.csv
        revision/C_variance_decomp.csv
        revision/C_bootstrap_beta.csv   (if per-participant file exists)
        revision/C_multiverse_summary.md
"""

from __future__ import annotations
import math
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))
from util_contrast_stats import jzs_bf10  # noqa: E402

DATA_CSV       = REPO / "results" / "unified_contrasts_ds004117.csv"
LIMITS_CSV     = REPO / "results" / "design_limits_n22.csv"
PER_PART_CSV   = REPO / "revision" / "roi_per_participant.csv"
OUT_FOREST     = REPO / "revision" / "C_forest_data.csv"
OUT_VARIANCE   = REPO / "revision" / "C_variance_decomp.csv"
OUT_BOOTSTRAP  = REPO / "revision" / "C_bootstrap_beta.csv"
OUT_MD         = REPO / "revision" / "C_multiverse_summary.md"

N             = 22
BF01_CEILING  = 1.0 / jzs_bf10(0.0, N)           # 4.485 at r = sqrt(2)/2
SE_DZ         = 1.0 / math.sqrt(N)                 # ≈ 0.2132 (approx SE of dz under H0)
BOOTSTRAP_SEED = 1
N_BOOT        = 2_000
BANDS_ORDER   = ["theta", "beta", "beta_low", "beta_high"]

# ── load ROI rows ─────────────────────────────────────────────────────────────
df = pd.read_csv(DATA_CSV)
roi = df[df["level"] == "roi"].copy()
assert len(roi) == 48

# ── (1) Forest-plot data with Holm correction ─────────────────────────────────
p_vals = roi["p"].values.astype(float)          # 48 parametric p values
# Holm step-down
order = np.argsort(p_vals)
holm_p = np.empty(len(p_vals))
cum_max = 0.0
for rank, idx in enumerate(order):
    corrected = p_vals[idx] * (len(p_vals) - rank)
    cum_max = max(cum_max, corrected)
    holm_p[idx] = min(cum_max, 1.0)

forest = roi[["band", "measure", "decomposition", "reference",
              "dz", "dz_ci_lo", "dz_ci_hi", "t", "p",
              "p_perm", "bf10", "bf01", "p_tost", "n"]].copy()
forest["bf01_ceiling"] = round(BF01_CEILING, 4)
forest["p_holm"] = holm_p
forest["holm_sig"] = (holm_p < 0.05).astype(int)
forest.to_csv(OUT_FOREST, index=False, float_format="%.6g")

# ── (2) Variance decomposition ────────────────────────────────────────────────
# 12 pipelines = 3 measures × 2 decomps × 2 refs; 4 bands
records_var = []

for band in BANDS_ORDER:
    sub = roi[roi["band"] == band].copy()
    assert len(sub) == 12, f"expected 12 pipelines for {band}, got {len(sub)}"

    dz = sub["dz"].values.astype(float)
    var_across = float(np.var(dz, ddof=1))
    sd_across  = float(np.std(dz, ddof=1))
    ratio_sd_to_se = sd_across / SE_DZ

    # Main effects
    for factor, col in [("measure", "measure"),
                        ("decomposition", "decomposition"),
                        ("reference", "reference")]:
        groups = sub.groupby(col)["dz"].mean()
        for level, mean_dz in groups.items():
            records_var.append(dict(
                band=band,
                factor=factor,
                level=level,
                mean_dz=float(mean_dz),
                sd_dz_across_pipelines=round(sd_across, 6),
                var_dz_across_pipelines=round(var_across, 6),
                se_dz_approx=round(SE_DZ, 6),
                sd_over_se=round(ratio_sd_to_se, 4),
            ))

var_df = pd.DataFrame(records_var)
var_df.to_csv(OUT_VARIANCE, index=False, float_format="%.6g")

# ── (3) Participant bootstrap (if per-participant file exists) ─────────────────
bootstrap_done = False
if PER_PART_CSV.exists():
    pp = pd.read_csv(PER_PART_CSV)
    # expected columns: participant, band, measure, decomposition, reference, roi_value
    beta_pp = pp[pp["band"] == "beta"].copy()
    pipeline_cols = ["measure", "decomposition", "reference"]
    pipelines = beta_pp.groupby(pipeline_cols)

    rng = np.random.default_rng(BOOTSTRAP_SEED)
    participants = beta_pp["participant"].unique()
    n_part = len(participants)

    boot_range_dz  = np.empty(N_BOOT)
    boot_range_bf01 = np.empty(N_BOOT)

    pip_keys = [(m, d, r)
                for (m, d, r), _ in beta_pp.groupby(pipeline_cols)]

    for b in range(N_BOOT):
        idx = rng.integers(0, n_part, size=n_part)
        subs = participants[idx]
        dz_vals = []
        bf01_vals = []
        for key in pip_keys:
            m, d, ref = key
            vals = beta_pp[(beta_pp["measure"] == m) &
                           (beta_pp["decomposition"] == d) &
                           (beta_pp["reference"] == ref)]["roi_value"].values
            sample = vals[idx]
            t_b = float(sample.mean() / (sample.std(ddof=1) / math.sqrt(n_part)))
            dz_b = t_b / math.sqrt(n_part)
            dz_vals.append(dz_b)
            bf01_vals.append(1.0 / jzs_bf10(t_b, n_part))
        boot_range_dz[b]  = max(dz_vals) - min(dz_vals)
        boot_range_bf01[b] = max(bf01_vals) - min(bf01_vals)

    boot_records = []
    for metric, arr in [("dz_range", boot_range_dz),
                        ("bf01_range", boot_range_bf01)]:
        lo, hi = np.percentile(arr, [2.5, 97.5])
        boot_records.append(dict(
            band="beta", metric=metric, seed=BOOTSTRAP_SEED, n_boot=N_BOOT,
            mean=float(arr.mean()), median=float(np.median(arr)),
            ci_lo_95=float(lo), ci_hi_95=float(hi),
        ))
    pd.DataFrame(boot_records).to_csv(OUT_BOOTSTRAP, index=False, float_format="%.6g")
    bootstrap_done = True

# ── markdown summary ──────────────────────────────────────────────────────────
lines = []
lines.append("# C. Multiverse summary (group-level Python portion)")
lines.append("")
lines.append(f"N = {N}; 48 cells (12 pipelines × 4 bands); BF₀₁ ceiling = {BF01_CEILING:.4f}.")
lines.append(f"Sampling SE of dz (approx, under H0) = 1/√N = {SE_DZ:.4f}.")
lines.append("")
lines.append("## Holm correction across 48 cells")
lines.append(f"Significant cells (Holm p < 0.05): {int(forest['holm_sig'].sum())}")
lines.append("")
lines.append("## SD of dz across 12 pipelines, by band")
for band in BANDS_ORDER:
    sub = roi[roi["band"] == band]["dz"].values.astype(float)
    sd = float(np.std(sub, ddof=1))
    ratio = sd / SE_DZ
    lines.append(f"- {band:10s}: SD = {sd:.4f}  (SD/SE = {ratio:.2f})")
lines.append("")
lines.append("## Main effects: mean dz by factor level, beta band")
beta_roi = roi[roi["band"] == "beta"]
for factor, col in [("measure", "measure"),
                    ("decomposition", "decomposition"),
                    ("reference", "reference")]:
    groups = beta_roi.groupby(col)["dz"].mean()
    lines.append(f"  {factor}: " + ", ".join(f"{k}={v:.4f}" for k, v in groups.items()))
lines.append("")
if bootstrap_done:
    bt = pd.read_csv(OUT_BOOTSTRAP)
    for _, r in bt.iterrows():
        lines.append(f"Bootstrap {r['metric']} (beta, {N_BOOT} resamples, seed {BOOTSTRAP_SEED}): "
                     f"median = {r['median']:.4f}, 95% CI [{r['ci_lo_95']:.4f}, {r['ci_hi_95']:.4f}]")
else:
    lines.append("## Participant bootstrap")
    lines.append(f"SKIPPED — {PER_PART_CSV.name} not found. Run revision/C_save_roi_per_participant.m first.")
lines.append("")
lines.append(f"Outputs: {OUT_FOREST.name}, {OUT_VARIANCE.name}")

OUT_MD.write_text("\n".join(lines), encoding="utf-8")
print(f"Wrote {OUT_FOREST}")
print(f"Wrote {OUT_VARIANCE}")
if bootstrap_done:
    print(f"Wrote {OUT_BOOTSTRAP}")
print(f"Wrote {OUT_MD}")
