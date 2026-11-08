"""
D_narrowband_check.py
Narrowband check: fractional bandwidth and |dz(Hilbert) - dz(wavelet)|
per measure × reference × band.  Table only, no statistical test.

Band definitions (from §2.3):
  theta:     4–7 Hz    (wavelet center freqs 4–7 Hz in 1-Hz steps)
  beta:     13–30 Hz   (wavelet center freqs 13–29 Hz in 2-Hz steps)
  beta_low: 13–20 Hz   (wavelet center freqs 13–19 Hz in 2-Hz steps)
  beta_high:21–30 Hz   (wavelet center freqs 21–29 Hz in 2-Hz steps)

Fractional bandwidth = (f_high - f_low) / f_center
  where f_center = (f_high + f_low) / 2.

Reads:  results/unified_contrasts_ds004117.csv (ROI rows only)
Writes: revision/D_narrowband_check.csv
        revision/D_narrowband_check.md
"""

from __future__ import annotations
import pathlib

import numpy as np
import pandas as pd

REPO     = pathlib.Path(__file__).resolve().parent.parent
DATA_CSV = REPO / "results" / "unified_contrasts_ds004117.csv"
OUT_CSV  = REPO / "revision" / "D_narrowband_check.csv"
OUT_MD   = REPO / "revision" / "D_narrowband_check.md"

# Band edge definitions (Hz) — from manuscript §2.3
BAND_EDGES = {
    "theta":     (4, 7),
    "beta":      (13, 30),
    "beta_low":  (13, 20),
    "beta_high": (21, 30),
}

# ── load ROI rows ─────────────────────────────────────────────────────────────
df = pd.read_csv(DATA_CSV)
roi = df[df["level"] == "roi"].copy()
assert len(roi) == 48

# ── compute fractional bandwidth ─────────────────────────────────────────────
def frac_bw(f_lo, f_hi):
    fc = (f_lo + f_hi) / 2.0
    return (f_hi - f_lo) / fc

# ── compute |dz(Hilbert) - dz(wavelet)| per measure × reference × band ────────
records = []
for band, (f_lo, f_hi) in BAND_EDGES.items():
    fb = frac_bw(f_lo, f_hi)
    band_df = roi[roi["band"] == band]
    for measure in ("PLV", "dwPLI", "ImC"):
        for reference in ("average", "csd"):
            row_h = band_df[(band_df["measure"] == measure) &
                            (band_df["decomposition"] == "hilbert") &
                            (band_df["reference"] == reference)]
            row_w = band_df[(band_df["measure"] == measure) &
                            (band_df["decomposition"] == "wavelet") &
                            (band_df["reference"] == reference)]
            assert len(row_h) == 1, f"no unique hilbert row: {band} {measure} {reference}"
            assert len(row_w) == 1, f"no unique wavelet row: {band} {measure} {reference}"
            dz_h = float(row_h["dz"].iloc[0])
            dz_w = float(row_w["dz"].iloc[0])
            records.append(dict(
                band=band,
                f_lo=f_lo,
                f_hi=f_hi,
                f_center=(f_lo + f_hi) / 2.0,
                fractional_bandwidth=round(fb, 4),
                measure=measure,
                reference=reference,
                dz_hilbert=round(dz_h, 6),
                dz_wavelet=round(dz_w, 6),
                abs_diff_dz=round(abs(dz_h - dz_w), 6),
            ))

out = pd.DataFrame(records)
out.to_csv(OUT_CSV, index=False, float_format="%.6g")

# ── markdown summary ──────────────────────────────────────────────────────────
lines = []
lines.append("# D. Narrowband check")
lines.append("")
lines.append("Fractional bandwidth = (f_high − f_low) / f_center.")
lines.append("")

lines.append("## Fractional bandwidth by band")
for band, (f_lo, f_hi) in BAND_EDGES.items():
    fb = frac_bw(f_lo, f_hi)
    lines.append(f"- {band:12s}: {f_lo}–{f_hi} Hz, f_c = {(f_lo+f_hi)/2:.1f} Hz, "
                 f"fractional BW = {fb:.3f}")
lines.append("")

lines.append("## |dz(Hilbert) − dz(wavelet)| per band (mean ± max across 6 measure×ref cells)")
for band in BAND_EDGES:
    sub = out[out["band"] == band]["abs_diff_dz"]
    lines.append(f"- {band:12s}: mean = {sub.mean():.4f},  max = {sub.max():.4f}")
lines.append("")

lines.append("## Full table (n=4 bands; see D_narrowband_check.csv)")
header = f"{'band':12s} {'frac_BW':8s} {'measure':7s} {'ref':8s} "
header += f"{'dz_H':7s}  {'dz_W':7s}  {'|diff|':7s}"
lines.append(header)
for _, r in out.iterrows():
    lines.append(f"{r['band']:12s} {r['fractional_bandwidth']:8.3f} "
                 f"{r['measure']:7s} {r['reference']:8s} "
                 f"{r['dz_hilbert']:+7.4f}  {r['dz_wavelet']:+7.4f}  "
                 f"{r['abs_diff_dz']:7.4f}")
lines.append("")
lines.append(f"Output: {OUT_CSV.name}")

OUT_MD.write_text("\n".join(lines), encoding="utf-8")
print(f"Wrote {OUT_CSV}")
print(f"Wrote {OUT_MD}")
