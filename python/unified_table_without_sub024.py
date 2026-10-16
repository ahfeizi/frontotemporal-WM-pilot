"""
unified_table_without_sub024.py -- ROI-level battery for ds004117 without sub-024 (N = 21).

Applies the computation of unified_table.py to battery_persubject.csv with sub-024
removed, ROI level only (no pair rows, so no max-|t| FWER), for Supplementary Table S2.
Loading, the load contrast (high - low), the ROI average over pairs and the inference
(util_contrast_stats.contrast_stats, exact sign-flip since n <= 22) are the functions of
unified_table.py itself, imported rather than copied.

Writes <results>/unified_contrasts_ds004117_without_sub024.csv with the columns of
unified_contrasts_ds004117.csv (p_fwer is empty, as it is on every ROI row there).

Usage
  python unified_table_without_sub024.py [--results <dir>]
Default --results is ../results relative to this file.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import unified_table as ut
import util_contrast_stats as ucs

EXCLUDE = ["sub-024"]
COLUMNS = ["dataset", "contrast", "level", "band", "measure", "decomposition", "reference",
           "pair", "n", "dz", "dz_ci_lo", "dz_ci_hi", "t", "df", "p", "p_perm", "p_fwer",
           "perm_method", "bf10", "bf01", "p_tost"]
DEFAULT_RESULTS = Path(__file__).resolve().parent.parent / "results"


def roi_contrasts(d: pd.DataFrame) -> pd.DataFrame:
    """ROI branch of unified_table.contrasts, unchanged."""
    rows = []
    for cell, g in d.groupby(ut.CELL, sort=True):
        base = dict(zip(ut.CELL, cell), dataset=ut.SPEC["P1"]["dataset"],
                    contrast=ut.SPEC["P1"]["contrast"])
        roi = g.groupby("subject")["diff"].mean()
        npairs = g.groupby("subject").pair.nunique()
        if npairs.nunique() != 1:
            raise ValueError(f"P1 {cell}: unequal pair counts across subjects")
        rows.append({**base, "level": "roi", "pair": f"ROI({npairs.iloc[0]})",
                     **ucs.contrast_stats(roi.values)})
    return pd.DataFrame(rows).reindex(columns=COLUMNS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    a = ap.parse_args()

    v, _ = ut.load(a.results, EXCLUDE)
    if set(EXCLUDE) & set(v.subject):
        raise ValueError("exclusion failed")
    tab = roi_contrasts(ut.diff_wide(v, "P1", flip=False))
    out = a.results / "unified_contrasts_ds004117_without_sub024.csv"
    tab.to_csv(out, index=False, float_format="%.6g")

    def cell(band, measure, dec, ref):
        return tab[(tab.band == band) & (tab.measure == measure)
                   & (tab.decomposition == dec) & (tab.reference == ref)].iloc[0]

    b_h = cell("beta", "dwPLI", "hilbert", "average")
    b_w = cell("beta", "dwPLI", "wavelet", "average")
    th = cell("theta", "PLV", "wavelet", "csd")
    print(f"wrote {out} ({len(tab)} rows)")
    print(f"n = 21 in every row: {bool((tab.n == 21).all())} (n values: {sorted(tab.n.unique())})")
    print(f"perm_method: {sorted(tab.perm_method.unique())}")
    print(f"beta dwPLI Hilbert average BF01 = {b_h.bf01:.5f} (expected 3.67)")
    print(f"beta dwPLI wavelet average BF01 = {b_w.bf01:.5f} (expected 2.03)")
    print(f"theta PLV wavelet CSD BF10 = {th.bf10:.5f} (manuscript 3.49), "
          f"parametric p = {th.p:.5f}, exact p = {th.p_perm:.5f} (manuscript 0.016)")
    print(f"ceiling BF01 at N = 21 = {1.0 / ucs.jzs_bf10(0.0, 21):.5f} (expected 4.40)")


if __name__ == "__main__":
    main()
