"""
power_exact.py -- frontal-midline power load effect (high - low) for one band.

Reads conn_results.csv (s03_connectivity.m), takes fm_power per subject x load for
--band, and runs util_contrast_stats.contrast_stats on log(high) - log(low), for all
subjects and again without sub-024. With N <= 22 the sign-flip p is exact (all 2^N
flips).

Usage
  python power_exact.py --band beta  [--csv <path to conn_results.csv>]
  python power_exact.py --band theta
Default --csv is ../results/conn_results.csv relative to this file.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from util_contrast_stats import contrast_stats

DEFAULT_CSV = Path(__file__).resolve().parent.parent / "results" / "conn_results.csv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", required=True, choices=["theta", "beta"])
    ap.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    a = ap.parse_args()

    c = pd.read_csv(a.csv)
    b = c[c["band"] == a.band].pivot(index="subject", columns="load", values="fm_power")
    if b.empty:
        raise ValueError(f"no rows for band {a.band!r} in {a.csv}")
    for excl in ([], ["sub-024"]):
        x = b.drop(index=excl)
        r = contrast_stats((np.log(x["high"]) - np.log(x["low"])).values)
        print(
            a.band,
            excl or "all",
            {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()},
        )


if __name__ == "__main__":
    main()
