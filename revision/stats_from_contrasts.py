"""
stats_from_contrasts.py  (parallel version)
Helper called by MATLAB revision scripts via system().

Usage:
  python stats_from_contrasts.py <input_csv> <output_csv>

Input CSV: one row per cell, columns:
  cell_key, n, x1, x2, ..., x_N

Output CSV: one row per cell with full inference.
Identical output to the original, but uses all CPU cores.
"""
from __future__ import annotations
import pathlib
import sys
import multiprocessing as mp

import pandas as pd
import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))
from util_contrast_stats import contrast_stats  # noqa: E402

KEY_COL = "cell_key"
N_COL   = "n"


def _process_row(args):
    key, x_vals = args
    x = np.asarray(x_vals, dtype=float)
    x = x[np.isfinite(x)]
    try:
        s = contrast_stats(x)
    except Exception as e:
        s = dict(n=len(x), mean=float("nan"), sd=float("nan"), dz=float("nan"),
                 t=float("nan"), df=len(x)-1, p=float("nan"), p_perm=float("nan"),
                 perm_method="error", bf10=float("nan"), bf01=float("nan"),
                 p_tost=float("nan"), dz_ci_lo=float("nan"), dz_ci_hi=float("nan"))
        print(f"  WARNING: {key}: {e}", file=sys.stderr)
    rec = {KEY_COL: key}
    rec.update(s)
    return rec


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: python {sys.argv[0]} <input_csv> <output_csv>")

    in_path  = pathlib.Path(sys.argv[1])
    out_path = pathlib.Path(sys.argv[2])

    df = pd.read_csv(in_path)
    x_cols = [c for c in df.columns if c not in (KEY_COL, N_COL)]

    tasks = [(row[KEY_COL], row[x_cols].values.tolist()) for _, row in df.iterrows()]

    n_workers = max(1, mp.cpu_count() - 1)
    print(f"Processing {len(tasks)} rows with {n_workers} workers ...", flush=True)

    with mp.Pool(processes=n_workers) as pool:
        records = pool.map(_process_row, tasks, chunksize=50)

    out = pd.DataFrame(records)
    out.to_csv(out_path, index=False, float_format="%.8g")
    print(f"Wrote {out_path} ({len(out)} rows)")
