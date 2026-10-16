"""s17_v2_checks.py -- additional checks for preprint v2 (ds004117, N = 22).

Usage:  python s17_v2_checks.py <results-dir>

Reads only files in results/ and writes four CSVs there:
  sensitivity_without_sub024.csv  key contrasts without sub-024, who scored
                                  below chance (36% / 44%) in both conditions
  emg_ratio_analysed.csv          temporal/frontal 30-45 Hz power ratio on the
                                  22 analysed participants (s15 also averaged
                                  over sub-005, who was excluded from analysis)
  trial_counts_by_subject.csv     Supplementary Table S2 (analysed correct
                                  trials per participant and load)
  design_limits_n22.csv           largest attainable BF01 and smallest
                                  attainable TOST p (bound d_z = 0.4) at N = 22

p-values here are exact: all 2^n sign-flip permutations are enumerated.
Bayes factors: one-sample JZS, r = 0.707 (Rouder et al., 2009), as in s12; the
design limits use python/util_contrast_stats.py (r = sqrt(2)/2).
Power is log-transformed before differencing, as in s04 (cfg.log_power).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import integrate, stats

RES = Path(sys.argv[1])


def jzs_bf10(t, n, r=0.707):
    nu = n - 1
    f = lambda g: ((1 + n * g * r**2) ** -0.5
                   * (1 + t**2 / ((1 + n * g * r**2) * nu)) ** (-(nu + 1) / 2)
                   * (2 * np.pi) ** -0.5 * g ** -1.5 * np.exp(-1 / (2 * g)))
    num, _ = integrate.quad(f, 0, np.inf, limit=500)
    return num / (1 + t**2 / nu) ** (-(nu + 1) / 2)


def exact_signflip_p(x):
    n = len(x)
    signs = ((np.arange(2**n)[:, None] >> np.arange(n)) & 1) * 2 - 1
    sx = signs * x
    tn = sx.mean(1) / (sx.std(1, ddof=1) / np.sqrt(n))
    to = x.mean() / (x.std(ddof=1) / np.sqrt(n))
    return float((np.abs(tn) >= abs(to) - 1e-12).mean())


def tost_p(dz, n, bound=0.4):
    t = dz * np.sqrt(n); d = bound * np.sqrt(n); df = n - 1
    return max(stats.t.sf(t + d, df), stats.t.cdf(t - d, df))


conn = pd.read_csv(RES / "conn_results.csv")

# 1. Sensitivity without sub-024
rows = []
d = conn[conn.subject != "sub-024"]
for band, metric in [("beta", "fm_power"), ("theta", "dwpli"), ("beta", "dwpli"),
                     ("theta", "plv"), ("beta", "plv")]:
    w = d[d.band == band].pivot(index="subject", columns="load", values=metric)
    if metric == "fm_power":
        w = np.log(w)
    x = (w["high"] - w["low"]).to_numpy()
    n = len(x); t = x.mean() / (x.std(ddof=1) / np.sqrt(n)); bf = jzs_bf10(t, n)
    rows.append(dict(band=band, metric=metric, n=n, dz=x.mean() / x.std(ddof=1), t=t,
                     p_exact=exact_signflip_p(x), BF10=bf, BF01=1 / bf))
pd.DataFrame(rows).to_csv(RES / "sensitivity_without_sub024.csv", index=False)

# 2. EMG ratio on analysed participants
emg = pd.read_csv(RES / "emg_band_power.csv")
analysed = set(conn.subject)
e = emg[emg.subject.isin(analysed)]
pd.DataFrame([dict(n=len(e), mean_temporal_frontal_ratio=e.temporal_frontal_ratio.mean(),
                   n_all_rows=len(emg),
                   mean_ratio_all_rows=emg.temporal_frontal_ratio.mean())]
             ).to_csv(RES / "emg_ratio_analysed.csv", index=False)

# 3. Supplementary Table S2
s2 = (conn[conn.band == "theta"].pivot(index="subject", columns="load", values="ntrials")
      .rename(columns={"low": "ntrials_low", "high": "ntrials_high"})
      [["ntrials_low", "ntrials_high"]])
s2.to_csv(RES / "trial_counts_by_subject.csv")

# 4. Design limits at N = 22, from the shared inference module (r = sqrt(2)/2, as in
#    the battery tables); the local r = 0.707 gave a ceiling of 4.4849 instead of 4.4854
sys.path.insert(0, str(Path(__file__).resolve().parent / "python"))
import util_contrast_stats as ucs  # noqa: E402

n = conn.subject.nunique()
pd.DataFrame([dict(n=n, max_BF01=1 / ucs.jzs_bf10(0.0, n),
                   min_TOST_p_bound_0p4=ucs.tost_p(0.0, n))]
             ).to_csv(RES / "design_limits_n22.csv", index=False)

print(pd.read_csv(RES / "sensitivity_without_sub024.csv").round(3).to_string(index=False))
print(pd.read_csv(RES / "emg_ratio_analysed.csv").round(3).to_string(index=False))
print(s2.describe().round(1))
print(pd.read_csv(RES / "design_limits_n22.csv").round(3).to_string(index=False))
