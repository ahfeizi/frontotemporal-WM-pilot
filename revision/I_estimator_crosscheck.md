# I. Estimator cross-check, Dataset 2

Subject: sub-001
Band: theta (4-7 Hz), Hilbert phase, average reference.
Trials: low = 30, high = 30.

Maximum absolute difference (repo vs FieldTrip):
  PLV:   3.838e-02
  dwPLI: 0.000e+00
  ImC:   0.000e+00

## Interpretation
dwPLI and ImC: difference = 0 (machine epsilon). Formulas are identical.

PLV: difference = 3.838e-02. This reflects a legitimate formula difference,
not a bug. util_conn_pairs computes PLV as:
  mean_over_time( abs(mean_over_trials(X/|X|)) )
The FieldTrip arm in this script uses:
  abs( mean_over_all_time_and_trials(X/|X|) )
The repository implementation (util_conn_pairs) is consistent with the
PLV definition used throughout the main battery (s03_connectivity.m).
The cross-check confirms dwPLI and ImC are numerically identical to
the FieldTrip equivalents; PLV uses a time-averaging convention that
differs from a simple pooled-trials formula.

Output: I_estimator_crosscheck.csv (n_pairs x 2 conds x 3 measures = 180 rows).
