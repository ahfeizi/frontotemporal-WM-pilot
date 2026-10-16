# Running the convergent-validation battery (ds004117)

## Before you run

1. In `config.m`, check the **BATTERY ADDITIONS** block at the end. Set `cfg.fieldtrip_dir` to your FieldTrip install. It is a relative placeholder, like `cfg.eeglab_dir`.
2. These existing outputs must already be in `cfg.deriv_dir`. Nothing is recomputed.
   - `manifest.mat`, `<key>_clean.set`/`.fdt` and `<key>_trials.mat` (from `s01_preprocess.m`, `s02_epoch.m`)
   - For the check only: `conn_results.csv` (`s03`), `conn_results_betasplit.csv` (`s09`), `conn_results_wavelet.csv` (`s16`)
3. Set MATLAB's current folder to this repository. That way `config.m`, `util_bandpass.m`, `util_chanidx.m` and `util_conn_pairs.m` all come from this repo.

## Toolboxes on the path

| Toolbox | Used for | How it is added |
|---|---|---|
| EEGLAB | `pop_loadset`, `eeglab2fieldtrip` | `addpath(cfg.eeglab_dir); eeglab nogui` (inside the scripts) |
| FieldTrip | `ft_scalpcurrentdensity` (CSD, method `'spline'`) | `addpath(cfg.fieldtrip_dir); ft_defaults` (inside the scripts) |
| Signal Processing Toolbox | `fir1`, `filtfilt`, `hilbert` | MATLAB installation |

## Run order

1. `s18_battery.m` writes these files to `cfg.deriv_dir`:
   - `battery_persubject.csv`
   - `battery_phasediff.csv`
   - `battery_phasediff_bins.csv`
   - `battery_done_subjects.txt`
2. `s19_battery_check.m` prints PASS/FAIL for each check, then an `OVERALL` line.

## Resuming and restarting

`s18_battery.m` appends one subject at a time. After a crash, run it again:
- Subjects listed in `battery_done_subjects.txt` are skipped.
- Rows from a subject that was interrupted partway are removed before that subject is recomputed.

To start from scratch, delete the three `battery_*.csv` files and `battery_done_subjects.txt`.

## Electrode positions (CSD)

CSD uses the 3-D positions that `s01_preprocess.m` loaded from the BIDS `*_electrodes.tsv` files (`pop_importbids`, `'bidschanloc','on'`). If any channel of a run has no X/Y/Z, `s18_battery.m` stops with error `s18_battery:no_positions` and names the channels. It never makes up coordinates.

## Expected runtime (estimate, not measured)

- **Per subject:** about 1–3 min (roughly 4 run files at 250 Hz). The script prints the actual minutes for each subject.
- **Where the time goes:** for each run file, loading the `.set` and computing the CSD once; then, for each reference, 4 Hilbert bands and 13 Morlet bins on the 11 ROI channels.
- **All 22 subjects:** about 0.5–1 h.
- **Check script:** under a minute.
- **Output size:** `battery_persubject.csv` has about 63,000 rows (22 × 2 × 4 × 3 × 2 × 2 × 30).
