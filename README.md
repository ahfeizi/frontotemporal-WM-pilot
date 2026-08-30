# Frontotemporal Connectivity Pilot — ds004117 (Sternberg verbal WM)

A modular MATLAB/EEGLAB pipeline that estimates **frontotemporal phase
connectivity** (debiased wPLI, with PLV as a secondary metric) in the **theta**
and **beta** bands during the **maintenance** period of a verbal working-memory
task, contrasted across **memory load**.

It is built to serve two purposes from one run:
1. a clean **preliminary-data figure** for the MSCA proposal (feasibility), and
2. a **preprint-ready** analysis skeleton (publication-grade stats from the start).

---

## Run order

```
s00_inventory       % FIRST. confirm events/channels -> fill TODOs in config.m
s01_preprocess      % continuous cleaned data per subject  -> *_clean.set
s02_epoch           % maintenance-window trial list        -> *_trials.mat
s03_connectivity    % dwPLI + PLV over frontal x temporal   -> conn_results.mat/.csv
s04_stats           % sign-flip permutation, high vs low    -> stats_results.mat
s05_figures         % 300-dpi headline figures               -> figures/
s06_cluster_conn    % FWER max-stat test over full 60-pair x band space
                     %   -> conn_cluster_results.csv, fig3_conn_tmap.png
s07_grandavg_ersp   % grand-average ERSP at FCz (pipeline validation figure)
                     %   -> grandavg_ersp.mat, fig_grandavg_ersp.png
```

Everything reads parameters from **`config.m`** — the only file you normally edit.

## Dependencies

- MATLAB + **Signal Processing Toolbox** (for `fir1`, `filtfilt`, `hilbert`)
- **EEGLAB** with plugins: `bids-matlab-tools` (BIDS import), `clean_rawdata`
  (ASR), `ICLabel`, `cleanline` (only if `cfg.line_hz` is set)
- `s07_grandavg_ersp.m` additionally requires EEGLAB's `newtimef` (bundled with
  EEGLAB; no extra plugin needed)

## Resolved from the s00 inventory (already set in config.m)

- `cfg.ev_maint_lock = 'WM'` — the `show_dash` / work_memory event marks the
  start of the retention/maintenance delay.
- `cfg.load_field = 'memory_cond'`, `load_low = {'3'}`, `load_high = {'7'}` —
  memory_cond is the set size (number of to-remember letters). Set size 5, if
  present, is treated as medium and excluded from the low-vs-high contrast.
- `cfg.maint_window = [0.5 2.5]` s — the delay is ~3.58 s (dash -> probe), so
  this window ends ~1 s before the probe and avoids probe-evoked contamination.
- `cfg.correct_only = true` — only correct-response trials enter the analysis.
- ROI labels are matched case-insensitively (the montage uses UPPERCASE).

## Subjects vs datasets (important)

ds004117 imports as ~85 **datasets** (~24 subjects x several runs each), not 85
people. The pipeline therefore writes a `manifest` mapping each dataset to its
subject; `s03` pools trials **within subject across runs** and emits one row per
subject. **N for statistics = number of subjects.** Use `cfg.subset_n` to process
a handful of datasets first (fast preliminary figure), then set it to `[]` for the
full run (ICA over ~85 datasets can take hours).

Final analysis N = **22 subjects** (one subject dropped for insufficient
good/pooled trials; see `s03_connectivity.m`'s `ntr < 5` skip condition and
`conn_results.csv` for the included subject list).

## Method notes (the parts that matter for review)

- **dwPLI** (Vinck et al., 2011) is primary: robust to volume conduction and to
  sample-size bias — essential for scalp frontotemporal pairs. **PLV** is reported
  alongside as a sensitivity check, not as the main inference.
- Connectivity is computed **across trials at each time sample**, then averaged
  over the maintenance window — the standard estimator, not raw pooling.
- Band-pass + Hilbert run on **continuous** data before windowing, so there are
  no within-epoch edge artifacts. ASR is used in **correction** mode so the data
  stays continuous.
- Inference is **non-parametric** (5000 sign-flip permutations), within-subject
  high-vs-low load, with Cohen's dz reported.
- `s06_cluster_conn.m` additionally runs a **max-statistic FWER-corrected**
  test across the full 60-cell frontal×temporal×band space, so a pair-averaged
  null in `s04` cannot hide a genuinely localized effect.

## Figures

- `fig1_fm_theta.png` / `fig1_fm_power.png` — frontal-midline theta and beta
  power × load. Beta reproduces the headline local load effect; theta serves
  as a non-significant control. This is your **sanity check**: if the beta
  effect is absent, fix preprocessing/epoching before interpreting connectivity.
- `fig2_ft_dwpli.png` — frontotemporal dwPLI in theta and beta × load.
- `fig3_conn_tmap.png` (from `s06`) — full frontal×temporal t-map per band,
  with FWER-significant cells marked. This is the figure that answers "is
  there any localized effect hiding in the pair-average."
- `fig_grandavg_ersp.png` (from `s07`) — grand-average ERSP at FCz locked to
  maintenance onset, collapsed across load. This is the **pipeline validation**
  figure: it should reproduce canonical WM-related spectral dynamics
  independent of the load contrast.
- `fig4_topography.png` (from `s10`) — per-channel t-value topography
  (High−Low) for theta and beta, restricted to the analyzed ROI channels.
  Exploratory/visual only — the inferential claims come from `s04`/`s06`.

## Proposal ↔ preprint mapping

| Output | Proposal use | Preprint use |
|---|---|---|
| fig_grandavg_ersp (s07) | shows method reproduces known physiology | §3.2 Pipeline validation |
| fig1 (power, s05) | local load effect (beta) | §3.3 Results |
| fig2 (dwPLI × load, s05) | feasibility: the connectivity contrast you propose to manipulate with tACS is measurable | §3.4 Main connectivity result |
| fig3_conn_tmap + conn_cluster_results.csv (s06) | rules out a diluted localized effect | §3.5 / Supplementary Table S1 |
| stats_results | quantitative feasibility (effect size powers the proposed design) | Results |
| this pipeline | "analysis tools already in hand" → Implementation | Methods + shared code |

## Honest limitations to state plainly

- ds004117 is **verbal but visually presented** (letters on screen), not auditory.
  It supports the verbal arm of "verbal/auditory"; do not overclaim the auditory
  modality from it.
- Scalp temporal electrodes are an imperfect proxy for temporal-lobe sources;
  dwPLI mitigates but does not eliminate this. Frame the pilot as **scalp-level
  feasibility**, with source-level and intracranial corroboration (e.g. ds004752)
  positioned as future/clinical work.
- The pooled-beta (13–30 Hz) frontotemporal dwPLI null is **not treated as a
  securely confirmed null result**: a wavelet-based robustness re-analysis
  changed its Bayesian support (BF01 dropped from ~4.5 to ~2.2). The claim is
  formally withdrawn in the associated preprint in favor of the sub-band-level
  result. *(Note: the wavelet re-analysis and sub-band split are not yet part of
  this pipeline — see the open items below.)*

## Status

This pipeline was executed against the full ds004117 dataset; final statistics
reflect **N = 22 subjects** after exclusion of subjects with too few pooled
good trials (see `results/` for the exact numerical outputs reported in the
associated preprint). Treat the first `s00`/`s01` pass as a structure-confirmation
step and adjust I/O if your EEGLAB version's BIDS import differs. The
connectivity and statistics core (`s01`–`s07`) is the part designed to stay fixed.

## Sub-band, sensitivity, and robustness checks (s08–s16)

These extend `s01`–`s07` with the additional analyses reported in the revised
preprint. Independent of one another **except** `s12`, which requires `s09`'s
output. Run after the core pipeline (`s00`–`s07`) has completed.

```
s08_site_robustness          % beta power at FCz vs Fz vs Cz              -> fcz_robustness_check.csv
s09_beta_split_robustness    % full pipeline re-run for low-/high-beta    -> conn_results_betasplit.*, stats_results_betasplit.mat, conn_cluster_results_betasplit.csv
s10_topography                % per-channel t-value topography (theta/beta) -> fig4_topography.png
s11_accuracy                  % behavioral accuracy, low vs high load       -> accuracy_by_subject_load.csv, accuracy_stats.mat
s12_bayes_tost                % BF10/BF01 (JZS) + TOST, main bands + sub-bands (needs s09) -> bayes_tost_expanded.csv
s13_power_conn_dissociation   % formal test: dz(power) vs dz(connectivity)  -> power_conn_dissociation_test.mat
s14_amplitude_phase_confound  % power/connectivity-change correlation + amplitude-normalized dwPLI sensitivity check -> amplitude_phase_correlation.mat, conn_results_ampnorm.csv, ampnorm_stats.mat
s15_emg_sensitivity           % 30-45 Hz (truncated) EMG-band proxy vs beta dwPLI -> emg_band_power.csv, emg_sensitivity_results.mat
s16_wavelet_beta              % Morlet wavelet re-analysis (13-29 Hz, 9 freqs), replaces Hilbert beta estimate; own BF/TOST -> conn_results_wavelet.*, stats_results_wavelet.mat, bayes_tost_wavelet.csv
```

**Cross-check before trusting these numbers:** the pooled-beta dwPLI Bayes
factor should read BF01 ≈ 4.5 from the primary Hilbert-based pipeline
(`s04`/`s12`) and BF01 ≈ 2.2 from the wavelet re-analysis (`s16`) — matching
the values already reported in the preprint. If these don't match, something
in the re-implementation has drifted from what the preprint actually reports.

**Known caveats already documented inline in these scripts (kept here for
visibility):**
- `s15`: the EMG-band proxy is truncated to 30–45 Hz because `cfg.lp_hz = 45`
  removes everything above that; true myogenic power mostly lives higher.
  Muscle-IC rejection counts are not recoverable without a full `s01` rerun.
- `s14` Part 2: amplitude normalization is applied per trial window (not to
  the continuous signal), so this sensitivity check is more exposed to edge
  effects than the primary pipeline.
