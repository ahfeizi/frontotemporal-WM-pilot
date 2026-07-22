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
s00_inventory      % FIRST. confirm events/channels -> fill TODOs in config.m
s01_preprocess     % continuous cleaned data per subject  -> *_clean.set
s02_epoch          % maintenance-window trial list        -> *_trials.mat
s03_connectivity   % dwPLI + PLV over frontal x temporal   -> conn_results.mat/.csv
s04_stats          % sign-flip permutation, high vs low    -> stats_results.mat
s05_figures        % 300-dpi figures                       -> figures/
```

Everything reads parameters from **`config.m`** — the only file you normally edit.

## Dependencies

- MATLAB + **Signal Processing Toolbox** (for `fir1`, `filtfilt`, `hilbert`)
- **EEGLAB** with plugins: `bids-matlab-tools` (BIDS import), `clean_rawdata`
  (ASR), `ICLabel`, `cleanline` (only if `cfg.line_hz` is set)

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
- **Recommended cross-check:** before trusting the numbers, validate the dwPLI/PLV
  implementation against FieldTrip's `ft_connectivityanalysis`
  (`method = 'wpli_debiased'` / `'plv'`) on one subject. Agreement is a cheap,
  convincing robustness line for a methods section.

## Figures

- `fig1_fm_theta.png` — frontal-midline theta power × load. This reproduces the
  established load effect (Onton, Delorme & Makeig, 2005) and is your **sanity
  check**: if this is absent, fix preprocessing/epoching before interpreting
  connectivity.
- `fig2_ft_dwpli.png` — frontotemporal dwPLI in theta and beta × load. This is
  the **headline** preliminary result and the figure most likely to anchor both
  the proposal's feasibility section and the preprint.

## Proposal ↔ preprint mapping

| Output | Proposal use | Preprint use |
|---|---|---|
| fig1 (sanity) | shows method reproduces known physiology | Methods/validation |
| fig2 (dwPLI × load) | feasibility: the connectivity contrast you propose to manipulate with tACS is measurable | Main result |
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

## Status

These scripts were written for your local Windows MATLAB/EEGLAB and have **not**
been executed against the dataset here. Treat the first `s00`/`s01` pass as a
structure-confirmation step and adjust I/O if your EEGLAB version's BIDS import
differs. The connectivity and statistics core is the part designed to stay fixed.
