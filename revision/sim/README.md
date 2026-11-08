# R3 — Dataset 2 injection power curve

The second half of the known-source simulation. A forward-projected **lagged
F3–T7 source pair** is injected into the **real preprocessed high-load trials**
at five coupling strengths, and the battery is asked how often it notices.

N = 22, each subject's real trial counts, each subject's runs pooled exactly as
`s03_connectivity.m` pools them.

This is a **power analysis within high load**, not a re-run of the high–low
contrast. The scenario simulation (S0–S4, both datasets) lives in P2's
`revision/sim/`.

## How the two arms differ, and only how

Both arms receive the **same real trials**, the **same two sources** at the
**same amplitude**. Only the coupling differs, and source 2's mix is variance
preserving:

```
arm A :  s2 = n2A                                   (coupling 0)
arm B :  s2 = kappa * delayed(s1) + sqrt(1-kappa^2) * n2B
```

so the injected power is identical in both arms at every coupling level and
`kappa = 0` is a genuine **false-positive row** — the arms then differ only by
two independent noise realisations of source 2. Because both arms share one
real recording, the paired comparison is as sensitive as the design allows,
which is what a power analysis wants.

Injection amplitude is `cfg.inject_frac` of **each subject's own** observed ROI
amplitude in the target band, measured from that subject's real data.

The **beta band carries no injection** and is therefore the specificity check.

## Nothing here re-implements the battery

| What | Comes from |
|---|---|
| PLV / dwPLI / ImC | `util_conn_pairs.m` |
| band-pass + Hilbert | `util_bandpass.m` + `hilbert`, as `s03`/`s18` |
| Morlet, 5 cycles | the `s18_battery.m` lines, from `s16_wavelet_beta.m` |
| CSD | the operator cached by P2's `sim_util_csd.m`, verified by its `sim_csd_check.m` |
| dz, 95 % CI, sign-flip p, JZS BF, TOST | `util_contrast_stats.contrast_stats` |

`s18`'s local functions are **copied verbatim** into `sim_inject_battery.m`
because MATLAB local functions cannot be shared between files.

The **forward model is built from each recording's own `chanlocs`**, converted
by FieldTrip's `eeglab2fieldtrip`. No coordinate convention is assumed, and the
68-versus-69-channel question raised by the supplied BIDS table does not arise
here at all.

## Run order

```matlab
cd <P1>/revision/sim

sim_inject_probe     % 1. Measures the cost and projects the whole stage.
sim_inject_power     % 2. The injection. Resumable, one subject at a time.
sim_inject_curve     % 3. Inference + the power curve.
```

Edit the `% EDIT` lines of `sim_inject_config.m` first — in particular
`cfg.shared_sim_dir`, which must point at P2's `revision/sim/`, where the shared
utilities live (`sim_util_path.m`, `sim_util_csd.m`, `sim_util_stats.m`,
`sim_stats_runner.py`).

P2's `sim_csd_check.m` should have passed before you trust any CSD row here:
the CSD operator caching it validates is the same mechanism this stage uses.

`sim_inject_power` is resumable: each subject is appended to
`sim_inject_contrasts.csv` and only then listed in `sim_inject_done.txt`. To
start over, delete both.

## Outputs (`results/`)

| File | Contents |
|---|---|
| `sim_inject_timing.md` | measured per-stage cost and the projected total |
| `sim_inject_contrasts.csv` | one row per subject × kappa × rep × cell, the raw contrasts |
| `sim_inject_percell.csv` | the same cells with full inference: dz, 95 % CI, t, p, sign-flip p, BF₁₀, BF₀₁, N, seed |
| `sim_inject_curve.csv` | detection probability, P(BF₁₀>3), P(TOST passes), mean dz per kappa × measure × decomposition × reference |
| `sim_inject_power.md` | the summary |

## Cost, and the shortcut not taken

Every arm at every coupling level is filtered and Hilberted from scratch. Since
the average reference, the CSD operator, `filtfilt` and `hilbert` are all
linear, the analytic signal of (real + injection) equals that of the real data
plus that of the injection, so the real data's transform could be computed once
per subject and reused across all six arms of every repetition — an exact
algebraic identity worth roughly an order of magnitude.

That is a **second** caching shortcut beyond the CSD operator that was
approved, so it is deliberately **not** taken. Ask for it and it goes in, with a
check script of the same shape as `sim_csd_check.m`.

`sim_inject_probe.m` gives the real number before you commit to the full run.
The knobs that scale it linearly are `cfg.n_rep_hilbert` and
`cfg.n_rep_wavelet`.
