# manuscript/

Figure scripts for the merged manuscript (Dataset 1 = ds003702, Dataset 2 = ds004117).
They read released results files only and run no new inference beyond
`python/util_contrast_stats.py`.

| Script | Makes |
|---|---|
| `make_main_figures.py` | Figures 1–4 (PDF and PNG) and `Figure3_source_values.csv` |
| `make_figS1_prisma.py` | Supplementary Figure S1, the PRISMA-ScR flow diagram of the literature audit |

## Running

```
python manuscript/make_main_figures.py --p2 <frontoparietal-gamma-PLV-dwPLI-pilot at v1.2.1> --p1 <this repo> --out <dir>
python manuscript/make_figS1_prisma.py
```

`--p2` and `--p1` are repository roots: Figures 1–3 read
`results/unified_contrasts_ds003702.csv`, `results/battery_phasediff.csv` and
`results/imc_diag_persubject.csv` from
[ahfeizi/frontoparietal-gamma-PLV-dwPLI-pilot](https://github.com/ahfeizi/frontoparietal-gamma-PLV-dwPLI-pilot)
v1.2.1, and Figure 4 reads `results/unified_contrasts_ds004117.csv` from this repository.
`make_figS1_prisma.py` writes `FigureS1_PRISMA.png/.pdf` to the current directory.

Requirements: numpy, pandas, scipy, matplotlib (Python 3.12; see `python/requirements.txt`
for the pinned numpy, scipy and pandas versions).

The PRISMA counts in `make_figS1_prisma.py` are hard-coded from
wm-phase-connectivity-audit v1.0.0 (doi:10.5281/zenodo.23426096): `stage0_output/queries.json`,
`screening_sheet_with_D12.csv` and deviations D9 and D12.
