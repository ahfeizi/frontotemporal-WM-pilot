"""Main-text Figures 1-4 for the merged manuscript.

Inputs (released results files; pass the two repository roots):
  P2 = ahfeizi/frontoparietal-gamma-PLV-dwPLI-pilot  (Dataset 1, ds003702)
  P1 = ahfeizi/frontotemporal-WM-pilot              (Dataset 2, ds004117)
Usage: python make_main_figures.py --p2 <P2 root> --p1 <P1 root> --out <dir>
Statistics come from each repo's python/util_contrast_stats.py (no new inference).
"""
import argparse, importlib.util
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ap = argparse.ArgumentParser()
ap.add_argument("--p2", type=Path, required=True)
ap.add_argument("--p1", type=Path, required=True)
ap.add_argument("--out", type=Path, default=Path("."))
a = ap.parse_args()
a.out.mkdir(parents=True, exist_ok=True)

spec = importlib.util.spec_from_file_location("ucs", a.p1 / "python" / "util_contrast_stats.py")
ucs = importlib.util.module_from_spec(spec); spec.loader.exec_module(ucs)

# validated categorical palette (dataviz reference, light mode); marker shape is the
# secondary encoding for every color, so identity never rests on color alone
COL = {"PLV": "#2a78d6", "dwPLI": "#eb6834", "ImC": "#1baf7a"}
MRK = {"PLV": "o", "dwPLI": "s", "ImC": "^"}
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3df"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
DEC = {"hilbert": "Hilbert", "wavelet": "Wavelet"}
REF = {"average": "Average", "csd": "CSD"}

def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(a.out / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

def fmt_bf(r):
    return f"BF$_{{10}}$ = {r.bf10:.1f}" if r.bf10 >= 1 else f"BF$_{{01}}$ = {r.bf01:.2f}"

# ---------------------------------------------------------------- Figure 1
d1 = pd.read_csv(a.p2 / "results" / "unified_contrasts_ds003702.csv")
roi1 = d1[d1.level == "roi"]
order = [(m, dc, rf) for rf in ("average", "csd") for m in ("PLV", "dwPLI", "ImC")
         for dc in ("hilbert", "wavelet")]
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.9), sharey=True)
for ax, band, title in zip(axes, ("gamma_high", "gamma_low"),
                           ("A  High gamma (55–80 Hz)", "B  Low gamma (30–45 Hz)")):
    ax.axvspan(-0.4, 0.4, color=GRID, zorder=0, lw=0)
    ax.axvline(0, color=MUTED, lw=0.8, zorder=1)
    for i, (m, dc, rf) in enumerate(order):
        r = roi1[(roi1.band == band) & (roi1.measure == m) & (roi1.decomposition == dc)
                 & (roi1.reference == rf)].iloc[0]
        y = len(order) - 1 - i - (0.6 if rf == "csd" else 0)
        ax.plot([r.dz_ci_lo, r.dz_ci_hi], [y, y], color=COL[m], lw=2, solid_capstyle="round", zorder=2)
        ax.plot(r.dz, y, MRK[m], ms=6, color=COL[m], mec="white", mew=1, zorder=3)
        ax.text(1.02, y, fmt_bf(r), transform=ax.get_yaxis_transform(), va="center",
                fontsize=6.5, color=INK)
    ax.set_xlim(-1.0, 0.75)
    ax.set_xlabel("Validity effect, $d_z$ (valid − invalid) with 95% CI")
    ax.set_title(title, loc="left", fontsize=8.5, color=INK, fontweight="bold")
    ax.axhline(5.7, color=MUTED, lw=0.5, ls=(0, (2, 2)))
labels = [f"{m}, {DEC[dc]}" for (m, dc, rf) in order]
ys = [len(order) - 1 - i - (0.6 if rf == "csd" else 0) for i, (_, _, rf) in enumerate(order)]
axes[0].set_yticks(ys); axes[0].set_yticklabels(labels)
axes[0].text(-0.98, 11.45, "Average reference", fontsize=7, color=MUTED, va="bottom")
axes[0].text(-0.98, 4.75, "Current source density", fontsize=7, color=MUTED, va="bottom")
fig.subplots_adjust(wspace=0.55)
save(fig, "Figure1_dataset1_battery")

# ---------------------------------------------------------------- Figure 2
ph = pd.read_csv(a.p2 / "results" / "battery_phasediff.csv")
f2 = ph[(ph.band == "gamma_high") & (ph.pair == "F2-P4")]
fig = plt.figure(figsize=(6.0, 6.2))
bins = np.deg2rad(np.arange(-180, 181, 10))
k = 0
for rf in ("average", "csd"):
    for dc in ("hilbert", "wavelet"):
        k += 1
        ax = fig.add_subplot(2, 2, k, projection="polar")
        x = f2[(f2.reference == rf) & (f2.decomposition == dc)]
        ang = np.deg2rad(((x.mean_angle_deg + 180) % 360) - 180)
        cnt, _ = np.histogram(ang, bins=bins)
        for w0 in (0, np.pi):  # zero-lag windows +-20 deg
            ax.bar(w0, cnt.max() * 1.08, width=np.deg2rad(40), color=GRID, lw=0, zorder=0)
        ax.bar(bins[:-1] + np.deg2rad(5), cnt, width=np.deg2rad(9), color=COL["PLV"], lw=0, zorder=2)
        zl = int(((np.abs(np.rad2deg(ang)) <= 20) | (np.abs(np.rad2deg(ang)) >= 160)).sum())
        ax.set_title(f"{'ABCD'[k-1]}  {REF[rf]} reference, {DEC[dc]}", fontsize=8, color=INK, fontweight="bold", pad=14)
        ax.text(0.5, -0.16, f"median R = {x.resultant_len.median():.2f}; {zl}/{len(x)} cells within ±20° of 0° or 180°",
                transform=ax.transAxes, ha="center", fontsize=6.5, color=MUTED)
        ax.set_yticklabels([]); ax.set_xticks(np.deg2rad([0, 90, 180, 270]))
        ax.set_xticklabels(["0°", "90°", "180°", "−90°"], fontsize=6.5)
        ax.grid(color=GRID, lw=0.5)
fig.subplots_adjust(hspace=0.55, wspace=0.35)
save(fig, "Figure2_F2P4_phase_difference")

# ---------------------------------------------------------------- Figure 3
dg = pd.read_csv(a.p2 / "results" / "imc_diag_persubject.csv")
dg["valid"] = dg.condition.str.endswith("_valid")
ab = dg[dg.measure.str.startswith("imc")].copy()
ab["measure"] = ab.measure.str.replace("imc", "absimc"); ab["value"] = ab.value.abs()
dg = pd.concat([dg, ab])
rows = []
for (m, dc, rf), x in dg.groupby(["measure", "decomposition", "reference"]):
    roi = x.groupby(["subject", "condition", "valid"]).value.mean().reset_index()
    c = roi.groupby(["subject", "valid"]).value.mean().unstack()
    s = ucs.contrast_stats((c[True] - c[False]).values)
    rows.append(dict(measure=m, dec=dc, ref=rf, **{k_: s[k_] for k_ in ("dz", "dz_ci_lo", "dz_ci_hi", "bf10", "bf01")}))
st = pd.DataFrame(rows)
st.to_csv(a.out / "Figure3_source_values.csv", index=False)
mlab = [("imc_joint", "Signed ImC, joint", "ImC", True), ("imc_trial", "Signed ImC, per sample", "ImC", True),
        ("absimc_joint", "|ImC|, joint", "ImC", False), ("absimc_trial", "|ImC|, per sample", "ImC", False),
        ("dwpli_trial", "dwPLI, per sample", "dwPLI", True), ("dwpli_joint", "dwPLI, joint", "dwPLI", True)]
fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True)
for ax, rf in zip(axes, ("average", "csd")):
    ax.axvspan(-0.4, 0.4, color=GRID, zorder=0, lw=0); ax.axvline(0, color=MUTED, lw=0.8)
    for i, (key, lab, fam, filled) in enumerate(mlab):
        for j, dc in enumerate(("hilbert", "wavelet")):
            r = st[(st.measure == key) & (st.dec == dc) & (st.ref == rf)].iloc[0]
            y = len(mlab) - 1 - i + (0.15 if dc == "hilbert" else -0.15)
            ax.plot([r.dz_ci_lo, r.dz_ci_hi], [y, y], color=COL[fam], lw=1.6, solid_capstyle="round")
            ax.plot(r.dz, y, "o" if dc == "hilbert" else "D", ms=5, color=COL[fam] if filled else "white",
                    mec=COL[fam], mew=1.2)
    ax.set_title(("A  Average reference" if rf == "average" else "B  Current source density"),
                 loc="left", fontsize=8.5, fontweight="bold", color=INK)
    ax.set_xlim(-0.85, 0.6); ax.set_xlabel("Validity effect, $d_z$ with 95% CI")
axes[0].set_yticks(range(len(mlab))); axes[0].set_yticklabels([m[1] for m in mlab][::-1])
axes[1].legend(handles=[Line2D([], [], marker="o", ls="", color=MUTED, label="Hilbert"),
                        Line2D([], [], marker="D", ls="", color=MUTED, label="Wavelet")],
               frameon=False, loc="lower right", fontsize=7)
save(fig, "Figure3_imc_diagnostic")

# ---------------------------------------------------------------- Figure 4
d2 = pd.read_csv(a.p1 / "results" / "unified_contrasts_ds004117.csv")
roi2 = d2[d2.level == "roi"].copy()
ceil22 = 1 / ucs.jzs_bf10(0.0, 22)
fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw={"width_ratios": [1.6, 1]})
bands = [("theta", "Theta\n4–7 Hz"), ("beta", "Beta\n13–30 Hz"), ("beta_low", "Low beta\n13–20 Hz"),
         ("beta_high", "High beta\n21–30 Hz")]
off = {(dc, rf): o for (dc, rf), o in zip([("hilbert", "average"), ("wavelet", "average"),
                                           ("hilbert", "csd"), ("wavelet", "csd")], (-0.27, -0.09, 0.09, 0.27))}
mo = {"PLV": -0.04, "dwPLI": 0.0, "ImC": 0.04}
for xb, (b, lab) in enumerate(bands):
    for _, r in roi2[roi2.band == b].iterrows():
        x = xb + off[(r.decomposition, r.reference)] + mo[r.measure]
        axA.plot(x, r.bf01, MRK[r.measure], ms=5, color=COL[r.measure] if r.reference == "average" else "white",
                 mec=COL[r.measure], mew=1.1)
axA.set_yscale("log"); axA.set_ylim(0.3, 6)
axA.axhline(ceil22, color=INK, lw=0.8, ls=(0, (4, 2)))
axA.text(3.45, ceil22 * 1.04, f"ceiling {ceil22:.2f}", ha="right", va="bottom", fontsize=6.5, color=INK)
for yv, t in ((1, "BF = 1"), (3, "BF$_{01}$ = 3")):
    axA.axhline(yv, color=GRID, lw=0.8, zorder=0)
    axA.text(-0.45, yv * 1.04, t, fontsize=6.5, color=MUTED, va="bottom")
axA.set_xticks(range(4)); axA.set_xticklabels([b[1] for b in bands]); axA.set_xlim(-0.5, 3.5)
from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator
axA.yaxis.set_major_locator(FixedLocator([0.5, 1, 2, 3, 4])); axA.yaxis.set_minor_locator(NullLocator())
axA.set_yticklabels(["0.5", "1", "2", "3", "4"])
axA.set_ylabel("BF$_{01}$ for load effect (log scale)")
axA.set_title("A  Evidence for the null across 12 pipelines", loc="left", fontsize=8.5, fontweight="bold", color=INK)
axA.legend(handles=[Line2D([], [], marker=MRK[m], ls="", color=COL[m], label=m) for m in ("PLV", "dwPLI", "ImC")] +
           [Line2D([], [], marker="o", ls="", color=MUTED, label="filled: average ref."),
            Line2D([], [], marker="o", ls="", mfc="white", mec=MUTED, label="open: CSD")],
           frameon=False, fontsize=6.5, loc="lower right", ncol=2)
h = roi2[roi2.decomposition == "hilbert"].set_index(["band", "measure", "reference"])
w = roi2[roi2.decomposition == "wavelet"].set_index(["band", "measure", "reference"])
lim = (-0.45, 0.6)
axB.plot(lim, lim, color=MUTED, lw=0.8)
for idx in h.index:
    axB.plot(h.loc[idx, "dz"], w.loc[idx, "dz"], MRK[idx[1]], ms=5,
             color=COL[idx[1]] if idx[2] == "average" else "white", mec=COL[idx[1]], mew=1.1)
axB.set_xlim(lim); axB.set_ylim(lim); axB.set_aspect("equal")
axB.set_xlabel("$d_z$, Hilbert phase"); axB.set_ylabel("$d_z$, wavelet phase")
axB.set_title("B  Decomposition dependence", loc="left", fontsize=8.5, fontweight="bold", color=INK)
fig.subplots_adjust(wspace=0.35)
save(fig, "Figure4_dataset2_null_strength")
print("ok; ceiling N=22:", round(ceil22, 4))
