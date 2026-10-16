"""Supplementary Figure S1: PRISMA-ScR flow diagram for the literature audit.

Counts from release_status.md (17 Oct 2026), cross-checked against
wm-phase-connectivity-audit v1.0.0 (screening_sheet_with_D12.csv, deviations D9, D12).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

N = dict(pubmed=546, epmc=310, epmc_recovered=4, dup=308, screened=548,
         excl_ta=418, fulltext=130, excl_ft=29, included=101, sensor=83, source=18)
assert N["pubmed"] + N["epmc"] - N["dup"] == N["screened"]
assert N["screened"] - N["excl_ta"] == N["fulltext"]
assert N["fulltext"] - N["excl_ft"] == N["included"]
assert N["sensor"] + N["source"] == N["included"]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
fig, ax = plt.subplots(figsize=(7.2, 7.6))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def box(x, y, w, h, text, bold=False, fill="white"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="square,pad=0",
                                linewidth=0.9, edgecolor="black", facecolor=fill))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontweight="bold" if bold else "normal", linespacing=1.35)

def arrow(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", lw=0.9, color="black",
                                shrinkA=0, shrinkB=0, mutation_scale=10))

# stage labels
for y, h, lab in [(80, 16, "Identification"), (53, 22, "Screening"), (10, 38, "Included")]:
    ax.add_patch(FancyBboxPatch((0, y), 6, h, boxstyle="square,pad=0", lw=0, facecolor="#d9d9d9"))
    ax.text(3, y + h / 2, lab, rotation=90, ha="center", va="center", fontweight="bold")

MX, MW = 10, 46   # main column
SX, SW = 64, 34   # side column

box(MX, 80, MW, 16,
    f"Records identified\nPubMed (n = {N['pubmed']})\n"
    f"Europe PMC (n = {N['epmc']}, including\n{N['epmc_recovered']} recovered after a retrieval error)")
box(SX, 82, SW, 12, f"Duplicate records removed\n(n = {N['dup']})")
arrow(MX + MW, 88, SX, 88)

box(MX, 64, MW, 11, f"Records screened on title and abstract\n(n = {N['screened']})")
box(SX, 64, SW, 11, f"Records excluded\n(n = {N['excl_ta']})")
arrow(MX + MW / 2, 80, MX + MW / 2, 75)
arrow(MX + MW, 69.5, SX, 69.5)

box(MX, 50, MW, 10, f"Full texts assessed for eligibility\n(n ≥ {N['fulltext']})")
box(SX, 50, SW, 10, f"Excluded at full text\nor at coding\n(n ≥ {N['excl_ft']})")
arrow(MX + MW / 2, 64, MX + MW / 2, 60)
arrow(MX + MW, 55, SX, 55)

box(MX, 33, MW, 11, f"Studies included\n(n = {N['included']})", bold=True)
arrow(MX + MW / 2, 50, MX + MW / 2, 44)

box(MX, 10, 22, 14, f"Sensor space\n(n = {N['sensor']})\nprimary\ndenominator")
box(MX + 24, 10, 22, 14, f"Source space\nonly (n = {N['source']})\nreported\nseparately")
arrow(MX + 11, 33, MX + 11, 24)
arrow(MX + 35, 33, MX + 35, 24)

for ext in ("png", "pdf"):
    fig.savefig(f"FigureS1_PRISMA.{ext}", dpi=300, bbox_inches="tight")
print("ok")
