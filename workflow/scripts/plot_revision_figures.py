"""
plot_revision_figures.py
========================
Manuscript figures for the revised submission, built from the finished pipeline
outputs so they cannot drift from the tables.

  Figure 1  adaptive-site rate, H2H-capable against spillover-only, by family
  Figure 2  position of differential sites along the protein, by functional role

The role map is imported from analyze_convergence rather than restated here. An
earlier version of this script hardcoded its own, put F_protein in Entry instead
of Entry_Helper, and produced a figure showing 11 entry sites where the pipeline
reported 9.

Usage:
    python workflow/scripts/plot_revision_figures.py <results_dir> <out_dir>
"""
import csv, glob, os, sys, collections, yaml
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import fisher_exact

V   = sys.argv[1]
OUT = sys.argv[2]
BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, GRID = "#1a1a19", "#5c5b55", "#dcdbd4"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.edgecolor": GRID, "axes.linewidth": 0.8,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "savefig.dpi": 300, "savefig.bbox": "tight", "figure.facecolor": "white",
})

cfg = yaml.safe_load(open("config/config.yaml"))
FAM = {g: cfg["virus_groups"][g]["family"] for g in cfg["virus_groups"]}
def klass(g):
    cats = {cfg["viruses"][v]["category"] for v in cfg["virus_groups"][g]["viruses"]}
    return "H2H-capable" if "H2H" in cats else "Spillover-only"

# ---------------------------------------------------------------- data
sites = collections.defaultdict(lambda: [0, 0])
pos   = collections.defaultdict(list)
sys.path.insert(0, "workflow/scripts")
from analyze_convergence import PROTEIN_ROLE as ROLE   # single source of truth
for f in sorted(glob.glob(f"{V}/04_selection/*/*_all_sites.tsv")):
    g = os.path.basename(os.path.dirname(f))
    p = os.path.basename(f)[:-len("_all_sites.tsv")]
    rs = list(csv.DictReader(open(f), delimiter="\t"))
    n  = max(int(r["site"]) for r in rs)
    sel = [r for r in rs if str(r.get("consensus_pos", "")).lower() == "true"]
    sites[g][0] += len(sel); sites[g][1] += len(rs)
    for r in sel:
        pos[ROLE.get(p, "?")].append((int(r["site"]) - 1) / max(n - 1, 1))

FAM_ORDER = ["Hantaviridae", "Paramyxoviridae", "Filoviridae"]
groups = sorted(sites, key=lambda g: (FAM_ORDER.index(FAM[g]), klass(g), g))
rate   = {g: 1000 * sites[g][0] / sites[g][1] for g in groups}

# ---------------------------------------------------------------- Figure 1
fig, ax = plt.subplots(figsize=(6.70, 3.284))   # matches Figure 1 extent
x = np.arange(len(groups)); w = 0.62
cols = [BLUE if klass(g) == "H2H-capable" else ORANGE for g in groups]
bars = ax.bar(x, [rate[g] for g in groups], w, color=cols, linewidth=0, zorder=3)
for xi, g in zip(x, groups):
    ax.annotate(f"{sites[g][0]}", (xi, rate[g]), textcoords="offset points",
                xytext=(0, 4), ha="center", fontsize=8, color=INK2)
ax.set_xticks(x)
ax.set_xticklabels([g.replace("_virus", "").replace("_", " ") for g in groups],
                   rotation=30, ha="right", fontsize=8.5, color=INK)
ax.set_ylabel("Adaptive sites per 1000 codons", color=INK, fontsize=9)
ax.yaxis.grid(True, color=GRID, linewidth=0.7, zorder=0)
ax.set_axisbelow(True)
for s_ in ("top", "right", "left"): ax.spines[s_].set_visible(False)

# Headroom for the family band and the annotations, so nothing overlaps a bar.
top = max(rate.values()) * 1.55
ax.set_ylim(0, top)

bounds, prev = [], None
for i, g in enumerate(groups):
    if prev is not None and FAM[g] != prev: bounds.append(i - 0.5)
    prev = FAM[g]
for b in bounds:
    ax.axvline(b, color=GRID, linewidth=0.9, zorder=1)
seg = [0] + [int(b + 0.5) for b in bounds] + [len(groups)]
for a_, b_, fam in zip(seg[:-1], seg[1:], FAM_ORDER):
    ax.annotate(fam, ((a_ + b_ - 1) / 2, top * 0.965), ha="center", va="top",
                fontsize=8.5, color=INK2, style="italic")
ax.axhline(top * 0.90, color=GRID, linewidth=0.8, zorder=1)

a = [sum(sites[g][i] for g in groups if klass(g) == "H2H-capable") for i in (0, 1)]
b = [sum(sites[g][i] for g in groups if klass(g) == "Spillover-only") for i in (0, 1)]
orr, pv = fisher_exact([[a[0], a[1] - a[0]], [b[0], b[1] - b[0]]])
ax.annotate(f"H2H-capable {1000*a[0]/a[1]:.2f}  vs  spillover-only {1000*b[0]/b[1]:.2f} per 1000 codons\n"
            f"Fisher exact OR = {orr:.2f}, p = {pv:.3f}",
            (0.985, 0.855), xycoords="axes fraction", ha="right", va="top",
            fontsize=8, color=INK2)
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=BLUE),
                   plt.Rectangle((0, 0), 1, 1, color=ORANGE)],
          labels=["H2H-capable", "Spillover-only"], frameon=False, fontsize=8.5,
          loc="upper left", bbox_to_anchor=(0.0, 0.885), labelcolor=INK,
          handlelength=1.1, handleheight=1.1, borderpad=0.0)
fig.savefig(f"{OUT}/Figure1_H2H_vs_spillover_rate.png")
plt.close(fig)

# ---------------------------------------------------------------- Figure 2
fig, ax = plt.subplots(figsize=(6.70, 3.304))   # matches Figure 4 extent
rows = [("Entry", BLUE), ("Replication", ORANGE)]
rng = np.random.default_rng(42)
for i, (role, col) in enumerate(rows):
    v = np.array(pos[role]); y = np.full(len(v), i) + rng.uniform(-0.13, 0.13, len(v))
    ax.scatter(v, y, s=34, color=col, alpha=0.85, linewidth=0.8,
               edgecolor="white", zorder=3)
    m = v.mean()
    ax.plot([m, m], [i - 0.28, i + 0.28], color=col, linewidth=2.4, zorder=4,
            solid_capstyle="round")
    ax.annotate(f"mean {m:.3f}", (m, i + 0.34), ha="center", fontsize=8, color=INK2)
ax.axvline(0.5, color=INK2, linewidth=1.0, linestyle=(0, (4, 3)), zorder=2)
ax.annotate("expected 0.500", (0.5, 1.62), ha="center", fontsize=8, color=INK2)
ax.set_yticks([0, 1])
qv0 = {}
with open(f"{V}/06_statistics/positional_bias_test.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        qv0[r["role"]] = float(r["q_value_across_roles"])
ax.set_yticklabels([f"Entry\n{len(pos['Entry'])} sites\nq = {qv0['Entry']:.3f}",
                    f"Replication\n{len(pos['Replication'])} sites\nq = {qv0['Replication']:.3f}"],
                                          fontsize=8.5, color=INK)
ax.set_ylim(-0.55, 1.75); ax.set_xlim(-0.02, 1.02)
ax.set_xlabel("Relative position along protein  (0 = N-terminus, 1 = C-terminus)",
              color=INK, fontsize=9)
ax.xaxis.grid(True, color=GRID, linewidth=0.7, zorder=0); ax.set_axisbelow(True)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
nh = len(pos.get("Entry_Helper", []))
if nh:
    ax.annotate(f"Entry_Helper (F protein): {nh} site(s), too few to test",
                (0.99, -0.42), ha="right", fontsize=7.5, color=INK2)
fig.savefig(f"{OUT}/Figure2_positional_bias.png")
plt.close(fig)
print("Figure 1 dan 2 ditulis ke", OUT)
for role in ("Entry","Replication"):
    print(f"   {role:12s} n={len(pos[role]):2d} mean={np.mean(pos[role]):.3f}")


# ---------------------------------------------------------------- Figure PRIME
# Replaces the manuscript's Figure 2: the physicochemical selection profile at
# NiV-B L-protein site 210. Lambda and p come from the PRIME JSON directly, so
# they cannot drift from the table.
import json
pj = json.load(open(f"{V}/04_selection/Nipah_NiVB/L_protein_prime.json"))
row = pj["MLE"]["content"]["0"][209]          # zero-based index of codon 210
props = [("Hydrophobicity", row[12], row[13]),
         ("Isoelectric point", row[15], row[16]),
         ("Volume", row[18], row[19])]
omni_p, omni_q = row[9], row[10]

fig, ax = plt.subplots(figsize=(6.70, 3.634))
y = np.arange(len(props))[::-1]
for yi, (name, lam, pv) in zip(y, props):
    sig = pv < 0.05
    col = (BLUE if lam > 0 else ORANGE) if sig else "#b9b8b0"
    ax.barh(yi, lam, height=0.5, color=col, linewidth=0, zorder=3)
    # Always annotate to the right of zero. A label trailing off a negative bar
    # runs into the property name on the axis.
    x = lam + 0.25 if lam > 0 else 0.25
    ax.annotate(f"lambda = {lam:+.2f}   p = {pv:.3f}", (x, yi), va="center",
                ha="left", fontsize=8.5, color=INK if sig else INK2)
ax.axvline(0, color=INK2, linewidth=1.0, zorder=4)
ax.set_yticks(y); ax.set_yticklabels([p[0] for p in props], fontsize=9, color=INK)
ax.set_xlabel("PRIME lambda   (positive = property conserved, negative = property diversifying)",
              fontsize=9, color=INK)
lim = max(abs(p[1]) for p in props) * 2.05
ax.set_xlim(-lim, lim); ax.set_ylim(-0.6, len(props) - 0.4)
ax.xaxis.grid(True, color=GRID, linewidth=0.7, zorder=0); ax.set_axisbelow(True)
for s_ in ("top", "right", "left"): ax.spines[s_].set_visible(False)
ax.annotate(f"NiV-B L-protein site 210   omnibus p = {omni_p:.3f}, q = {omni_q:.2f}",
            (0.5, 1.02), xycoords="axes fraction", ha="center", va="bottom",
            fontsize=8.5, color=INK2)
ax.annotate("grey = not significant", (0.99, 0.02), xycoords="axes fraction",
            ha="right", fontsize=7.5, color=INK2)
fig.savefig(f"{OUT}/Figure2_PRIME_L210.png")
plt.close(fig)
print("Figure PRIME (situs 210) ditulis")
for n, l, pv in props:
    print(f"   {n:20s} lambda={l:+.4f}  p={pv:.6f}")
