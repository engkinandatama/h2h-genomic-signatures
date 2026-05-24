"""
visualize_selection.py
======================
Generates a publication-quality Manhattan plot of positive selection signals
across the protein sequence for a given virus_group × protein combination.

Plot layers (top to bottom):
  1. -log10(MEME p-value)      — blue bars (episodic selection)
  2. -log10(FEL p-value)       — orange bars (pervasive selection)
  3. FUBAR posterior prob       — secondary y-axis, purple line
  4. Contrast-FEL significance  — red diamonds on top of significant sites
  5. Significance threshold     — horizontal dashed line at p=0.05 / PP=0.90
  6. PRIME property labels      — annotated text for sites with property bias

Usage:
  python visualize_selection.py \
      --all-sites  results/.../GnGc_all_sites.tsv \
      --out-png    results/.../GnGc_manhattan.png \
      --virus-group Andes_virus \
      --protein GnGc
"""

import argparse
import os
import sys
import csv
import math

try:
    import matplotlib
    matplotlib.use("Agg")  # headless rendering for server
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_sites(path):
    """Load the master per-site TSV from parse_additional_selection.py."""
    sites = []
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return sites
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sites.append(row)
    return sites


def safe_float(val, default=1.0):
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


def neg_log10(p, floor=1e-10):
    """Compute -log10(p) with floor to avoid inf."""
    p = max(float(p), floor)
    return -math.log10(p)


# ---------------------------------------------------------------------------
# Plot generation
# ---------------------------------------------------------------------------

def make_manhattan(sites, virus_group, protein, out_png):
    """Generate and save the Manhattan plot."""
    if not sites:
        print("Warning: no site data, skipping plot.", file=sys.stderr)
        return

    positions  = [int(s["site"]) for s in sites]
    meme_vals  = [neg_log10(safe_float(s.get("meme_p", 1))) for s in sites]
    fel_vals   = [neg_log10(safe_float(s.get("fel_p",  1))) for s in sites]
    fubar_vals = [safe_float(s.get("fubar_pp", 0), default=0.0) for s in sites]

    # Contrast-FEL significant sites
    cfel_sig_pos  = [int(s["site"]) for s in sites if s.get("cfel_sig",  "").lower() == "true"]
    cfel_sig_meme = [neg_log10(safe_float(s.get("meme_p", 1)))
                     for s in sites if s.get("cfel_sig", "").lower() == "true"]

    # Consensus positive selection
    consensus_pos = [int(s["site"]) for s in sites
                     if s.get("consensus_pos", "").lower() == "true"]

    # PRIME significant sites (for annotation)
    prime_sites = [s for s in sites if s.get("prime_sig_properties", "None") != "None"
                   and s.get("prime_sig_properties", "") != ""]

    # Recombination warning
    gard_warn = any(s.get("recombination_warning", "").lower() == "true" for s in sites)

    # -----------------------------------------------------------------------
    # Figure layout
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(
        3, 1,
        figsize=(max(14, len(positions) // 20 + 8), 12),
        gridspec_kw={"height_ratios": [3, 1.2, 1]},
        sharex=True
    )
    fig.patch.set_facecolor("#0f172a")
    for ax in axes:
        ax.set_facecolor("#1e293b")
        ax.spines["bottom"].set_color("#475569")
        ax.spines["left"].set_color("#475569")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(colors="#94a3b8", labelsize=9)

    # -----------------------------------------------------------------------
    # Panel 1: Manhattan (MEME + FEL -log10 p)
    # -----------------------------------------------------------------------
    ax1 = axes[0]
    bar_width = max(0.4, min(2.0, len(positions) // 300 + 0.5))

    ax1.bar(positions, meme_vals, width=bar_width,
            color="#3b82f6", alpha=0.75, label="MEME −log₁₀(p)", zorder=2)
    ax1.bar(positions, fel_vals,  width=bar_width,
            color="#f59e0b", alpha=0.55, label="FEL −log₁₀(p)",  zorder=2)

    # Significance threshold line (p = 0.05 → -log10 = 1.301)
    ax1.axhline(y=-math.log10(0.05), color="#ef4444", linestyle="--",
                linewidth=1.2, alpha=0.8, label="p = 0.05", zorder=3)

    # Mark consensus positive selection sites
    if consensus_pos:
        ax1.scatter(consensus_pos,
                    [max(meme_vals[p - 1], fel_vals[p - 1]) + 0.15 for p in consensus_pos],
                    marker="*", color="#fbbf24", s=60, zorder=4, label="Consensus (+sel)")

    # Mark Contrast-FEL significant sites (H2H > Reservoir)
    if cfel_sig_pos:
        ax1.scatter(cfel_sig_pos, [v + 0.3 for v in cfel_sig_meme],
                    marker="D", color="#f43f5e", s=55, zorder=5,
                    label="Contrast-FEL sig (H2H > Reservoir)")

    # Annotate PRIME sites (only top-N to avoid clutter)
    shown = 0
    for s in prime_sites:
        if shown >= 10:
            break
        site = int(s["site"])
        y    = neg_log10(safe_float(s.get("meme_p", 1)))
        props = s.get("prime_sig_properties", "")
        ax1.annotate(
            props[:12],
            xy=(site, y), xytext=(site, y + 0.6),
            fontsize=6, color="#a78bfa",
            arrowprops=dict(arrowstyle="-", color="#a78bfa", lw=0.6),
        )
        shown += 1

    ax1.set_ylabel("−log₁₀(p-value)", color="#94a3b8", fontsize=10)
    ax1.legend(loc="upper right", fontsize=8, facecolor="#0f172a",
               labelcolor="#cbd5e1", framealpha=0.8)

    title = f"Selection Landscape: {virus_group} — {protein}"
    if gard_warn:
        title += "  ⚠️ Recombination detected (GARD)"
    ax1.set_title(title, color="#f8fafc", fontsize=13, fontweight="bold", pad=12)

    # -----------------------------------------------------------------------
    # Panel 2: FUBAR posterior probability
    # -----------------------------------------------------------------------
    ax2 = axes[1]
    ax2.fill_between(positions, fubar_vals, color="#8b5cf6", alpha=0.4, label="FUBAR PP")
    ax2.plot(positions, fubar_vals, color="#a78bfa", linewidth=0.8, alpha=0.9)
    ax2.axhline(y=0.90, color="#ef4444", linestyle="--", linewidth=1.0, alpha=0.7,
                label="PP = 0.90")
    ax2.set_ylim(0, 1.05)
    ax2.set_ylabel("FUBAR Posterior", color="#94a3b8", fontsize=9)
    ax2.legend(loc="upper right", fontsize=8, facecolor="#0f172a",
               labelcolor="#cbd5e1", framealpha=0.8)

    # -----------------------------------------------------------------------
    # Panel 3: SLAC dN−dS
    # -----------------------------------------------------------------------
    ax3 = axes[2]
    slac_dndS = []
    for s in sites:
        v = s.get("slac_dndS", "NA")
        try:
            slac_dndS.append(float(v))
        except (ValueError, TypeError):
            slac_dndS.append(0.0)

    colors = ["#ef4444" if v > 0 else "#3b82f6" for v in slac_dndS]
    ax3.bar(positions, slac_dndS, width=bar_width, color=colors, alpha=0.7)
    ax3.axhline(y=0, color="#475569", linewidth=0.8)
    ax3.set_ylabel("SLAC dN−dS", color="#94a3b8", fontsize=9)
    ax3.set_xlabel("Codon Position", color="#94a3b8", fontsize=10)

    pos_patch = mpatches.Patch(color="#ef4444", alpha=0.7, label="dN > dS (positive)")
    neg_patch = mpatches.Patch(color="#3b82f6", alpha=0.7, label="dS > dN (purifying)")
    ax3.legend(handles=[pos_patch, neg_patch], loc="upper right",
               fontsize=8, facecolor="#0f172a", labelcolor="#cbd5e1", framealpha=0.8)

    # -----------------------------------------------------------------------
    # Save
    # -----------------------------------------------------------------------
    plt.tight_layout(h_pad=0.5)
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    plt.savefig(out_png, dpi=200, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"Manhattan plot saved to: {out_png}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Generate per-protein dN/dS selection Manhattan plot"
    )
    p.add_argument("--all-sites",   required=True, dest="all_sites",
                   help="Master site TSV from parse_additional_selection.py")
    p.add_argument("--out-png",     required=True, dest="out_png",
                   help="Output PNG file path")
    p.add_argument("--virus-group", required=True, dest="virus_group")
    p.add_argument("--protein",     required=True)
    return p.parse_args()


def main():
    args = parse_args()

    if not HAS_MATPLOTLIB:
        print("Error: matplotlib is not installed. Cannot generate plot.", file=sys.stderr)
        # Create an empty placeholder so Snakemake does not fail the DAG
        os.makedirs(os.path.dirname(args.out_png) or ".", exist_ok=True)
        open(args.out_png, "w").close()
        sys.exit(0)

    sites = load_sites(args.all_sites)
    make_manhattan(sites, args.virus_group, args.protein, args.out_png)


if __name__ == "__main__":
    main()
