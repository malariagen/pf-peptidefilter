"""
Supplementary figure: how far each benchmark antigen sits from each Pf-PeptideFilter
threshold.

Replaces a retained/excluded matrix, which for these eight genes is uniformly
"excluded" and therefore carries no information. Distance to the requirement does:
it separates an antigen that narrowly misses from one that misses by a wide margin.

Encoding notes
  - panels A, C and D plot a *slider setting* -- the highest "% of the gene's peptides
    that must pass" the gene still survives -- not the exact fraction of its peptides that
    pass. The app's gene rule tests an interpolated percentile of the peptide metric, so
    the two differ by a few points in either direction (e.g. CelTOS conservation: 40 on
    the slider, 50.0% of peptides actually at >=0.99). Panel B is exact to the 0.01 grid.
  - position carries the value; the dashed rule is the requirement. Colour encodes
    nothing, so the figure survives greyscale printing and any form of colour blindness.
  - a filled marker passes that filter, a hollow marker fails it - a second,
    non-colour channel, so status is never carried by colour alone.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
from matplotlib.lines import Line2D      # noqa: E402
import pandas as pd                      # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_sensitivity import OUT          # noqa: E402

# tokens (contrast against the surface computed, not eyeballed:
# series 4.30:1, muted ink 3.50:1, secondary ink 7.73:1)
SURFACE   = "#fcfcfb"
SERIES    = "#2a78d6"
INK       = "#0b0b0b"
INK_2     = "#52514e"
MUTED     = "#898781"
GRID      = "#e1e0d9"
BASELINE  = "#c3c2b7"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.facecolor": SURFACE, "figure.facecolor": SURFACE,
    "axes.edgecolor": BASELINE, "axes.labelcolor": INK_2,
    "xtick.color": MUTED, "ytick.color": INK_2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False,
})


def dot_panel(ax, labels, values, requirement, title, xlabel, xlim, passes,
              fmt="{:.0f}", req_label=None, label_side="right"):
    y = range(len(labels))
    ax.axvline(requirement, color=MUTED, lw=1.2, ls=(0, (4, 3)), zorder=1)
    for yi, (val, ok) in enumerate(zip(values, passes)):
        if pd.isna(val):
            ax.annotate("below the lowest setting the app offers", (xlim[0], yi),
                        xytext=(4, 0), textcoords="offset points", va="center",
                        fontsize=6.5, color=MUTED, style="italic")
            continue
        ax.plot([xlim[0], val], [yi, yi], color=SERIES, lw=2, alpha=0.35,
                solid_capstyle="round", zorder=2)
        ax.plot([val], [yi], marker="o", ms=8, zorder=3,
                color=SERIES if ok else SURFACE,
                markeredgecolor=SERIES, markeredgewidth=1.8)
        dx, ha = (11, "left") if label_side == "right" else (-11, "right")
        ax.annotate(fmt.format(val), (val, yi), xytext=(dx, 0),
                    textcoords="offset points", va="center", ha=ha,
                    fontsize=7.5, color=INK_2)
    ax.set_yticks(list(y)); ax.set_yticklabels(labels, fontsize=8)
    ax.set_ylim(-0.7, len(labels) - 0.3); ax.invert_yaxis()
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel, fontsize=7.5)
    ax.set_title(title, fontsize=8.5, color=INK, loc="left", pad=16, fontweight="bold")
    ax.annotate(req_label or f"requirement {fmt.format(requirement)}",
                (requirement, -0.72), xytext=(0, 3), textcoords="offset points",
                ha="center", va="bottom", fontsize=6.8, color=MUTED,
                annotation_clip=False)
    ax.xaxis.grid(True, color=GRID, lw=0.7); ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)


def main():
    t = pd.read_csv(os.path.join(OUT, "supplementary_table_B_antigens.csv"))
    t = t.sort_values("conservation_gene_rule_pct_tolerated", ascending=False)
    lab = t.antigen.tolist()

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.4))

    dot_panel(axes[0][0], lab, t.conservation_gene_rule_pct_tolerated, 95,
              "A  Strain conservation",
              "highest \"% of peptides must pass\" setting still met, at ≥99% frequency",
              (0, 108), t.conservation.eq("PASS"), req_label="requirement 95%")

    dot_panel(axes[0][1], lab, t.conservation_max_threshold_tolerated, 0.99,
              "B  Conservation threshold tolerated",
              "highest haplotype frequency the gene still meets", (0, 1.13),
              t.conservation.eq("PASS"), fmt="{:.2f}", req_label="requirement 0.99",
              label_side="left")

    dot_panel(axes[1][0], lab, t.human_gene_rule_pct_tolerated, 75,
              "C  Human similarity",
              "highest \"% of peptides must pass\" setting still met, identity and length",
              (0, 108), t.human.eq("PASS"), req_label="requirement 75%")

    dot_panel(axes[1][1], lab, t.indel_gene_rule_pct_tolerated, 100,
              "D  Indel frequency",
              "highest \"% of peptides must pass\" setting still met, at ≤5% indel frequency",
              (0, 108), t.indel.eq("PASS"), req_label="requirement 100%",
              label_side="left")

    handles = [
        Line2D([], [], marker="o", ls="", ms=8, color=SERIES,
               markeredgecolor=SERIES, markeredgewidth=1.8, label="passes this filter"),
        Line2D([], [], marker="o", ls="", ms=8, color=SURFACE,
               markeredgecolor=SERIES, markeredgewidth=1.8, label="fails this filter"),
        Line2D([], [], color=MUTED, lw=1.2, ls=(0, (4, 3)),
               label="requirement under the published configuration"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=7.5, bbox_to_anchor=(0.5, -0.012), labelcolor=INK_2)

    fig.suptitle("Behaviour of eight established malaria vaccine antigens under "
                 "Pf-PeptideFilter's criteria",
                 fontsize=10, fontweight="bold", color=INK, x=0.011, ha="left", y=1.035)
    fig.text(0.011, 0.952,
             "All eight are excluded by the published configuration. Conservation is the "
             "filter responsible in every case (panels A and B).\n"
             "Panels A, C and D are settings of the app's own gene-rule slider, not exact "
             "peptide fractions: the rule tests an interpolated percentile, so a value here can "
             "sit a few points either side of the true share of peptides passing.",
             fontsize=7.8, color=INK_2, ha="left", va="top", linespacing=1.5)

    fig.tight_layout(rect=[0, 0.035, 1, 0.88])
    fig.subplots_adjust(hspace=0.62, wspace=0.42)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"supplementary_figure_antigen_thresholds.{ext}"),
                    dpi=300, facecolor=SURFACE, bbox_inches="tight")
    print("wrote supplementary_figure_antigen_thresholds.png / .pdf")


if __name__ == "__main__":
    main()
