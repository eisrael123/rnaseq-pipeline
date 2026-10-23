"""Manuscript-ready style for the pipeline's matplotlib figures (PI review, October 2026).

- Arial everywhere: every journal accepts it. The Docker image gets it from conda-forge's
  mscorefonts; a missing Arial is an error rather than a silent fallback.
- Text about twice matplotlib's defaults, so it stays legible once a figure is shrunk to a
  manuscript panel.
- 1 pt black axis lines and ticks, and 1 pt black outlines around bars.
- Colours from the lab's figures in Nguyen et al., "A BRRF1-CCR4-NOT axis underlies conserved
  transcriptome-wide loss of splicing fidelity during gammaherpesvirus reactivation" (bioRxiv
  2026, doi:10.64898/2026.05.20.726682), sampled from Figures 1B and 6.

deseq2_analysis_ercc.R applies the same rules to the R figures (PCA, ERCC).
"""

# Figure 1B: Neg IncDiff / Pos IncDiff bars.
STEEL_BLUE = "#8CBBCE"
PALE_BLUE = "#CFDBE7"
# Figure 6 bars.
SEA_GREEN = "#5BA493"
PEACH = "#FFC99D"
PINK = "#FFD8E0"
LIGHT_GREY = "#D9D9D9"

FONT = "Arial"
FONT_SIZE = 20       # matplotlib's default is 10
LINE_WIDTH = 1.0     # points
DPI = 300


def apply():
    """Set the style for every figure drawn afterwards in this process."""
    import matplotlib
    from matplotlib import font_manager

    # Raises if Arial is not installed, instead of quietly drawing in DejaVu Sans.
    font_manager.findfont(FONT, fallback_to_default=False)
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [FONT],
        "font.size": FONT_SIZE,
        "axes.titlesize": "large",
        "axes.labelsize": "medium",
        "xtick.labelsize": "medium",
        "ytick.labelsize": "medium",
        "legend.fontsize": "medium",
        "legend.frameon": False,
        "axes.linewidth": LINE_WIDTH,
        "axes.edgecolor": "black",
        "axes.labelcolor": "black",
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": "black",
        "ytick.color": "black",
        "xtick.major.width": LINE_WIDTH,
        "ytick.major.width": LINE_WIDTH,
        "patch.linewidth": LINE_WIDTH,
        "patch.edgecolor": "black",
        "patch.force_edgecolor": True,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        # Keep text as text if a figure is ever saved as PDF/SVG for a journal.
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })
