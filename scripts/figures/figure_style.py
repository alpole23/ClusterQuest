"""Shared look for the paper figures, and deterministic output.

Every figure script imports this before matplotlib.pyplot, so the salt pinning
in utils/plotting.py is in force and two runs produce byte-identical SVGs.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.plotting import SVG_METADATA, canonicalise_svg  # noqa: E402  (pins the salt)

import matplotlib.pyplot as plt  # noqa: E402

# Colourblind-safe. OFF/BEFORE is the muted one, ON/AFTER the saturated one,
# so the figure reads correctly in greyscale as well.
BEFORE = '#9aa5b1'
AFTER = '#1b5e7e'
ACCENT = '#8d3a2c'
GOOD = '#2f7d5d'
INK = '#15181d'
FAINT = '#6d7683'

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 9,
    'axes.titlesize': 10,
    'axes.titleweight': 'bold',
    'axes.labelsize': 9,
    'axes.edgecolor': '#c9ced6',
    'axes.linewidth': 0.8,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.color': INK,
    'ytick.color': INK,
    'xtick.labelsize': 8.5,
    'ytick.labelsize': 8.5,
    'legend.frameon': False,
    'legend.fontsize': 8.5,
    'figure.dpi': 150,
})


def save(fig, outdir, stem):
    """Write <stem>.svg and <stem>.png, both timestamp-free and reproducible."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    svg, png = outdir / f'{stem}.svg', outdir / f'{stem}.png'
    fig.savefig(svg, format='svg', bbox_inches='tight', metadata=SVG_METADATA)
    canonicalise_svg(svg)
    fig.savefig(png, format='png', bbox_inches='tight', dpi=300)
    print(f'wrote {svg}\nwrote {png}')
    return svg, png


def panel_label(ax, letter, dx=-0.16, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=12,
            fontweight='bold', va='top', ha='left', color=INK)


def bar_values(ax, bars, fmt='{:.0f}', dy=0.01, **kw):
    """Print each bar's value just above it, in axis-relative offset."""
    span = ax.get_ylim()[1] - ax.get_ylim()[0]
    for b in bars:
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + dy * span,
                fmt.format(b.get_height()), ha='center', va='bottom',
                fontsize=7.8, color=FAINT, **kw)
