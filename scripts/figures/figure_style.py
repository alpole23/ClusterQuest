"""Shared look for the paper figures, and deterministic output.

Every figure script imports this before matplotlib.pyplot, so the salt pinning
in utils/plotting.py is in force and two runs produce byte-identical SVGs.

SVG text is written as outlines by default, which is what makes a figure render
identically on a machine that does not have the font. Set
``CQ_SVG_EDITABLE_TEXT=1`` to emit live <text> elements instead: labels then
stay selectable and editable in Illustrator or Inkscape, which is what a journal
usually wants, at the cost of depending on the font being present. Both are real
vector output -- this only changes whether the glyphs are paths or characters.
"""
import os
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

if os.environ.get('CQ_SVG_EDITABLE_TEXT') not in (None, '', '0'):
    # 'none' leaves glyphs as characters and names the font in the SVG; the
    # default 'path' outlines them. Helvetica/Arial are substituted for DejaVu
    # so a designer opening the file gets a face they actually have.
    plt.rcParams['svg.fonttype'] = 'none'

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


def strip_doctype(path):
    """Drop matplotlib's SVG 1.1 DOCTYPE.

    It is legacy boilerplate no browser validates against, and it makes the file
    unusable wherever DTD machinery is refused -- an XML parser configured
    against external entities rejects the document outright rather than ignoring
    the declaration. Removing it changes nothing about how the figure renders.
    """
    path = Path(path)
    text = path.read_text(encoding='utf-8')
    start = text.find('<!DOCTYPE')
    if start == -1:
        return
    end = text.find('>', start)
    if end == -1:
        return
    path.write_text(text[:start] + text[end + 1:].lstrip('\n'), encoding='utf-8')


def save(fig, outdir, stem):
    """Write <stem>.svg and <stem>.png, both timestamp-free and reproducible."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    svg, png = outdir / f'{stem}.svg', outdir / f'{stem}.png'
    fig.savefig(svg, format='svg', bbox_inches='tight', metadata=SVG_METADATA)
    canonicalise_svg(svg)
    strip_doctype(svg)
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
