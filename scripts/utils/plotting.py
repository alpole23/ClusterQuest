"""Deterministic matplotlib output — import before saving any figure.

Matplotlib gives SVGs two run-varying fields: a `<dc:date>` creation timestamp
and per-element `id` attributes drawn from a random salt. Both leak into
`bgc_report.html`, which embeds several SVGs as base64, so two runs over
identical data produced byte-different reports. Importing this module pins the
salt; pass `metadata=SVG_METADATA` to `savefig` to drop the timestamp.
"""

import re
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt'] = 'clusterquest'

# savefig(..., metadata=SVG_METADATA) — omits the creation timestamp.
# SVG only; the PNG writer takes a different metadata dict.
SVG_METADATA = {'Date': None}

# matplotlib clip-path ids look like id="p1a2b3c4d5e" / url(#p1a2b3c4d5e)
_CLIP_ID = re.compile(r'(id="|url\(#)p([0-9a-f]{10})')


def canonicalise_svg(path):
    """Rewrite matplotlib's clip-path ids to stable, appearance-ordered names.

    `svg.hashsalt` pins marker and hatch ids but NOT clip-path ids: matplotlib
    keys those on `(id(clippath), str(clippath_trans))` (backend_svg.py), and
    `id()` is a CPython memory address. Addresses do repeat when the allocation
    sequence is identical — which is why simple figures compare equal by luck —
    but any change in allocation order shifts them, so it is not determinism to
    rely on. Renaming in order of first appearance makes the bytes stable
    regardless.

    Rendering is untouched: these ids are internal document references, and both
    the definition and every reference are rewritten together.
    """
    path = Path(path)
    try:
        text = path.read_text()
    except (OSError, UnicodeDecodeError):
        return path

    mapping = {}
    for _, digest in _CLIP_ID.findall(text):
        if digest not in mapping:
            mapping[digest] = f'{len(mapping):010d}'
    if not mapping:
        return path

    path.write_text(_CLIP_ID.sub(lambda m: f'{m.group(1)}p{mapping[m.group(2)]}', text))
    return path
