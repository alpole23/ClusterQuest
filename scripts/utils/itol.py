"""Writers for iTOL (https://itol.embl.de) annotation datasets.

Every iTOL dataset file has the same shape:

    DATASET_<TYPE>
    SEPARATOR TAB
    DATASET_LABEL   <label>
    COLOR           <hex>
    <options>        e.g. STRIP_WIDTH, WIDTH, SHOW_INTERNAL, SIZE_FACTOR
    <field lines>    DATASET_BINARY only: FIELD_SHAPES / FIELD_LABELS / FIELD_COLORS
    <legend lines>   LEGEND_TITLE / LEGEND_SHAPES / LEGEND_COLORS / LEGEND_LABELS
    DATA
    <one row per leaf>

write_dataset() emits exactly that; the helpers below cover the dataset types the
BGC scripts use. Rows are pre-formatted strings (without trailing newline), since
each dataset type has its own column layout.
"""

import os


def write_dataset(path, dataset_type, dataset_label, color, rows,
                  options=None, fields=None, legend=None):
    """Write one iTOL dataset file.

    options : [(key, value), ...]      extra header lines, in order
    fields  : [(shape, label, color), ...]  DATASET_BINARY field definitions
    legend  : (title, [(label, color, shape), ...])
    rows    : iterable of tab-joined data lines
    """
    with open(path, 'w') as f:
        f.write(f'{dataset_type}\n')
        f.write('SEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\t{dataset_label}\n')
        f.write(f'COLOR\t{color}\n')

        for key, value in (options or []):
            f.write(f'{key}\t{value}\n')

        if fields:
            f.write('FIELD_SHAPES\t' + '\t'.join(str(s) for s, _, _ in fields) + '\n')
            f.write('FIELD_LABELS\t' + '\t'.join(l for _, l, _ in fields) + '\n')
            f.write('FIELD_COLORS\t' + '\t'.join(c for _, _, c in fields) + '\n')

        if legend:
            title, items = legend
            f.write(f'LEGEND_TITLE\t{title}\n')
            f.write('LEGEND_SHAPES\t' + '\t'.join(str(s) for _, _, s in items) + '\n')
            f.write('LEGEND_COLORS\t' + '\t'.join(c for _, c, _ in items) + '\n')
            f.write('LEGEND_LABELS\t' + '\t'.join(l for l, _, _ in items) + '\n')

        f.write('DATA\n')
        for row in rows:
            f.write(f'{row}\n')


def write_colorstrip(path, dataset_label, entries, legend=None, color='#777777',
                     options=None):
    """DATASET_COLORSTRIP: one colored band per leaf.

    entries : [(leaf_label, hex_color, display_label), ...]
    """
    write_dataset(
        path, 'DATASET_COLORSTRIP', dataset_label, color,
        rows=(f'{leaf}\t{col}\t{disp}' for leaf, col, disp in entries),
        options=options,
        legend=legend,
    )


def write_binary(path, dataset_label, fields, entries, legend=None, color='#333333'):
    """DATASET_BINARY: presence/absence squares.

    fields  : [(shape, label, color), ...] — one per column
    entries : [(leaf_label, [0/1, ...]), ...]
    """
    write_dataset(
        path, 'DATASET_BINARY', dataset_label, color,
        rows=(f'{leaf}\t' + '\t'.join(str(v) for v in values) for leaf, values in entries),
        fields=fields,
        legend=legend,
    )


def write_simplebar(path, dataset_label, entries, color='#5b5ea6', width=200,
                    show_internal=0):
    """DATASET_SIMPLEBAR: one bar per leaf. entries: [(leaf_label, value), ...]"""
    write_dataset(
        path, 'DATASET_SIMPLEBAR', dataset_label, color,
        rows=(f'{leaf}\t{value}' for leaf, value in entries),
        options=[('WIDTH', width), ('SHOW_INTERNAL', show_internal)],
    )


def write_text(path, dataset_label, entries, color='#333333', options=None,
               position=1, text_color='#333333', style='normal', size_factor=1):
    """DATASET_TEXT: free text beside each leaf. entries: [(leaf_label, text), ...]

    Rows carry iTOL's per-row overrides: position (-1 before / 1 after the leaf),
    color, style and size factor.
    """
    write_dataset(
        path, 'DATASET_TEXT', dataset_label, color,
        rows=(f'{leaf}\t{text}\t{position}\t{text_color}\t{style}\t{size_factor}'
              for leaf, text in entries if text),
        options=options,
    )


def simple_legend(items):
    """Build legend items with iTOL shape 1 (square): [(label, color), ...]."""
    return [(label, color, 1) for label, color in items]


def dataset_path(outdir, bgc_type, name):
    """Conventional output path: {outdir}/{bgc_type}_itol_{name}.txt"""
    return os.path.join(outdir, f'{bgc_type}_itol_{name}.txt')
