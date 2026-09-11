"""GCF rarefaction curve — how GCF discovery saturates as genomes are added."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# pins svg.hashsalt and selects the Agg backend
from utils.plotting import SVG_METADATA, canonicalise_svg


# region_counts.tsv records carry the genome-file extension (Foo.gbff) because
# count_regions.py names them from the antiSMASH input file; the BiG-SCAPE gbk
# path yields the bare directory name (Foo). Strip only known genbank suffixes —
# Path.stem would eat the version off accessions like GCA_963520565.1.
_GENOME_SUFFIXES = ('.gbff', '.gbk', '.gb', '.genbank')


def strip_genome_suffix(name):
    """Genome name without its genbank file extension, if it has one."""
    lowered = name.lower()
    for suffix in _GENOME_SUFFIXES:
        if lowered.endswith(suffix):
            return name[:-len(suffix)]
    return name


def analysed_genomes(counts_file):
    """Every genome the run analysed, read from region_counts.tsv.

    The BiG-SCAPE database only knows genomes that produced at least one region
    GBK, so it cannot supply this: a genome with no phosphonate BGC is absent
    from `gbk` entirely. Sampling only those genomes answers "as I add more
    BGC-carrying genomes...", not "as I sequence more of this clade..." — the
    question the report's axis label claims to answer.

    Returns a set of genome names, or None if the file is unusable.
    """
    if not counts_file:
        return None
    path = Path(counts_file)
    if not path.exists():
        return None
    try:
        import csv
        with path.open() as f:
            # count_regions.py writes a '#' provenance line above the header
            rows = [ln for ln in f if not ln.startswith('#')]
        names = set()
        for row in csv.DictReader(rows, delimiter='\t'):
            record = (row.get('record') or '').strip()
            if not record:
                continue
            # count_per_contig=true emits "genome|contig"; collapse to the genome,
            # then drop the .gbff so it matches the BiG-SCAPE-derived name
            names.add(strip_genome_suffix(record.split('|')[0]))
        return names or None
    except Exception as e:
        print(f"Warning: could not read genome list from {counts_file}: {e}")
        return None


def chao2(genome_gcfs):
    """Bias-corrected Chao2 richness estimate for incidence data.

    GCFs are "species", genomes are sampling units, and each genome records
    presence/absence — which is exactly the design Chao2 is built for. The
    estimate leans on how much of the diversity was seen only once or twice:

        S_est = S_obs + ((m-1)/m) * Q1^2 / (2*Q2)

    with Q1/Q2 the GCFs found in exactly one/two genomes. Coverage is then
    S_obs/S_est — the fraction of estimated diversity actually recovered.

    Returns dict with s_obs, s_est, coverage (%), q1, q2.
    """
    incidence = {}
    for gcfs in genome_gcfs.values():
        for fam in gcfs:
            incidence[fam] = incidence.get(fam, 0) + 1

    s_obs = len(incidence)
    q1 = sum(1 for c in incidence.values() if c == 1)
    q2 = sum(1 for c in incidence.values() if c == 2)
    m = len(genome_gcfs)

    if s_obs == 0 or m < 2:
        s_est = float(s_obs)
    elif q2 > 0:
        s_est = s_obs + ((m - 1) / m) * (q1 * q1) / (2 * q2)
    elif q1 > 1:
        # no doubletons: fall back to the Q2=0 form rather than dividing by zero
        s_est = s_obs + ((m - 1) / m) * q1 * (q1 - 1) / 2
    else:
        s_est = float(s_obs)

    return {
        's_obs': s_obs,
        's_est': s_est,
        'coverage': (100.0 * s_obs / s_est) if s_est > 0 else 100.0,
        'q1': q1,
        'q2': q2,
    }


def generate_rarefaction_curve(bigscape_db_path, outdir, taxon, n_iterations=50, seed=0,
                               counts_file=None, cutoff=0.3):
    """
    Generate GCF rarefaction curve from BiG-SCAPE database.

    Shows how the number of unique Gene Cluster Families (GCFs) discovered
    increases as more genomes are sampled.

    The genome order is resampled `n_iterations` times from a local RNG seeded
    with `seed`, so the curve and its confidence band are reproducible across
    runs. Pass seed=None for a nondeterministic draw.

    `counts_file` (region_counts.tsv) supplies the full analysed genome set so
    that BGC-free genomes appear on the x-axis as zero-yield draws. Without it
    the curve falls back to BGC-positive genomes only and says so in the plot.

    Returns:
        dict with rarefaction statistics, or None if generation failed.
        'chao2' holds the asymptotic richness estimate and coverage; prefer it
        over 'saturation', which is an uncalibrated slope ratio.
    """
    import os
    import sqlite3
    from collections import defaultdict
    import random

    if not bigscape_db_path or not os.path.exists(bigscape_db_path):
        return None

    def extract_genome_name(path):
        """Extract genome name from BiG-SCAPE GBK path."""
        parts = path.split('/')
        for i, part in enumerate(parts):
            if part == 'antismash_input' and i + 1 < len(parts):
                return parts[i + 1]
        return Path(path).parent.name

    try:
        conn = sqlite3.connect(bigscape_db_path)
        cursor = conn.cursor()

        # Genome -> GCF mapping. Filtered to region records at one cutoff, matching
        # bgc_gcf_tree.py, bgc_all_bgcs_tree.py and viz/report_sections.py — otherwise
        # a multi-value `bigscape_cutoffs` merges family ids from every cutoff into one
        # set and inflates the GCF count. (Measured on Pantoea 2026-08-25 the filters
        # are a no-op: BiG-SCAPE assigns families only to region records, and the run
        # used a single cutoff. They matter only for multi-cutoff runs.)
        cursor.execute("""
            SELECT g.path, bf.family_id
            FROM gbk g
            JOIN bgc_record b ON g.id = b.gbk_id
            JOIN bgc_record_family bf ON b.id = bf.record_id
            JOIN family f ON f.id = bf.family_id
            WHERE b.record_type = 'region' AND f.cutoff = ?
        """, (cutoff,))

        genome_gcfs = defaultdict(set)
        for path, family_id in cursor.fetchall():
            genome = extract_genome_name(path)
            genome_gcfs[genome].add(family_id)

        # Get BGC type counts for top types
        cursor.execute("""
            SELECT product, COUNT(*) as count
            FROM bgc_record
            GROUP BY product
            HAVING count >= 50
            ORDER BY count DESC
            LIMIT 6
        """)
        top_types = [(row[0], row[1]) for row in cursor.fetchall()]

        conn.close()

        n_bgc_positive = len(genome_gcfs)

        # Genomes with no phosphonate BGC never reach the BiG-SCAPE DB. Add them
        # back as empty draws so the x-axis is the analysed set, not the hit set.
        all_genomes = analysed_genomes(counts_file)
        if all_genomes and not (all_genomes & set(genome_gcfs)):
            # Nothing lines up: the two sides name genomes differently. Padding here
            # would silently add every genome twice and inflate the axis, so refuse.
            print(f"Warning: no genome names shared between the BiG-SCAPE database and "
                  f"{counts_file} ({len(genome_gcfs)} vs {len(all_genomes)} names); "
                  f"falling back to BGC-positive genomes only")
            all_genomes = None
        if all_genomes:
            for name in all_genomes:
                genome_gcfs.setdefault(name, set())
            denominator = 'analysed'
        else:
            denominator = 'bgc_positive'

        genomes = sorted(genome_gcfs)   # sorted so the seeded shuffle is stable
        n_genomes = len(genomes)
        total_gcfs = len(set().union(*genome_gcfs.values())) if genome_gcfs else 0

        if n_genomes < 10:
            return None

        richness = chao2(genome_gcfs)

        # Calculate rarefaction curve
        n_points = min(50, n_genomes)
        x_values = sorted(set(np.linspace(1, n_genomes, n_points).astype(int)))

        rng = random.Random(seed)

        all_curves = []
        for _ in range(n_iterations):
            shuffled = genomes.copy()
            rng.shuffle(shuffled)
            seen_gcfs = set()
            curve = []
            for i, genome in enumerate(shuffled, 1):
                seen_gcfs.update(genome_gcfs[genome])
                if i in x_values:
                    curve.append(len(seen_gcfs))
            all_curves.append(curve)

        all_curves = np.array(all_curves)
        mean_gcfs = np.mean(all_curves, axis=0)
        lower_ci = np.percentile(all_curves, 2.5, axis=0)
        upper_ci = np.percentile(all_curves, 97.5, axis=0)

        # Calculate saturation
        if len(mean_gcfs) >= 10:
            early_rate = (mean_gcfs[len(mean_gcfs)//10] - mean_gcfs[0]) / (x_values[len(x_values)//10] - x_values[0]) if x_values[len(x_values)//10] != x_values[0] else 0
            late_rate = (mean_gcfs[-1] - mean_gcfs[-len(mean_gcfs)//10]) / (x_values[-1] - x_values[-len(x_values)//10]) if x_values[-1] != x_values[-len(x_values)//10] else 0
            saturation = (1 - (late_rate / early_rate)) * 100 if early_rate > 0 else 100
        else:
            saturation = 0

        # Plot
        fig, ax = plt.subplots(figsize=(10, 6), facecolor='white')
        ax.set_facecolor('white')

        ax.plot(x_values, mean_gcfs, color='#2c5aa0', linewidth=2.5, label='All BGC types')
        ax.fill_between(x_values, lower_ci, upper_ci, color='#2c5aa0', alpha=0.2)

        ax.set_xlabel('Number of Genomes Sampled'
                      + ('' if denominator == 'analysed' else ' (BGC-positive only)'),
                      fontsize=12)
        ax.set_ylabel('Unique Gene Cluster Families (GCFs)', fontsize=12)
        ax.set_title(f'GCF Rarefaction Curve - {taxon}', fontsize=14)
        ax.grid(True, alpha=0.3)

        # Add annotation. Chao2 coverage leads; the slope-ratio 'saturation' is
        # kept for continuity with earlier reports but is an optimistic heuristic.
        axis_note = (f'{n_bgc_positive} of {n_genomes} carry a BGC'
                     if denominator == 'analysed'
                     else 'BGC-positive genomes only (no counts file)')
        ax.text(0.02, 0.98,
                f'Coverage (Chao2): {richness["coverage"]:.0f}%\n'
                f'{total_gcfs} GCFs observed, {richness["s_est"]:.0f} estimated\n'
                f'{n_genomes} genomes — {axis_note}\n'
                f'Slope-ratio saturation: {saturation:.0f}%',
                transform=ax.transAxes, fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='#e8f4e8', alpha=0.8, edgecolor='#4a9'))

        # Add interpretation guide
        ax.text(0.98, 0.02,
                'Plateau = diversity saturated\nRising = more diversity to discover',
                transform=ax.transAxes, fontsize=9, verticalalignment='bottom',
                ha='right', color='#666')

        plt.tight_layout()
        output_path = Path(outdir) / 'rarefaction_curve.png'
        svg_path    = Path(outdir) / 'rarefaction_curve.svg'
        plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.savefig(svg_path,              bbox_inches='tight', facecolor='white', metadata=SVG_METADATA)
        canonicalise_svg(svg_path)
        plt.close()

        # Embed SVG as base64 for self-contained HTML
        import base64
        with open(svg_path, 'rb') as f:
            rarefaction_svg_b64 = base64.b64encode(f.read()).decode('ascii')

        return {
            'generated': True,
            'n_genomes': n_genomes,
            'n_bgc_positive': n_bgc_positive,
            'denominator': denominator,
            'total_gcfs': total_gcfs,
            'saturation': saturation,
            'chao2': richness,
            'top_types': top_types,
            'svg_b64': rarefaction_svg_b64,
        }

    except Exception as e:
        print(f"Warning: Could not generate rarefaction curve: {e}")
        return None
