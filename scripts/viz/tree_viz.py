#!/usr/bin/env python3
"""Prune a GTDB-Tk phylogenetic tree to the genomes analysed in this run.

This file used to be 1,137 lines: a hand-rolled Newick parser, three circular-tree
plotters, and this function. The plotters were dead — exported from ``viz/__init__`` but
called from nowhere — and ``plot_circular_phylogenetic_tree`` timed out after two minutes
on a real GTDB-Tk tree, because it drove the hand-rolled parser rather than Bio.Phylo.
The parser itself survived only as an ``except`` fallback here, unreachable in practice
since biopython is a hard dependency of every conda environment that runs this code.

What remains is the Bio.Phylo path, which prunes the same tree in about two seconds.
Its ``newick`` payload is not rendered anywhere — the report has no tree renderer — but
the pruned Newick it writes to ``pruned_phylo_tree.nwk`` is a published output, used with
iTOL, FigTree or Dendroscope.
"""

import sys
from io import StringIO
from pathlib import Path

import pandas as pd
from Bio import Phylo

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


def prepare_phylo_tree_for_js(newick_file, gtdbtk_summary, counts_file, outdir, outgroup=None):
    """Prune the GTDB-Tk tree to the genomes analysed in this run.

    Writes ``pruned_phylo_tree.nwk`` beside the report, which is the useful output: it
    opens in iTOL, FigTree or Dendroscope. The returned dict is also handed to the report
    generator, but nothing there renders it — the report has no tree renderer.

    Args:
        newick_file: Path to Newick tree file from GTDB-Tk
        gtdbtk_summary: Path to GTDB-Tk summary TSV
        counts_file: Path to region_counts.tsv for BGC data
        outdir: Output directory
        outgroup: Optional outgroup taxon pattern (e.g., "g__Escherichia") to keep in pruned tree

    Returns:
        dict with 'newick' (pruned tree string) and 'metadata' (genome info), or None on failure
    """
    newick_path = Path(newick_file)
    if not newick_path.exists():
        print(f"Newick file not found: {newick_file}")
        return None

    # Read user genomes from GTDB-Tk summary
    user_genomes = set()
    genome_metadata = {}
    if gtdbtk_summary and Path(gtdbtk_summary).exists():
        try:
            gtdbtk_df = pd.read_csv(gtdbtk_summary, sep='\t')
            if 'user_genome' in gtdbtk_df.columns:
                for _, row in gtdbtk_df.iterrows():
                    user_genome = row['user_genome']
                    user_genomes.add(user_genome)
                    genome_metadata[user_genome] = {
                        'classification': row.get('classification', ''),
                    }
            print(f"Found {len(user_genomes)} user genomes for tree pruning")
        except Exception as e:
            print(f"Warning: Could not read GTDB-Tk summary: {e}")
            return None

    if not user_genomes:
        print("No user genomes found - cannot create tree")
        return None

    # Read BGC counts
    if counts_file and Path(counts_file).exists():
        try:
            df = pd.read_csv(counts_file, sep='\t', skiprows=lambda i: i == 0)
            for _, row in df.iterrows():
                genome_name = row.get('genome', row.get('file', ''))
                if genome_name:
                    clean_name = Path(genome_name).stem if '/' in str(genome_name) else genome_name
                    clean_name = clean_name.replace('.gbff', '').replace('.fna', '')
                    if clean_name in genome_metadata:
                        genome_metadata[clean_name]['total_bgcs'] = int(row.get('total_count', 0))
        except Exception as e:
            print(f"Warning: Could not read counts file: {e}")

    # Try to use Bio.Phylo for efficient parsing
    from Bio import Phylo
    import sys

    # Increase recursion limit for large trees (GTDB reference trees can be very deep)
    old_limit = sys.getrecursionlimit()
    sys.setrecursionlimit(max(old_limit, 15000))

    print(f"Reading phylogenetic tree from {newick_file} using Bio.Phylo...")

    # Parse tree with Bio.Phylo (much faster than custom parser)
    tree = Phylo.read(newick_file, 'newick')

    # Get all terminal names (leaves)
    terminals = tree.get_terminals()
    print(f"Tree has {len(terminals)} leaves")

    # Find user genomes in tree
    user_terminals = [t for t in terminals if t.name in user_genomes]
    print(f"Found {len(user_terminals)} user genomes in tree")

    if not user_terminals:
        print("No user genomes found in tree - check genome name matching")
        # Try partial matching
        for t in terminals[:10]:
            print(f"  Sample terminal: {t.name}")
        return None

    # Find outgroup terminal if specified
    outgroup_terminal = None
    if outgroup:
        # Search for a terminal matching the outgroup pattern
        # Pattern can be taxonomy rank (e.g., "g__Escherichia") or partial name
        for t in terminals:
            if t.name and outgroup in t.name:
                outgroup_terminal = t
                print(f"Found outgroup: {t.name} (matching pattern '{outgroup}')")
                break
        if not outgroup_terminal:
            print(f"Warning: No outgroup found matching pattern '{outgroup}'")

    # Prune tree to user genomes (and outgroup if found)
    terminals_to_keep = len(user_terminals) + (1 if outgroup_terminal else 0)
    print(f"Pruning tree to {terminals_to_keep} terminals...")

    # For small numbers of genomes, build a simple tree with just their relationships
    if len(user_terminals) == 1 and not outgroup_terminal:
        # Single genome without outgroup - create simple tree
        name = user_terminals[0].name
        bl = user_terminals[0].branch_length or 0.0
        pruned_newick = f"({name}:{bl});"
    else:
        # Use Bio.Phylo's built-in pruning - O(n) instead of O(n²) distance matrix
        # Strategy: remove all terminals except user genomes and outgroup
        terminals_to_keep = set(t.name for t in user_terminals)
        if outgroup_terminal:
            terminals_to_keep.add(outgroup_terminal.name)
        non_user_terminals = [t for t in terminals if t.name not in terminals_to_keep]

        print(f"Removing {len(non_user_terminals)} terminals from tree...")

        # Remove non-user terminals in batches with progress reporting
        removed_count = 0
        total_to_remove = len(non_user_terminals)
        report_interval = max(1, total_to_remove // 20)  # Report ~20 times

        for terminal in non_user_terminals:
            try:
                tree.prune(terminal)
                removed_count += 1
                if removed_count % report_interval == 0:
                    print(f"  Pruning progress: {removed_count}/{total_to_remove} ({100*removed_count//total_to_remove}%)")
            except Exception as e:
                # Terminal may already be removed if it was part of a collapsed branch
                pass

        print(f"Removed {removed_count} terminals")

        # Collapse single-child internal nodes to clean up the tree
        def collapse_single_child_clades(clade):
            """Recursively collapse internal nodes with single children."""
            if clade.is_terminal():
                return clade

            # Process children first
            new_clades = []
            for child in clade.clades:
                collapsed_child = collapse_single_child_clades(child)
                if collapsed_child is not None:
                    new_clades.append(collapsed_child)

            clade.clades = new_clades

            # If this node has only one child, merge branch lengths
            if len(clade.clades) == 1:
                child = clade.clades[0]
                # Add this node's branch length to child
                if clade.branch_length and child.branch_length:
                    child.branch_length += clade.branch_length
                elif clade.branch_length:
                    child.branch_length = clade.branch_length
                return child

            # If no children left, return None
            if len(clade.clades) == 0:
                return None

            return clade

        print("Collapsing single-child internal nodes...")
        tree.root = collapse_single_child_clades(tree.root)

        # Write to Newick
        output = StringIO()
        Phylo.write(tree, output, 'newick')
        pruned_newick = output.getvalue().strip()

    print(f"Pruned tree newick length: {len(pruned_newick)} chars")

    # Save pruned Newick to file
    pruned_newick_path = Path(outdir) / 'pruned_phylo_tree.nwk'
    with open(pruned_newick_path, 'w') as f:
        f.write(pruned_newick)
    print(f"Saved pruned Newick to {pruned_newick_path}")

    # Count leaves in pruned tree
    pruned_tree = Phylo.read(StringIO(pruned_newick), 'newick')
    pruned_leaf_count = len(pruned_tree.get_terminals())

    return {
        'newick': pruned_newick,
        'metadata': genome_metadata,
        'leaf_count': pruned_leaf_count
    }

