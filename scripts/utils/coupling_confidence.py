"""Reference support for coupling enzyme class assignments.

The class itself is assigned from antiSMASH SMCOG/domain markers, which are broad by
design — a virtue here, because characterised phosphonate coupling enzymes are scarce
(Transaminase and Reductase have exactly one reference each). A narrow, reference-driven
classifier would only find enzymes resembling the handful already known, which is the
opposite of the goal.

What this adds is a *support* value: how similar the assigned enzyme is to the nearest
characterised reference of that class. It is advisory. It does not filter, reorder or
downweight anything.

Interpreting a low value needs care, and the ambiguity is irreducible with so few
references:

  - the protein may not be a coupling enzyme of that class at all; or
  - it may be a genuinely novel variant, divergent from the one characterised example.

Both warrant manual inspection, which is the point. Treating low support as "probably
wrong" would discard precisely the candidates worth investigating.

Support is reported against *every* class, not only the assigned one, so an assignment
that is nearly as close to a different class is visible rather than hidden.
"""

from Bio import Align
from Bio.Align import substitution_matrices

# Global alignment, BLOSUM62, affine gaps. Support is reported as percent identity
# rather than a normalised alignment score: unrelated proteins score negative under
# BLOSUM62, which clamps to zero and throws away the low end — exactly the range where
# a divergent, possibly novel enzyme would sit. Percent identity keeps that resolution
# and is directly interpretable ("25% identical to PnaA"). Reconstructing the alignment
# costs ~3.7 ms per pair, so scoring every BGC against every reference is a few seconds.
_ALIGNER = None


def _aligner():
    global _ALIGNER
    if _ALIGNER is None:
        al = Align.PairwiseAligner()
        al.mode = 'global'
        al.open_gap_score = -11
        al.extend_gap_score = -1
        al.substitution_matrix = substitution_matrices.load("BLOSUM62")
        _ALIGNER = al
    return _ALIGNER


def percent_identity(a, b):
    """Percent identity over aligned columns, ignoring gap-only positions."""
    aln = _aligner().align(str(a), str(b))[0]
    top, bottom = aln[0], aln[1]
    ident = sum(1 for x, y in zip(top, bottom) if x == y and x != '-')
    cols = sum(1 for x, y in zip(top, bottom) if x != '-' and y != '-')
    return 100.0 * ident / cols if cols else 0.0


def class_of_reference(description):
    """Coupling class a reference FASTA record belongs to, from its description."""
    d = description.lower()
    if 'malate synthase' in d:
        return 'Synthase'
    if 'decarboxylase' in d:
        return 'Decarboxylase'
    if 'transaminase' in d:
        return 'Transaminase'
    if 'reductase' in d:
        return 'Reductase'
    return None


def genus_species(name):
    """`Streptomyces durhamensis NRRL B-3309` -> `Streptomyces durhamensis`.

    Keeps the binomial and drops the strain designation. The source organism belongs
    beside the score: 23.9% against a *Streptomyces* reference means something quite
    different from 23.9% against a same-genus one, and most characterised phosphonate
    enzymes come from *Streptomyces* while the analysed genomes may not.
    """
    parts = (name or '').split()
    if len(parts) < 2 or not parts[0][:1].isalpha():
        return name or ''
    return f'{parts[0]} {parts[1]}'


def load_references(fasta_path):
    """{class_id: [(name, sequence, organism), ...]} from the reference FASTA.

    Headers are `accession|protein_id|name|description|organism`; `rec.id` stops at the
    first space, so the organism has to come from the full description.
    """
    from Bio import SeqIO
    out = {}
    for rec in SeqIO.parse(str(fasta_path), 'fasta'):
        cls = class_of_reference(rec.description)
        if not cls:
            continue
        fields = rec.description.split('|')
        name = fields[2].strip() if len(fields) > 2 else rec.id
        organism = genus_species(fields[4].strip() if len(fields) > 4 else '')
        out.setdefault(cls, []).append((name, str(rec.seq), organism))
    return out


def support(query_seq, references):
    """Similarity of `query_seq` to every reference class.

    Returns {class_id: {'best_ref', 'pct_id', 'n_refs'}} where `pct_id` is percent
    identity to the closest reference of that class. On the Pantoea genus run a true
    orthologue scored 93.8-100% and unrelated members of the same superfamily scored
    21.8-30.8%, so those are the empirical poles — but with one reference for some
    classes, a low value cannot distinguish "wrong class" from "novel variant".
    """
    out = {}
    for cls, refs in references.items():
        best_name, best_org, best = None, '', 0.0
        for name, seq, organism in refs:
            pid = percent_identity(query_seq, seq)
            if best_name is None or pid > best:
                best_name, best_org, best = name, organism, pid
        out[cls] = {'best_ref': best_name, 'best_ref_organism': best_org,
                    'pct_id': round(best, 1), 'n_refs': len(refs)}
    return out


# Characterised references of *different* classes score 26.7-29.7% against each other
# (HvrC vs Reductase 26.7, PnaA vs Decarboxylase 29.7, VlpB vs Transaminase 28.3,
# Fom2 vs Transaminase 29.1). At or below that, an identity carries no class
# information — it is what unrelated members of the superfamily look like. This is the
# only boundary the current reference set can actually justify; above it there are too
# few references, drawn from too few genera, to calibrate anything.
BACKGROUND_CEILING_PCT = 30.0
