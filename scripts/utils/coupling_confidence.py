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


# Classes that share another class's references. The CDP-activating variant is still a
# phosphonopyruvate decarboxylase — it just additionally encodes a cytidylyltransferase —
# so it is scored against the same DhpF/Fom2/Ppd set rather than against nothing.
_SHARED_REFS = {'Decarboxylase-Nucleotidyltransferase': 'Decarboxylase'}


def load_references(fasta_path):
    """{class_id: [(name, sequence), ...]} from the coupling enzyme reference FASTA."""
    from Bio import SeqIO
    out = {}
    for rec in SeqIO.parse(str(fasta_path), 'fasta'):
        cls = class_of_reference(rec.description)
        if not cls:
            continue
        name = rec.id.split('|')[2] if rec.id.count('|') >= 2 else rec.id
        out.setdefault(cls, []).append((name, str(rec.seq)))
    for derived, source in _SHARED_REFS.items():
        if source in out:
            out.setdefault(derived, list(out[source]))
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
        best_name, best = None, 0.0
        for name, seq in refs:
            pid = percent_identity(query_seq, seq)
            if best_name is None or pid > best:
                best_name, best = name, pid
        out[cls] = {'best_ref': best_name, 'pct_id': round(best, 1),
                    'n_refs': len(refs)}
    return out
