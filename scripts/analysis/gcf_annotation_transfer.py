#!/usr/bin/env python3
"""Complete each BGC's annotation from its better-annotated relatives in the same GCF.

Half this dataset is unreadable at the gene level. Measured on Erwiniaceae: only
41.5% of CDS carry an informative product, and 170 of 333 regions carry *none* —
not because the clusters are bare but because the assemblies are GenBank-only, with
no functional annotation at all. Any gene-content metric built on those products is
therefore comparing NCBI annotation pipelines rather than biology.

The members of a gene cluster family are homologous by construction, and annotation
quality is distributed very unevenly across them: the typical family here has a
*median* member at 0% and a *best* member at 88-100%. One RefSeq-quality genome
carries the whole family. So orthologues are grouped inside each family and a
product name found on any member is propagated to the rest, taking coverage from
41.5% to ~83% and making 16 of 19 families comparable.

What this is not: an observation. A transferred product is an inference from a
homologue, and every row records where it came from — source genome, identity to
that source, how many members of the group carried the same name, and how many
disagreed. The failure mode is error propagation: one mis-annotated RefSeq gene
becomes N mis-annotated genes, and the agreement count looks reassuring because
they all descend from the same source. `n_sources` is the column that exposes that
— a call supported by one genome is a very different thing from one supported by
twelve independent assemblies.

Three families here gain nothing: singletons with no annotated relative to learn
from (GCF-8, 13, 16 on Erwiniaceae). They are reported with transferred == original
so that downstream code can exclude them rather than silently treating an absent
annotation as an absent gene.
"""
import argparse
import collections
import csv
import json
import re
import sqlite3
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils import domain_functions

# A product that names a function. BiG-SCAPE stores no products at all, so these
# come from the region GenBanks antiSMASH wrote, where an unannotated assembly
# leaves either an empty string or a placeholder.
UNINFORMATIVE = re.compile(
    r'^\s*$|hypothetical|unknown|uncharacteri[sz]ed|^-$|^putative protein$|'
    r'^conserved protein$|^DUF\d+', re.I)


def informative(product):
    return bool(product) and not UNINFORMATIVE.match(product.strip())


def family_members(db_path, cutoff):
    """{family_id: [(genome, region_basename)]} at the given BiG-SCAPE cutoff.

    `family.id` is an AUTOINCREMENT column recording write order, not a stable
    biological identity — it shifts between runs. It is used here only to group
    members within this run, never persisted as a label.
    """
    db = sqlite3.connect(db_path)
    rows = db.execute(
        """SELECT bf.family_id, g.path
           FROM bgc_record_family bf
           JOIN bgc_record r ON r.id = bf.record_id
           JOIN gbk g        ON g.id = r.gbk_id
           JOIN family f     ON f.id = bf.family_id
           WHERE f.cutoff = ?""", (cutoff,)).fetchall()
    db.close()
    fams = collections.defaultdict(list)
    for fam_id, path in rows:
        p = Path(path)
        fams[fam_id].append((p.parent.name, p.name))
    return fams


def read_region(gbk_path):
    """[(locus_tag, product, translation, [pfam_acc])] for one region GenBank.

    Domains come from antiSMASH's own clusterhmmer scan (`PFAM_domain` features), so
    unlike the product they are present whether or not NCBI annotated the assembly.
    That is what makes them usable for functional comparison across families; see
    utils/domain_functions.py.
    """
    from Bio import SeqIO
    out, doms = [], collections.defaultdict(list)
    for rec in SeqIO.parse(str(gbk_path), 'genbank'):
        for feat in rec.features:
            if feat.type != 'PFAM_domain':
                continue
            q = feat.qualifiers
            tag = (q.get('locus_tag') or ['?'])[0]
            for ref in q.get('db_xref', []):
                if ref.startswith('PF'):
                    doms[tag].append(ref.split('.')[0])
    for rec in SeqIO.parse(str(gbk_path), 'genbank'):
        for feat in rec.features:
            if feat.type != 'CDS':
                continue
            q = feat.qualifiers
            tag = (q.get('locus_tag') or q.get('gene') or ['?'])[0]
            prod = (q.get('product') or [''])[0]
            seq = (q.get('translation') or [''])[0]
            if seq:
                out.append((tag, prod, seq, doms.get(tag, [])))
    return out


def diamond_edges(binary, fasta, workdir, min_id, min_cov, threads):
    """All-vs-all within one family; yields (qi, si, pident) above the cuts.

    Mutual coverage, not one-sided: a short fragment aligning cleanly inside a long
    multidomain protein is not an orthologue, and one-sided coverage would accept it.
    """
    dbp = workdir / 'fam'
    subprocess.run([binary, 'makedb', '--in', str(fasta), '-d', str(dbp), '--quiet'],
                   check=True, capture_output=True)
    proc = subprocess.run(
        [binary, 'blastp', '-d', str(dbp), '-q', str(fasta), '--quiet',
         '--threads', str(threads), '--max-target-seqs', '0',
         '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen'],
        capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f'diamond blastp failed:\n{proc.stderr[:2000]}')
    for line in proc.stdout.splitlines():
        q, s, pid, ln, ql, sl = line.split('\t')
        if q == s:
            continue
        pid, ln, ql, sl = float(pid), int(ln), int(ql), int(sl)
        if pid < min_id or ln / ql < min_cov or ln / sl < min_cov:
            continue
        yield int(q), int(s), pid


def group_orthologues(n, edges):
    """Single-linkage components. Returns {member_index: group_id}."""
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    best_id = {}
    for a, b, pid in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
        key = (min(a, b), max(a, b))
        best_id[key] = max(best_id.get(key, 0.0), pid)
    return {i: find(i) for i in range(n)}, best_id


def consensus(names):
    """Majority product among annotated members; ties broken by the longer name.

    The longer name is the more specific one far more often than not
    ("phosphonopyruvate decarboxylase" over "decarboxylase"), and a tie means the
    counts gave us nothing to go on.
    """
    counts = collections.Counter(names)
    top = max(counts.values())
    return sorted((n for n, c in counts.items() if c == top), key=len)[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--db', type=Path, required=True, help='BiG-SCAPE sqlite database')
    ap.add_argument('--antismash', type=Path, required=True,
                    help='antiSMASH results dir holding <genome>/<region>.gbk')
    ap.add_argument('--outdir', type=Path, default=Path('.'))
    ap.add_argument('--cutoff', type=float, default=0.3)
    ap.add_argument('--min_identity', type=float, default=50.0,
                    help='percent identity floor for calling two CDS orthologous')
    ap.add_argument('--min_coverage', type=float, default=0.70,
                    help='mutual alignment coverage floor')
    ap.add_argument('--diamond', default='diamond')
    ap.add_argument('--threads', type=int, default=4)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    fams = family_members(args.db, args.cutoff)
    if not fams:
        sys.exit(f'no families at cutoff {args.cutoff} in {args.db}')

    per_cds, per_group, summary = [], [], {}
    work = args.outdir / '_transfer'
    work.mkdir(exist_ok=True)

    for fam_id in sorted(fams):
        members = fams[fam_id]
        cds = []                      # (region_label, genome, tag, product, seq)
        for genome, region_file in members:
            path = args.antismash / genome / region_file
            if not path.exists():
                print(f'  warn: missing {path}', file=sys.stderr)
                continue
            label = region_file[:-4] if region_file.endswith('.gbk') else region_file
            for tag, prod, seq, dm in read_region(path):
                cds.append((label, genome, tag, prod, seq, dm))
        if not cds:
            continue

        before = sum(1 for c in cds if informative(c[3]))

        fasta = work / f'fam{fam_id}.faa'
        with fasta.open('w') as fh:
            for i, (_, _, _, _, seq, _) in enumerate(cds):
                fh.write(f'>{i}\n{seq}\n')

        # A family with one member has no relative to learn from; skip the search
        # rather than paying for an all-vs-all that can only find self-hits.
        if len(members) > 1:
            edges = list(diamond_edges(args.diamond, fasta, work,
                                       args.min_identity, args.min_coverage,
                                       args.threads))
        else:
            edges = []
        gid, best_id = group_orthologues(len(cds), edges)

        groups = collections.defaultdict(list)
        for i, g in gid.items():
            groups[g].append(i)

        after = 0
        for g, idxs in groups.items():
            named = [(cds[i][3].strip(), cds[i][1]) for i in idxs
                     if informative(cds[i][3])]
            cons = consensus([n for n, _ in named]) if named else ''
            sources = sorted({src for n, src in named if n == cons})
            agree = sum(1 for n, _ in named if n == cons)
            disagree = len(named) - agree
            # Domains are per-group, not per-CDS: a group is one gene seen across
            # members, so take the domains any member carries. Unlike the product this
            # needs no transfer — antiSMASH scanned every member.
            gdoms = collections.Counter()
            for i in idxs:
                for d in cds[i][5]:
                    gdoms[d] += 1
            top_doms = [d for d, _ in gdoms.most_common(6)]
            roles = {domain_functions.category(d) for d in top_doms} - {'other'}
            per_group.append(dict(
                family=fam_id, group=g,
                consensus_product=cons or '(unnamed)',
                domains=';'.join(domain_functions.name(d) for d in top_doms),
                domain_accessions=';'.join(top_doms),
                role=sorted(roles)[0] if roles else 'other',
                group_size=len(idxs), n_annotated=len(named),
                n_sources=len(sources), n_agree=agree, n_disagree=disagree,
                prevalence=round(len(idxs) / len(members), 3)))
            for i in idxs:
                label, genome, tag, prod, _, _ = cds[i]
                orig_ok = informative(prod)
                transferred = prod.strip() if orig_ok else cons
                if transferred:
                    after += 1
                per_cds.append(dict(
                    family=fam_id, region=label, genome=genome, locus_tag=tag,
                    original_product=prod.strip(),
                    transferred_product=transferred,
                    origin='observed' if orig_ok else ('transferred' if cons else 'none'),
                    source_genomes=';'.join(sources[:5]) if not orig_ok and cons else '',
                    n_sources=len(sources) if not orig_ok and cons else '',
                    n_agree=agree if not orig_ok and cons else '',
                    n_disagree=disagree if not orig_ok and cons else '',
                    group_size=len(idxs)))

        n = len(cds)
        summary[str(fam_id)] = dict(
            members=len(members), cds=n, groups=len(groups),
            annotated_before=before, annotated_after=after,
            pct_before=round(100 * before / n, 1), pct_after=round(100 * after / n, 1))
        print(f'  GCF-{fam_id}: {len(members)} members, {n} CDS, {len(groups)} groups, '
              f'{100*before/n:.1f}% -> {100*after/n:.1f}% annotated')

    for name, rows, fields in (
        ('gcf_annotation_transfer.tsv', per_cds,
         ['family', 'region', 'genome', 'locus_tag', 'original_product',
          'transferred_product', 'origin', 'source_genomes', 'n_sources',
          'n_agree', 'n_disagree', 'group_size']),
        ('gcf_consensus_clusters.tsv', per_group,
         ['family', 'group', 'consensus_product', 'role', 'domains',
          'domain_accessions', 'prevalence', 'group_size', 'n_annotated',
          'n_sources', 'n_agree', 'n_disagree'])):
        with (args.outdir / name).open('w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
            w.writeheader()
            w.writerows(rows)

    tot_cds = sum(v['cds'] for v in summary.values())
    tot_b = sum(v['annotated_before'] for v in summary.values())
    tot_a = sum(v['annotated_after'] for v in summary.values())
    overall = dict(
        cutoff=args.cutoff, min_identity=args.min_identity,
        min_coverage=args.min_coverage, families=len(summary), cds=tot_cds,
        pct_before=round(100 * tot_b / tot_cds, 1) if tot_cds else 0.0,
        pct_after=round(100 * tot_a / tot_cds, 1) if tot_cds else 0.0,
        per_family=summary)
    (args.outdir / 'gcf_annotation_transfer.json').write_text(
        json.dumps(overall, indent=2))

    print(f'\n{tot_cds} CDS across {len(summary)} families: '
          f'{overall["pct_before"]}% -> {overall["pct_after"]}% annotated')

    for p in work.glob('fam*'):
        p.unlink()
    try:
        work.rmdir()
    except OSError:
        pass


if __name__ == '__main__':
    main()
