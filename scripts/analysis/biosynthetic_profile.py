#!/usr/bin/env python3
"""Per-family biosynthetic domain profile, and which family splits are chemistry.

BiG-SCAPE compares whole regions, and a region is a rule core plus a symmetric
flank. At the 10 kb neighbourhood this pipeline uses, that flank reaches ~9 kb
past the core in both directions and reliably catches chromosomal neighbours.
Two groups of genomes with different neighbours therefore land in different
families even when their biosynthesis is identical.

That is not hypothetical, and the *Erwiniaceae* run has a measured example:
GCF 7 == GCF 21 and GCF 8 == GCF 19 are identical on filtered content, each
pairing a real family with a singleton. A third case is stronger still. Within
the 29-member half of the pantaphos family, 11 clusters (all *P. agglomerans*)
stay separate from the other 18 even after their regions are trimmed to the same
size and their extra genes removed -- neither group carries the ATP-grasp, so
what divides them is the ordinary flanking content of one species versus
another.

Note the 186/29 split ITSELF is not such a case, though an earlier version of
this docstring used it as the headline example. Ablating PF13535 and
re-clustering merges 18 of the 29 into the large family, so that split is
substantially biosynthetic. See `docs/comparisons/pantaphos_family_split/`.

So this reports each family's domain content **filtered to the categories that
are about making a molecule**: core, tailoring, lipid and transport.
`primary` and `mobile` are dropped by construction, which is what
`utils/domain_functions` exists for. It is deliberately NOT the antiSMASH
`proto_core`: that is 3.3 kb here -- pepM, a monooxygenase and a homoaconitate
synthase -- and would discard LeuC/LeuD, the methyltransferase, SanS and the
GNAT acetyltransferase, which is most of the tailoring chemistry anyone would
want to compare.

The actionable output is the last column. Two families whose filtered profiles
are **identical** were separated by something other than their chemistry, and
that is worth knowing before a family count is read as biological resolution.

Two families that are merely *close* are a different claim, and the distinction
is the point rather than a nicety. On this run the pantaphos pair sits at Jaccard
0.909 and differs by exactly one filtered domain — PF13535 ATP-grasp_4, an
amide-bond ligase, in 96.8% of the 186-member family and 3.4% of the 29-member
one. That is the single most interesting difference in the dataset, and an
earlier version of this script printed "chemically indistinguishable — split is
not biosynthetic" over it, because one threshold covered both cases. A near miss
is now reported with the differing domains named, and never as indistinguishable.

**`other` is not evidence of absence.** The domain map is curated and covers
~93% of observed hits; everything else returns `other` and is excluded here by
ignorance rather than by classification. A family whose chemistry lives in
unmapped domains will look empty. Read a small profile as "nothing recognised",
never as "nothing there".

    python biosynthetic_profile.py --db run.db --out biosynthetic_profile.tsv
"""
import argparse
import collections
import csv
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.domain_functions import category, name  # noqa: E402

# Categories that describe making a molecule. `primary` and `mobile` are the
# two the region boundary sweeps in, and `other` means unmapped, not absent.
BIOSYNTHETIC = ('core', 'tailoring', 'lipid', 'transport')

# A domain is part of a family's profile when at least this share of its members
# carry it. Below it the domain is accessory -- often a neighbour one member's
# boundary happened to catch -- and including it would reintroduce the noise
# this script exists to remove.
PREVALENCE = 0.5

# Filtered profiles this similar are reported as a possible context-driven split.
# Jaccard over domain sets. Only 1.0 means chemically indistinguishable; between
# SIMILAR and 1.0 the families differ by a domain or two and the difference is
# named rather than dismissed. Collapsing those two cases into one verdict is a
# mistake this script made and got caught on -- see `verdict_for`.
SIMILAR = 0.9


def load(db, cutoff):
    """{family_id: {record_id: {domain accessions}}} for region records."""
    con = sqlite3.connect(db)
    con.row_factory = sqlite3.Row
    fam = {}
    for r in con.execute(
            """SELECT br.id rid, br.gbk_id gid, f.id fid FROM bgc_record br
               JOIN bgc_record_family brf ON brf.record_id = br.id
               JOIN family f ON f.id = brf.family_id
               WHERE br.record_type = 'region' AND f.cutoff = ?""", (cutoff,)):
        fam[r['rid']] = (r['fid'], r['gid'])
    if not fam:
        con.close()
        sys.exit(f'no region records at cutoff {cutoff}; check --cutoff')

    by_gbk = collections.defaultdict(set)
    for r in con.execute("""SELECT cds.gbk_id gid, hsp.accession acc
                            FROM cds JOIN hsp ON hsp.cds_id = cds.id"""):
        by_gbk[r['gid']].add(r['acc'].split('.')[0])
    con.close()

    out = collections.defaultdict(dict)
    for rid, (fid, gid) in fam.items():
        out[fid][rid] = by_gbk.get(gid, set())
    return out


def profile(members):
    """{domain: prevalence} over the biosynthetic categories only."""
    n = len(members)
    counts = collections.Counter()
    for doms in members.values():
        for d in doms:
            if category(d) in BIOSYNTHETIC:
                counts[d] += 1
    return {d: c / n for d, c in counts.items() if c / n >= PREVALENCE}


def jaccard(a, b):
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b) if (a | b) else 1.0


def verdict_for(sim, other, diff_txt, n_diff):
    """What a nearest-neighbour similarity actually licenses saying.

    Only an EMPTY symmetric difference licenses "indistinguishable". An earlier
    version applied that sentence to everything at Jaccard >= 0.9, and on the
    Erwiniaceae run that put it on the pantaphos pair, which differs by exactly
    one filtered domain: PF13535 ATP-grasp_4, in 96.8% of the 186-member family
    and 3.4% of the 29-member one. An ATP-grasp is an amide-bond ligase, so the
    one domain the threshold waved through is the one that would separate a
    monopeptide product from a dipeptide -- the tool argued against the most
    interesting hypothesis in the run on a rounding margin (10/11 = 0.909).
    """
    if not other or sim < SIMILAR:
        return ''
    if not n_diff:
        return (f'chemically indistinguishable from GCF {other} '
                '— split is not biosynthetic')
    return (f'differs from GCF {other} by {n_diff} biosynthetic domain'
            f'{"s" if n_diff > 1 else ""} only ({diff_txt}) — near-identical, '
            'but NOT indistinguishable; check the difference before calling '
            'this split contextual')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--db', required=True, help='BiG-SCAPE SQLite database')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--cutoff', type=float, default=0.3)
    ap.add_argument('--prevalence', type=float, default=PREVALENCE)
    args = ap.parse_args()

    fams = load(args.db, args.cutoff)
    globals()['PREVALENCE'] = args.prevalence

    profiles = {fid: profile(m) for fid, m in fams.items()}
    sets = {fid: set(p) for fid, p in profiles.items()}

    # Nearest other family by filtered content. A family that is chemically
    # indistinguishable from another was split by something else.
    nearest = {}
    for fid in sets:
        best = max(((jaccard(sets[fid], sets[o]), o) for o in sets if o != fid),
                   default=(0.0, ''))
        nearest[fid] = best

    rows = []
    for fid in sorted(fams, key=lambda f: -len(fams[f])):
        p = profiles[fid]
        cats = collections.Counter(category(d) for d in p)
        sim, other = nearest[fid]
        diff = (sets[fid] ^ sets[other]) if other != '' else set()
        diff_txt = ';'.join(f'{d}({name(d) or "?"})' for d in sorted(diff))
        rows.append({
            'gcf': fid,
            'members': len(fams[fid]),
            'biosynthetic_domains': len(p),
            'core': cats.get('core', 0),
            'tailoring': cats.get('tailoring', 0),
            'lipid': cats.get('lipid', 0),
            'transport': cats.get('transport', 0),
            'domains': ';'.join(f'{d}({name(d) or "?"})'
                                for d in sorted(p, key=lambda x: (-p[x], x))),
            'nearest_gcf': other,
            'nearest_similarity': f'{sim:.3f}',
            'differing_domains': diff_txt if sim >= SIMILAR else '',
            'verdict': verdict_for(sim, other, diff_txt, len(diff)),
        })

    with args.out.open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    same = [r for r in rows if r['verdict'] and not r['differing_domains']]
    near = [r for r in rows if r['verdict'] and r['differing_domains']]
    print(f'{len(rows)} families at cutoff {args.cutoff} '
          f'(Jaccard on {"/".join(BIOSYNTHETIC)} domains)')
    print(f'{len(same)} chemically indistinguishable from another family:')
    for r in same:
        print(f"  GCF {r['gcf']:>3} ({r['members']:>3} members) == "
              f"GCF {r['nearest_gcf']:>3}  J={r['nearest_similarity']}")
    print(f'{len(near)} near-identical but distinguishable — read the domain, '
          f'not the number:')
    for r in near:
        print(f"  GCF {r['gcf']:>3} ({r['members']:>3} members) ~ "
              f"GCF {r['nearest_gcf']:>3}  J={r['nearest_similarity']}  "
              f"differs by {r['differing_domains']}")
    empty = [r for r in rows if r['biosynthetic_domains'] == 0]
    if empty:
        print(f'{len(empty)} families have no recognised biosynthetic domain — '
              f'unmapped, not empty; see the docstring')
    return 0


if __name__ == '__main__':
    sys.exit(main())
