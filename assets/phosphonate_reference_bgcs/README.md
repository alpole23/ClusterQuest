# Phosphonate reference BGCs

Characterised phosphonate clusters that **MIBiG does not carry**, supplied to BiG-SCAPE
via `--reference-dir` so every run measures a real cluster-level distance to them.

## Why this exists

antiSMASH's KnownClusterBlast cannot see this chemistry. Measured on the Erwiniaceae run
(334 regions):

| | |
|---|---:|
| regions with any KCB ranking | 248 |
| best-hit similarity — min / median / max | 6 / 15 / **15** |
| above the pipeline's floor (`sim > 15`) | **0** |
| hits to any of MIBiG's 13 pepM-bearing clusters | **1** (dehydrophos, similarity 11) |

The most common best hit was *luminmycin/glidobactin* on 236 regions. KCB scores shared
genes plus synteny, and an Erwiniaceae phosphonate cluster shares only its 2-3 core genes
with an actinomycete one — below the noise floor of unrelated hits.

BiG-SCAPE measures the same thing far better: domain content plus adjacency, all-pairs,
which is already the basis of every GCF assignment in this pipeline. Putting known
clusters into that matrix gives a distance that means something.

## What MIBiG already has

Searching MIBiG 4.0's proteins with the 7 curated pepM references finds **13** clusters
carrying a pepM homologue — five more than a keyword search finds:

    BGC0000383  luminmycin/glidobactin (Photorhabdus)   BGC0000937  fosfazinomycin A/B
    BGC0000806  phosphonoglycans (Glycomyces)           BGC0000938  fosfomycin
    BGC0000807  phosphonoglycans (Stackebrandtia)       BGC0001411  polysaccharide B (B. fragilis)
    BGC0000897  dehydrophos                             BGC0001739  phosphonoacetic acid
    BGC0000904  FR-900098                               BGC0001859  fosfomycin
    BGC0000926  rhizocticin A                           BGC0002036  dehydrofosmidomycin
    BGC0002670  fosfonochlorin

Set `--bigscape_mibig_version 4.0` to bring those in. Cost on a 334-BGC run: 2,766 BGCs
total, ~2.1 GB peak RAM against a 32 GB allocation, ~62 CPU-min.

**Note it changes clustering.** MIBiG entries join the same distance matrix, so family
composition can shift. That is why the default is off rather than on.

## What MIBiG lacks — this directory

| file | product | organism | source |
|---|---|---|---|
| `pantaphos_LMG5342.region001.gbk` | pantaphos | *Pantoea ananatis* LMG 5342 | curated HiVir cluster, Sanger-corrected (this work); see below |
| `argolaphos.region001.gbk` | argolaphos | *Streptomyces monomycini* NRRL B-24309 | MZ612424.1 |
| `bialaphos.region001.gbk` | bialaphos | *Streptomyces hygroscopicus* ATCC 21705 | KP026916.1 |
| `phosphinothricin_PTT.region001.gbk` | phosphinothricin tripeptide | *Streptomyces viridochromogenes* Tü494 | X65195.2 |
| `phosphonothrixin.region001.gbk` | phosphonothrixin | *Saccharothrix* sp. ST-888 | AB863705.1 |

All four of the fetched clusters carry a pepM (43-70% identity to the curated
references), so they are unambiguously phosphonate clusters — **but antiSMASH detects no
region in any of them**, at every strictness setting including `loose`. Its phosphonate
rule needs pepM *plus* a partner domain from a fixed list, and these deposits do not
supply one. `scripts/genome/make_reference_bgc.py` declares the region instead, which is
a statement of fact for a curated deposit: the whole record is the cluster.

Pantaphos matters most of the three gaps: it is a characterised bioactive small molecule
from an organism in the Erwiniaceae dataset, and GCF-1's 215 members *are* this cluster.
Its nearest MIBiG relative is BGC0000383 at 75.2% on the synthase alone — a different
cluster.

### Pantaphos: the curated HiVir cluster, Sanger-corrected

The file is no longer antiSMASH's region from the 2012 deposit. It is the curated HiVir
cluster — 12,531 bp, 12 CDS — generated from the corrected record by
`scripts/genome/make_reference_bgc.py`. **Re-run that script if the curated record
changes**; a raw export has no region feature and BiG-SCAPE ignores it silently.

Aligned against HE617160, the corrected sequence is the deposit's 801,910-814,435 plus
**five single-base insertions and nothing else**, each one lengthening a 7-base
homopolymer (the classic assembly error of that era). They are marked in the file as
`misc_difference`, confirmed via Sanger sequencing.

Two of those insertions repair genes RefSeq calls pseudo: `RS26500` (MFS transporter)
and `RS26505` (hypothetical) now translate as intact ORFs, and the MFS contributes
PF07690 to the reference's domain content.

Two edits were made to the curated record before that script ran. The first: Geneious
exports leave the record-level `SOURCE`/`ORGANISM` empty, and BiG-SCAPE stores that
verbatim — the reference arrived as organism `.`, taxonomy `Unknown`, where the other
four carry proper names. Both were filled from the record's own `source` feature and the
deposit's lineage.

The second, `RS26480` (2-phosphonomethylmaleate dehydratase small subunit), needed its
**end moved from 5086 to 5109**, and this is not cosmetic: the insertion at 5066-5073 falls inside
the gene, so the deposit's end leaves a 529 bp CDS. That is not a multiple of 3, has no
fuzzy start or end to trim against, and **BiG-SCAPE therefore discards the CDS with only
a log warning** — the reference would silently lose a gene and its PF00694. Reading on
from the annotated TTG gives 183 aa; 195 of the 212 homologues across the family share
the resulting C-terminus, so the read-through is right. 191 of them are 177 aa, i.e.
from the ATG at 4576, which is the shorter start model if one is preferred.

**The old file was a duplicate.** It was byte-identical to this pipeline's own antiSMASH
output for LMG 5342, and BiG-SCAPE deduplicates input on the sha256 of file bytes — so in
every run that included LMG 5342 it was dropped ("Skipping duplicate", INFO level) and
contributed nothing. The corrected sequence is not a duplicate of anything.

**Delimitation is now the limiting factor, not sequence.** This reference is the cluster
proper; antiSMASH's query regions are a fixed 13,338 bp window (rule core ± neighbourhood)
that carries ~5 kb of upstream flank and stops before the cluster's last four genes — the
MFS transporter, the hypothetical, the FMN reductase and the second ATP-grasp are absent
from all 215 regions of the family. Measured consequences, on the Erwiniaceae run:

| comparison | distance |
|---|---:|
| reference vs LMG 5342's own antiSMASH region | **0.364** (jaccard 0.62) |
| reference vs the 23 contig-edge members (`auto` → glocal, trims to the shared part) | 0.000-0.199, all ≤ 0.30 |
| reference vs the 192 complete members (compared end to end) | median 0.364, **none** ≤ 0.30 |

That is why every reference here carries **`contig_edge=True`** — set by
`make_reference_bgc.py`, with a `note` qualifier in each file saying it is not a claim
about the DNA. Under `auto`, BiG-SCAPE compares a pair by its shared part whenever either
record is on a contig edge, which is the right question for a curated cluster against a
window: *is this cluster's content present here*. The edge parameter is stored per run, so
the flag recomputes nothing. Measured on a 20-BGC set: BGCs within 0.30 of pantaphos went
from 3 to 11, LMG 5342's own region from 0.364 to 0.000, while genuinely different
clusters stayed put (0.946 and 0.953 either way).

The region boundaries themselves are a separate fix: see
`antismash_phosphonate_neighbourhood` in `nextflow.config`, which widens the flank
antiSMASH keeps from 5 kb to 10 kb so the whole cluster lands inside the region.

Already covered by `--bigscape_mibig_version`: dehydrofosmidomycin (BGC0002036),
rhizocticin A (BGC0000926), fosfomycin, FR-900098, dehydrophos, the two phosphonoglycans
and B. fragilis polysaccharide B.

Still missing, and worth adding:

- **phosphonoalamide** — *Streptomyces* sp. NRRL B-2790 (references PnaA/PnaD are in
  `reference_sequences/`; nearest MIBiG entry BGC0000434 at 34.7% is not the same cluster)
- **valinophos** — *S. durhamensis* NRRL B-3309 (VlpA/VlpB; nearest BGC0002039 at 34.2%)
- **fosmidomycin** proper — only *dehydro*fosmidomycin is deposited as a cluster
- **plumbemycin** — no nuccore record found; it is a *Bacillus* rhizocticin-family
  phosphonate, and rhizocticin A (BGC0000926) is the closest available stand-in

## Adding a cluster

Files must be **antiSMASH-processed** GenBank — BiG-SCAPE's `--reference-dir` requires
region records, not raw deposits. So:

1. Fetch the source GenBank (NCBI nuccore).
2. Run `scripts/genome/make_reference_bgc.py --in <deposit> --out <here>/<name>.region001.gbk`.
   If antiSMASH *does* call a region on the deposit, prefer its output instead.
3. Add a row to the table above.

BiG-SCAPE's AS5 reader walks a four-level hierarchy — region → cand_cluster →
protocluster → proto_core — and rejects the file at each missing level in turn, with a
different error each time. It also picks its parser from
`structured_comment['antiSMASH-Data']['Version']` and silently falls back to the
antiSMASH-4 reader (which wants `cluster`, not `region`) when that key is absent. The
script writes all of it; do not hand-edit a GenBank and expect it to load.

Filenames must be unique across this directory — BiG-SCAPE warns about duplicates.
