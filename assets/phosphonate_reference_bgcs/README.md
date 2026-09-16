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
| `pantaphos_LMG5342.region001.gbk` | pantaphos | *Pantoea ananatis* LMG 5342 | this pipeline's own antiSMASH output, HE617160 region 1 |

Pantaphos matters most of the three gaps: it is a characterised bioactive small molecule
from an organism in the Erwiniaceae dataset, and GCF-1's 215 members *are* this cluster.
Its nearest MIBiG relative is BGC0000383 at 75.2% on the synthase alone — a different
cluster.

Still missing, and worth adding:

- **phosphonoalamide** — *Streptomyces* sp. NRRL B-2790 (references PnaA/PnaD are in
  `reference_sequences/`; nearest MIBiG entry BGC0000434 at 34.7% is not the same cluster)
- **valinophos** — *S. durhamensis* NRRL B-3309 (VlpA/VlpB; nearest BGC0002039 at 34.2%)

## Adding a cluster

Files must be **antiSMASH-processed** GenBank — BiG-SCAPE's `--reference-dir` requires
region records, not raw deposits. So:

1. Fetch the source GenBank (NCBI nuccore).
2. Run antiSMASH on it, without `--hmmdetection-limit-to-rule-names` if the cluster is
   not a phosphonate by antiSMASH's rules.
3. Copy the resulting `*.regionNNN.gbk` here and add a row to the table above.

Filenames must be unique across this directory — BiG-SCAPE warns about duplicates.
