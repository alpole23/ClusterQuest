# How ClusterQuest works

A plain-language walkthrough, for readers who want to know what the pipeline does
without reading Nextflow or Python. Figure S1 (`docs/figures/figS1_architecture.png`)
is this document as one picture.

## What it is looking for

Most natural products are built from carbon, nitrogen and oxygen. A **phosphonate** has
something rarer: a direct bond between a carbon atom and a phosphorus atom. That bond is
chemically stubborn — it resists the enzymes that break ordinary phosphate esters — which
is why phosphonates make useful drugs and herbicides. Fosfomycin, still used for urinary
infections, is one.

The genes a bacterium needs to build one of these compounds usually sit side by side on
the chromosome, in what is called a **biosynthetic gene cluster**. Every phosphonate
cluster known begins the same way, with an enzyme called **PEP mutase** that performs the
single chemical step creating the C–P bond. That shared starting enzyme is the hook this
pipeline hangs everything on: a genome with no PEP mutase cannot contain what we are
looking for.

ClusterQuest takes the name of a bacterial group — a species, a genus, an entire order —
and answers three questions about it: which genomes carry a phosphonate cluster, which of
those clusters are variations on the same theme, and what each theme probably makes.

## The shape of the problem

Almost nothing matches. On the full GenBank set for the order *Enterobacterales*:

> **1,349,504 genomes** examined → **1,530** carry a PEP mutase → **1,476** clusters
> found → **112** families.

That is 99.89% discarded. Every design decision below follows from it. The expensive
analysis tools are good but slow, so the pipeline's real job is to make sure they only
ever see the 0.11% that could possibly be interesting — and to throw the rest away as
early and as cheaply as it can.

## The steps

### 1. Ask NCBI what exists, without downloading it

The public archive is asked for a *list* of genomes in the chosen group, plus a line of
metadata about each — organism name, strain, assembly quality. No sequence yet. A listing
entry is about 3.7 kB against roughly 8.7 MB for the genome itself, so this step is
essentially free even for a million genomes.

### 2. Decide what to call everything

Archive accessions like `GCA_963520565.1` are unreadable, and two deposits of the same
organism must not collide. Each genome is given a name built from its organism and
strain, with numbered suffixes where that would produce a duplicate.

### 3. Build the search set

A small searchable database is built from seven PEP mutase proteins taken from clusters
whose chemistry has been worked out in a laboratory. This is the reference the screen in
the next step compares against.

### 4. Fetch, screen and discard — all in one step

**This is the central design decision and about 91% of the computing time.**

Genomes are handled twenty-five at a time. For each batch the pipeline downloads the
sequences, renames them, searches each one for a PEP mutase, keeps the few that have one,
and **deletes the rest before the next batch starts**.

The ordering is what matters. An obvious design would download everything first and
filter afterwards — but for this order that means writing 1.3 million genomes to disk,
around 12 terabytes, to keep 1,530 of them. Screening inside the download means a
rejected genome is never written at all. Disk use stays at a few hundred megabytes
regardless of whether the group has a thousand genomes or a million.

The screen itself is deliberately crude and fast: roughly 0.9 seconds of computing per
genome, against 41 seconds for the full analysis. It is tuned to never miss — across
every validation set tried, including one deliberately chosen to share no reference
sequences with the search set, it found **every** genome that a full analysis later
confirmed. It does let a few through that turn out to have nothing, which costs a little
time and loses no data.

### 5. Look up where each survivor sits in the tree of life

The surviving genomes get their formal classification from NCBI's taxonomy — kingdom
down to species — so that later results can be read by clade.

### 6. Fill in missing gene annotations

Public genome deposits vary enormously in quality. Many list the organism's DNA but label
only some of its genes, and a few label none. This is a real problem rather than a
cosmetic one: the detection tool used in the next step assumes that if a deposit labels
*any* genes, it labels all of them, so it will not look for the missing ones.

The consequence was found in a cluster whose chemistry is known from laboratory work.
Its 2012 deposit labels 6 of its 15 genes; the other 9, including two that matter for
identifying the product, were invisible to every measurement built on top.

So before detection, the pipeline pools the proteins from every retained genome, aligns
that pool back against each genome's DNA, and recovers the genes that were never
labelled. These are handed to the detection step as genuine genes — not as a footnote —
which means a recovered gene can itself trigger detection and extend the boundaries of a
cluster.

### 7. Detect the clusters

**antiSMASH**, the standard tool for this task, examines each surviving genome and marks
the region around any phosphonate signature it finds. ClusterQuest widens the window
antiSMASH uses by default, from 5 kb to 10 kb either side, because the narrower setting
was measured to cut four genes off the end of a known cluster.

### 8. Count what was found

The clusters are counted per genome and tabulated with their coordinates, their genes and
the organism they came from. This is the raw result table everything downstream reads.

### 9. Group the clusters into families

**BiG-SCAPE** compares every cluster against every other one, scoring them on which genes
they share and in what order, and groups the similar ones into **gene cluster families**.
A family is a set of clusters that are recognisably variations on one biosynthetic theme,
usually meaning they make the same compound or close relatives of it.

This step is quadratic — doubling the clusters quadruples the work — which would be a
serious constraint if the funnel were not so narrow. Having only 1,476 clusters to
compare makes it under 0.1% of the run.

### 10. Measure the distance to known chemistry

Seventeen clusters whose products have been characterised are compared against every
family, so each family can be reported as *this is a known compound*, *this is a distant
relative of one*, or *nothing published resembles this*.

The comparison runs on a **copy** of the results database, and this is not fussiness: a
reference similar enough to a real cluster would otherwise join its family and be counted
as though a laboratory strain were part of the dataset.

### 11. Place the genomes on the bacterial tree

**GTDB-Tk** assigns each cluster-carrying genome a position in a standardised
phylogeny — more consistent than NCBI's names, which carry historical baggage. This is
what allows a statement like "this family is confined to one genus" to be checked rather
than assumed.

### 12. Work out what each family makes

Several analyses run over the families, each answering a different question:

- **Which enzyme acts next?** The step immediately after PEP mutase decides which branch
  of phosphonate chemistry the cluster is on. The pipeline identifies that enzyme from
  the cluster's gene content and reports the branch, along with how good the evidence is.
- **Which intermediate does it go through?** Two compounds, 2-AEP and 2-HEP, account for
  most phosphonate pathways; profile-based searching calls which one applies, or says
  neither.
- **What does a typical member look like?** Rather than show one arbitrary cluster, each
  family is summarised by which genes appear in what fraction of its members, so a reader
  can tell the defining genes from the incidental neighbours.
- **What are the unlabelled genes?** Where one well-annotated member of a family names a
  gene its relatives leave blank, that name is carried across — raising usable gene
  labels from 40% to 80% in one measured run. Every transferred label records how many
  genomes it rests on, because a name resting on one genome is one genome's opinion.

### 13. Write the report

Everything is assembled into a single self-contained HTML file: summary statistics, every
cluster found, the families with their consensus gene content and predicted chemistry,
the taxonomic distribution, and a page per genome. Alongside it the pipeline writes plain
tab-separated tables and standard-format trees, so the results can be taken into other
software rather than only read.

## Two things worth knowing about the result

**Cost runs opposite to intuition.** The expensive-sounding analysis is cheap here and
the plumbing is expensive. Downloading and filtering is 91% of the computing; the cluster
detection everybody thinks of as the hard part is 4%; the step that is mathematically
quadratic is 0.07%. At this rarity, moving and discarding data costs more than analysing
it — which is why the screen is the architecture rather than a speed-up bolted onto it.

**The pipeline reports evidence, not verdicts.** Predicted chemistry is inference from
gene content, and the set of phosphonate enzymes anyone has characterised in a laboratory
is small. Family groupings are reproducible and solid; what a family *makes* is a
hypothesis the output is designed to let you check, with every call carrying how well it
is evidenced and what it rests on.
