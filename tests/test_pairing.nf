/*
 * Pairing genomes with their recovered-ORF GFF3s.
 *
 * antismash_analysis.nf joins <genome>.gbff to <genome>.gff3 on the basename. The
 * obvious key, Nextflow's `simpleName`, strips EVERY extension -- so a genome named
 * `..._GCA_963520565.1` keys as `..._GCA_963520565`, fails to join, and is dropped
 * from antiSMASH without a word. 652 of 2,771 Erwiniaceae genome names contain a dot.
 *
 * This asserts the join keeps dotted names. It is a regression test, not a
 * hypothetical: the same class of bug already cost this pipeline a silent 25-of-333
 * mismatch in the novelty score.
 */
workflow {
    names = Channel.of('Pantoea_GCA_963520565.1', 'Simple_name', 'A.b.c', 'Buchnera_B.tra')
    genomes = names.map { n -> tuple(n, file("/nonexistent/${n}.gbff")) }
    gffs    = names.map { n -> file("/nonexistent/${n}.gff3") }

    paired = genomes
        .join(gffs.map { f -> tuple(f.name.replaceAll(/\.gff3$/, ''), f) })
    // Delimited, not newline-terminated: Nextflow view() output can arrive
    // on a single line, and a $-anchored grep then matches only the last entry.
    paired.view { k, g, f -> "PAIRED<${k}>" }
}
