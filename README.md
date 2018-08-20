# motifAnalysis
identification of shared motifs on a subset of genes, with an enrichment test against all genes in a genome

Algorithm (pseudocode)

Input: * a full set of transcripts (e.g. genes, promoters, 5' UTRs, 3' UTRs) and a subset
the subset could be a set of genes forming a pathway, or a set of genes co-expressed
to ensure a meaningful significance analysis, use a subset with at least 10 transcripts
* an annotation set used for annotating the motifs (e.g. a list of annotated miRNAs or transcription factors)

Output: a set of motifs significantly enriched or depleted relative to the full set.
For each motif, its potential annotation and the list of genes sharing that motif are presented.

The analysis investigates all possible windows with length between Wmin and Wmax
windows with low sequence complexity (e.g. AAAAA, ATATAT) are excluded

the number of occurrences of a motif in the subset and the full set are used as input for a Fisher exact test

the significant motifs are annotated (on both positive and negative strand).
