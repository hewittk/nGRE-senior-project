## Process for retrieval of gene sequences

Chromosome locations of the genes, adjusted for the proper amount of base pairs upstream and downstream, are needed to get the proper gene sequences. To get information about the genes' transcription start sequences, two files containing information about the upregulated and downregulated genes, respectively, were downloaded from the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables). A small subset of genes were not included in these files due to the names not being recognized by the UCSC database. The locations of these genes, 

`[regulation]_genes.txt`: genes in DEG file that were upregulated (fold change > 1) or downregulated (fold_change < 1)

 `ucsc_[regulation]_output`: information about list of genes from UCSC table browser including transciption start and end points
 `ucsc_[regulation]_unknown`: genes that were pronounced unknown to the UCSC table browser that aren't currently included in the `ucsc_[regulation]_output` datasets

 `[regulation]_gene_ENMUST_ids`: list of gene names in original DEG files translated to their ENMUST ids
