# Deriving lists of upregulated/downregulated genes

### Files derived from table browser

`upregulated_output_with_gene_names.tsv`: information about list of upregulated genes from given treatment from UCSC table browser
`downregulated_output_with_gene_names.tsv`: information about list of downregulated genes of given treatment from UCSC table browser


# Process for retrieval of gene sequences

Chromosome locations of the genes, adjusted for the proper amount of base pairs upstream and downstream, are needed to get the proper gene sequences. To get information about the genes' transcription start sequences, two files containing information about the upregulated and downregulated genes, respectively, were downloaded from the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables).

A small subset of genes were not included in these files due to the names not being recognized by the UCSC database. The list of these genes were documented in each treatment's `ucsc_<regulation>_unknown.txt` files. Many of the genes that weren't recognized are currently not in the UCSC database because they haven't yet been officially annotated but have sequences that are known to be a relative gene site, and some other genes weren't annotated under their name given in the original dataset but may have other names that they are annotated under. A few other names that consisted of a single or low double digit number followed by a dash and a month abbreviation, such as `7-Mar` and `8-Sep`, were errors in the transfer of the original spreadsheet that can be omitted.

### Files derived from original datasets

`upregulated.txt`: list of genes in original dataset that were upregulated (fold change > 1)
`downregulated.txt`: list of genes in original dataset that were downregulated (fold change < 1)

### Columns in `<regulation>_output_with_gene_names.tsv`

 `mm10.knownGene.name`: Ensembl Stable ID of genetic transcript

 `mm10.knownGene.chrom`: Chromosome that genetic transcript exists on

 `mm10.knownGene.strand`: DNA strand (+ or -) that given gene transcript is found on

 `mm10.knownGene.txStart`: Start coordinate of gene's transcription on its chromosome

 `mm10.knownGene.txEnd`: End coordinate of gene's transcription on its chromosome

 `mm10.knownGene.cdsStart`: Start coordinate of gene's transcription on its chromosome

 `mm10.knownGene.cdsEnd`: End coordinate of gene's transcription on its chromosome

 `mm10.knownGene.proteinID`: "UniProt display ID, UniProt accession, or RefSeq protein ID" (according to [Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables))

 `mm10.knownGene.alignID`: "Unique identifier (GENCODE transcript ID for GENCODE Basic)" (according to [Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables))

 `mm10.kgXref.geneSymbol`: Common gene name that is same as name given in originral `upregulated.txt`/`downregulated.txt` files



 `ucsc_[regulation]_output`: information about list of genes from UCSC table browser including transciption start and end points

 `ucsc_[regulation]_unknown`: genes that were pronounced unknown to the UCSC table browser that aren't currently included in the `ucsc_[regulation]_output` datasets
