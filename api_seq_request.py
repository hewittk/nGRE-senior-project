import requests
import pandas as pd

genes_df = pd.read_csv("annotated_gene_datasets/ADX_F_DexvsADX_F_Veh/downregulated_genes/downregulated_output_with_gene_names.tsv", sep = "\t")

def geneCoordinates(dataframe):
    print("In geneCoordinates")
    chromosome = dataframe["mm10.knownGene.chrom"][1]
    transcriptionStart = dataframe["mm10.knownGene.txStart"][1]
    transcriptionEnd = dataframe["mm10.knownGene.txEnd"][1]

    return chromosome, transcriptionStart, transcriptionEnd

def method():
    # get gene name and chromosome, txStart, txEnd for each gene


    # request sequence for gene
    response = requests.get("http://api.genome.ucsc.edu/getData/sequence?genome=mm10;chrom=chr1;start=20542952;end=20618064")
    print(response.status_code)
    print(response.json())

    # mine json for response

    # save to file named after gene


def main():
    print("In main")
    print(geneCoordinates(genes_df))


if __name__ == "__main__":
    main()
