import requests
import pandas as pd

genes_df = pd.read_csv("annotated_gene_datasets/ADX_F_DexvsADX_F_Veh/downregulated_genes/downregulated_output_with_gene_names.tsv", sep = "\t")

def geneCoordinates(column):
    print("In geneCoordinates")
    chromosome = genes_df["mm10.knownGene.chrom"][column]
    transcriptionStart = genes_df["mm10.knownGene.txStart"][column]
    transcriptionEnd = genes_df["mm10.knownGene.txEnd"][column]

    return chromosome, transcriptionStart, transcriptionEnd

def request(chromosome, startCoordinate, endCoordinate):
    response = requests.get("http://api.genome.ucsc.edu/getData/sequence?genome=mm10;chrom=" + str(chromosome) + ";start=" + str(startCoordinate) + ";end=" + str(endCoordinate))
    print(response.status_code)

def method():
    # get gene name and chromosome, txStart, txEnd for each gene


    # request sequence for gene
    response = requests.get("http://api.genome.ucsc.edu/getData/sequence?genome=mm10;chrom=chr1;start=20542952;end=20618064")
    # print(response.json())

    # mine json for response

    # save to file named after gene


def main():
    print("In main")
    for i in range(len(genes_df)):
        chromosome, startCoordinate, endCoordinate = geneCoordinates(i)
        request(chromosome, startCoordinate, endCoordinate)

if __name__ == "__main__":
    main()
