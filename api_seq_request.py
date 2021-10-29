import requests
import os
import json
import pandas as pd

genes_df = pd.read_csv("annotated_gene_datasets/ADX_F_DexvsADX_F_Veh/downregulated_genes/downregulated_output_with_gene_names.tsv", sep = "\t")

def geneCoordinates(column):
    chromosome = genes_df["mm10.knownGene.chrom"][column]
    transcriptionStart = int(genes_df["mm10.knownGene.txStart"][column]) - 30000
    transcriptionEnd = int(genes_df["mm10.knownGene.txEnd"][column]) + 5000

    return chromosome, transcriptionStart, transcriptionEnd

def request(chromosome, startCoordinate, endCoordinate):
    response = requests.get("http://api.genome.ucsc.edu/getData/sequence?genome=mm10;chrom=" + str(chromosome) + ";start=" + str(startCoordinate) + ";end=" + str(endCoordinate))

    json_response = response.json()

    return json_response['dna']

def write_to_file(gene_name, string):

    # make file for gene in gene sequences directory
    os.mkdir("gene_sequences/" + gene_name)
    file_name = gene_name + ".txt"
    f = open("gene_sequences/" + gene_name + "/" + file_name, "x")
    f.close()

    # write gene's sequence to file
    with open("gene_sequences/" + gene_name + "/" + file_name, "w") as file:
        file.write(string)
        file.close()

def main():
    """
    for i in range(len(genes_df)):
        # chromosome, startCoordinate, endCoordinate = geneCoordinates(i)
        # json_response = request(chromosome, startCoordinate, endCoordinate)
    """
    chromosome, startCoordinate, endCoordinate = geneCoordinates(0)
    gene_json = request(chromosome, startCoordinate, endCoordinate)
    write_to_file("Genee", gene_json)

if __name__ == "__main__":
    main()
