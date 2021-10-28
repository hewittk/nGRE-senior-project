import requests
import os
import json
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

    json_sequence = response.json()
    print(json_sequence)
    # print(json_sequence['dna'])

    return json_sequence['dna']
    # print(response.status_code)

def write_to_file(gene_name, string):
    # save to file named after gene
    os.mkdir("gene_sequences/" + gene_name)
    file_name = gene_name + ".txt"

    f = open("gene_sequences/" + gene_name + "/" + file_name, "x")
    f.close()

    with open("gene_sequences/" + gene_name + "/" + file_name, "w") as file:
        file.write(string)
        file.close()


def main():
    print("In main")
    """
    for i in range(len(genes_df)):
        # chromosome, startCoordinate, endCoordinate = geneCoordinates(i)
        # json_response = request(chromosome, startCoordinate, endCoordinate)
    """
    chromosome, startCoordinate, endCoordinate = geneCoordinates(0)
    json_response = request(chromosome, startCoordinate, endCoordinate)
    write_to_file("Genee", json_response)

if __name__ == "__main__":
    main()
