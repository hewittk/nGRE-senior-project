import requests
import os
import json
import time
import pandas as pd

genes_df = pd.read_csv("annotated_gene_datasets/ADX_F_DexvsADX_F_Veh/downregulated_genes/downregulated_output_with_gene_names.tsv", sep = "\t")

def geneCoordinates(column):
    chromosome = genes_df["mm10.knownGene.chrom"][column]
    transcriptionStart = int(genes_df["mm10.knownGene.txStart"][column]) - 30000
    transcriptionEnd = int(genes_df["mm10.knownGene.txEnd"][column]) + 5000
    geneName = genes_df["mm10.kgXref.geneSymbol"][column] + "_" + genes_df["#mm10.knownGene.name"][column]

    return chromosome, transcriptionStart, transcriptionEnd, geneName

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
    start = time.time()
    for i in range(len(genes_df)):
        chromosome, startCoordinate, endCoordinate, geneName = geneCoordinates(i)
        gene_json = request(chromosome, startCoordinate, endCoordinate)
        write_to_file(geneName, gene_json)
    end = time.time()
    print("Total runtime: ", end - start)

if __name__ == "__main__":
    main()
