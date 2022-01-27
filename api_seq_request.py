import requests
import os
import json
import time
import pandas as pd

genes_df = pd.read_csv("annotated_gene_datasets/OVX_ADX_F_DexvsOVX_ADX_F_Veh/upregulated_genes/upregulated_output_with_gene_names.tsv", sep = "\t")

def geneCoordinates(column):
    """Retrieve relevant gene name/location information necessary to make and store information from API request."""

    # retrieve gene chromosome and relevant chromosome coordinates to start and end search for gene's nGREs
    chromosome = genes_df["mm10.knownGene.chrom"][column]
    transcriptionStart = int(genes_df["mm10.knownGene.txStart"][column]) - 30000
    transcriptionEnd = int(genes_df["mm10.knownGene.txEnd"][column]) + 5000

    # make name to name gene's file by combining the common gene symbol and the ENMUST ID
    geneName = genes_df["mm10.kgXref.geneSymbol"][column] + "_" + genes_df["#mm10.knownGene.name"][column]
    if geneName[len(geneName) - 1] == 0:
        geneName = geneName[:len(geneName) - 1]
    print("geneName: ", geneName)

    return chromosome, transcriptionStart, transcriptionEnd, geneName

def request(chromosome, startCoordinate, endCoordinate):
    """Send API request for relevant sequence and store response."""

    response = requests.get("http://api.genome.ucsc.edu/getData/sequence?genome=mm10;chrom=" + str(chromosome) + ";start=" + str(startCoordinate) + ";end=" + str(endCoordinate))
    print("Status Code: ", response.status_code)
    json_response = response.json()

    return json_response['dna']

def write_to_file(gene_name, string):
    """Write gene sequence to file."""

    # make file for gene in gene sequences directory
    print("Writing to file: ", gene_name)
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
        if not(os.path.isdir("gene_sequences/" + geneName)):
            print("Adding new gene sequence file")
            gene_json = request(chromosome, startCoordinate, endCoordinate)
            write_to_file(geneName, gene_json)
        else:
            print("Gene sequence file already exists")
        print()
    end = time.time()
    print("Total runtime: ", end - start)

if __name__ == "__main__":
    main()
