"""Make subset csv of all potential nGREs file that only contains rows of nGREs that have up to a maximum number of mutations."""

import csv
import time
import pandas as pd
import regex

treatment = ""
regulation = ""

def fewer_mutations_parse():
    """Parse nGRE sequence matches csv file and return all that have zero or one mutations."""

    global treatment
    global regulation

    # dataset containing information about each gene/ENMUST ID
    genes_df = pd.read_csv("annotated_gene_datasets/" + treatment + "/" + regulation + "_genes/" + regulation + "_output_with_gene_names_known_sequences.tsv", sep = "\t", index_col=False)
    genes_df = genes_df.reset_index(drop=True)
    gene_name_list = genes_df["mm10.kgXref.geneSymbol"].to_list()
    gene_name_set = set(gene_name_list)

    # dataset containing information about nGREs found
    promoter_region_df = pd.read_csv("nGRE_parse_output/" + treatment + "/" + regulation + "_promoter_gene_output.csv", index_col=False) # nGREs in promoter region of gene
    all_region_df = pd.read_csv("nGRE_parse_output/" + treatment + "/" + regulation + "_gene_output.csv", index_col=False) # nGREs in all regions of gene

    # numbers of genes with zero or one mutations
    zero_mutation_subsequences = 0
    one_mutation_subsequences = 0

    with open('nGRE_parse_output/promoter_fewer_mutations_parse_statistics.txt', 'a') as output_file:
        print("In the " + treatment + " " + regulation + " dataset (length " + str(len(gene_name_set)) + " genes):", file = output_file)

        # amount of genes with zero mutations
        zero_mutations = promoter_region_df.loc[promoter_region_df['mutations'] == 0]
        number_zero_mutations, percent_zero_mutations = find_nGRE_amounts(zero_mutations, genes_df, gene_name_set)
        print("Number of potential nGREs with zero mutations found in all regions: ", number_zero_mutations, file = output_file)
        print(number_zero_mutations)
        print("Percent of genes containing nGREs with zero mutations in all regions: ", str(percent_zero_mutations), "%", file = output_file)

        # amount of genes with one mutation
        one_mutations = promoter_region_df.loc[promoter_region_df['mutations'] == 1]
        number_one_mutations, percent_one_mutations = find_nGRE_amounts(one_mutations, genes_df, gene_name_set)
        print("Number of potential nGREs with one mutation found in all regions: ", number_one_mutations, file = output_file)
        print("Percent of genes containing nGREs with one mutation in all regions: ", str(percent_one_mutations), "%", file = output_file)
        print(" ", file = output_file)

def find_nGRE_amounts(nGRE_df, all_genes_df, gene_name_set):
    """Find total numbers of nGREs and proportion of genes in dataset containing nGREs."""

    global treatment
    global regulation

    # initialize dictionary of lists of nGRE locations for each gene
    locations_dict = {}
    nGRE_df = nGRE_df.reset_index()

    # initialize dictionary of genes containing the highest numbers of nGREs
    highest_nGRE_genes = {}


    # parse each gene for nGREs
    total_nGREs = 0 # nGRE counter
    for column in range(len(nGRE_df)):

        # derive gene name and enmust id
        gene_id = nGRE_df["gene_id"][column]
        gene_name = regex.search(r'.+_', gene_id)
        gene_name = str(gene_name.group(0)[:(len(gene_name.group(0)) - 1)])
        enmust_id = regex.search(r'_\w+\.?\d+', gene_id)
        enmust_id = str(enmust_id.group(0)[1:])

        # derive given nGRE's start site
        relative_nGRE_start = int(nGRE_df["start_site"][column])

        # derive enmust id's start/end sites on chromosome
        enmust_id_info = all_genes_df.loc[all_genes_df['#mm10.knownGene.name'] == enmust_id]
        txStart_data = enmust_id_info['mm10.knownGene.txStart']
        txStart = txStart_data.values[0]
        cdsStart_data = enmust_id_info['mm10.knownGene.cdsStart']
        cdsStart = cdsStart_data.values[0]
        cdsEnd_data = enmust_id_info['mm10.knownGene.cdsEnd']
        cdsEnd = cdsEnd_data.values[0]
        txEnd_data = enmust_id_info['mm10.knownGene.txEnd']
        txEnd = txEnd_data.values[0]

        # find non-relative, chromosomal start site of given nGRE
        chromosomal_nGRE_start = txStart + (relative_nGRE_start - 30000)

        # add each nGRE found to location_dict
        if gene_name in locations_dict:
            if chromosomal_nGRE_start not in locations_dict[gene_name]:
                locations_dict[gene_name].append(chromosomal_nGRE_start)
                total_nGREs += 1
        else:
            locations_dict[gene_name] = [chromosomal_nGRE_start]

        if (gene_id == "Lipg" or gene_id == "Pstpip2" or gene_id == "Banf1" or gene_id == "Yif1a"):
            print(gene_id)
            print(locations_dict)

        # check if gene is eligible to be added to top 100 genes containing the most nGREs
        if (len(highest_nGRE_genes) < 100): # automatically add first 100 genes
            highest_nGRE_genes[gene_name] = len(locations_dict[gene_name])
        else:
            for standing_gene in highest_nGRE_genes:
                if(len(locations_dict[gene_name]) > highest_nGRE_genes[standing_gene]):
                    # remove standing gene in highest_nGRE_genes that has less
                    highest_nGRE_genes.pop(standing_gene)
                    # add current gene containing more nGREs
                    highest_nGRE_genes[gene_name] = len(locations_dict[gene_name])
                    break

    # write genes found with the most nGREs to a file
    with open('nGRE_parse_output/' + treatment + '/' + regulation + '_promoter_top100_nGRE_genes.txt', 'w') as output_file:
        for gene in highest_nGRE_genes:
            print(gene + ": " + str(highest_nGRE_genes[gene]), file = output_file)

    with open('nGRE_parse_output/' + treatment + '/' + regulation + '_list_promoter_top100_nGRE_genes.txt', 'w') as output_file:
        for gene in highest_nGRE_genes:
            print(gene, file = output_file)

    # divide total number of nGREs found by amount of genes in dataset
    nGRE_genes = 0
    for gene in gene_name_set:
        if gene in locations_dict:
            nGRE_genes += 1
            # print("Gene: " + gene + " nGRE genes: " + str(nGRE_genes))
    percent_nGRE_genes = round((nGRE_genes/len(gene_name_set)*100),2)

    return nGRE_genes, percent_nGRE_genes


def main():

    global treatment
    global regulation

    treatment = "ADX_F_DexvsADX_F_Veh"
    regulation = "downregulated"

    start = time.time()
    fewer_mutations_parse()
    end = time.time()

    runtime = end - start
    print("Total runtime: {} seconds".format(runtime))


if __name__ == "__main__":
    main()
