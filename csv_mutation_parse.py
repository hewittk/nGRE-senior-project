"""Make subset csv of all potential nGREs file that only contains rows of nGREs that have up to a maximum number of mutations."""

import csv
import time
import pandas as pd

def fewer_mutations_parse(treatment, regulation):
    """Parse nGRE sequence matches csv file and return all that have up to maximum number of mutations."""

    all_regions_df = pd.read_csv("nGRE_parse_output/" + treatment + "/" + regulation + "_gene_output.csv")
    promoter_region_df = pd.read_csv("nGRE_parse_output/" + treatment + "/" + regulation + "_promoter_gene_output.csv")

    # numbers of genes with zero or one mutations
    zero_mutation_subsequences = 0
    one_mutation_subsequences = 0

    with open('nGRE_parse_output/promoter_fewer_mutations_parse_statistics.txt', 'a') as output_file:
        print("In the " + treatment + " " + regulation + " dataset (length " + str(len(all_regions_df)) + " genes):", file = output_file)
        zero_mutations = promoter_region_df.loc[promoter_region_df['mutations'] == 0]
        print("Number of potential nGREs with zero mutations found in promoter region: ", len(zero_mutations), file = output_file)
        percent_zero_mutations = round(((len(zero_mutations)/len(all_regions_df)) * 100),2)
        print("Percent of genes containing nGREs with zero mutations in promoter region: ", str(percent_zero_mutations), "%", file = output_file)

        one_mutation = promoter_region_df.loc[promoter_region_df['mutations'] == 1]
        print("Number of potential nGREs with one mutation found in promoter region: ", len(one_mutation), file = output_file)
        percent_one_mutations = round(((len(one_mutation)/len(all_regions_df)) * 100),2)
        print("Percent of genes containing nGREs with one mutation in promoter region: ", str(percent_one_mutations), "%", file = output_file)
        print(" ", file = output_file)


def main():

	start = time.time()
	fewer_mutations_parse("ADX_F_DexvsADX_F_Veh", "downregulated")
	end = time.time()

	runtime = end - start
	print("Total runtime: {} seconds".format(runtime))


if __name__ == "__main__":
    main()
