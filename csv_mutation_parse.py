"""Make subset csv of all potential nGREs file that only contains rows of nGREs that have up to a maximum number of mutations."""

import csv
import time
import pandas as pd

def fewer_mutations_parse(treatment, regulation):
    """Parse nGRE sequence matches csv file and return all that have up to maximum number of mutations."""

    df = pd.read_csv("nGRE_parse_output/" + treatment + "/" + regulation + "_gene_output.csv")

    # numbers of genes with zero or one mutations
    zero_mutation_subsequences = 0
    one_mutation_subsequences = 0

    with open('nGRE_parse_output/fewer_mutations_parse_statistics.txt', 'a') as output_file:
        print("In the " + treatment + " " + regulation + " dataset (length " + str(len(df)) + " genes):", file = output_file)
        zero_mutations = df.loc[df['mutations'] == 0]
        print("Number of potential nGREs with zero mutations found: ", len(zero_mutations), file = output_file)
        percent_zero_mutations = round(((len(zero_mutations)/len(df)) * 100),2)
        print("Percent of genes containing nGREs with zero mutations: ", str(percent_zero_mutations), "%", file = output_file)

        one_mutation = df.loc[df['mutations'] == 1]
        print("Number of potential nGREs with one mutation found: ", len(one_mutation), file = output_file)
        percent_one_mutations = round(((len(one_mutation)/len(df)) * 100),2)
        print("Percent of genes containing nGREs with one mutation: ", str(percent_one_mutations), "%", file = output_file)
        print(" ", file = output_file)


def main():

	start = time.time()
	fewer_mutations_parse("ADX_F_DexvsADX_F_Veh", "downregulated")
	end = time.time()

	runtime = end - start
	print("Total runtime: {} seconds".format(runtime))


if __name__ == "__main__":
    main()
