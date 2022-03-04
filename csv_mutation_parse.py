"""Make subset csv of all potential nGREs file that only contains rows of nGREs that have up to a maximum number of mutations."""

import csv
import time
import pandas as pd

def fewer_mutations_parse(treatment, regulation):
    """Parse potential nGRE csv file and return all that have up to maximum number of mutations."""

    df = pd.read_csv("nGRE_parse_output/" + treatment + "/" + regulation + "_gene_output.csv")

    # trackers for how many genes in subset have zero or one mutations
    zero_mutation_subsequences = 0
    one_mutation_subsequences = 0

    # dataset containing potential nGREs with fewer mutations
    fewer_mutations_table = pd.DataFrame(columns = ["gene_id", "chromosome", "nGRE_sequence", "start_site", "end_site", "mutations", "mismatch_mutations", "insertion_mutations", "deletion_mutations"])

    with open('nGRE_parse_output/fewer_mutations_parse_statistics.txt', 'a') as output_file:
        print("In the " + treatment + " " + regulation + " dataset (length " + str(len(df)) + " genes):", file = output_file)
        zero_mutations = df.loc[df['mutations'] == 0]
        fewer_mutations_table = fewer_mutations_table.append(zero_mutations, ignore_index = True)
        print("Number of potential nGREs with zero mutations found: ", len(zero_mutations), file = output_file)

        one_mutation = df.loc[df['mutations'] == 1]
        fewer_mutations_table = fewer_mutations_table.append(one_mutation, ignore_index = True)
        print("Number of potential nGREs with one mutation found: ", len(one_mutation), file = output_file)
        print(" ", file = output_file)

    print(fewer_mutations_table)

    fewer_mutations_table.to_csv("nGRE_parse_output/" + treatment + "/zero_to_one_mutations_" + regulation + "_gene_output.csv")

def main():

	start = time.time()
	fewer_mutations_parse("OVX_ADX_F_DexvsOVX_ADX_F_Veh", "downregulated")
	end = time.time()

	runtime = end - start
	print("Total runtime: {} seconds".format(runtime))


if __name__ == "__main__":
    main()
