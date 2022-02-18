import csv
import time
import pandas as pd

def fewer_mutations_parse(treatment, regulation, maximum_mutations):
    df = pd.read_csv("nGRE_parse_output/" + treatment + "/" + regulation + "_gene_output.csv")

    two_mutations = 0
    minimal_mutations = 0

    fewer_mutations_table = pd.DataFrame(columns = ["gene_id", "chromosome", "nGRE_sequence", "start_site", "end_site", "mutations", "mismatch_mutations", "insertion_mutations", "deletion_mutations"])

    zero_mutations = df.loc[df['mutations'] == 0]
    fewer_mutations_table = fewer_mutations_table.append(zero_mutations, ignore_index = True)
    one_mutation = df.loc[df['mutations'] == 1]
    fewer_mutations_table = fewer_mutations_table.append(one_mutation, ignore_index = True)

    print(fewer_mutations_table)

    fewer_mutations_table.to_csv("nGRE_parse_output/" + treatment + "/zero_to_one_mutations_" + regulation + "_gene_output.csv")

def main():

	start = time.time()
	fewer_mutations_parse("ADX_F_DexvsADX_F_Veh", "upregulated", 1)
	end = time.time()

	runtime = end - start
	print("Total runtime: {} seconds".format(runtime))


if __name__ == "__main__":
    main()
