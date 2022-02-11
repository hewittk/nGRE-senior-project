import csv
import pandas as pd

df = pd.read_csv('nGRE_parse_output/ADX_F_DexvsADX_F_Veh/downregulated_gene_output.csv')

two_mutations = 0
minimal_mutations = 0

fewer_mutations_table = pd.DataFrame(columns = ["gene_id", "chromosome", "nGRE_sequence", "start_site", "end_site", "mutations", "mismatch_mutations", "insertion_mutations", "deletion_mutations"])

zero_mutations = df.loc[df['mutations'] == 0]
fewer_mutations_table = fewer_mutations_table.append(zero_mutations, ignore_index = True)
one_mutation = df.loc[df['mutations'] == 1]
fewer_mutations_table = fewer_mutations_table.append(one_mutation, ignore_index = True)

print(fewer_mutations_table)

fewer_mutations_table.to_csv("nGRE_parse_output/ADX_F_DexvsADX_F_Veh/zero_to_one_mutations_downregulated_gene_output.csv")
