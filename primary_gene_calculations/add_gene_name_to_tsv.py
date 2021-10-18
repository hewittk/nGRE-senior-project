import pandas as pd

ucsc_downregulated_data = pd.read_csv("ADX_F_DexvsADX_F_Veh/downregulated_genes/ucsc_downregulated_output.tsv", sep = "\t")
names_to_enmust_ids = pd.read_csv("ADX_F_DexvsADX_F_Veh/downregulated_genes/downregulated_gene_ENMUST_IDs.tsv", sep = "\t")

name_dict = {}
# make dictionary of ENMUST IDs mapped to gene names
for index, row in names_to_enmust_ids.iterrows():
    name_dict[row['To']] = row['From']

printed = False
found_count = 0
not_found_count = 0
for index, row in ucsc_downregulated_data.iterrows():

    # get non-decimal ENMUST ID
    ENMUST_ID_split = row['ENMUST_ID'].split(".")
    ENMUST_ID = ENMUST_ID_split[0]
    if ENMUST_ID in name_dict:
        print("Found: ", ENMUST_ID)
        found_count += 1
    else:
        print("Couldn't find: ", ENMUST_ID)
        not_found_count += 1

print("Percent not found: ", not_found_count/found_count)


"""
# make ordered list of gene names in ucsc data
gene_names = []
for index, row in ucsc_downregulated_data.iterrows():
    ENMUST_ID = int(row['ENMUST_ID'])
    if gene_names.has_key(ENMUST_ID):
        print("Found: ", ENMUST_ID)
    else:
        print("Couldn't find: ", ENMUST_ID)
"""
