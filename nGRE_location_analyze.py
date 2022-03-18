import pandas as pd
import csv
import regex

# retrieve gene's chr, trx start site, coding start site, coding end site, and trx end site
def retrieve_gene_sites(gene_id, treatment, regulation, nGRE_start, nGRE_end):
    genes_df = pd.read_csv("annotated_gene_datasets/" + treatment + "/" + regulation + "_genes/" + regulation + "_output_with_gene_names_known_sequences.tsv", sep = "\t")

    separated_id = regex.search(r'_\w+\.?\d+', gene_id)
    enmust_id = separated_id.group(0)[1:]

    # retrieve start sites within genes
    # print("enmust_id: ", enmust_id)
    enmust_id_info = genes_df.loc[genes_df['#mm10.knownGene.name'] == enmust_id]
    start = enmust_id_info['mm10.knownGene.txStart']
    print(start.values)
    # print("gene_start: " + gene_start)

    # print("Test")

# retreive nGRE's chromosome and start site from nGRE parse output file
def retrieve_nGRE_sites(gene_id, treatment, regulation, nGRE_start, nGRE_end):
    print("Test")

# classify what regions of gene nGREs were found in
def classify_sites(gene_id, treatment, regulation, nGRE_start, nGRE_end):
    genes_df = pd.read_csv("annotated_gene_datasets/" + treatment + "/" + regulation + "_genes/" + regulation + "_output_with_gene_names_known_sequences.tsv", sep = "\t")

    separated_id = regex.search(r'_\w+\.?\d', gene_id)
    enmust_id = separated_id.group(0)[1:]

    print(enmust_id)


def main():
    treatment = "ADX_F_DexvsADX_F_Veh"
    regulation = "upregulated"
    potential_nGRE_file = "nGRE_parse_output/" + treatment + "/" + regulation + "_gene_output.csv"

    previous_gene_id = ""
    with open(potential_nGRE_file, 'r') as potential_nGREs:
        data = csv.reader(potential_nGREs)
        for row in data:
            gene_id = row[0]
            nGRE_start = row[3] # nGRE start site within retrieved gene
            nGRE_end = row[4] # nGRE end site within retrieved gene

            if(gene_id == previous_gene_id):
                retrieve_gene_sites(gene_id, treatment, regulation, nGRE_start, nGRE_end)
            # else:
            # count genes in each region for current gene_id
            # record and switch to new gene_id if get to next gene_id

            previous_gene_id = gene_id

if __name__ == "__main__":
    main()
