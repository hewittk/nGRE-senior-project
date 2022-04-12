import pandas as pd
import csv

def csv_parse(treatment, regulation):
    potential_nGRE_file = "nGRE_parse_output/" + treatment + "/" + regulation + "_gene_output.csv"

    promoter_nGRE_sites = pd.DataFrame(columns = ["gene_id", "chromosome", "nGRE_sequence", "start_site", "end_site", "mutations", "mismatch_mutations", "insertion_mutations", "deletion_mutations"])
    nGRE_dict = {}

    with open(potential_nGRE_file, 'r') as potential_nGREs:
        data = csv.reader(potential_nGREs)
        for row in data:
            if(row[0] != "gene_id"): # if isn't first row
                nGRE_start = row[3] # nGRE start site within retrieved gene
                if (int(nGRE_start) < 35000):
                    print(row)
                    nGRE_dict["gene_id"] = row[2]
                    nGRE_dict["chromosome"] = row[3]
                    nGRE_dict["nGRE_sequence"] = row[4]
                    nGRE_dict["start_site"] = row[0]
                    nGRE_dict["end_site"] = row[1]
                    nGRE_dict["mutations"] = row[5]
                    nGRE_dict["mismatch_mutations"] = row[6]
                    nGRE_dict["insertion_mutations"] = row[7]
                    nGRE_dict["deletion_mutations"] = row[8]

                    previous_gene_id = row[2]

                    promoter_nGRE_sites = promoter_nGRE_sites.append(nGRE_dict, ignore_index = True)

    return promoter_nGRE_sites

def main():
    treatment = "ADX_F_DexvsADX_F_Veh"
    regulation = "downregulated"
    potential_nGREs = csv_parse(treatment, regulation)

    potential_nGREs.to_csv("nGRE_parse_output/" + treatment + "/" + regulation + "_promoter_gene_output.csv", mode="w", index=False, header=False)


if __name__ == "__main__":
    main()
