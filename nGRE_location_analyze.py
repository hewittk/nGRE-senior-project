import pandas as pd
import csv
import regex

nGRE_locations = {"Pre transcription start site":0, "Pre coding start site":0, "During coding site":0, "After coding site":0, "After transcription site":0}

def retrieve_gene_sites(gene_id, treatment, regulation):
    """Return gene's transcription and coding start/end sites on chromosome."""

    genes_df = pd.read_csv("annotated_gene_datasets/" + treatment + "/" + regulation + "_genes/" + regulation + "_output_with_gene_names_known_sequences.tsv", sep = "\t")

    # derive gene's enmust id to search for in csv file
    print("retrieve_gene_sites gene_id: " + gene_id)
    separated_id = regex.search(r'_\w+\.?\d+', gene_id)
    enmust_id = separated_id.group(0)[1:]

    # retrieve start and end sites of gene
    enmust_id_info = genes_df.loc[genes_df['#mm10.knownGene.name'] == enmust_id]
    txStart_data = enmust_id_info['mm10.knownGene.txStart']
    txStart = txStart_data.values[0]
    cdsStart_data = enmust_id_info['mm10.knownGene.cdsStart']
    cdsStart = cdsStart_data.values[0]
    cdsEnd_data = enmust_id_info['mm10.knownGene.cdsEnd']
    cdsEnd = cdsEnd_data.values[0]
    txEnd_data = enmust_id_info['mm10.knownGene.txEnd']
    txEnd = txEnd_data.values[0]

    return txStart, cdsStart, cdsEnd, txEnd

def classify_sites(gene_id, treatment, regulation, gene_txStart, gene_cdsStart, gene_cdsEnd, gene_txEnd):
    """Classify relative start/end landmark sites of gene within returned gene site instead of within chromosome."""

    relative_txStart = 30000
    relative_cdStart = 30000 + (gene_cdsStart - gene_txStart)
    relative_cdEnd = relative_cdStart + (gene_cdsEnd - gene_cdsStart)
    relative_txEnd = relative_cdEnd + (gene_txEnd - gene_cdsEnd)

    print("Gene: " + gene_id)
    print("relative_txStart: " + str(relative_txStart))
    print("relative_cdStart: " + str(relative_cdStart))
    print("relative_cdEnd: " + str(relative_cdEnd))
    print("relative_txEnd: " + str(relative_txEnd))

    return relative_txStart, relative_cdStart, relative_cdEnd, relative_txEnd

def cumulate_sites(nGRE_start, relative_txStart, relative_cdStart, relative_cdEnd, relative_txEnd):
    """Add nGRE site location to running total of nGRE site locations found in dataset."""

    if (int(nGRE_start) < relative_txStart):
        nGRE_locations["Pre transcription start site"] += 1
    elif (int(nGRE_start) < relative_cdStart):
        nGRE_locations["Pre coding start site"] += 1
    elif (int(nGRE_start) < relative_cdEnd):
        nGRE_locations["During coding site"] += 1
    elif ((int(nGRE_start) > relative_txStart) and (int(nGRE_start) < relative_txEnd)):
        nGRE_locations["After coding site"] += 1
    elif ((int(nGRE_start) > relative_txEnd) and (int(nGRE_start) < (relative_txEnd + 5000))):
        nGRE_locations["After transcription site"] += 1
    else:
        print("Error: nGRE site of " + str(nGRE_start) + " doesn't appear to be in range 0 to " + str((relative_txEnd + 5000)))

def main():
    treatment = "OVX_ADX_F_DexvsOVX_ADX_F_Veh"
    regulation = "downregulated"
    potential_nGRE_file = "nGRE_parse_output/" + treatment + "/" + regulation + "_gene_output.csv"

    previous_gene_id = ""
    gene_id = ""
    with open(potential_nGRE_file, 'r') as potential_nGREs:
        data = csv.reader(potential_nGREs)
        for row in data:
            if(row[0] != "gene_id"): # if isn't first row
                gene_id = row[0]
                nGRE_start = row[3] # nGRE start site within retrieved gene
                nGRE_end = row[4] # nGRE end site within retrieved gene
                nGRE_mutations = row[5]

                if(gene_id != previous_gene_id):
                    txStart, cdsStart, cdsEnd, txEnd = retrieve_gene_sites(gene_id, treatment, regulation)
                    relative_txStart, relative_cdStart, relative_cdEnd, relative_txEnd = classify_sites(gene_id, treatment, regulation, txStart, cdsStart, cdsEnd, txEnd)

                cumulate_sites(nGRE_start, relative_txStart, relative_cdStart, relative_cdEnd, relative_txEnd)

            previous_gene_id = gene_id

    locations_file = open("nGRE_parse_output/nGRE_location_parse_statistics.txt", "a")
    locations_file.write(treatment + " " + regulation + " genes \n")
    locations_file.write(("nGRE_locations: " + str(nGRE_locations) + "\n\n"))
    locations_file.close()


if __name__ == "__main__":
    main()
