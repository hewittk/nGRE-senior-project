"""Calculate amounts of genes downregulated and upregulated in initial DEG datasets."""

import pandas as pd
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt

minimum_p_value = 0.001

dex_vs_vehicle_df = pd.read_csv("primary_DEG_datasets/ADX_F_DexvsADX_F_Veh_results_STAR_FeatureCount_DESeq2_annotation_p.01.xlsx - ADX_F_DexvsADX_F_Veh_results_ST.csv")
ovx_dex_vs_vehicle_df = pd.read_csv("primary_DEG_datasets/OVX_ADX_F_DexvsOVX_ADX_F_Veh_results_STAR_FeatureCount_DESeq2_annotation_p.01.xlsx - OVX_ADX_F_DexvsOVX_ADX_F_Veh_re.csv")

def up_down_amounts(DEG_dataframe):
    """Return list of genes that are upregulated/downregulated based on fold change."""

    upregulated_genes = []
    downregulated_genes = []

    for index, row in DEG_dataframe.iterrows():
        if(row['pvalue']) < minimum_p_value:
            if(row['Fold_Change']) > 1:
                upregulated_genes.append(row['Gene_Name'])
            elif(row['Fold_Change']) < 1:
                downregulated_genes.append(row['Gene_Name'])

    return upregulated_genes, downregulated_genes


def write_to_file(treatment, upregulated_genes, downregulated_genes):
    """Write given lists of genes that are upregulated/downregulated to files."""

    with open("primary_gene_calculations/" + treatment + "/upregulated.txt", "w") as file:
        for gene_name in upregulated_genes:
            file.write(gene_name)
            file.write("\n")
        file.close()
    with open("primary_gene_calculations/" + treatment + "/downregulated.txt", "w") as file:
        for gene_name in downregulated_genes:
            file.write(gene_name)
            file.write("\n")
        file.close()

def find_mutual_genes(list1, list2):
    """Return list of genes common between two given gene lists."""

    mutual_genes = []

    for list1_gene in list1:
        for list2_gene in list2:
            if(list1_gene == list2_gene):
                mutual_genes.append(list1_gene)

    return mutual_genes

def find_unique_genes(list1, list2):
    """Return list of unique genes between two given gene lists."""
    unique_list1_genes = []

    for list1_gene in list1:
        if not(list1_gene in list2):
            unique_list1_genes.append(list1_gene)

    return unique_list1_genes


def main():
    # write upregulated and downregulated gene names to files
    dex_upregulated_genes, dex_downregulated_genes = up_down_amounts(dex_vs_vehicle_df)
    write_to_file("ADX_F_DexvsADX_F_Veh", dex_upregulated_genes, dex_downregulated_genes)

    ovx_dex_upregulated_genes, ovx_dex_downregulated_genes = up_down_amounts(ovx_dex_vs_vehicle_df)
    write_to_file("OVX_ADX_F_DexvsOVX_ADX_F_Veh", ovx_dex_upregulated_genes, ovx_dex_downregulated_genes)

    # Report amount of differential expression before/after estrogen depletion
    print("Before estrogen depletion, total genes upregulated by dexamethasone: ", len(dex_upregulated_genes))
    print("Before estrogen depletion, total genes downregulated by dexamethasone: ", len(dex_downregulated_genes))
    print()

    print("After estrogen depletion, total genes upregulated by dexamethasone: ", len(ovx_dex_upregulated_genes))
    print("After estrogen depletion, total genes downregulated by dexamethasone: ", len(ovx_dex_downregulated_genes))
    print()

    # Report numbers of genes upregulated and downregulated in each group
    mutual_upregulated_genes = find_mutual_genes(dex_upregulated_genes, ovx_dex_upregulated_genes)
    mutual_downregulated_genes = find_mutual_genes(dex_downregulated_genes, ovx_dex_downregulated_genes)

    print("Genes upregulated by dexamethasone both before and after estrogen depletion: ", len(mutual_upregulated_genes))
    print("Genes downregulated by dexamethasone both before and after estrogen depletion: ", len(mutual_downregulated_genes))
    print()

    ovx_upregulated_genes = find_unique_genes(ovx_dex_upregulated_genes, dex_upregulated_genes)
    ovx_downregulated_genes = find_unique_genes(ovx_dex_downregulated_genes, dex_downregulated_genes)
    non_ovx_upregulated_genes = find_unique_genes(dex_upregulated_genes, ovx_dex_upregulated_genes)
    non_ovx_downregulated_genes = find_unique_genes(dex_downregulated_genes, ovx_dex_downregulated_genes)

    print("Genes only upregulated by dexamethasone after estrogen depletion: ", len(ovx_upregulated_genes))
    print("Genes only downregulated by dexamethasone after estrogen depletion:  ", len(ovx_downregulated_genes))
    print("Genes only upregulated by dexamethasone before estrogen depletion: ", len(non_ovx_upregulated_genes))
    print("Genes only downregulated by dexamethasone before estrogen depletion: ", len(non_ovx_downregulated_genes))

    # create Venn diagrams
    venn2(subsets = (len(ovx_upregulated_genes), len(non_ovx_upregulated_genes), len(mutual_upregulated_genes)), set_labels = ('OVX_ADX', 'ADX'), set_colors=('purple', 'skyblue'), alpha = 0.5)
    plt.savefig("images/upregulated.png")
    plt.show()

    venn2(subsets = (len(ovx_downregulated_genes), len(non_ovx_downregulated_genes), len(mutual_downregulated_genes)), set_labels = ('OVX_ADX', 'ADX'), set_colors=('purple', 'skyblue'), alpha = 0.5)
    plt.savefig("images/downregulated.png")
    plt.show()


if __name__ == "__main__":
    main()
