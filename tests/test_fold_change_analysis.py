import pytest
import fold_change_analysis
import pandas as pd

def test_find_unique_genes():
    gene_list1 = ["cyx", "abx", "abc", "sat1", "lyrc2", "125"]
    gene_list2 = ["zyc", "sat1", "lyrc2", "abx", "xyz"]

    assert fold_change_analysis.find_unique_genes(gene_list1, gene_list2) == ["cyx", "abc", "125"]

def test_up_down_amounts():
    sample_gene_data = [["cyx", 0.235, 0.01], ["agre", 0.001, 0], ["ipjg", 0.9999, .00005], ["nvid", 1.001, 3], ["yentl", 100, 0.00003], ["tral2", -5, 0.02], ["tal3", 1.35, 0.00000035], ["pront2", 0.75, 0.0002], ["cral", 0.5, 0.0001]]

    sample_df = pd.DataFrame(sample_gene_data, columns = ['Gene_Name', 'Fold_Change', 'pvalue'])

    genes_returned = fold_change_analysis.up_down_amounts(sample_df)
    assert genes_returned[0] == ["yentl", "tal3"]
    assert genes_returned[1] == ["agre", "ipjg", "pront2", "cral"]

def test_find_mutual_genes():
    list1 = ["cyx", "agre", "ipjg", "tal3"]
    list2 = ["pront2", "cyx", "cral", "agre"]

    assert fold_change_analysis.find_mutual_genes(list1, list2) == ["cyx", "agre"]
