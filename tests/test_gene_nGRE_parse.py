import pytest
import gene_nGRE_parse

def test_nGRE_to_regex():
    TSLP_file = open("test_gene_sequences/TSLP.txt")
    TSLP_sequence = TSLP_file.read()
    TSLP_df = gene_nGRE_parse.regexSearch(TSLP_sequence, "TSLP", "chr18", "+")
    print(TSLP_df.values)
    assert "CTCCAGGAGA" in TSLP_df.values
