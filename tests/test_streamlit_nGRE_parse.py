import pytest
import streamlit_site.nGRE_parse as nGRE_parse

def test_nGRE_info_dataframe():
    """
    sample_sequence = "TAGTAGCATGGGATACAGCTCCGGGAGATAGCCTGATCATGGGCTCCAAGGAGAATGATCCAGGAGA"
    sample_sequence_df = nGRE_parse.nGRE_info_dataframe(sample_sequence, "sample", "chr18", "+")
    assert "CTCCGGGAGA" in sample_sequence_df.values
    assert "CTCCAAGGAGA" in sample_sequence_df.values
    assert "ATCCAGGAGA" in sample_sequence_df.values
    """
    print("Test")
