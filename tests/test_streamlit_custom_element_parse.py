import pytest
import streamlit_site.custom_element_parse as custom_element_parse

def test_nucleotide_bracket():
    nGRE_halfsite1 = "CTCC"
    assert custom_element_parse.nucleotide_bracket(nGRE_halfsite1) == "[Cc][Tt][Cc][Cc]"
    nGRE_halfsite2 = "GGAGA"
    assert custom_element_parse.nucleotide_bracket(nGRE_halfsite2) == "[Gg][Gg][Aa][Gg][Aa]"

    GRE_halfsite1 = "GGTACA"
    assert custom_element_parse.nucleotide_bracket(GRE_halfsite1) == "[Gg][Gg][Tt][Aa][Cc][Aa]"
    GRE_halfsite2 = "TGTTCT"
    assert custom_element_parse.nucleotide_bracket(GRE_halfsite2) == "[Tt][Gg][Tt][Tt][Cc][Tt]"

    sample_element = "AcagtacccAAAT"
    assert custom_element_parse.nucleotide_bracket(sample_element) == "[Aa][Cc][Aa][Gg][Tt][Aa][Cc][Cc][Cc][Aa][Aa][Aa][Tt]"

def test_element_to_regex():
    nGRE_sequence = "CTCC(N)(0-2)GGAGA"
    nGRE_sequence_regex = str(custom_element_parse.element_to_regex(nGRE_sequence))
    assert nGRE_sequence_regex == "[Cc][Tt][Cc][Cc]([AaCcTtGg]){0,2}[Gg][Gg][Aa][Gg][Aa]"

    GRE_sequence = "GGTACA(N)(3)TGTTCT"
    GRE_sequence_regex = str(custom_element_parse.element_to_regex(GRE_sequence))
    assert GRE_sequence_regex == "[Gg][Gg][Tt][Aa][Cc][Aa]([AaCcTtGg]){3}[Tt][Gg][Tt][Tt][Cc][Tt]"

    sample_seq = "ATAG(CA)(4-6)AYG"
    sample_seq_regex = str(custom_element_parse.element_to_regex(sample_seq))
    assert sample_seq_regex == "[Aa][Tt][Aa][Gg]([Cc][Aa]){4,6}[Aa][CcTt][Gg]"

def test_sequence_search():
    sample_nGRE_search = custom_element_parse.sequence_search("TAGTAGCATGGGATACAGCTCCGGGAGATAGCCTGATCATGGGCTCCAAGGAGAATGATCCAGGAGA", "[Cc][Tt][Cc][Cc]([AaCcTtGg]){0,2}[Gg][Gg][Aa][Gg][Aa]", 0)
    assert "CTCCGGGAGA" in sample_nGRE_search
    assert "CTCCAAGGAGA" in sample_nGRE_search
    assert "ATCCAGGAGA" not in sample_nGRE_search
