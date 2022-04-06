import pytest
import streamlit_site.custom_element_parse as custom_element_parse

def test_nucleotide_bracket():
    """Test casting of each nucleotide/letter in a string into bracket containing it in upper/lower case letters."""

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
    """Test conversion of an element containing nucleotides/IUPAC codes/numbers into its regular expression equivalent."""

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
    """Test that sequence_search properly finds a given regular expression with up to a certain amount of errors within a sequence."""

    m0_nGRE_search = custom_element_parse.sequence_search("TAGTAGCATGGGATACAGCTCCGGGAGATAGCCTGATCATGGGCTCCAAGGAGAATGATCCAGGAGA", "[Cc][Tt][Cc][Cc]([AaCcTtGg]){0,2}[Gg][Gg][Aa][Gg][Aa]", 0)
    assert "CTCCGGGAGA" in m0_nGRE_search
    assert "CTCCAAGGAGA" in m0_nGRE_search
    assert "ATCCAGGAGA" not in m0_nGRE_search # assert doesn't return mutated nGRE in sequence

    m1_nGRE_search = custom_element_parse.sequence_search("TAGTAGCATGGGATACAGCTCCGGGAGATAGCCTGATCATGGGCTCCAAGGAGAATGATCCAGGAGA", "[Cc][Tt][Cc][Cc]([AaCcTtGg]){0,2}[Gg][Gg][Aa][Gg][Aa]", 1)
    assert "CTCCGGGAGA" in m1_nGRE_search
    assert "CTCCAAGGAGA" in m1_nGRE_search
    assert "ATCCAGGAGA" in m1_nGRE_search

def test_iupac_handling():
    """Test that iupac_handling replaces IUPAC codes with their nucleotide equivalents in the appropriate position of a gene sequence."""

    sample_purine_seq = custom_element_parse.iupac_handling("R",9,"[Aa][Cc][Rr]{1,3}[Gg]")
    assert sample_purine_seq == "[Aa][Cc][AaGg]{1,3}[Gg]"

    sample_n_seq = custom_element_parse.iupac_handling("N", 5, "[Cc][Nn][Aa][Gg][Gg][Gg]")
    assert sample_n_seq == "[Cc][AaCcTtGg][Aa][Gg][Gg][Gg]"

def test_count_length():
    """Test that count_length properly returns the amount of non-optional characters in a regular expression."""

    assert custom_element_parse.count_length("[Cc][Tt][Cc][Cc]([AaCcTtGg]){0,2}[Gg][Gg][Aa][Gg][Aa]") == 9
    assert custom_element_parse.count_length("[Aa]?[Gg][Tt][Cc][Aa][Cc]?[Aa][Gg]?") == 5
    assert custom_element_parse.count_length("[Aa][Gg][Cc][Aa]") == 4
