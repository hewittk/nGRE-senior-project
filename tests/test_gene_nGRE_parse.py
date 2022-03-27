import pytest
import gene_nGRE_parse

def test_regexSearch():
    """Test regexSearch's ability to find and return matches to the nGRE consensus sequence within a gene sequence."""
    # test on generic sequence containg nGREs
    sample_sequence = "TAGTAGCATGGGATACAGCTCCGGGAGATAGCCTGATCATGGGCTCCAAGGAGAATGATCCAGGAGA"
    sample_sequence_df = gene_nGRE_parse.regexSearch(sample_sequence, "sample", "chr18", "+")
    print(sample_sequence_df.values)
    assert "CTCCGGGAGA" or "ctccgggaga" in sample_sequence_df.values
    assert "CTCCAAGGAGA" or "ctccaaggaga" in sample_sequence_df.values
    assert "ATCCAGGAGA" or "atccaggaga" in sample_sequence_df.values

    # test tool's ability to find known nGRE in TSLP sequence
    TSLP_file = open("test_gene_sequences/TSLP.txt")
    TSLP_sequence = TSLP_file.read()
    TSLP_df = gene_nGRE_parse.regexSearch(TSLP_sequence, "TSLP", "chr18", "+")
    print(TSLP_df.values)
    assert "CTCCAGGAGA" or "ctccaggaga" in TSLP_df.values

    # test tool's ability to find known nGRE in STAT1 sequence
    STAT1_file = open("test_gene_sequences/STAT1.txt")
    STAT1_sequence = STAT1_file.read()
    STAT1_df = gene_nGRE_parse.regexSearch(STAT1_sequence, "Stat1", "chr1", "+")
    print(STAT1_df.values)
    assert "CACCTGGAGA" or "cacctggaga" in STAT1_df.values
    assert "CTCCAGGACA" or "ctccaggaca" in STAT1_df.values


def test_nGRE_to_regex():
    """Test nGRE_to_regex's casting gene sequence's nucleotides into regex format."""

    # nGREs containing different spacer nucleotides
    IR0_nGRE = "CTCCGGAGA"
    assert gene_nGRE_parse.nGRE_to_regex(IR0_nGRE) == "[Cc][Tt][Cc][Cc][Gg][Gg][Aa][Gg][Aa]"
    IR1_nGRE = "CTCCTGGAGA"
    assert gene_nGRE_parse.nGRE_to_regex(IR1_nGRE) == "[Cc][Tt][Cc][Cc][Tt][Gg][Gg][Aa][Gg][Aa]"
    IR2_nGRE = "CTCCGAGGAGA"
    assert gene_nGRE_parse.nGRE_to_regex(IR2_nGRE) == "[Cc][Tt][Cc][Cc][Gg][Aa][Gg][Gg][Aa][Gg][Aa]"

    # nGREs containing different mutations
    mismatch_mutation_nGRE = "CTACGGAGA"
    assert gene_nGRE_parse.nGRE_to_regex(mismatch_mutation_nGRE) == "[Cc][Tt][Aa][Cc][Gg][Gg][Aa][Gg][Aa]"
    insertion_mutation_nGRE = "CTCCGGAGCA"
    assert gene_nGRE_parse.nGRE_to_regex(insertion_mutation_nGRE) == "[Cc][Tt][Cc][Cc][Gg][Gg][Aa][Gg][Cc][Aa]"
    deletion_mutation_nGRE = "CTCGGAGA"
    assert gene_nGRE_parse.nGRE_to_regex(deletion_mutation_nGRE) == "[Cc][Tt][Cc][Gg][Gg][Aa][Gg][Aa]"
