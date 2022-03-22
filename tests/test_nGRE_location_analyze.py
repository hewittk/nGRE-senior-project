import pytest
import nGRE_location_analyze

def test_retrieve_gene_sites():
    Pcmdt1_txStart, Pcmdt1_cdsStart, Pcmdt1_cdsEnd, Pcmdt1_txEnd = nGRE_location_analyze.retrieve_gene_sites("Pcmtd1_ENSMUST00000182306.1", "ADX_F_DexvsADX_F_Veh", "upregulated")
    assert Pcmdt1_txStart == 7089103
    assert Pcmdt1_cdsStart == 7089103
    assert Pcmdt1_cdsEnd == 7089103
    assert Pcmdt1_txEnd == 7160902

    Mybl1_txStart, Mybl1_cdsStart, Mybl1_cdsEnd, Mybl1_txEnd = nGRE_location_analyze.retrieve_gene_sites("Mybl1_ENSMUST00000115468.8", "OVX_ADX_F_DexvsOVX_ADX_F_Veh", "downregulated")
    assert Mybl1_txStart == 9669536
    assert Mybl1_cdsStart == 9669883
    assert Mybl1_cdsEnd == 9699954
    assert Mybl1_txEnd == 9700024

def test_classify_sites():
    gene_id = "test"
    gene_txStart = 6245
    gene_cdStart = 6250
    gene_cdEnd = 7540
    gene_txEnd = 8000

    relative_gene_txStart, relative_gene_cdStart, relative_gene_cdEnd, relative_gene_txEnd = nGRE_location_analyze.classify_sites("test", "upregulated", "OVX_ADX_F_DexvsOVX_ADX_F_Veh", gene_txStart, gene_cdStart, gene_cdEnd, gene_txEnd)

    assert relative_gene_txStart == 30000
    assert relative_gene_cdStart == 30005
    assert relative_gene_cdEnd == 31295
    assert relative_gene_txEnd == 31755

def test_cumulate_sites():
    nGRE_location_analyze.cumulate_sites(1000, 30000, 30000, 50000, 50006)
    print(nGRE_location_analyze.nGRE_locations["Pre transcription start site"])
    assert nGRE_location_analyze.nGRE_locations["Pre transcription start site"] == 1

    nGRE_location_analyze.cumulate_sites(32000, 30000, 31000, 46000, 46000)
    assert nGRE_location_analyze.nGRE_locations["After coding site"] == 1

    nGRE_location_analyze.cumulate_sites(48000, 30000, 31000, 46000, 46000)
    nGRE_location_analyze.cumulate_sites(30005, 30000, 31000, 46000, 46000)
    nGRE_location_analyze.cumulate_sites(29000, 30000, 31000, 46000, 46000)
    nGRE_location_analyze.cumulate_sites(48000, 30000, 31000, 46000, 46000)
    nGRE_location_analyze.cumulate_sites(80000, 30000, 31000, 46000, 46000)
    nGRE_location_analyze.cumulate_sites(45999, 30000, 31000, 46000, 46000)
    nGRE_location_analyze.cumulate_sites(25000, 30000, 31000, 46000, 46000)

    assert nGRE_location_analyze.nGRE_locations["Pre transcription start site"] == 3
