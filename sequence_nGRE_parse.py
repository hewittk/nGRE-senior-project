import Bio
from Bio import pairwise2
from Bio.Seq import Seq

file = open("Glyr1_datasets/ncbi_dataset/data/gene.fna")
gene_sequence = file.read()

nGRE_consensus_sequence = "CTCCGGAGA"

myAlignments = pairwise2.align.globalmd(gene_sequence, nGRE_consensus_sequence)

for thisAlignment in myAlignments:
	print(thisAlignment)
