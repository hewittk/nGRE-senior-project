import time
import Bio
from Bio import pairwise2
from Bio.Seq import Seq


nGRE_consensus_sequence = "CTCCGGAGA"

def localAlignment(gene_sequence, gene_id):
	# local alignment with match score of 3, mismatch penalty of -1.5, gap opening penalty of -.5, gap extension penalty of -5
	myAlignments = pairwise2.align.localms(gene_sequence, nGRE_consensus_sequence, 3, -1.5, -.5, -.5)

	with open("nGRE_sample_parse_output/parse_" + gene_id + ".txt", "w") as file:
		for nGRE_alignment in myAlignments:
			file.write(nGRE_alignment)


def main():
	start = time.time()

	file = open("gene_sequences/Anapc16_ENSMUST00000182912.1/Anapc16_ENSMUST00000182912.1.txt")
	gene_id = "Anapc16_ENSMUST00000182912.1"
	interest_gene = file.read()
	localAlignment(interest_gene, gene_id)

	end = time.time()
	runtime = end - start
	print("Total runtime: {} seconds".format(runtime))

if __name__ == "__main__":
    main()
