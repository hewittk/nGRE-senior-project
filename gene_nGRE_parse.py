import time
import regex
import Bio
from Bio import pairwise2
from Bio.Seq import Seq

import pandas as pd



nGRE_consensus_sequence = "CTCCGGAGA"

def localAlignment(gene_sequence, gene_id):
	"""Perform local alignment on gene sequences."""

	alignment_start = time.time()
	# local alignment with match score of 3, mismatch penalty of -1.5, gap opening penalty of -.2, gap extension penalty of -.2
	myAlignments = pairwise2.align.localms(gene_sequence, nGRE_consensus_sequence, 3, -1, -.2, -.2)

	with open("nGRE_sample_parse_output/parse_" + gene_id + ".txt", "a") as file:
		for nGRE_alignment in myAlignments:
			file.write(str(nGRE_alignment))

	alignment_end = time.time()
	alignment_runtime = alignment_end - alignment_start
	print("Local alignment finished, total runtime: {} seconds".format(alignment_runtime))


def regexSearch(gene_sequence, gene_id):
	"""Find sequences and positions of relative matches to nGRE consensus sequence in gene."""

	regex_start = time.time()
	# find matches
	possible_matches = regex.findall("([Cc][Tt][Cc][Cc][TAGCtagc]?[TAGCtagc]?[TAGCtagc]?[Gg][Gg][Aa][Gg][Aa]){e<=3}", gene_sequence)

	# print(possible_matches)

	# obtain information about location and number of errors in each nGRE
	for nGRE_sequence in possible_matches:

		# obtain potential nGRE sequence's location in gene sequence
		pos_information = regex.search(nGRE_sequence, gene_sequence)
		print(pos_information)

		# retrieve and store start/end coordinates of nGRE
		start_coordinate_regex = regex.findall(r"\(\d+,", str(pos_information))
		start_coordinate_element = start_coordinate_regex[0]
		# print("start coordinate: ", start_coordinate_element[1:6])

		end_coordinate_regex = regex.findall(r" \d+\)", str(pos_information))
		end_coordinate_element = end_coordinate_regex[0]
		# print("end coordinate: ", end_coordinate_element[1:6])

		# compare to perfect nGRE for error counts
		comparison = (regex.search(r"([Cc][Tt][Cc][Cc][TAGCtagc]?[TAGCtagc]?[TAGCtagc]?[Gg][Gg][Aa][Gg][Aa]){e<=3}", nGRE_sequence))
		# print(comparison.fuzzy_counts)
		mutation_counts = comparison.fuzzy_counts
		mismatch_count = mutation_counts[0]
		insertion_count = mutation_counts[1]
		deletion_count = mutation_counts[2]
		# print(mismatch_count, insertion_count, deletion_count)

		# print()

	regex_end = time.time()
	regex_runtime = regex_end - regex_start
	print("Regex search completed, total runtime: {} seconds".format(regex_runtime))

def csv_parse(treatment, regulation):
	genes_df = pd.read_csv("annotated_gene_datasets/" + treatment + "/" + regulation + "_genes/" + regulation + "_output_with_gene_names.tsv", sep = "\t")
	print("Test")

	for column in range(len(genes_df)):
		print(column)
		gene_id = genes_df["mm10.kgXref.geneSymbol"][column] + "_" + genes_df["#mm10.knownGene.name"][column]
		print(gene_id)

		file = open("gene_sequences/" + gene_id + "/" + gene_id + ".txt")
		interest_gene_sequence = file.read()

		regexSearch(interest_gene_sequence, gene_id)

def main():
	start = time.time()

	file = open("gene_sequences/Ccnd2_ENSMUST00000201637.1/Ccnd2_ENSMUST00000201637.1.txt")
	gene_id = "Ccnd2_ENSMUST00000201637.1"
	interest_gene = file.read()

	# localAlignment(interest_gene, gene_id)
	regexSearch(interest_gene, gene_id)

	csv_parse("ADX_F_DexvsADX_F_Veh", "upregulated")

	end = time.time()
	runtime = end - start
	print("Total runtime: {} seconds".format(runtime))

if __name__ == "__main__":
    main()
