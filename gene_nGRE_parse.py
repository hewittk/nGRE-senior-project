import time
import regex
import Bio
from Bio import pairwise2
from Bio.Seq import Seq

import pandas as pd
<<<<<<< HEAD

nGRE_consensus_sequence = "CTCCGGAGA"


=======



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


>>>>>>> 2186d65a55acd86563de23c580c4affbb6f97760
def regexSearch(gene_sequence, gene_id, chromosome, strand):
	"""Find sequences and positions of relative matches to nGRE consensus sequence in gene."""

	regex_start = time.time()
	# find matches
	print("Searching for nGREs in ", gene_id)
	possible_matches = regex.findall("([Cc][Tt][Cc][Cc][TAGCtagc]?[TAGCtagc]?[TAGCtagc]?[Gg][Gg][Aa][Gg][Aa]){e<=2}", gene_sequence)
<<<<<<< HEAD
=======
	# print(possible_matches)

	# print(possible_matches)
>>>>>>> 2186d65a55acd86563de23c580c4affbb6f97760

	# obtain information about location and number of errors in each nGRE
	nGRE_sites = pd.DataFrame(columns = ["gene_id", "chromosome", "nGRE_sequence", "start_site", "end_site", "mismatch_mutations", "insertion_mutations", "deletion_mutations"])

	nGRE_list = []

	for nGRE_sequence in possible_matches:
		nGRE_information_dict = {}

		# obtain potential nGRE sequence's location in gene sequence
		pos_information = regex.search(nGRE_sequence, gene_sequence)

		# retrieve and store start/end coordinates of nGRE
		start_coordinate_regex = regex.findall(r"\(\d+,", str(pos_information))[0]
		start_coordinate = int (regex.findall("\d+", str(start_coordinate_regex))[0])
<<<<<<< HEAD
		end_coordinate_regex = regex.findall(r" \d+\)", str(pos_information))[0]
		end_coordinate = int (regex.findall("\d+", str(end_coordinate_regex))[0])

		# compare to perfect nGRE for error counts
		comparison = (regex.search(r"((?e)[Cc][Tt][Cc][Cc][TAGCtagc]?[TAGCtagc]?[TAGCtagc]?[Gg][Gg][Aa][Gg][Aa]){e<=2}", nGRE_sequence))
=======

		# start_coordinate = start_coordinate_element[1:len(start_coordinate_element)]

		end_coordinate_regex = regex.findall(r" \d+\)", str(pos_information))[0]
		end_coordinate = int (regex.findall("\d+", str(end_coordinate_regex))[0])

		# compare to perfect nGRE for error counts
		comparison = (regex.search(r"([Cc][Tt][Cc][Cc][TAGCtagc]?[TAGCtagc]?[TAGCtagc]?[Gg][Gg][Aa][Gg][Aa]){e<=2}", nGRE_sequence))
		# print(comparison.fuzzy_counts)
>>>>>>> 2186d65a55acd86563de23c580c4affbb6f97760
		mutation_counts = comparison.fuzzy_counts
		mismatch_count = mutation_counts[0]
		insertion_count = mutation_counts[1]
		deletion_count = mutation_counts[2]
<<<<<<< HEAD
		total_mutations = mismatch_count + insertion_count + deletion_count
=======
		# print(mismatch_count, insertion_count, deletion_count)
>>>>>>> 2186d65a55acd86563de23c580c4affbb6f97760

		# assemble dictionary of nGRE information
		nGRE_information_dict["gene_id"] = gene_id
		nGRE_information_dict["chromosome"] = chromosome
		nGRE_information_dict["nGRE_sequence"] = nGRE_sequence
		nGRE_information_dict["start_site"] = start_coordinate
		nGRE_information_dict["end_site"] = end_coordinate
<<<<<<< HEAD
		nGRE_information_dict["mutations"] = total_mutations
=======
>>>>>>> 2186d65a55acd86563de23c580c4affbb6f97760
		nGRE_information_dict["mismatch_mutations"] = mismatch_count
		nGRE_information_dict["insertion_mutations"] = insertion_count
		nGRE_information_dict["deletion_mutations"] = deletion_count

		# nGRE_sites = nGRE_sites.append(nGRE_information_dict, ignore_index = True)
<<<<<<< HEAD
=======
		# print(nGRE_information_dict)
		# print(nGRE_sites)
		# print(nGRE_information_dict)
		print()
>>>>>>> 2186d65a55acd86563de23c580c4affbb6f97760
		nGRE_list.append(nGRE_information_dict)
		# nGRE_information_dict.clear()

	regex_end = time.time()
	regex_runtime = regex_end - regex_start
	print("Regex search completed, total runtime: {} seconds".format(regex_runtime))
<<<<<<< HEAD

	nGRE_table = pd.DataFrame(columns = ["gene_id", "chromosome", "nGRE_sequence", "start_site", "end_site", "mutations", "mismatch_mutations", "insertion_mutations", "deletion_mutations"])
	for nGRE in nGRE_list:
		nGRE_table = nGRE_table.append(nGRE, ignore_index = True)

	return nGRE_table
=======
	# print(nGRE_sites)
	print()

	return nGRE_list
>>>>>>> 2186d65a55acd86563de23c580c4affbb6f97760


def csv_parse(treatment, regulation):
	genes_df = pd.read_csv("annotated_gene_datasets/" + treatment + "/" + regulation + "_genes/" + regulation + "_output_with_gene_names_known_sequences.tsv", sep = "\t")
<<<<<<< HEAD


	for column in range(len(genes_df)):
		gene_id = genes_df["mm10.kgXref.geneSymbol"][column] + "_" + genes_df["#mm10.knownGene.name"][column]
=======
	print("Test")

	nGRE_table = pd.DataFrame(columns = ["gene_id", "chromosome", "nGRE_sequence", "start_site", "end_site", "mismatch_mutations", "insertion_mutations", "deletion_mutations"])

	for column in range(0, 10):
		# print(column)
		gene_id = genes_df["mm10.kgXref.geneSymbol"][column] + "_" + genes_df["#mm10.knownGene.name"][column]
		# print(gene_id)
>>>>>>> 2186d65a55acd86563de23c580c4affbb6f97760
		gene_chromosome = genes_df["mm10.knownGene.chrom"][column]
		gene_strand = genes_df["mm10.knownGene.strand"][column]

		file = open("gene_sequences/" + gene_id + "/" + gene_id + ".txt")
		interest_gene_sequence = file.read()

<<<<<<< HEAD
		potential_nGREs = regexSearch(interest_gene_sequence, gene_id, gene_chromosome, gene_strand)
		potential_nGREs.to_csv("test.csv", mode="a", index=False, header=False)
=======
		"""
		potential_nGRE_dataframe = regexSearch(interest_gene_sequence, gene_id, gene_chromosome, gene_strand)
		print(potential_nGRE_dataframe)
		nGRE_table.append(potential_nGRE_dataframe, ignore_index=True)
		"""
		# print(nGRE_table)
		# print()
		potential_nGREs = regexSearch(interest_gene_sequence, gene_id, gene_chromosome, gene_strand)
		for potential_nGRE in potential_nGREs:
			nGRE_table = nGRE_table.append(potential_nGRE, ignore_index = True)

	nGRE_table.to_csv("test.csv")
>>>>>>> 2186d65a55acd86563de23c580c4affbb6f97760


def main():
	start = time.time()

	file = open("gene_sequences/Ccnd2_ENSMUST00000201637.1/Ccnd2_ENSMUST00000201637.1.txt")
	gene_id = "Ccnd2_ENSMUST00000201637.1"
	interest_gene = file.read()

<<<<<<< HEAD
=======
	# localAlignment(interest_gene, gene_id)
	# regexSearch(interest_gene, gene_id)

>>>>>>> 2186d65a55acd86563de23c580c4affbb6f97760
	csv_parse("ADX_F_DexvsADX_F_Veh", "upregulated")

	end = time.time()
	runtime = end - start
	print("Total runtime: {} seconds".format(runtime))

if __name__ == "__main__":
    main()
