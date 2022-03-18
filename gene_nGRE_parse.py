import time
import regex
import Bio
from Bio import pairwise2
from Bio.Seq import Seq

import pandas as pd

nGRE_consensus_sequence = "CTCCGGAGA"


def regexSearch(gene_sequence, gene_id, chromosome, strand):
	"""Find sequences and positions of relative matches to nGRE consensus sequence in gene."""

	regex_start = time.time()
	# find matches
	print("Searching for nGREs in ", gene_id)
	possible_matches = regex.findall("((?e)[Cc][Tt][Cc][Cc][TAGCtagc]?[TAGCtagc]?[TAGCtagc]?[Gg][Gg][Aa][Gg][Aa]){e<=1}", gene_sequence)

	# remove duplicates in possible matches set
	nonrepeat_matches = []
	for match in possible_matches:
		if match not in nonrepeat_matches:
		    nonrepeat_matches.append(match)

	# obtain information about location and number of errors in each nGRE
	nGRE_sites = pd.DataFrame(columns = ["gene_id", "chromosome", "nGRE_sequence", "start_site", "end_site", "mutations", "mismatch_mutations", "insertion_mutations", "deletion_mutations"])
	nGRE_list = []
	nGRE_information_dict = {}
	for nGRE_sequence in nonrepeat_matches:

		# compare to perfect nGRE for error counts
		comparison = (regex.search(r"((?e)[Cc][Tt][Cc][Cc][TAGCtagc]?[TAGCtagc]?[TAGCtagc]?[Gg][Gg][Aa][Gg][Aa]){e<=1}", nGRE_sequence))
		mutation_counts = comparison.fuzzy_counts
		mismatch_count = int(mutation_counts[0])
		insertion_count = int(mutation_counts[1])
		deletion_count = int(mutation_counts[2])
		total_mutations = int(int(mismatch_count) + int(insertion_count) + int(deletion_count))

		# find all start/end locations of element
		nGRE_locations = []
		match_element_regex = regex.compile(nGRE_to_regex(nGRE_sequence))
		for match_location in match_element_regex.finditer(gene_sequence):
		    nGRE_locations.append([match_location.start(), match_location.end()])

		for location in nGRE_locations:
			nGRE_information_dict.clear()
			nGRE_information_dict["start_site"] = location[0]
			nGRE_information_dict["end_site"] = location[1]
			nGRE_information_dict["gene_id"] = gene_id
			nGRE_information_dict["chromosome"] = chromosome
			nGRE_information_dict["nGRE_sequence"] = nGRE_sequence
			nGRE_information_dict["mutations"] = int(total_mutations)
			nGRE_information_dict["mismatch_mutations"] = int(mismatch_count)
			nGRE_information_dict["insertion_mutations"] = int(insertion_count)
			nGRE_information_dict["deletion_mutations"] = int(deletion_count)
			nGRE_sites = nGRE_sites.append(nGRE_information_dict, ignore_index = True)

	regex_end = time.time()
	regex_runtime = regex_end - regex_start
	print("Regex search completed, total runtime: {} seconds".format(regex_runtime))


	# commented out while trying to replace nGRE_list
	# nGRE_table = pd.DataFrame(columns = ["gene_id", "chromosome", "nGRE_sequence", "start_site", "end_site", "mutations", "mismatch_mutations", "insertion_mutations", "deletion_mutations"])
	# for nGRE in nGRE_list:
		# nGRE_table = nGRE_table.append(nGRE, ignore_index = True)

	return nGRE_sites


def csv_parse(treatment, regulation):
	"""Parse each gene for nGREs and record results in csv."""

	genes_df = pd.read_csv("annotated_gene_datasets/" + treatment + "/" + regulation + "_genes/" + regulation + "_output_with_gene_names_known_sequences.tsv", sep = "\t")

	for column in range(len(genes_df)):
		gene_id = genes_df["mm10.kgXref.geneSymbol"][column] + "_" + genes_df["#mm10.knownGene.name"][column]
		gene_chromosome = genes_df["mm10.knownGene.chrom"][column]
		gene_strand = genes_df["mm10.knownGene.strand"][column]

		# read sequence of given gene
		file = open("gene_sequences/" + gene_id + "/" + gene_id + ".txt")
		interest_gene_sequence = file.read()

		# record potential nGREs found in given gene sequence
		potential_nGREs = regexSearch(interest_gene_sequence, gene_id, gene_chromosome, gene_strand)
		potential_nGREs.to_csv("nGRE_parse_output/" + treatment + "/" + regulation + "_gene_output.csv", mode="a", index=False, header=False)

def nGRE_to_regex(nGRE_sequence):
	nGRE_sequence_regex = ""
	for nucleotide in nGRE_sequence:
		nGRE_sequence_regex += "["
		nGRE_sequence_regex += nucleotide.upper() + nucleotide.lower()
		nGRE_sequence_regex += "]"

	return nGRE_sequence_regex

def main():
	start = time.time()
	csv_parse("ADX_F_DexvsADX_F_Veh", "downregulated")
	end = time.time()

	runtime = end - start
	print("Total runtime: {} seconds".format(runtime))


if __name__ == "__main__":
    main()
