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
	possible_matches = regex.findall("([Cc][Tt][Cc][Cc][TAGCtagc]?[TAGCtagc]?[TAGCtagc]?[Gg][Gg][Aa][Gg][Aa]){e<=2}", gene_sequence)

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
		end_coordinate_regex = regex.findall(r" \d+\)", str(pos_information))[0]
		end_coordinate = int (regex.findall("\d+", str(end_coordinate_regex))[0])

		# compare to perfect nGRE for error counts
		comparison = (regex.search(r"((?e)[Cc][Tt][Cc][Cc][TAGCtagc]?[TAGCtagc]?[TAGCtagc]?[Gg][Gg][Aa][Gg][Aa]){e<=2}", nGRE_sequence))
		mutation_counts = comparison.fuzzy_counts
		mismatch_count = mutation_counts[0]
		insertion_count = mutation_counts[1]
		deletion_count = mutation_counts[2]
		total_mutations = mismatch_count + insertion_count + deletion_count

		# assemble dictionary of nGRE information
		nGRE_information_dict["gene_id"] = gene_id
		nGRE_information_dict["chromosome"] = chromosome
		nGRE_information_dict["nGRE_sequence"] = nGRE_sequence
		nGRE_information_dict["start_site"] = start_coordinate
		nGRE_information_dict["end_site"] = end_coordinate
		nGRE_information_dict["mutations"] = total_mutations
		nGRE_information_dict["mismatch_mutations"] = mismatch_count
		nGRE_information_dict["insertion_mutations"] = insertion_count
		nGRE_information_dict["deletion_mutations"] = deletion_count

		# nGRE_sites = nGRE_sites.append(nGRE_information_dict, ignore_index = True)
		nGRE_list.append(nGRE_information_dict)
		# nGRE_information_dict.clear()

	regex_end = time.time()
	regex_runtime = regex_end - regex_start
	print("Regex search completed, total runtime: {} seconds".format(regex_runtime))

	nGRE_table = pd.DataFrame(columns = ["gene_id", "chromosome", "nGRE_sequence", "start_site", "end_site", "mutations", "mismatch_mutations", "insertion_mutations", "deletion_mutations"])
	for nGRE in nGRE_list:
		nGRE_table = nGRE_table.append(nGRE, ignore_index = True)

	return nGRE_table


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
		potential_nGREs.to_csv("test.csv", mode="a", index=False, header=False)


def main():

	start = time.time()
	csv_parse("ADX_F_DexvsADX_F_Veh", "upregulated")
	end = time.time()

	runtime = end - start
	print("Total runtime: {} seconds".format(runtime))

if __name__ == "__main__":
    main()
