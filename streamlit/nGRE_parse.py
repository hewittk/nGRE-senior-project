import regex
import streamlit as st
import pandas as pd
maximum_mutations = 0

def main():
    global maximum_mutations
    st.write("Parse gene sequences for nGREs")

    gene_sequence = st.text_area("Insert DNA sequence to parse for nGREs")
    st.write(gene_sequence)

    maximum_mutations = st.selectbox(
         'What is the maximum number of mutations in the nGRE consensus sequence that you want to tolerate?',
         ('0', '1', '2', '3'))
    print("Maximum mutations in potential nGREs: ", maximum_mutations)

    mismatch_penalty = st.selectbox(
         'What do you want the mismatch penalty (number of points subtracted from a subsequence\'s score for every mismatch mutation) to be?',
         ('1', '2', '3', '0'))

    insertion_penalty = st.selectbox(
         'What do you want the insertion penalty (number of points subtracted from a subsequence\'s score for every insertion mutation) to be?',
         ('1', '2', '3', '0'))

    deletion_penalty = st.selectbox(
         'What do you want the deletion penalty (number of points subtracted from a subsequence\'s score for every deletion mutation) to be?',
         ('1', '2', '3', '0'))

    search_string = "((?e)[Cc][Tt][Cc][Cc][TAGCtagc]?[TAGCtagc]?[TAGCtagc]?[Gg][Gg][Aa][Gg][Aa]){e<=" + maximum_mutations + "}"
    possible_nGREs = regex.findall(search_string, gene_sequence)
    st.write(possible_nGREs)

    nGRE_info_dataframe(possible_nGREs, gene_sequence)


def nGRE_info_dataframe(possible_nGREs, gene_sequence):
    nGRE_table = pd.DataFrame(columns = ["sequence", "start", "end", "mutations"])

    nGRE_information = {}
    for nGRE_sequence in possible_nGREs:

        nGRE_information["sequence"] = nGRE_sequence
        pos_information = regex.search(nGRE_sequence, gene_sequence)

        # find and extract start coordinate from pos_information using regex
        start_coordinate_regex = regex.findall(r"\(\d+,", str(pos_information))[0]
        start_coordinate = int(regex.findall("\d+", str(start_coordinate_regex))[0])
        nGRE_information["start"] = start_coordinate

        # find and extract end coordinate from pos_information using regex
        end_coordinate_regex = regex.findall(r" \d+\)", str(pos_information))[0]
        end_coordinate = int(regex.findall("\d+", str(end_coordinate_regex))[0])
        nGRE_information["end"] = end_coordinate

        comparison = (regex.search(r"((?e)[Cc][Tt][Cc][Cc][TAGCtagc]?[TAGCtagc]?[TAGCtagc]?[Gg][Gg][Aa][Gg][Aa]){e<=" + maximum_mutations + "}", nGRE_sequence))
        print(comparison)
        mutation_counts = comparison.fuzzy_counts
        mismatch_count = mutation_counts[0]
        insertion_count = mutation_counts[1]
        deletion_count = mutation_counts[2]
        total_mutations = mismatch_count + insertion_count + deletion_count
        nGRE_information["mutations"] = total_mutations

        nGRE_table = nGRE_table.append(nGRE_information, ignore_index = True)

    st.dataframe(nGRE_table)
