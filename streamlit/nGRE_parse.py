import regex
import streamlit as st
import pandas as pd

def main():
    st.write("Parse gene sequences for nGREs")

    maximum_mutations = st.selectbox(
         'What is the maximum number of mutations in the nGRE consensus sequence that you want to tolerate?',
         ('0', '1', '2', '3'))
    print("Maximum mutations in potential nGREs: ", maximum_mutations)

    gene_sequence = st.text_area("Insert DNA sequence to parse for nGREs")
    st.write(gene_sequence)

    search_string = "((?e)[Cc][Tt][Cc][Cc][TAGCtagc]?[TAGCtagc]?[TAGCtagc]?[Gg][Gg][Aa][Gg][Aa]){e<=" + maximum_mutations + "}"
    possible_nGREs = regex.findall(search_string, gene_sequence)
    st.write(possible_nGREs)

    nGRE_info_dataframe(possible_nGREs, gene_sequence)


def nGRE_info_dataframe(possible_nGREs, gene_sequence):
    nGRE_table = pd.DataFrame(columns = ["sequence", "start", "end"])

    nGRE_information = {}
    for nGRE_sequence in possible_nGREs:
        nGRE_information["sequence"] = nGRE_sequence

        pos_information = regex.search(nGRE_sequence, gene_sequence)

        start_coordinate_regex = regex.findall(r"\(\d+,", str(pos_information))[0]
        st.write("start_coordinate_regex: " + start_coordinate_regex)
        start_coordinate = int(regex.findall("\d+", str(start_coordinate_regex))[0])
        st.write("start_coordinate: " + str(start_coordinate))
        nGRE_information["start"] = start_coordinate

        end_coordinate_regex = regex.findall(r" \d+\)", str(pos_information))[0]
        st.write("end_coordinate_regex: " + end_coordinate_regex)
        end_coordinate = int(regex.findall("\d+", str(end_coordinate_regex))[0])
        st.write("end_coordinate: " + str(end_coordinate))
        nGRE_information["end"] = end_coordinate

        nGRE_table = nGRE_table.append(nGRE_information, ignore_index = True)

    st.dataframe(nGRE_table)
