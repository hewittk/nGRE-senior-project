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
