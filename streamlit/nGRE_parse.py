import regex
import streamlit as st
import pandas as pd

def main():
    st.write("My First Streamlit Web App")

    maximum_mutations = st.selectbox(
         'What is the maximum number of mutations in the nGRE consensus sequence that you want to tolerate?',
         ('0', '1', '2', '3'))
    print("Maximum mutations in potential nGREs: ", maximum_mutations)

    user_input = st.text_area("Insert DNA sequence to parse for nGREs")
    st.write(user_input)

    search_string = "([Cc][Tt][Cc][Cc][TAGCtagc]?[TAGCtagc]?[TAGCtagc]?[Gg][Gg][Aa][Gg][Aa]){e<=2}"
    possible_nGREs = regex.findall(search_string, user_input)
    st.write(possible_nGREs)
