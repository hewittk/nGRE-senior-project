import regex
import streamlit as st
import pandas as pd

def main():
    st.write("Parse gene sequences for custom nuclear hormone response element or other subsequence")

    target_element = st.text_area("Insert target element to parse for")
    st.write(target_element)
    element_to_regex(target_element)


    gene_sequence = st.text_area("Insert DNA sequence to search for target element in")
    st.write(gene_sequence)

    maximum_mutations = st.selectbox(
        'What is the maximum number of mutations in the nGRE consensus sequence that you want to tolerate?',
        ('0', '1', '2', '3'))

def element_to_regex(element):
    """Transform user-inputted subsequence into its regex form."""

    # separate subsequence into list elements by parantheses
    subcomponent = ""
    element_components = []

    # append element in groups by parantheses
    for character in element:
        print("character: " + character)
        if character == "(":
            if subcomponent:
                element_components.append(subcomponent)
                subcomponent = ""
            print()
        elif character == ")":
            element_components.append(subcomponent)
            subcomponent = ""
            print()
        if character != "(" and character != ")":
            subcomponent += character
        print(subcomponent)

    # append last subcomponent if not in parantheses
    if (element[len(element)-1]) != ")":
        element_components.append(subcomponent)
    print("Subcomponent list: ", element_components)
