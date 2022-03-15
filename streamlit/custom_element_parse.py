import regex
import streamlit as st
import pandas as pd

def main():
    st.write("Parse gene sequences for custom nuclear hormone response element or other subsequence")

    regex_element = ""
    target_element = st.text_area("Insert target element to parse for")
    if(target_element):
        st.write(target_element)
        regex_element = element_to_regex(target_element)

    gene_sequence = st.text_area("Insert DNA sequence to search for target element in")
    st.write(gene_sequence)

    maximum_mutations = st.selectbox(
        'What is the maximum number of mutations in the nGRE consensus sequence that you want to tolerate?',
        ('0', '1', '2', '3'))

    sequence_search(gene_sequence, regex_element, maximum_mutations)

def sequence_search(gene_sequence, regex_element, maximum_mutations):
    """Search given sequence for given regex subsequence pattern and return any matches."""

    # put regex element into python regex package's processing format
    if(int(maximum_mutations) > 0): # process any potential mutations
        regex_element = "((?e)" + regex_element + "){e<=" + str(maximum_mutations) + "}"
    else:
        regex_element = "(" + regex_element + ")"

    element_matches = regex.findall(str(regex_element), gene_sequence)
    print(element_matches)
    st.write(element_matches)

def iupac_handling(iupac_code, position, str):
    """Convert IUPAC codes to nucleotide possibilities that they symbolize."""

    if(iupac_code == "Y"):
        return str[:(position)] + "CcTt" + str[(position+2):]
    if(iupac_code == "S"):
        return str[:(position)] + "GgCc" + str[(position+2):]
    if(iupac_code == "W"):
        return str[:(position)] + "AaTt" + str[(position+2):]
    if(iupac_code == "R"):
        return str[:(position)] + "AaGg" + str[(position+2):]
    if(iupac_code == "M"):
        return str[:(position)] + "AaCc" + str[(position+2):]
    if(iupac_code == "K"):
        return str[:(position)] + "GgTt" + str[(position+2):]
    if(iupac_code == "B"):
        return str[:(position)] + "CcGgTt" + str[(position+2):]
    if(iupac_code == "D"):
        return str[:(position)] + "AaGgTt" + str[(position+2):]
    if(iupac_code == "H"):
        return str[:(position)] + "AaCcTt" + str[(position+2):]
    if(iupac_code == "V"):
        return str[:(position)] + "AaGgCc" + str[(position+2):]
    if(iupac_code == "N"):
        return str[:(position)] + "AaCcTtGg" + str[(position+2):]

    return "None" # if letter submitted not found


def nucleotide_bracket(sequence_component):
    """Put each nucleotide letter into bracket with letter of other case."""

    print("sequence_component: ", sequence_component)
    expanded_sequence_component = ""
    expansion = ""

    for character in sequence_component:
        if character.isalpha():
            if character == character.upper():
                expansion = "[" + character + character.lower() + "]"
            if character == character.lower():
                expansion = "[" + character.upper() + character + "]"
            expanded_sequence_component += expansion
        else:
            expanded_sequence_component += character

    print("expanded_sequence_component: ", expanded_sequence_component)

    return expanded_sequence_component

def element_to_regex(element):
    """Transform user-inputted subsequence into its regex form."""

    # separate subsequence of interest into list elements by parantheses
    subcomponent = ""
    element_components = []

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

    if (element[len(element)-1]) != ")": # append last subcomponent if not in parantheses
        element_components.append(subcomponent)
    print("Subcomponent list: ", element_components)

    # cast subsequence into its regular expression
    previous_component = element_components[0]
    regex_string = ""
    for i in range(1, len(element_components)):

        if element_components[i][0].isnumeric():
            print("Numeric component found: " + element_components[i])
            if previous_component.isalpha(): # check if previous string is all letters
                if "-" in element_components[i]:
                    component_split = element_components[i].split("-")
                    print("Component split: ", component_split)
                    previous_component = "(" + nucleotide_bracket(previous_component) + ")"
                    previous_component += "{" + component_split[0] + "," + component_split[1] + "}"
                elif (element_components[i].isdigit()):
                    previous_component = "(" + nucleotide_bracket(previous_component) + ")"
                    previous_component += "{" + element_components[i] + "}"
                else:
                    st.write("Error: Repitition number or number range in parantheses is not preceded by nucleotides")
        else:
            previous_component = nucleotide_bracket(previous_component)

        if (bool(regex.search('[a-zA-Z]', previous_component))): # if previous string contains nucleotides
            regex_string += previous_component

        previous_component = element_components[i] # increment

    if(element_components[len(element_components)-1].isalpha()):
        regex_string += nucleotide_bracket(element_components[len(element_components)-1])

    print()
    print(regex_string)


    # handle any IUPAC codes in sequence
    iupac_codes = ["N", "S", "W", "Y", "R", "M", "K", "B", "D", "V", "H"]
    i = 0 # incrementer
    prev_length = len(regex_string)
    while (i < len(regex_string)):
        if(regex_string[i] in iupac_codes):
            regex_string = iupac_handling(regex_string[i], i, regex_string)
            i += ((len(regex_string) - prev_length) + 1)
        else:
            i += 1

    print(regex_string)
    print("--------------")
    print()

    return regex_string

def count_length(regex_sequence):
    """Count all mandatory nucleotides to find and return minimum length of sequence."""
    length = 0
    for i in range(len(regex_sequence)):
        if(regex_sequence[i] == "]"):
            if((not(i+1 > len(regex_sequence))) and (regex_sequence[i+1] != "?")):
                length += 1
    return length
