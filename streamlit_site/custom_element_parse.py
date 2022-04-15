import regex
import streamlit as st
import pandas as pd

def main():

    st.write("Parse gene sequences for custom nuclear hormone response element or other subsequence")

    target_element = st.text_area("Insert target element to parse for")
    regex_element = "" # regex equivalent of target element

    gene_sequence = st.text_area("Insert DNA sequence to search for target element in")

    maximum_mutations = int(st.selectbox(
        'What is the maximum number of mutations in the nGRE consensus sequence that you want to tolerate?',
        ('0', '1', '2', '3')))

    mismatch_penalty = int(st.selectbox(
         'What do you want the mismatch penalty (number of points subtracted from a subsequence\'s score for every mismatch mutation) to be?',
         ('1', '2', '3', '0')))

    insertion_penalty = int(st.selectbox(
         'What do you want the insertion penalty (number of points subtracted from a subsequence\'s score for every insertion mutation) to be?',
         ('1', '2', '3', '0')))

    deletion_penalty = int(st.selectbox(
         'What do you want the deletion penalty (number of points subtracted from a subsequence\'s score for every deletion mutation) to be?',
         ('1', '2', '3', '0')))

    if(target_element):
        regex_element = element_to_regex(target_element)

        print("Regex element in main: " + regex_element)

        matches = sequence_search(gene_sequence, regex_element, maximum_mutations)

        matches_df(matches, target_element, regex_element, gene_sequence, maximum_mutations, mismatch_penalty, insertion_penalty, deletion_penalty)

def matches_df(matches, target_element, regex_element, gene_sequence, maximum_mutations, mismatch_penalty, insertion_penalty, deletion_penalty):
    """Score matches to custom element based on amount of mutations and generate dataframe of scored sequences."""

    # obtain minimum length of target element to use as base score
    length = count_length(regex_element)
    print("Length: ", length)

    # remove duplicates in matches list
    nonrepeat_matches = []
    for match in matches:
        if match not in nonrepeat_matches:
            nonrepeat_matches.append(match)
    print("nonrepeat_matches: " + str(nonrepeat_matches))

    # score and add each element to dataframe
    element_table = pd.DataFrame(columns = ["sequence", "start", "end", "score", "mutations", "mismatch", "insertion", "deletion"])

    element_information = {}
    for element in nonrepeat_matches:

        # find number of mutations in sequence
        comparison = (regex.search(r"((?e)" + regex_element + "){e<=" + str(maximum_mutations) + "}", str(element)))
        print(comparison)
        mutation_counts = comparison.fuzzy_counts
        mismatch_count = mutation_counts[0]
        insertion_count = mutation_counts[1]
        deletion_count = mutation_counts[2]
        total_mutations = mismatch_count + insertion_count + deletion_count

        # calculate score
        score = length - insertion_penalty*(insertion_count) - deletion_penalty*(deletion_count) - mismatch_penalty*(mismatch_count)

        # find all start/end locations of element
        element_locations = []
        p = regex.compile(element_to_regex(element))
        for match_location in p.finditer(gene_sequence):
            element_locations.append([match_location.start(), match_location.end()])

        # create match location's row in output dataframe
        for location in element_locations:
            element_information.clear()
            element_information["sequence"] = element
            element_information["start"] = location[0]
            element_information["end"] = location[1]
            element_information["score"] = score
            element_information["mutations"] = total_mutations
            element_information["mismatch"] = mismatch_count
            element_information["insertion"] = insertion_count
            element_information["deletion"] = deletion_count
            element_table = element_table.append(element_information, ignore_index = True)

    st.dataframe(element_table)


def sequence_search(gene_sequence, regex_element, maximum_mutations):
    """Search given sequence for given regex subsequence pattern and return any matches."""

    print("Regex element in sequence search: " + regex_element)

    # process any potential mutations
    if(int(maximum_mutations) > 0):
        regex_element = "((?e)" + regex_element + "){e<=" + str(maximum_mutations) + "}"
    else:
        regex_element = "(" + regex_element + ")"
    print("Regex element after processing " + regex_element)

    # only take first element if regex returns list with second element being spacer nucleotide
    if(type(regex_element) == list):
        regex_element = str(regex_element[0])

    # find and report matches to element within gene sequence
    element_matches = regex.findall(str(regex_element), gene_sequence)
    print("Element matches: " + str(element_matches))

    # only take first element if regex returns tuple with second element being spacer nucleotide
    for index in range(len(element_matches)):
        if type(element_matches[index]) == tuple:
            element_matches[index] = str(element_matches[index][0])

    return element_matches

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

    print("element: " + element)
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
        print("Counted element: " + regex_sequence[i])
        # add to length for each element group contained in brackets
        if(regex_sequence[i] == "]"):
            if((not(i+1 >= len(regex_sequence))) and (regex_sequence[i+1] != "?")):
                length += 1
            elif (i+1 >= len(regex_sequence)):
                length += 1

    # account for any zeros in regular expression
    elements_before_zero = 0
    total_elements_before_zeroes = 0
    for i in range(len(regex_sequence)):
        if(regex_sequence[i] == "0"):
            print("Zero element found, index: " + str(i))
            for j in range(i, 0, -1):
                print("regex_sequence[j]: " + regex_sequence[j])
                if(regex_sequence[j] == "]"):
                    elements_before_zero += 1
                if(regex_sequence[j] == "("):
                    break
        total_elements_before_zeroes += elements_before_zero
        elements_before_zero = 0
    print("Total elements before zeroes: " + str(total_elements_before_zeroes))
    length -= total_elements_before_zeroes

    return length
