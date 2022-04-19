import streamlit as st

def main():
    st.title("About Tool")

    st.write("This tool aims to allow biologists to detect matches to the nGRE consensus sequence or other subsequences containing variability within gene sequences.")

    st.header("Mutations and mutation penalties")
    st.write("A user can set the maximum number of mutations that they would like to return matches to their interest sequence containing. This number is set by the drop-down menu following the question \"What is the maximum number of mutations in the nGRE consensus sequence/element sequence that you want to tolerate?\"")

    st.write("Penalties for different types of mutations specified by the user are then used to determine a score for each match to the element of interest. An element's base `score` is the minimum length of the sequence, such as 9 for CTCC(n)(0-2)GGAGA, 6 for a custom element of GTG(A)(0-1)TGC, and 7 for a custom element of GGA(C)(1-2)TTG. The `penalty` for a given type of mutation (mismatch, insertion, deletion) is then subtracted from a match's score each time that the given mutation appears in a match to the element. Users determine a penalty for a given mutation in the drop-down menu following the question \"What do you want the [mutation type] penalty (number of points subtracted from a subsequence's score for every [mutation type] mutation) to be?\"")

    st.header("Representing custom elements")
    st.write("On the custom element parsing page, users can input custom subsequences to parse for within given gene sequences. These custom elements can have variability in the form of variable length components and in the form of IUPAC codes representing a base that could consist of one of two or more nucleotides.")

    st.write("Variable length components can be represented using parantheses containing the nucleotide(s) having variable numbers of repeats followed by parantheses indicating the range of times that component could repeat. Examples of elements containing variable length components: ")
    st.write("* CTA(C)(0-3)GT")
    st.write("* GG(A)(2-4)CTC")
    st.write("* ATG(CA)(1-4)GGA")
    st.write("* GA(N)(0-2)CTC")
    st.write("* C(Y)(3-5)ATG")
    st.write("* ACGA(GK)(3-5)")

    st.write("") # spacer

    st.image("IUPAC.png")
    st.write("Bases that can be represented by one of two or more nucleotides can by represented in custom elements in the form of International Union of Pure and Applied Chemistry (IUPAC) codes as listed above. This image was retrieved from the George Mason University paper [Effective Automated Feature Construction and Selection for Classification of Biological Sequences](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4102475/) in a reformatted version of the original IUPAC codes table published by the [Nomenclature Committee of the International Union of Biochemistry](https://pubmed.ncbi.nlm.nih.gov/2417239/). A cropped version of the original image was redistributed on this page via the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/)")
